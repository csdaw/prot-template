library(curl)
library(httr)
library(magrittr)
library(rentrez)
## also uses functions from:
# Biostrings
# data.table

# define the urls we need
BASE <- "https://www.uniprot.org"
TOOL_ENDPOINT <- "/uploadlists/" # REST endpoint
KB_ENDPOINT <- "/uniprot/" # UniProt website search endpoint
UPARC_ENDPOINT <- "/uniparc/" # UniParc website search endpoint

## MaxQuant, Max Plank Institute of Biochemistry (link not always working)
## See http://coxdocs.org/doku.php?id=maxquant:start_downloads.htm&s[]=contaminants

## This FASTA is also relatively old and not annotated properly. Therefore, we 
## want to get the proteins from this original MaxQuant contaminants FASTA and 
## create a new FASTA with the protein sequences from the most recent UniProt 
## release.

# download original MaxQuant contaminants FASTA if it does not already exist
if (!file.exists("fasta/mq_contaminants_original.fasta")) {
  curl::curl_download(
    url = "http://lotus1.gwdg.de/mpg/mmbc/maxquant_input.nsf/7994124a4298328fc125748d0048fee2/$FILE/contaminants.fasta", 
    destfile = "fasta/mq_contaminants_original.fasta"
  )
}

# load MaxQuant contaminants FASTA
mq_crap <- Biostrings::readAAStringSet("fasta/mq_contaminants_original.fasta")

# extract all accessions (UniProt or 'other') from FASTA headers
accessions <- regexec("^[^ ]+", names(mq_crap)) %>% 
  regmatches(names(mq_crap), .) %>% 
  unlist()

# extract UniProt and non-Uniprot accessions
accessions_up <- grep("^[QPOA][A-Z,0-9]{5}", accessions, value = TRUE)
accessions_other <- grep("^[QPOA][A-Z,0-9]{5}", accessions, invert = TRUE, value = TRUE)

## Here we map canonical (non-isoform) UniProt accessions (209 entries) from the 
## original MaxQuant fasta to modern UniProt accessions.

# extract canonical UniProt accessions
up_can_accessions <- grep("-[2-9]", accessions_up, invert = TRUE, value = TRUE)

payload <- list(
  query = paste(up_can_accessions, collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "tab"
)

response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload)

if (response$status_code == 200) {
  up_can_mapping <- data.table::fread(content(response)) %>% 
    `colnames<-`(c("old_accession", "new_accession"))
  print(paste("Input length:", length(up_can_accessions)))
  print(paste("Output length:", length(unique(up_can_mapping$new_accession))))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

## Note that the mapping is not 1:1 (input: 209 entries, output: 201 entries) as 
## the original contaminants FASTA contains old accessions that now map to the 
## same, modern accession as 'secondary' accessions. IMPORTANTLY, the FASTA
## headers will now include both the UniProt accession and sequence version 
## number (SV).

# obtain FASTA for modern, canonical (non-isoform) UniProt accessions
payload <- list(
  query = paste(unique(up_can_mapping$new_accession), collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  up_can_fasta <- Biostrings::readAAStringSet(tmp)
  print(paste("Input length:", length(unique(up_can_mapping$new_accession))))
  print(paste("Output FASTA length:", length(up_can_fasta)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

## Note that there are some accessions (4 entries) that don't map to a FASTA,
## because they are obsolete and have been removed from UniProt.
## We will leave these sequences out of the final FASTA. 

# identify obsolete accessions
extracted_accessions <- regexec(
  "(?<=sp\\||tr\\|)[QPOA][A-Z,0-9,-]{5,8}", 
  names(up_can_fasta), 
  perl = TRUE
) %>% 
  regmatches(names(up_can_fasta), .) %>% 
  unlist()

obsolete <- up_can_mapping$new_accession[!up_can_mapping$new_accession %in% 
                                           extracted_accessions]

## Next we map the ISOFORM UniProt accessions from the original MaxQuant FASTA 
## to modern UniProt accessions.

# extract isoform UniProt accessions
up_iso_accessions <- grep("-[2-9]", accessions_up, value = TRUE)

# obtain FASTA for modern, canonical (non-isoform) UniProt accessions
payload <- list(
  query = paste(up_iso_accessions, collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  up_iso_fasta <- Biostrings::readAAStringSet(tmp)
  print(paste("Input length:", length(up_iso_accessions)))
  print(paste("Output FASTA length:", length(up_iso_fasta)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

## Now the non-UniProt accessions are dealt with.

# extract refseq accessions
refseq_accessions <- regexec("(?<=REFSEQ:)[A-Z,0-9,_]+", accessions_other, perl = TRUE) %>% 
  regmatches(accessions_other, .) %>% 
  unlist()

## Try rentrez.
# get summary for each refseq accession (aka reference sequence)
bbb <- entrez_summary(db = "protein", id = refseq_accessions)
names(bbb) <- refseq_accessions

# see the status of each record
refseq_statuses <- sapply(bbb, "[", "status") %>% unlist()

## We can see that the records for 4 out of 6 refseq accessions have been 
## "removed as a result of standard genome annotation processing". Therefore 
## we will not include their sequences in our final updated FASTA.

# get suppressed refseq accessions
bad_refseq <- names(refseq_statuses)[grep("suppressed", refseq_statuses)] %>% 
  gsub(".status", "", .)

# remove them from dasdasdasdlasmdlasmdlasdlamsdlamlsdmlasdmalsdmalsdmlasmdmasldmasldasld
bbb <- bbb[-grep("suppressed", refseq_statuses)]

## Also we can see that sequences referred to by the remaining 2 accessions 
## have been replaced. Therefore we will use the modern protein sequences in
## our final updated FASTA.

# get modern refseq accessions
good_refseq <- sapply(bbb, "[", "replacedby") %>% unlist()

## Now we can map these to Uniprot accessions.
payload <- list(
  query = paste(good_refseq, collapse = " "),
  from = "P_REFSEQ_AC",
  to = "ACC",
  format = "tab"
)

response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload)

if (response$status_code == 200) {
  refseq_mapping <- data.table::fread(content(response)) %>% 
    `colnames<-`(c("old_accession", "new_accession"))
  print(paste("Input length:", length(good_refseq)))
  print(paste("Output length:", length(unique(refseq_mapping$new_accession))))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

## And lastly obtain the FASTA sequences for these 2 proteins
payload <- list(
  query = paste(refseq_mapping$new_accession, collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  refseq_fasta <- Biostrings::readAAStringSet(tmp)
  print(paste("Input length:", length(unique(refseq_mapping$new_accession))))
  print(paste("Output FASTA length:", length(refseq_fasta)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

# H-InvDB (an old database of human genes and transcripts
# which hasn't been updated recently, not since 2015!)
# See http://h-invitational.jp/hinv/ahg-db/index.jsp

# extract H-InvDB accessions
hinv_accessions <- regexec("(?<=H-INV:)[A-Z,0-9]+", accessions_other, perl = TRUE) %>% 
  regmatches(accessions_other, .) %>% 
  unlist()

# map them to UniParc accessions
payload <- list(
  query = paste(hinv_accessions, collapse = " OR "),
  format = "tab"
)

response <- GET(url = paste0(BASE, UPARC_ENDPOINT), query = payload)

if (response$status_code == 200) {
  hinv_mapping <- data.table::fread(content(response))
  print(paste("Input length:", length(hinv_accessions)))
  print(paste("Output length:", nrow(hinv_mapping)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

## The 3 H-InvDB accessions have mapped to 6 potential UniParc sequences.
## Lets see exactly which sequences were in the original MaxQuant FASTA.
hinv_lengths <- lengths(mq_crap[grep("H-INV", names(mq_crap))])

hinv_mapping <- hinv_mapping[hinv_mapping$Length %in% hinv_lengths]

## We can see that these sequences have been removed from UniProt. Therefore
## we will not include them in the final updated FASTA.

## Moving onto the Ensembl accessions.
# extract Ensembl accessions
ens_accessions <- regexec("(?<=ENSEMBL:)[A-Z,0-9,_]+", accessions_other, perl = TRUE) %>% 
  regmatches(accessions_other, .) %>% 
  unlist()

# map them to modern UniProt accessions
payload <- list(
  query = paste(ens_accessions, collapse = " "),
  from = "ENSEMBL_PRO_ID",
  to = "ACC",
  format = "tab"
)

response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload)

if (response$status_code == 200) {
  ens_mapping <- data.table::fread(content(response)) %>% 
    `colnames<-`(c("old_accession", "new_accession"))
  print(paste("Input length:", length(ens_accessions)))
  print(paste("Output length:", length(ens_mapping$new_accession)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

# get the fasta for Ensembl accessions that do map
payload <- list(
  query = paste(ens_mapping$new_accession, collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  ens_fasta <- Biostrings::readAAStringSet(tmp)
  print(paste("Input length:", length(unique(ens_mapping$new_accession))))
  print(paste("Output FASTA length:", length(ens_fasta)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

# extract the Ensembl accessions that don't map
ens_bad <- ens_accessions[!ens_accessions %in% ens_mapping$old_accession]

# convert these to UniParc accessions
payload <- list(
  query = paste(ens_bad, collapse = " OR "),
  format = "tab"
)

response <- GET(url = paste0(BASE, UPARC_ENDPOINT), query = payload)

if (response$status_code == 200) {
  ens_bad_uparc <- data.table::fread(content(response))
  print(paste("Input length:", length(ens_bad)))
  print(paste("Output length:", nrow(ens_bad_uparc)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

# then convert these UniParc accessions to UniProtKB
payload <- list(
  query = paste(ens_bad_uparc$Entry, collapse = " "),
  from = "UPARC",
  to = "ACC",
  format = "tab"
)

response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload)

if (response$status_code == 200) {
  ens_bad_mapping <- data.table::fread(content(response)) %>% 
    `colnames<-`(c("old_accession", "new_accession"))
  print(paste("Input length:", nrow(ens_bad_uparc)))
  print(paste("Output length:", length(ens_bad_mapping$new_accession)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

## Only 3 UniProtKB accessions are output. Lets get their sequences.
payload <- list(
  query = paste(ens_bad_mapping$new_accession, collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  ens_bad_fasta <- Biostrings::readAAStringSet(tmp)
  print(paste("Input length:", length(unique(ens_bad_mapping$new_accession))))
  print(paste("Output FASTA length:", length(ens_bad_fasta)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

## Two of the sequences are identical.
ens_bad_fasta[[2]] == ens_bad_fasta[[3]]

## Therefore we will just drop the sequence from the incorrect cow proteome
ens_bad_fasta <- ens_bad_fasta[-2]

## Lastly, we will deal with any other proteins in the original FASTA
# extract anything else non-UniProt
man_accessions <- grep(
  "(ENSEMBL|REFSEQ|H-INV)", 
  accessions_other, 
  invert = TRUE, value = TRUE, perl = TRUE
)

# it is streptavidin which we will just retrieve manually
payload <- list(
  query = "P22629",
  from = "ACC+ID",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  man_fasta <- Biostrings::readAAStringSet(tmp)
  print(paste("Input length:", length("P22629")))
  print(paste("Output FASTA length:", length(man_fasta)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

## Combine all the FASTAs
output <- c(
  up_can_fasta,
  up_iso_fasta,
  refseq_fasta,
  ens_fasta,
  ens_bad_fasta,
  man_fasta
)
