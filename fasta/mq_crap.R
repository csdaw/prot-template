library(curl)
library(httr)
library(magrittr)
## also uses functions from:
# Biostrings
# data.table

# define the urls we need
BASE <- "https://www.uniprot.org"
TOOL_ENDPOINT <- "/uploadlists/" # REST endpoint
KB_ENDPOINT <- "/uniprot/" # UniProt website search endpoint
UPARC_ENDPOINT <- "/uniparc/" # UniParc website search endpoint

# MaxQuant, Max Plank Institute of Biochemistry (link not always working)
# See http://coxdocs.org/doku.php?id=maxquant:start_downloads.htm&s[]=contaminants

if (!file.exists("fasta/mq_crap_original.fasta")) {
  curl::curl_download(
    url = "http://lotus1.gwdg.de/mpg/mmbc/maxquant_input.nsf/7994124a4298328fc125748d0048fee2/$FILE/contaminants.fasta", 
    destfile = "fasta/mq_crap_original.fasta"
  )
}

# load MaxQuant cRAP fasta
crap_mq <- Biostrings::readAAStringSet("fasta/mq_crap_original.fasta")

# extract all accessions (UniProt or 'other') from fasta headers
accessions <- regexec("^[^ ]+", names(crap_mq)) %>% 
  regmatches(names(crap_mq), .) %>% 
  unlist()

# extract UniProt and non-Uniprot accessions
accessions_up <- grep("^[QPOA][A-Z,0-9]{5}", accessions, value = TRUE)
accessions_other <- grep("^[QPOA][A-Z,0-9]{5}", accessions, invert = TRUE, value = TRUE)

## Here we map canonical (non-isoform) UniProt accessions from the original 
## MaxQuant fasta to modern UniProt accessions. Note that the mapping is not
## 1:1 as the original fasta contains secondary and obsolete accessions.

# extract canonical UniProt accessions
accessions_up_can <- grep("-[2-9]", accessions_up, invert = TRUE, value = TRUE)

payload <- list(
  query = paste(accessions_up_can, collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "tab"
)

response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload)

if (response$status_code == 200) {
  mapping_up_can <- data.table::fread(content(response))
  print(paste("Input length:", length(accessions_up_can)))
  print(paste("Output length:", length(mapping_up_can$To)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

# obtain FASTA for modern, canonical (non-isoform) UniProt accessions
payload <- list(
  query = paste(unique(mapping_up_can$To), collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  up_can_fasta <- Biostrings::readAAStringSet(tmp)
  print(paste("Input length:", length(unique(mapping_up_can$To))))
  print(paste("Output FASTA length:", length(up_can_fasta)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

## Note that there are some accessions that don't map to a FASTA,
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

obsolete <- mapping_up_can$To[!mapping_up_can$To %in% extracted_accessions]

## Now we map isoform UniProt accessions from the original 
## MaxQuant fasta to modern UniProt accessions.

# extract isoform UniProt accessions
accessions_up_iso <- grep("-[2-9]", accessions_up, value = TRUE)

# obtain FASTA for modern, canonical (non-isoform) UniProt accessions
payload <- list(
  query = paste(accessions_up_iso, collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  up_iso_fasta <- Biostrings::readAAStringSet(tmp)
  print(paste("Input length:", length(accessions_up_iso)))
  print(paste("Output FASTA length:", length(up_iso_fasta)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}


## Now the non-UniProt accessions are dealt with.
# extract refseq accessions
accessions_refseq <- regexec("(?<=REFSEQ:)[A-Z,0-9,_]+", accessions_other, perl = TRUE) %>% 
  regmatches(accessions_other, .) %>% 
  unlist()

# search UniParc using the refseq accessions to obtain protein fasta
payload <- list(
  query = paste(accessions_refseq, collapse = " OR "),
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, UPARC_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  refseq_fasta <- Biostrings::readAAStringSet(tmp)
  print(paste("Input length:", length(accessions_refseq)))
  print(paste("Output FASTA length:", length(refseq_fasta)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

## The mapping of refseq accessions to UniParc sequences was not 1:1.
## Therefore we compare these sequences to the original MaxQuant fasta to get
## the UniParc sequences of interest.

# keep only sequences of interest
refseq_fasta <- Biostrings::intersect(refseq_fasta, crap_mq)

# extract h-inv accessions
accessions_hinv <- regexec("(?<=H-INV:)[A-Z,0-9]+", accessions_other, perl = TRUE) %>% 
  regmatches(accessions_other, .) %>% 
  unlist()

# extract ensembl accessions
accessions_ens <- regexec("(?<=ENSEMBL:)[A-Z,0-9,_]+", accessions_other, perl = TRUE) %>% 
  regmatches(accessions_other, .) %>% 
  unlist()

# extract anything else non-UniProt
accessions_manual <- grep("(ENSEMBL|REFSEQ|H-INV)", accessions_other, invert = TRUE, perl = TRUE)
