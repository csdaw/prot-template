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

## Try out the NCBI Entrez E-utilities
ENTREZ_BASE <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
ENTREZ_SEARCH <- "esearch.fcgi"

payload <- list(
  db = "protein",
  term = paste(refseq_accessions, collapse = ","),
  idtype = "acc"
)

response <- GET(url = paste0(ENTREZ_BASE, ENTREZ_SEARCH), query = payload)
response$headers
response$status_code
XML::xmlParse(content(response)) 

ENTREZ_FETCH <- "efetch.fcgi"

payload <- list(
  db = "protein",
  id = paste(refseq_accessions, collapse = ","),
  retmode = "xml",
  parsed = "true"
)

response <- GET(url = paste0(ENTREZ_BASE, ENTREZ_FETCH), query = payload)
response$headers
response$status_code
xxx <- content(response)
class(xxx)
yyy <- XML::xmlParse(content(response)) 
zzz <- XML::xmlToDataFrame(yyy)

XML::xpathSApply(xxx, "\\LineageEx", XML::xmlValue)

ENTREZ_LINK <- "elink.fcgi?dbfrom=protein"

# extract h-inv accessions
hinv_accessions <- regexec("(?<=H-INV:)[A-Z,0-9]+", accessions_other, perl = TRUE) %>% 
  regmatches(accessions_other, .) %>% 
  unlist()

# extract ensembl accessions
ens_accessions <- regexec("(?<=ENSEMBL:)[A-Z,0-9,_]+", accessions_other, perl = TRUE) %>% 
  regmatches(accessions_other, .) %>% 
  unlist()

# extract anything else non-UniProt
man_accessions <- grep("(ENSEMBL|REFSEQ|H-INV)", accessions_other, invert = TRUE, perl = TRUE)
