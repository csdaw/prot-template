###-------------------------###
#### What this script does ####
###-------------------------###

## Author: Charlotte Dawson (csdaw@outlook.com)

## Original FASTA source: MaxQuant, Max Plank Institute of Biochemistry
## http://coxdocs.org/doku.php?id=maxquant:start_downloads.htm&s[]=contaminants

## Description of MaxQuant contaminants fixed FASTA (x sequences).
## Many of these sequences are old/obsolete. Therefore, I have written this 
## script to generate a FASTA file with the most recent sequences for these 
## proteins (i.e. from the latest UniProt release)

## Input: mq_contaminants_fixed.fasta

## Output: mq_contaminants_XXXX_XX.fasta where XXXX_XX is the name of the 
## latest UniProt release.

###----------------###
#### Dependencies ####
###----------------###

library(httr)
library(magrittr)

## This script also calls functions from the following packages:
# Biostrings

# define some urls for GET requests
BASE <- "https://www.uniprot.org"
TOOL_ENDPOINT <- "/uploadlists/" # REST endpoint
KB_ENDPOINT <- "/uniprot/" # UniProt website search endpoint

###------------------------------------------###
#### Load MaxQuant contaminants fixed FASTA ####
###------------------------------------------###

# load my MaxQuant contaminants fixed FASTA
crap_fixed <- Biostrings::readAAStringSet("fasta/mq_contaminants_fixed.fasta")

###----------------------------------------------###
#### Extract accessions from outdated sequences ####
###----------------------------------------------###
# extract sequences with UniParc headers (i.e. outdated sequences)
crap_uniprot <- crap_fixed[grep("^(sp|tr)", names(crap_fixed))]
crap_uniparc <- crap_fixed[grep("^(sp|tr)", names(crap_fixed), invert = TRUE)]

# extract accessions from the headers of UniParc headers
accessions <- regexec("(?<=formerly\\=).+", names(crap_uniparc), perl = TRUE) %>% 
  regmatches(names(crap_uniparc), .) %>% 
  unlist()

# extract UniProt and non-Uniprot accessions
# UniProt accessions (length: 56)
accessions_up <- grep(
  "^[QPOA][A-Z,0-9]{5}", 
  accessions, 
  value = TRUE
)

# Non-UniProt accessions (length: 29)
accessions_non_up <- grep(
  "^[QPOA][A-Z,0-9]{5}", 
  accessions, 
  invert = TRUE, 
  value = TRUE
)

###--------------------------------------------------------###
#### Deal with outdated sequences with UniProt accessions ####
###--------------------------------------------------------###

# see which accessions map to a modern sequence
payload <- list(
  query = paste(accessions_up, collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  crap_uniparc_up_fasta <- Biostrings::readAAStringSet(tmp)
  message(paste("Input length:", length(accessions_up)))
  message(paste("Output FASTA length:", length(crap_uniparc_up_fasta)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

length(crap_uniparc_up_fasta) # 52/56

# which are the 4/56 left behind?
payload <- list(
  query = paste(accessions_up, collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "tab",
  columms = paste(c("entry_name", "entry"), collapse = ",")
)

response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload)

if (response$status_code == 200) {
  crap_uniparc_up_mapping <- data.table::fread(content(response, encoding = "UTF-8"))
  message(paste("Input length:", length(accessions_up)))
  message(paste("Output table length:", nrow(crap_uniparc_up_mapping)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

# 2/4 are duplicates and have been merged
any(duplicated(crap_uniparc_up_mapping$To))
crap_uniparc_up_mapping[duplicated(crap_uniparc_up_mapping$To) | 
                             duplicated(crap_uniparc_up_mapping$To, fromLast = TRUE), ]

# 2/4 have been deleted
accessions_deleted <- regexec("(?<=sp\\||tr\\|)[A-Z0-9]{6}", names(crap_uniparc_up_fasta), perl = TRUE) %>% 
  regmatches(names(crap_uniparc_up_fasta), .) %>% 
  unlist() %>% 
  setdiff(crap_uniparc_up_mapping$To, .)

###------------------------------------------------------------###
#### Deal with outdated sequences with non-UniProt accessions ####
###------------------------------------------------------------###
