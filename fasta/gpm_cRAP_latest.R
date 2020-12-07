###-------------------------###
#### What this script does ####
###-------------------------###

## Author: Charlotte Dawson (csdaw@outlook.com)

## Original FASTA source: Ron Beavis group, the Global Proteome Machine (GPM)
## https://www.thegpm.org/crap/ and ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta

## Description of GPM cRAP fixed FASTA (x sequences).
## Many of these sequences are old/obselete. Therefore, I have written this 
## script to generate a FASTA file with the most recent sequences for these 
## proteins (i.e. from the latest UniProt release)

## Input: gpm_cRAP_fixed.fasta

## Output: gpm_XXXX_XX.fasta where XXXX_XX is the name of the latest UniProt 
## release.

###----------------###
#### Dependencies ####
###----------------###

library(httr)
library(magrittr)

## This script also calls functions from the following packages:
# Biostrings

## Define some urls for GET requests
BASE <- "https://www.uniprot.org"
KB_ENDPOINT <- "/uniprot/" # UniProt website search endpoint

###---------------------------###
#### Download GPM cRAP FASTA ####
###---------------------------###

# load my GPM cRAP fixed FASTA
crap_fixed <- Biostrings::readAAStringSet("fasta/gpm_crap_fixed.fasta")

# separate the "good" sequences from the "bad"
crap_good <- crap_fixed[grep("^sp", names(crap_fixed))]
crap_bad <- crap_fixed[grep("^sp", names(crap_fixed), invert = TRUE)]

# sort crap_bad by entry names
sort_order <- regexec(
  "(?<=formerly=)[A-Z,0-9]+_[A-Z]+", 
  names(crap_bad), 
  perl = TRUE
) %>%
  regmatches(names(crap_bad), .) %>% 
  unlist() %>% 
  order()

crap_bad <- crap_bad[sort_order]

###--------------------------------------------------###
#### Extract UniProt entry names from FASTA headers ####
###--------------------------------------------------###

## First, we extract the entry names (aka mnemonics) from the FASTA headers.
## Note: it is not good to use entry names as they are not stable identifiers.
## Rather it is important always to use UniProt accessions AND sequence versions
## in FASTA headers (or UniParc accessions).
entry_names <- regexec("[A-Z,0-9]+_[A-Z]+", names(crap_bad)) %>% 
  regmatches(names(crap_bad), .) %>% 
  unlist()

payload <- list(
  query = paste(entry_names, collapse = " OR "),
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, KB_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  crap_bad_new <- Biostrings::readAAStringSet(tmp)
  message(paste("Input length:", length(crap_bad)))
  message(paste("Output FASTA length:", length(crap_bad_new)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

# sort crap_bad_new alphabetically by entry names
sort_order <- regexec(
  "(?<=sp\\|[A-Z,0-9]{6}\\|)[A-Z,0-9]+_[A-Z]+", 
  names(crap_bad_new), 
  perl = TRUE
) %>%
  regmatches(names(crap_bad_new), .) %>% 
  unlist() %>% 
  order()

crap_bad_new <- crap_bad_new[sort_order]

# compare crap_bad and crap_bad_new headers
# names(crap_bad)
# names(crap_bad_new)

###----------###
#### Output ####
###----------###

# combine crap_good and crap_bad_new
output <- c(
  crap_good, 
  crap_bad_new
)

# sort alphabetically by entry names
sort_order <- regexec(
  "(?<=sp\\|[A-Z,0-9]{6}\\|)[A-Z,0-9]+_[A-Z]+", 
  names(output), 
  perl = TRUE
) %>%
  regmatches(names(output), .) %>% 
  unlist() %>% 
  order()

output <- output[sort_order]

# get the current UniProt release
cur_release <- response$headers$`x-uniprot-release`

# write updated GPM cRAP FASTA file
Biostrings::writeXStringSet(
  output, 
  filepath = paste0("fasta/gpm_cRAP_", cur_release, ".fasta")
)
