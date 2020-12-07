###-------------------------###
#### What this script does ####
###-------------------------###

## Author: Charlotte Dawson (csdaw@outlook.com)

## Description of GPM FASTA (x sequences).
## This FASTA is relatively old and the headers contain several incorrect and 
## out-of-date accessions. 

## I have written this script to generate a FASTA file with (almost) identical
## sequences to the original GPM cRAP.fasta. However, these 
## sequences are annotated with headers that contain modern, working
## accessions.

## Original FASTA source: Ron Beavis group, the Global Proteome Machine (GPM)
## https://www.thegpm.org/crap/ and ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta

###----------------###
#### Dependencies ####
###----------------###

library(httr)
library(magrittr)
library(rvest)

## This script also calls functions from the following packages:
# Biostrings
# data.table
# curl

## Define some urls for GET requests
BASE <- "https://www.uniprot.org"
TOOL_ENDPOINT <- "/uploadlists/" # REST endpoint
KB_ENDPOINT <- "/uniprot/" # UniProt website search endpoint
UPARC_ENDPOINT <- "/uniparc/" # UniParc website search endpoint

###---------------------------###
#### Download GPM cRAP FASTA ####
###---------------------------###

# download original GPM cRAP FASTA if it does not already exist
if (!file.exists("fasta/gpm_cRAP_original.fasta")) {
  curl::curl_download(
    url = "ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta", 
    destfile = "fasta/gpm_cRAP_original.fasta"
  )
}

# load GPM cRAP fasta
gpm_crap <- Biostrings::readAAStringSet("fasta/gpm_crap_original.fasta")

###--------------------------------------------------###
#### Extract UniProt entry names from FASTA headers ####
###--------------------------------------------------###

## First, we extract the entry names (aka mnemonics) from the FASTA headers.
## Note: it is not good to use entry names as they are not stable identifiers.
## Rather it is important always to use UniProt accessions AND sequence versions
## in FASTA headers (or UniParc accessions).
entry_names <- regexec("[A-Z,0-9]+_[A-Z]+", names(gpm_crap)) %>% 
  regmatches(names(gpm_crap), .) %>% 
  unlist()

###------------------------------------------------###
#### Extract UniProt accessions from cRAP website ####
###------------------------------------------------###

## Next we use the {rvest} package to scrape the corresponding UniProt 
## accessions from the GPM website
gpm_url <- "https://www.thegpm.org/crap/"

# parse the tables on the website into a list of data.frames
gpm_crap_list <- read_html(gpm_url) %>% 
  html_nodes("table") %>% 
  html_table(fill = TRUE)
 
# combine the separate tables into one
gpm_crap_table <- do.call(rbind, gpm_crap_list[3:7]) %>% 
  `colnames<-`(c("n", "id", "description", "reason")) %>% 
  dplyr::filter(!is.na(n))

# extract UniProt accessions from the description column
gpm_crap_table$accession <- regexec(
  "(?<=\\()[A-Z,0-9]{6}", 
  gpm_crap_table$description, 
  perl = TRUE
) %>% 
  regmatches(gpm_crap_table$description, .) %>% 
  unlist()

## For some reason the GPM website cRAP table (115 entries) is shorter 
## than the GPM cRAP fasta (116 entries). Lets find the extra protein.

# check GPM cRAP table against entry_names we extracted earlier
extra_protein <- entry_names[!entry_names %in% gpm_crap_table$id]

## It is some E. coli protein, which we can quickly get the accession for and
## append it to the gpm_crap_table.
payload <- list(
  query = extra_protein,
  from = "ACC+ID",
  to = "ACC",
  format = "list"
)

response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload)

# overwrite entry_name of extra protein with UniProt accession
extra_protein <- strsplit(content(response), split = "\n") %>% 
  unlist()

###----------------------------------------###
#### Obtain FASTA with updated accessions ####
###----------------------------------------###

# Use UniProt accessions to obtain sequences from UniParc
payload <- list(
  query = paste(c(gpm_crap_table$accession, extra_protein), collapse = " OR "),
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, UPARC_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  gpm_uparc <- Biostrings::readAAStringSet(tmp)
  message(paste("Input length:", length(gpm_crap_table$accession)))
  message(paste("Output FASTA length:", length(gpm_uparc)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

# compare new UniParc FASTA to old cRAP FASTA (should be 116/116)
length(Biostrings::intersect(gpm_uparc, gpm_crap)) == 116

## 115/116 are good because the UniParc sequences match the protein
## sequences in the cRAP FASTA.
gpm_uparc_good <- Biostrings::intersect(gpm_uparc, gpm_crap)
length(gpm_uparc_good)

## 1/116 are bad because the UniParc sequences do not match the protein
## sequences in the cRAP FASTA.
gpm_uparc_bad <- Biostrings::setdiff(gpm_crap, gpm_uparc)
length(gpm_uparc_bad)

## Which protein is this? It is P00883 (UniParc accession UPI000016C534).
gpm_uparc_bad

## The reason the sequence doesn't match is because in the gpm_crap FASTA it 
## is missing the final 4 amino acids for some reason. Therefore the new UniParc
## sequence is the correct one.
gpm_uparc <- c(
  gpm_uparc_good,
  gpm_uparc[180] # manually add P0083
)

## This is why it is important that FASTA files retain the full UniProt headers
## with the accession and sequence version number (SV) intact! Because the 
## original GPM FASTA did not include the sequence versions we can't easily
## replicate it.

# use UniParc accessions to get some modern UniProtKB headers if available
payload <- list(
  query = paste(
    gsub("(?<=UPI[0-9,A-Z]{10}).*", "", names(gpm_uparc), perl = TRUE), 
    collapse = " "
  ),
  from = "UPARC",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  gpm_new <- Biostrings::readAAStringSet(tmp)
  message(paste("Input length:", length(gpm_uparc)))
  message(paste("Output FASTA length:", length(gpm_new)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

# keep only sequences from: HUMAN, MOUSE, BOVIN, PIG, SHEEP, HORSE, CHICK, 
# RABIT, HEVBR, AEQVI, LYSEN, ECOLX, STAAU, YEAST, GRIFR, SCVLA, ECOLI
gpm_new <- gpm_new[grep(
  paste(c("HUMAN", "MOUSE", "BOVIN", "PIG", "SHEEP", "HORSE", "CHICK",
          "RABIT", "YEAST", "ECOLI", "ECOLX", "HEVBR", "AEQVI", "LYSEN", 
          "STAAU", "GRIFR", "SCVLA"), collapse = "|"),
  names(gpm_new)
)]

# combine the 116 sequences with updated headers into a single FASTA
gpm_new_fasta <- c(
  Biostrings::intersect(gpm_new, gpm_crap), 
  Biostrings::setdiff(gpm_uparc, gpm_new),
  gpm_new[133] # manually add P0083
)

# length of new FASTA should equal length of cRAP FASTA
length(gpm_new_fasta) == length(gpm_crap)

# check intersection of sequences, should be 115/116
length(Biostrings::intersect(gpm_new_fasta, gpm_crap)) == 115

# sequence that doesn't match should be P00883 (UniParc accession UPI000016C534)
Biostrings::setdiff(gpm_new_fasta, gpm_crap)

###----------------------------###
#### Compare input and output ####
###----------------------------###

# sort input and output FASTAs by sequence lengths
gpm_crap <- gpm_crap[order(lengths(gpm_crap))]
gpm_new_fasta <- gpm_new_fasta[order(lengths(gpm_new_fasta))]

# create a data.frame to compare the input and output
df <- data.frame(
  input_headers = names(gpm_crap),
  input_lengths = lengths(gpm_crap),
  output_headers = names(gpm_new_fasta),
  output_lengths = lengths(gpm_new_fasta)
)

# manually fix order of some things
df[c(21, 20, 27, 26, 30, 29, 64, 63, 
     71, 70, 76, 75, 80, 79, 82, 81), 3] <- df[c(20, 21, 26, 27, 29, 30, 63, 64, 
                                                 70, 71, 75, 76, 79, 80, 81, 82), 3]

# append former headers to UniParc headers
former_headers <- regexec("[A-Z,0-9]+_[A-Z]+", df$input_headers) %>% 
  regmatches(df$input_headers, .) %>% 
  unlist()

df[grep("UPI", df$output_headers), 3] <- paste0(
  df[grep("UPI", df$output_headers), 3],
  " formerly=",
  former_headers[grep("UPI", df$output_headers)]
)

# update headers in new FASTA
names(gpm_new_fasta) <- df$output_headers

###----------###
#### Output ####
###----------###

## Finally, we write our updated FASTA to a file with the current UniProt
## release name appended.

# sort according to entry name (alphabetical)
sort_order <- regexec(
  "(?<=sp\\|[A-Z,0-9]{6}\\|)[A-Z,0-9]+_[A-Z]+|UPI[0-9,A-Z]{10}", 
  names(gpm_new_fasta), 
  perl = TRUE
) %>%
  regmatches(names(gpm_new_fasta), .) %>% 
  unlist() %>% 
  order()

gpm_new_fasta <- gpm_new_fasta[sort_order]

# get the current UniProt release
cur_release <- response$headers$`x-uniprot-release`

# write updated GPM cRAP FASTA file
Biostrings::writeXStringSet(
  gpm_new_fasta, 
  filepath = "fasta/gpm_cRAP_fixed.fasta"
)
