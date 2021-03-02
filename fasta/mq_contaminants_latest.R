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

# obtain sequences in FASTA format
crap_uniparc_non_up <- crap_uniparc[
  grep(paste(accessions_non_up, collapse = "|"), 
       names(crap_uniparc))
]

## One can observe that many of the sequences are "inactive" which means they've 
## been deleted from UniProt and only exist as UniParc sequences now.
names(crap_uniparc_non_up)

## Here we'll BLAST each sequence individually by hand and see if it not already
## represented in the fixed FASTA.
df_non_up <- data.frame(
  headers = names(crap_uniparc_non_up),
  decision = c(
    "replace with UniProtKB Streptavidin P22629",
    "discard, 93.9% identity with existing FASTA sequence P05787",
    "replace, 78.2% identity with P05783",
    "replace, 100% identity with G3N188",
    "replace, 96.8% identity with Q7SIH1",
    "replace, 100% identity with E1BJK2",
    "replace, 100% identity with E1BF81",
    "discard, 98.5% identity with existing FASTA sequence A2I7N3",
    "replace, 89.2% identity with P35445",
    "replace, 98.6% identity with A0A3Q1NGH3 A0A3Q1MCP0 F1MJK3",
    "discard, 97.6% identity with existing FASTA sequence Q3MHN5",
    "replace, 99.1% identity with E1BNR0",
    "replace, 97.3% identity with E1BKY4",
    "replace, 98.9% identity with F1N076",
    "replace, 95.3% identity with A0A3Q1LSF0",
    "discard, 91.7% identity with G3N188",
    "replace, 89.9% identity with E1BCW0",
    "discard, 77.9% identity with A0A3Q1M3L6",
    "replace, 99.2% identity with A0A3Q1M3L6",
    "replace, 94.0% identity with G3N0S9",
    "replace, 70.9% identity with A0A3Q1MHR7",
    "discard, 85.4% identity with existing FASTA sequence P05787",
    "discard, 50.5% identity with E1BIL2",
    "discard, 95.5% identity with existing FASTA sequence P05787",
    "discard, 96.4% identity with existing FASTA sequence A2I7N3",
    "replace, 99.4% identity with F1MVK1",
    "replace, 100% identity with Q61897",
    "replace, 98.7% identity with Q925H3",
    "replace, 100% identity with P13646"
  )
)
df_non_up$accession <- sub('^.* ([[:alnum:]]+)$', '\\1', df_non_up$decision)

# get FASTA sequences to replace the "inactive" sequences
payload <- list(
  query = paste(
    df_non_up[grep("replace", df_non_up$decision), "accession"],
    collapse = " "
  ),
  from = "ACC+ID",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  crap_uniparc_non_up_fasta <- Biostrings::readAAStringSet(tmp)
  message(paste("Input length:", length(df_non_up[grep("replace", df_non_up$decision), "accession"])))
  message(paste("Output FASTA length:", length(crap_uniparc_non_up_fasta)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

###-------------------------------###
#### Combine all FASTAs into one ####
###-------------------------------###

crap_combined <- c(
  crap_uniprot,
  crap_uniparc_up_fasta,
  crap_uniparc_non_up_fasta
)

# check length
length(crap_combined) # 231 sequences

# any duplicated sequences?
any(duplicated(crap_combined)) # yes

# remove duplicates
crap_combined_unique <- crap_combined[!duplicated(crap_combined)]

# check length
length(crap_combined_unique) # 227 sequences

###----------###
#### Output ####
###----------###

# sort alphabetically by entry names
sort_order <- regexec(
  "(?<=\\|)[A-Z,0-9]+_[A-Z]+|UPI[0-9,A-Z]{10}", 
  names(crap_combined_unique), 
  perl = TRUE
) %>%
  regmatches(names(crap_combined_unique), .) %>% 
  unlist() %>% 
  order()

output <- crap_combined_unique[sort_order]

# get the current UniProt release
cur_release <- response$headers$`x-uniprot-release`

# write updated MaxQuant contaminants FASTA file
Biostrings::writeXStringSet(
  output, 
  filepath = paste0("fasta/mq_contaminants_", cur_release, ".fasta")
)
