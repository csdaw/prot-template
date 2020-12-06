###-------------------------###
#### What this script does ####
###-------------------------###

## Author: Charlotte Dawson (csdaw@outlook.com)

## MaxQuant comes with a contaminants FASTA file (245 sequences) that is used 
## during peptide searching using the Andromeda search engine. This FASTA is 
## relatively old and the headers contain several incorrect and out-of-date 
## accessions. 

## I have written this script to generate a FASTA file with (almost) identical
## sequences to the original MaxQuant contaminants.fasta. However, these 
## sequences are annotated with headers that contain modern, working
## accessions.

## Original FASTA source: MaxQuant, Max Plank Institute of Biochemistry
## http://coxdocs.org/doku.php?id=maxquant:start_downloads.htm&s[]=contaminants

###----------------###
#### Dependencies ####
###----------------###

library(httr)
library(magrittr)
library(rentrez)

## This script also calls functions from the following packages:
# Biostrings
# curl
# 

## Define some urls for GET requests
BASE <- "https://www.uniprot.org"
TOOL_ENDPOINT <- "/uploadlists/" # REST endpoint
KB_ENDPOINT <- "/uniprot/" # UniProt website search endpoint
UPARC_ENDPOINT <- "/uniparc/" # UniParc website search endpoint

###----------------------------------------###
#### Download MaxQuant contaminants FASTA ####
###----------------------------------------###

# download original MaxQuant contaminants FASTA if it does not already exist
if (!file.exists("fasta/mq_contaminants_original.fasta")) {
  curl::curl_download(
    url = "http://lotus1.gwdg.de/mpg/mmbc/maxquant_input.nsf/7994124a4298328fc125748d0048fee2/$FILE/contaminants.fasta", 
    destfile = "fasta/mq_contaminants_original.fasta"
  )
}

# load MaxQuant contaminants FASTA
mq_crap <- Biostrings::readAAStringSet("fasta/mq_contaminants_original.fasta")

###-----------------------------------------###
#### Extract accessions from FASTA headers ####
###-----------------------------------------###

# extract all accessions (UniProt or non-UniProt) from FASTA headers
accessions <- regexec("^[^ ]+", names(mq_crap)) %>% 
  regmatches(names(mq_crap), .) %>% 
  unlist()

# extract UniProt and non-Uniprot accessions

# UniProt accessions (length: 210)
accessions_up <- grep(
  "^[QPOA][A-Z,0-9]{5}", 
  accessions, 
  value = TRUE
)

# Non-UniProt accessions (length: 35)
accessions_non_up <- grep(
  "^[QPOA][A-Z,0-9]{5}", 
  accessions, 
  invert = TRUE, 
  value = TRUE
)

###----------------------------###
#### UniProt accessions (210) ####
###----------------------------###

## Here we separate the canonical (non-isoform) UniProt accessions (209 entries)
## from the isoform UniProt accessions (1 entry).

# extract canonical or isoform UniProt accessions
up_can_accessions <- grep("-[2-9]", accessions_up, invert = TRUE, value = TRUE)
up_iso_accessions <- grep("-[2-9]", accessions_up, value = TRUE)

###----------------------------------###
#### UniProt isoform accessions (1) ####
###----------------------------------###

## First we will deal with the isoform UniProt accessions (1 entry)

# extract isoform UniProt accessions
up_iso_accessions <- grep("-[2-9]", accessions_up, value = TRUE)

payload <- list(
  query = "UPI00001D9675", # MaxQuant sequence doesn't match header
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, UPARC_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  up_iso_fasta <- Biostrings::readAAStringSet(tmp)
  print(paste("Input length:", length(up_iso_accessions)))
  print(paste("Output FASTA length:", length(up_iso_fasta)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

# check intersection of sequences, should be 1/1
length(Biostrings::intersect(up_iso_fasta, mq_crap)) == 1

###--------------------------------------###
#### UniProt canonical accessions (209) ####
###--------------------------------------###

## Then we will deal with the canonical accessions (209 entries)

payload <- list(
  query = paste(up_can_accessions, collapse = " "),
  from = "ACC+ID",
  to = "UPARC",
  format = "tab"
)

response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload)

if (response$status_code == 200) {
  up_can_mapping <- data.table::fread(content(response, encoding = "UTF-8")) %>% 
    `colnames<-`(c("old_accession", "uparc_accession"))
  message(paste("Input length:", length(up_can_accessions)))
  message(paste("Output length:", length(unique(up_can_mapping$uparc_accession))))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

## 16 MaxQuant FASTA accessions don't map 1:1 to UniParc accessions.
## We will deal with these "duplicated" accessions separately.

# get sequences for old accessions that don't map 1:1
up_can_dups <- up_can_mapping[duplicated(up_can_mapping$uparc_accession) | 
                                duplicated(up_can_mapping$uparc_accession, fromLast = TRUE), ]

idx <- sapply(
  up_can_dups$old_accession,
  function(x) grep(x, names(mq_crap))
)

up_can_dups_mq <- mq_crap[idx]

# manually fix an incorrect accession before obtaining sequences
up_can_dups[12, 1] <- "Q6E0U4"

# obtain sequences with UniParc accession headers using old accessions 
# (manually add extra UniParc accession)
payload <- list(
  query = paste(c(up_can_dups$old_accession, "UPI0000EFFEE3"), collapse = " OR "),
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, UPARC_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  up_can_dups_uparc <- Biostrings::readAAStringSet(tmp)
  message(paste("Input length:", length(up_can_dups$old_accession)))
  message(paste("Output FASTA length:", length(up_can_dups_uparc)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

# compare sequences from new FASTA to old MaxQuant FASTA (should be 16/16)
up_can_dups_uparc <- Biostrings::intersect(up_can_dups_uparc, up_can_dups_mq)

length(up_can_dups_uparc) == 16

# use UniParc accessions to get some modern UniProtKB headers if available
payload <- list(
  query = paste(gsub("(?<=UPI[0-9,A-Z]{10}).*", "", names(up_can_dups_uparc), perl = TRUE), collapse = " "),
  from = "UPARC",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  up_can_dups_new <- Biostrings::readAAStringSet(tmp)
  message(paste("Input length:", length(up_can_dups_uparc)))
  message(paste("Output FASTA length:", length(up_can_dups_new)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

# remove a duplicate sequence in the "duplicated" sequences
up_can_dups_new <- up_can_dups_new[!Biostrings::duplicated(up_can_dups_new)]

# combine the 16 "duplicate" sequences into a single FASTA
up_can_dups_fasta <- c(
  up_can_dups_new, 
  Biostrings::setdiff(up_can_dups_uparc, up_can_dups_new)
)

# check intersection of sequences, should be 16/16
length(Biostrings::intersect(up_can_dups_fasta, mq_crap)) == 16

## 193 MaxQuant FASTA accessions do map 1:1 to UniParc accessions.
## We will deal with these "unique" accessions separately.

# get sequences for old accessions that do map 1:1
up_can_unique <- up_can_mapping[!up_can_mapping$uparc_accession %in% 
                                  up_can_dups$uparc_accession, ]

# obtain sequences with UniParc accession headers using old accessions 
payload <- list(
  query = paste(up_can_unique$old_accession, collapse = " OR "),
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, UPARC_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  up_can_unique_uparc <- Biostrings::readAAStringSet(tmp)
  message(paste("Input length:", length(up_can_unique$old_accession)))
  message(paste("Output FASTA length:", length(up_can_unique_uparc)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

idx <- sapply(
  up_can_unique$old_accession,
  function(x) grep(x, names(mq_crap))
)

up_can_unique_mq <- mq_crap[idx]

# compare sequences from new FASTA to old MaxQuant FASTA (should be 193/193)
length(Biostrings::intersect(up_can_unique_uparc, up_can_unique_mq)) == 193

## 162/193 are good because MaxQuant headers are correct and do match the protein
## sequences in the MaxQuant FASTA.
up_can_unique_good <- Biostrings::intersect(up_can_unique_uparc, up_can_unique_mq)
length(up_can_unique_good)

## 31/193 are bad because MaxQuant headers are wrong and don't match the protein
## sequences in the MaxQuant FASTA.
## I manually found new headers that match sequences by BLAST searching the
## sequences against the UniParc database.
up_can_unique_bad <- Biostrings::setdiff(up_can_unique_mq, up_can_unique_uparc)
length(up_can_unique_bad)

manual_accessions <- c(
  "Q86Y46", # UPI000000DCB8
  "UPI0000D9FD95",
  "Q9R4J4", # UPI0001D8971A MaxQuant sequence not exact match, is fragment
  "UPI0000F01707",
  "UPI000061528B",
  "UPI0000614433",
  "UPI0000EBCB5C",
  "A0A4W2IJA9", # UPI0000615295,
  "A5PKC2", # UPI0000EBDC90
  "F1N3A1", # UPI00006160D9
  "UPI0000EBDFDE",
  "Q68RU0", # UPI0000418D31
  "A0A4W2II28", # UPI0000EBD983
  "E1BH06", # UPI0000EBE2E7,
  "UPI000061462D",
  "UPI000061587D",
  "UPI0000EBD5A3",
  "A5D7R6", # UPI0000EBD43D
  "UPI000061326A",
  "Q08E14", # UPI000003BB1F
  "UPI00001A4820",
  "UPI000245108C",
  "UPI0000174051",
  "UPI0001D7C7DE",
  "P31096", # UPI000016C265 MaxQuant sequence not exact match, is fragment
  "UPI0000111542", 
  "E9Q1Y9", # UPI0000429C2F
  "Q9EQD6", # UPI00000E6123
  "Q6IFX4", # UPI00001C4C3F
  "Q0IIN1", # UPI0000DD7597
  "UPI00001C1108"
)

# obtain sequences with UniParc accession headers using old accessions 
payload <- list(
  query = paste(manual_accessions, collapse = " OR "),
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, UPARC_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  up_can_unique_bad_fixed <- Biostrings::readAAStringSet(tmp)
  message(paste("Input length:", length(manual_accessions)))
  message(paste("Output FASTA length:", length(up_can_unique_bad_fixed)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

# combine all the UniParc accessions obtained for the various sequences into
# one FASTA
up_can_unique_uparc <- c(
  up_can_unique_good,
  Biostrings::intersect(up_can_unique_bad_fixed, mq_crap),
  up_can_unique_bad_fixed[c(33, 7)]
)

# check that the length is still 193 sequences
length(up_can_unique_uparc) == 193

# check intersection of sequences with MaxQuant fasta, should be 191/193
length(Biostrings::intersect(up_can_unique_uparc, mq_crap)) == 191

# use UniParc accessions to get some modern UniProtKB headers if available
payload <- list(
  query = paste(gsub("(?<=UPI[0-9,A-Z]{10}).*", "", names(up_can_unique_uparc), perl = TRUE), collapse = " "),
  from = "UPARC",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  up_can_unique_new <- Biostrings::readAAStringSet(tmp)
  message(paste("Input length:", length(up_can_unique_uparc)))
  message(paste("Output FASTA length:", length(up_can_unique_new)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

# combine 193 "unique" into single FASTA
up_can_unique_fasta <- c(
  Biostrings::intersect(up_can_unique_new, mq_crap),
  Biostrings::setdiff(up_can_unique_uparc, up_can_unique_new),
  Biostrings::setdiff(up_can_unique_new, mq_crap)
)

# check that the length is still 193 sequences
length(up_can_unique_fasta) == 193

# check intersection of sequences, should be 191/193
length(Biostrings::intersect(up_can_unique_fasta, mq_crap)) == 191

###---------------------###
#### RefSeq accessions ####
###---------------------###

# extract RefSeq accessions
refseq_accessions <- regexec("(?<=REFSEQ:)[A-Z,0-9,_]+", accessions_non_up, perl = TRUE) %>% 
  regmatches(accessions_non_up, .) %>% 
  unlist()

# get Entrez summary for each RefSeq accession
refseq_summary <- entrez_summary(db = "protein", id = refseq_accessions)
names(refseq_summary) <- refseq_accessions

# see the status of each record
refseq_statuses <- sapply(refseq_summary, "[", "status") %>% unlist()

## We can see that the records for 4 out of 6 refseq accessions have been 
## "removed as a result of standard genome annotation processing".
## We will deal with these separately.

# get suppressed refseq accessions
refseq_bad <- names(refseq_statuses)[grep("suppressed", refseq_statuses)] %>% 
  gsub(".status", "", .)

# remove them from the summary for now
refseq_summary <- refseq_summary[-grep("suppressed", refseq_statuses)]

# obtain sequences with UniParc accession headers using bad refseq accessions 
payload <- list(
  query = paste(refseq_bad, collapse = " OR "),
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, UPARC_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  refseq_bad_uparc <- Biostrings::readAAStringSet(tmp)
  print(paste("Input length:", length(refseq_summary)))
  print(paste("Output FASTA length:", length(refseq_bad_fasta)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

idx <- sapply(
  refseq_bad,
  function(x) grep(x, names(mq_crap))
)

# check intersection of sequences, should be 4/4
refseq_bad_uparc <- Biostrings::intersect(refseq_bad_uparc, mq_crap[idx])
length(refseq_bad_uparc) == 4

## Also we can see that sequences referred to by the remaining 2 accessions 
## have been replaced. Therefore we will use the modern protein sequences in
## our final updated FASTA.

refseq_good <- names(refseq_summary)

payload <- list(
  query = paste(refseq_good, collapse = " OR "),
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, UPARC_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  refseq_good_uparc <- Biostrings::readAAStringSet(tmp)
  print(paste("Input length:", length(refseq_summary)))
  print(paste("Output FASTA length:", length(refseq_good_uparc)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

idx <- sapply(
  refseq_good,
  function(x) grep(x, names(mq_crap))
)

# check intersection of sequences, should be 2/2
refseq_good_uparc <- Biostrings::intersect(refseq_good_uparc, mq_crap[idx])
length(refseq_good_uparc) == 2

## Now we combine the RefSeq sequences with their new UniParc headers and see
## if we can map UniProtKB headers to any of them.
refseq_uparc <- c(refseq_good_uparc, refseq_bad_uparc)

# use UniParc accessions to get some modern UniProtKB headers if available
payload <- list(
  query = paste(gsub("(?<=UPI[0-9,A-Z]{10}).*", "", names(refseq_uparc), perl = TRUE), collapse = " "),
  from = "UPARC",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  refseq_new <- Biostrings::readAAStringSet(tmp)
  message(paste("Input length:", length(refseq_uparc)))
  message(paste("Output FASTA length:", length(refseq_new)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

# remove a duplicate sequence in refseq_new
refseq_new <- refseq_new[!Biostrings::duplicated(refseq_new)]

refseq_fasta <- c(
  Biostrings::setdiff(refseq_uparc, refseq_new),
  refseq_new
)

# check that the length is still 6 sequences
length(refseq_fasta) == 6

# check intersection of sequences, should be 6/6
length(Biostrings::intersect(refseq_fasta, mq_crap)) == 6

###----------------------###
#### H-InvDb accessions ####
###----------------------###

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

###----------------------###
#### Ensembl accessions ####
###----------------------###

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

###---------###
#### Other ####
###---------###

## Lastly, we will deal with any other proteins in the original FASTA
# extract anything else non-UniProt
man_accessions <- grep(
  "(ENSEMBL|REFSEQ|H-INV)", 
  accessions_other, 
  invert = TRUE, value = TRUE, perl = TRUE
)

# it is streptavidin which we will just retrieve manually UPI000002B867
payload <- list(
  query = "P22629",
  from = "ACC+ID",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  oth_fasta <- Biostrings::readAAStringSet(tmp)
  print(paste("Input length:", length("P22629")))
  print(paste("Output FASTA length:", length(oth_fasta)))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

###----------###
#### Output ####
###----------###

## Combine all the FASTAs
output <- c(
  up_can_unique_fasta,
  up_can_dups_fasta,
  up_iso_fasta,
  refseq_fasta,
  ens_fasta,
  oth_fasta
)
