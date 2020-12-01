library(curl)
library(httr)
library(magrittr)
## also uses functions from:
# data.table
# Biostrings

# Ron Beavis group, the Global Proteome Machine (GPM) (link not always working)
# See https://www.thegpm.org/crap/

# Note that this FASTA is relatively old and some of the sequences do not 
# precisely match those in the current UniProt database e.g. ADH1 S. cerevisiae
# so we will not use this FASTA directly.
if (!file.exists("fasta/gpm_crap_original.fasta")) {
  curl::curl_download(
    url = "ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta", 
    destfile = "fasta/gpm_crap_original.fasta"
  )
}

# load GPM cRAP fasta
crap_gpm <- Biostrings::readAAStringSet("fasta/gpm_crap_original.fasta")

# Extract protein entry names (aka mnemonics)
# note: it is not good to use entry names as they are not stable identifiers,
# rather UniProt accessions should always be used. Therefore we will use the
# entry names to obtain accessions.
entry_names <- regexec("[A-Z,0-9]+_[A-Z]+", names(crap_gpm)) %>% 
  regmatches(names(crap_gpm), .) %>% 
  unlist()

# define the urls we need
BASE <- "https://www.uniprot.org"
TOOL_ENDPOINT <- "/uploadlists/" # REST endpoint
KB_ENDPOINT <- "/uniprot/" # UniProt website search endpoint

# check which mnemonics are still valid and which are old
payload <- list(
  query = paste(entry_names, collapse = " "),
  from = "ACC+ID",
  to = "ID",
  format = "list"
)

response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload)

valid_mnemonics <- strsplit(content(response), split = "\n") %>% 
  unlist()

old_mnemonics <- entry_names[!entry_names %in% valid_mnemonics]

# obtain FASTA with UniProt accessions for valid mnemonics
payload <- list(
  query = paste(valid_mnemonics, collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  valid_mnemonics_fasta <- Biostrings::readAAStringSet(tmp)
  print(length(valid_mnemonics_fasta))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

# obtain FASTA with UniProt accession for old mnemonics
payload <- list(
  query = paste0("mnemonic:", paste(old_mnemonics, collapse = " OR ")),
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, KB_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  old_mnemonics_fasta <- Biostrings::readAAStringSet(tmp)
  print(length(old_mnemonics_fasta))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

# combine valid and old mnemonic FASTAs
crap_gpm_updated <- Biostrings::union(valid_mnemonics_fasta, old_mnemonics_fasta)

# sort according to entry name (alphabetical)
sort_order <- regexec("(?<=sp\\|[A-Z,0-9]{6}\\|)[A-Z,0-9]+_[A-Z]+", names(crap_gpm_updated), perl = TRUE) %>%
  regmatches(names(crap_gpm_updated), .) %>% 
  unlist() %>% 
  order()

crap_gpm_updated <- crap_gpm_updated[sort_order]

# get the current UniProt release
cur_release <- response$headers$`x-uniprot-release`

# write updated GPM cRAP FASTA file
Biostrings::writeXStringSet(crap_gpm_updated, filepath = paste0("fasta/gpm_cRAP_", cur_release, ".fasta"))
