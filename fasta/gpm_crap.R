library(curl)
library(httr)
library(magrittr)
library(rvest)
## also uses functions from:
# data.table
# Biostrings

BASE <- "https://www.uniprot.org"
TOOL_ENDPOINT <- "/uploadlists/" # REST endpoint

## Ron Beavis group, the Global Proteome Machine (GPM) (link not always working)
## See https://www.thegpm.org/crap/

## This FASTA is relatively old and not well annotated. Therefore, we want to 
## get the proteins from this original GPM cRAP FASTA and create a new FASTA
## with the protein sequences from the most recent UniProt release.

# download original GPM cRAP FASTA if it does not already exist
if (!file.exists("fasta/gpm_cRAP_original.fasta")) {
  curl::curl_download(
    url = "ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta", 
    destfile = "fasta/gpm_cRAP_original.fasta"
  )
}

# load GPM cRAP fasta
gpm_crap <- Biostrings::readAAStringSet("fasta/gpm_crap_original.fasta")

## First, we extract the entry names (aka mnemonics) from the FASTA headers.
## Note: it is not good to use entry names as they are not stable identifiers.
## Rather it is important always to use UniProt accessions AND sequence versions
## in FASTA headers.
entry_names <- regexec("[A-Z,0-9]+_[A-Z]+", names(gpm_crap)) %>% 
  regmatches(names(gpm_crap), .) %>% 
  unlist()

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

# overwrite entry_name of extra protein with accession
extra_protein <- strsplit(content(response), split = "\n") %>% 
  unlist()

## Now we have the UniProt accessions for the GPM cRAP FASTA, we can obtain a
## new FASTA with the modern UniProt headers for each protein. These headers
## include the most important information to retain with the protein sequences,
## specifically their UniProt accessions and sequence version numbers.
payload <- list(
  query = paste(c(gpm_crap_table$accession, extra_protein), collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "tab",
  columns = paste(c("entry_name", "id"), collapse = ",")
)

response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload)

# extract table of old accessions mapped to modern accessions and entry names
gpm_mapping <- data.table::fread(content(response)) %>% 
  `colnames<-`(c("new_entry_name", "new_accession", "old_accession"))

# add old entry names too
gpm_mapping$old_entry_name <- camprotR::match_id_(
  gpm_mapping$old_accession, 
  gpm_crap_table, 
  "accession", 
  "id"
  ) %>% 
  unlist()

# fill in the NA for the E. coli protein, in the old_entry_name column
source("R/fill_na.R")

gpm_mapping <- fill_na(gpm_mapping, "old_entry_name", "new_entry_name")

## The 116 proteins from gpm_crap have mapped to 119 modern UniProt entries.
## Here we manually remove the undesired duplicates.

# find duplicates
dup <- gpm_mapping[duplicated(gpm_mapping$old_accession) | 
                     duplicated(gpm_mapping$old_accession, fromLast = TRUE), ]

# remove duplicates
gpm_mapping <- gpm_mapping[!gpm_mapping$new_entry_name %in% 
                             c("SSPA_STAA8", "AMY1B_HUMAN", "AMY1C_HUMAN")]

## Finally we obtain a FASTA with modern, full UniProtKB headers
payload <- list(
  query = paste(gpm_mapping$new_accession, collapse = " "),
  from = "ACC",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(
  url = paste0(BASE, TOOL_ENDPOINT), 
  query = payload, 
  config = write_disk(tmp)
)

if (response$status_code == 200) {
  gpm_crap_new <- Biostrings::readAAStringSet(tmp)
  print(length(gpm_crap_new))
} else {
  print("Something went wrong. Status code: ", response$status_code)
}

## Note that the sequences of the original GPM fasta and our new GPM fasta
## do not completely overlap. This is due to updated in protein sequences 
## between now and then (e.g. PEPA_PIG in the original FASTA was v2 of the
## sequence whereas the sequence we use in the updated FASTA is v3). 

## This is why it is important that FASTA files retain the full UniProt headers
## with the accession and sequence version number (SV) intact! Because the 
## original GPM FASTA did not include the sequence versions we can't easily
## replicate it.
length(Biostrings::intersect(gpm_crap_new, gpm_crap))

## Finally, we write our updated FASTA to a file with the current UniProt
## release name appended.

# sort according to entry name (alphabetical)
sort_order <- regexec(
  "(?<=sp\\|[A-Z,0-9]{6}\\|)[A-Z,0-9]+_[A-Z]+", 
  names(gpm_crap_new), 
  perl = TRUE
) %>%
  regmatches(names(gpm_crap_new), .) %>% 
  unlist() %>% 
  order()

gpm_crap_new <- gpm_crap_new[sort_order]

# get the current UniProt release
cur_release <- response$headers$`x-uniprot-release`

# write updated GPM cRAP FASTA file
Biostrings::writeXStringSet(
  gpm_crap_new, 
  filepath = paste0("fasta/gpm_cRAP_", cur_release, ".fasta")
)
