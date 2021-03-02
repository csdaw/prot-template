## Dependencies
library(httr)
library(magrittr)

# define some urls for GET requests
BASE <- "https://www.uniprot.org"
TOOL_ENDPOINT <- "/uploadlists/" # REST endpoint

## This script also calls functions from the following packages:
# Biostrings

## These are additional RNase contaminants present in OOPS experiments
crap_oops <- c(
  "P13717", # Nuclease, GN=nucA
  "P61823", # Ribonuclease pancreatic, GN = RNASE1
  "P00651" # Guanyl-specific ribonuclease T1, GN=rntA
)

# get FASTA sequences
payload <- list(
  query = paste(crap_oops, collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  output <- Biostrings::readAAStringSet(tmp)
  message(paste("Input length:", length(crap_oops)))
  message(paste("Output FASTA length:", length(output)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

# get the current UniProt release
cur_release <- response$headers$`x-uniprot-release`

# load csd_cRAP FASTA
crap_csd <- Biostrings::readAAStringSet(paste0("fasta/csd_cRAP_", cur_release, ".fasta"))

# combine with OOPS contaminants
combined <- c(output, crap_csd)

# sort alphabetically
sort_order <- regexec(
  "(?<=\\|)[A-Z,0-9]+_[A-Z]+|UPI[0-9,A-Z]{10}", 
  names(combined), 
  perl = TRUE
) %>%
  regmatches(names(combined), .) %>% 
  unlist() %>% 
  order()

output <- combined[c(sort_order, length(combined)-1, length(combined))]

# write updated csd_cRAP FASTA file with OOPS contaminants
Biostrings::writeXStringSet(
  output, 
  filepath = paste0("fasta/OOPS_cRAP_", cur_release, ".fasta")
)
