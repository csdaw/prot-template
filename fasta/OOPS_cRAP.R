## Dependencies
library(httr)

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

# write updated MaxQuant contaminants FASTA file
Biostrings::writeXStringSet(
  output, 
  filepath = paste0("fasta/OOPS_cRAP_", cur_release, ".fasta")
)
