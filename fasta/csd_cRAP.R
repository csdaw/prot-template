## Dependencies
library(httr)
library(magrittr)

# define some urls for GET requests
BASE <- "https://www.uniprot.org"
TOOL_ENDPOINT <- "/uploadlists/" # REST endpoint

## This script also calls functions from the following packages:
# Biostrings

# First, add some extra proteins from an internal FASTA
extra_accessions <- c(
  "P00921", # Carbonic anhydrase 2, GN=CA2
  "P01012", # Ovalbumin, GN=SERPINB14
  "P81054", # Peptidyl-Lys metalloendopeptidase, GN=MEP
  "P0AEX9", # Maltose/maltodextrin-binding periplasmic protein, GN=malE
  "P00772" # Chymotrypsin-like elastase family member 1, GN=CELA1
)

# get their FASTA sequences
payload <- list(
  query = paste(extra_accessions, collapse = " "),
  from = "ACC+ID",
  to = "ACC",
  format = "fasta"
)

tmp <- tempfile()
response <- GET(url = paste0(BASE, TOOL_ENDPOINT), query = payload, config = write_disk(tmp))

if (response$status_code == 200) {
  extra_seqs <- Biostrings::readAAStringSet(tmp)
  message(paste("Input length:", length(extra_accessions)))
  message(paste("Output FASTA length:", length(extra_seqs)))
} else {
  stop("Something went wrong. Status code: ", response$status_code)
}

# get the current UniProt release
cur_release <- response$headers$`x-uniprot-release`

# combine GPM and MaxQuant FASTAs from the current release
gpm <- Biostrings::readAAStringSet(paste0("fasta/gpm_cRAP_", cur_release, ".fasta"))
mq <- Biostrings::readAAStringSet(paste0("fasta/mq_contaminants_", cur_release, ".fasta"))

combined <- Biostrings::union(gpm, mq)

# add sequence for NEB Endoproteinase GluC (P8100S)
gluc <- Biostrings::AAStringSet("AGYRDGFGASGSCEVDAVCATQSGTRAYDNATAAVAKMVFTSSADGGSYICTGTLLNNGNSPKRQLFWSAAHCIEDQATAATLQTIWFYNTTQCYGDASTINQSVTVLTGGANILHRDAKRDTLLLELKRTPPAGVFYQGWSATPIANGSLGHDIHHPRGDAKKYSQGNVSAVGVTYDGHTALTRVDWPSAVVEGGSSGSGLLTVAGDGSYQLRGGLYGGPSYCGAPTSQRNDYFSDFSGVYSQISRYFAPHQHQHQHQHQ")
names(gluc) <- "P8100S|NEB Endoproteinase GluC"

# add sequence for Promega recombinant Lys-C (V1671)
lysc <- Biostrings::AAStringSet("VILPNNDRHQITDTTNGHYAPVTYIQVEAPTGTFIASGVVVGKDTLLTNKHVVDATHGDPHALKAFPSAINQDNYPNGGFTAEQITKYSGEGDLAIVKFSPNEQNKHIGEVVKPATMSNNAETQVNQNITVTGYPGDKPVATMWESKGKITYLKGEAMQYDLSTTGGNSGSPVFNEKNEVIGIHWGGVPNEFNGAVFINENVRNFLKQNIEDIHFANDDQPNNPDNPDNPNNPDNPNNPDEPNNPDNPNNPDNPDNGDNNNSDNPDAAHHHHHH")
names(lysc) <- "V1671|Promega recombinant Lys-C"

# add the extra UniProt and commercial protease sequences
combined <- c(combined, extra_seqs, gluc, lysc)

## Output
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

# write FASTA
Biostrings::writeXStringSet(
  output, 
  filepath = paste0("fasta/csd_cRAP_", cur_release, ".fasta")
)
