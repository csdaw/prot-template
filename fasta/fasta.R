library(curl)
library(magrittr)
## also uses functions from:
# data.table
# Biostrings

# UniProt reference proteome url
BASE <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes"

# Find out what proteomes are available by parsing the README
readme_url <- paste(BASE, "README", sep = "/")

# Open a connection to the ftp server and grab the README
connection <- curl(url = readme_url, open = "r")
readme <- readLines(connection)
close(connection)

# Which release are we downloading from?
release <- readme[base::which(grepl("Release", readme))][1] %>% 
  sub("Release ", "", .) %>% 
  strsplit(", ") %>% 
  unlist()

# How many species are available in the FTP reference_proteomes folder
total_species <- readme[which(grepl("Total of species", readme))][1] %>% 
  sub(".* ([0-9]+).*", "\\1", .) %>% 
  as.integer()

# Get the start and end line numbers for the table of interest
idx1 <- which(grepl("Proteome_ID", readme))
idx2 <- idx1 + total_species

# Extract the table of interest
readme_subset <- paste(readme[idx1:idx2], collapse = "\n")

proteomes_tbl <- data.table::fread(
  readme_subset,
  col.names = c(
    "proteome_id", 
    "tax_id", 
    "oscode", 
    "superregnum", 
    "n_canonical", 
    "n_isoform", 
    "n_gene2acc", 
    "species_name"
  )
)

# superregnum: one of Eukaryota, Bacteria, Archaea, Viruses
SUPERREGNUM <- "Eukaryota" # case sensitive!

# reference proteome of interest (Homo sapiens, canonical FASTA)
REF_PROTEOME <- "UP000005640_9606.fasta.gz"

url <- paste(BASE, SUPERREGNUM, REF_PROTEOME, sep = "/")
destfile <- paste("fasta", REF_PROTEOME, sep = "/")
curl::curl_download(url = url, destfile = destfile)

xxx <- Biostrings::fasta.index("fasta/UP000005640_9606.fasta.gz")

source("R/regex.R")

hhh <- lapply(
  xxx$desc,
  parse_uniprot_header
) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

colnames(hhh) <- c(
  "uniprot_db",
  "accession",
  "entry_name",
  "protein_name",
  "organism_name",
  "organism_id",
  "gene_name",
  "protein_existence",
  "sequence_version"
)

xxx2 <- do.call(cbind, c(xxx, hhh)) %>% 
  as.data.frame()
