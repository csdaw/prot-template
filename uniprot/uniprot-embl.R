library(httr)

# simple example which mirrors the python one
BASE <- "https://www.uniprot.org"
KB_ENDPOINT <- "/uniprot/" # REST endpoint

payload <- list(
  query = 'name:"polymerase alpha" AND taxonomy:mus AND reviewed:yes',
  format = "list"
)

response <- GET(paste0(BASE, KB_ENDPOINT), query = payload)

if (response$status_code < 400) {
  print(strsplit(content(response, "text"), split = "\n")[[1]])
} else {
  print("Something went wrong: ", response$status_code)
}

response$headers

# a real life example
# try a larger query payload (a real life example, Homo sapiens)
# change the output format
payload <- list(
  query = "proteome:UP000005640 AND reviewed:yes",
  format = "tab",
  columns = paste(
    c(
      "entry name", "go", "genes(ORF)", "genes(PREFERRED)", "genes", "protein names",
      "length", "feature(DOMAIN EXTENT)", "comment(DOMAIN)",
      "feature(BINDING SITE)", "features", "sequence", "reviewed"
    ),
    collapse = ","
  )
)

response <- GET(paste0(BASE, KB_ENDPOINT), query = payload)

response$headers

response$status_code == 200 # check if the request is 'ok'

# r_content should be character vector length 1. It is always reencoded to UTF-8.
# If the re-encoding fails then r_content will just be NA.
r_content <- content(response)
xxx <- readr::read_tsv(r_content)
xxx <- data.table::fread(r_content)
