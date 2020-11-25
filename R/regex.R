parse_uniprot_header <- function(x) {
  # UniProtKB database (swissprot or trembl), Accession, Entry
  start <- unlist(strsplit(x, "\\|"))
  
  # Extract values from remainder
  remainder <- regexec(
    "^([^\\s]+) (.*) OS=(.*) OX=([0-9]{2,7}) GN=(.*) PE=([1-5]) SV=([0-9])", 
    start[3], 
    perl = TRUE
  ) %>% 
    regmatches(start[3], .) %>% 
    unlist()
  
  # Create named vector as output
  result <- c(start[1:2], remainder[-1])
  result
}
