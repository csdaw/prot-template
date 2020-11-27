import requests

requests.__version__

BASE = "https://www.uniprot.org"
KB_ENDPOINT = "/uniprot/" # REST endpoint

## There are several different outputs potentially available 
## (not necessarily for all REST endpoints):
# list: list of accession numbers
# txt: full uniprot entries in plain text format (i.e. a long string)
# xml: full uniprot entries in xml format
# tab: selected fields in tabular format (must specify columns of interest)
# fasta: sequences in FASTA format (one sequence per entry (canonical) unless
#        include=yes which will include isoforms)
# gff: annotated sequence features in GFF format
# rdf: full uniprot entries in rdf format

#### LIST OUTPUT ####

# searching for proteins by protein name using a query string
payload = {
  "query": 'name:"polymerase alpha" AND taxonomy:mus AND reviewed:yes',
  "format": "list"
}

result = requests.get(BASE + KB_ENDPOINT, params = payload)

if result.ok:
    print(result.text)
else:
    print("Something went wrong: ", result.status_code)

# result is an object encapsulating the response that the UniProt server
# sent back when we asked for data

# check the header for more info about the reponse
for key, value in result.headers.items():
    print("{}: {}".format(key, value))

#### TAB OUTPUT ####
# see here for columns: https://www.uniprot.org/help/uniprotkb%5Fcolumn%5Fnames

# change the output format
payload = {
  "query": "proteome:UP000000559",
  "format": "tab",
  "columns": "entry name,genes(ORF)"
}

result = requests.get(BASE + KB_ENDPOINT, params = payload)
for key, value in result.headers.items():
    print("{}: {}".format(key, value))
    

if result.ok:
    print(result.text[:50])
else:
    print("Something went wrong: ", result.status_code)

# Need to get the tab output into a useable format
# as it is currently just one giant long string.
import pandas as pd
import numpy as np

# Several ways:
a = np.loadtxt(result.text, delimiter="\n", dtype=str)
# 1. Convert tab separated string to list of lists, then into array or df
result_list = [i.split("\t") for i in result.text.rstrip().split("\n")]

a = np.array(result_list) # convert to numpy array
a.shape

df = pd.DataFrame(result_list) # convert to pandas dataframe
df.shape

# 2. Use io module and numpy
from io import StringIO

a = np.loadtxt(
  StringIO(result.text),
  dtype = str,
  delimiter = "\t"
)

# 3. Use io module and pandas

df = pd.read_csv(
  StringIO(result.text),
  sep = "\t"
)


# try a larger query payload (a real life example, Homo sapiens)
# change the output format
payload = {
  "query": "proteome:UP000005640 AND reviewed:yes",
  "format": "tab",
  "include": "no",
  "columns": ",".join(
    [
      "entry name", "go", "genes(ORF)", "genes(PREFERRED)", "genes", "protein names",
      "length", "feature(DOMAIN EXTENT)", "comment(DOMAIN)", 
      "feature(BINDING SITE)", "features", "sequence", "reviewed"
    ]
  )
}

result = requests.get(BASE + KB_ENDPOINT, params = payload)

if result.ok:
    df = pd.read_csv(
  StringIO(result.text),
  sep = "\t"
)
    print(df.shape)
else:
    print("Something went wrong: ", result.status_code)

payload = {
  "query": "proteome:UP000005640",
  "format": "list",
  "include": "yes"
}

result = requests.get(BASE + KB_ENDPOINT, params = payload)

if result.ok:
    len(result.text.split("\n"))
else:
    print("Something went wrong: ", result.status_code)
    
    
#### Use urllib instead of requests ####
import urllib, urllib3

BASE = "https://www.uniprot.org"
KB_ENDPOINT = "/uniprot/" # REST endpoint

## There are several different outputs potentially available 
## (not necessarily for all REST endpoints):
# list: list of accession numbers
# txt: full uniprot entries in plain text format (i.e. a long string)
# xml: full uniprot entries in xml format
# tab: selected fields in tabular format (must specify columns of interest)
# fasta: sequences in FASTA format (one sequence per entry (canonical) unless
#        include=yes which will include isoforms)
# gff: annotated sequence features in GFF format
# rdf: full uniprot entries in rdf format

#### LIST OUTPUT ####

# searching for proteins by protein name using a query string
params = {
  "query": 'name:"polymerase alpha" AND taxonomy:mus AND reviewed:yes',
  "format": "list"
}
data = urllib.parse.urlencode(params).encode("utf-8")
request = urllib.request.Request(BASE + KB_ENDPOINT)
contact = "csdaw@outlook.com"
request.add_header("User Agent", "Python %s" % contact)

# this step seems to be relatively error prone
# probably why people prefer to use the requests module
response = urllib.request.urlopen(request, data = data)

print(response.headers)

# get the encoding of the data and decode the response
content = response.read().decode(response.headers.get_content_charset())
