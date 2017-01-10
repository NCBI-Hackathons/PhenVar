# PhenVar
PhenVar is designed to take one or more rsids and generate a list of PubMed IDs to query and generate novel associations between publications. It utilizes a local sqlite cache in a configurable location to keep a local copy of relevant srids, pmids, and the abstrat blobs for the pmids.  

### settings.py
settings.py is a file containing a single dictionary with configuration options that is imported into all of the other tools.  Configuration options you need can easily be added or changed.
### ncbiutils.py
#### get_pmids 
Takes a single rsid input (string) or several rsid inputs (list of strings) and returns a dictionary from Entrez.read/esearch.  The results returned are all pmids which explicitly cited the rsids given.  
#### get_abstracts
Expects the results from get_pmids as an input.  Will return a list of abstracts.  Each item in the list is an abstract to a pmid from the get_pmids search results.  Can return a list with only one item.  
