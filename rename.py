#!/usr/bin/env python2

from Bio import Entrez
import time # We need this if we want to filter based on age
from urllib2 import HTTPError

# Configuratoin settings to be moved externally later
email = "jon.demasi@colorado.edu"


""" 
Take a search term, intended to be an rs#, 
and return a list of pmids
"""
def get_pmids(interm):
    #print("This is doing stuff")
    Entrez.email = email     # Always tell NCBI who you are
    search_results = Entrez.read(Entrez.esearch(db="pubmed", term=interm, usehistory="y"))
    #print(record)
    #print(record["IdList"])
    return search_results


def get_abstracts(results):
    Entrez.email = email 
    #handle = Entrez.efetch(db="pubmed", id="26999119", rettype="abstract", retmode="text")
    #print(handle.readline().strip())
    #handle.close()
    #fetch_handle = Entrez.efetch(db="pubmed",rettype="abstract",retmode="text",webenv=results["WebEnv"],query_key=results["QueryKey"])
    #data=fetch_handle.read()
    #print(data)
    #fetch_handle.close()
    for start in range(0,int(results["Count"]),1):
        end = min(int(results["Count"]), start+1)
        print("Going to download record %i to %i" % (start+1, end))
        fetch_handle = Entrez.efetch(db="pubmed",rettype="abstract",
                                         retmode="text",retstart=start,
                                         retmax=1,
                                         webenv=results["WebEnv"],
                                         query_key=results["QueryKey"])
        data = fetch_handle.read()
        fetch_handle.close()
        print(data)

# Just using this to test get_pmids function

toquery=get_pmids("rs328")
print(toquery["Count"])
get_abstracts(toquery)

