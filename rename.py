#!/usr/bin/env python2

from Bio import Entrez


""" 
Take a search term, intended to be an rs#, 
and return a list of pmids
"""
def get_pmids(interm):
    print("This is doing stuff")
    Entrez.email = "jon.demasi@colorado.edu"     # Always tell NCBI who you are
    handle = Entrez.esearch(db="pubmed", term=interm)
    record = Entrez.read(handle)
    #print(record)
    #print(record["IdList"])
    return(record["IdList"])


def get_abstract(pmid):
    return
# Just using this to test get_pmids function
list = get_pmids("rs328")
print(list)
