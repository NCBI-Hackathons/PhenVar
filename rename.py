#!/usr/bin/env python2

from Bio import Entrez

def get_from_search_term(interm):
    print("This is doing stuff")
    Entrez.email = "jon.demasi@colorado.edu"     # Always tell NCBI who you are
    handle = Entrez.esearch(db="pubmed", term=interm)
    record = Entrez.read(handle)
    print(record)
    print(record["IdList"])
    return

get_from_search_term("rs328")
