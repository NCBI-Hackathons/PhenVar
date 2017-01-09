#!/usr/bin/env python2

from Bio import Entrez
import time # We need this if we want to filter based on age
from urllib2 import HTTPError
import xml.etree.ElementTree as ET

# Configuration settings to be moved externally later
email = "jon.demasi@colorado.edu"


""" 
Take a search term, intended to be an rs#, 
and return a list of pmids
"""
def get_pmids(interm):
    if isinstance(interm, str):
        interm = interm + " AND pubmed_snp_cited[sb]"
        Entrez.email = email     # Always tell NCBI who you are
        search_results = Entrez.read(Entrez.esearch(db="pubmed", term=interm, usehistory="y"))
        return search_results
    elif isinstance(interm, list):
        searchstring = " OR ".join(interm)
        searchstring = "(" + searchstring + ") AND pubmed_snp_cited[sb]"
        Entrez.email = email     # Always tell NCBI who you are
        search_results = Entrez.read(Entrez.esearch(db="pubmed", term=searchstring, usehistory="y"))
        return search_results



""" 
Takes the saved list of results from get_pmids and retrieves
the raw XML to parse.  Currently goes article by article instead of 
doing them in bulk.  Could go either way. Unsure which way
is better / more efficient.
""" 
def get_abstracts(results):
    abstracts_list=[]
    Entrez.email = email 
    for start in range(0,int(results["Count"]),1):
        end = min(int(results["Count"]), start+1)
        #print("Going to download record %i to %i" % (start+1, end))
        fetch_handle = Entrez.efetch(db="pubmed",rettype="abstract",
                                         retmode="xml",retstart=start,
                                         retmax=1,
                                         webenv=results["WebEnv"],
                                         query_key=results["QueryKey"])
        data = fetch_handle.read()
        fetch_handle.close()
        root=ET.fromstring(data)
        for abst in root.iter('Abstract'):
            for sec in abst.iter('AbstractText'):
                abstracts_list.append(sec.text)
    return abstracts_list


# Just using this to test get_pmids function

toquery=get_pmids(["rs328", "rs360"])

# List for abstracts, ready for nltk
abstracts = get_abstracts(toquery)
print(abstracts)
