from Bio import Entrez
import time  # We need this if we want to filter based on age
from urllib2 import HTTPError
import xml.etree.ElementTree as ET
from settings import configuration
email = configuration["email"]


""" 
Generate a list of rsids that are explicitly cited
in a pubmed paper somewhere :-) 
"""

def get_complete_rsids():
    Entrez.email = email
    rsidlist = []
    numresults = 0 
    retstart = 0
    search_string = "snp_pubmed_cited[sb]"
    search_results = Entrez.read(Entrez.esearch(db="snp", term=search_string,
                                                retmax=100000, retstart=retstart, usehistory="y"))
    print("Found a total of " +
        search_results["Count"] + " results using search string '" + search_string + "'")
    numresults = search_results["Count"]
    rsidlist = rsidlist + search_results["IdList"]
    additional_queries = int(int(numresults)/100000)
    while additional_queries != 0:
        retstart = retstart + 100000
        search_results = Entrez.read(Entrez.esearch(db="snp", term=search_string,
                                                    retmax=100000, retstart=retstart, usehistory="y"))
        rsidlist = rsidlist + search_results["IdList"]
        additional_queries = additional_queries - 1
    return rsidlist

        
"""
Take a search term, intended to be an rs#,
and return a dictionary that can be used to 
efetch later
"""


def get_pmids(interm):
    # This is obsolete now, essentially, but
    # allows a user to pass a single string
    # which can be nice.
    if isinstance(interm, str):
        interm = interm + " AND pubmed_snp_cited[sb]"
        Entrez.email = email     # Always tell NCBI who you are
        search_results = Entrez.read(Entrez.esearch(db="pubmed", term=interm,
                                                    retmax=100000,
                                                    usehistory="y"))
        print("Found a total of " +
              search_results["Count"] + " results using search string '" + interm + "'")
        return search_results

    elif isinstance(interm, list):
        searchstring = " OR ".join(interm)
        searchstring = "(" + searchstring + ") AND pubmed_snp_cited[sb]"
        Entrez.email = email     # Always tell NCBI who you are
        search_results = Entrez.read(Entrez.esearch(db="pubmed",
                                                    term=searchstring,
                                                    retmax=100000,
                                                    usehistory="y"))
        print("Found a total of " +
              search_results["Count"] + " results using search string'" + searchstring + "'")
        return search_results


"""
Takes the saved list of results from get_pmids and retrieves
the raw XML to parse.  Currently goes article by article instead of
doing them in bulk.  Could go either way. Unsure which way
is better / more efficient.
"""


def get_abstracts(results):
    abstracts_list = []
    Entrez.email = email
    for start in range(0, int(results["Count"]), 1):
        end = min(int(results["Count"]), start + 1)
        # print("Going to download record %i to %i" % (start+1, end))
        fetch_handle = Entrez.efetch(db="pubmed", rettype="abstract",
                                        retmode="xml", retstart=start,
                                        retmax=1,
                                        webenv=results["WebEnv"],
                                        query_key=results["QueryKey"])
        data = fetch_handle.read()
        fetch_handle.close()
        root = ET.fromstring(data)
        for abst in root.iter('Abstract'):
            for sec in abst.iter('AbstractText'):
                abstracts_list.append(sec.text)
    return abstracts_list


def get_abstracts_from_list(pmids_list):
    #    abstracts_list = []
    pmids_abstracts_dict = {}
    Entrez.email = email
    for each_pmid in pmids_list:
        fetch_handle = Entrez.efetch(db="pubmed", id=each_pmid, retmode='xml')
        data = fetch_handle.read()
        fetch_handle.close()
        root = ET.fromstring(data)
        for abst in root.iter('Abstract'):
            for sec in abst.iter('AbstractText'):
                # abstracts_list.append(sec.text)
                #                pmids_abstracts_dict[each_pmid] = sec.text
                #    print abstracts_list
                # return abstracts_list
                pmids_abstracts_dict[each_pmid] = sec.text
    #print pmids_abstracts_dict
    return pmids_abstracts_dict

list = get_complete_rsids()
spot = len(list) - 1
while spot != 0:
    dict = get_pmids("rs"+list[spot])
    #print(dict)
    for y in dict["IdList"]:
        ylist=[]
        ylist.append(y)
        newdict = get_abstracts_from_list(ylist)
        print("ID: " + y)
        print("Abstract: ")
        print(newdict[y])
    #newlisty = get_abstracts(dict)
    time.sleep(5)
    spot = spot - 1
