from Bio import Entrez
import time  # We need this if we want to filter based on age
from urllib2 import HTTPError
import xml.etree.ElementTree as ET
from settings import configuration
email = configuration["email"]


"""
Take a search term, intended to be an rs#,
and return a list of pmids
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
        #print("Found a total of " +
              #search_results["Count"] + " results using search string '" + interm + "'")
        return search_results

    elif isinstance(interm, list):
        searchstring = " OR ".join(interm)
        searchstring = "(" + searchstring + ") AND pubmed_snp_cited[sb]"
        Entrez.email = email     # Always tell NCBI who you are
        search_results = Entrez.read(Entrez.esearch(db="pubmed",
                                                    term=searchstring,
                                                    retmax=100000,
                                                    usehistory="y"))
        #print("Found a total of " +
              #search_results["Count"] + " results using search string'" + searchstring + "'")
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
        end = min(int(results["Count"]), start+1)
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
    """ return a dict where key is pmid and value is abstract """
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
                #abstracts_list.append(sec.text)
#                pmids_abstracts_dict[each_pmid] = sec.text
#    print abstracts_list
    #return abstracts_list
                pmids_abstracts_dict[each_pmid] = sec.text
#    print pmids_abstracts_dict
    return pmids_abstracts_dict
