import nltk
from Bio import Entrez
import time # We need this if we want to filter based on age
from urllib2 import HTTPError
import xml.etree.ElementTree as ET
from collections import Counter
import numpy as np
import operator
from os import path
import random

# Configuratoin settings to be moved externally later
email = "jon.demasi@colorado.edu"

"""
Take a search term, intended to be an rs#,
and return a list of pmids
"""
def get_pmids(interm):
    interm = interm + " AND pubmed_snp_cited[sb]"
    Entrez.email = email     # Always tell NCBI who you are
    search_results = Entrez.read(Entrez.esearch(db="pubmed", term=interm, usehistory="y"))
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

def tokenize_abstracts(abstracts):

    """ Takes a list of abstracts and breaks up each abstract into tokens """

    tokenized_abstracts_list = []
    for abstract in abstracts:
        tokens = nltk.word_tokenize(abstract)
        tokenized_abstracts_list.append(tokens)
    return tokenized_abstracts_list

def tagged_abstracts(tokenized_abstracts_list):

    """ Takes a list of tokenized abstracts
    and tags them using the NLTK module for Natural Language Entities"""

    tagged_abstracts_list = []
    for tokenized_abstract in tokenized_abstracts_list:
        tagged = nltk.pos_tag(tokenized_abstract)
        tagged_abstracts_list.append(tagged)
    return tagged_abstracts_list

def extract_nouns(tagged_abstracts_list, def_tags_per_abs = 0.0):

    """Takes a list of tuples of the form (word, tag) and returns a dictionary of counts for each
    word with tag "NN", "NNS", "NNP" or "NNPS" """

    noun_counter = []
    all_abstract_noun_counts = []
    normalized_all_counts = {}


    for tags in tagged_abstracts_list:

        per_abstract_noun_counts = []

        for tag in tags:

            if tag[1] == "NN" or tag[1] == "NNS" or tag[1] == "NNP" or tag[1] == "NNPS":


                per_abstract_noun_counts.append(str(tag[0].encode('ascii', 'ignore')))

                noun_counter.append(str(tag[0].encode('ascii', 'ignore')))


        all_abstract_noun_counts.append(dict(Counter(per_abstract_noun_counts)))

    all_counts = dict(Counter(noun_counter))

    num_abstracts = float(len(tagged_abstracts_list))

    for key in all_counts.keys():

        total_occurrences = float(all_counts[key])

        for each_abstract in all_abstract_noun_counts:

            if key in each_abstract:

                single_abstract_count = float(each_abstract[key])

                if (single_abstract_count/total_occurrences) < def_tags_per_abs:

                    normalized_all_counts[key] = float(all_counts[key])/num_abstracts

    return normalized_all_counts

#test implementation

toquery=get_pmids("rs328")

abstracts = get_abstracts(toquery)

t = tokenize_abstracts(abstracts)
q = tagged_abstracts(t)
n = extract_nouns(q)

print create_wordcloud(n)
