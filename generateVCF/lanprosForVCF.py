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

def tokenize_abstracts(RS_pmids_abstracts_dict):

    """ Takes in a dictionary where key is the RSID, value is a dictionary where key is the pmid and value is abstract. Breaks up each abstract into tokens. Return: dictionary where key is RSID, value is a dictionary where key is the pmid and value is a list of tokenized abstracts """
    RS_pmids_tokenizedabstracts_dict = {}
    for each_RS in RS_pmids_abstracts_dict:
        pmids_tokenizedabstracts = {}
        pmids_abstracts = RS_pmids_abstracts_dict[each_RS]
        for pmid in pmids_abstracts:
            tokenizedabtracts_list = []
            tokens = nltk.word_tokenize(pmids_abstracts[pmid])
            tokenizedabtracts_list.append(tokens)
            pmids_tokenizedabstracts[pmid] = tokenizedabtracts_list
        RS_pmids_tokenizedabstracts_dict[each_RS] = pmids_tokenizedabstracts
    return RS_pmids_tokenizedabstracts_dict

def tagged_abstracts(RS_pmids_tokenizedabstracts_dict):

    """ Takes a dict of tokenized abstracts
    and tags them using the NLTK module for Natural Language Entities.
    Input dictionary: key is the RS ID, value is a dictionary where key is the pmid and value is a list of tokens"""
    RS_pmids_taggedabstracts_dict = {}
    for each_RS in RS_pmids_tokenizedabstracts_dict:
        pmids_taggedabstracts = {}
        pmids_tokenizedabstracts = RS_pmids_tokenizedabstracts_dict[each_RS]
        for pmid in pmids_tokenizedabstracts:
            taggedabstracts_list = []
            for token in pmids_tokenizedabstracts[pmid]:
                tagged = nltk.pos_tag(token)
                taggedabstracts_list.append(tagged)
            pmids_taggedabstracts[pmid] = taggedabstracts_list
        RS_pmids_taggedabstracts_dict[each_RS] = pmids_taggedabstracts
    return RS_pmids_taggedabstracts_dict

def extract_nouns(RS_pmids_taggedabstracts_dict):

    """Takes a dict where key is the RS ID, values is a dict where key is the pmid and value is list of tuples of the form (word, tag).
    Return a dictionary of counts for each
    word with tag "NN", "NNS", "NNP" or "NNPS" """

    # noun_counter = []
    # all_abstract_noun_counts = []
    # normalized_all_counts = {}

    RS_pmids_nounsabstracts_dict = {}
    for each_RS in RS_pmids_taggedabstracts_dict:
        pmids_nounsabstracts = {}
        pmids_taggedabstracts = RS_pmids_taggedabstracts_dict[each_RS]
        for pmid in pmids_taggedabstracts:
            nounsabstracts = []
            for tags in pmids_taggedabstracts[pmid]:
                for tag in tags:
                    if tag[1] == "NN" or tag[1] == "NNS" or tag[1] == "NNP" or tag[1] == "NNPS":
                        nounsabstracts.append(str(tag[0].encode('ascii', 'ignore')))
            nounsabstracts_count = dict(Counter(nounsabstracts))
            nouns = sorted(nounsabstracts_count, key=nounsabstracts_count.__getitem__, reverse=True)

            #nouns = nounsabstracts_count.keys()
            # pmids_nounsabstracts[pmid] = nounsabstracts_count
            pmids_nounsabstracts[pmid] = nouns
        RS_pmids_nounsabstracts_dict[each_RS] = pmids_nounsabstracts

    return RS_pmids_nounsabstracts_dict

# def obtain_all_abtracts_counts(RS_pmids_nounsabstracts_dict):
#     all_abstract_noun_counts = {}
#     for each_RS in RS_pmids_nounsabstracts_dict:
#         pmids_nounsabstracts_dict = RS_pmids_nounsabstracts_dict[each_RS]
#         for pmid in pmids_nounsabstracts_dict:
#             for each_noun in pmids_nounsabstracts_dict[pmid]:
#                 if each_noun not in all_abstract_noun_counts.keys():
#                     all_abstract_noun_counts[each_noun] = pmids_nounsabstracts_dict[pmid][each_noun]
#                 else:
#                     all_abstract_noun_counts[each_noun] += pmids_nounsabstracts_dict[pmid][each_noun]
#     print len(all_abstract_noun_counts)



    # for tags in taggedabstracts_list:

    #     per_abstract_noun_counts = []

    #     for tag in tags:

    #         if tag[1] == "NN" or tag[1] == "NNS" or tag[1] == "NNP" or tag[1] == "NNPS":


    #             per_abstract_noun_counts.append(str(tag[0].encode('ascii', 'ignore')))

    #             noun_counter.append(str(tag[0].encode('ascii', 'ignore')))


    #     all_abstract_noun_counts.append(dict(Counter(per_abstract_noun_counts)))

    # all_counts = dict(Counter(noun_counter))

    # num_abstracts = float(len(taggedabstracts_list))

    # for key in all_counts.keys():

    #     total_occurrences = float(all_counts[key])

    #     for each_abstract in all_abstract_noun_counts:

    #         if key in each_abstract:

    #             single_abstract_count = float(each_abstract[key])

    #             # if def_tags_per_abs != 0:

    #             #     if (single_abstract_count/total_occurrences) < def_tags_per_abs:

    #             #         normalized_all_counts[key] = float(all_counts[key])/num_abstracts

    #             # else:
                    
    #             #     normalized_all_counts[key] = float(all_counts[key])/num_abstracts


    # return normalized_all_counts
