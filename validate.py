import nltk
from Bio import Entrez
import time  # We need this if we want to filter based on age
from urllib2 import HTTPError
import xml.etree.ElementTree as ET
from collections import Counter
import numpy as np
import operator
from os import path
import random

import settings
import ncbiutilsTanya as ncbiutils
import lanpros
import argparse
import wordcloudfornouns
import argparse
import csv


def parse_args():
    """
    Parse command-line arguments
    """
    parser = argparse.ArgumentParser(
        description="This script parses command-line arguments and will run PhenVar pipeline: ncbiutils, lanpros, and wordcloudfornouns")

    parser.add_argument(
        "--RS_ID", nargs="+", required=True,
        help="REQUIRED. RS ID of interest. Could input one RS ID or more than one RS IDs. If there are more than 1 RS ID, just separate by space")
    parser.add_argument(
        "--prior_dict_filename", required=True,
        help="REQUIRED. Pre-filtering: Name of the output file for the output dictionary - prefilter")
    parser.add_argument(
        "--prior_plot_filename", required=True,
        help="REQUIRED. Pre-filtering: Name of the output file for the plot - prefilter")
    parser.add_argument(
        "--post_dict_filename", required=True,
        help="REQUIRED. Post-filtering: Name of the output file for the output dictionary - postfilter")
    parser.add_argument(
        "--post_plot_filename", required=True,
        help="REQUIRED. Post-filtering: Name of the output file for the plot - postfilter")
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    RS_ID = args.RS_ID
#    print "Getting pmids from RS ID..."
#    pmids = ncbiutils.get_pmids(RS_ID)
    pmids = [15146469, 15060124, 8834249]
    print "Getting abstracts from pmids..."
    abstracts = ncbiutils.get_abstracts_from_list(pmids)
    print "Analyze using language process script..."
    tokens = lanpros.tokenize_abstracts(abstracts)
    tagged_abstracts = lanpros.tagged_abstracts(tokens)

    print "Getting results PRIOR to any filtering of key words"
    prior_results = lanpros.extract_nouns(tagged_abstracts)
    with open(args.prior_dict_filename, "w") as f:
        output_list = []
        for key in sorted(prior_results):
            output_list.append([key, prior_results[key]])
        w = csv.writer(f)
        w.writerows(output_list)

    wordcloudfornouns.create_wordcloud(prior_results, args.prior_plot_filename)

    print "Getting results POST filtering of key words"
    post_results = lanpros.extract_nouns_filter(tagged_abstracts)
    with open(args.post_dict_filename, "w") as f:
        output_list = []
        for key in sorted(post_results):
            output_list.append([key, post_results[key]])
        w = csv.writer(f)
        w.writerows(output_list)
    wordcloudfornouns.create_wordcloud(post_results, args.post_plot_filename)
main()


# RS_ID = "rs328"
# pmids = ncbiutils.get_pmids(RS_ID)
# #clinvar_pmids = [11389159,15695382,18607349,20104584,21520273,22703879,24348212,24728327,25637381,25741868,26586665,8896551,9971877]
# #abstracts = ncbiutils.get_abstracts(pmids)

# abstracts = ncbiutils.get_abstracts(pmids)
# tokens = lanpros.tokenize_abstracts(abstracts)
# tagged_abstracts = lanpros.tagged_abstracts(tokens)
# #results = lanpros.extract_nouns_filter(tagged_abstracts)
# results = lanpros.extract_nouns(tagged_abstracts)
# wordcloudfornouns.create_wordcloud(results, "figure.png")

# for k, v in results.iteritems():
#    print k, v
# main()
