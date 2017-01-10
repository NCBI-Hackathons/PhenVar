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

import settings
import ncbiutils
import lanpros
import argparse

# def parse_args():
#     """
#     Parse command-line arguments
#     """
#     parser = argparse.ArgumentParser(description="This script parses command-line arguments and will run the language processing script.")

#     parser.add_argument(
#             "--RS_ID", required=True,
#             help="REQUIRED. RS ID of interest. Could input one RS ID or more than one RS IDs")
#     args = parser.parse_args()
#     return args

# def main():
RS_ID = "rs328"
pmids = ncbiutils.get_pmids(RS_ID)
#clinvar_pmids = [11389159,15695382,18607349,20104584,21520273,22703879,24348212,24728327,25637381,25741868,26586665,8896551,9971877]
#abstracts = ncbiutils.get_abstracts(pmids)

abstracts = ncbiutils.get_abstracts_from_list(pmids)
#tokens = lanpros.tokenize_abstracts(abstracts)
#tagged_abstracts = lanpros.tagged_abstracts(tokens)
#results = lanpros.extract_nouns(tagged_abstracts)


#for k, v in results.iteritems():
#    print k, v
# main()
