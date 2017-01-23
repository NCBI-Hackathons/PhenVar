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
from settings import configuration
email = configuration["email"]

class PhenVar_WordCloud(object):
    """Takes a single or list of rsids and gives a wordcloud for all rsids
    or each rsid or both"""

    def __init__(self,interm):
        self.interm = interm

    def get_complete_rsids(self):

        """
        Generate a list of rsids that are explicitly cited
        in a pubmed paper somewhere :-)
        """
        Entrez.email = email
        rsidlist = []
        numresults = 0
        retstart = 0
        search_string = "snp_pubmed_cited[sb]"
        #search_results = Entrez.read(Entrez.esearch(db="snp", term=search_string,retmax=100000, retstart=retstart, usehistory="y"))
        search_results = Entrez.read(Entrez.esearch(db="snp", term=self.interm,retmax=100000, retstart=retstart, usehistory="y"))
        print("Found a total of " +
              search_results["Count"] + " results using search string '" + search_string + "'")
        numresults = search_results["Count"]
        rsidlist = rsidlist + search_results["IdList"]
        additional_queries = int(int(numresults) / 100000)
        while additional_queries != 0:
            retstart = retstart + 100000
            search_results = Entrez.read(Entrez.esearch(db="snp", term=search_string,
                                                        retmax=100000, retstart=retstart, usehistory="y"))
            rsidlist = rsidlist + search_results["IdList"]
            additional_queries = additional_queries - 1

        self.rsidlist = rsidlist
        return

    def get_pmids(self):

        """
        Take a search term, intended to be an rs#,
        and return a dictionary that can be used to
        efetch later
        """

        # This is obsolete now, essentially, but
        # allows a user to pass a single string
        # which can be nice.
        if isinstance(self.interm, str):
            interm = interm + " AND pubmed_snp_cited[sb]"
            Entrez.email = email     # Always tell NCBI who you are
            search_results = Entrez.read(Entrez.esearch(db="pubmed", term=self.interm,
                                                        retmax=100000,
                                                        usehistory="y"))
            print("Found a total of " +
                  search_results["Count"] + " results using search string '" + interm + "'")

            self.search_results = search_results
            return

        elif isinstance(self.interm, list):
            searchstring = " OR ".join(interm)
            searchstring = "(" + searchstring + ") AND pubmed_snp_cited[sb]"
            Entrez.email = email     # Always tell NCBI who you are
            search_results = Entrez.read(Entrez.esearch(db="pubmed",
                                                        term=searchstring,
                                                        retmax=100000,
                                                        usehistory="y"))
            print("Found a total of " +
                  search_results["Count"] + " results using search string'" + searchstring + "'")

            self.search_results = search_results

            return

    def get_abstracts(self):
        """
        Takes the saved list of results from get_pmids and retrieves
        the raw XML to parse.  Currently goes article by article instead of
        doing them in bulk.  Could go either way. Unsure which way
        is better / more efficient.
        """
        abstracts_list = []
        Entrez.email = email

        for start in range(0, int(self.search_results["Count"]), 1):
            end = min(int(self.search_results["Count"]), start + 1)
            # print("Going to download record %i to %i" % (start+1, end))
            fetch_handle = Entrez.efetch(db="pubmed", rettype="abstract",
                                            retmode="xml", retstart=start,
                                            retmax=1,
                                            webenv=self.search_results["WebEnv"],
                                            query_key=self.search_results["QueryKey"])
            data = fetch_handle.read()
            fetch_handle.close()
            root = ET.fromstring(data)
            for abst in root.iter('Abstract'):
                for sec in abst.iter('AbstractText'):
                    abstracts_list.append(sec.text)

            self.abstracts_list = abstracts_list
        return


    def get_abstracts_from_list(self):
        abstracts_list = []
        pmids_abstracts_dict = {}
        Entrez.email = email
        #for each_pmid in pmids_list:
        for each_pmid in self.search_results:
            fetch_handle = Entrez.efetch(db="pubmed", id=each_pmid, retmode='xml')
            data = fetch_handle.read()
            fetch_handle.close()
            root = ET.fromstring(data)
            for abst in root.iter('Abstract'):
                for sec in abst.iter('AbstractText'):
                    abstracts_list.append(sec.text)

        self.abstracts_list = abstracts_list
        return


    def tokenize_abstracts(self):
        """ Takes a list of abstracts and breaks up each abstract into tokens """

        tokenized_abstracts_list = []
        #for abstract in abstracts:
        for abstract in self.abstracts_list:
            tokens = nltk.word_tokenize(abstract)
            tokenized_abstracts_list.append(tokens)

        self.tokenized_abstracts_list = tokenized_abstracts_list
        return


    def tagged_abstracts(self):
        """ Takes a list of tokenized abstracts
        and tags them using the NLTK module for Natural Language Entities"""

        tagged_abstracts_list = []
        for tokenized_abstract in self.tokenized_abstracts_list:
            tagged = nltk.pos_tag(tokenized_abstract)
            tagged_abstracts_list.append(tagged)

        self.tagged_abstracts_list = tagged_abstracts_list
        return


    def extract_nouns(self):
        """Takes a list of tuples of the form (word, tag) and returns a dictionary of counts for each
        word with tag "NN", "NNS", "NNP" or "NNPS" """

        noun_counter = []
        all_abstract_noun_counts = []
        normalized_all_counts = {}

        for tags in self.tagged_abstracts_list:

            per_abstract_noun_counts = []

            for tag in tags:

                if tag[1] == "NN" or tag[1] == "NNS" or tag[1] == "NNP" or tag[1] == "NNPS":

                    per_abstract_noun_counts.append(
                        str(tag[0].encode('ascii', 'ignore')))

                    noun_counter.append(str(tag[0].encode('ascii', 'ignore')))

            all_abstract_noun_counts.append(
                dict(Counter(per_abstract_noun_counts)))

        all_counts = dict(Counter(noun_counter))

        num_abstracts = float(len(tagged_abstracts_list))

        for key in all_counts.keys():

            total_occurrences = float(all_counts[key])

            for each_abstract in all_abstract_noun_counts:

                if key in each_abstract:

                    single_abstract_count = float(each_abstract[key])

                    if def_tags_per_abs != 0:

                        if (single_abstract_count / total_occurrences) < def_tags_per_abs:

                            normalized_all_counts[key] = float(
                                all_counts[key]) / num_abstracts

                    else:

                        normalized_all_counts[key] = float(
                            all_counts[key]) / num_abstracts

        self.normalized_all_counts = normalized_all_counts
        return


    def extract_nouns_filter(self):
        """Takes a list of tuples of the form (word, tag) and returns a
        dictionary of counts for each
        word with tag "NN", "NNS", "NNP" or "NNPS.
        This function is different from the function extract_nouns in that
        this function is doing some external filtering.
        Right now, it will filter using a predefined list of termed defined
        from using textools (TODO: add more info here)" """

        filter_list = ['polymorphism', 'polymorphisms', 'nucleotide', 'nucleotides', 'snp', 'snps', 'allele', 'alleles', 'gene', 'genes',
                       'genotype', 'genotypes', 'genotyped',  'single', 'singles', 'genetic', 'genetics', 'study', 'studies', 'variant', 'variants', 'analysis', 'analyses']
        noun_counter = []
        all_abstract_noun_counts = []
        normalized_all_counts = {}

        for tags in self.tagged_abstracts_list:

            per_abstract_noun_counts = []

            for tag in tags:

                if tag[1] == "NN" or tag[1] == "NNS" or tag[1] == "NNP" or tag[1] == "NNPS":
                    if tag[0].lower() not in filter_list:

                        per_abstract_noun_counts.append(
                            str(tag[0].encode('ascii', 'ignore')))

                        noun_counter.append(str(tag[0].encode('ascii', 'ignore')))

            all_abstract_noun_counts.append(
                dict(Counter(per_abstract_noun_counts)))

        all_counts = dict(Counter(noun_counter))

        num_abstracts = float(len(tagged_abstracts_list))

        for key in all_counts.keys():

            total_occurrences = float(all_counts[key])

            for each_abstract in all_abstract_noun_counts:

                if key in each_abstract:

                    single_abstract_count = float(each_abstract[key])

                    if def_tags_per_abs != 0:

                        if (single_abstract_count / total_occurrences) < def_tags_per_abs:

                            normalized_all_counts[key] = float(
                                all_counts[key]) / num_abstracts

                    else:

                        normalized_all_counts[key] = float(
                            all_counts[key]) / num_abstracts

        self.normalized_all_counts = normalized_all_counts

        return 
