#!/usr/bin/env python2

from Bio import Entrez
import time # We need this if we want to filter based on age
from urllib2 import HTTPError
import xml.etree.ElementTree as ET
import argparse
import csv

def parse_args():
    """
    Parse command-line arguments
    """
    parser = argparse.ArgumentParser(description="This script grabs information from esummary. Currently this script will grab GLOBAL_MAF, CLINICAL_SIGNIFICANCE, CHR, FXN_CLASS")

    parser.add_argument(
            "--RS_ID", required=True,
            help="REQUIRED. RS ID of interest")
    parser.add_argument(
            "--email_address", required=True,
            help="REQUIRED. Email address is needed to call eutils")
    parser.add_argument("--outfile", required=True, 
            help="REQUIRED. Name of output file.")
    args = parser.parse_args()
    return args

def annotate_RS_ID(RS_ID, email):
    Entrez.email = email
    handle = Entrez.esummary(db="snp", id=RS_ID)
    record = handle.read()
    handle.close()
    root = ET.fromstring(record)
    result = []
    for item in root.iter('Item'):
        # print x.attrib, x.text
        attribute = item.attrib
        if attribute['Name'] == 'GLOBAL_MAF':
            result.append(item.text)
        if attribute['Name'] == 'CLINNICAL_SIGNIFICANCE':
            result.append(item.text)
        if attribute['Name'] == 'CHR':
            result.append(item.text)
        if attribute['Name'] == 'FXN_CLASS':
            result.append(item.text)
    return result

def main():
    args = parse_args()
    RS_ID = str(args.RS_ID)
    email = str(args.email_address)
    result = annotate_RS_ID(RS_ID, email)
    with open(args.outfile,"w") as f:
        w = csv.writer(f)
        w.writerows(result)
main()


