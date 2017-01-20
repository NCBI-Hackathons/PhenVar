
from __future__ import print_function
#!/usr/bin/env python2
from Bio import Entrez
import time  # We need this if we want to filter based on age
from urllib2 import HTTPError
import xml.etree.ElementTree as ET
import os
import sys


def annotate_RS_ID(RS_ID, email):
    """This script returns annotation from eutils esummary associated with each RS ID. Currently it will return CLINICAL_SIGNIFICANCE, FXN_CLASS, VALIDATED.
    TODO: would be better to not hardcode these options but to be user-defined."""
    Entrez.email = email
    handle = Entrez.esummary(db="snp", id=RS_ID)
    record = handle.read()
    handle.close()
    root = ET.fromstring(record)
    result = []
    for item in root.iter('Item'):
        attribute = item.attrib
        if attribute['Name'] == 'CLINICAL_SIGNIFICANCE':
            result.append(item.text)
        if attribute['Name'] == 'FXN_CLASS':
            result.append(item.text)
        if attribute['Name'] == 'VALIDATED':
            result.append(item.text)
    return result

# Incorporate information from esummary and also from the file PubMedAlleleFreqByRs1000G_All.txt obtained from Lon Phan. This file has allele frequency information for (1) global maf, (2) AFR, (3) AMR, (4) EAS, (5) EUR, and (6) SAS

email='tnphung@ucla.edu'
output = []
header = ['#chr', 'pos', 'rs_id', 'GMAF', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'clinical_sig', 'fxn_class', 'validated']
output.append(header)
# Load the file PubMedAlleleFreqByRs1000G_All_tab.txt from command line
with open(sys.argv[1], "r") as file:
	for line in file:
		line = line.rstrip("\n")
		line = line.split("\t")
		each_RS_esummary = annotate_RS_ID(int(line[2]), email)
		output.append([line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], each_RS_esummary[0], each_RS_esummary[1], each_RS_esummary[2]])

for i in output:
	print(*i, sep=' ', end='\n')


# RS_ID_list = [328, 45468592, 193922659, 199474738]
# results = annotate_RS_ID(RS_ID_list, "tnphung@ucla.edu")
# header = ["RS_ID", "Alt_AF", "clinical_sig", "chr", "func"]
# print "\t".join(item for item in header)
# for k, v in results.iteritems():
#     print k, v[0], v[1], v[2], v[3]
    # toprint = [k, v[0], v[1], v[2], v[3]]
    # if toprint[0] != 'None':
    #     print "\t".join(item for item in toprint)
    # else:
    #     toprint = ['NA', v[0], v[1], v[2], v[3]]
    #     print "\t".join(item for item in toprint)
