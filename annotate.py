#!/usr/bin/env python2

from Bio import Entrez
import time  # We need this if we want to filter based on age
from urllib2 import HTTPError
import xml.etree.ElementTree as ET
#from settings import configuration

def annotate_RS_ID(RS_ID, email):
    """This script should return annotation asspciate with each RS ID. Can incorporate more information later"""
    Entrez.email = email
    toreturn = {}
    for each_RS_ID in RS_ID:
        handle = Entrez.esummary(db="snp", id=each_RS_ID)
        record = handle.read()
        handle.close()
        root = ET.fromstring(record)
        result = []
        for item in root.iter('Item'):
            # print x.attrib, x.text
            attribute = item.attrib
            if attribute['Name'] == 'GLOBAL_MAF':
                result.append(item.text)
            if attribute['Name'] == 'CLINICAL_SIGNIFICANCE':
                result.append(item.text)
            if attribute['Name'] == 'CHR':
                result.append(item.text)
            if attribute['Name'] == 'FXN_CLASS':
                result.append(item.text)
        toreturn[each_RS_ID] = result
    return toreturn

RS_ID_list = [328, 45468592, 193922659, 199474738]
results = annotate_RS_ID(RS_ID_list, "tnphung@ucla.edu")
header = ["RS_ID", "Alt_AF", "clinical_sig", "chr", "func"]
print "\t".join(item for item in header)
for k, v in results.iteritems():
    print k, v[0], v[1], v[2], v[3]
    # toprint = [k, v[0], v[1], v[2], v[3]]
    # if toprint[0] != 'None':
    #     print "\t".join(item for item in toprint)
    # else:
    #     toprint = ['NA', v[0], v[1], v[2], v[3]]
    #     print "\t".join(item for item in toprint)