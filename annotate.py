#!/usr/bin/env python2

from Bio import Entrez
import time  # We need this if we want to filter based on age
from urllib2 import HTTPError
import xml.etree.ElementTree as ET
from settings import configuration
email = configuration["email"]


def annotate_RS_ID(RS_ID, email):
    """This script should return annotation asspciate with each RS ID. Can incorporate more information later"""
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
# TODO: add header to output file

