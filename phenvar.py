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
