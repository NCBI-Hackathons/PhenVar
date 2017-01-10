import sqlite3
import ncbiutils
import time
from settings import configuration

db_location = configuration["dbloc"]
db_exists = False

# Check if db already exists 
def checkdb(location):
    if os.path.isfile(location):
        db_exists = True

def createcache(location):
    # Since this doesn't exist, creates it
    cachedb = sqlite3.connect(location)
    c = cachedb.cursor()
    # Table for rsid and associated pmids
    c.execute('''CREATE TABLE rsids
             (rsid text, pmids text)''')
    # A table to store update history in case we'd need it?
    c.execute('''CREATE TABLE updatehist
             (date text, added integer)''')
    # Table for pmid + abstract 
    c.execute('''CREATE TABLE pmidabs
             (pmid integer, abstract text)''')
    cachedb.commit()
    cachdb.close()

""" Check if there are updates to our current
state of cachedb.  If so, download and append
to relevant tables. """
def check_updates(location):
    pass

def main():
    if db_exists(db_location):
        check_updates(db_location)
    else:
        createcache(db_location)

main()
