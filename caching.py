import sqlite3
import ncbiutils
import time
from settings import configuration

db_location = configuration["dbloc"]
db_exists = False

# Dedicated to inserting date in update table on creation/update
# pulled out in its own function in case we want to store metdata
# about updates later.
def insert_date(location, records):
    pass

# Check if db already exists 
def check_db(location):
    if os.path.isfile(location):
        db_exists = True

# Called if cache does not yet exist
def create_cache(location):
    numadded = 0
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
    c.execute('''CREATE TABLE pminfo
             (pmid integer, abstract text)''')
    insert_date(location, numadded)
    cachedb.commit()
    cachdb.close()

""" Check if there are updates to our current
state of cachedb.  If so, download and append
to relevant tables. """
def check_updates(cursor):
    pass

""" For a given rsid, get the list of pmids that
explicitly cite it. """
def get_pmids(rsid):
    t = (rsid,)
    c.execute('SELECT pmids FROM rsids WHERE rsid=?', t)
    return pmids

""" For a given pmid, return its abstract. If multiple
pmids are given, return a list of all abstracts. """
def get_abstracts:
    pass

def main():
    if db_exists(db_location):
        check_updates(db_location)
    else:
        createcache(db_location)

main()
