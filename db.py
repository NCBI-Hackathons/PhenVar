import sqlite3
import ncbiutils
import time
from settings import configuration

db_location = configuration["dbloc"]
db_exists = False


""" Establishes a connection to an sql database that may or may not
previous exist and returns a list including the connection as well
as the cursor """
def connect(db_location):
    conn = sqlite3.connect(db_location)
    cur = conn.cursur()
    return(conn, cur)

""" Commits all changes that may not have actually been written yet
and closes the connection.  Should be used at the end of all db 
operations """
def disconnect(conn):
    conn.commit()
    conn.close()
    
""" Dedicated to inserting date in update table on creation/update
pulled out in its own function in case we want to store metdata
about updates later. """
def insert_date(db, cursor, records):
    now = time.strftime("%Y-%m-%d")
    cursor.execute("""INSERT INTO updatehist VALUES (?, ?)""", (now, records))
    db.commit()

# Check if db already exists 
def check_db(location):
    if os.path.isfile(location):
        db_exists = True

# Called if cache does not yet exist
def create_cache(location):
    numadded = 0
    # Since this doesn't exist, creates it
    
    # Table for rsid and associated pmids
    c.execute('''CREATE TABLE rsids
             (rsid text, pmids text)''')
    # A table to store update history in case we'd need it?
    c.execute('''CREATE TABLE updatehist
             (date text, added integer)''')
    # Table for pmid + abstract 
    c.execute('''CREATE TABLE pminfo
             (pmid integer, abstract text)''')
    # Get list of rsids cited in pubmed
    # For each rsid in list, get pmids citing them
    # Update our date table with the records we just added
    insert_date(location, numadded)
    # Replace the following two with a dedicated function
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

""" This function is intended to be run to establish a connection,
as well as a cursor.  It will also check if the DB exists or needs
to be updated (can be configured to ignore updates).  Returns a tuple
of (connection, cursor). """
def initdb():
    # Regardless of which of the following conditions is true 
    # we're going to need to open db and cursor, so let's just
    # go ahead and do it.  
    if db_exists(db_location):
        check_updates(db_location)
    else:
        createcache(db_location)

