import sqlite3
import ncbiutils
import time
from settings import configuration
import os

db_location = configuration["dbloc"]
db_exists = False


""" Establishes a connection to an sql database that may or may not
previous exist and returns a list including the connection as well
as the cursor """
def connect(db_location):
    conn = sqlite3.connect(db_location)
    cur = conn.cursor()
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
    print("Adding entry to update history: " + str(records)  + " records on " + now)
    cursor.execute("""INSERT INTO updatehist VALUES (?, ?)""", (now, records))
    db.commit()

""" This function can be used to effectively print the entire update history
for a given dbs update table """ 
def print_update_history(db, cursor):
    cursor.execute("""SELECT * FROM updatehist""")
    print(cursor.fetchall())
    pass

""" Check if db already exists.  This should always be run 
prior to running updates, because they do not currently check
on their own """
def check_db(location):
    global db_exists
    if os.path.isfile(location):
        db_exists = True

# Called if cache does not yet exist
def create_cache(conn, cursor):
    numadded = 0
    # Opens db and returns list of conn, cursur
    # Table for rsid and associated pmids
    cursor.execute('''CREATE TABLE rsids
             (rsid text, pmids text)''')
    # A table to store update history in case we'd need it?
    cursor.execute('''CREATE TABLE updatehist
             (date text, added integer)''')
    # Table for pmid + abstract 
    cursor.execute('''CREATE TABLE pminfo
             (pmid integer, abstract text)''')
    # Get list of rsids cited in pubmed
    # For each rsid in list, get pmids citing them
    # Update our date table with the records we just added
    list = ncbiutils.get_complete_rsids()
    spot = len(list) - 1
    stop = spot - 100
    while spot != stop:
        dict = ncbiutils.get_pmids("rs"+list[spot])
        cursor.execute("""INSERT INTO rsids VALUES(?, ?)""", (list[spot], ",".join(dict["IdList"])))
        for y in dict["IdList"]:
            ylist=[]
            ylist.append(y)
            newdict = ncbiutils.get_abstracts_from_list(ylist)
            print("ID: " + y)
            #print("Abstract: ")
            #print(newdict[y])
            if newdict:
                cursor.execute("""INSERT INTO pminfo VALUES (?, ?)""", (y, newdict[y]))
        time.sleep(5)
        spot = spot - 1

""" Check if there are updates to our current
state of cachedb.  If so, download and append
to relevant tables. """
def check_updates(conn, cursor):
    pass

""" For a given rsid, get the list of pmids that
explicitly cite it. """
def get_pmids(rsid):
    t = (rsid,)
    c.execute('SELECT pmids FROM rsids WHERE rsid=?', t)
    return pmids

""" For a given pmid, return its abstract. If multiple
pmids are given, return a list of all abstracts. """
def get_abstracts():
    pass

""" This function is intended to be run to establish a connection,
as well as a cursor.  It will also check if the DB exists or needs
to be updated (can be configured to ignore updates).  Returns a tuple
of (connection, cursor). 
CURRENTLY BEING USED AS A TEST FUNCTION"""
def initdb():
    check_db(db_location)
    if db_exists is True:
        print("Database already exists. Checking for updates")
        sqlinfo = connect(db_location)    
        check_updates(sqlinfo[0], sqlinfo[1])
        disconnect(sqlinfo[0])
    elif db_exists is False:
        print("Database does not exist. Creating now... This may take awhile")
        sqlinfo = connect(db_location)
        create_cache(sqlinfo[0], sqlinfo[1])
        disconnect(sqlinfo[0])


initdb()
