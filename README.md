# PhenVar
PhenVar is designed to take one or more rsids and generate a list of PubMed IDs to query and generate novel associations between publications. It utilizes a local sqlite cache in a configurable location to keep a local copy of relevant srids, pmids, and the abstrat blobs for the pmids.  
## Installation Notes
At the time of writing, the PyQt4 and SIP packages were broken in pip, so we had to compile them from source.  They make break a pip -r requirements.txt until that issue is resolved.  Otherwise, these scripts are meant to basically be standalone utilities that can be mixed and matched together to work in whatever fashion you are hoping for.  The validate.py can be run to test that all of the functionality is working, you should get two tables and two word clouds as output.  Please utilize the issue tracking on GitHub if you experience any issues or have feature requests.  
### settings.py
settings.py is a file containing a single dictionary with configuration options that is imported into all of the other tools.  Configuration options you need can easily be added or changed.
### ncbiutils.py
#### get_pmids 
Takes a single rsid input (string) or several rsid inputs (list of strings) and returns a dictionary from Entrez.read/esearch.  The results returned are all pmids which explicitly cited the rsids given.  
#### get_abstracts
Expects the results from get_pmids as an input.  Will return a list of abstracts.  Each item in the list is an abstract to a pmid from the get_pmids search results.  Can return a list with only one item.  
#### get_abstracts_from_list
Instead of relying on a previous, saved search, this function will return the abstracts from a pre-defined list of pmids
### db.py
#### connect
Plain and simply connects to an sqlite3 database and opens a cursor. Function returns a list where index 0 is the connection name and index 1 is the cursor object.
#### disconnect
commits any in-RAM changes and then closes the sqlite3 database
#### insert_date
Assumes that the updatehist table is created and inserts a record with the supplied number of "records."  This is meant to be run upon database creation and updates to keep an update history.  Can certainly be improved to contain more pertinent information
#### print_update_history
Prints a mostly human-readable form of the update table
#### check_db
Checks if a given database location exists -- useful to use before running an update hook or upon initial loading of program.  
#### create_cache
A little convoluted right now.  Assumes the db/tables don't exist (so use check_db first) and creates them.  Inserts all values from esearch for rsid/pmids and pmids/abstracts.   Not currently primary key based, so duplicate rows will exist primarily in the pmid/abstracts table.  
#### check_updates
Placeholder function that currently doesn't work.  Will check if there are new additions for any of the tables after the initial cache is already created.  
#### get_pmids
Given an rsid input, return all pmids that explicitly cite that rsid in a python list
#### get abstracts
Placeholder.  Intended purpose is to pass in a list of pmids and get back a list of abstracts. 
#### initdb
Used to automatically check_db and then either run init or update, then disconnect.
