#!/usr/bin/env python3
from database import initialize_database, load_all_data, load_data_not_present, update_data, build_complete_noun_mapping,\
    build_rsid_list
from sys import argv


help_text = """USAGE:
./manage.py <command>
 commands:
    - initialize    creates tables in database, should only be run once
    - load          loads all snp cited articles from pubmed
    - resume        if load fails partway through, loads any articles that were not added to the database
    - update        loads any articles that have been added/modified since latest update
    - build_json    builds json noun cache from articles in database
"""

def help():
    print(help_text)


def build_json():
    build_complete_noun_mapping()
    build_rsid_list()


commands = {
    "--help": help,
    "-h": help,
    "initialize": initialize_database,
    "load": load_all_data,
    "resume": load_data_not_present,
    "update": update_data,
    "build_json": build_json,
}

if len(argv) == 2 and argv[1] in commands:
    commands[argv[1]]()
else:
    help()
