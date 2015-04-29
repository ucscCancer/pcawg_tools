#!/usr/bin/env python
"""
Tracker for job allocation for workloads distributed amoung multiple compute centers

Configuration:
All command lines need the following defined:
 - table_id : The synapse ID for the table the holds the registry
 - primary_col : The primary id column for the assignment: this must be unique
 - assignee_col : The column to hold the name of the assignee (value is your Synapse username)
 - state_col : Column to hold the state of the job

SynQueue will search for the file '.synqueue' up the directory hiearchy, starting
in the current directory, and stopping when it finds a file
ie, /home/myname/project_folder/test1/.synqueue -> /home/myname/project_folder/.synqueue ->
/home/myname/.synqueue -> /home/.synqueue

Example .synqueue file:
{
    "table_id" : "syn3498886",
    "primary_col" : "Donor_ID",
    "assignee_col" : "Assignee",
    "state_col" : "Processing State"
}


Usage:

Register for 5 assignments
./synqueue.py register -c 5

List your assignments
./synqueue.py list

Set Assignment state
./synqueue.py set live 178b28cd-99c3-48dc-8d09-1ef71b4cee80


"""
import argparse
import os
import subprocess
import sys
import synapseclient
import json
from synapseclient.exceptions import *
import itertools
import collections
import uuid
import csv

#HACK: fixing weirdness on older installs
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()


def listAssignments(syn, list_all, table_id, primary_col, assignee_col, state_col, debug):
    table = syn.get(table_id)
    if table.entityType != "org.sagebionetworks.repo.model.table.TableEntity":
        return
    username = syn.getUserProfile()['userName']
    results = syn.tableQuery('select "%s", "%s", "%s" from %s' % (primary_col, assignee_col, state_col, table.id))
    total = 0
    for row in results:
        if row[3] == username or list_all:
            print "\t".join( list( str(a) for a in row[2:]) )
            total += 1
    print "Total:", total

def registerAssignments(syn, count, table_id, primary_col, assignee_col, state_col, debug):
    table = syn.get(table_id)
    if table.entityType != "org.sagebionetworks.repo.model.table.TableEntity":
        return
    username = syn.getUserProfile()['userName']

    results = syn.tableQuery('select * from %s where "%s" is NULL limit %s' % (table.id, assignee_col, count))
    df = results.asDataFrame()
    for row in df.index:
        df.loc[row,assignee_col] = username
    syn.store(synapseclient.Table(table, df, etag=results.etag))

def setStates(syn, state, ids, table_id, primary_col, assignee_col, state_col, debug):
    table = syn.get(table_id)
    if table.entityType != "org.sagebionetworks.repo.model.table.TableEntity":
        return
    username = syn.getUserProfile()['userName']
    query = 'select * from %s where "%s" = \'%s\'' % (table.id, assignee_col, username)
    results = syn.tableQuery(query)
    df = results.asDataFrame()
    for row in df.index:
        if df.loc[row,primary_col] in ids:
            df.loc[row,state_col] = state
    syn.store(synapseclient.Table(table, df, etag=results.etag))

def add_config_arguments(parser):
    parser.add_argument("-t", "--table_id")
    parser.add_argument("--primary_col", default="id")
    parser.add_argument("--assignee_col", default="assignee")
    parser.add_argument("--state_col", default="state")


def build_parser():
    """Builds the argument parser and returns the result."""

    parser = argparse.ArgumentParser(
            description='Interfaces with the Synapse repository for sample tracking')
    parser.add_argument('--debug', dest='debug', action='store_true')

    subparsers = parser.add_subparsers( title='commands',
                                        description='The following commands are available:',
                                        help='For additional help: "synqueue <COMMAND> -h"')

    parser_register = subparsers.add_parser('register',
                                       help='Returns set of new assignments')
    parser_register.add_argument("-c", "--count", help="Number of assignments to get", type=int, default=1)
    add_config_arguments(parser_register)
    parser_register.set_defaults(func=registerAssignments)

    parser_list = subparsers.add_parser('list')
    parser_list.add_argument("-a", "--list_all", action="store_true", default=False)
    add_config_arguments(parser_list)
    parser_list.set_defaults(func=listAssignments)

    parser_set = subparsers.add_parser('set')
    add_config_arguments(parser_set)
    parser_set.add_argument("state")
    parser_set.add_argument("ids", nargs="+")
    parser_set.set_defaults(func=setStates)

    return parser



if __name__ == "__main__":
    args = build_parser().parse_args()
    syn = synapseclient.Synapse(debug=args.debug, skip_checks=True)
    syn.login(silent=True)
    if 'func' in args:
        kwds = dict(vars(args))
        del kwds['func']
        cur_dir = os.path.abspath(os.getcwd())
        while True:
            cur_path = os.path.join(cur_dir,".synqueue")
            print "Checking", cur_path
            if os.path.exists(cur_path):
                with open(cur_path) as handle:
                    meta = json.loads(handle.read())
                for k,v in meta.items():
                    kwds[k] = v
                break
            cur_dir = os.path.dirname(cur_dir)
            if cur_dir == "/":
                break

        for i in ['table_id', 'primary_col', 'assignee_col', 'state_col']:
            if kwds[i] is None:
                raise Exception("Please define --%s" % (i))
        sys.exit(args.func(syn=syn, **kwds))
