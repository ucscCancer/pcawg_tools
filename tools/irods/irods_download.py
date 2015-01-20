#!/usr/bin/env python

import os
import argparse
import yaml
import json
import subprocess
from irods.session import iRODSSession
from irods.models import Collection, DataObject, DataObjectMeta
from irods.exception import CollectionDoesNotExist


def which(cmd):
    cmd = ["which",cmd]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    res = p.stdout.readline().rstrip()
    if len(res) == 0: return None
    return res

def get_env(config):
    env = {
        'irodsHost' : config['irods_host'],
        'irodsPort' : str(config['irods_port']),
        'irodsUserName': config['irods_user'],
        'irodsZone' : config['irods_zone'],
        'irodsAuthFileName' : os.path.join(os.getcwd(), "AUTH")
    }
    return env

def call_iinit(config):
    iinit = which("iinit")
    FNULL = open(os.devnull, 'w')
    subprocess.check_call("echo %s | %s" % (config['irods_pass'], iinit), env=get_env(config), shell=True, stdout=FNULL, stderr=FNULL)
    FNULL.close()
    #proc = subprocess.Popen(iinit, env=get_env(config), stdin=subprocess.PIPE)
    #proc.communicate(config['irods_pass'] + "\n")
    

def call_iget(config, src, dst):
    iget = which("iget")
    subprocess.check_call([iget, "-f", src, dst], env=get_env(config))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    parser.add_argument("--query")

    parser.add_argument("config")
    parser.add_argument("outfile")
    parser.add_argument("--unique", action="store_true", default=False)
    parser.add_argument("--test", action="store_true", default=False)
    parser.add_argument("--dataset_id", default=None)

    args = parser.parse_args()
    with open(args.config) as handle:
        txt = handle.read()
        config = yaml.load(txt)

    if args.path is not None:
        call_iinit(config)
        call_iget(config, args.path, args.outfile)
    else:
        sess = iRODSSession(host=config['irods_host'], port=config['irods_port'], user=config['irods_user'], password=config['irods_pass'], zone=config['irods_zone'])

        q = sess.query(Collection, DataObject)

        with open(args.query) as handle:
            query = yaml.load(handle.read())

        elm_map = {
            'meta' : {
                'name' : DataObjectMeta.name,
                'value' : DataObjectMeta.value
            },
            'collection' : {
                'value': Collection.name
            },
            'name' : {
                'value' : DataObject.name
            }
        }

        for filt in query:
            filter_args = []
            for field_name, field in elm_map[filt['method']].items():
                p = field == filt[field_name]
                filter_args.append(p)
            q = q.filter(*filter_args)

        #print q.criteria
        #print q.first() #this is broken right now
        val = None
        for row in q.all():
            if val is None:
                val = row
            else:
                if args.unique:
                    path1 = os.path.join(val[Collection.name], val[DataObject.name])
                    path2 = os.path.join(row[Collection.name], row[DataObject.name])
                    raise Exception("Non-unique file: %s and %s" % (path1, path2))

        if val is not None:
            path = os.path.join(val[Collection.name], val[DataObject.name])
            if not args.test:
                call_iinit(config)
                call_iget(config, path, args.outfile)
                if args.dataset_id is not None:
                    json_file = open( 'galaxy.json', 'w' )
                    info = dict(
                        type = 'dataset',
                        name = os.path.basename(path),
                        dataset_id = args.dataset_id
                    )
                    json_file.write( json.dumps( info ) )
                    json_file.close()
            else:
                print "Found", path

        else:
            raise Exception("File not found")
