#!/usr/bin/env python

import os
import argparse
import yaml
import subprocess
from irods.session import iRODSSession
from irods.models import Collection, DataObject
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
    parser.add_argument("--dataset_id", default=None)

    args = parser.parse_args()
    with open(args.config) as handle:
        txt = handle.read()
        config = yaml.load(txt)

    call_iinit(config)
    if args.path is not None:
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

        for k, v in query.items():
            args = []
            for field_name, field in elm_map[k].items():
                p = field == v[field_name]
                args.append(p)
            q = q.filter(*args)

        val = q.first()
        if val is not None:
            path = os.path.join(val[Collection.name], val[DataObject.name])
            call_iget(config, path, args.outfile)
        else:
            raise Exception("File not found")
