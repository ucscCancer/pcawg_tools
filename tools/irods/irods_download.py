#!/usr/bin/env python


import argparse
import yaml
from irods.session import iRODSSession
from irods.models import Collection, DataObject
from irods.exception import CollectionDoesNotExist



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    parser.add_argument("config")
    parser.add_argument("outfile")
    parser.add_argument("--dataset_id", default=None)

    args = parser.parse_args()
    with open(args.config) as handle:
        txt = handle.read()
        config = yaml.load(txt)

    sess = iRODSSession(host=config['irods_host'], port=config['irods_port'], user=config['irods_user'], password=config['irods_pass'], zone=config['irods_zone'])

    obj = sess.data_objects.get(args.path)

    with open(args.outfile, "wb") as out:
        with obj.open('r') as f:
            for block in f.read(2048):
                if block:
                    out.write(block)
