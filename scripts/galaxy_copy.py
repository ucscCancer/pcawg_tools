#!/usr/bin/env python

import re
import argparse
import requests
import nebula.docstore

from nebula.warpdrive import RemoteGalaxy
from nebula.target import Target

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("docstore")
    parser.add_argument("host")
    parser.add_argument("--apikey", default="admin")
    parser.add_argument("-f", "--filter", default=None)
    parser.add_argument("-l", "--list", action="store_true", default=False)
    

    args = parser.parse_args()

    docstore = nebula.docstore.from_url(args.docstore)

    rg = RemoteGalaxy(args.host, args.apikey)
    for hist in rg.history_list():
        for ds in rg.get_history_contents(hist['id']):
            copy = True           
            if not ds['visible']:
                copy = False
            if ds['name'] == "Unnamed dataset": #these datasets are inputs to the pipeline (misnamed due to bug)
                copy = False
            if ds['state'] != 'ok' and not args.list:
                copy = False
            if args.filter is not None and ds['name'] not in args.filter:
                copy = False
            if copy:
                if args.list:
                    print "Found", ds['name'], ds['state']
                else:
                    print "Copying", ds['name']
                    meta = rg.get_dataset(ds['id'], 'hda')
                    prov = rg.get_provenance(meta['history_id'], ds['id'])
                    meta['provenance'] = prov
                    meta['job'] = rg.get_job(prov['job_id'])
                    meta['id'] = meta['uuid'] #use the glocal id

                    hda = Target(uuid=meta['id'])
                    docstore.put(meta['id'], meta)

                    docstore.create(hda)
                    path = docstore.get_filename(hda)
                    rg.download(meta['download_url'], path)
                    docstore.update_from_file(hda)

