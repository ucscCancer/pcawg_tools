#!/usr/bin/env python

import re
import argparse
import requests
import nebula.docstore
from nebula.warpdrive import RemoteGalaxy
from nebula.target import Target

def report(ds, tags, all):
    if all or ds['visible']:
        if tags:
            print host, ds['name'], ",".join(ds.get('tags', []))
        else:
            print host, ds['name']

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("hostlist")
    parser.add_argument("--port", default="19999")
    parser.add_argument("--apikey", default="admin")
    parser.add_argument("-d", "--docstore", default=None)
    parser.add_argument("-t", "--tags", action="store_true", default=False)
    parser.add_argument("-a", "--all", action="store_true", default=False)
    parser.add_argument("-f", "--filter", default=None)

    parser.add_argument("cmd", nargs="*")
    args = parser.parse_args()

    hosts = []
    with open(args.hostlist) as handle:
        for line in handle:
            hosts.append(line.rstrip())

    docstore = None
    if args.docstore:
        docstore = nebula.docstore.from_url(args.docstore)

    if len(args.cmd):
        for host in hosts:
            rg = RemoteGalaxy("http://%s:%s" % (host, args.port), args.apikey)
            try:
                for hist in rg.history_list():
                    donor = None
                    out = []
                    if isinstance(hist, dict):
                        for ds in rg.get_history_contents(hist['id']):
                            visible = False
                            for t in ds.get('tags', []):
                                tmp = t.split(":")
                                if tmp[0] == 'donor':
                                    donor = tmp[1]
                                
                            if ds['visible'] or args.all:
                                visible = True
                            if visible and args.filter:
                                visible = False
                                if re.search(args.filter, ds['name']):
                                    visible = True
                            if visible:
                                if args.cmd[0] == 'error':
                                    if ds['state'] == 'error':
                                        out.append(ds)
                                elif args.cmd[0] == 'ok':
                                    if ds['state'] == 'ok':
                                        out.append(ds)
                                        if docstore:
                                            if ds['visible']:
                                                t = Target(uuid=ds['uuid'])
                                                if not docstore.exists(t):
                                                    prov = rg.get_provenance(ds['history_id'], ds['id'])
                                                    ds['provenance'] = prov
                                                    ds['id'] = ds['uuid']
                                                    ds['job'] = rg.get_job(prov['job_id'])

                                                    docstore.put(ds['uuid'], ds)
                                                    docstore.create(t)
                                                    path = docstore.get_filename(t)
                                                    print "Downoading %s" % (ds['download_url'])
                                                    rg.download(ds['download_url'], path)
                                                    docstore.update_from_file(t)
                                elif args.cmd[0] == 'running':
                                    if ds['state'] == 'running':
                                        out.append(ds)
                                elif args.cmd[0] == "all":
                                    out.append(ds)
                                else:
                                    print "unknown", ds['state']
                        for ds in out:
                            #report(ds, args.tags, args.all)
                            print host, donor, ds['name']
                    else:
                        print "host error", hist
            except requests.exceptions.ConnectionError:
                pass
