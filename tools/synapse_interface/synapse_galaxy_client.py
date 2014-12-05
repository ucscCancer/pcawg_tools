#!/usr/bin/env python

import synapseclient
import sys
import json
import os
import shutil
import csv

def dict_list_to_tsv(data, handle):
    headers = {}
    for row in data:
        for col in row:
            if col not in headers:
                headers[col] = len(headers)

    head = headers.keys()
    head.sort(key=lambda x : headers[x])

    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(head)

    for row in data:
        out = []
        for c in head:
            out.append(row.get(c, ""))
        writer.writerow(out)


if __name__ == "__main__":
    cmd = sys.argv[1]
    info_file = sys.argv[2]

    handle = open(info_file)
    mode = handle.readline().rstrip("\n\r")
    syn = None
    if mode == 'password':
        username = handle.readline().rstrip("\n\r")
        password = handle.readline().rstrip("\n\r")
        syn = synapseclient.Synapse()
        syn.login(username, password, silent=True)
    elif mode == 'apikey':
        username = handle.readline().rstrip("\n\r")
        apikey = handle.readline().rstrip("\n\r")
        syn = synapseclient.Synapse()
        syn.login(username, apiKey=apikey, silent=True)
    elif mode == 'token':
        token = handle.readline().rstrip("\n\r")
        syn = synapseclient.Synapse()
        syn.login(sessionToken=token, silent=True)
    handle.close()

    if syn is None:
        sys.stderr.write("No login info\n")
        sys.exit(1)

    if cmd == "query":
        querypath = sys.argv[3]
        outpath = sys.argv[4]
        handle = open(querypath)
        query = handle.read()
        handle.close()
        print "Query:", query
        ohandle = open(outpath, "w")
        out = syn.query(query)
        if 'results' in out:
            dict_list_to_tsv(out['results'], ohandle)
        else:
            sys.stderr.write("Bad Result %s\n" % (out))
            sys.exit(1)
        ohandle.close()


    if cmd == "get":
        synid = sys.argv[3]
        outpath = sys.argv[4]
        dataset_id = int(sys.argv[5])

        ent = syn.downloadEntity(synid)
        src_path = os.path.join(ent['cacheDir'], ent['files'][0])
        shutil.copy(src_path, outpath)
        import uuid
        import binascii
        new_uuid = uuid.uuid4()
        if 'md5' in ent:
            new_uuid = uuid.UUID(bytes=binascii.unhexlify(ent['md5'])[:16], version=3)
            print "Setting uuid", new_uuid
        elif 'etag' in ent:
            new_uuid = uuid(ent['etag'])

        handle = open("galaxy.json", "w")
        handle.write(json.dumps({
            'type' : 'dataset',
            'dataset_id' : dataset_id,
            'ext' : 'auto',
            'stdout' : '',
            'name' : ent.name,
            'uuid' : str(new_uuid)
        }) + "\n")
        handle.close()


    if cmd == "upload":
        etype = sys.argv[3]
        properties_file = sys.argv[4]
        annotations_file = sys.argv[5]
        outfile = sys.argv[6]
        attach_file = sys.argv[7]
        attach_file_name = os.path.basename(sys.argv[8])
        synid = sys.argv[9]

        handle = open(properties_file)
        props = {}
        for line in handle:
            tmp = line.rstrip("\r\n").split("\t")
            props[tmp[0]] = tmp[1]
        handle.close()

        handle = open(annotations_file)
        annon = {}
        for line in handle:
            tmp = line.rstrip("\r\n").split("\t")
            if tmp[0] in annon:
                annon[tmp[0]].append( tmp[1] )
            else:
                annon[tmp[0]] = [tmp[1]]
        handle.close()


        if synid is not None and synid.startswith("syn"):
            entity = syn.getEntity(synid)
        else:
            entityData = { u'entityType': u'org.sagebionetworks.repo.model.Data' }
            if etype == 'folder':
                entityData = { u'entityType': u'org.sagebionetworks.repo.model.Folder' }
            entityData['name'] = attach_file_name
            entityData['parentId'] = props['parentId']
            entity = syn.createEntity(entityData)
            synid = entity['id']

        if len(props):
            for p in props:
                if p != 'parentId':
                    entity[p] = props[p]
            entity = syn.updateEntity(entity)

        if len(annon):
            ann = syn.getAnnotations(entity)
            for a in annon:
                ann[a] = annon[a]
            syn.setAnnotations(entity, ann)

        if attach_file != '-' and attach_file_name != '-':
            if not os.path.exists(attach_file_name):
                os.symlink(attach_file, attach_file_name)
                attach_file = attach_file_name
            print "upload", attach_file
            entity = syn.getEntity(synid)
            syn.uploadFile(entity, attach_file)

        handle = open(outfile, "w")
        handle.write(json.dumps(dict(entity)))
        handle.close()
