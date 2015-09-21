#!/bin/bash

BASE=$(dirname "$(cd $(dirname $0) && pwd)")
SUDO="sudo" #set to 'sudo' if need to sudo to access docker
SCRIPT="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"

$SUDO docker run --rm -v $BASE:$BASE -u `id -u` -w $BASE pancancer_upload_download /bin/bash $*
