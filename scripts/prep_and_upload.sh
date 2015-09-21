#!/bin/bash

DIR=$1

if [ ! -e $DIR/prep.sh.ready ]; then
    bash $DIR/prep.sh
fi

if [ ! -e $DIR/upload.sh.submitted ]; then 
    bash $DIR/upload.sh
fi
