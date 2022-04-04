#!/bin/bash

# primitive bash script
# add csv suffix to files in directory
# based on https://unix.stackexchange.com/questions/19654/how-do-i-change-the-extension-of-multiple-files

DIR=$1
for f in "$DIR"/*; do # (!) not working
    echo $f
    mv -- "$f" "${f}.csv"
done