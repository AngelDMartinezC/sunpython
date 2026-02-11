#!/bin/bash

dir="../bash/"
cd $dir

files=(
    ../sunpython/hmi/hmi_disambig.py
    ../sunpython/map_cut/map_cut.py
    ../sunpython/coordinate_conversion/*2*
    ../sunpython/hmi/limbcorrect.py
    ../sunpython/map_projection/*2*
)

for file in "${files[@]}"
do
    name="${file##*/}"
    rm -rf "$name"
    ln -s "$file" "$name"
done
