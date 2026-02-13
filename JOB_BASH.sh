#!/bin/bash

target_dir="../Angel/files_python/bash/"
current_dir="$(pwd)"

files=(
    hmi/hmi_disambig.py
    hmi/limbcorrect.py
    map_cut/map_cut.py
    map_projection/los2los.py
    map_projection/los2postel.py
    map_projection/postel2los.py
    coordinate_conversion/hpxy2lonlat.py
    coordinate_conversion/hpxy2xy.py
    coordinate_conversion/lonlat2hpxy.py
    coordinate_conversion/lonlat2xy.py
    coordinate_conversion/xy2hpxy.py
    coordinate_conversion/xy2lonlat.py
)

for file in "${files[@]}"
do
    name="${file##*/}"
    # rm -rf "$name"
    # rm -rf $dir/$name
    ln -sF ${current_dir}/$file ${target_dir}/$name
done
