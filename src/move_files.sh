#!/bin/bash

for dir in ~/PycharmProjects/Benchmark-Datasets/data/*; do
    if [ -d "$dir" ]; then
        dir_name=$(basename "$dir")
        cp -a "$dir/clean_nx_json/." "$dir/$dir_name/data"
#        touch "$dir/$dir_name/README.md"
#        echo "Downloaded from [Graph Layout Benchmark Datasets Repository on OSF](https://osf.io/j7ucv/)" >> "$dir/$dir_name/README.md"
#        echo "" >> "$dir/$dir_name/README.md"
#        if [ -f "$dir/clean/metadata.txt" ]; then
#            cat "$dir/clean/metadata.txt" >> "$dir/$dir_name/README.md"
#        fi
    fi
done