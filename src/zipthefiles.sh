#!/bin/bash

for dir in ~/PycharmProjects/Benchmark-Datasets/data/*; do
    if [ -d "$dir" ]; then
        dir_name=$(basename "$dir")
        if [ -n "$(find "$dir/${dir_name}_gml/data" -maxdepth 1 -type f -print -quit)" ]; then
            tar -czf "$dir/${dir_name}_gml.tar.gz" -C "$dir/${dir_name}_gml" .
            mv "$dir/${dir_name}_gml.tar.gz" "$dir/../../zipfiles_gml/${dir_name}_gml.tar.gz"
        fi
    fi
done