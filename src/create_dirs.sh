#!/bin/bash

for dir in ~/PycharmProjects/Benchmark-Datasets/data/*; do
    if [ -d "$dir" ]; then
        dir_name=$(basename "$dir")
        mkdir "$dir/${dir_name%/}_gml"
        cp "$dir/$dir_name/README.md" "$dir/${dir_name%/}_gml"
        mkdir "$dir/${dir_name%/}_gml/data"
    fi
done