#!/bin/bash

for d in *; do
    if [ -d "$d" ]; then
        cd $d;
        echo $d;
        for dd in *; do
            if [ -d "$dd" ]; then
                cd $dd;
                cp ../*.py .;
                cp ../*.sh .;
                ./prepare_data.sh
                cd ..;
            fi
        done
        cd ..;
    fi
done
