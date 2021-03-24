#!/bin/bash

topdir="."
for d in "$topdir"/PYRIDINE-[MCU][OHN]*;
do
  if [ -d "$d" ]; then
     echo "$d"
     cd $d
     cp ../prepare_data.sh .
     ./prepare_data.sh
     cd ..
     for dd in "$d"/*;
     do
       if [ -d "$dd" ]; then
          echo "$dd"
          cd $dd
          cp ../prepare_data.sh .
          ./prepare_data.sh
          cd ../..
       fi
     done
  fi
done
