#!/bin/bash

topdir="."
for d in "$topdir"/PYRIDINE-[MCU][OHN]*;
do
  if [ -d "$d" ]; then
     echo "$d"
     cd $d
     sbatch geoopt-submit
      cd ..
     for dd in "$d"/*;
     do
       if [ -d "$dd" ]; then
          echo "$dd"
          cd $dd
          sbatch geoopt-submit
         cd ../..
       fi
     done
  fi
done
