#!/bin/bash

topdir="."
for p in "$topdir"/POL-*;
do
  if [ -d "$p" ]; then
     echo "$p"
     cd $p
     for d in PYRIDINE-[CMU][HON]*;
     do
       if [ -d "$d" ]; then
          echo "$d"
          cd $d
# grep 'DELTA\|INTENSITY\|POLARISATION' GEO-OPT.inp
sbatch geoopt-submit
          cd ..
          for dd in "$d"/*;
          do
            if [ -d "$dd" ]; then
               echo "$dd"
               cd $dd
# grep 'DELTA\|INTENSITY\|POLARISATION' GEO-OPT.inp
sbatch geoopt-submit
               cd ../..
            fi
          done
       fi
     done
     cd ..
  fi
done
