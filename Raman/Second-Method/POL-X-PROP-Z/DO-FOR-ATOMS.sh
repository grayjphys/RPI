#!/bin/bash

topdir="."
for d in "$topdir"/PYRIDINE-[M][O]*;
do
  if [ -d "$d" ]; then
     echo "$d"
     cd $d
#     rm *.out
#     rm *wfn*
#     rm *.cube
#     rm *.pkl
#     rm *MATRIX
#     rm MOLog
#     rm *.MOLog
#     rm *.csv
#     cp ../change-lines-print.py .
#     cp ../GEO-OPT.inp .
#     python change-lines-print.py
     sbatch geoopt-submit
     #./prepare_data.sh
#     sed -i 's/SCF_GUESS ATOMIC/SCF_GUESS RESTART/g' GEO-OPT.inp
     cd ..
     for dd in "$d"/*;
     do
       if [ -d "$dd" ]; then
          echo "$dd"
          cd $dd
#     rm *wfn*
#     rm *.cube
#     rm *.pkl
#     rm *MATRIX
#     rm MOLog
#     rm *.MOLog
#     rm *.csv
#     cp ../change-lines-print.py .
#     cp ../GEO-OPT.inp .
#     python change-lines-print.py
     sbatch geoopt-submit
     #./prepare_data.sh
#     sed -i 's/SCF_GUESS ATOMIC/SCF_GUESS RESTART/g' GEO-OPT.inp
          cd ../..
       fi
     done
  fi
done
