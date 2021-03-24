#!/bin/bash

topdir="."
for d in "$topdir"/PYRIDINE-[UC][NH]*;
do
  if [ -d "$d" ]; then
     echo "$d"
     cd $d
     rm *.out
     rm *wfn*
     rm *.cube
     rm *.pkl
     rm *MATRIX
     rm MOLog
     rm *.MOLog
     rm *.csv
     cp ../change-lines-print.py .
#     cp ../GEO-OPT.inp .
#     python change-lines-print.py
     sbatch geoopt-submit
     #./prepare_data.sh
     cd ..
     for dd in "$d"/*;
     do
       if [ -d "$dd" ]; then
          echo "$dd"
          cd $dd
          rm *.out
          rm *wfn*
          rm *.cube
          rm *.pkl
          rm *MATRIX
          rm MOLog
          rm *.MOLog
          rm *.csv
          cp ../change-lines-print.py .
#          cp ../GEO-OPT.inp .
#          python change-lines-print.py
          sbatch geoopt-submit
          #./prepare_data.sh
          cd ../..
       fi
     done
  fi
done

