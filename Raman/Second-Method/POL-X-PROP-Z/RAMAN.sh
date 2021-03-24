#!/bin/bash

topdir="."
for d in "$topdir"/PYRIDINE-[MCU][OHN]*;
do
  if [ -d "$d" ]; then
     echo "$d"
     cd $d
     cp ../prepare_data.sh .
     cp ../start-matrix-files.sh .
     cp ../Get-CORE-Matrix.py .
     cp ../Get-Density-Matrix.py .
     ./prepare_data.sh
     cd ..
     for dd in "$d"/*;
     do
       if [ -d "$dd" ]; then
          echo "$dd"
          cd $dd
          cp ../prepare_data.sh .
          cp ../start-matrix-files.sh .
          cp ../Get-CORE-Matrix.py .
          cp ../Get-Density-Matrix.py .
          ./prepare_data.sh
          cd ../..
       fi
     done
  fi
done

python3 First-Derivative-External-Potential.py
python3 Second-Derivative-Density.py
python3 calc-intensity.py
