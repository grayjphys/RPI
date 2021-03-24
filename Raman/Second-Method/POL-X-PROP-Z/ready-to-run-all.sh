#!/bin/bash

echo "did you change the sample input?"
nano Pyridine-Sample-TDDFT.inp 
echo "done editing sample input"

read -p "Press any key to resume ..."

cp ../PYRIDINE-VIBRATIONS-1.mol .
python shift-positions-along-normal-modes.py
./move-displaced-geometries.sh

cp Pyridine-Sample-TDDFT.inp ../RTTDDFT-POL-Y
for d in */; do cp Pyridine-Sample-TDDFT.inp $d; done
for d in */; do rm "$d"PYRIDINE-*; done
./run-change-files.sh
for d in */; do cd $d; sbatch geoopt-submit; cd ..; done

read -p "Press any key to resume ..."

cd ../RTTDDFT-POL-Y/

cp ../PYRIDINE-VIBRATIONS-1.mol .
python shift-positions-along-normal-modes.py
./move-displaced-geometries.sh

for d in */; do cp Pyridine-Sample-TDDFT.inp $d; done
for d in */; do rm "$d"PYRIDINE-*; done
./run-change-files.sh
for d in */; do cd $d; sbatch geoopt-submit; cd ..; done
