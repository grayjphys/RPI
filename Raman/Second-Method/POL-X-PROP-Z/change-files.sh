#!/bin/bash

# 1 is number, 2 is PLUS, 3 is plus, 4 is Plus
cp ../GEO-OPT.inp  .
cp ../geoopt-submit .
mv GEO-OPT.inp  Pyridine-TDDFT-$4-Mode-$1.inp
sed -i 's/PYRIDINE-SAMPLE/PYRIDINE-TDDFT-'$2'-MODE-'$1'/g' Pyridine-TDDFT-$4-Mode-$1.inp
sed -i 's/pyridine.xyz/pyridine-'$3'-mode-'$1'.xyz/g' Pyridine-TDDFT-$4-Mode-$1.inp
sed -i 's/Pyridine-Sample-TDDFT.out/Pyridine-TDDFT-'$4'-Mode-'$1'.out/g' geoopt-submit
sed -i 's/Pyridine-Sample-TDDFT.inp/Pyridine-TDDFT-'$4'-Mode-'$1'.inp/g' geoopt-submit
sed -i 's/tddft/'$1'-Pyridine-TDDFT-'$4'-Mode/g' geoopt-submit
