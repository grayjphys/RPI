#!/bin/bash

#python3 update_com_input.py;
#  sed -i 's/      &CENTER_COORDINATES//g' PYRIDINE-UNPERTURBED/GEO-OPT.inp;
#  sed -i 's/      &END CENTER_COORDINATES//g' PYRIDINE-UNPERTURBED/GEO-OPT.inp;

python3 Displace-Along-Modes.py;

for d in */;do
  cd $d;
  cp ../PYRIDINE-UNPERTURBED/GEO-OPT.inp .;
  rm Vibration-Pyridine.out;
  rm *.cube;
  sed -i 's/INTENSITY 1.0E-2/INTENSITY 4.0E-2/g' GEO-OPT.inp
  sbatch geoopt-submit;
  cd ..;
done
