#!/bin/bash

for d in */; do cd $d; sbatch geoopt-submit; cd ..; done

cd ../RTTDDFT-Y/

for d in */; do cd $d; sbatch geoopt-submit; cd ..; done
