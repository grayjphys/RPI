#!/bin/bash

topdir="."
for d in "$topdir"/PYRIDINE-[CUM][HNO]*;
do
  if [ -d "$d" ]; then
     echo "$d"
     cd $d
#    cp ../pyridine.xyz .
     cp ../PYRIDINE-UNPERTURBED/GEO-OPT.inp .
     sed -i 's/1.0E10/1.0E11/g' GEO-OPT.inp
#      sbatch geoopt-submit
    sed -i 's/8/16/g' geoopt-submit
      cd ..
     for dd in "$d"/*;
     do
       if [ -d "$dd" ]; then
          echo "$dd"
          cd $dd
#      cp ../pyridine.xyz
     cp ../../PYRIDINE-UNPERTURBED/GEO-OPT.inp .
     sed -i 's/1.0E10/1.0E11/g' GEO-OPT.inp
#        sbatch geoopt-submit
sed -i 's/8/16/g' geoopt-submit
         cd ../..
       fi
     done
  fi
done
