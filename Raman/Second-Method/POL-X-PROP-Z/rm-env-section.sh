#!/bin/bash

sed -i 's/&CONSTANT_ENV//g' GEO-OPT.inp;
sed -i 's/&END CONSTANT_ENV//g' GEO-OPT.inp;
sed -i 's/START_STEP 1//g' GEO-OPT.inp;
sed -i 's/END_STEP 2//g' GEO-OPT.inp;

