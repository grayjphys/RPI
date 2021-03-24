#!/bin/bash

for i in {1..27}; do cd PYRIDINE-MINUS-MODE-$i; sed -i "s/mode-1/mode-{$i}/g" GEO-OPT.inp; sed -i 's/{//g' GEO-OPT.inp; sed -i 's/}//g' GEO-OPT.inp; cd ..; done
for i in {1..27}; do cd PYRIDINE-PLUS-MODE-$i; sed -i "s/minus-mode-\$i/plus-mode-{$i}/g" GEO-OPT.inp; sed -i 's/{//g' GEO-OPT.inp; sed -i 's/}//g' GEO-OPT.inp; cd ..; done
for i in {1..27}; do cd PYRIDINE-MINUS-MODE-$i; sed -i "s/MODE-1/MODE-{$i}/g" GEO-OPT.inp; sed -i 's/{//g' GEO-OPT.inp; sed -i 's/}//g' GEO-OPT.inp; cd ..; done
for i in {1..27}; do cd PYRIDINE-PLUS-MODE-$i; sed -i "s/MINUS-MODE-1/PLUS-MODE-{$i}/g" GEO-OPT.inp; sed -i 's/{//g' GEO-OPT.inp; sed -i 's/}//g' GEO-OPT.inp; cd ..; done

