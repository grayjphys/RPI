#!/bin/bash

for d in */; do cp change-files.sh $d; done
for i in {1..27}; do cd PYRIDINE-PLUS-MODE-$i; pwd; ./change-files.sh $i PLUS plus Plus; cd ..; done;
for i in {1..27}; do cd PYRIDINE-MINUS-MODE-$i; pwd; ./change-files.sh $i MINUS minus Minus; cd ..; done;
