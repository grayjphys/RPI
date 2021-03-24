#!/bin/bash

for d in */; do cd $d; rm *.cube; rm *.restart*; rm *wfn*; rm *bak*; cd ..; done
