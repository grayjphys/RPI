#!/bin/bash

grep -A 3250 "MO EIGEN" *.MOLog | tail -3249 > MOLog
sed -i 's/PYRIDINE-TDDFT-cartesian-mos-1_0.MOLog-//g' MOLog
