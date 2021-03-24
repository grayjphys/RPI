#!/bin/bash

grep -A 2990 KOHN-SHAM Vibration-Pyridine.out | tail -2990 > KS-MATRIX
grep -A 2990 "DENSITY MATRIX" Vibration-Pyridine.out | tail -2990 > D-MATRIX

