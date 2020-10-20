#!/bin/bash
grep -A 5980 "DENSITY MATRIX" Vibration-Pyridine.out | tail -5980 > D-MATRIX
#grep -A 5980 "KOHN-SHAM" Vibration-Pyridine.out | tail -5980 > KS-MATRIX
#grep -A 5980 "POTENTIAL ENERGY MATRIX" Vibration-Pyridine.out | tail -5980 > PE-MATRIX
#grep -A 5980 "Matrix VXC" Vibration-Pyridine.out | tail -5980 > VXC-MATRIX
#grep -A 5980 "KINETIC ENERGY MATRIX" Vibration-Pyridine.out | tail -5980 > KE-MATRIX
grep -A 5980 "CORE HAMILTONIAN MATRIX" Vibration-Pyridine.out | tail -5980 > CORE-MATRIX

