#!/bin/bash
grep -A 14 "ATOMIC FORCES in" Vibration-Pyridine.out | tail -14
grep -A 12 "i =" Pyridine-pos-1.xyz | tail -12
