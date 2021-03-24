#!/bin/bash
for d in POL-*/; do cd $d; rm CALCU*; cp ../CALCULATE_RAMAN_AT_COM.py .; python3 CALCULATE_RAMAN_AT_COM.py; cd ..; done; python3 TOTAL_RAMAN.py
