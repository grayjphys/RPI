#!/bin/bash
for d in POL-*/; do cd $d; rm CALCU*; cp ../CALCULATE_RAMAN.py .; ~/anaconda3/bin/python CALCULATE_RAMAN.py; cd ..; done; ~/anaconda3/bin/python TOTAL_RAMAN.py
