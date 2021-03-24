#!/bin/bash

#./clip-file.sh;
./start-matrix-files.sh;
python3 Get-Density-Matrix.py;
#python3 Get-Efield.py;
#python3 Get-Kohn-Sham-Matrix.py;
#python3 Get-States.py;
#python3 Get-VXC-Matrix.py;
#python3 Get-PE-Matrix.py;
#python3 Get-KE-Matrix.py;
python3 Get-CORE-Matrix.py;
#python3 check-data-for-error.py;
