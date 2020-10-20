#!/bin/bash

./start-matrix-files.sh;
python3 Get-Density-Matrix.py;
python3 Get-CORE-Matrix.py;
