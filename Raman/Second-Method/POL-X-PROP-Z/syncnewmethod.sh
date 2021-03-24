#!/bin/bash

rsync -rav -e ssh --exclude="*.cube" --exclude="*.sh" --exclude="*.py" --exclude="*.pkl"\
 --exclude="*MATRIX" --exclude="*MOLog" grayj6@edge.phys.rpi.edu:/home/grayj6/data/CP2K-Pyridine/VIBRATION-better/GEO-OPT/New-Method ~/Desktop/EDGE/CP2K-Pyridine/VIBRATION-better/GEO-OPT

