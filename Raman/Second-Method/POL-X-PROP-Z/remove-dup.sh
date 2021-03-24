#!/bin/bash

sed -i '/^      DELTA_PULSE_DIRECTION 0 0 1/{N;/\n      DELTA_PULSE_DIRECTION 0 0 1\s*$/D}' GEO-OPT.inp 
sed -i '/^      DELTA_PULSE_DIRECTION 0 0 1/{N;/\n     DELTA_PULSE_DIRECTION 0 0 1\s*$/D}' GEO-OPT.inp 

