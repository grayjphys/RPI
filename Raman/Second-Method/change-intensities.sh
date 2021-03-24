cd POL-X-PROP-Z/;
cp ../../CP2K-Pyridine-1E11-Intensity/POLX-PROP-Z/DO-FOR-ALL-DIRS.sh .;
./DO-FOR-ALL-DIRS.sh;
sed -i 's/1E10/1E11/g' create-efield-distortion.py;
python3 create-efield-distortion.py;
cp ../../CP2K-Pyridine-1E11-Intensity/POLX-PROP-Z/RUN-ALL.sh .;
./RUN-ALL.sh;

cd ../POL-X-PROP-Y/;
cp ../POL-X-PROP-Z/DO-FOR-ALL-DIRS.sh .;
./DO-FOR-ALL-DIRS.sh;
sed -i 's/1E10/1E11/g' create-efield-distortion.py;
python3 create-efield-distortion.py;
cp ../POLX-PROP-Z/RUN-ALL.sh .;
./RUN-ALL.sh;

cd ../POL-Z-PROP-Y/;
cp ../POL-X-PROP-Z/DO-FOR-ALL-DIRS.sh .;
./DO-FOR-ALL-DIRS.sh;
sed -i 's/1E10/1E11/g' create-efield-distortion.py;
python3 create-efield-distortion.py;
cp ../POLX-PROP-Z/RUN-ALL.sh .;
./RUN-ALL.sh;

cd ../POL-Z-PROP-X/;
cp ../POL-X-PROP-Z/DO-FOR-ALL-DIRS.sh .;
./DO-FOR-ALL-DIRS.sh;
sed -i 's/1E10/1E11/g' create-efield-distortion.py;
python3 create-efield-distortion.py;
cp ../POLX-PROP-Z/RUN-ALL.sh .;
./RUN-ALL.sh;

cd ../POL-Y-PROP-X/;
cp ../POL-X-PROP-Z/DO-FOR-ALL-DIRS.sh .;
./DO-FOR-ALL-DIRS.sh;
sed -i 's/1E10/1E11/g' create-efield-distortion.py;
python3 create-efield-distortion.py;
cp ../POLX-PROP-Z/RUN-ALL.sh .;
./RUN-ALL.sh;

cd ../POL-Y-PROP-Z/;
cp ../POL-X-PROP-Z/DO-FOR-ALL-DIRS.sh .;
./DO-FOR-ALL-DIRS.sh;
sed -i 's/1E10/1E11/g' create-efield-distortion.py;
python3 create-efield-distortion.py;
cp ../POLX-PROP-Z/RUN-ALL.sh .;
./RUN-ALL.sh;
