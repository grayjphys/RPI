This is the newest iteration of the code to calculate the Raman spectrum of a molecule using dft results from CP2K.
This project looks at the spectrum for the Pyridine molecule.

POL-X-PROP-Y means that incident light is polarized in the X direction and propagates in the Y direction. All combinations of X,Y,Z 
are required for the spectrum except instances where the polarization and propagation are in the same direction. This means there should be a total of 
6 directories.

For each directory above, there should be two directories (PLUS, MINUS) for each vibrational mode of the molecule. Included in the POL-X-PROP-Y directory is the 
PYRIDINE-MINUS-MODE-01 directory. This means that the vibrational mode is subtracted (MINUS) from the positions of the atoms in the lowest energy structure of 
the molecule. There is both a minus and a plus directory for each mode because for each mode, derivatives need to be calculated and I chose the central difference 
derivative for higher accuracy. 

RUN-ALL.sh should be edited based on what you want to change across all files in a POL-[]-PROP-[] directory. It should then be copied to all of the POL-[]-PROP-[] 
directories and ran to submit to the edge cluster. To copy and run one simply needs to type: 

for d in */; do cp RUN-ALL.sh $d; cd $d; ./RUN-ALL.sh; cd ..; done

The geoopt-submit file is the submit script to be edited depending on your cluster specifications. This can be found in the PYRIDINE-[]-MODE-[] directories.

Once your jobs complete, run the ./CALC-RAMAN-AT-COM.sh command. This goes through and takes the electric fields and dipole moments and fourier transforms 
them into frequency space, then reads in the vibrational modes (To be calculated with the VIBRATIONAL_ANALYSIS routine in CP2K), and returns the Raman intensities
for each mode and each pol-[]-prop-[], then it averages over the pol-[]-prop-[] results to get a Raman spectrum.
