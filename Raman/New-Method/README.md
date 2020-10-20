#Simplified instructions for how to compute the Raman spectrum using CP2K.

To prepare data for DFT calculations:
1. Calculate the vibrational spectrum and the optimal geometry.
2. (a. b. c. in any order)
   a. For the optimal geometry you will need to calculate the Density Matrix and Core Hamiltonian (Kinetic plus Potential Energy) in the atomic orbital basis.
   b. i. Create folders where the geometry has an atom displaced a direction for every direction (x,y,z) and every atom (N atoms), with a small value. 
      Should have 2*3N folders. Remember to have positive and negative displacements. Should have (Atom 1 +x, Atom 1 -x, Atom 1, +y, Atom 1 -y, Atom 1 +z, 
      Atom 1 -z, Atom 2 +x, etc.).
      ii. Run displace-geometry.py
   c. i. Create folders where the electric field components are varied up to the second derivative. Should have x, y, z first derivatives with positive and negative       changes as well as xx, yy, zz, xy, xz, and yz mixed derivatives, with positive and negative changes. All: (+E_x, -E_x, +E_y, -E_y, +E_z, -E_z, 
      {+E_x then +E_x}, {+E_x then +E_y}, {+E_x then -E_y}, {+E_x then +E_z}, {+E_x then -E_z}, {-E_x then -E_x}, {-E_x then +E_y}, {-E_x then -E_y}, 
      {-E_x then +E_z}, {-E_x then -E_z}, {+E_y then +E_y}, {+E_y then +E_z}, {+E_y then -E_z}, {-E_y then -E_y}, {-E_y then +E_z}, {-E_y then -E_z},
      {+E_z then +E_z}, {-E_z then -E_z}).
      ii. Run create-efield-distortion.py
      
Note on naming convention: for the folders pertaining to displacements call your file MOLECULE-MOVE-ATOM-01-X-POS or MOLECULE-MOVE-ATOM-11-Z-NEG as examples. for the folders pertaining to single changes to electric field components call your file MOLECULE-CHANGE-X-FIELD-POSITIVE or MOLECULE-CHANGE-Y-FIELD-NEGATIVE as examples. For folders pertaining to double changes to electric field components, place the folder within the single change folder corresponding to the correct order of change, then call the folder MOLECULE-CHANGE-XX-FIELD-POS-POS or MOLECULE-CHANGE-YZ-NEG-POS as examples. MOLECULE can be named anything.

To calculate spectrum run ./RAMAN.sh which:
1. Copies the prepare-data.sh, start-matrix-files.sh, Get-Density-Matrix.py, and Get-CORE-Matrix.py scripts to each subdirectory, and runs it.
    a. prepare_data.sh creates the Density Matrix (D-MATRIX) and Core Hamiltonian Matrix (CORE-MATRIX) files by reading the output file with start-matrix-files.sh.
    b. prepare_data.sh then runs python scripts called Get-Density-Matrix.py and Get-CORE-Matrix.py to create pkl files of numpy arrays of the matrices so 
       later python scripts can quickly read in the matrices.
2. Runs First-Derivative-External-Potential.py which reads in the Core Hamiltonian Matrix and calculates the first derivative of the matrix with respect to changes
   to the atomic coordinates. 
3. Runs Second-Derivative-Density.py which reads in the Density Matrix and calculates the second derivative of the matrx with respect to changes to the components of
   the electric field.
4. Runs calc-intensity.py which takes the first derivative of the Core Hamiltonian Matrix and the second derivative of the Density Matrix along with the eigenvalues
   and eigenmodes of vibration for a molecule, then it calculates the Raman spectrum and plots the results.
   
The math behind the code can be found in the paper:
First-Principles Calculation of Vibrational Raman Spectra in Large Systems: Signature of Small Rings in Crystalline SiO2
Michele Lazzeri and Francesco Mauri
Phys. Rev. Lett. 90, 036401 – Published 23 January 2003



NOTE: YOU WILL HAVE TO CHANGE THE PYTHON SCRIPTS DEPENDING ON YOUR MOLECULE BECAUSE THE NAMES OF THE FOLEDERS WILL CHANGE. There will be sections in the code with d.name[27] or whatever and you will have to change the number/range to match the file name so that it finds the number of the atom, "P" or "N" for positive or negative, and "X", "Y", "Z" for displacement directions.
