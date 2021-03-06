# 5D Flavourful Axion

Quality 5D Axion flavour-changing couplings to quarks and leptons.

## Package

`WarpFlavourAxion` includes all constants and functions.  

Examples of function usages are included in `sample.nb`.

## Notebooks

### `pmns.nb`
calculates the Wolfenstein paramters of the PMNS matrix. 

### `yukawa_generator.nb` 
randomly generates pairs of Y<sub>u</sub>, Y<sub>d</sub> (Y<sub>n</sub>, Y<sub>e</sub>) then systematically eliminates through 3 constraints: 
1. Whether the pairs can generate the right rho and eta parameters of the CKM (PMNS) matrix.
2. Whether the pairs result in a product of the fermion effective profile less than 1 (such that there exists corresponding profile parameters c<sub>L</sub>, c<sub>R</sub>).
3. Whether the approximate A matrices (from the SVD decomposition of 4d Yukawa matrices) are unitary.

The raw data, i.e. surviving pairs after each constraints, are saved in `data/`. The final result is saved in `output/`.

### `aff_offdiag_couling.nb` 
calculates off-diagonal axion-fermion-fermion couplings and produce corresponding figures. 
