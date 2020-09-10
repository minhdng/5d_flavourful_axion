# 5d_flavourful_axion

## Package

'WarpFlavourAxion' includes all constants and functions  

## Notebooks

'yukawa_generator.nb' randomly generates pairs of Y_u, Y_d (Y_n, Y_e) then select through 3 constraints
1. Whether the pairs can generate the right rho and eta parameters of the CKM (PMNS) matrix 
2. Whether the pairs result in a product of the fermion effective profile less than 1 (such that there exists corresponding profile parameters c_L, c_R)
3. Whether the approximate A matrices (from the SVD decomposition of 4d Yukawa matrices) are unitary

The raw data, i.e. surviving pairs after each constraints, are saved in 'data/'. The final result is saved in 'output/'


