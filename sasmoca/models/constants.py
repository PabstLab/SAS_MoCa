############### GLOBAL VARIABLES ###############

""" Scattering model """
Nd = 12 # spacing of the dD_C array (chain thickness) for polydispersity integration
Nt = 45 # spacing of the t-array (cos(gamma)) for orintation integration

###########################################################################

""" Water parameters """
V_PW	= 24.5          # (AA^3) Perkins 2001 (Hydration shell in proteins)
b_HW	= 2.8179E-04    # (AA)
d_shl	= 3.1           # (AA) Perkins 2001 (Hydration shell in proteins) 

""" Buffer molecules """

# TRIS buffer
b_tris	= 1.860E-03 # scattering length (AA)
# EDTA
b_EDTA	= 4.340E-03 # scattering length (AA)

###########################################################################

""" ############################################# Lipids """

""" X-ray scattering length of acyl chain groups (AA) """
b_CH    = 1.97256E-04 
b_CH2   = 2.25435E-04 
b_CH3   = 2.53615E-04 

""" Phosphatidylcholine (PC) """

""" X-ray scattering length of head-related quasimolecular groups (AA)"""
b_PC	= 2.73340E-03
b_CG	= 1.88802E-03
b_PCN	= 1.97256E-03
b_Chol	= 7.60844E-04

""" Lipid-volume temperature-dependencies a0 + a1*T (AA^3)"""
""" POPC """
a0_V_POPC =	1228.1
a1_V_POPC =	0.935
""" DLPC """
a0_V_DLPC = 970.3
a1_V_DLPC = 0.670

""" Phosphatidylglycerol (PG) """

""" Lipid-head volume (AA^3) """
V_H_POPG = 289   
V_H_DLPG = 289   

""" X-ray scattering length of head-related quasimolecular groups (AA)"""
b_PG	= 2.47979E-03 
b_PG1	= 1.32443E-03 
b_PG2	= 1.15536E-03 

""" Lipid-volume temperature-dependencies a0 + a1*T (AA^3)"""
""" POPG """
a0_V_POPG =	1178.8
a1_V_POPG =	1.08
""" DLPG """
a0_V_DLPG = 934.1
a1_V_DLPG = 0.60

""" ############################################# OmpLA """

""" OmpLA parameters """
rho_ompla	= 1.202E-05 ## (AA^-2)
H_ompla 	= 44.0 # (AA) estimation obtained via PyMol 

""" Equivalent radius and volume of OmpLA"""
R_ompla_monomer = 17.2  # (AA) from PyMol  primary structure (1qd5)
V_ompla_monomer = 38906 # (AA^3) primary structure (1qd5)

R_ompla_dimer   = 20.3  # (AA) from PyMol  primary structure (1qd6)
V_ompla_dimer   = 71932 # (AA^3) primary structure (1qd6)