#!/usr/bin/python

import numpy as np
import pandas as pd
from scipy.special import erf

from models.shared import (mu4,
                           PDF_normal,
                           FTreal_erf, FTreal_gauss,
                           Slab, Gauss,
                           Water_volume)

from models.constants import *

 
###################################################################
###################################################################
###############		DMPC/DMPG 95/5 mol/mol LUVs		############### 
###############		Suspension in water 			############### 
###################################################################
###################################################################
 
#	##############	DMPC and DMPG	############### 
#	1,2-dimyristoyl-sn-glycero-3-phosphocholine 
#	1,2-dimyristoyl-sn-glycero-3-phosphoglycerol

############### GLOBAL VARIABLES ###############

# POPG molar ratio
x_PG  = 0.05

# number of chain groups
n_CH  = 0
n_CH2 = 24
n_CH3 = 2  		

################################################################################################################## 
################################################################################################################## 

#########################################################
def Lipid_volume(T) :
	"""Return the average lipid volume including a molar fraction of DMPG"""
	return (1-x_PG) * (a0_V_DMPC +T * a1_V_DMPC) + x_PG * (a0_V_DMPG + T * a1_V_DMPG) 

#########################################################
def Volumes_DMPC(V_H_DMPC, T, r32, r_CG, r_PCN) :	
	"""Return volumes of each quasi-molecular group"""

	# Volumes 
	V_L = Lipid_volume(T)
	V_HC = V_L - ( (1-x_PG) * V_H_DMPC + x_PG * V_H_PG )

	# Quasi-molecular volumes 
	V_CH2	= V_HC / ( n_CH2 + n_CH3*r32 )	# Volume of CH2 groups  
	V_CH3	= V_CH2	*	r32 							# Volume of CH3 groups

	V_CG	= V_H_DMPC	*	r_CG				# Volume of CG group 
	V_PCN	= V_H_DMPC	* 	r_PCN 				# Volume of PCN group 
	V_Chol	= V_H_DMPC	*   (1 - r_PCN - r_CG)	# Volume of CholCH3 group  

	V_PG1	= V_H_PG	*	0.16 				# Pan et al. BBA 2012
	V_PG2	= V_H_PG	*	( 1 - 0.51 - 0.16) 	# Pan et al. BBA 2012
	
	return V_L, V_HC, V_CH2, V_CH3, V_CG, V_PCN, V_Chol, V_PG1, V_PG2

#########################################################
def SLDs_DMPC(T, V_Chol, V_PCN, V_CG, V_PG2, V_PG1, V_CH2, V_CH3, V_BW) :
	"""Return SLDs of solvent and quasi-molecular lipi groups (AA^-2)"""
	
	############### X-ray scattering length densities (nm^-2)
	rho_sol		= b_HW / Water_volume(T) 
	drho_Chol	= ( (1-x_PG)*b_Chol/V_Chol	+ x_PG*b_PG2/V_PG2 )	- rho_sol 
	drho_PCN	= ( (1-x_PG)*b_PCN/V_PCN	+ x_PG*b_P/V_PG1 )		- rho_sol 
	drho_CG		= b_CG  / V_CG  									- rho_sol 
	drho_CH2	= b_CH2 / V_CH2  									- rho_sol 
	drho_CH3	= b_CH3 / V_CH3  									- rho_sol 
	drho_BW		= b_HW	/ V_BW										- rho_sol 	
	
	return rho_sol, drho_Chol, drho_PCN, drho_CG, drho_CH2, drho_CH3, drho_BW

#########################################		VESICLE SCATTERING MODEL	######################################### 

#########################################		SDP MODELLING			######################################### 
####################################		SEPARATED FORM FACTOR			#################################### 
##################################################################################################################

class LUV_DMPC:

##################
	def __init__(self, q, PAR, _) :
		self.q = q
		[self.Norm, self.nv, 
		self.Rm, self.Z, 
		self.d_Chol, self.s_Chol, self.d_PCN, self.s_PCN, self.d_CG, self.s_CG, 
		self.dD_C, self.r_D_C, 
		self.s_CH2, self.s_CH3, 
		self.r_PCN, self.r_CG, self.r32, 
		self.T, self.V_H, self.V_BW, 
		self.Con] = PAR

		### polydispersity value 		
		alp = self.Rm/(self.Z+1)

		### Volumes 
		[_, 
		V_HC, V_CH2, self.V_CH3, 
		self.V_CG, self.V_PCN, self.V_Chol, 
		self.V_PG1, self.V_PG2] = Volumes_DMPC(self.V_H, self.T, self.r32, self.r_CG, self.r_PCN)
			
		### Reduce to single lipid thickness
		self.D_C = 0.5*self.dD_C
		### Get mean area per lipid (A_L)
		self.A_L = V_HC / self.D_C
		
		### X-ray scattering length densities (nm^-2) 	
		[_, self.drho_Chol, self.drho_PCN, self.drho_CG,
   	 	self.drho_CH2, self.drho_CH3, 
	 	self.drho_BW] = SLDs_DMPC(	self.T,
							 		self.V_Chol, self.V_PCN, self.V_CG, 
	 								self.V_PG2, self.V_PG1, 
	 								V_CH2, self.V_CH3, self.V_BW )
			
		### D_C polydispersity
		self.s_D_C = self.r_D_C*self.D_C
		HC_array = np.linspace(self.D_C-3*self.s_D_C, self.D_C+3*self.s_D_C, Nd)
		Normal = PDF_normal(HC_array, self.D_C, self.s_D_C)

		### c-prefactors
		c_CH3 = np.zeros(Nd, dtype=float)

		c_Chol	= ( (1-x_PG)*self.V_Chol	+ x_PG*self.V_PG2 )	/ self.A_L 
		c_PCN	= ( (1-x_PG)*self.V_PCN		+ x_PG*self.V_PG1 )	/ self.A_L
		c_CG	= self.V_CG / self.A_L 

		c_CH3	= self.V_CH3	* n_CH3	/ (V_HC / HC_array ) 
	
		### Assembling membrane scattering amplitude		
		Am = np.zeros([self.q.size, Nd], dtype=float)

		self.d_BW = self.d_CG + self.d_PCN + self.d_Chol + d_shl	

		# Adding hydrocarbon-chain envelope 
		Am += 2 * self.drho_CH2 * HC_array[None,:] * FTreal_erf(self.q[:,None], 0, 2*HC_array[None,:], self.s_CH2)
		# Adding CH and CH3 groups 
		Am += 2 * (self.drho_CH3 - self.drho_CH2) * c_CH3[None,:] * FTreal_gauss(self.q[:,None], 0,         self.s_CH3) 
		# Adding hydration-water envelope 
		Am += 4 * self.drho_BW * self.d_BW/2. * FTreal_erf(self.q[:,None], (HC_array[None,:]+self.d_BW/2.), self.d_BW, self.s_CH2)
		# Adding CG, PCN and CholCH3 groups 
		Am += 2 * (self.drho_CG		- self.drho_BW) * c_CG	* FTreal_gauss(self.q[:,None], (HC_array[None,:] + self.d_CG),	self.s_CG) 
		Am += 2 * (self.drho_PCN	- self.drho_BW) * c_PCN	* FTreal_gauss(self.q[:,None], (HC_array[None,:] + self.d_CG + self.d_PCN), self.s_PCN) 
		Am += 2 * (self.drho_Chol	- self.drho_BW) * c_Chol* FTreal_gauss(self.q[:,None], (HC_array[None,:] + self.d_CG + self.d_PCN + self.d_Chol), self.s_Chol)

		# Calculating vesicle scattering intensity from OmpLA-vesicle interactions; HC-dependent
		self.Iv = Am[:,:]**2 * 16*np.pi**2*mu4(self.q[:,None],self.Z,alp)
		
		# Initializing final intensity array
		self.Ipoly = np.empty([self.q.size], dtype=float)	

		### Calculating average over D_C PDF		
		self.Iv[:,] *= Normal				
		self.Ipoly = np.trapz(self.Iv, HC_array, 6*self.s_D_C/(Nd-1), axis=1)

##################			
	def intensity(self):
		return ( self.Norm * self.nv*1e5 ) * self.Ipoly + self.Con*( 0.99*(1./(1+np.exp(-8*(self.q-0.1)))) + 0.01 )	
		
##################			
	def negative_water(self):

		# Calculate the Luzzati thickness
		D_B = 2*Lipid_volume(self.T)/self.A_L

		# Set z arrasy to plt SDP profile...
		z_array = np.linspace(0., D_B*1.1, int(D_B*1.1)*10, endpoint=False)

		self.check = False

		CG		= Gauss(z_array, self.V_CG,			self.D_C+self.d_CG,							self.s_CG,		self.A_L)
		PCN		= Gauss(z_array, self.V_PCN,		self.D_C+self.d_CG+self.d_PCN,				self.s_PCN,		self.A_L)
		Chol	= Gauss(z_array, self.V_Chol,		self.D_C+self.d_CG+self.d_PCN+self.d_Chol,	self.s_Chol,	self.A_L)
		BW		= Slab(z_array,	self.D_C+self.d_BW/2., self.d_BW, self.s_CH2) - CG - PCN - Chol

		if (BW<0).any(): self.check = True

		return self.check
		
##################			
	def SDP_profile(self):

		# Calculate the Luzzati thickness
		D_B = 2*Lipid_volume(self.T)/self.A_L

		# Set z arrasy to plt SDP profile...
		z_array = np.linspace(0., D_B*1.1, int(D_B*1.1)*10, endpoint=False)
		
		CH2		= Slab(z_array, 0, 2*self.D_C, self.s_CH2)
		BW		= Slab(z_array,	self.D_C+self.d_BW/2., self.d_BW, self.s_CH2) 

		CH3		= Gauss(z_array, 2*n_CH3*self.V_CH3, 0, self.s_CH3, self.A_L)
		CH2		= CH2 - CH3

		CG		= Gauss(z_array, self.V_CG,			self.D_C+self.d_CG,							self.s_CG,		self.A_L)
		PCN		= Gauss(z_array, self.V_PCN,		self.D_C+self.d_CG+self.d_PCN,				self.s_PCN,		self.A_L)
		Chol	= Gauss(z_array, self.V_Chol,		self.D_C+self.d_CG+self.d_PCN+self.d_Chol,	self.s_Chol,	self.A_L)
		BW		= BW - CG - PCN - Chol

		# ...and DeltaSLD profile
		SLD = np.zeros_like(z_array)
		SLD += CH3 * self.drho_CH3 + CH2 * self.drho_CH2
		SLD += CG * self.drho_CG + PCN * self.drho_PCN + Chol * self.drho_Chol
		SLD += BW * self.drho_BW

		# Calculate peak-to-peak distance D_pp (proxy for phosphate-to-phosphate distance)
		SLD_peak = np.argmax(SLD)
		D_pp = 2*z_array[SLD_peak]

		# Calculate the number of water molecules per lipid headgroup
		n_W =  ( self.A_L * self.d_BW*erf(self.d_BW/2./(np.sqrt(2)*self.s_CH2)) - self.V_H ) / self.V_BW
		

		SDP_matrix = pd.DataFrame(np.column_stack((z_array, CH3, CH2, CG, PCN, Chol, BW, SLD)), 
								  columns=["z", "CH3", "CH2", "CG", "PCN", "Chol", "BW", "SLD"])

		return 	SDP_matrix, n_W, D_pp, D_B, self.D_C, self.A_L
				
##################################################################################################################
##################################################################################################################
