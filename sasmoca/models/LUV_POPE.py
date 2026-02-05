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
################################################################### 
###############		POPE/POPG 90/10 mol/mol LUVs	############### 
###############		Suspension in water 			############### 
################################################################### 
################################################################### 
################################################################### 

#	##############	POPE and POPG	############### 
#	1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine 
#	1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoglycerol 

############### GLOBAL VARIABLES ###############

# POPG molar ratio
x_PG  = 0.10

# number of chain groups
n_CH  = 2 
n_CH2 = 28
n_CH3 = 2  		

################################################################################################################## 
################################################################################################################## 


#########################################################
def Lipid_volume(T) :
	"""Return the average lipid volume including a molar fraction of POPG"""
	return (1-x_PG) * (a0_V_POPE + T * a1_V_POPE) + x_PG * (a0_V_POPG + T * a1_V_POPG) 

#########################################################
def Volumes_POPE(V_H_POPE, T, r12, r32, r_CG, r_P) :	
	"""Return volumes of each quasi-molecular group"""

	# Volumes 
	V_L = Lipid_volume(T)
	V_HC = V_L - ( (1-x_PG) * V_H_POPE + x_PG * V_H_PG )

	# Quasi-molecular volumes 
	V_CH2	= V_HC / ( n_CH2 + n_CH*r12 + n_CH3*r32 )	# Volume of CH2 groups  
	V_CH	= V_CH2	*	r12 							# Volume of CH groups 
	V_CH3	= V_CH2	*	r32 							# Volume of CH3 groups

	V_CG	= V_H_POPE	*	r_CG				# Volume of CG group 
	V_P		= V_H_POPE	* 	r_P 				# Volume of P group 
	V_ENX	= V_H_POPE	*   (1 - r_CG - r_P)	# Volume of ENX group  

	V_PG1	= V_H_PG	*	0.16 				# Pan et al. BBA 2012
	V_PG2	= V_H_PG	*	( 1 - 0.51 - 0.16) 	# Pan et al. BBA 2012
	
	return V_L, V_HC, V_CH2, V_CH, V_CH3, V_CG, V_P, V_ENX, V_PG1, V_PG2

#########################################################
def SLDs_POPE(T, V_ENX, V_P, V_CG, V_PG2, V_PG1, V_CH, V_CH2, V_CH3, V_BW) :
	"""Return SLDs of solvent and quasi-molecular lipi groups (AA^-2)"""
	
	############### X-ray scattering length densities (nm^-2)
	rho_sol		= b_HW / Water_volume(T) 
	drho_ENX	= ( (1-x_PG)*b_ENX/V_ENX	+ x_PG*b_PG2/V_PG2 )    - rho_sol 
	drho_P		= b_P*( (1-x_PG)/V_P		+ x_PG/V_PG1 )	        - rho_sol 
	drho_CG		= b_CG  / V_CG  									- rho_sol 
	drho_CH		= b_CH  / V_CH  									- rho_sol 
	drho_CH2	= b_CH2 / V_CH2  									- rho_sol 
	drho_CH3	= b_CH3 / V_CH3  									- rho_sol 
	drho_BW		= b_HW	/ V_BW										- rho_sol 	
	
	return rho_sol, drho_ENX, drho_P, drho_CG, drho_CH, drho_CH2, drho_CH3, drho_BW

#########################################		VESICLE SCATTERING MODEL	######################################### 

#########################################		SDP MODELLING			######################################### 
####################################		SEPARATED FORM FACTOR			#################################### 
##################################################################################################################

class LUV_POPE:

##################
	def __init__(self, q, PAR, _) :
		self.q = q
		[self.Norm, self.nv, 
		self.Rm, self.Z, 
		self.d_ENX, self.s_ENX, self.d_P, self.s_P, self.d_CG, self.s_CG, 
		self.dD_C, self.r_D_C, 
		self.s_CH2, self.d_CH, self.s_CH, self.s_CH3, 
		self.r_P, self.r_CG, self.r12, self.r32, 
		self.T, self.V_H, self.V_BW, 
		self.Con] = PAR

		### polydispersity value 		
		alp = self.Rm/(self.Z+1)

		### Volumes 
		[_, 
		V_HC, V_CH2, self.V_CH, self.V_CH3, 
		self.V_CG, self.V_P, self.V_ENX, 
		self.V_PG1, self.V_PG2] = Volumes_POPE(self.V_H, self.T, self.r12, self.r32, self.r_CG, self.r_P)
			
		### Reduce to single lipid thickness
		self.D_C = 0.5*self.dD_C
		### Get mean area per lipid (A_L)
		self.A_L = V_HC / self.D_C
		
		### X-ray scattering length densities (nm^-2) 	
		[_, self.drho_ENX, self.drho_P, self.drho_CG,
   	 	self.drho_CH, self.drho_CH2, self.drho_CH3, 
	 	self.drho_BW] = SLDs_POPE(	self.T,
							 		self.V_ENX, self.V_P, self.V_CG, 
	 								self.V_PG2, self.V_PG1, 
	 								self.V_CH, V_CH2, self.V_CH3, self.V_BW )
			
		### D_C polydispersity
		self.s_D_C = self.r_D_C*self.D_C
		HC_array = np.linspace(self.D_C-3*self.s_D_C, self.D_C+3*self.s_D_C, Nd)
		Normal = PDF_normal(HC_array, self.D_C, self.s_D_C)

		### c-prefactors
		c_CH  = np.zeros(Nd, dtype=float)
		c_CH3 = np.zeros(Nd, dtype=float)

		c_ENX	= ( (1-x_PG)*self.V_ENX	+ x_PG*self.V_PG2 )	/ self.A_L 
		c_P		= ( (1-x_PG)*self.V_P	+ x_PG*self.V_PG1 )	/ self.A_L
		c_CG	= self.V_CG / self.A_L 

		c_CH	= self.V_CH		* n_CH	/ (V_HC / HC_array ) 
		c_CH3	= self.V_CH3	* n_CH3	/ (V_HC / HC_array ) 
	
		### Assembling membrane scattering amplitude		
		Am = np.zeros([self.q.size, Nd], dtype=float)
		
		self.d_BW = self.d_CG + self.d_P + self.d_ENX + d_shl

		Am += 2 * self.drho_CH2 * HC_array[None,:] * FTreal_erf(self.q[:,None], 0, 2*HC_array[None,:], self.s_CH2)
		# Adding CH and CH3 groups 
		Am += 2 * (self.drho_CH  - self.drho_CH2) * c_CH[None,:]  * FTreal_gauss(self.q[:,None], self.d_CH, self.s_CH) 
		Am += 2 * (self.drho_CH3 - self.drho_CH2) * c_CH3[None,:] * FTreal_gauss(self.q[:,None], 0,         self.s_CH3) 
		# Adding hydration-water envelope 
		Am += 4 * self.drho_BW * self.d_BW/2. * FTreal_erf(self.q[:,None], (HC_array[None,:]+self.d_BW/2.), self.d_BW, self.s_CH2)
		# Adding CG, PCN and CholCH3 groups 
		Am += 2 * (self.drho_CG		- self.drho_BW) * c_CG	* FTreal_gauss(self.q[:,None], (HC_array[None,:] + self.d_CG),	self.s_CG) 
		Am += 2 * (self.drho_P		- self.drho_BW) * c_P	* FTreal_gauss(self.q[:,None], (HC_array[None,:] + self.d_CG + self.d_P),	self.s_P) 
		Am += 2 * (self.drho_ENX	- self.drho_BW) * c_ENX * FTreal_gauss(self.q[:,None], (HC_array[None,:] + self.d_CG + self.d_P + self.d_ENX), self.s_ENX)

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

		self.check = 0

		CG  = Gauss(z_array, self.V_CG,		self.D_C+self.d_CG,						self.s_CG,	self.A_L)
		P   = Gauss(z_array, self.V_P,		self.D_C+self.d_CG+self.d_P,			self.s_P,	self.A_L)
		ENX = Gauss(z_array, self.V_ENX,	self.D_C+self.d_CG+self.d_P+self.d_ENX,	self.s_ENX,	self.A_L)
		BW  = Slab(z_array,	self.D_C+self.d_BW/2., self.d_BW, self.s_CH2) - CG - P - ENX

		for i in(BW) : 
			if i <-0.001 : self.check+= 0   

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
		CH		= Gauss(z_array, n_CH*self.V_CH, self.d_CH, self.s_CH, self.A_L)
		CH2		= CH2 - CH3 - CH
		
		CG		= Gauss(z_array, self.V_CG,     self.D_C+self.d_CG,                     self.s_CG,  self.A_L)
		P		= Gauss(z_array, self.V_P,		self.D_C+self.d_CG+self.d_P,			self.s_P,	self.A_L)
		ENX 	= Gauss(z_array, self.V_ENX,	self.D_C+self.d_CG+self.d_P+self.d_ENX,	self.s_ENX,	self.A_L)
		BW		= BW - CG - P - ENX
		
        # ...and DeltaSLD profile
		SLD = np.zeros_like(z_array)
		SLD += CH3 * self.drho_CH3 + CH2 * self.drho_CH2 + CH * self.drho_CH 
		SLD += CG * self.drho_CG + P * self.drho_P + ENX * self.drho_ENX
		SLD += BW * self.drho_BW

		# Calculate peak-to-peak distance D_pp (proxy for phosphate-to-phosphate distance)
		SLD_peak = np.argmax(SLD)
		D_pp = 2*z_array[SLD_peak]

		# Calculate the number of water molecules per lipid headgroup
		n_W =  ( self.A_L * self.d_BW*erf(self.d_BW/2./(np.sqrt(2)*self.s_CH2)) - self.V_H ) / self.V_BW
		
		SDP_matrix = pd.DataFrame(np.column_stack((z_array, CH3, CH, CH2, CG, P, ENX, BW, SLD)), 
								  columns=["z", "CH3", "CH", "CH2", "CG", "P", "ENX", "BW", "SLD"])

		return 	SDP_matrix, n_W, D_pp, D_B, self.D_C, self.A_L

				
##################################################################################################################
##################################################################################################################