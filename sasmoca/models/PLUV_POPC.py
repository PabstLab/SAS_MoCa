#!/usr/bin/python

import os
import numpy as np
import pandas as pd
from scipy.special import erf

from models.shared import (mu0, mu2, mu4,
                           PDF_normal,
                           FTreal_erf, FTreal_gauss,
                           Slab, Gauss,
                           Beta_t_np, virtual,
                           Water_volume, RecBuf_mol_ratio)

from models.constants import *

###################################################################################################### 
###################################################################################################### 
###################################################################################################### 
###############		POPC/POPG 95/5 mol/mol LUVs						############### 
###############		Reconstitution buffer: TRIS 20 mM, EDTA 2 mM 	############### 
###############		Buffer used for PLUV exp. in DESY in 2017		############### 
###################################################################################################### 
###################################################################################################### 
###################################################################################################### 

#	##############	POPC and POPG	############### 
#	1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine 
#	1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoglycerol 

############### GLOBAL VARIABLES ###############

# POPG molar ratio
x_PG  = 0.05

# number of chain groups
n_CH  = 2 
n_CH2 = 28
n_CH3 = 2  		

################################################################################################################## 
################################################################################################################## 


#########################################################
def Lipid_volume(T) :
	"""Return the average lipid volume including a molar fraction of POPG"""
	return (1-x_PG) * (a0_V_POPC +T * a1_V_POPC) + x_PG * (a0_V_POPG + T * a1_V_POPG) 

#########################################################
def Volumes_POPC(V_H_POPC, T, r12, r32, r_CG, r_PCN) :	
	"""Return volumes of each quasi-molecular group"""

	# Volumes 
	V_L = Lipid_volume(T)
	V_HC = V_L - ( (1-x_PG) * V_H_POPC + x_PG * V_H_POPG )

	# Quasi-molecular volumes 
	V_CH2	= V_HC / ( n_CH2 + n_CH*r12 + n_CH3*r32 )	# Volume of CH2 groups  
	V_CH	= V_CH2	*	r12 							# Volume of CH groups 
	V_CH3	= V_CH2	*	r32 							# Volume of CH3 groups

	V_CG	= V_H_POPC	*	r_CG				# Volume of CG group 
	V_PCN	= V_H_POPC	* 	r_PCN 				# Volume of PCN group 
	V_Chol	= V_H_POPC	*   (1 - r_PCN - r_CG)	# Volume of CholCH3 group  

	V_PG1	= V_H_POPG	*	0.16 				# Pan et al. BBA 2012
	V_PG2	= V_H_POPG	*	( 1 - 0.51 - 0.16) 	# Pan et al. BBA 2012
	
	return V_L, V_HC, V_CH2, V_CH, V_CH3, V_CG, V_PCN, V_Chol, V_PG1, V_PG2

#########################################################
def SLDs_POPC(T, xtris, xEDTA, V_Chol, V_PCN, V_CG, V_PG2, V_PG1, V_CH, V_CH2, V_CH3, V_BW) :
	"""Return SLDs of solvent and quasi-molecular lipi groups (AA^-2)"""
	
	############### X-ray scattering length densities (nm^-2)
	rho_sol		= ( b_HW + xtris*b_tris + xEDTA*b_EDTA ) / Water_volume(T) 
	drho_Chol	= ( (1-x_PG)*b_Chol/V_Chol	+ x_PG*b_PG2/V_PG2 )	- rho_sol 
	drho_PCN	= ( (1-x_PG)*b_PCN/V_PCN	+ x_PG*b_PG1/V_PG1 )	- rho_sol 
	drho_CG		= b_CG  / V_CG  									- rho_sol 
	drho_CH		= b_CH  / V_CH  									- rho_sol 
	drho_CH2	= b_CH2 / V_CH2  									- rho_sol 
	drho_CH3	= b_CH3 / V_CH3  									- rho_sol 
	drho_BW		= b_HW	/ V_BW										- rho_sol 	
	
	return rho_sol, drho_Chol, drho_PCN, drho_CG, drho_CH, drho_CH2, drho_CH3, drho_BW

#########################################		PROTEOLIPOSOME MODEL	######################################### 

#########################################		SDP MODELLING			######################################### 
####################################		SEPARATED FORM FACTOR			#################################### 
#################################		HOMOGENEOUSLY DISTRIBUTED PROTEINS			################################# 
###############################		DESCRIBED BY SYMMETRIC STACKS OF DISKS			############################### 
##################################################################################################################

class pLUV_POPC_OmpLA_RecBuf:

##################
	def __init__(self, q, PAR, ompla_state) :
		self.q = q
		[self.Norm, self.nv, 
		self.Rm, self.Z, 
		self.xp,
		self.d_Chol, self.s_Chol, self.d_PCN, self.s_PCN, self.d_CG, self.s_CG, 
		self.dD_C, self.r_D_C, 
		self.s_CH2, self.d_CH, self.s_CH, self.s_CH3, 
		self.r_PCN, self.r_CG, self.r12, self.r32, 
		self.T, self.V_H, self.V_BW, 
		self.Con] = PAR

		# load either monomer or dimer virtual parameters
		if ompla_state == 'monomer':
			R_om, H_om = np.load(str(os.path.dirname(os.path.realpath(__file__)))+"/ompla_virtual_parameters_monomer.npy")
			R_ompla		= R_ompla_monomer
			V_ompla		= V_ompla_monomer

		elif ompla_state == 'dimer':
			R_om, H_om = np.load(str(os.path.dirname(os.path.realpath(__file__)))+"/ompla_virtual_parameters_dimer.npy")
			R_ompla		= R_ompla_dimer
			V_ompla		= V_ompla_dimer

		# Assign carrot size
		R_crt	= R_ompla

		# renormalization of protein-to-lipid ratio
		self.xp*=1e-2

		### polydispersity value 		
		alp = self.Rm/(self.Z+1)

		### Volumes 
		[V_L, 
		V_HC, V_CH2, self.V_CH, self.V_CH3, 
		self.V_CG, self.V_PCN, self.V_Chol, 
		self.V_PG1, self.V_PG2] = Volumes_POPC(self.V_H, self.T, self.r12, self.r32, self.r_CG, self.r_PCN)
			
		### Reduce to single lipid thickness
		self.D_C = 0.5*self.dD_C
		### Get mean area per lipid (A_L)
		self.A_L = V_HC / self.D_C

		### Reconstitution buffer: molar ratios		
		xtris, xEDTA = RecBuf_mol_ratio(self.T)
		
		### X-ray scattering length densities (nm^-2) 	
		[rho_sol, self.drho_Chol, self.drho_PCN, self.drho_CG,
   	 	self.drho_CH, self.drho_CH2, self.drho_CH3, 
	 	self.drho_BW] = SLDs_POPC(	self.T, xtris, xEDTA,
							 		self.V_Chol, self.V_PCN, self.V_CG, 
	 								self.V_PG2, self.V_PG1, 
	 								self.V_CH, V_CH2, self.V_CH3, self.V_BW )
			
		### D_C polydispersity
		self.s_D_C = self.r_D_C*self.D_C
		HC_array = np.linspace(self.D_C-3*self.s_D_C, self.D_C+3*self.s_D_C, Nd)
		Normal = PDF_normal(HC_array, self.D_C, self.s_D_C)

		### t array for OmpLA orientation
		t = np.linspace(0, 1, Nt)
		
		### c-prefactors
		c_CH  = np.zeros(Nd, dtype=float)
		c_CH3 = np.zeros(Nd, dtype=float)

		c_Chol	= ( (1-x_PG)*self.V_Chol	+ x_PG*self.V_PG2 )	/ self.A_L 
		c_PCN	= ( (1-x_PG)*self.V_PCN		+ x_PG*self.V_PG1 )	/ self.A_L
		c_CG	= self.V_CG / self.A_L 

		c_CH	= self.V_CH		* n_CH	/ (V_HC / HC_array ) 
		c_CH3	= self.V_CH3	* n_CH3	/ (V_HC / HC_array ) 
	
		### Assembling membrane scattering amplitude		
		Am = np.zeros([self.q.size, Nd, Nt], dtype=float)
		
		self.d_BW = self.d_CG + self.d_PCN + self.d_Chol + d_shl

		# Adding hydrocarbon-chain envelope 
		Am += 2 * self.drho_CH2 * HC_array[None,:,None] * FTreal_erf(self.q[:,None,None]*t[None,None,:], 0, 2*HC_array[None,:,None], self.s_CH2)
		# Adding CH and CH3 groups 
		Am += 2 * (self.drho_CH  - self.drho_CH2) * c_CH[None,:,None]  * FTreal_gauss(self.q[:,None,None]*t[None,None,:], self.d_CH, self.s_CH) 
		Am += 2 * (self.drho_CH3 - self.drho_CH2) * c_CH3[None,:,None] * FTreal_gauss(self.q[:,None,None]*t[None,None,:], 0,         self.s_CH3) 
		# Adding hydration-water envelope 
		Am += 4 * self.drho_BW * self.d_BW/2. * FTreal_erf(self.q[:,None,None]*t[None,None,:], (HC_array[None,:,None]+self.d_BW/2.), self.d_BW, self.s_CH2)
		# Adding CG, PCN and CholCH3 groups 
		Am += 2 * (self.drho_CG		- self.drho_BW) * c_CG	* FTreal_gauss(self.q[:,None,None]*t[None,None,:], (HC_array[None,:,None] + self.d_CG),	self.s_CG) 
		Am += 2 * (self.drho_PCN	- self.drho_BW) * c_PCN	* FTreal_gauss(self.q[:,None,None]*t[None,None,:], (HC_array[None,:,None] + self.d_CG + self.d_PCN), self.s_PCN) 
		Am += 2 * (self.drho_Chol	- self.drho_BW) * c_Chol* FTreal_gauss(self.q[:,None,None]*t[None,None,:], (HC_array[None,:,None] + self.d_CG + self.d_PCN + self.d_Chol), self.s_Chol)

		# Calculating vesicle scattering intensity from OmpLA-vesicle interactions; HC-dependent
		self.Iv = Am[:,:,Nt-1]**2 * 16*np.pi**2*mu4(self.q[:,None], self.Z, alp)
		
		# Initializing final intensity array
		self.Ipoly = np.empty([self.q.size], dtype=float)	

		# Calculating overall intensity in the presence of OmpLA
		if self.xp != 0 :

			### Calculating number of OmpLA per vesicle
			N_L = np.empty(Nd, dtype=float)	
			N_L = 8*np.pi* ( self.Rm**2*(self.Z+2)/(self.Z+1) )  / ( self.A_L + 2*np.pi*R_crt**2*self.xp ) 

			Np = self.xp*N_L

			# Getting OmpLA scattering amplitude, scaled by OmpLA volume and SLD
			A_virtual = virtual(q, R_om, H_om, t, R_ompla, V_ompla, rho_sol)

			# Calculating the "carrot" scattering amplitude
			Bcrt	= Beta_t_np(self.q[:,None], R_crt, t[None,:])  	

			A_crt	= Bcrt[:,None,:]*Am	
			# Effective OmpLA amplitude
			#A_ompla = A_virtual[:,None,:] * np.cos(self.q[:,None,None]*t[None,None,:]*self.delta) - A_crt
			A_ompla = A_virtual[:,None,:] - A_crt
			# Effective OmpLA intensity
			#I_ompla = A_virtual[:,None,:]**2 + A_crt**2 - 2*A_virtual[:,None,:]*A_crt*np.cos(self.q[:,None,None]*t[None,None,:]*self.delta)
			I_ompla = A_virtual[:,None,:]**2 + A_crt**2 - 2*A_virtual[:,None,:]*A_crt
			
			# Calculating beta-factor
			beta = ( np.trapz(A_ompla, t, 1/(Nt-1), axis=2) )**2
			beta/= np.trapz( I_ompla, t, 1/(Nt-1), axis=2)

			# Assembling Spp scattering intensity from OmpLA-OmpLA interactions; HC-dependent & t-dependent
			self.Spp = Np*I_ompla * ( 1 + (Np-1) * mu0(self.q[:,None,None],self.Z,alp) * beta[:,:,None] ) # Epsilon(q)~beta(q) 

			# Assembling Spv scattering intensity from OmpLA-vesicle interactions; HC-dependent & t-dependent
			Am3 = Am[:,:,Nt-1]
			self.Spv = Am3[:,:,None] * A_ompla 
			self.Spv *= Np*8*np.pi * mu2(self.q[:,None,None],self.Z,alp)	
			
			### Calculating total scattering intensity, HC-dependent & t-dependent
			Itot = self.Iv[:,:,None] + self.Spp + self.Spv
			
			### Calculating orientation average; output is HC-dependent
			Iori = np.empty([self.q.size, Nd], dtype=float)
			Iori = np.trapz(Itot, t, 1/(Nt-1), axis=2)

			### Calculating average over D_C PDF		
			Iori[:,] *= Normal				
			self.Ipoly = np.trapz(Iori, HC_array, 6*self.s_D_C/(Nd-1), axis=1)

		# Calculating overall intensity in the absence of OmpLA
		elif self.xp == 0 :

			### Calculating average over D_C PDF		
			self.Iv[:,] *= Normal				
			self.Ipoly = np.trapz(self.Iv, HC_array, 6*self.s_D_C/(Nd-1), axis=1)

				
##################			
	def intensity(self):
		return ( self.Norm * self.nv*1e5 ) * self.Ipoly + self.Con*( 0.99*(1./(1+np.exp(-8*(self.q-0.1)))) + 0.01 )	
		
##################			
	def intensity_separated(self):
	
		### D_C polydispersity
		self.s_D_C = self.r_D_C*self.D_C
		HC_array = np.linspace(self.D_C-3*self.s_D_C, self.D_C+3*self.s_D_C, Nd)
		Normal = PDF_normal(HC_array, self.D_C, self.s_D_C)

		### t array for OmpLA orientation
		t = np.linspace(0, 1, Nt)
	
		### Calculating orientation average; output is HC-dependent
		Spp_ori = np.empty([self.q.size, Nd], dtype=float)
		Spv_ori = np.empty([self.q.size, Nd], dtype=float)

		Spp_ori = np.trapz(self.Spp, t, 1/(Nt-1), axis=2)
		Spv_ori = np.trapz(self.Spv, t, 1/(Nt-1), axis=2)
				
		### Calculating average over D_C PDF		
		Ivv_ave = np.empty([self.q.size], dtype=float)
		Spp_ave = np.empty([self.q.size], dtype=float)
		Spv_ave = np.empty([self.q.size], dtype=float)

		Ivv_ave = np.trapz(self.Iv*Normal, HC_array, 6*self.s_D_C/(Nd-1), axis=1)	
		Spp_ave = np.trapz(Spp_ori*Normal, HC_array, 6*self.s_D_C/(Nd-1), axis=1)	
		Spv_ave = np.trapz(Spv_ori*Normal, HC_array, 6*self.s_D_C/(Nd-1), axis=1)	
						
		return ( self.Norm * self.nv*1e5 ) * Ivv_ave, ( self.Norm * self.nv*1e5 ) * Spp_ave, ( self.Norm * self.nv*1e5 ) * Spv_ave 		

##################			
	def negative_water(self):

		# Calculate the Luzzati thickness
		D_B = 2*Lipid_volume(self.T)/self.A_L

		# Set z arrasy to plt SDP profile...
		z_array = np.linspace(0., D_B*1.1, int(D_B*1.1)*10, endpoint=False)

		self.check = 0
		
		CG		= Gauss(z_array, self.V_CG,			self.D_C+self.d_CG,							self.s_CG,		self.A_L)
		PCN		= Gauss(z_array, self.V_PCN,		self.D_C+self.d_CG+self.d_PCN,				self.s_PCN,		self.A_L)
		Chol	= Gauss(z_array, self.V_Chol,		self.D_C+self.d_CG+self.d_PCN+self.d_Chol,	self.s_Chol,	self.A_L)
		BW		= Slab(z_array,	self.D_C+self.d_BW/2., self.d_BW, self.s_CH2) - CG - PCN - Chol

		for i in(BW) : 
			if i <-0.001 : self.check+= 1   

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

		CG		= Gauss(z_array, self.V_CG,			self.D_C+self.d_CG,							self.s_CG,		self.A_L)
		PCN		= Gauss(z_array, self.V_PCN,		self.D_C+self.d_CG+self.d_PCN,				self.s_PCN,		self.A_L)
		Chol	= Gauss(z_array, self.V_Chol,		self.D_C+self.d_CG+self.d_PCN+self.d_Chol,	self.s_Chol,	self.A_L)
		BW		= BW - CG - PCN - Chol

		# ...and DeltaSLD profile
		SLD = np.zeros_like(z_array)
		SLD += CH3 * self.drho_CH3 + CH2 * self.drho_CH2 + CH * self.drho_CH 
		SLD += CG * self.drho_CG + PCN * self.drho_PCN + Chol * self.drho_Chol
		SLD += BW * self.drho_BW

		# Calculate peak-to-peak distance D_pp (proxy for phosphate-to-phosphate distance)
		SLD_peak = np.argmax(SLD)
		D_pp = 2*z_array[SLD_peak]

		# Calculate the number of water molecules per lipid headgroup
		n_W =  ( self.A_L * self.d_BW*erf(self.d_BW/2./(np.sqrt(2)*self.s_CH2)) - self.V_H ) / self.V_BW
		

		SDP_matrix = pd.DataFrame(np.column_stack((z_array, CH3, CH, CH2, CG, PCN, Chol, BW, SLD)), 
								  columns=["z", "CH3", "CH", "CH2", "CG", "PCN", "Chol", "BW", "SLD"])

		return 	SDP_matrix, n_W, D_pp, D_B, self.D_C, self.A_L
				
##################################################################################################################
##################################################################################################################
