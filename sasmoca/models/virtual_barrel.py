#!/usr/bin/python

import numpy as np
from scipy.special import spherical_jn, jv

H_ompla = 44.0 # AA from PyMol
R_ompla = 17.2 # AA from PyMol
R_ompla_dimer = 20.3 # AA from PyMol

###################################################################################################### 
###################################################################################################### 
###################################################################################################### 
###############		DUMMY PROTEIN MODEL IN VACUO   														############### 
###############		z-SYMMETIC STACK OF DISKS (ROTATION SYMMETRY AROUN AN AXIS - CYLINDRICAL GEOMETRY) 	############### 
###################################################################################################### 
###################################################################################################### 
###################################################################################################### 


#########################################################	
def Lambda_t_np(q, d, t) :
	return d * spherical_jn( 0, q * t * d/2.) 

def Beta_t_np(q, R, t) :
	return np.where( t == 1, np.pi*R**2, np.pi*R**2 * 2 * jv(1, q*R*np.sqrt(1-t**2)) / (q*R*np.sqrt(1-t**2)) )

##################################################################################################################
##################################################################################################################

###################################################################

class Stacked_disks_singleSLD:
	"""OmpLA monomer virtual model
	
	Symmetric n-stacked disks with equal SLD
	"""

##################
	def __init__(self, q, PAR, _) :
		[self.Norm, self.n,
		dR0, dR1, dR2, dR3, dR4, dR5, dR6,
		dH0, dH1, dH2, dH3, dH4, dH5, dH6,
		self.Con] = PAR

		R = R_ompla*np.array([dR0, dR1, dR2, dR3, dR4, dR5, dR6])
		N = R.shape[0] # number of of stacks (N-1)*2+1
		H_frac = H_ompla/((N-1)*2+1)
		H = H_frac*np.array([dH0, dH1, dH2, dH3, dH4, dH5, dH6])

		### t array for dummy-object orientation
		Nt = 45
		t = np.linspace(0, 1, Nt)

		B = np.empty([N,q.shape[0],Nt])
		L = np.empty([N,q.shape[0],Nt])
		crs = np.ones([N,q.shape[0],Nt])
		dist = np.zeros([N])

		B = Beta_t_np(q[None,:,None], R[:,None,None], t[None,None,:])
		L = Lambda_t_np(q[None,:,None], H[:,None,None], t[None,None,:])   	
		dist[1:] = np.cumsum(H[1:])+(H[0]-H[1:])/2 #distance from center along z
		crs[1:] = 2 * np.cos( q[None,:,None]*t[None,None,:] * dist[1:,None,None] )

		# Assembling OmpLA scattering amplitude
		I_dummy = np.zeros([q.shape[0],Nt])
		I_dummy = (np.sum(B*L*crs,axis=0))**2
		# Orietation average
		self.Iori = np.empty([q.shape[0]],dtype=float)
		self.Iori = np.trapz(I_dummy, t, 1/(Nt-1), axis=1)
		# Normalising by disk' volumes
		V = R[0]**2*H[0]
		for i in range(N): V+=2*R[i]**2*H[i]
		V*=2*np.pi
		self.Iori/=V**2

##################			
	def intensity(self):
		return  (self.Norm * self.n)* self.Iori + self.Con
		
##################			
	def negative_water(self):

		self.check = 0

		return self.check
				
##################################################################################################################
##################################################################################################################

###################################################################

class Stacked_disks_singleSLD_dimer:
	"""OmpLA dimer virtual model
	
	Symmetric n-stacked disks with equal SLD
	"""
	
##################
	def __init__(self, q, PAR, _) :
		self.q = q
		[self.Norm, self.n,
		dR0, dR1, dR2, dR3, dR4,
		dH0, dH1, dH2, dH3, dH4, 
		self.Con] = PAR

		R = R_ompla_dimer*np.array([dR0, dR1, dR2, dR3, dR4])
		N = R.shape[0] # number of of stacks (N-1)*2+1
		H_frac = H_ompla/((N-1)*2+1)
		H = H_frac*np.array([dH0, dH1, dH2, dH3, dH4])

		### t array for dummy-object orientation
		Nt = 45
		t = np.linspace(0, 1, Nt)

		B = np.empty([N,q.shape[0],Nt])
		L = np.empty([N,q.shape[0],Nt])
		crs = np.ones([N,q.shape[0],Nt])
		dist = np.zeros([N])

		B = Beta_t_np(q[None,:,None], R[:,None,None], t[None,None,:])
		L = Lambda_t_np(q[None,:,None], H[:,None,None], t[None,None,:])   	
		dist[1:] = np.cumsum(H[1:])+(H[0]-H[1:])/2 #distance from center along z
		crs[1:] = 2 * np.cos( q[None,:,None]*t[None,None,:] * dist[1:,None,None] )

		# Assembling OmpLA scattering amplitude
		I_dummy = np.zeros([q.shape[0],Nt])
		I_dummy = (np.sum(B*L*crs,axis=0))**2
		# Orietation average
		self.Iori = np.empty([q.shape[0]],dtype=float)
		self.Iori = np.trapz(I_dummy, t, 1/(Nt-1), axis=1)
		# Normalising by disk' volumes
		V = R[0]**2*H[0]
		for i in range(N): V+=2*R[i]**2*H[i]
		V*=2*np.pi
		self.Iori/=V**2

##################			
	def intensity(self):
		return  (self.Norm * self.n)* self.Iori + self.Con
		
##################			
	def negative_water(self):

		self.check = 0

		return self.check
				
##################################################################################################################
##################################################################################################################
