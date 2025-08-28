"""
Shared function for (proteo)liposome models
"""

import numpy as np
from scipy.special import (spherical_jn, jv, erf)

from models.constants import *

#####################################################################################################################
""" 0th, 2nd and 4th moment (Schulz PDF) of  j_0(x) """

def mu0(q, Z, a) :
	x = 2*q*a
	return np.where(q==0, (Z-1)/Z, ( 1 - (1+x**2)**(-(Z-1)/2.) * np.cos((Z-1)*np.arctan(x)) ) / ( Z*(Z-1)*(x)**2/2 ) )

def mu2(q, Z, a) :
	x = 2*q*a
	return np.where(q==0, a*a*(Z+1)*(Z+1), ( 1 - (1+x**2)**(-(Z+1)/2.) * np.cos((Z+1)*np.arctan(x)) ) / ( 2*q**2 ) )
		
def mu4(q, Z, a) :
	x = 2*q*a
	return np.where(q==0, a**4*(Z+1)*(Z+2)*(Z+3)**2,  a**2 * (Z+2)*(Z+1) * ( 1 - (1+x**2)**(-(Z+3)/2.) * np.cos((Z+3)*np.arctan(x)) ) / ( 2*q**2 ) )	

#####################################################################################################################
""" Auxiliary functions:
    Real and inverse space description of SLD profile, slaba and Gauss curves
"""

def PDF_normal(x, mu, sig) :
    """ Normal PDF """
    return np.exp(-(x-mu)**2 / (2*sig**2) ) /( sig*np.sqrt(2*np.pi))

def FTreal_erf(q, mu, d, sig) :
    """ Fourier trasnform of a smooth slab profile (error-function-based) """
    return spherical_jn(0, q*d/2.) * np.exp(-(q*sig)**2/2.) * np.cos(q*mu) 
 
def FTreal_gauss(q, mu, sig) :
    """ Fourier trasnform of a Gaussian profile """
    return np.exp(-(q*sig)**2/2.) * np.cos(q*mu)

def Slab(x, mu, L, sig) :
    """ Smooth slab profile (error-function-based) (SLD profile)"""
    return 0.5 * ( erf( (x - (mu-L/2.))/(np.sqrt(2)*sig) ) - erf( (x - (mu+L/2))/(np.sqrt(2)*sig) ) )

def Gauss( x, V, mu, sig, A_L ) :
    """ Gaussian profile (SLD profile)"""
    return V * PDF_normal(x, mu, sig) / A_L

#####################################################################################################################
""" Auxiliary functions:
    Radial and axial components of the scattering amplitude of a cylinder
"""

def Lambda_t_np(q, d, t) :
    """ Normal PDF """
    return d * spherical_jn( 0, q * t * d/2.) 

def Beta_t_np(q, R, t) :
    """ Normal PDF """
    return np.where( t == 1, np.pi*R**2, np.pi*R**2 * 2 * jv(1, q*R*np.sqrt(1-t**2)) / (q*R*np.sqrt(1-t**2)) )

#####################################################################################################################

def virtual(q, dR, dH, t, R_base, V_ompla, rho_sol) :
	"""Normalized scattering amplitude of OmpLA virtual model
	
	SLDs of the single disks are not included (assuming homogeneous SLD)
	"""
	# initialize radii and thicknesses of stacks
	R = R_base*dR
	Nst = R.shape[0] # number of of stacks (N-1)*2+1
	H_frac = H_ompla/((Nst-1)*2+1)
	H = H_frac*dH
	# stacks' radial term
	B = np.empty([Nst, q.shape[0], t.shape[0]])
	B = Beta_t_np(q[None,:,None], R[:,None,None], t[None,None,:])
	# stacks' director term
	L = np.empty([Nst, q.shape[0], t.shape[0]])
	L = Lambda_t_np(q[None,:,None], H[:,None,None], t[None,None,:])   	
	# stacks' cross-term depending on distances
	dist = np.zeros([Nst])
	dist[1:] = np.cumsum(H[1:])+(H[0]-H[1:])/2 #distance from center along z
	crs = np.ones([Nst,q.shape[0],t.shape[0]])
	crs[1:] = 2 * np.cos( q[None,:,None]*t[None,None,:] * dist[1:,None,None] )
	# Assembling OmpLA scattering amplitude
	A_virtual = np.zeros([q.shape[0],t.shape[0]])
	A_virtual = np.sum(B*L*crs,axis=0)
	A_virtual*= rho_ompla-rho_sol
	# Normalising by disks' volumes
	V_disks = np.pi*(R[0]**2*H[0]+2*np.sum(R[1:]**2*H[1:]))
	A_virtual/=V_disks
	A_virtual*=V_ompla
	
	# Add top and bottom bound-water layers
	B_W = Beta_t_np(q[:,None], R_base, t[None,:])
	L_W = Lambda_t_np(q[:,None], d_shl, t[None,:])   	
	dist_W = H_ompla/2 + d_shl/2 #distance from center along z
	crs_W = 2 * np.cos( q[:,None]*t[None,:] * dist_W )
	A_virtual+= ((b_HW/V_PW)-rho_sol) * B_W*L_W*crs_W

	return A_virtual

#####################################################################################################################

def Water_volume(T) :
	"""Return the volume of water molecules (AA) as a function of temperature"""
	# polynome coefficient for T-dependency of bulk-water-molecule volume (V_HW) 
	# Units in degree Celsius 
	p0_VW	=	29.9218322593344
	p1_VW	=	-0.00225941007461549
	p2_VW	=	0.000256750019262826
	p3_VW	=	-1.69660959036021e-06
	p4_VW	=	6.52029089103885e-09
	return p0_VW + p1_VW*T + p2_VW*T**2 + p3_VW*T**3 + p4_VW*T**4	

def RecBuf_mol_ratio(T) :
    """ Reconstitution buffer
    Return mole fraction of TRIS and EDTA in buffer as a function of temperature
    """
      
    # Composition of the Reconstitution buffer (M) 
    ctris = 0.02
    cEDTA = 0.002
    
    # polynome coefficient for T-dependency of bulk-water molar concentration (Cw) 
    # Units in degree Celsius 
    p0_Cw	= 55.5052     
    p1_Cw	= 0.00131894     
    p2_Cw	= -0.000334396     
    p3_Cw	= 9.10861e-07  

    Cw = p0_Cw + p1_Cw*T + p2_Cw*T**2 + p3_Cw*T**3       
    xtris = ctris / Cw # mole fraction of free TRIS in bulk 	
    xEDTA = cEDTA / Cw # mole fraction of free EDTA in bulk 
    return xtris, xEDTA

#####################################################################################################################
