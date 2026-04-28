#!/usr/bin/python

"""THERMODYNAMIC SIMULATED ANNEALING ALGORITHM"""

from __future__ import print_function
import time
import random
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from moca import FitTools

#suppress warnings
import warnings
warnings.filterwarnings('ignore')

###########################################################################################

def X2function(Q, IDATA, I, ERR, N):
	'X^2 calculation'
	X2 = np.sum(  ( ( IDATA-I ) / ERR )**2 )
	return X2 / (Q.shape[0]-(N-1)) 
	
###########################################################################################
###########################################################################################
def Min_likelyhood (prior_ratio) : 
	min = 1 
	if prior_ratio < 1 :
		min = prior_ratio
	return min
    
#------------------------- Partial ratio between priors
def PriorRatio (start, check, rel, prev) :
	new = np.exp(-(start-check)**2/(2*(rel*start)**2))
	old = np.exp(-(start-prev)**2/(2*(rel*start)**2))
	return Min_likelyhood(new/old)

###########################################################################################
###########################################################################################
def SimAnnealing ( input_list ):
	'Simulated Annealing Algorithm'
	
	# initilize input
	(data, NAME, PAR, FREE, PRIOR,
     L_LIM, H_LIM,
     function, T0, thermo, X2_trg, state, prt_progress,
     maxcount, conv_threshold, neg_water_scale) = input_list
	# initialize parameter dictionary
	par_dict = {'name': NAME, 'value': PAR, 'free': FREE, 'prior': PRIOR, 'low_l': L_LIM, 'high_l': H_LIM} 

	#####################################################
	#------------------------ Define Controls

 	# initialize counters and limits for convergence
	loop=0
	good = 0
	stim = 0

 	# initialize Boltzmann-like probability
	prBol = 0 
	
	# alpha: amplitude factor for random perturbations
	alpha		= 0.15 
	alphamax	= 0.33 
	alphamin	= 0.02

	# initialize temperature and cumulative entropy and enthalpy variables		
	T = T0
	cumDEnt	= 0
	cumDX2 	= 0
	
	# initialize stop condition
	stop = "CONTINUE"

	# set number of free parameters
	N_Free = 0
	for v in range(len(par_dict['value'])):
		if par_dict['free'][v]: N_Free+=1
		
	# initialize starting and best parameter sets
	par_start	= par_dict['value'].copy()
	par_MIN 	= par_dict['value'].copy()

	# initialize trace and likelihood array to track convergence
	trace = [[] for i in range(3+len(par_dict['value']))]
	lkhd = []

	#####################################################
	#------------------------ Evaluate first I(q) and X^2
	intensity_init = function(data[:,0], par_start, state)
	I = intensity_init.intensity()
	X2 = X2function(data[:,0], data[:,1], I, data[:,2], N_Free)
	X2_min = X2

	#####################################################
	#------------------------ Start Loop
	t0=time.time()	
	
	while stop == "CONTINUE":
		
		#####################################################
		#------------------------ Create temporary parameter set via random perturbations
		#------------------------ Check for FIX parameters & boundarys
		par_check = par_dict['value'].copy() 
		
		for p in range(len(par_dict['value'])):
			if par_dict['free'][p]:
				par_check[p]= par_dict['value'][p] + np.random.normal( 0, alpha * ( par_dict['high_l'][p] - par_dict['low_l'][p] )/6. )
				if par_check[p] < par_dict['low_l'][p] : 
					par_check[p] = par_dict['value'][p] - (par_dict['value'][p]-par_dict['low_l'][p]) / 10.
				elif par_check[p] > par_dict['high_l'][p] : 
					par_check[p] =  par_dict['value'][p] + (par_dict['high_l'][p]-par_dict['value'][p]) / 10.				

		#####################################################
		#------------------------ Evaluate temporary parameter set
		#------------------------ Add penalty in the case of negative water		
		intensity_init = function(data[:,0],par_check,state)
		check_negH2O = intensity_init.negative_water()	
		I_TEMP = intensity_init.intensity()
		X2_TEMP = X2function(data[:,0], data[:,1], I_TEMP, data[:,2], N_Free)
		if check_negH2O : X2_TEMP*= neg_water_scale

		#####################################################
		#------------------------ Calculate cumulative prior ratio between temporary and best parameter set
		prior = 1
		for v in range(len(par_dict['value'])):
			if not math.isnan(par_dict['prior'][v]):
				prior*= PriorRatio( par_start[v], par_check[v], par_dict['prior'][v], par_dict['value'][v] )		
		
		#------------------------ Calculate Boltzmann-like probability
		prBol = np.exp(-(X2_TEMP-X2)/T)
		prBol*= prior

		#####################################################
		#------------------------ Check temporary X^2 value
		#------------------------ Update configuration
		
		#----- Excellent instance
		# accept temporary set
		# update parameter set 
		# update cumulative enthapy
		if X2_TEMP < X2 :
			par_dict['value'] = par_check.copy()
			I = I_TEMP.copy()
			cumDX2+= ( X2_TEMP - X2 )			
			X2 = X2_TEMP
			good+= 1
			# Update best parameter set
			if X2 < X2_min :
				par_MIN = par_dict['value'].copy()
				X2_min = X2
			if X2_min < X2_trg :
				X2_trg = X2_min
			# Update trace
			if good % 10 == 0 : 
				trace[0].append(good)
				trace[1].append(T)
				trace[2].append(X2_TEMP)
				for p in range(len(par_check)) : 
					trace[p+3].append(par_check[p])
				lkhd.append(np.exp(-(X2_TEMP-X2_trg)/X2_trg))
				stim = np.cumsum(np.array(lkhd))[-1]

		#----- Good instance
		# accept temporary set
		# update parameter set 
		# update cumulative enthapy
		elif random.random() <= prBol :
			par_dict['value'] = par_check.copy()
			I = I_TEMP.copy()
			cumDX2+= ( X2_TEMP - X2 )			
			X2 = X2_TEMP
			good += 1
			##### Update trace
			if good % 10 == 0 : 
				trace[0].append(good)
				trace[1].append(T)
				trace[2].append(X2_TEMP)
				for p in range(len(par_check)) : 
					trace[p+3].append(par_check[p])
				lkhd.append(np.exp(-(X2_TEMP-X2_trg)/X2_trg))
				stim = np.cumsum(np.array(lkhd))[-1]

		#----- Bad instance
		# reject temporary set
		# update cumulative entropy
		else :
			cumDEnt+= -( X2_TEMP - X2 ) / T	

		#####################################################
		#------------------------ Update loop counter
		loop+= 1

		#------------------------ Update temperature scheme
		T_min = 1 + np.log10(X2_min) #1 #0.05*X2_min
		if cumDEnt == 0 or cumDX2 >= 0 : 
			T = T0
		else :
			T = thermo * cumDX2 / cumDEnt
			if T < T_min : 
				T = T_min
 
 		#------------------------ Update alpha value
		if		(good+1)/(loop+1) <= 0.62 : alpha-= alpha*0.1
		elif	(good+1)/(loop+1) >= 0.68 : alpha+= alpha*0.1
        
		if 		alpha < alphamin : alpha = alphamin
		elif 	alpha > alphamax : alpha = alphamax
	
		#####################################################
		#------------------------ Check for convergence: Plots
#		if loop % 1000 == 0 and loop > 1001 :
#			fig1, (figA) = plt.subplots(3, 1, figsize=(7.0, 4.0*3))
#			figA[0].scatter(trace[0],trace[2], marker='o', facecolors='none', color='b')
#			figA[0].axhline(y=X2_trg, color='black')
#			fig1.tight_layout()
#			figA[1].plot(trace[0],lkhd, color='b')
#			figA[1].axhline(y=0, color='black')
#			figA[2].plot(trace[0],np.cumsum(np.array(lkhd)), color='r', linewidth=3)
#			figA[2].axhline(y=0, color='black')
#			figA[2].axhline(y=1, color='black')
#			plt.show()
				
		#####################################################
		#------------------------ Time counting
		t1=time.time()-t0
		(minutes,seconds)=FitTools.TimeCount(t1)
		success = (good+1)/(loop+1)

		#####################################################
		#------------------------ Print Loop Results	
		if prt_progress == 1 :
			if loop%50 == 0 :
				out="T = %0.2f \ a = %0.2f \ min.X²= %0.5s/%0.5s \ N. %0.3d/%0.3d (%0.3f) \ %.2f/%.2f \ %2d:%2d s \ " %(T, alpha, X2_min, X2_trg, good, loop, success, stim, conv_threshold, minutes, seconds)
				print (out, end='\r')	
			
		#####################################################
		#------------------------ Set end-loop conditions
		if stim > conv_threshold or loop > maxcount:
				stop="STOP"

	#------------------------ End Loop	
	#####################################################

	# Update parameter dictionary
	par_dict['value'] = par_MIN.copy()

	return [ pd.DataFrame.from_dict(par_dict), X2_min]
	
###########################################################################################
###########################################################################################





