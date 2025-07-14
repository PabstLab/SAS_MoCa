#!/usr/bin/python

import pandas as pd
import numpy as np

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

class Load_data:
	'Read data file'
	
	def __init__ ( self , datafiles, qrange, bing) :
		
		# Load datafile 
		file = pd.read_csv(datafiles, sep='\s+', header=None, names=['q', 'I', 'error', 'Dq'], engine='python', skiprows=lambda r: r % bing)
		# Restrict data to selected q-range
		file = file[(file['q'] >= qrange[0]) & (file['q'] <= qrange[1])]
		# Remove Dq column if empty
		if file['Dq'].all():
			file = file[['q', 'I', 'error']]

		# Convert pandas dataframes to numpy array
		self.file = file.to_numpy()
		
###########################################################################################
###########################################################################################

	def convert(self, err_mul):
	
		# Set minimum relative error
		min_rel = 0.02

		###### Reshape error values
		#----- Set minimum relative error to min_rel
		#----- If err_mul == 1 leave the errorbars untouched
		self.file[:,2] = self.file[:,2]*err_mul
		if err_mul != 1 :
			self.file[:,2] = np.where( self.file[:,2] < self.file[:,1]*min_rel, self.file[:,1]*min_rel, self.file[:,2] )

		###### Fix minimum relative error in 2-columns only datasets
		self.file[:,2] = np.where( np.isnan(self.file[:,2]), self.file[:,1]*min_rel, self.file[:,2] )

		return self.file