#!/usr/bin/python

import os
import sys
import math
import pandas as pd
import yaml
import numpy as np

###############
############### Miscellaneous

def file_exists(arg):
	if os.path.exists(arg):
		pass
	else:
		raise IOError(arg,"The file does not exist")
##############

def is_float_positive(val):
	if isinstance(val, (int, float)) and isinstance(val, bool)==False:
		if float(val) >= 0:
			return float(val)
		else:
			raise Exception(val," must be a positive float number")
	else:
		raise ValueError(val," must be a float number")

###############

def is_int_positive(val):
	if isinstance(val, int) and isinstance(val, bool)==False:
		if int(val) >= 0:
			return int(val)
		else:
			raise Exception(val," must be a positive integer number")
	else:
		raise ValueError(val," must be an integer number")

###############

def is_within_01(val):
	if isinstance(val, (int, float)) and isinstance(val, bool)==False:
		if float(val) >= 0 and float(val) <= 1:
			return float(val)
		else:
			raise Exception(val," must be a within 0 and 1")
	else:
		raise ValueError(val," must be a float number number")		

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

class Load_input:
	'Read the Options & Parameters file'
	
	def __init__ ( self , ARGV, settings_file):

		if len(ARGV)>3:
			sys.exit("--- Too many arguments")
		elif len(ARGV)==1:
			sys.exit("--- Specify the parameters file")

		file_exists(ARGV[1])
		file_exists(settings_file)

		# Reading settings from the YAML file
		with open(settings_file, 'r', encoding='utf8') as file:
			self.settings = yaml.safe_load(file)

		# Reading cofigurations from the YAML file
		with open(ARGV[1], 'r', encoding='utf8') as file:
			par_file = yaml.safe_load(file)

		self.cfg = par_file['config']
	
		self.param = pd.DataFrame.from_dict(par_file['parameters'],
																			  orient='index',
																			  columns=['value', 'free', 'prior', 'low_l', 'high_l']).reset_index()
		self.param = self.param.rename(columns={"index": "name"})

###########################################################################################
	def get_settings ( self ):
		""" Load settings	"""

		# Check numerical or Boolean values

		is_float_positive(self.settings['convergence']['thermo'])
		is_int_positive(self.settings['convergence']['maxcount'])
		is_float_positive(self.settings['convergence']['conv_threshold'])
		is_float_positive(self.settings['convergence']['neg_water_scale'])

		if isinstance(self.settings['statistics']['scale_MAD_to_std'], bool):
				pass
		else:
			raise Exception(self.settings['statistics']['scale_MAD_to_std'], "must be Boolean")

		is_int_positive(self.settings['statistics']['hist_mincount'])

		# print settings  recap
		print("# Settings")
		print("convergence:")
		lmax = max(len(key) for key, value in self.settings['convergence'].items())
		for key, value in self.settings['convergence'].items():
			key=key+str(":")
			print("\t", key.ljust(lmax+1), value)
		print("statistics:")
		lmax = max(len(key) for key, value in self.settings['statistics'].items())
		for key, value in self.settings['statistics'].items():
			key=key+str(":")
			print("\t", key.ljust(lmax+1), value)
		print("-----------------------------------")

		return self.settings

###########################################################################################
	def get_config ( self ):
		""" Load configuration	"""

		# Check datafile
		file_exists(self.cfg['datafile'])

		# Check and create if not existing destination folder
		self.cfg['save-folder'] = str("RES_")+self.cfg['save-folder']
		os.makedirs("./"+self.cfg['save-folder'], exist_ok=True)

		# Check numerical or Boolean values
		for el in self.cfg['qrange']:
			is_float_positive(el)
		if self.cfg['qrange'][0] < self.cfg['qrange'][1]:
				pass
		else:
			raise Exception("q_min must be lower than q_max")
	
		is_float_positive(self.cfg['temperature-init'])
		#is_float_positive(self.cfg['temperature-gain'])
		is_float_positive(self.cfg['target-X2'])

		is_int_positive(self.cfg['data-binning'])

		is_within_01(self.cfg['error-scale'])

		is_int_positive(self.cfg['iterations'])
		is_int_positive(self.cfg['processes'])

		if isinstance(self.cfg['plot-only'], bool):
				pass
		else:
			raise Exception(self.cfg['plot-only'], "must be Boolean")

		# Check if 'state' parameter (for pLUVs) exists 
		if 'state' in self.cfg:
			if self.cfg['state']=="monomer" or self.cfg['state']=="dimer":
					pass
			else:
				raise Exception(self.cfg['state'], "must be either (monomer) or (dimer)")
		else:
			pass
		
		# print configuration recap
		print("# Configuration")
		lmax = max(len(key) for key, value in self.cfg.items())
		for key, value in self.cfg.items():
			key=key+str(":")
			print(key.ljust(lmax+1), value)
		print("-----------------------------------")

		return self.cfg
	
###########################################################################################
	def get_parameters ( self ):
		'Load parameter matrix'

		for index, row in self.param.iterrows():

			# Check if the parameter initialization is a number
			if not isinstance(row['value'], (int, float)):
				raise ValueError("Initialization value must be a float number (parameter initialization)")

			# Check if the free/fixed entry is boolean
			if not isinstance(row['free'], (bool)):
				raise ValueError("Free/fix entry must be a boolean value (on/off or True/False)")

			# Check if the prior entry is a number or a nan value
			if not isinstance(row['prior'], (float)):
				if row['prior'] is not None:
						raise ValueError("The relative sigma of the prior must be either a float value or Null (NAN)")
			# Check if the prior entry is a number or a nan value
			if row['prior'] is not None and row['prior']<=0:
				raise ValueError("The relative sigma of the prior must be a positive float value")

			# Check if the parameter initialization is a number
			if not isinstance(row['high_l'], (int, float)):
				raise ValueError("High boundary must be either a float value or Null (NAN)")
			
			# Check if the baudary entries are numbers or a nan values, and if they are consistent
			if not isinstance(row['low_l'], (float)):
				raise ValueError("Low boundary must be either a float value or Null (NAN)")
			#elif not isinstance(row['high_l'], (float)):
			#	raise ValueError("High boundary must be either a float value or Null (NAN)")
			elif (not math.isnan(row['low_l']) and math.isnan(row['high_l'])) or (math.isnan(row['low_l']) and not math.isnan(row['high_l'])):
				raise ValueError("Both low and high limits must be either float values or Null (NAN)")

			# Check if at leat one among boundaries or priors are set
			if math.isnan(row['low_l']) and math.isnan(row['high_l']) and math.isnan(row['prior']):
				raise ValueError("set at least one between prior and low-high bondaries entries")

		# force 'prior' column to numeric in the case of all Null entries
		self.param['prior'] = pd.to_numeric(self.param['prior'], errors='coerce')

		for index, row in self.param.iterrows():
			# Overwrite lw and high limits is the prior is set
			if not math.isnan(row['prior']):
				self.param.at[index, 'low_l'] =  row['value']*(1-5*row['prior'])
				self.param.at[index, 'high_l'] = row['value']*(1+5*row['prior'])


		print("\n# Parameter matrix:\n", self.param)

		return self.param
	
###########################################################################################	
			
