#!/usr/bin/python

import os
import sys
import time
import math
import pandas as pd
import yaml

###############
############### Miscellaneous

def file_exists(arg, name):
	if os.path.exists(arg):
		pass
		#print(name,"\t",arg)
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
	
	def __init__ ( self , ARGV ):

		if len(ARGV)>3:
			sys.exit("--- Too many arguments")
		elif len(ARGV)==1:
			sys.exit("--- Specify the parameters file")

		print("-----------------------------------------------------")
		localtime = time.asctime( time.localtime(time.time()) )
		print(localtime)

		file_exists(ARGV[1], "\n-- Parameters File")

		# Reading cofigurations from the YAML file
		with open(ARGV[1], 'r', encoding='utf8') as file:
			par_file = yaml.safe_load(file)

		self.cfg = par_file['config']
	
		self.param = pd.DataFrame.from_dict(par_file['parameters'], orient='index', columns=['value', 'free', 'prior', 'low_l', 'high_l']).reset_index()
		self.param = self.param.rename(columns={"index": "name"})

###########################################################################################
	def get_config ( self ):
		""" Load configuration	"""

		# Check datafile
		file_exists(self.cfg['datafile'], "\n- datafile:\t")

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
		is_float_positive(self.cfg['temperature-gain'])
		is_float_positive(self.cfg['target-X2'])

		is_int_positive(self.cfg['data-binning'])

		is_within_01(self.cfg['error-scale'])

		is_int_positive(self.cfg['iterations'])
		is_int_positive(self.cfg['processes'])

		if isinstance(self.cfg['plot-only'], bool):
				pass
		else:
			raise Exception(self.cfg['plot-only'], "must be Boolean")

		if self.cfg['state']=="monomer" or self.cfg['state']=="dimer":
				pass
		else:
			raise Exception(self.cfg['state'], "must be either (monomer) or (dimer)")
		
		# print configuration recap
		lmax = max(len(key) for key, value in self.cfg.items())
		for key, value in self.cfg.items():
			key=key+str(":")
			print(key.ljust(lmax+1), value)

		return self.cfg
	
###########################################################################################
	def get_parameters ( self ):
		'Load parameter matrix'	

		for index, row in self.param.iterrows():
			if not math.isnan(row['prior']):
				self.param.at[index, 'low_l'] =  row['value']*(1-5*row['prior'])
				self.param.at[index, 'high_l'] = row['value']*(1+5*row['prior'])
			#else:
				#self.param.at[index, 'prior'] =  0

		print("\n-- Parameter matrix:\n", self.param)

		return self.param
	
###########################################################################################	
			
