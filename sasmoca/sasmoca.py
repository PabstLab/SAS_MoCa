#!/usr/bin/python

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import multiprocessing
from tqdm import tqdm
import time

from in_out import (Load_input, Load_data,
					PlotData, PlotStat, PlotSDP_profile)

from models.models_list import ChooseFunction
from moca.TSA_algorithm import SimAnnealing
from moca.TSA_algorithm import X2function
from version import __version__

############################################################
############################################################
############################################################

# Protect entry point for multiprocessing
if __name__ == "__main__":

	print("-----------------------------------")
	print("# SAS_MoCa Version", __version__)
	localtime = time.asctime( time.localtime(time.time()) )
	print("#", localtime)
	print("-----------------------------------")

	# Set the start method for new processes
	multiprocessing.set_start_method('spawn')

	# Tune pandas dataframe printing options
	pd.options.display.max_colwidth = 100

	############################################################ Read configuration options & model parameters

	############### Read File & Options
	# Load configuration and parameter file
	input = Load_input(sys.argv)
	# Load configuration details
	config = input.get_config()
	# Load parameter matrix
	parameters = input.get_parameters()

	# Read model function and plotting options
	(function, proteo, LUV) = ChooseFunction( config['model'] )

	# Initialize control to print progression for each iteration
	if config['processes'] == 1 :	prt_progress = 1
	else :							prt_progress = 0

	# Print 'state' of the virtual proteins in pLUVs
	# or initialize a fictive 0 state for LUVs
	if proteo:
		print('############################', config['state'])
	else:
		config['state'] = 0

	############################################################ Read & Select Data
	data_init = Load_data( config['datafile'], config['qrange'], config['data-binning'] )
	data = data_init.convert( config['error-scale'] )
			
	############################################################ Fit Routine
	collection = []
	collection_X2 = []

	#### Plot only
	#------ show data and simulated model
	if config['plot-only']:

		par_plot = []

		for v in range(len(parameters['value'])):
			par_plot.append(parameters.iloc[v,1])

		I_plot = function(data[:,0], par_plot, config['state']).intensity()

		PlotData(config['qrange'], config['save-folder']).plot_fit(data, I_plot)	

	#### Fit data
	#------ run the minimization routine
	elif not config['plot-only']:

		#### One process only: 
		#----------- serial computation of all iterations
		#----------- display minimization stats and progress
		if config['processes'] == 1 :

			for it in range(config['iterations']) :
				print("\n---- Iteration N.", it+1,"/", config['iterations'],"-----")
				# List of inputs for minimization algorighm
				fit_inputs = [	data,
								parameters['name'].to_list(), parameters['value'].to_list(), parameters['free'].to_list(), parameters['prior'].to_list(), parameters['low_l'].to_list(), parameters['high_l'].to_list(),
								function, config['temperature-init'], config['temperature-gain'], config['target-X2'], config['state'], prt_progress]	

				# Run the minimization routine		
				(par_res, X2_min) = SimAnnealing( fit_inputs )	

				# Collect stored results
				collection.append(par_res)	
				collection_X2.append(X2_min)	

		#### Multiple processes: 
		#----------- parallel computation of all iterations
		#----------- display iterations progress
		elif config['processes'] > 1 :

			# Check OS name
			# 'nt' for Windows
			# 'posix' for Linux and Mac
			print("\n---- OS standard: ", os.name,"\n---- Platform: ", sys.platform,"\n---- N. CPU: ", multiprocessing.cpu_count())
			# Set max number of processes in Windows
			if os.name == 'nt':
				if config['processes'] > 61:
					config['processes'] = 61
					print("---- Set new max number of processes: ", config['processes'])
			else:
				print("---- Current number of processes: ", config['processes'])	

			# List of inputs for minimization routine
			fit_inputs = [	(data,
							parameters['name'].to_list(), parameters['value'].to_list(), parameters['free'].to_list(), parameters['prior'].to_list(), parameters['low_l'].to_list(), parameters['high_l'].to_list(),
							function, config['temperature-init'], config['temperature-gain'], config['target-X2'], config['state'], prt_progress)]			
				
			# Create a pool of processes with the number of processors
			pool = multiprocessing.Pool(processes=config['processes'])

			# Calculate the total number of iterations for tqdm
			total_iterations = len(fit_inputs) * config['iterations']

			# Use tqdm to display a progress bar
			with tqdm(total=total_iterations) as pbar:
				# Initialize list for starage
				results_list = []
				# Define a function to update progress and store results
				def update_progress(result):
					results_list.append(result)
					pbar.update()
				# Use imap to get run the minimization routine asynchronously and update progress
				results = pool.imap_unordered(SimAnnealing, fit_inputs * config['iterations'])
				for result in results:
					update_progress(result)

			# Close the pool to free up resources
			pool.close()
			pool.join()

			# Collect stored results
			for result in results_list:
				collection.append(result[0])    
				collection_X2.append(result[1])

	######################################################################

	#### Build statistics and save results
	if not config['plot-only']:

		#------ Prepare the metadata header to be attached to the following output data files.
		#------ This is a summary of the configuration file.
		def add_metadata(file):
			file.writelines("# -----------------------------------\n")
			file.writelines("# SAS_MoCa Version "+ str(__version__)+"\n")
			file.writelines("# "+ str(localtime)+"\n")
			file.writelines("# -----------------------------------\n")
			for c in range(len(config)):
				file.writelines( str(list(config.keys())[c]) + ":\t" + str(list(config.values())[c]) + "\n" )		
			file.writelines("# -----------------------------------\n")

		#------ Fill collected results into the dataframe "results_collection"
		dict = {}
		for p in range(len(parameters['name'])):
			tmp = []
			for clt in (collection):
				tmp.append(clt.iloc[p,1])
			dict[p] = tmp
		del tmp
		results_collection=pd.DataFrame.from_dict(dict,orient='index').transpose()
		results_collection.columns = parameters['name']

		#------ get list of names of free parameters
		names=[]
		for i,p in enumerate(parameters['name']):
			if parameters['free'].iloc[i]:
				names.append(parameters['name'].iloc[i])   
		
		#------ Calculate Pearson correlation coefficient matrix for inters>=10
		if config['iterations'] >= 10:
			pearson_correlation = np.empty([len(names),len(names)])
			# Create Pearson correlation coefficient matrix
			for i_n, n in enumerate(names,0):
				for i_m, m in enumerate(names,0):
					pearson_correlation[i_n][i_m], _  =  stats.pearsonr(results_collection[n], results_collection[m])
			pearson_correlation = pd.DataFrame(pearson_correlation, index=names, columns=names)
			# Save Pearson correlation coefficient matrix
			with open("./"+config['save-folder']+"/correlations.dat", 'w') as fl:
				add_metadata(fl)
				fl.writelines("# Pearson correlation coefficient matrix,\n")
				fl.writelines("# see scipy.stats.pearsonr() docs for details.\n")
				fl.writelines("# -----------------------------------\n")
			pearson_correlation.to_csv("./"+config['save-folder']+"/correlations.dat", sep='\t', na_rep='-', header=True, index=True, mode='a')

		#------ Add X2 column to "results_collection"
		results_collection['X2'] = collection_X2

		#------ Get median values and standard deviation for each parameter
		median = []
		MAD = []
		for p in range(len(parameters['name'])):
			tmp = []
			for clt in (collection):
				tmp.append(clt.iloc[p,1])
			median.append( np.median(np.array(tmp)) )
			MAD.append( stats.median_abs_deviation(np.array(tmp), nan_policy='propagate', scale=1.0) )
		del tmp
		parameters['median'] = median
		parameters['MAD'] = MAD

		#------ print parameters recap and results
		print("\n----- Results recap -----\n")
		print(parameters)

		#------ initilize object to calculate final intensity and extra parameters
		par_plot = []
		for p in range(len(parameters['name'])):
			par_plot.append(parameters['median'].iloc[p])
		res_function = function(data[:,0], par_plot, config['state'])	
		#------ plot data and fitted model
		I_plot = res_function.intensity()
		PlotData(config['qrange'], config['save-folder']).plot_fit( data, I_plot )	
		
		if LUV:
			#------ calculate extra parameters and real space description
			SDP_matrix, n_W, D_pp, D_B, D_C, A_L = res_function.SDP_profile()
			#------ plot SDP and SLD profiles
			SDP = PlotSDP_profile(SDP_matrix, D_B, D_C, config['save-folder'])
			SDP.plot_sdp()
			#------ save SDP and SLD profiles
			with open("./"+config['save-folder']+"/plot_SDP-SLD.dat", 'w') as fl:
				add_metadata(fl)
				fl.writelines("# Probabilities of finding a quasimolecular group at a position z (Angstrom) (distance from the bilayer center).\n")
				fl.writelines("# SLD column: SLD contrast (difference from suspension medium) as a function of z; unit Angstrom^-2.\n")
				fl.writelines("# -----------------------------------\n")
			SDP_matrix.to_csv("./"+config['save-folder']+"/plot_SDP-SLD.dat", sep='\t', na_rep='-', header=True, index=False, mode='a', index_label=False)

		#------ compute equivalent X^2 from the set of median results
		N_Free = 0
		for v in range(len(parameters['value'])):
			if parameters.iloc[v,2]: N_Free+=1
		X2_median = X2function(data[:,0], data[:,1], I_plot, data[:,2], N_Free)

		#------ plot and save statistics 
		Plots=PlotStat(results_collection, parameters)
		#------ parameters and chi-squared histograms
		if config['iterations'] >= 10: 
			Plots.histograms("./"+config['save-folder']+"/plot_histograms.png")
			Plots.histogram_X2("./"+config['save-folder']+"/plot_histogram_X2.png")
			Plots.correlations(pearson_correlation, "./"+config['save-folder']+"/plot_correlations.png")

		#------------------ Saving results
		if LUV:
			# Calculate extra physical quantities and add it to parameters dataframe
			rel_err_dD_C = float(parameters.loc[parameters['name']=='dD_C','MAD']/parameters.loc[parameters['name']=='dD_C','median'])
			new_row = pd.DataFrame({"name": ["A_L*", "D_B*", "D_pp*", "n_W*"],
									"median": [A_L, D_B, D_pp, n_W],
									"MAD": [A_L*rel_err_dD_C, D_B*rel_err_dD_C, "-", "-"]})
			parameters = pd.concat([parameters, new_row], ignore_index=True)
		
		#------ Recap of results. It includes parameter initialization
		with open("./"+config['save-folder']+"/results_recap.dat", 'w') as fl:
			add_metadata(fl)
			fl.writelines("# Table of results (median and MAD) along with parameter initialization.\n")
			fl.writelines("# * Calculated values\n")
			fl.writelines("# X^2 from median values: "+str(X2_median)+"\n")
			fl.writelines("# -----------------------------------\n")
		parameters.to_csv("./"+config['save-folder']+"/results_recap.dat", sep='\t', na_rep='-', header=True, index=False, mode='a', index_label=False)
		
		#------ Collection of the result parameters (medians only) for each iteration
		with open("./"+config['save-folder']+"/iterations_collection.dat", 'w') as fl:
			add_metadata(fl)
			fl.writelines("# Collection of the result parameters (medians only) for each iteration.\n")
			fl.writelines("# -----------------------------------\n")
		results_collection.to_csv("./"+config['save-folder']+"/iterations_collection.dat", sep='\t', na_rep='-', header=True, index=True, mode='a')

	#------ Save data and model scattering intensity ( 0: q, 1: I_data, 2: I_err, 3: I_fit)
	results_intensity= pd.DataFrame(np.column_stack((data,I_plot)), columns=['q', 'I_data', 'I_err', 'I_fit'])
	results_intensity=results_intensity.set_index('q')
	if not config['plot-only']:
		with open("./"+config['save-folder']+"/plot_intensity.dat", 'w') as fl:
			add_metadata(fl)
	results_intensity.to_csv("./"+config['save-folder']+"/plot_intensity.dat", sep='\t', na_rep='-', header=True, index=True, mode='a')
