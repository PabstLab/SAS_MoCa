#!/usr/bin/python

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import multiprocessing
from tqdm import tqdm

from in_out import (Load_input, Load_data,
					PlotData, PlotStat)

from Model_list import ChooseFunction
import TSA_algorithm


############################################################
############################################################
############################################################

def X2function(Q,IDATA,I,ERR,N):
	'X^2 calculation'
	X2 = np.sum(  ( ( IDATA-I ) / ERR )**2 )
	return X2/(Q.shape[0]-(N-1))

############################################################
############################################################
############################################################

# Protect entry point for multiprocessing
if __name__ == "__main__":

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
	(function, _) = ChooseFunction( config['model'] )

	# Initialize control to print progression for each iteration
	if config['processes'] == 1 :	prt_progress = 1
	else :							prt_progress = 0

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

		print('############################', config['state'])
		I_plot = function(data[:,0],par_plot,config['state']).intensity()

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
				(par_res, X2_min) = TSA_algorithm.SimAnnealing( fit_inputs )	

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
							function, config['temperature-init'], config['temperature-gain'], config['target-X2'], config['state'], prt_progress)	]			
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
				results = pool.imap_unordered(TSA_algorithm.SimAnnealing, fit_inputs * config['iterations'])
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
			print("\n## Pearson correlation coefficients\n",pearson_correlation)
			pearson_correlation.to_csv("./"+config['save-folder']+"/Results_pearsonr.dat", sep='\t', na_rep='-', header=True, index=True, mode='w')

		#------ Add X2 column to "results_collection"
		results_collection['X2'] = collection_X2

		#------ Get mean values and standard deviation for each parameter
		mean = []
		stdv = []
		for p in range(len(parameters['name'])):
			tmp = []
			for clt in (collection):
				tmp.append(clt.iloc[p,1])
			mean.append( np.mean(np.array(tmp)) )
			stdv.append( np.std(np.array(tmp), ddof=1) )
		del tmp
		parameters['mean'] = mean
		parameters['stdev'] = stdv

		#------ print parameters recap and results
		print("\n----- Results -----\n")
		print(parameters)

		#------ plot data and fitted model
		par_plot = []
		for p in range(len(parameters['name'])):
			par_plot.append(parameters['mean'].iloc[p])
		I_plot = function(data[:,0],par_plot,config['state']).intensity()
		PlotData(config['qrange'], config['save-folder']).plot_fit( data, I_plot )	

		#------ plot and save histograms
		if config['iterations'] >= 10:
			PlotStat(results_collection, parameters, "./"+config['save-folder']+"/Plot_histograms.png").histograms()

		#------ compute equivalent X^2 from the set of mean results
		N_Free = 0
		for v in range(len(parameters['value'])):
			if parameters.iloc[v,2]!="f" : N_Free+=1
		X2_mean = X2function(data[:,0],data[:,1],I_plot,data[:,2],N_Free)
		X2_mean_to_print = np.empty(1,dtype=float)
		X2_mean_to_print[0] = X2_mean		
		
		#------ plot and save X^2-histogram
		if config['iterations'] >= 2: PlotStat(results_collection, parameters, "./"+config['save-folder']+"/Plot_histogram_X2.png").histogram_X2(X2_mean)
		
		#------ save results
		
		# Global X^2 (from mean values)
		np.savetxt("./"+config['save-folder']+"/Results_X2_mean.dat", X2_mean_to_print, header='X2 from mean values')		
		# Parameter set and info: from start values to results
		parameters.to_csv("./"+config['save-folder']+"/Results_recap.dat", sep='\t', na_rep='-', header=True, index=False, mode='w', index_label=False)
		# Collection of results from each instance
		results_collection.to_csv("./"+config['save-folder']+"/Results_collection.dat", sep='\t', na_rep='-', header=True, index=True, mode='w')

	# Metadata (config dictionary)
	with open("./"+config['save-folder']+"/Results_metadata.dat", 'w') as fl:
		for c in range(len(config)):
			fl.writelines( str(list(config.keys())[c]) + "\t" + str(list(config.values())[c]) + "\n" )

	# Save data and model scattering intensity ( 0: q, 1: I_data, 2: I_err, 3: I_fit)
	np.savetxt("./"+config['save-folder']+"/Results_intensity.dat", np.column_stack((data,I_plot)), header='0: q, 1: I_data, 2: I_err, 3: I_fit')	