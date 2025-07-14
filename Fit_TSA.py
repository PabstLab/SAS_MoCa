#!/usr/bin/python

import sys
import numpy as np
import pandas as pd
from scipy import stats

from toremove_Get_input import Load_input
from toremove_Get_data import Load_data
from Model_list import ChooseFunction
import TSA_algorithm
import toremove_Show_output

############################################################
############################################################
############################################################

def X2function(Q,IDATA,I,ERR,N):
	'X^2 calculation'
	X2 = np.sum(  ( ( IDATA-I ) / ERR )**2 )
	return X2/(Q.shape[0]-(N-1))

############################################################

# Tune pandas dataframe printing options
pd.options.display.max_colwidth = 100
# Print X^2 progression for each iteration
prt_progress = 1

############################################################ Read configuration options & model parameters

############### Read File & Options
# Load configuration and parameter file
input = Load_input(sys.argv)
# Load configuration details
config = input.get_config()
# Load parameter matrix
parameters = input.get_parameters()

############### Read model function and plotting options
"-- Reading Function"
(function, pltOptions) = ChooseFunction( config['function'] )

############################################################ Read & Select Data
data_init = Load_data( config['datafile'], config['qrange'], config['bing'] )
data = data_init.convert( config['err_mul'] )

############################################################ Fit Routine
collection = []
collection_X2 = []

#### Plot only
#------ show data and simulated model
if config['plot'] == 1:

	par_plot = []

	for v in range(len(parameters['value'])):
		par_plot.append(parameters.iloc[v,1])

	I_plot = function(data[:,0],par_plot,config['state']).intensity()

	plot_init = toremove_Show_output.PlotData( config['label'], config['qrange'], config['save'] )	
	plot_init.plot_fit( data, I_plot )

#### Fit data
#------ run the minimization routine
#### One process only: 
#----------- serial computation of all iterations
#----------- display minimization stats and progress
elif config['plot'] == 0:
	for it in range(config['iterations']) :
		print("\n---- Iteration N.", it+1," -----\n")
			
		# List of inputs for minimization algorighm
		fit_inputs = [	data,
				 		parameters['name'].to_list(), parameters['value'].to_list(), parameters['fix'].to_list(), parameters['prior'].to_list(), parameters['low_l'].to_list(), parameters['high_l'].to_list(),
						function, config['temp'], config['thermo'], config['target'], config['state'], prt_progress]			
		# Run the minimization routine		
		(par_res, X2_min) = TSA_algorithm.SimAnnealing( fit_inputs )	

		# Collect stored results
		collection.append(par_res)	
		collection_X2.append(X2_min)	

######################################################################

#### Build statistics and save results
if config['plot'] == 0:

	#------ Fill collected results into a ordered dataframe 
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
		if parameters['fix'].iloc[i] != "f":
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
		pearson_correlation.to_csv("./"+config['save']+"/Results_pearsonr.dat", sep='\t', na_rep='-', header=True, index=True, mode='w')

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
	toremove_Show_output.PlotData( config['label'], config['qrange'], config['save'] ).plot_fit( data, I_plot )	

	#------ plot and save histograms
	if config['iterations'] >= 10:
		toremove_Show_output.PlotStat(results_collection, parameters, "./"+config['save']+"/Plot_histograms.png").histograms()

	#------ compute equivalent X^2 from the set of mean results
	N_Free = 0
	for v in range(len(parameters['value'])):
		if parameters.iloc[v,2]!="f" : N_Free+=1
	X2_mean = X2function(data[:,0],data[:,1],I_plot,data[:,2],N_Free)
	X2_mean_to_print = np.empty(1,dtype=float)
	X2_mean_to_print[0] = X2_mean		

	#------ plot and save X^2-histogram
	if config['iterations'] >= 2: toremove_Show_output.PlotStat(results_collection, parameters, "./"+config['save']+"/Plot_histogram_X2.png").histogram_X2(X2_mean)

	#------ save results

	# Global X^2 (from mean values)
	np.savetxt("./"+config['save']+"/Results_X2_mean.dat", X2_mean_to_print, header='X2 from mean values')		
	# Parameter set and info: from start values to results
	parameters.to_csv("./"+config['save']+"/Results_recap.dat", sep='\t', na_rep='-', header=True, index=False, mode='w', index_label=False)
	# Collection of results from each instance
	results_collection.to_csv("./"+config['save']+"/Results_collection.dat", sep='\t', na_rep='-', header=True, index=True, mode='w')

# Metadata (config dictionary)
with open("./"+config['save']+"/Results_metadata.dat", 'w') as fl:
	for c in range(len(config)):
		fl.writelines( str(list(config.keys())[c]) + "\t" + str(list(config.values())[c]) + "\n" )

# Save data and model scattering intensity ( 0: q, 1: I_data, 2: I_err, 3: I_fit)
np.savetxt("./"+config['save']+"/Results_intensity.dat", np.column_stack((data,I_plot)), header='0: q, 1: I_data, 2: I_err, 3: I_fit')	