#!/usr/bin/python

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

###########################################################################################
###########################################################################################

def reshape (base, c, r):
	'reshape array into matrix'
	new = []
	i=0    
	last = 0
	while i < r :
		row = []
		for ele in base[last:last+c]:
			row.append(ele)
		new.append(row)
		i+=1     
		last = i*c
	return new    


###########################################################################################
###########################################################################################

class PlotData:
	'Plot Results'

	def __init__ ( self, qrange, save ):

		self.fig1, self.axes = plt.subplots(2, 1, figsize=(9.0, 12), facecolor='#2E3436')

		bk_color = '#222222' #'#2E3436'
		label_color = '#D3D7CF'
		
		self.axes[0].set_yscale('log', base=10)
		self.axes[0].set_title(save, color=label_color, fontsize=16)

		for g in range(len(self.axes)):
			self.axes[g].set_xscale('log', base=10)

			self.axes[g].grid()		
			self.axes[g].set_facecolor(bk_color)			

			self.axes[g].tick_params(axis='both', which='major', labelsize=14, colors=label_color)
			self.axes[g].tick_params(axis='both', which='minor', labelsize=8, colors=label_color)

			self.axes[g].set_xlim([qrange[0],qrange[1]])
			self.axes[g].set_xlabel('$q$ ($\\AA^{-1}$)', fontsize=14, color=label_color)
			
		self.axes[0].set_ylabel('$I(q)$ (mm$^{-1}$)', fontsize=14, color=label_color)
		self.axes[1].set_ylabel('relative deviation $|I(q)-f(q)|/I(q)$', fontsize=14, color=label_color)
	
		self.save_folder = "./"+str(save)+"/Plot.eps"
		self.save_folder2 = "./"+str(save)+"/Plot.png"

###########################################################################################

	def plot_data ( self, data ):	

		point_color = ['orange', 'magenta', 'lightblue']
		point_edge = '#551000'		

		self.axes[0].errorbar(data[0][:,0], data[0][:,1], yerr=data[0][:,2], fmt='o', color=point_color[0], markeredgecolor=point_edge, markersize=5, linewidth=3.0, label='SAXS data', zorder=0)		

		self.fig1.tight_layout()

		plt.show()

###########################################################################################

	def plot_fit ( self, data, I_plot ):	

		point_color = ['orange', 'magenta', 'blue', 'green', 'orange']
		point_edge = '#551000'		
		line_color = ['red', 'violet', 'lightblue', 'lightgreen', 'red']

		self.axes[0].errorbar(data[:,0], data[:,1], yerr=data[:,2], fmt='o', color=point_color[0], markeredgecolor=point_edge, markersize=5, linewidth=1.0, label='SAXS data', zorder=0)		
		self.axes[0].plot(data[:,0], I_plot, linewidth=2.0, color=line_color[0], ls='-', label='Best fit', zorder=1)

		self.axes[1].axhline(y=0, linewidth=2)		
		self.axes[1].scatter(data[:,0], data[:,2]/data[:,1], marker='o', color=point_color[0], edgecolor=point_edge, s=18, linewidth=1.2, label='Exp. error', zorder=0)		
		self.axes[1].plot(data[:,0], np.abs(I_plot-data[:,1])/data[:,1], linewidth=2.0, color=line_color[0], ls='-', label='Best fit', zorder=1)

		self.axes[0].set_ylim([0.1*np.min(data[:,1]),10*np.max(data[:,1])])

		for g in range(len(self.axes)):
			self.axes[g].legend()
		self.fig1.tight_layout()

		#plt.savefig(self.save_folder, transparent=False, dpi=150, format='eps',
        #metadata={'Creator': "SAS-MoCa_proteo"}, 
		#bbox_inches='tight', facecolor='auto', edgecolor='auto')

		plt.savefig(self.save_folder2, transparent=False, dpi=150, format='png',
        metadata={'Creator': "SAS-MoCa_proteo"}, 
		bbox_inches='tight', facecolor='auto', edgecolor='auto')

		#plt.show()
		plt.close()

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

class PlotStat:
	'Plot Statistics'

	def __init__ ( self , results, parameters, save ):

		# restrict the list to free parameteres only
		self.name=[]
		self.start=[]
		self.mean=[]
		self.stdev=[]
		self.low_l=[]
		self.high_l=[]
		self.prior=[]
		for i, p in enumerate(parameters['name']):
			if parameters['free'].iloc[i]:
				self.name.append(parameters['name'].iloc[i])        
				self.start.append(parameters['value'].iloc[i])       
				self.prior.append(parameters['prior'].iloc[i])    
				self.low_l.append(parameters['low_l'].iloc[i])        
				self.high_l.append(parameters['high_l'].iloc[i])      
				self.mean.append(parameters['mean'].iloc[i])        
				self.stdev.append(parameters['stdev'].iloc[i])        

		self.results = results
		self.save = save

###########################################################################################

	def histograms(self):
	
		# reshape the free parametr list to a i,j-matrix to plot
		cols=5
		rows=int(len(self.name)/cols)+1
		left = rows*cols-len(self.name)

		bk_color = '#222222' #'#2E3436'
		label_color = '#D3D7CF'
		
		hist_color = '#BB86FC'
		hist_edge = '#000000' 	
		kde_color = '#F0E4FD'
		prior_color = '#F8F042'
		res_color = '#FF0000'

		name_reshaped = reshape(self.name, cols, rows)
		start_reshaped = reshape(self.start, cols, rows)
		prior_reshaped = reshape(self.prior, cols, rows)
		low_l_reshaped = reshape(self.low_l, cols, rows)
		high_l_reshaped = reshape(self.high_l, cols, rows)
		mean_reshaped = reshape(self.mean, cols, rows)
		stdev_reshaped = reshape(self.stdev, cols, rows)

		# initialize plot
		fig, axs = plt.subplots(rows, cols, figsize=(4.*cols, 4.*rows), facecolor='#2E3436')

		for i in range(len(name_reshaped)):
			if i!=rows-1 :
				end=len(name_reshaped[0])
			else: 
				end=len(name_reshaped[0])-left
			for j in range(end):

				# plot options	
				axs[i][j].set_facecolor(bk_color)	

				axs[i][j].tick_params(axis='both', which='major', labelsize=14, colors=label_color)
				axs[i][j].tick_params(axis='both', which='minor', labelsize=8, colors=label_color)		

				axs[i][j].set_xlim(low_l_reshaped[i][j],high_l_reshaped[i][j])
				axs[i][j].set_xlabel(name_reshaped[i][j], color=label_color)
				
				# set x-range within boundaries
				x = np.arange(low_l_reshaped[i][j],high_l_reshaped[i][j],(high_l_reshaped[i][j]-low_l_reshaped[i][j])/100)
				# get the y-array of the calculated histogram
				y, _, _ = plt.hist(self.results[name_reshaped[i][j]], density=True)

				# generate and plot histogram
				axs[i][j].hist(self.results[name_reshaped[i][j]], density=True, facecolor=hist_color, edgecolor=hist_edge, linewidth=2.5, zorder=0)
				# generate and plot histogram kernel density
				kde = stats.gaussian_kde(self.results[name_reshaped[i][j]])
				axs[i][j].plot(x,kde(x), ls='-', lw=3.5, color=kde_color, alpha=1.0, zorder=1)
				
				# generate and plot priors
				if prior_reshaped[i][j]!=0 : 
					prior_pdf = y.max()*np.exp( -(x-start_reshaped[i][j])**2 / (2*(start_reshaped[i][j]*prior_reshaped[i][j])**2) )
					axs[i][j].plot(x, prior_pdf, ls='-', lw=3.5, color=prior_color, alpha=1.0, zorder=2)

				# generate mean and standard deviation			
				axs[i][j].errorbar( x=mean_reshaped[i][j], y=y.max()/2, xerr=stdev_reshaped[i][j], marker='D', ms=8, markerfacecolor=res_color, markeredgecolor='None', 
					   				ecolor=res_color, elinewidth=3.0, capsize=8, zorder=4)

		fig.tight_layout()

		plt.savefig(self.save, transparent=False, dpi=100, format='png',
       				metadata={'Creator': "SAS-MoCa_proteo"}, 
					bbox_inches='tight', facecolor='auto', edgecolor='auto')

		#plt.show()
		plt.close()

###########################################################################################

	def histogram_X2(self, X2_mean):
	
		# reshape the free parametr list to a i,j-matrix to plot
		bk_color = '#222222'
		label_color = '#D3D7CF'
		
		hist_color = '#BB86FC'
		hist_edge = '#000000' 	
		kde_color = '#F0E4FD'

		# initialize plot
		fig, axs = plt.subplots(1, 1, figsize=(4., 4.), facecolor='#2E3436')

		# plot options	
		axs.set_facecolor(bk_color)	
		axs.tick_params(axis='both', which='major', labelsize=10, colors=label_color)
		axs.tick_params(axis='both', which='minor', labelsize=6, colors=label_color)		
				
		# set x-range within boundaries
		_, x, _ = plt.hist(self.results['X2'].to_list(), density=True)
		x = np.arange(x.min()*0.9,x.max()*1.1,(x.max()*1.1-x.min()*0.9)/100)

		# generate and plot histogram
		axs.hist(self.results['X2'].to_list(), density=True, facecolor=hist_color, edgecolor=hist_edge, linewidth=2.5)#, zorder=0)
		# generate and plot histogram kernel density
		kde = stats.gaussian_kde(self.results['X2'].to_list())
		axs.plot(x,kde(x), ls='-', lw=3.5, color=kde_color, alpha=1.0, zorder=1)
		axs.set_xlabel('$\chi^2$', fontsize=14, color=label_color)

		#axs.axvline(x=X2_mean, ls=':', lw=3.5, color=kde_color, alpha=1.0, zorder=1)

		fig.tight_layout()

		plt.savefig(self.save, transparent=False, dpi=150, format='png',
      				metadata={'Creator': "SAS-MoCa_proteo"}, 
					bbox_inches='tight', facecolor='auto', edgecolor='auto')

		#plt.show()
		plt.close()

###########################################################################################