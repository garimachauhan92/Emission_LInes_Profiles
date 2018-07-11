import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import os
import math as math
import time
from scipy import integrate
from scipy import interpolate 
import random as random
import sys

import pandas as pd

#################################################################################################################################################################


with np.errstate(divide = 'ignore', invalid = 'ignore'):

	exec(open("/home/garima/pleiades_icrar/Lightcones/Codes/light_cones_reading_file.py").read())
	exec(open("/home/garima/pleiades_icrar/Lightcones/Codes/loading_data_file.py").read())






#path_save is in all the functions - use the below addresses to specify what and where


path_central	= 	'/home/garima/pleiades_icrar/Lightcones/Plots/Plots_centrals/Pattern/'
path_sub		=	'/home/garima/pleiades_icrar/Lightcones/Plots/Plots_centrals_sub/Pattern/'
path_20			=	'/home/garima/Desktop/'

#***************************************************************************************************************************************************************


mass_galaxies = mass_galaxies[index_central]


def plot_hist(central,satellite,property_hist):

	central = central[~np.isnan(central)]
	satellite = satellite[~np.isnan(satellite)]

	central = central[central != 0]
	satellite = satellite[satellite != 0]
	
	bins_host = np.linspace(np.log10(min(central)),np.log10(max(central)),num = 50)
	bins_sat = np.linspace(np.log10(min(satellite)),np.log10(max(satellite)),num = 50)
	

	plt.figure(figsize = (12,8))
	plt.hist(np.log10(central), bins = bins_host,   lw = 1.5, ec = 'b', histtype = 'bar', alpha = 0.5)
	#plt.hist(np.log10(satellite), bins = bins_sat,  lw = 1.5, ec = 'r', histtype = 'bar', alpha =0.5)
	#plt.yscale('log')
	#plt.xscale('log')
	#plt.legend(('Central','Satellite'))
	plt.xlabel('$log_{10}$(%s)' %property_hist)
	plt.ylabel('$Num_{galaxies}$')
	plt.title('%s' %(property_hist))
	plt.show()









# Plotting Line-functions ###

def plot_all(velocity,width, name, name_width, path,log) :

	path_save = path

	# Normal - W_50 for all

	plt.figure(figsize = (12,8))

	plt.scatter(velocity, width, s=2, alpha = 0.5, c = mass_galaxies)
	plt.plot(velocity, velocity, linestyle = ':', linewidth=0.2)
	plt.nipy_spectral()
	cbar = plt.colorbar()
	cbar.ax.set_title('$(M_{*}$  $M_{0})$', size=10)

	#plt.xlim(0, 100)
	#plt.ylim(0,100)
	#plt.xlim(0.5*10**2,0.5*10**3)
	#plt.ylim(0.5*10**2,0.5*10**3)

	if log == 1:
		plt.xscale('log')
		plt.yscale('log')
		plt.xlabel('$Velocity_{%s}$ km/s' %(name))
		plt.ylabel('$%s/2$ km/s' %(name_width))
		plt.title('$%s/2$ vs $Velocity_{%s}$' %(name_width,name))
		#plt.savefig(path_save + name + name_width + 'log' + '.png')
		plt.show()
	
	else:		

		plt.xlabel('$Velocity_{%s}$ km/s' %(name))
		plt.ylabel('$%s/2$ km/s' %(name_width))
		plt.title('$%s/2$ vs $Velocity_{%s}$' %(name_width,name))
		#plt.savefig(path_save + name + name_width + '.png')
		plt.show()


def plot_trial(velocity,width, name, name_width, path) :

	path_save = path

	# Normal - W_50 for all

	plt.figure(figsize = (12,8))

	plt.scatter(velocity, width, s=2, alpha = 0.5, c = 'k')
	plt.plot(velocity, velocity, linestyle = ':', linewidth=0.2)
	
	#plt.xlim(min(velocity), max(velocity))
	#plt.ylim(min(width), max(width))
	#plt.xlim(0.5*10**2,0.5*10**3)
	#plt.ylim(0.5*10**2,0.5*10**3)

	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('$Velocity_{%s}$ km/s' %(name))
	plt.ylabel('$%s/2$ km/s' %(name_width))
	plt.title('$%s/2$ vs $Velocity_{%s}$' %(name_width,name))
	plt.savefig(path_save + name + name_width + '.png')
	plt.show()


def plot_mass_BT():

	plt.figure(figsize = (12,8))
	plt.scatter(B_T_central[B_T_central > 0], np.log10(mass_gas[B_T_central > 0]), s=2, alpha = 0.5, c = mass_galaxies[B_T_central > 0])
	plt.nipy_spectral()
	cbar = plt.colorbar()
	cbar.ax.set_title('$M_{*}$ $M_{0}$', size = 10)
	plt.show()






def plotting_BT():

	yo_7 = np.poly1d(np.polyfit(Median_stars_central_7["B_T"], Median_stars_central_7["Ratio"],2))

	yo_8 = np.poly1d(np.polyfit(Median_stars_central_8["B_T"], Median_stars_central_8["Ratio"],2))

	yo_9 = np.poly1d(np.polyfit(Median_stars_central_9["B_T"], Median_stars_central_9["Ratio"],2))

	yo_10 = np.poly1d(np.polyfit(Median_stars_central_10["B_T"], Median_stars_central_10["Ratio"],2))



	plt.plot(Median_stars_central_7["B_T"], np.log10(yo_7(Median_stars_central_7["B_T"])))
	plt.plot(Median_stars_central_8["B_T"], np.log10(yo_8(Median_stars_central_8["B_T"])))
	plt.plot(Median_stars_central_9["B_T"], np.log10(yo_9(Median_stars_central_9["B_T"])))
	plt.plot(Median_stars_central_10["B_T"], np.log10(yo_10(Median_stars_central_10["B_T"])))
	plt.legend(('$M_{*}$: < 1e7$M_{o}$', '$M_{*}$: 1e8:1e9$M_{o}$', '$M_{*}$: 7e9:3e10$M_{o}$', '$M_{*}$: > 3e10 $M_{o}$'))
	plt.title('B/T graph')
	plt.xlabel('B/T ratio')
	plt.ylabel('$M_{HI}/M_{*}$')
	plt.show()









def plot_velocity(k):

	name = 'circ_flat_%d.png' %k
	
	plt.figure(figsize = (12,8))
	plt.plot(r_x[k], V_circ_all[k],  '--b', linewidth = 0.5)
	plt.plot(r_x_flat[k], V_flat_all[k],  '--r', linewidth = 0.5)
	
	plt.xlabel('Radius (Mpc) km/s')
	plt.ylabel('$V_{circ}$ km/s')
	plt.title('Comparison of $V_{circ}$ and $V_{flat}$')
	plt.savefig(path_save_comparison + name)
	plt.show()



def plot_comparisons(width_1, width_2, name_1, name_2, path, log):

	path_save = path

	plt.figure(figsize = (12,8))

	plt.scatter(width_1, width_2, s=2, alpha = 0.5, c = 'r' )#c = mass_galaxies, cmap = 'nipy_spectral')
	
	plt.plot(width_1,width_1, linewidth=0.5, linestyle='--', c='b')
	#plt.legend(('W50', 'Galaxies'))

	#cbar = plt.colorbar()
	#cbar.ax.set_title('$(M_{*}$  $M_{0})$', size=10)


	if log == 1:

		plt.xscale('log')
		plt.yscale('log')
		plt.ylabel('%s km/s' %(name_2))
		plt.xlabel('%s km/s' %(name_1))
		plt.title('%s vs %s' %(name_2, name_1))
		plt.savefig(path_save + name_1 + '_' + name_2 + 'log.png')
		plt.show()
	else:
		plt.ylabel('%s km/s' %(name_2))
		plt.xlabel('%s km/s' %(name_1))
		plt.title('%s vs %s' %(name_2, name_1))
		plt.savefig(path_save + name_1 + '_' + name_2 + '.png')
		plt.show()



def plot_comparison_lines(k):

	name = 'Number_%d.png' %k

	plt.figure(figsize = (12,8))

	edge = plt.plot(v_x[k],s_flux[k] , '--b', label = "Edge-on")
	random = plt.plot(v_x[k+1],s_flux[k+1] , '--m', label = "Random Orientation")

	plt.axhline(y = s_50[k], color = 'r', linestyle = '-', linewidth = 0.5)
	#plt.axhline(y = s_20[k], color = 'g', linestyle = '-', linewidth = 0.5)
	
	plt.axhline(y = s_50[k+1], color = 'r', linestyle = '--', linewidth = 0.5)
	#plt.axhline(y = s_20_r[k], color = 'g', linestyle = '--', linewidth = 0.5)
	
	plt.annotate('$W_{50}^{Random}$', xy = (0,s_50[k+1]), fontsize = 10)
	#plt.annotate('$W_{20}^{Random}$', xy = (0,s_20_r[k]), fontsize = 10)
	plt.annotate('$W_{50}^{Edge}$', xy = (0,s_50[k]), fontsize = 10)
	#plt.annotate('$W_{20}^{Edge}$', xy = (0,s_20[k]), fontsize = 10)

	plt.legend(('Edge-on','Random orientation'), loc = "best")#, ['Edge-on', 'Random'],loc = 2)
	
	plt.title('H1 emission line')
	plt.xlabel('$V_{obs}$ - km/s')
	plt.ylabel('$ \Psi_{H1} $')
	#plt.savefig(path_save_lines + name)
	plt.show()

	print(M_stars_disk[k])
	print(W_50[k])
	print(np.arcsin(theta_r[k])*180/np.pi)

	#plt.close()



def plot_emission_line_single(k):

	plt.figure(figsize = (12,8))

	plt.plot(v_x[k],s_flux[k] , '--b', label = "Edge-on")
	#random = plt.plot(v_x_r[k],s_final_r[k] , '--m', label = "Random Orientation")

	plt.axhline(y = s_50[k], color = 'r', linestyle = '-', linewidth = 0.5)
	plt.axhline(y = s_20[k], color = 'g', linestyle = '-', linewidth = 0.5)
	plt.axhline(y = s_high[k], color = 'm', linestyle = '-', linewidth = 0.5)
	#plt.axhline(y = s_50_r[k], color = 'r', linestyle = '--', linewidth = 0.5)
	#plt.axhline(y = s_20_r[k], color = 'g', linestyle = '--', linewidth = 0.5)
	
	#plt.annotate('$W_{50}^{Random}$', xy = (0,s_50_r[k]), fontsize = 10)
	#plt.annotate('$W_{20}^{Random}$', xy = (0,s_20_r[k]), fontsize = 10)
	plt.annotate('$W_{50}$', xy = (0,s_50[k]), fontsize = 10)
	plt.annotate('$W_{20}$', xy = (0,s_20[k]), fontsize = 10)
	plt.annotate('$W_{peak}$', xy = (0,s_high[k]), fontsize = 10)
	#plt.legend(('Edge-on','Random orientation'), loc = "best")#, ['Edge-on', 'Random'],loc = 2)
	
	plt.title('H1 emission line')
	plt.xlabel('$V_{obs}$ - km/s')
	plt.ylabel('$ \Psi_{H1} $')
	#plt.savefig(path_save_lines + name)
	plt.show()
	#plt.close()





def plot_BT_thing(log):

	if log == 1:
		f,axarr = plt.subplots(4,4,sharex='all',sharey='all')
		axarr[3,2].set_ylim(0.00001,5)
		axarr[3,2].set_xlim(-400,400)
		axarr[3,2].set_yscale("log")

	else:
		f,axarr = plt.subplots(4,4,sharex='col',sharey='row')
		

	
	#axarr[0,0].plot(v_x_7[0],s_final_7[0]*np.array(flux_c_7[0]),'--b')
	axarr[0,0].plot(v_x_s_7[0],s_final_s_7[0]*np.array(flux_s_7[0]),'--r')
	axarr[0,0].set_ylabel('B/T - 0.7 < ')

	axarr[0,1].plot(v_x_7[1],s_final_7[1]*np.array(flux_c_7[1]),'--b')
	axarr[0,1].plot(v_x_s_7[1],s_final_s_7[1]*np.array(flux_s_7[1]),'--r')
	
	axarr[0,2].plot(v_x_7[2],s_final_7[2]*np.array(flux_c_7[2]),'--b')
	axarr[0,2].plot(v_x_s_7[2],s_final_s_7[2]*np.array(flux_s_7[2]),'--r')

	axarr[0,3].plot(v_x_7[3],s_final_7[3]*np.array(flux_c_7[3]),'--b')
	axarr[0,3].plot(v_x_s_7[3],s_final_s_7[3]*np.array(flux_s_7[3]),'--r')




	axarr[1,0].plot(v_x_5[0],s_final_5[0]*np.array(flux_c_5[0]),'--b')
	axarr[1,0].plot(v_x_s_5[0],s_final_s_5[0]*np.array(flux_s_5[0]),'--r')
	axarr[1,0].set_ylabel('B/T - 0.50:0.70 ')

	axarr[1,1].plot(v_x_5[1],s_final_5[1]*np.array(flux_c_5[1]),'--b')
	axarr[1,1].plot(v_x_s_5[1],s_final_s_5[1]*np.array(flux_s_5[1]),'--r')
	
	axarr[1,2].plot(v_x_5[2],s_final_5[2]*np.array(flux_c_5[2]),'--b')
	axarr[1,2].plot(v_x_s_5[2],s_final_s_5[2]*np.array(flux_s_5[2]),'--r')

	axarr[1,3].plot(v_x_5[3],s_final_5[3]*np.array(flux_c_5[3]),'--b')
	axarr[1,3].plot(v_x_s_5[3],s_final_s_5[3]*np.array(flux_s_5[3]),'--r')




	axarr[2,0].plot(v_x_2[0],s_final_2[0]*np.array(flux_c_2[0]),'--b')
	axarr[2,0].plot(v_x_s_2[0],s_final_s_2[0]*np.array(flux_s_2[0]),'--r')
	axarr[2,0].set_ylabel('B/T - 0.20:0.35 ')

	axarr[2,1].plot(v_x_2[1],s_final_2[1]*np.array(flux_c_2[1]),'--b')
	axarr[2,1].plot(v_x_s_2[1],s_final_s_2[1]*np.array(flux_s_2[1]),'--r')
	
	axarr[2,2].plot(v_x_2[2],s_final_2[2]*np.array(flux_c_2[2]),'--b')
	axarr[2,2].plot(v_x_s_2[2],s_final_s_2[2]*np.array(flux_s_2[2]),'--r')

	axarr[2,3].plot(v_x_2[3],s_final_2[3]*np.array(flux_c_2[3]),'--b')
	axarr[2,3].plot(v_x_s_2[3],s_final_s_2[3]*np.array(flux_s_2[3]),'--r')




	axarr[3,0].plot(v_x_0[0],s_final_0[0]*np.array(flux_c_0[0]),'--b')
	axarr[3,0].plot(v_x_s_0[0],s_final_s_0[0]*np.array(flux_s_0[0]),'--r')
	axarr[3,0].set_xlabel('Mass - 1e6:1e7')
	axarr[3,0].set_ylabel('B/T - 0:0.15')

	axarr[3,1].plot(v_x_0[1],s_final_0[1]*np.array(flux_c_0[1]),'--b')
	axarr[3,1].plot(v_x_s_0[1],s_final_s_0[1]*np.array(flux_s_0[1]),'--r')
	axarr[3,1].set_xlabel('Mass - 1e8:1e9')

	axarr[3,2].plot(v_x_0[2],s_final_0[2]*np.array(flux_c_0[2]),'--b')
	axarr[3,2].plot(v_x_s_0[2],s_final_s_0[2]*np.array(flux_s_0[2]),'--r')
	axarr[3,2].set_xlabel('Mass - 7e9:3e10')


	#axarr[3,3].plot(v_x[14886],s_flux[14886],'--b')
	#axarr[3,3].plot(v_x_s[14886],s_flux_s[14886],'--r')
	axarr[3,3].set_xlabel('Mass - 4e10 < ')

	

	plt.show()
	




def plot_BT_satellite():

	f,axarr = plt.subplots(4,4,sharex='col',sharey='row')
	
	#axarr[0,0].plot(v_x_7[0],s_final_7[0]*np.array(flux_c_7[0]),'--b')
	axarr[0,0].plot(v_x_s_7[0],s_final_s_7[0]*np.array(flux_s_7[0]),'--r')
	axarr[0,0].set_ylabel('B/T - 0.7 < ')

	#axarr[0,1].plot(v_x_7[1],s_final_7[1]*np.array(flux_c_7[1]),'--b')
	axarr[0,1].plot(v_x_s_7[1],s_final_s_7[1]*np.array(flux_s_7[1]),'--r')
	
	#axarr[0,2].plot(v_x_7[2],s_final_7[2]*np.array(flux_c_7[2]),'--b')
	axarr[0,2].plot(v_x_s_7[2],s_final_s_7[2]*np.array(flux_s_7[2]),'--r')

	#axarr[0,3].plot(v_x_7[3],s_final_7[3]*np.array(flux_c_7[3]),'--b')
	axarr[0,3].plot(v_x_s_7[3],s_final_s_7[3]*np.array(flux_s_7[3]),'--r')




	#axarr[1,0].plot(v_x_5[0],s_final_5[0]*np.array(flux_c_5[0]),'--b')
	axarr[1,0].plot(v_x_s_5[0],s_final_s_5[0]*np.array(flux_s_5[0]),'--r')
	axarr[1,0].set_ylabel('B/T - 0.50:0.70 ')

	#axarr[1,1].plot(v_x_5[1],s_final_5[1]*np.array(flux_c_5[1]),'--b')
	axarr[1,1].plot(v_x_s_5[1],s_final_s_5[1]*np.array(flux_s_5[1]),'--r')
	
	#axarr[1,2].plot(v_x_5[2],s_final_5[2]*np.array(flux_c_5[2]),'--b')
	axarr[1,2].plot(v_x_s_5[2],s_final_s_5[2]*np.array(flux_s_5[2]),'--r')

	#axarr[1,3].plot(v_x_5[3],s_final_5[3]*np.array(flux_c_5[3]),'--b')
	axarr[1,3].plot(v_x_s_5[3],s_final_s_5[3]*np.array(flux_s_5[3]),'--r')




	#axarr[2,0].plot(v_x_2[0],s_final_2[0]*np.array(flux_c_2[0]),'--b')
	axarr[2,0].plot(v_x_s_2[0],s_final_s_2[0]*np.array(flux_s_2[0]),'--r')
	axarr[2,0].set_ylabel('B/T - 0.20:0.35 ')

	#axarr[2,1].plot(v_x_2[1],s_final_2[1]*np.array(flux_c_2[1]),'--b')
	axarr[2,1].plot(v_x_s_2[1],s_final_s_2[1]*np.array(flux_s_2[1]),'--r')
	
	#axarr[2,2].plot(v_x_2[2],s_final_2[2]*np.array(flux_c_2[2]),'--b')
	axarr[2,2].plot(v_x_s_2[2],s_final_s_2[2]*np.array(flux_s_2[2]),'--r')

	#axarr[2,3].plot(v_x_2[3],s_final_2[3]*np.array(flux_c_2[3]),'--b')
	axarr[2,3].plot(v_x_s_2[3],s_final_s_2[3]*np.array(flux_s_2[3]),'--r')




	#axarr[3,0].plot(v_x_0[0],s_final_0[0]*np.array(flux_c_0[0]),'--b')
	axarr[3,0].plot(v_x_s_0[0],s_final_s_0[0]*np.array(flux_s_0[0]),'--r')
	axarr[3,0].set_xlabel('Mass - 1e6:1e7')
	axarr[3,0].set_ylabel('B/T - 0:0.15')

	#axarr[3,1].plot(v_x_0[1],s_final_0[1]*np.array(flux_c_0[1]),'--b')
	axarr[3,1].plot(v_x_s_0[1],s_final_s_0[1]*np.array(flux_s_0[1]),'--r')
	axarr[3,1].set_xlabel('Mass - 1e8:1e9')

	#axarr[3,2].plot(v_x_0[2],s_final_0[2]*np.array(flux_c_0[2]),'--b')
	axarr[3,2].plot(v_x_s_0[2],s_final_s_0[2]*np.array(flux_s_0[2]),'--r')
	axarr[3,2].set_xlabel('Mass - 7e9:3e10')


	#axarr[3,3].plot(v_x[14886],s_flux[14886],'--b')
	#axarr[3,3].plot(v_x_s[14886],s_flux_s[14886],'--r')
	axarr[3,3].set_xlabel('Mass - 4e10 < ')

	axarr[3,2].set_ylim(0.00001,5)
	axarr[3,2].set_xlim(-400,400)
	axarr[3,2].set_yscale("log")


	plt.show()
	

