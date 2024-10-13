import pycali
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

# AJF initialize parsers
parser = argparse.ArgumentParser(description='uses PyCali to intercalibrate lightcurves')
parser.add_argument('openpath', metavar = 'openpath', help = 'path to open cont_combined files from, probably ../pycali/')  
parser.add_argument('season_num', metavar='season_num', help='total number of seasons, i.e. "5"')
arg = parser.parse_args()
openpath = arg.openpath

# create new directory to store all calibrated files in after pycali
if os.path.isdir(f'{openpath}for_create_cont')==False:
	os.mkdir(f'{openpath}for_create_cont')
	

names = []
for i in range(1, (int(arg.season_num)+1)):
	j = str(i)
	name = 'season'+(j)
	names.append(name) 

"""
for i in range(5):
	print('\n')
print(f'You are going to re-intercalibrate {names} and NOT re-intercalibrate seasons {skip}.')
for i in range(5):
	print('\n')

# copy cont_combined.txt_cali file to new directory if good intercalibration
for i in range(len(skip)):
	os.system(f'cp {openpath}season{skip[i]}/pycali_results/season{skip[i]}_cont_combined.txt_cali {openpath}for_create_cont/season{skip[i]}_cont_combined_new.txt_cali')
	print(f'Since season{skip[i]} does not need to be recalibrated, it has been moved to for_create_cont folder')
for i in range(5):
	print('\n')    
"""
#######################################################
# setup configurations, there are two ways:
# 1) load from a param file
#    cfg = pycali.Config("param.txt")
# 2) direct call setup()
redone = []
not_redone = []
for name in names:
	print(f'Removing outliers for {name} and writing it to files: {arg.openpath}{name}/pycali_results/{name}_cont_combined_new.txt and {arg.openpath}{name}/{name}_cont_combined_new.txt')
	pycali.remove_outliers(arg.openpath+name+"/pycali_results/"+name+"_cont_combined.txt", dev=5, doplot=True)
	os.system(f'cp ../pycali/{name}/pycali_results/{name}_cont_combined_new.txt ../pycali/{name}/{name}_cont_combined_new.txt')
	os.system(f'evince ../pycali/{name}/pycali_results/{name}_intercalibrated_lightcurves_only.pdf &')
	check = input(f'\nAfter seeing these plots, would you still like to re-intercalibrate {name}? If the pdf shows there is only one dataset or if the outlier plot looks good, type n or no\n')
	if check in ['y', 'yes', 'Y', 'Yes']:
		redone.append(name)
		os.mkdir(f'{openpath}{name}/redo_pycali_results')
        	# now, redo pycali with _new file. create redo_pycali_results directories in each recalibrated season
		cfg = pycali.Config()

		# except for the argument "fcont", the rest arguments are optional.
		# e.g.,  cfg.setup(fcont="data/ngc5548_cont.txt")

		cfg.setup(fcont=(arg.openpath+name+'/'+name+"_cont_combined_new.txt"),     # fcont is a string
          		nmcmc=10000, ptol=0.1,
          		scale_range_low=0.5, scale_range_up=2.0,
          		shift_range_low=-1.0, shift_range_up=1.0,
          		syserr_range_low=0.0, syserr_range_up=0.2,
          		errscale_range_low=0.5, errscale_range_up=2.0,
          		sigma_range_low=1.0e-4, sigma_range_up=1.0,
          		tau_range_low=1.0, tau_range_up=1.0e4,
          		fixed_scale=False, fixed_shift=False,
          		fixed_syserr=True, fixed_error_scale=True,
          		fixed_codes=[],
          		fixed_scalecodes=[],
          		flag_norm=True,
          		)	
		cfg.print_cfg()

		######################################################
		# do intercalibration
		#
		cali = pycali.Cali(cfg)  # create an instance
		cali.mcmc()              # do mcmc
		cali.get_best_params()   # calculate the best parameters
		cali.output()            # print output
		cali.recon()             # do reconstruction
		print('\n')
		print(f'You are re-intercalibrating {name}')
		print('\n')
		# plot results to PyCALI_results.pdf
		pycali.plot_results(cfg)

		# a simple plot
		#pycali.simple_plot(cfg)
        
		# AJF move output files to season folder
		os.rename("./PyCALI_results.pdf", "../pycali/"+name+"/"+"redo_pycali_results/"+"PyCALI_results_redo_"+name+".pdf")
		os.rename("./data/factor.txt", "../pycali/"+name+"/"+"redo_pycali_results/"+"factor_redo_"+name+".txt")
		os.rename("./data/levels.txt", "../pycali/"+name+"/"+"redo_pycali_results/"+"levels_redo_"+name+".txt")	
		os.rename("./data/param_input", "../pycali/"+name+"/"+"redo_pycali_results/"+"param_input_redo_"+name+".txt")
		os.rename("./data/posterior_sample_info.txt", "../pycali/"+name+"/"+"redo_pycali_results/"+"posterior_sample_info_redo_"+name+".txt")
		os.rename("./data/posterior_sample.txt", "../pycali/"+name+"/"+"redo_pycali_results/"+"posterior_sample_redo_"+name+".txt")
		os.rename("./data/PyCALI_output.txt", "../pycali/"+name+"/"+"redo_pycali_results/"+"PyCALI_output_redo_"+name+".txt")
		os.rename("./data/sample_info.txt", "../pycali/"+name+"/"+"redo_pycali_results/"+"sample_info_redo_"+name+".txt")
		os.rename("./data/sampler_state.txt", "../pycali/"+name+"/"+"redo_pycali_results/"+"sampler_state_redo_"+name+".txt")
		os.rename("./data/sample.txt", "../pycali/"+name+"/"+"redo_pycali_results/"+"sample_redo_"+name+".txt")
		os.rename("../pycali/"+name+"/"+name+"_cont_combined_new.txt_cali", "../pycali/"+name+"/"+"redo_pycali_results/"+name+"_cont_combined_new.txt_cali")
		os.rename("../pycali/"+name+"/"+name+"_cont_combined_new.txt_recon", "../pycali/"+name+"/"+"redo_pycali_results/"+name+"_cont_combined_new.txt_recon")
		os.rename("../pycali/"+name+"/"+name+"_cont_combined_new.txt_sort", "../pycali/"+name+"/"+"redo_pycali_results/"+name+"_cont_combined_new.txt_sort")
		os.system(f'cp ../pycali/{name}/redo_pycali_results/{name}_cont_combined_new.txt_cali {arg.openpath}for_create_cont')      
		os.rmdir("./data")
		print(f'{name} has been recalibrated and its _new.txt_cali file is now stored in for_create_cont')
	else:
		print(f'You have chosen to NOT re-intercalibrate {name} after viewing the outliers.')
		print(f'Thus, the original pycali file {name}_cont_combined.txt_cali has been moved to for_create_cont.')
		os.system(f'cp ../pycali/{name}/pycali_results/{name}_cont_combined.txt_cali {arg.openpath}for_create_cont')
		not_redone.append(name)           




for name in redone:
	print(f'Showing re-intercalibrated data for {name}')
	data_cali = np.loadtxt(arg.openpath+name+"/redo_pycali_results/"+name+"_cont_combined_new.txt_cali", usecols=(0, 1, 2))
	code = np.loadtxt(arg.openpath+name+"/redo_pycali_results/"+name+"_cont_combined_new.txt_cali", usecols=(3), dtype=str)
	fig = plt.figure(figsize=(10, 6))
	ax = fig.add_subplot(111)
	for c in np.unique(code):
    		idx = np.where(code == c)[0]
    		ax.errorbar(data_cali[idx, 0],  data_cali[idx, 1], yerr=data_cali[idx, 2], ls='none', marker='o', markersize=3, label=c)
	ax.legend()
	ax.set_title("Re-Intercalibrated data"+" "+name)
	plt.savefig(arg.openpath+name+"/redo_pycali_results/"+name+"_intercalibrated_lightcurves_only_new.pdf", format = 'pdf')
	plt.show()
        
for name in not_redone:
	"""
        print(f'Showing singly-intercalibrated data for {name}')
	data_cali = np.loadtxt(arg.openpath+name+"/pycali_results/"+name+"_cont_combined.txt_cali", usecols=(0, 1, 2))
	code = np.loadtxt(arg.openpath+name+"/pycali_results/"+name+"_cont_combined.txt_cali", usecols=(3), dtype=str)
	fig = plt.figure(figsize=(10, 6))
	ax = fig.add_subplot(111)
	for c in np.unique(code):
    		idx = np.where(code == c)[0]
    		ax.errorbar(data_cali[idx, 0],  data_cali[idx, 1], yerr=data_cali[idx, 2], ls='none', marker='o', markersize=3, label=c)
	ax.legend()
	ax.set_title("Singly-Intercalibrated data"+" "+name)
	plt.savefig(arg.openpath+name+"/pycali_results/"+name+"_intercalibrated_lightcurves_only_singly.pdf", format = 'pdf')
	plt.show()
	"""
	os.system(f'evince ../pycali/{name}/pycali_results/{name}_intercalibrated_lightcurves_only.pdf &')	


# do this for all_season data below

print(f'Removing outliers for all_seasons and writing it to file: {arg.openpath}all_seasons/pycali_results/all_seasons_cont_combined_new.txt and {arg.openpath}all_seasons/')
pycali.remove_outliers(arg.openpath + 'all_seasons/pycali_results/all_seasons_cont_combined.txt', dev=5, doplot=True)
print(f'\nView the outlier data. After seeing each seasons data, if one of them needed to be re-intercalibrated, you should probably reintercalibrate the all_season data. ')
redo_as = input('Would you like to? I would choose yes if the outlier plot for all_seasons looked bad or if at least one season had to be re-intercalibrated.\n')
if redo_as in ['yes', 'y', 'Yes', 'Y']:
	os.system(f'cp {openpath}all_seasons/pycali_results/all_seasons_cont_combined_new.txt {openpath}all_seasons/')
	os.mkdir(f'{openpath}all_seasons/redo_pycali_results')
	print(f'You are re-intercalibrating the all_seasons data.')
	print('\n')
	cfg = pycali.Config()
	cfg.setup(fcont=(arg.openpath+"all_seasons/all_seasons_cont_combined_new.txt"),     # fcont is a string
       		nmcmc=10000, ptol=0.1,
       		scale_range_low=0.5, scale_range_up=2.0,
       		shift_range_low=-1.0, shift_range_up=1.0,
       		syserr_range_low=0.0, syserr_range_up=0.2,
       		errscale_range_low=0.5, errscale_range_up=2.0,
      		sigma_range_low=1.0e-4, sigma_range_up=1.0,
       		tau_range_low=1.0, tau_range_up=1.0e4,
       		fixed_scale=False, fixed_shift=False,
       		fixed_syserr=True, fixed_error_scale=True,
       		fixed_codes=[],
       		fixed_scalecodes=[],
       		flag_norm=True,
       		)
	cfg.print_cfg()
	######################################################
	# do intercalibration
	#
	cali = pycali.Cali(cfg)  # create an instance
	cali.mcmc()              # do mcmc
	cali.get_best_params()   # calculate the best parameters
	cali.output()            # print output
	cali.recon()             # do reconstruction
	# plot results to PyCALI_results.pdf
	pycali.plot_results(cfg)
	# a simple plot
	#pycali.simple_plot(cfg)
        
	# AJF move output files to all_season folder
	# rename all files so they aren't re-written by loop 
	os.rename("./PyCALI_results.pdf", "../pycali/all_seasons/redo_pycali_results/PyCALI_results_redo_all_seasons.pdf")
	os.rename("./data/factor.txt", "../pycali/all_seasons/redo_pycali_results/factor_redo_all_seasons.txt")
	os.rename("./data/levels.txt", "../pycali/all_seasons/redo_pycali_results/levels_redo_all_seasons.txt")	
	os.rename("./data/param_input", "../pycali/all_seasons/redo_pycali_results/param_redo_input_all_seasons.txt")
	os.rename("./data/posterior_sample_info.txt", "../pycali/all_seasons/redo_pycali_results/posterior_sample_info_redo_all_seasons.txt")
	os.rename("./data/posterior_sample.txt", "../pycali/all_seasons/redo_pycali_results/posterior_sample_redo_all_seasons.txt")
	os.rename("./data/PyCALI_output.txt", "../pycali/all_seasons/redo_pycali_results/PyCALI_output_redo_all_seasons.txt")
	os.rename("./data/sample_info.txt", "../pycali/all_seasons/redo_pycali_results/sample_info_redo_all_seasons.txt")
	os.rename("./data/sampler_state.txt", "../pycali/all_seasons/redo_pycali_results/sampler_state_redo_all_seasons.txt")
	os.rename("./data/sample.txt", "../pycali/all_seasons/redo_pycali_results/sample_redo_all_seasons.txt")
	os.rename("../pycali/all_seasons/all_seasons_cont_combined_new.txt_cali", "../pycali/all_seasons/redo_pycali_results/all_seasons_cont_combined_new.txt_cali")
	os.rename("../pycali/all_seasons/all_seasons_cont_combined_new.txt_recon", "../pycali/all_seasons/redo_pycali_results/all_seasons_cont_combined_new.txt_recon")
	os.rename("../pycali/all_seasons/all_seasons_cont_combined_new.txt_sort", "../pycali/all_seasons/redo_pycali_results/all_seasons_cont_combined_new.txt_sort")
	os.system(f'cp {openpath}all_seasons/redo_pycali_results/all_seasons_cont_combined_new.txt_cali {openpath}for_create_cont/')      
	os.rmdir("./data")
	print(f'all_seasons has been recalibrated and its _new.txt_cali file is now stored in for_create_cont')


	print(f'Showing re-intercalibrated data for all_seasons')
	data_cali = np.loadtxt(arg.openpath + 'all_seasons/redo_pycali_results/all_seasons_cont_combined_new.txt_cali', usecols=(0, 1, 2))
	code = np.loadtxt(arg.openpath + 'all_seasons/redo_pycali_results/all_seasons_cont_combined_new.txt_cali', usecols=(3), dtype=str)
	fig = plt.figure(figsize=(10, 6))
	ax = fig.add_subplot(111)
	for c in np.unique(code):
    		idx = np.where(code == c)[0]
    		ax.errorbar(data_cali[idx, 0],  data_cali[idx, 1], yerr=data_cali[idx, 2], ls='none', marker='o', markersize=3, label=c)
	ax.legend()
	ax.set_title("Re-Intercalibrated data for all_seasons")
	plt.savefig(arg.openpath + 'all_seasons/redo_pycali_results/all_seasons_intercalibrated_lightcurves_only.pdf', format = 'pdf')
	plt.show()
else:
	print(f'You chose not to re-intercalibrate the all_seasons data. The singly-intercalibrated data is stored in for_create_cont!')
	os.system(f'cp {openpath}all_seasons/pycali_results/all_seasons_cont_combined.txt_cali {openpath}for_create_cont/all_seasons_cont_combined_new.txt_cali')
	os.system(f'evince ../pycali/all_seasons/pycali_results/all_seasons_intercalibrated_lightcurves_only.pdf &')
