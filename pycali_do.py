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
# AJF create list of season names

names = []
for i in range(1, (int(arg.season_num)+1)):
	name = 'season'+(str(i))
	names.append(name)
for i in range(3):
	print('\n')
print(f'You are going to intercalibrate {names}.')
for i in range(3):
	print('\n')      
#######################################################
# setup configurations, there are two ways:
# 1) load from a param file
#    cfg = pycali.Config("param.txt")
# 2) direct call setup()

for name in names:


	if os.path.isdir(f'{openpath}{name}/pycali_results'): 
		break                
	else:
		os.mkdir(f'{openpath}{name}/pycali_results')
	cfg = pycali.Config()

	# except for the argument "fcont", the rest arguments are optional.
	# e.g.,  cfg.setup(fcont="data/ngc5548_cont.txt")

	cfg.setup(fcont=(arg.openpath+name+'/'+name+"_cont_combined.txt"),     # fcont is a string
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
	print(f'You are intercalibrating {names}.')
	print('\n')
	# plot results to PyCALI_results.pdf
	pycali.plot_results(cfg)

	# a simple plot
	#pycali.simple_plot(cfg)
        
	# AJF move output files to season folder
	# rename all files so they aren't re-written by loop 
	os.rename("./PyCALI_results.pdf", "../pycali/"+name+"/"+"pycali_results/"+"PyCALI_results_"+name+".pdf")
	os.rename("./data/factor.txt", "../pycali/"+name+"/"+"pycali_results/"+"factor_"+name+".txt")
	os.rename("./data/levels.txt", "../pycali/"+name+"/"+"pycali_results/"+"levels_"+name+".txt")	
	os.rename("./data/param_input", "../pycali/"+name+"/"+"pycali_results/"+"param_input_"+name+".txt")
	os.rename("./data/posterior_sample_info.txt", "../pycali/"+name+"/"+"pycali_results/"+"posterior_sample_info_"+name+".txt")
	os.rename("./data/posterior_sample.txt", "../pycali/"+name+"/"+"pycali_results/"+"posterior_sample_"+name+".txt")
	os.rename("./data/PyCALI_output.txt", "../pycali/"+name+"/"+"pycali_results/"+"PyCALI_output_"+name+".txt")
	os.rename("./data/sample_info.txt", "../pycali/"+name+"/"+"pycali_results/"+"sample_info_"+name+".txt")
	os.rename("./data/sampler_state.txt", "../pycali/"+name+"/"+"pycali_results/"+"sampler_state_"+name+".txt")
	os.rename("./data/sample.txt", "../pycali/"+name+"/"+"pycali_results/"+"sample_"+name+".txt")
	os.rename("../pycali/"+name+"/"+name+"_cont_combined.txt_cali", "../pycali/"+name+"/"+"pycali_results/"+name+"_cont_combined.txt_cali")
	os.rename("../pycali/"+name+"/"+name+"_cont_combined.txt_recon", "../pycali/"+name+"/"+"pycali_results/"+name+"_cont_combined.txt_recon")
	os.rename("../pycali/"+name+"/"+name+"_cont_combined.txt_sort", "../pycali/"+name+"/"+"pycali_results/"+name+"_cont_combined.txt_sort")
	os.system(f'cp -f ../pycali/{name}/{name}_cont_combined.txt ../pycali/{name}/pycali_results/')    
	os.rmdir("./data")

for name in names:
	print(f'Showing intercalibrated data for {name}')
	data_cali = np.loadtxt(arg.openpath+name+"/pycali_results/"+name+"_cont_combined.txt_cali", usecols=(0, 1, 2))
	code = np.loadtxt(arg.openpath+name+"/pycali_results/"+name+"_cont_combined.txt_cali", usecols=(3), dtype=str)
	fig = plt.figure(figsize=(10, 6))
	ax = fig.add_subplot(111)
	for c in np.unique(code):
    		idx = np.where(code == c)[0]
    		ax.errorbar(data_cali[idx, 0],  data_cali[idx, 1], yerr=data_cali[idx, 2], ls='none', marker='o', markersize=3, label=c)
	ax.legend()
	ax.set_title("Intercalibrated data"+" "+name)
	plt.savefig(arg.openpath+name+"/pycali_results/"+name+"_intercalibrated_lightcurves_only.pdf", format = 'pdf')
	plt.show()


# do this for all_season data below
os.mkdir(f'{openpath}all_seasons/pycali_results')
cfg = pycali.Config()

cfg.setup(fcont=(arg.openpath+"all_seasons/all_seasons_cont_combined.txt"),     # fcont is a string
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
os.rename("./PyCALI_results.pdf", "../pycali/all_seasons/pycali_results/PyCALI_results_all_seasons.pdf")
os.rename("./data/factor.txt", "../pycali/all_seasons/pycali_results/factor_all_seasons.txt")
os.rename("./data/levels.txt", "../pycali/all_seasons/pycali_results/levels_all_seasons.txt")	
os.rename("./data/param_input", "../pycali/all_seasons/pycali_results/param_input_all_seasons.txt")
os.rename("./data/posterior_sample_info.txt", "../pycali/all_seasons/pycali_results/posterior_sample_info_all_seasons.txt")
os.rename("./data/posterior_sample.txt", "../pycali/all_seasons/pycali_results/posterior_sample_all_seasons.txt")
os.rename("./data/PyCALI_output.txt", "../pycali/all_seasons/pycali_results/PyCALI_output_all_seasons.txt")
os.rename("./data/sample_info.txt", "../pycali/all_seasons/pycali_results/sample_info_all_seasons.txt")
os.rename("./data/sampler_state.txt", "../pycali/all_seasons/pycali_results/sampler_state_all_seasons.txt")
os.rename("./data/sample.txt", "../pycali/all_seasons/pycali_results/sample_all_seasons.txt")
os.rename("../pycali/all_seasons/all_seasons_cont_combined.txt_cali", "../pycali/all_seasons/pycali_results/all_seasons_cont_combined.txt_cali")
os.rename("../pycali/all_seasons/all_seasons_cont_combined.txt_recon", "../pycali/all_seasons/pycali_results/all_seasons_cont_combined.txt_recon")
os.rename("../pycali/all_seasons/all_seasons_cont_combined.txt_sort", "../pycali/all_seasons/pycali_results/all_seasons_cont_combined.txt_sort")
os.system(f'cp -f ../pycali/all_seasons/all_seasons_cont_combined.txt ../pycali/all_seasons/pycali_results/')        
os.rmdir("./data")


print(f'Showing intercalibrated data for all_seasons')
data_cali = np.loadtxt(arg.openpath + 'all_seasons/pycali_results/all_seasons_cont_combined.txt_cali', usecols=(0, 1, 2))
code = np.loadtxt(arg.openpath + 'all_seasons/pycali_results/all_seasons_cont_combined.txt_cali', usecols=(3), dtype=str)
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111)
for c in np.unique(code):
    	idx = np.where(code == c)[0]
    	ax.errorbar(data_cali[idx, 0],  data_cali[idx, 1], yerr=data_cali[idx, 2], ls='none', marker='o', markersize=3, label=c)
ax.legend()
ax.set_title("Intercalibrated data for all_seasons")
plt.savefig(arg.openpath + 'all_seasons/pycali_results/all_seasons_intercalibrated_lightcurves_only.pdf', format = 'pdf')
plt.show()
