# creator Aidan Ferguson on 2/26/24 for use with BRAINS
# updated 04/15/24 as per Yan-Rong's emails on 4/13/24 and earlier about only using cont. data from time delay
# updated 5/7/24 to have code search time delay to make sure there were enough data points for continuum reconstruction
# updated 6/27 for completion

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

# AJF initialize parsers
parser = argparse.ArgumentParser(description='create continuum files for BRAINS after using PyCali')
parser.add_argument('orpath', metavar='orpath', help='path the object folder (probably "original" directory) ') 
parser.add_argument('pycalipath', metavar='pycali_path', help='path to the pycali folder (probably "pycali" directory) ')
parser.add_argument('savepath', metavar='save_path', help='path to place you want to save the output spectra text file (probably "brains" directory')
parser.add_argument('season_num', metavar='season_num', help='total number of seasons, i.e. "5"')
parser.add_argument('-td', '--time_delay', metavar='time_delay', type = int, nargs = '*', help = 'approximate time delay of emission line from continuum')
arg = parser.parse_args()
td = (arg.time_delay)

# AJF create list of season names
names = []
for i in range(1, (int(arg.season_num)+1)):
	j = str(i)
	name = 'season'+(j)
	names.append(name)

out_cmd = open(arg.savepath+'/'+'cmdoutput_most_recent_run_create_cont_dynmod.txt', 'w')

# write continuum files from pycali
# write the "with surveys" files, no cross-checking with these

for name in names:
	# read in flux.lst_season#
	if os.path.isfile(arg.pycalipath+'for_create_cont/'+name+'_cont_combined_new.txt_cali'):
		data = np.loadtxt(arg.pycalipath+'for_create_cont/'+name+'_cont_combined_new.txt_cali', dtype = np.float64, usecols = (0,1,2))
		jd, flux, err = (data[:,0]), data[:,1], data[:,2]
	
		# print out these columns with one space between each column
		out = open(arg.savepath+name+'/'+'cont_dynmod_with_surveys_'+name+'.txt', 'w')
		for i in range(len(flux)):
			out.write(str(jd[i]) + ' ' + str(flux[i]) + ' ' + str(err[i]) + '\n')
		out.close()
	else:
		data = np.loadtxt(arg.pycalipath+'for_create_cont/'+name+'_cont_combined_.txt_cali', dtype = np.float64, usecols = (0,1,2))
		jd, flux, err = (data[:,0]), data[:,1], data[:,2]
	
		# print out these columns with one space between each column
		out = open(arg.savepath+name+'/'+'cont_dynmod_with_surveys_'+name+'.txt', 'w')
		for i in range(len(flux)):
			out.write(str(jd[i]) + ' ' + str(flux[i]) + ' ' + str(err[i]) + '\n')
		out.close()



# load in jd, flux, error arrays from all_season_cont_combined post-pycali using entire continuum across all seasons without caring about number of data points within time delay

  
for name in names:
	print(f'\nThis is the loop not checking to see if there is data within the time delay limits for {name}.')
	out_cmd.write('\nThis is the loop not checking to see if there is data within the time delay limits for '+name)
	if os.path.isfile(arg.pycalipath+'for_create_cont/'+name+'_cont_combined_new.txt_cali'): #using re-intercalibrated data
		print(f'Found _new file for {name}, so using that.')
		out_cmd.write('\nFound _new file for '+name+' so using that.')
		i = int(name[-1])
		data1 = np.loadtxt(arg.pycalipath+'for_create_cont/'+ name +'_cont_combined_new.txt_cali', dtype = np.float64, usecols = (0,1,2))
		jd1 = data1[:,0]
		jd1_val = int(jd1[0])
		data_all = np.loadtxt(arg.pycalipath+'for_create_cont/all_seasons_cont_combined_new.txt_cali', dtype = np.float64, usecols = (0,1,2))		
		jd, flux, err = data_all[:,0], data_all[:,1], data_all[:,2]
		data = np.array((jd, flux, err))       
		if i == 1: # season1 only               
			if int(jd[0]) < jd1_val - 1.4*td[i-1]:
				jd1_val = jd1_val - 1.4*td[i-1]
			else:
				jd1_val = int(jd[0])
		else: #other seasons
			jd1_val = jd1_val - 1.4*td[i-1]	
		season_lims = ((jd>jd1_val) & (jd<(jd1[-1]+2)))
		data_final= [vals[season_lims] for vals in data]		                
		jd_final = data_final[0]
		print(f'this is the final jd list with unchecked time data in {name}: {jd_final}\n')
		jd, flux, err = data_final[0], data_final[1], data_final[2]
		out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_unchecked_time_delay.txt', 'w')
		for j in range(len(flux)):
			out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
		out.close()

	else: # using singly-intercalibrated data
		print(f'\nNo _new file found, so using singly-intercalibrated data.')
		out_cmd.write('\nNo _new file found, so using singly-intercalibrated data.')
		i = int(name[-1])
		data1 = np.loadtxt(arg.pycalipath+'for_create_cont/'+ name +'_cont_combined.txt_cali', dtype = np.float64, usecols = (0,1,2))
		jd1 = data1[:,0]
		jd1_val = int(jd1[0])
		data_all = np.loadtxt(arg.pycalipath+'for_create_cont/all_seasons_cont_combined.txt_cali', dtype = np.float64, usecols = (0,1,2))		
		jd, flux, err = data_all[:,0], data_all[:,1], data_all[:,2]
		data = np.array((jd, flux, err))       
		if i == 1: # season1 only               
			if int(jd[0]) < jd1_val - 1.4*td[i-1]:
				jd1_val = jd1_val - 1.4*td[i-1]
			else:
				jd1_val = int(jd[0])
		else: #other seasons
			jd1_val = jd1_val - 1.4*td[i-1]	
		season_lims = ((jd>jd1_val) & (jd<(jd1[-1]+2)))
		data_final= [vals[season_lims] for vals in data]		                
		jd_final = data_final[0]
		print(f'this is the final jd list with unchecked time data in {name}: {jd_final}\n\n')
		jd, flux, err = data_final[0], data_final[1], data_final[2]
		out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_unchecked_time_delay.txt', 'w')
		for j in range(len(flux)):
			out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
		out.close()


# load in jd, flux, error arrays from all_season_cont_combined post-pycali using entire continuum across all seasons AND SEARCH TIME DELAY FOR DATA COUNT

  
for name in names:
	print(f'\nThis is the loop that IS checking to see if there is data within the time delay limits for {name}.')
	out_cmd.write('\n\nThis is the loop that IS checking to see if there is data within the time delay limits for '+name)
	if os.path.isfile(arg.pycalipath+'for_create_cont/'+name+'_cont_combined_new.txt_cali'): #using re-intercalibrated data
		print(f'\nFound _new file for {name}, so using that.')
		out_cmd.write('\n\nFound _new file for '+name+' so using that.')
		i = int(name[-1])
		print(f'\nthis is time delay for {name}: {td[i-1]}')
		out_cmd.write('\nthis is time delay for '+name+', so using that: '+str(td[i-1]))
		data1 = np.loadtxt(arg.pycalipath+'for_create_cont/'+ name +'_cont_combined_new.txt_cali', dtype = np.float64, usecols = (0,1,2))
		jd1 = data1[:,0]
		jd1_val = int(jd1[0])
		data_all = np.loadtxt(arg.pycalipath+'for_create_cont/all_seasons_cont_combined_new.txt_cali', dtype = np.float64, usecols = (0,1,2))		
		jd, flux, err = data_all[:,0], data_all[:,1], data_all[:,2]
		data = np.array((jd, flux, err))
		jd_check = jd1_val
		if i == 1: # season1 only               
			if int(jd[0]) < jd1_val - 1.4*td[i-1]:
				jd1_val = jd1_val - 1.4*td[i-1]
			else:
				jd1_val = int(jd[0])
			#jd1_val = jd1_val - 1.4*td[i-1]	
			count = 0
			while jd_check > jd1_val:  
				day_check = (jd<jd_check) & (jd>jd_check-1) # note above does not break until full range (or more) of time delay is checked
				data_check = [vals[day_check] for vals in data]
				#print(f'this is data_check: {data_check}')
				if len(data_check[0]) > 0:
					count = count+1
				jd_check = jd_check -1
				#print(f'this is count: {count}')
				#print(f'this is jd_check and jd1_val: {jd_check}, {jd1_val}.')
			if count<((td[i-1])/3):
				print(f'You may want to try to find more ZTF/ASAS-SN data,\nas only {count} days within {np.round(td[i-1],2)} days prior to the first WIRO continuum day have data which may make it difficult to reconstruct the continuum.')
				out_cmd.write('\nYou may want to try to find more ZTF/ASAS-SN data,\nas only '+str(count)+' days within '+str(np.round(td[i-1],2))+' days prior to the first WIRO continuum day have data which may make it difficult to reconstruct the continuum.\n')
			season_lims = ((jd>jd1_val) & (jd<(jd1[-1]+2)))
			data_final= [vals[season_lims] for vals in data]		                
			jd_final = data_final[0]
			print(f'this is the final jd list with time delay check data in {name}: {jd_final}\n')
			jd, flux, err = data_final[0], data_final[1], data_final[2]
			out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_time_delay_checked.txt', 'w')
			for j in range(len(flux)):
				out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
			out.close()
		else: #other seasons
			jd_check = jd1_val
			jd1_val = jd1_val - 1.4*td[i-1]	
			count = 0
			while count<((td[i-1])/3) or jd_check > jd1_val:  # checks to see if there is td/3 days that contain data and breaks when this condition is met
				day_check = (jd<jd_check) & (jd>jd_check-1) # note above does not break until full range (or more) of time delay is checked
				data_check = [vals[day_check] for vals in data]
				#print(f'this is data_check: {data_check}')
				if len(data_check[0]) > 0:
					count = count+1
				jd_check = jd_check -1
				#print(f'this is count: {count}')
				#print(f'this is jd_check and jd1_val: {jd_check}, {jd1_val}.')
			if jd_check+1>jd1_val: # if there are td/3 data points between season's first jd and the first jd minus the time delay
				season_lims = ((jd>jd1_val) & (jd<(jd1[-1]+2)))
				data_final= [vals[season_lims] for vals in data]		                
				jd_final = data_final[0]
				print(f'this is the final jd list with time delay check data in {name}: {jd_final}\n')
				jd, flux, err = data_final[0], data_final[1], data_final[2]
				out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_time_delay_checked.txt', 'w')
				for j in range(len(flux)):
					out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
				out.close()
			else: # if you have to search further past time delay to find enough data
				season_lims = ((jd>jd_check-1) & (jd<(jd1[-1]+2)))
				data_final= [vals[season_lims] for vals in data]		                
				jd_final = data_final[0]
				print(f'this is the final jd list with time delay check data in {name}: {jd_final}\n')
				print(f'Not enough data found in {np.round(1.4*int(td[i-1]),2)} days before first h_beta data for {name}, so searched an additional ~{np.round(jd1_val-jd_check,2)} days prior before finding total {np.round((td[i-1])/3,0)} days with data.\n')
				out_cmd.write('\nNot enough data found in '+str(np.round(1.4*int(td[i-1]),2))+' days before first h_beta data for '+name+', so searched an additional ~'+str(np.round(jd1_val-jd_check,2))+' days prior before finding total {np.round(td[i-1]/3,0)} days with data.\n')
				jd, flux, err = data_final[0], data_final[1], data_final[2]
				out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_time_delay_checked.txt', 'w')
				for j in range(len(flux)):
					out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
				out.close()
	else: # using singly-intercalibrated data
		print(f'\nNo _new file found, so using singly-intercalibrated data.')
		out_cmd.write('\nNo _new file found, so using singly-intercalibrated data.')
		i = int(name[-1])
		print(f'this is the time delay for {name}: {td[i-1]}')
		out_cmd.write('\nthis is the time delay for '+name+': '+str(td[i-1]))
		data1 = np.loadtxt(arg.pycalipath+'for_create_cont/'+ name +'_cont_combined.txt_cali', dtype = np.float64, usecols = (0,1,2))
		jd1 = data1[:,0]
		jd1_val = int(jd1[0])
		data_all = np.loadtxt(arg.pycalipath+'for_create_cont/all_seasons_cont_combined.txt_cali', dtype = np.float64, usecols = (0,1,2))		
		jd, flux, err = data_all[:,0], data_all[:,1], data_all[:,2]
		data = np.array((jd, flux, err))
		jd_check = jd1_val
		if i == 1: # season1 only               
			if int(jd[0]) < jd1_val - 1.4*td[i-1]:
				jd1_val = jd1_val - 1.4*td[i-1]
			else:
				jd1_val = int(jd[0])
			#jd1_val = jd1_val - 1.4*td[i-1]	
			count = 0
			while jd_check > jd1_val:  
				day_check = (jd<jd_check) & (jd>jd_check-1) # note above does not break until full range (or more) of time delay is checked
				data_check = [vals[day_check] for vals in data]
				#print(f'this is data_check: {data_check}')
				if len(data_check[0]) > 0:
					count = count+1
				jd_check = jd_check -1
				#print(f'this is count: {count}')
				#print(f'this is jd_check and jd1_val: {jd_check}, {jd1_val}.')
			if count<((td[i-1])/3):
				print(f'You may want to try to find more ZTF/ASAS-SN data,\nas only {count} days within {np.round(td[i-1],2)} days prior to the first WIRO continuum day have data which may make it difficult to reconstruct the continuum.')
				out_cmd.write('\nYou may want to try to find more ZTF/ASAS-SN data,\nas only '+str(count)+' days within '+str(np.round(td[i-1],2))+' days prior to the first WIRO continuum day have data which may make it difficult to reconstruct the continuum.')
			season_lims = ((jd>jd1_val) & (jd<(jd1[-1]+2)))
			data_final= [vals[season_lims] for vals in data]		                
			jd_final = data_final[0]
			print(f'this is the final jd list with data in {name}: {jd_final}\n')
			jd, flux, err = data_final[0], data_final[1], data_final[2]
			out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_time_delay_checked.txt', 'w')
			for j in range(len(flux)):
				out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
			out.close()
		else:
			jd_check = jd1_val
			jd1_val = jd1_val - 1.4*td[i-1]	
			count = 0
			while count<((td[i-1])/3) or jd_check > jd1_val:
				day_check = (jd<jd_check) & (jd>jd_check-1)
				data_check = [vals[day_check] for vals in data]
				if len(data_check) > 0:
					count = count+1
				jd_check = jd_check -1
			if jd_check+1>jd1_val:		
				season_lims = ((jd>jd1_val) & (jd<(jd1[-1]+2)))
				data_final= [vals[season_lims] for vals in data]		                
				jd_final = data_final[0]
				print(f'this is the final jd list with time delay check data in {name}: {jd_final}\n')
				jd, flux, err = data_final[0], data_final[1], data_final[2]
				out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_time_delay_checked.txt', 'w')
				for j in range(len(flux)):
					out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
				out.close()
			else:
				season_lims = ((jd>jd_check+1) & (jd<(jd1[-1]+2)))
				data_final= [vals[season_lims] for vals in data]		                
				jd_final = data_final[0]
				print(f'this is the final jd list with time delay check data in {name}: {jd_final}\n')
				print(f'Not enough data found in {np.round(1.4*int(td[i-1]),2)} days before first h_beta data for {name}, so searched an additional ~{np.round(jd1_val-jd_check,2)} days prior before finding 10 days with data.\n')
				out_cmd.write('\nNot enough data found in '+str(np.round(1.4*int(td[i-1]),2))+' days before first h_beta data for '+name+', so searched an additional ~'+str(np.round(jd1_val-jd_check,2))+' days prior before finding 10 days with data.\n')
				jd, flux, err = data_final[0], data_final[1], data_final[2]
				out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_time_delay_checked.txt', 'w')
				for j in range(len(flux)):
					out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
				out.close()


"""

# do same process, but for mean fractional variance-adjusted data

queue = str(input(f'Did you previously use the mean-fractional variance code (split_conts.py) to adjust the MFV of the continuum?\n')
if queue in ['y', 'yes', 'Y', 'Yes']:
	queue2 = list(input(f'\nFor what seasons did you do this for? List them as numbers, like "1 3 5", for example.\n')

	# AJF create list of season names for MFV adjusted
	names = []
	for i in range(len(queue2)):
		j = str(queue2[i])
		name = 'season'+(j)
		names.append(name)

# load in jd, flux, error arrays from all_season_cont_combined post-pycali using entire continuum across all seasons without caring about number of data points within time delay

	for name in names:
		print(f'\nThis is the loop not checking to see if there is data within the time delay limits for {name}.')
		if os.path.isfile(arg.pycalipath+'for_create_cont/'+name+'_cont_combined_new.txt_cali'): #using re-intercalibrated data
			print(f'Found _new file for {name}, so using that.')
			i = int(name[-1])
			data1 = np.loadtxt(arg.pycalipath+'for_create_cont/'+ name +'_cont_combined_new.txt_cali', dtype = np.float64, usecols = (0,1,2))
			jd1 = data1[:,0]
			jd1_val = int(jd1[0])
			data_all = np.loadtxt(arg.pycalipath+'for_create_cont/all_seasons_cont_combined_new.txt_cali', dtype = np.float64, usecols = (0,1,2))		
			jd, flux, err = data_all[:,0], data_all[:,1], data_all[:,2]
			data = np.array((jd, flux, err))       
			if i == 1: # season1 only               
				if int(jd[0]) < jd1_val - 1.4*td[i-1]:
					jd1_val = jd1_val - 1.4*td[i-1]
				else:
					jd1_val = int(jd[0])
			else: #other seasons
				jd1_val = jd1_val - 1.4*td[i-1]	
			season_lims = ((jd>jd1_val) & (jd<(jd1[-1]+2)))
			data_final= [vals[season_lims] for vals in data]		                
			jd_final = data_final[0]
			print(f'this is the final jd list with data in {name}: {jd_final}\n')
			jd, flux, err = data_final[0], data_final[1], data_final[2]
			out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_unchecked_data_count_in_time_delay_mfv_adj.txt', 'w')
			for j in range(len(flux)):
				out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
			out.close()

		else: # using singly-intercalibrated data
			print(f'\nNo _new file found, so using singly-intercalibrated data.')
			i = int(name[-1])
			data1 = np.loadtxt(arg.pycalipath+'for_create_cont/'+ name +'_cont_combined.txt_cali', dtype = np.float64, usecols = (0,1,2))
			jd1 = data1[:,0]
			jd1_val = int(jd1[0])
			data_all = np.loadtxt(arg.pycalipath+'for_create_cont/all_seasons_cont_combined.txt_cali', dtype = np.float64, usecols = (0,1,2))		
			jd, flux, err = data_all[:,0], data_all[:,1], data_all[:,2]
			data = np.array((jd, flux, err))       
			if i == 1: # season1 only               
				if int(jd[0]) < jd1_val - 1.4*td[i-1]:
					jd1_val = jd1_val - 1.4*td[i-1]
				else:
					jd1_val = int(jd[0])
			else: #other seasons
				jd1_val = jd1_val - 1.4*td[i-1]	
			season_lims = ((jd>jd1_val) & (jd<(jd1[-1]+2)))
			data_final= [vals[season_lims] for vals in data]		                
			jd_final = data_final[0]
			print(f'this is the final jd list with data in {name}: {jd_final}\n\n')
			jd, flux, err = data_final[0], data_final[1], data_final[2]
			out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_unchecked_data_count_in_time_delay_mfv_adj.txt', 'w')
			for j in range(len(flux)):
				out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
			out.close()


	# load in jd, flux, error arrays from all_season_cont_combined post-pycali using entire continuum across all seasons AND SEARCH TIME DELAY FOR DATA COUNT

  
	for name in names:
		print(f'\nThis is the loop that IS checking to see if there is data within the time delay limits for {name}.')
		if os.path.isfile(arg.pycalipath+'for_create_cont/'+name+'_cont_combined_new.txt_cali'): #using re-intercalibrated data
			print(f'Found _new file for {name}, so using that.')
			i = int(name[-1])
			print(f'this is time delay for {name}: {td[i-1]}')
			data1 = np.loadtxt(arg.pycalipath+'for_create_cont/'+ name +'_cont_combined_new.txt_cali', dtype = np.float64, usecols = (0,1,2))
			jd1 = data1[:,0]
			jd1_val = int(jd1[0])
			data_all = np.loadtxt(arg.pycalipath+'for_create_cont/all_seasons_cont_combined_new.txt_cali', dtype = np.float64, usecols = (0,1,2))		
			jd, flux, err = data_all[:,0], data_all[:,1], data_all[:,2]
			data = np.array((jd, flux, err))
			jd_check = jd1_val
			if i == 1: # season1 only               
				if int(jd[0]) < jd1_val - 1.4*td[i-1]:
					jd1_val = jd1_val - 1.4*td[i-1]
				else:
					jd1_val = int(jd[0])
				#jd1_val = jd1_val - 1.4*td[i-1]	
				count = 0
				while jd_check > jd1_val:  
					day_check = (jd<jd_check) & (jd>jd_check-1) # note above does not break until full range (or more) of time delay is checked
					data_check = [vals[day_check] for vals in data]
					#print(f'this is data_check: {data_check}')
					if len(data_check[0]) > 0:
						count = count+1
					jd_check = jd_check -1
					#print(f'this is count: {count}')
					#print(f'this is jd_check and jd1_val: {jd_check}, {jd1_val}.')
				if count<((td[i-1])/3):
					print(f'You may want to try to find more ZTF/ASAS-SN data,\nas only {count} days within {np.round(td[i-1],2)} days prior to the first WIRO continuum day have data which may make it difficult to reconstruct the continuum.')
				season_lims = ((jd>jd1_val) & (jd<(jd1[-1]+2)))
				data_final= [vals[season_lims] for vals in data]		                
				jd_final = data_final[0]
				print(f'this is the final jd list with data in {name}: {jd_final}\n')
				jd, flux, err = data_final[0], data_final[1], data_final[2]
				out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_mfv_adj.txt', 'w')
				for j in range(len(flux)):
					out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
				out.close()
			else: #other seasons
				jd_check = jd1_val
				jd1_val = jd1_val - 1.4*td[i-1]	
				count = 0
				while count<((td[i-1])/3) or jd_check > jd1_val:  # checks to see if there is td/3 days that contain data and breaks when this condition is met
					day_check = (jd<jd_check) & (jd>jd_check-1) # note above does not break until full range (or more) of time delay is checked
					data_check = [vals[day_check] for vals in data]
					#print(f'this is data_check: {data_check}')
					if len(data_check[0]) > 0:
						count = count+1
					jd_check = jd_check -1
					#print(f'this is count: {count}')
					#print(f'this is jd_check and jd1_val: {jd_check}, {jd1_val}.')
				if jd_check+1>jd1_val: # if there are td/3 data points between season's first jd and the first jd minus the time delay
					season_lims = ((jd>jd1_val) & (jd<(jd1[-1]+2)))
					data_final= [vals[season_lims] for vals in data]		                
					jd_final = data_final[0]
					print(f'this is the final jd list with data in {name}: {jd_final}\n')
					jd, flux, err = data_final[0], data_final[1], data_final[2]
					out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_mfv_adj.txt', 'w')
					for j in range(len(flux)):
						out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
					out.close()
				else: # if you have to search further past time delay to find enough data
					season_lims = ((jd>jd_check+1) & (jd<(jd1[-1]+2)))
					data_final= [vals[season_lims] for vals in data]		                
					jd_final = data_final[0]
					print(f'this is the final jd list with data in {name}: {jd_final}\n')
					print(f'Not enough data found in {np.round(1.4*int(td[i-1]),2)} days before first h_beta data for {name}, so searched an additional ~{np.round(jd1_val-jd_check,2)} days prior before finding 10 days with data.\n')
					jd, flux, err = data_final[0], data_final[1], data_final[2]
					out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_mfv_adj.txt', 'w')
					for j in range(len(flux)):
						out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
					out.close()
		else: # using singly-intercalibrated data
			print(f'\nNo _new file found, so using singly-intercalibrated data.')
			i = int(name[-1])
			print(f'this is the time delay for {name}: {td[i-1]}')
			data1 = np.loadtxt(arg.pycalipath+'for_create_cont/'+ name +'_cont_combined.txt_cali', dtype = np.float64, usecols = (0,1,2))
			jd1 = data1[:,0]
			jd1_val = int(jd1[0])
			data_all = np.loadtxt(arg.pycalipath+'for_create_cont/all_seasons_cont_combined.txt_cali', dtype = np.float64, usecols = (0,1,2))		
			jd, flux, err = data_all[:,0], data_all[:,1], data_all[:,2]
			data = np.array((jd, flux, err))
			jd_check = jd1_val
			if i == 1: # season1 only               
				if int(jd[0]) < jd1_val - 1.4*td[i-1]:
					jd1_val = jd1_val - 1.4*td[i-1]
				else:
					jd1_val = int(jd[0])
				#jd1_val = jd1_val - 1.4*td[i-1]	
				count = 0
				while jd_check > jd1_val:  
					day_check = (jd<jd_check) & (jd>jd_check-1) # note above does not break until full range (or more) of time delay is checked
					data_check = [vals[day_check] for vals in data]
					#print(f'this is data_check: {data_check}')
					if len(data_check[0]) > 0:
						count = count+1
					jd_check = jd_check -1
					#print(f'this is count: {count}')
					#print(f'this is jd_check and jd1_val: {jd_check}, {jd1_val}.')
				if count<((td[i-1])/3):
					print(f'You may want to try to find more ZTF/ASAS-SN data,\nas only {count} days within {np.round(td[i-1],2)} days prior to the first WIRO continuum day have data which may make it difficult to reconstruct the continuum.')
				season_lims = ((jd>jd1_val) & (jd<(jd1[-1]+2)))
				data_final= [vals[season_lims] for vals in data]		                
				jd_final = data_final[0]
				print(f'this is the final jd list with data in {name}: {jd_final}\n')
				jd, flux, err = data_final[0], data_final[1], data_final[2]
				out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_mfv_adj.txt', 'w')
				for j in range(len(flux)):
					out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
				out.close()
			else:
				jd_check = jd1_val
				jd1_val = jd1_val - 1.4*td[i-1]	
				count = 0
				while count<((td[i-1])/3) or jd_check > jd1_val:
					day_check = (jd<jd_check) & (jd>jd_check-1)
					data_check = [vals[day_check] for vals in data]
					if len(data_check) > 0:
						count = count+1
					jd_check = jd_check -1
				if jd_check+1>jd1_val:		
					season_lims = ((jd>jd1_val) & (jd<(jd1[-1]+2)))
					data_final= [vals[season_lims] for vals in data]		                
					jd_final = data_final[0]
					print(f'this is the final jd list with data in {name}: {jd_final}\n')
					jd, flux, err = data_final[0], data_final[1], data_final[2]
					out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_mfv_adj.txt', 'w')
					for j in range(len(flux)):
						out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
					out.close()
				else:
					season_lims = ((jd>jd_check+1) & (jd<(jd1[-1]+2)))
					data_final= [vals[season_lims] for vals in data]		                
					jd_final = data_final[0]
					print(f'this is the final jd list with data in {name}: {jd_final}\n')
					print(f'Not enough data found in {np.round(1.4*int(td[i-1]),2)} days before first h_beta data for {name}, so searched an additional ~{np.round(jd1_val-jd_check,2)} days prior before finding 10 days with data.\n')
					jd, flux, err = data_final[0], data_final[1], data_final[2]
					out = open(arg.savepath+name+'/'+'cont_dynmod_from_all_season_'+name+'_mfv_adj.txt', 'w')
					for j in range(len(flux)):
						out.write(str(jd[j]) + ' ' + str(flux[j]) + ' ' + str(err[j]) + '\n')
					out.close()

"""
