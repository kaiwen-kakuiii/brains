import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os.path

# updated 6/27 by AJF

# intialize data and input
parser = argparse.ArgumentParser(description='Adjust the variation amplitude of the continuum in reference to H-Beta')
parser.add_argument('path_cont', metavar = 'path_cont', help = 'path to brains continuum folder, probably ../brains/')
parser.add_argument('path_hb', metavar = 'path_hb', help = 'path to hb folder, probably ../integrated/')
parser.add_argument('season', metavar = 'season', help = 'season dir, i.e. season4, season2, etc.')
parser.add_argument('factor', metavar='factor', type = float, help='variance factor of continuum; i.e. "2" would mean continuum has mean fractional variance that is twice that of H_beta')
parser.add_argument('m_err', metavar='m_err', type = str, help = 'which error method for hbeta? m0,m1,m2,m3,m4')
arg = parser.parse_args()
factor = arg.factor
season = arg.season
path_cont = arg.path_cont
path_hb = arg.path_hb
m_err = arg.m_err

# check to see if files exist
names = []
if os.path.isfile(path_cont+season+'/cont_dynmod_from_all_season_'+season+'_unchecked_time_delay.txt'):
	names.append(str('cont_dynmod_from_all_season_{}_unchecked_time_delay.txt'.format(season)))
	print('checked1')
if os.path.isfile(path_cont+season+'/cont_dynmod_from_all_season_'+season+'_time_delay_checked.txt'):
	names.append(str('cont_dynmod_from_all_season_{}_time_delay_checked.txt'.format(season)))	
	print('checked2')


for name in names:
	out2 = open(path_cont+season+'/'+'cmdoutput_mfv_cont_adj_'+name + '_'+m_err+'_factor_'+str(factor), 'w')
	cont_data = np.loadtxt(path_cont+season+'/'+name)
	time, cont, conterr = cont_data[:,0], cont_data[:,1], cont_data[:,2]
	hb = np.loadtxt(f'{path_hb}{season}/integrated_{season}_hb.txt_{m_err}', delimiter = ',')
	timehb, hbf, hbferr = hb[:,0], hb[:,1], hb[:,2]
	len_before = len(cont)	
	og_len = len(cont)
	# do hb calcs
	hbf_av = np.average(hbf)
	sum_square = 0
	hbf_norm = hbf/np.average(hbf)
	for i in range(len(hbf)):
		sum_square = sum_square+(hbf[i]-hbf_av)**2
	hbf_var = (1/(len(hbf)-1))*sum_square
	sum_delta = 0
	for i in range(len(hbf)):
		sum_delta = sum_delta + (hbferr[i])**2
	delta_square = (1/len(hbf))*sum_delta
	hbf_mfv = np.sqrt(hbf_var-delta_square)/hbf_av
	print(f'\nFor {name} below:\n')
	print(f'H_Beta Mean Fractional Variation: {hbf_mfv}\n')
	out2.write('\nFor '+name+' below:\n')
	out2.write('H_Beta Mean Fractional Variation: '+str(hbf_mfv)+'\n')
	# do original continuum calcs
	cont_av = np.average(cont)
	cont_norm = cont/np.average(cont)
	cont_norm_1 = cont_norm
	sum_square_cont = 0
	for i in range(len(cont)):
		sum_square_cont = sum_square_cont+(cont[i]-cont_av)**2
	cont_var = (1/(len(cont)-1))*sum_square_cont
	sum_delta_cont = 0
	for i in range(len(cont)):
		sum_delta_cont = sum_delta_cont + (conterr[i])**2
	delta_square_cont = (1/len(cont))*sum_delta_cont
	old_mfv_cont = np.sqrt(cont_var-delta_square_cont)/cont_av
	print(f'this is cont_var: {cont_var} and delta_square: {delta_square_cont}\n')
	out2.write('this is cont_var: '+str(cont_var)+' and delta_square: '+str(delta_square_cont)+'\n')
	if cont_var - delta_square_cont < 0:
		print(f'There appears to be a continuum data point that deviates significantly from the others. Please check the plot.\n')
		out2.write('There appears to be a continuum data point that deviates significantly from the others. Please check the plot.\n')
		print(f'Take note of the y value of the point and round down to the nearest half (i.e. if point is at 3.75, round down to 3.5. If point is at 4.3, round down to 4).\n')
		out2.write('Take note of the y value of the point and round down to the nearest half (i.e. if point is at 3.75, round down to 3.5. If point is at 4.3, round down to 4).\n')
		fig = plt.figure(figsize = (12,10))
		plt.plot(time, cont_norm, 'k.')
		plt.xlabel('Time, JD')
		plt.ylabel(f'Normalized Continuum')
		plt.yticks(ticks = np.arange(0, np.round(max(cont_norm)+0.5, 0), 0.5))
		plt.grid()
		plt.show()
		#len_before = len_list[j]
		check = input(f'Do you want to remove this point? y or n\n')
		out2.write('Do you want to remove this point? y or n\n'+check+'\n')
		while check in ['yes', 'y', 'Yes', 'Y']:
			remove_count = 0
			check2 = input(f'What is the y value (rounded down as instructed before) of this point?\n')
			out2.write('What is the y value (rounded down as instructed before) of this point?\n'+check2+'\n')
			for i in range(len(cont_norm)):
				if cont_norm[i]>int(check2) and remove_count < 1:
					print(f'\n\n\nDELETING row {i} with jd {time[i]} and normalized continuum value of {cont[i]/np.average(cont)}\n\n\n')
					out2.write('\n\nDELETING row '+str(i)+' with jd '+str(time[i])+' and normalized continuum value of '+str(cont[i]/np.average(cont))+'\n\n')
					#print(f'this is time before: {time}')
					time = np.delete(time, i, 0)
					#print(f'this is time after: {time}')
					#print(f'this is cont before: {cont}')
					cont = np.delete(cont, i, 0)
					#print(f'this is cont after: {cont}')
					#print(f'this is conterr before: {conterr}')
					conterr = np.delete(conterr, i, 0)
					#print(f'this is conterr after: {conterr}')
					remove_count+=1
			#print(f'this is time outside loop: {time}')
			cont_av = np.average(cont)
			cont_norm = cont/np.average(cont)
			cont_norm_1 = cont_norm
			sum_square_cont = 0
			for i in range(len(cont)):
				sum_square_cont = sum_square_cont+(cont[i]-cont_av)**2
			cont_var = (1/(len(cont)-1))*sum_square_cont
			sum_delta_cont = 0
			for i in range(len(cont)):
				sum_delta_cont = sum_delta_cont + (conterr[i])**2
			delta_square_cont = (1/len(cont))*sum_delta_cont
			cont_norm = cont/np.average(cont)
			print(f'this is now cont_var: {cont_var} and delta_square: {delta_square_cont}\n')
			out2.write('this is now cont_var: '+str(cont_var)+'and delta_square: '+str(delta_square_cont)+'\n')
			print(f'The point has been removed if matching criteria; check the plot to see that you are satisfied.\n')
			out2.write('The point has been removed if matching criteria; check the plot to see that you are satisfied.\n')
			fig = plt.figure(figsize = (12,10))
			plt.xlabel('Time, JD')
			plt.ylabel(f'New Normalized Continuum')
			plt.yticks(ticks = np.arange(0, np.round(max(cont_norm)+0.5, 0), 0.5))
			plt.plot(time, cont_norm, 'k.')
			plt.grid()
			plt.show()
			#len_list[j] = int(len_list[j])-int(remove_count)
			print(f'New Num Cont. Points: {len(cont)}; Old Num Cont. Points: {len_before}')
			out2.write('New num cont. points: '+str(len(cont))+'; Old num cont. points: '+str(len_before)+'\n')
			len_before = len(cont)
			check = input(f'Do you want to remove another point? y or n\n') 
			out2.write('Do you want to remove another point? y or n\n'+check+'\n')                                       		                                                                               
	cont_mfv = np.sqrt((cont_var-delta_square_cont))/cont_av
	print(f'Continuum (Unadjusted) Mean Fractional Variation: {cont_mfv}\n')
	out2.write('Continuum (Unadjusted) Mean Fractional Variation: '+str(cont_mfv)+'\n')
	
   	#initialize plot
	fig = plt.figure(figsize = (12,10))
	ax1 = fig.add_subplot(311)
	ax2 = fig.add_subplot(312)
	ax3 = fig.add_subplot(313)
	
     	# plot original continuum
	ax1.plot(time, cont_norm, 'k.', label = 'cont')
	ax1.plot(timehb, hbf_norm, 'r.', label = 'hb')
	ax1.plot(time, cont_norm, 'k--', alpha = 0.1)
	ax1.plot(timehb, hbf_norm, 'r--', alpha = 0.1)
	ax1.set_title(f'Unadjusted, Num Points: {og_len}, for {name}', loc = 'left')
	ax1.set_xticklabels([])
	ax1.legend()

	# do adjusted cont calcs
	sub = ((factor*hbf_mfv*cont_av)-(np.sqrt((cont_var-delta_square_cont))))/(factor*hbf_mfv)
	cont_adj = cont-sub
	cont_av = np.average(cont_adj)
	sum_square_cont = 0
	for i in range(len(cont_adj)):
		sum_square_cont = sum_square_cont+(cont_adj[i]-cont_av)**2
	cont_var = (1/(len(cont_adj)-1))*sum_square_cont
	sum_delta_cont = 0
	for i in range(len(cont_adj)):
		sum_delta_cont = sum_delta_cont + (conterr[i])**2
	delta_square_cont = (1/len(cont_adj))*sum_delta_cont		
	cont_mfv = np.sqrt((cont_var-delta_square_cont))/cont_av
	print(f'Continuum (Adjusted) Mean Fractional Variation: {cont_mfv}\n')
	out2.write('Continuum (Adjusted) Mean Fractional Variation: '+str(cont_mfv)+'\n')

	#check math
	print(f'Constant Flux Subtracted: {sub:.4e}\n')
	sub_out2 = str("{:.4e}".format(sub))
	out2.write('Constant Flux Subtracted: '+sub_out2+'\n')
	print(f'Ratio factor you gave of Adjusted Continuum MFV / H_Beta MVF: {factor}\n')
	out2.write('Ratio factor you gave of Adjusted Continuum MFV / H_Beta MVF: '+str(factor)+'\n')
	print(f'Calculated ratio factor from calculated Adjusted Continuum MFV and HB MFV: {np.round(cont_mfv/hbf_mfv, 4)}\n\n\n\n\n\n')
	out2.write('Calculated ratio factor from calculated Adjusted Continuum MFV and HB MFV: '+str(np.round(cont_mfv/hbf_mfv, 4))+'\n\n\n')

	# plot adjusted cont
	hbf_norm = hbf/np.average(hbf)
	cont_norm_adj = cont_adj/np.average(cont_adj)
	ax2.plot(time, cont_norm_adj, 'b.', label = 'cont_adj')
	ax2.plot(timehb, hbf_norm, 'r.', label = 'hb')
	ax2.plot(time, cont_norm_adj, 'b--', alpha = 0.1)
	ax2.plot(timehb, hbf_norm, 'r--', alpha = 0.1)
	ax2.set_title(f'Adjusted by subtracting {sub:.4e}', loc = 'left')
	ax2.set_xticklabels([])
	ax2.legend()

	#plot all together
	ax3.plot(timehb, hbf_norm, 'r.', label = 'hb')
	ax3.plot(time, cont_norm_adj, 'b.', label = 'cont_adj')
	ax3.plot(time, cont_norm, 'k.', label = 'cont', alpha = 0.1)
	ax3.plot(timehb, hbf_norm, 'r--', alpha = 0.1)
	ax3.plot(time, cont_norm_adj, 'b--', alpha = 0.1)
	ax3.plot(time, cont_norm, 'k--', alpha = 0.1)
	ax3.set_title(f'Old and New Normalized Continuums and H-Beta\nOld MFV ratio: {np.round((old_mfv_cont/hbf_mfv),4)}  |  New MFV: {np.round((cont_mfv/hbf_mfv), 4)}', loc = 'left', fontsize = 11)
	ax3.legend(loc = 'lower right')
	plt.savefig(path_cont+season+'/plot_mfv_adjusted_'+name+'_'+m_err+'.png', format = 'png') 
	plt.show()
			
	# write MFV file			
	out = open(path_cont+season+'/'+'mfv_adj_factor_'+str(factor)+'_'+name +'_'+m_err, 'w')
	for i in range(len(cont_adj)):
		out.write(str(time[i]) + '   ' + str(cont_adj[i]) + '   ' + str(conterr[i]) + "\n")	
	out.close()

"""			
name_list2 = []
for name in name_list:
	name_list2.append(name[-2:])	
count2 = 0
for line3 in open(f'{path_cont}all_seasons/all_seasons_cont_combined.txt'):
	line3=line3.split()
	count2+=1
	if line3[0] == '#':
		if str(line3[1][-2:]) not in name_list2:
			print(f'\n\n\nWriting to file for {str(line3[1][-2:])} not in {name_list2}.\n\n\n')
			out = open(path_cont+'/all_seasons/all_seasons_mfv_adjusted_cont_combined.txt', 'a')
			skip_row = int(count2)
			length_data = int(line3[2])
			cont_data = np.loadtxt(f'{path_cont}all_seasons/all_seasons_cont_combined.txt', skiprows = skip_row, max_rows = length_data)
			if length_data!=0:
				time, cont, conterr = cont_data[:,0], cont_data[:,1], cont_data[:,2]
				out.write('# {} {}'.format(line3[1], len(cont)) +'\n')
				for i in range(len(cont)):
					out.write(str(time[i]) + '   ' + str(cont[i]) + '   ' + str(conterr[i]) + "\n")

out.close()
			                       
				

for line in open(f'{path_cont}{season}/{season}_cont_combined.txt'):
	line = line.split()
	check_hash = line[0]
	name = line[1]
	check_len = line[2]
	if check_hash == '#':
		hash_list.append(count)
		name_list.append(name)
		len_list.append(line[2])		
	count += 1


"""
