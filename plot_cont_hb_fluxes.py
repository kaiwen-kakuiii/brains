"""
Created on Thu Jan 18 12:51:46 2024

@author: aidan
"""
# finsihed editing on 2/21/24 at 7:52 am
# updated 6/27

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

# initialize parser (copied in part from Thea's code sub_nar_updated.py)
parser = argparse.ArgumentParser(description='plots integrated h-beta flux and continuum flux for RM visualization') 
parser.add_argument('objpath', metavar='objpath', help='path to the directory containing the continuum and h-beta integrated fluxes (probably "integrated" directory)')
parser.add_argument('savepath', metavar='savepath', help='path to save the plot to (probably the "integrated" directory)')
parser.add_argument('season', metavar='season', help='name of season, i.e. "season#"')
#parser.add_argument('sliml', metavar='sliml', type = int, help = 'lower season# modified jd value (probably in the 1000s)')
#parser.add_argument('slimu', metavar='slimu', type = int, help = 'upper season# modified jd value (probably in the 1000s)')
parser.add_argument('choose', metavar='choose', type=int, help='1 is using cont_dynmod_no_ztf.txt, 2 is using cont_dynmod_ztf.txt, and 3 is using cont_dynmod_from_all_season_time_delay_checked.txt')
parser.add_argument('-em', '--ems', metavar='ems', type = str, nargs = '+', help = 'error methods to be used')
arg = parser.parse_args()
#sliml = arg.sliml - 1
#slimu = arg.slimu + 1

if arg.choose==1:

	# set up fig, load in data

	fig = plt.figure(figsize = (12,10))
	ax1 = fig.add_subplot(311)
	ax2 = fig.add_subplot(312)
	ax3 = fig.add_subplot(313)
	os.system(f'cp ../brains/{arg.season}/cont_dynmod_just_wiro_{arg.season}.txt ../integrated/{arg.season}/')
	data1 = np.loadtxt(arg.objpath+arg.season+'/'+'cont_dynmod_just_wiro_'+arg.season+'.txt', usecols = (0,1))
	data2 = np.loadtxt(arg.objpath+arg.season+'/'+'integrated_'+arg.season+'_hb.txt', delimiter =',', usecols = (0,1))


	jd, flux_con = data1[:,0], data1[:,1]
	jd_hb, flux_hb, err_hb = data2[:,0], data2[:,1], data2[:,2]
	#season = (jd>sliml) & (jd<slimu)
	#jd = jd[season]
	#flux_con = flux_con[season]
              
	# plot lightcurves for cont and hbeta

	ax1.plot(jd, flux_con,'k.', label = 'continuum 5100')
	ax1.set_ylim(min(flux_con)-np.std(flux_con), max(flux_con)+np.std(flux_con))
	ax1.set_xlim(jd_hb[0]-15, jd_hb[-1]+15)
	ax1.legend(loc=2)
	ax1.set_title('From cont_dynmod_just_wiro.txt')
	ax2.plot(jd_hb, flux_hb,'r.', label = 'h_beta')
	ax2.set_ylabel('Flux')
	ax2.set_xlim(jd_hb[0]-15, jd_hb[-1]+15)
	ax2.set_ylim(min(flux_hb)-np.std(flux_hb), max(flux_hb)+np.std(flux_hb))
	ax2.legend(loc = 2)

	# do stacked plot with same x axis, different y axis

	twin_stacked = ax3.twiny().twinx()
	twin1 = ax3.plot(jd, flux_con, 'k.', label = 'continuum 5100')
	ax3.set_xlim(jd_hb[0]-15, jd_hb[-1]+15)
	ax3.set_ylim(min(flux_con)-np.std(flux_con), max(flux_con)+np.std(flux_con))
	ax3.set_xlabel('JD')
	ax3.set_ylabel('Continuum Flux')
	ticks = np.round((np.linspace(jd_hb[0]-10, jd_hb[-1]+10, 20)),0)
	twin2 = twin_stacked.plot(jd_hb, flux_hb, 'r.', label = 'h_beta')
	twin_stacked.set_xlim(jd_hb[0]-15, jd_hb[-1]+15)
	twin_stacked.set_ylim(min(flux_hb)-np.std(flux_hb), max(flux_hb)+np.std(flux_hb))
	twin_stacked.set_ylabel('h_beta Flux')
	twin_stacked.set_xticklabels([])
	twin_stacked.set_xticks(ticks)
	ax3.set_xticks(ticks)
	ax3.grid()
	lns = twin1+twin2
	labs = [l.get_label() for l in lns]
	ax3.legend(lns, labs, loc='best', framealpha = 0.5)
	plt.savefig(arg.savepath+arg.season+'/'+'lightcurve_plot_just_wiro_'+arg.season+'.pdf', format = 'pdf')
	plt.show()
if arg.choose==2:

	# set up fig, load in data

	fig = plt.figure(figsize = (12,10))
	ax1 = fig.add_subplot(311)
	ax2 = fig.add_subplot(312)
	ax3 = fig.add_subplot(313)
	os.system(f'cp ../brains/{arg.season}/cont_dynmod_with_surveys_{arg.season}.txt ../integrated/{arg.season}/')
	data1 = np.loadtxt(arg.objpath+arg.season+'/'+'cont_dynmod_with_surveys_'+arg.season+'.txt', usecols = (0,1))
	data2 = np.loadtxt(arg.objpath+arg.season+'/'+'integrated_'+arg.season+'_hb.txt', delimiter =',', usecols = (0,1))

	jd, flux_con = data1[:,0], data1[:,1]
	jd_hb, flux_hb = data2[:,0], data2[:,1]
	#season = (jd>sliml) & (jd<slimu)
	#jd = jd[season]
	#flux_con = flux_con[season]
              
	# plot lightcurves for cont and hbeta

	ax1.plot(jd, flux_con,'k.', label = 'continuum 5100')
	ax1.set_ylim(min(flux_con)-np.std(flux_con), max(flux_con)+np.std(flux_con))
	ax1.set_xlim(jd_hb[0]-15, jd_hb[-1]+15)
	ax1.legend(loc=2)
	ax1.set_title('From cont_dynmod_with_surveys.txt')
	ax2.plot(jd_hb, flux_hb, 'r.', label = 'h_beta')
	ax2.set_ylabel('Flux')
	ax2.set_xlim(jd_hb[0]-15, jd_hb[-1]+15)
	ax2.set_ylim(min(flux_hb)-np.std(flux_hb), max(flux_hb)+np.std(flux_hb))
	ax2.legend(loc = 2)

	# do stacked plot with same x axis, different y axis

	twin_stacked = ax3.twiny().twinx()
	twin1 = ax3.plot(jd, flux_con, 'k.', label = 'continuum 5100')
	ax3.set_xlim(jd_hb[0]-15, jd_hb[-1]+15)
	ax3.set_ylim(min(flux_con)-np.std(flux_con), max(flux_con)+np.std(flux_con))
	ax3.set_xlabel('JD')
	ax3.set_ylabel('Continuum Flux')
	ticks = np.round((np.linspace(jd_hb[0]-10, jd_hb[-1]+10, 20)),0)
	twin2 = twin_stacked.plot(jd_hb, flux_hb, 'r.', label = 'h_beta')
	twin_stacked.set_xlim(jd_hb[0]-15, jd_hb[-1]+15)
	twin_stacked.set_ylim(min(flux_hb)-np.std(flux_hb), max(flux_hb)+np.std(flux_hb))
	twin_stacked.set_ylabel('h_beta Flux')
	twin_stacked.set_xticklabels([])
	twin_stacked.set_xticks(ticks)
	ax3.set_xticks(ticks)
	ax3.grid()
	lns = twin1+twin2
	labs = [l.get_label() for l in lns]
	ax3.legend(lns, labs, loc='best', framealpha = 0.5)
	plt.savefig(arg.savepath+arg.season+'/'+'lightcurve_plot_with_surveys_'+arg.season+'.pdf', format = 'pdf')
	plt.show()
   
out_cmd = open(arg.savepath+arg.season+'/'+'cmdoutput_plot_cont_hb_'+arg.season+'_time_delay_checked', 'w')     
for em in arg.ems:
	if arg.choose==3:
		
		# set up fig, load in data
		out_cmd.write(em+'\n')
		fig = plt.figure(figsize = (12,24))
		ax1 = fig.add_subplot(411)
		ax2 = fig.add_subplot(412)
		ax3 = fig.add_subplot(413)
		ax4 = fig.add_subplot(414)
		os.system(f'cp ../brains/{arg.season}/cont_dynmod_from_all_season_{arg.season}_time_delay_checked.txt ../integrated/{arg.season}/')
		data1 = np.loadtxt(arg.objpath+arg.season+'/'+'cont_dynmod_from_all_season_'+arg.season+'_time_delay_checked.txt', usecols = (0,1,2))
		data2 = np.loadtxt(arg.objpath+arg.season+'/'+'integrated_'+arg.season+'_hb.txt_'+em, delimiter =',', usecols = (0,1,2))	
		jd, flux_con, con_err = data1[:,0], data1[:,1], data1[:,2]
		#print(f'this is jd for season4: {jd}')
		jd_hb, flux_hb, err_hb = data2[:,0], data2[:,1], data2[:,2]
		#season = (jd>sliml) & (jd<slimu)
		#jd = jd[season]
		#flux_con = flux_con[season]
              
		# plot lightcurves for cont and hbeta

		ax1.plot(jd, flux_con,'k.', label = 'continuum 5100')
		ax1.set_ylim(min(flux_con)-np.std(flux_con), max(flux_con)+np.std(flux_con))
		ax1.set_xlim(jd[0]-5, jd_hb[-1]+5)
		ax1.legend(loc=2)
		ax1.set_title(f'From cont_dynmod_from_all_season_{arg.season}_time_delay_checked.txt\nand using {em} error method')
		ax2.errorbar(jd_hb, flux_hb, yerr = err_hb, color = 'r', label = 'h_beta', ls = 'none', capsize = 2.0)
		ax2.scatter(jd_hb, flux_hb, c = 'r', s = 3)
		ax2.set_ylabel('Flux')
		ax2.set_xlim(jd[0]-5, jd_hb[-1]+5)
		ax2.set_ylim(min(flux_hb)-np.std(flux_hb), max(flux_hb)+np.std(flux_hb))
		ax2.legend(loc = 2)
	
		# do stacked plot with same x axis, different y axis
	
		twin_stacked = ax3.twiny().twinx()
		twin1 = ax3.plot(jd, flux_con, 'k.', label = 'continuum 5100')
		ax3.set_xlim(jd[0]-5, jd_hb[-1]+5)
		ax3.set_ylim(min(flux_con)-np.std(flux_con), max(flux_con)+np.std(flux_con))
		ax3.set_ylabel('Continuum Flux')
		ticks = np.round((np.linspace(jd[0]-5, jd_hb[-1]+5, 20)),0)
		twin2 = twin_stacked.plot(jd_hb, flux_hb, 'r.', label = 'h_beta')
		twin_stacked.set_xlim(jd[0]-5, jd_hb[-1]+5)
		twin_stacked.set_ylim(min(flux_hb)-np.std(flux_hb), max(flux_hb)+np.std(flux_hb))
		twin_stacked.set_ylabel('h_beta Flux')
		twin_stacked.set_xticklabels([])
		twin_stacked.set_xticks(ticks)
		ax3.set_xticks(ticks)
		ax3.grid()
		lns = twin1+twin2
		labs = [l.get_label() for l in lns]
		ax3.legend(lns, labs, loc='best', framealpha = 0.5)

		# do normalized plot

		#data = np.loadtxt('../brains/season4/cont_dynmod_from_all_season_season4.txt')
		time, cont, conterr = data1[:,0], data1[:,1], data1[:,2]
		#hb_data = np.loadtxt('../integrated/season4/integrated_season4_hb.txt', delimiter = ',')
		time_hb, hb, hberr = data2[:,0], data2[:,1], data2[:,2]

		cont_norm = cont/np.average(cont)
		hb_norm = hb/np.average(hb)

		ss_hb = 0
		for i in range(len(hb)):
			ss_hb = ss_hb + (hb[i]-np.average(hb))**2

		hb_var = (1/(len(hb)-1))*ss_hb
		sd_hb = 0
		for i in range(len(hb)):
			sd_hb = sd_hb + (hberr[i])**2
		ds_hb = (1/len(hb))*sd_hb
		hb_mfv = np.sqrt(hb_var-ds_hb)/(np.average(hb))
		print(f'for {em}:')
		print(f'This is MFV of hb: {hb_mfv}')
		out_cmd.write('\nThis is MFV of hb: '+str(hb_mfv))

		ss_cont = 0
		for i in range(len(cont)):
			ss_cont = ss_cont + (cont[i]-np.average(cont))**2

		cont_var = (1/(len(cont)-1))*ss_cont
		sd_cont = 0
		for i in range(len(cont)):
			sd_cont = sd_cont + (conterr[i])**2
		ds_cont = (1/len(cont))*sd_cont
		cont_mfv = np.sqrt(cont_var-ds_cont)/(np.average(cont))
		print(f'This is MFV of cont: {cont_mfv}')
		out_cmd.write('\nThis is MFV of cont: '+str(cont_mfv))
	
		ratio = cont_mfv/hb_mfv

		print(f'This is ratio of cont_mfv to hb_mfv: {ratio}.\n\n')
		out_cmd.write('\nThis is ratio of cont_mfv to hb_mfv: '+str(ratio))
		ax4.set_xlabel('JD')
		ax4.set_ylabel('Normalized Continuum Flux, no MFV Adj.')
		ax4.plot(time, cont_norm, 'k.', label = 'Cont.')
		ax4.plot(time, cont_norm, 'k--', alpha = 0.25)
		ax4.plot(time_hb, hb_norm, 'r.', label = 'hb')
		ax4.plot(time_hb, hb_norm, 'r--', alpha = 0.25)
		ax4.legend(loc=2, framealpha = 0.5)

		plt.savefig(arg.savepath+arg.season+'/'+'lightcurve_plot_all_season_'+arg.season+'_time_delay_checked_err_'+ em+'.pdf', format = 'pdf')
		plt.show()


