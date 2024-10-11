import numpy as np
import astropy.units as u
from astropy.time import Time
import argparse
import pandas as pd
from ztfquery import lightcurve, query
import sys
import os

# edited amd updated 02/29/24 by Aidan Ferguson
# added all_season stuff

def get_ztf(ra, dec, objpath, arcsearch):
	"""Uses ZTF query to download ZTF lightcurves from the web and filters based on input criteria.                                                                                               
    
	Parameters
	----------
	ra : :class:`float'
	RA of object in decimal degrees
                                                                                                   
	dec : :class:`float`
	Dec of object in decimal degrees

	r : :class:`integer`
	Radius to search in arcseconds. Default is 3
        
        objpath : :class'string'
        Path to object directory

 	Returns
 	-------
 	:class:`.txt file`
	A text file named "ztf" that has jd, flux, flux err and filter
        
        AJF NOTES:
        
        User does not need to modify this function.
        
	"""

	# TZ intitialize the query
	lcq = lightcurve.LCQuery()
	arcsearch = float(arcsearch)       
	result = lcq.from_position(ra, dec, arcsearch)
	data = result.data

	# TZ pull only g filter data
	datatable = data[(data.filtercode != 'zi') & (data.exptime <= 100)]

	# TZ pull jd 
	jd = np.array(datatable['hjd'])
	jd -= 2457000

	# TZ pull filter code
	filt = np.array(datatable['filtercode'])

	# TZ pull mag, magerr and convert from mags to cgs units
	flux0 = 10**(-0.4*(datatable['mag']-21))
	flux_err0 = flux0 * datatable['magerr'] / 1.1816
        
	# TZ assign Janskys as the unit
	flux1 = np.array(flux0)*u.Jy
	flux_err1 = np.array(flux_err0)*u.Jy

	#convert to cgs
	flux2 = flux1.to(u.erg/u.s/u.cm**2/u.Hz)
	flux_err2 = flux_err1.to(u.erg/u.s/u.cm**2/u.Hz)

	#convert to angstroms
	#define reference wavelength
	wave = 5100.* u.AA #angstroms

	flux3 = flux2.to(u.erg/u.s/u.cm**2/u.AA, u.spectral_density(wave))
	flux_err3 = flux_err2.to(u.erg/u.s/u.cm**2/u.AA, u.spectral_density(wave))

	flux = flux3.value
	flux_err = flux_err3.value

	#write output to a text file called ztf.txt
	out = open(objpath+"ztf.txt", "w")
	out.write('JD-2457000' + ',' + 'flux' + ',' + 'err' + ',' + 'filter' +'\n')
	for i in range(len(jd)):
		out.write(str(jd[i]) + ',' + str(flux[i]) + ',' + str(flux_err[i]) + ',' + str(filt[i]) + "\n")
	out.close()

	return jd, flux, flux_err, filt

def get_asas(objpath):
	"""
	uses ASAS-SN data to supplement WIRO continuum data
	
        objpath : :class'string'
	Path to object directory
	
	Returns:
        --------
        data arrays asas_data and asas_data_cont;
        asas_data contains extra       
                                
	"""
	# AJF read in ASAS-SN csv, check to see if units are in mJy       
	asas = pd.read_csv(objpath+'asas.csv', delimiter = ',')
	colnames = asas.columns.tolist()
	if any('mJy' in word for word in colnames) is not True:
		print('\n')
		print(colnames)
		que = input('\nFlux unit may not be in mJy (see above). Still proceed?\n')                        
		if que not in ['yes', 'y', 'Yes', 'Y']:
			sys.exit()
	if 'hjd' not in colnames:
		asas_hjd = asas['HJD']-2457000
	else:
		asas_hjd = asas['hjd']-2457000
	if 'flux (mJy)' and 'flux(mJy)' not in colnames:
		asas_flux = asas['flux']	
	else:
		if 'flux (mJy)' in colnames:
        		asas_flux = asas['flux (mJy)'] 
		else:
			asas_flux = asas['flux(mJy)']           
	if 'flux err' not in colnames:
		asas_err = asas['flux_err']
	else:
		if 'flux_err' not in colnames:
			asas_err = asas['flux err']
		else:
			asas_err = asas['fluxerr']	                
	if 'camera' not in colnames:
		asas_cam = asas['Camera']
	else:
		asas_cam = asas['camera']
	asas_data = np.array((asas_hjd, asas_flux, asas_err, asas_cam))        	
	asas_data_cont = np.array((asas_hjd, asas_flux, asas_err))	

        	        
	# convert to cgs and F_lamda        
	# TZ assign Janskys as the unit
	flux0 = asas_flux/1000
	flux_err0 = asas_err/1000
	flux1 = np.array(flux0)*u.Jy
	flux_err1 = np.array(flux_err0)*u.Jy

	#convert to cgs
	flux2 = flux1.to(u.erg/u.s/u.cm**2/u.Hz)
	flux_err2 = flux_err1.to(u.erg/u.s/u.cm**2/u.Hz)

	#convert to angstroms
	#define reference wavelength
	wave = 5100.* u.AA #angstroms

	flux3 = flux2.to(u.erg/u.s/u.cm**2/u.AA, u.spectral_density(wave))
	flux_err3 = flux_err2.to(u.erg/u.s/u.cm**2/u.AA, u.spectral_density(wave))

	flux_asas = flux3.value
	flux_err_asas = flux_err3.value
	#write output to a text file called asas.txt
	out = open(objpath+"asas.txt", "w")
	out.write('HJD-2457000' + ',' + 'flux' + ',' + 'err' + ',' + 'filter' +'\n')
	for i in range(len(asas_hjd)):
		out.write(str(asas_hjd[i]) + ',' + str(flux_asas[i]) + ',' + str(flux_err_asas[i]) + ',' + str(asas_cam[i]) + "\n")
	out.close()

        
	return asas_hjd, flux_asas, flux_err_asas, asas_cam
	        
	        
	        
	        
	        
        
def split_seasons(jd, flux, err, filtercode, objpath, season_num, hjd, flux_asas, err_asas, cam_asas):
	"""Takes text file with ztf data and splits it by filter and season> Reads in wiro data and splits it into seasons.                                                                    
                                                                                                   
	Parameters                                                                                     
	----------                                                                                     
	jd : :class:`float'                                                                          
	Reduced julian date of data point                                                             
                                                                                                   
	flux : :class:`float`                                                                           
	Flux of data point in erg/s/cm^2/AA                                                   
                                                                                                   
	err : :class:`integer`                                                                        
	Error in the flux measurement in erg/s/cm^2/AA

	filtercode : :class'string'
	Filter data was observed in. zg for g band and zr for r band
        
        objpath : :class'string'
        Path to object directory

	Returns                                                                                        
	-------                                                                                        
	g :class:`array-like`                                                                            
	An array of arrays of ztf g band data split by season

	r :class:`array-like`                                                                            
	An array of arrays of ztf r band data split by season

	wirodata :class:`array-like`                                                                            
 	An array of arrays of wiro data (cont and hb) split by season
        
        hb_med :class:'array-like'
        An array of median hb flux values for each season
        
        AJF NOTES:
        
        need to add/remove seasons in this function; below parameters are set for IRAS data, which currently has 5 seasons worth of data from WIRO
        
	"""
	#np.set_printoptions(threshold=sys.maxsize)
        
	# TZ put jd, flux, err into an array for code-compactness
	data = np.array((jd, flux, err, filtercode))

	# TZ define season bounds, indiv seasons
	season1 = (jd > 1338) & (jd < 1564)
	season2 = (jd > 1700) & (jd < 1980)      
	season3 = (jd > 2080) & (jd < 2335) 
	season4 = (jd > 2440) & (jd < 2698) 
	season5 = (jd > 2785) & (jd < 3071)        
        
	# TZ split data into seasons
	data1 = [vals[season1] for vals in data]
	data2 = [vals[season2] for vals in data]
	data3 = [vals[season3] for vals in data]
	data4 = [vals[season4] for vals in data]
	data5 = [vals[season5] for vals in data]

	# TZ make indices for g band data for each season
	g1 = data1[3] == 'zg'
	g2 = data2[3] == 'zg'
	g3 = data3[3] == 'zg'
	g4 = data4[3] == 'zg'
	g5 = data5[3] == 'zg'
       
	# TZ pull g band data for each season
	gdata1 = ([vals[g1] for vals in data1])
	gdata2 = ([vals[g2] for vals in data2])
	gdata3 = ([vals[g3] for vals in data3])
	gdata4 = ([vals[g4] for vals in data4])
	gdata5 = ([vals[g5] for vals in data5])
	g = np.array((gdata1, gdata2, gdata3, gdata4, gdata5))    
        
        
	"""
        # TZ/AJF make indices for r band data for each season
	r1 = data1[3] == 'zr'
	r2 = data2[3] == 'zr'
	r3 = data3[3] == 'zr'
	r4 = data4[3] == 'zr'
	r5 = data5[3] == 'zr'
        
        # TZ pull r band data for each season
	rdata1 = [vals[r1] for vals in data1]
	rdata2 = [vals[r2] for vals in data2]
	rdata3 = [vals[r3] for vals in data3]
	rdata4 = [vals[r4] for vals in data4]
	rdata5 = [vals[r5] for vals in data5]
	#r = np.array((rdata1, rdata2, rdata3, rdata4, rdata5))
	"""
	# ASAS-SN data below, indiv seasons
	adata = np.array((hjd, flux_asas, err_asas, cam_asas))     
	aseason1 = (hjd > 1338) & (hjd < 1564)
	aseason2 = (hjd > 1700) & (hjd < 1980)      
	aseason3 = (hjd > 2080) & (hjd < 2335) 
	aseason4 = (hjd > 2440) & (hjd < 2698) 
	aseason5 = (hjd > 2785) & (hjd < 3071) 
                
	# TZ split data into seasons
	adata1 = [vals[aseason1] for vals in adata]
	adata2 = [vals[aseason2] for vals in adata]
	adata3 = [vals[aseason3] for vals in adata]
	adata4 = [vals[aseason4] for vals in adata]
	adata5 = [vals[aseason5] for vals in adata]                 
	asas_all = np.array((adata1, adata2, adata3, adata4, adata5))
	                	                
	# TZ read in wiro data
	wiro = pd.read_csv(objpath+'flux.lst', delim_whitespace=True)
	wiro_jd = wiro['jd'] - 2457000
	wiro_flux = wiro['f5100']
	wiro_err = wiro['ef5100']
	hbflux = np.array(wiro['fhb'])
	hberr = np.array(wiro['efhb'])

	# TZ season bounds for hb data, same as continuum bounds
	wiro_data = np.array((wiro_jd, wiro_flux, wiro_err, hbflux, hberr))
	wiro_all_season_cont = np.array((wiro_jd, wiro_flux, wiro_err)) # needed for pycali all_seasons
	wiroseason1 = (wiro_jd > 1338) & (wiro_jd < 1564)
	wiroseason2 = (wiro_jd > 1708) & (wiro_jd < 1980)
	wiroseason3 = (wiro_jd > 2088) & (wiro_jd < 2335) 
	wiroseason4 = (wiro_jd > 2450) & (wiro_jd < 2698) 
	wiroseason5 = (wiro_jd > 2799) & (wiro_jd < 3071)
	
	# TZ split wiro data into seasons
	wiro1 = [vals[wiroseason1] for vals in wiro_data]
	wiro2 = [vals[wiroseason2] for vals in wiro_data]
	wiro3 = [vals[wiroseason3] for vals in wiro_data]
	wiro4 = [vals[wiroseason4] for vals in wiro_data]
	wiro5 = [vals[wiroseason5] for vals in wiro_data]
	wirodata = np.array((wiro1,wiro2,wiro3,wiro4,wiro5))
	
	################################################################################################################
	################################################################################################################
	################################################################################################################        
	# AJF create g and r band full set of data across all seasons if needed for pycali, created from indiv seasons
	empty = []
	g_all = []
	#r_all = []
        
        # g-band ztf below
	for i in range(1,5):
        	g_all.append(empty)
                
	for i in range(1,(int(season_num))+1):
		if len(g[i-1][0])!=0:
			for j in range(len(g[i-1])):
				g_all[j] = np.concatenate((g_all[j], ((g[i-1][j]))))
			print(f'Season{i} has g-band ZTF data')		
		else:
                	print('\n')
                	print(f'ALERT:')
                	print('\n')
                	print(f'Season{i} does not have g-band ZTF data')
                	print('\n')
                	continue                                					

	# r-band below  
	"""
        empty = []
	r_all = []
	for i in range(1,(int(season_num)+1)):
        	g_all.append(empty)
                
	for i in range(1,(int(season_num))+1):
		if len(r[i-1][0]!=0:
                	for j in range(len(g[i-1])):
                        	r_all[j] = np.concatenate((r_all[j], ((r[i-1][j]))))
                        	print(f'this is r_all: {r_all}')		
		else:
                	print(f'Season{i} does not have r-band ZTF data.')
			continue                                					
		continue
	"""

	# filter different cameras asassn data
	cameras = []
	dic = {}
	cam_dic = {}
	cam_list = []
	cs_dic = {}
	cs_list = []
	cameras.append(adata[3][0])
	for i in range(1, len(adata[3])):
		if adata[3][i] not in cameras:
			cameras.append(adata[3][i])
	print(f'These are all the different cameras used by ASAS-SN across all seasons: {cameras}')        	
	for i in range(len(asas_all)):
		for j in range(len(cameras)):
			check = asas_all[i][3] == cameras[j]
			if not isinstance(check, np.ndarray):
				check2 = []
				check = asas_all[i][3] == cameras[j]
				check2.append(check)
				check = check2
				check2 = [] 	
			if any(check) is True:
				key_i_j = '{}_{}'.format(str(int(i)+1), cameras[j])                        
				dic[key_i_j] = ([vals[check] for vals in asas_all[i]])
				cam_list.append(cameras[j])
				cs_list.append(key_i_j)	        
		key = (int(i)+1)
		cam_dic[key] = cam_list
		cs_dic[key] = cs_list
		cam_list = []
		cs_list = []	        
	dict_keys = list(dic.keys())
 
        
        
        # create all_season data for asassn below from indiv seasons
	all_dic = {}
	initial = []
	empty = np.array([])
	try1=[]
	for i in range(4):
		initial.append(empty)
		try1.append(empty)  
	for i in range(len(cameras)):
		all_dic[cameras[i]]=initial                
	for j in range(len(cameras)):
		check = str(cameras[j])
		results = [key for key in dict_keys if check in key]
		for k in range(len(results)):
			for w in range(4):
				if any(check in word for word in results):					
					dic_array = (dic.get(results[k]))                                        
					try1[w] = np.concatenate( (try1[w], dic_array[w]) )                       
		all_dic.update({str(check): try1})
		try1=[]
		for i in range(4):
			try1.append(empty)                                                                                                 
	# rename and print season data
	asas_all_dic = all_dic # all_seasons, catagorized by camera                              
	asas_cs_dic = dic # season uncatagorized dictionary (catgorized by camera)
       
	for i in range(int(season_num)):	
		if len(g[i][0])!=0 and len(cam_dic.get(i+1))==0:
			print(f'Season{i+1} has g-band ZTF data, but no ASAS-SN data')
			print('\n')                        
		elif len(g[i][0])==0 and len(cam_dic.get(i+1))==0:
			print(f'Season{i+1} has no g-band ZTF data and no ASAS-SN data')	                               
			print('\n') 
		elif len(g[i][0])==0 and len(cam_dic.get(i+1))!=0:
			print(f'Season{i+1} has ASAS-SN data, but no g-band ZTF data')	  
			print('\n') 
		else:
			print(f'Season{i+1} has g-band ZTF data and ASAS-SN data')	                               
			print('\n') 	  


	################################################################################################################
	################################################################################################################
	################################################################################################################
        # create all_season data across full range; season1 first date to last season last date
                        
	# AJF/TZ define seasonall bounds
	seasonall = (jd > 1338) & (jd < 3072)        
        
	# AJF/TZ split data into seasons
	dataall = [vals[seasonall] for vals in data]

	# AJF/TZ make indices for g band data for each season
	gall = dataall[3] == 'zg'

       
	# AJF/TZ pull g band data for each season
	gdataall = ([vals[gall] for vals in dataall])
	g_all_meta = np.array((gdataall))    
        
        
	"""
        # TZ/AJF make indices for r band data for each season
	rall = dataall[3] == 'zr'

        
        # TZ pull r band data for each season
	rdataall = [vals[rall] for vals in dataall]
	#r = np.array((rdataall))
	"""
	# ASAS-SN data below, all seasons combined
    
	aseasonall = (hjd > 1338) & (hjd < 3072) 
                
	# TZ split data into seasons
	adataall = [vals[aseasonall] for vals in adata]                
	asas_all_meta = np.array((adataall))

	# filter different cameras asassn data
	dic_meta = {}
	cam_dic_meta = {}
	cam_list_meta = []
	cs_dic_meta = {}
	cs_list_meta = []        	

	for j in range(len(cameras)):
		check = asas_all_meta[3] == cameras[j]
		if not isinstance(check, np.ndarray):
			check2 = []
			check = asas_all_meta[3] == cameras[j]
			check2.append(check)
			check = check2
			check2 = [] 	
		if any(check) is True:
			key_i_j_meta = '{}_{}'.format(str(1), cameras[j])                        
			dic_meta[key_i_j_meta] = ([vals[check] for vals in asas_all_meta])
			cam_list_meta.append(cameras[j])
			cs_list_meta.append(key_i_j_meta)	        
	key_meta = (int(1))
	cam_dic_meta[key_meta] = cam_list_meta
	cs_dic_meta[key_meta] = cs_list_meta	        
	dict_keys_meta = list(dic_meta.keys())
 
        
        # create all_season data for asassn below from meta-season
	all_dic_meta = {}
	initial = []
	empty = np.array([])
	try1=[]
	for i in range(4):
		initial.append(empty)
		try1.append(empty)  
	for i in range(len(cameras)):
		all_dic_meta[cameras[i]]=initial                
	for j in range(len(cameras)):
		check = str(cameras[j])
		results = [key_meta for key_meta in dict_keys_meta if check in key_meta]
		for k in range(len(results)):
			for w in range(4):
				if any(check in word for word in results):					
					dic_array_meta = (dic_meta.get(results[k]))                                        
					try1[w] = np.concatenate( (try1[w], dic_array_meta[w]) )                       
		all_dic_meta.update({str(check): try1})
		try1=[]
		for i in range(4):
			try1.append(empty)  

		         
	# rename and print meta-season data
	asas_all_dic_meta = all_dic_meta # all_seasons, catagorized by camera                              	  

	# end function	                      
	return g, asas_cs_dic, wirodata, g_all, asas_all_dic, wiro_all_season_cont, cs_dic, g_all_meta, asas_all_dic_meta  


def write_cont(wiro, g, asas_cs_dic, num_seasons, seasonnames, savepath, g_all, wiro_all_season_cont, asas_all_dic, cs_dic, g_all_meta, asas_all_dic_meta):
	"""Writes sepearate seasons of continuum data to a text file in the PyCALI format

	Parameters
    	----------
    	telescope : :class:`array-like`
    	List of data from the different telescopes. Each telescope is an array containing arrays for each season.

    	num_seasons : :class:`integer
    	Number of seasons

    	seasonnames : :class:`array-like'
    	List of directory names for each season
        
        hb_median : :class:'array-like'
        List of median hbeta flux value for each season
        
        objpath : :class'string'
        Path to object directory

    	Returns                                                                
    	-------                                                                                        
    	:class:`.txt file`
    	Outputs are text files named with cont_combined_season_ for each season with the season number filling in the blank
        
        AJF NOTES:
        
        User changes this function to get correct number of seasons and correct bands (if r is included)        
        
    	"""

	# make array
	num = np.arange(1,int(num_seasons)+1,1) 
        
	# open a text file for each season
	for i in range(int(num_seasons)):
		letters = ['WIRO', 'G_band_ztf'] 
		telescope = [wiro, g]
		out = open(savepath+seasonnames[i]+'/'+seasonnames[i]+'_cont_combined.txt', 'w')	
		#loop through letters for wiro and g-band-ztf to write pycali format
		for j, label in enumerate(letters):           
			num_points = len(telescope[j][i][0])
			if num_points != 0:
				out.write('# {} '.format(label) + str(num_points) +'\n')
				jd_season = np.array(telescope[j][i][0])
				flux_season = np.array(telescope[j][i][1])
				err_season = np.array(telescope[j][i][2])

				#now write values for that telescope
				for k in range(num_points):
					out.write(str(jd_season[k]) + '   ' + str(flux_season[k]) + '   ' + str(err_season[k]) + "\n")
		letters = cs_dic.get(i+1)
		telescope = []	# create new letters for ASAS-SN data
		for w in range(len(letters)):
			telescope.append(asas_cs_dic.get(letters[w]))	
		# now, loop through new letters for ASAS-SN data		
		for j, label in enumerate(letters):           
			num_points = len(telescope[j][0])
			if num_points != 0:
				out.write('# {} '.format(label) + str(num_points) +'\n')

				jd_season = np.array(telescope[j][0])
				flux_season = np.array(telescope[j][1])
				err_season = np.array(telescope[j][2])

				#now write values for that telescope
				for k in range(num_points):
					out.write(str(jd_season[k]) + '   ' + str(flux_season[k]) + '   ' + str(err_season[k]) + "\n")	                   	
		out.close()
	
	# all season cont file made below
	out=open(savepath+'all_seasons/all_seasons_cont_combined.txt', 'w')
	out.write('# WIRO ' + str(len(wiro_all_season_cont[0]))+'\n')
	for i in range(len(wiro_all_season_cont[0])):
        	out.write(str(wiro_all_season_cont[0][i])+ '   ' + str(wiro_all_season_cont[1][i])+ '   ' +str(wiro_all_season_cont[2][i])+'\n') 

	########### below pound commented out uses ztf g-band data from indiv seasons, not from season1 start all the way to last season last day

	#out.write('# G_band_ztf ' + str(len(g_all[0]))+'\n')
	#for i in range(len(g_all[0])):
        	#out.write(str(g_all[0][i])+ '   ' + str(g_all[1][i]) + '   ' + str(g_all[2][i]) + '\n')
	#keys = list(asas_all_dic.keys())
	#for key in keys:
		#cam = asas_all_dic.get(key)
		#out.write('# ASAS_{} '.format(key) + str(len(cam[0]))+'\n')
		#for i in range(len(cam[0])):
			#out.write(str(cam[0][i]) + '   ' + str(cam[1][i]) + '   ' + str(cam[2][i]) + '\n')                
	"""out.write('# C ' + str(len(r_all[0]))+'\n')
	for i in range(len(r_all[0])):
        	out.write(str(r_all[0][i])+ '   ' + str(r_all[1][i]) + '   ' + str(r_all[2][i]) + '\n')"""

	######## below creates all_seasons_combined.txt using ztf and ASAS-SN data across entire time period from day1 of season1 to last day of last season

	out.write('# G_band_ztf ' + str(len(g_all_meta[0]))+'\n')
	for i in range(len(g_all_meta[0])):
        	out.write(str(g_all_meta[0][i])+ '   ' + str(g_all_meta[1][i]) + '   ' + str(g_all_meta[2][i]) + '\n')		
	keys = list(asas_all_dic_meta.keys())
	for key in keys:
		cam = asas_all_dic_meta.get(key)
		out.write('# ASAS_{} '.format(key) + str(len(cam[0]))+'\n')
		for i in range(len(cam[0])):
			out.write(str(cam[0][i]) + '   ' + str(cam[1][i]) + '   ' + str(cam[2][i]) + '\n')
	out.close()
        
        
def writehb(filenames, jd, flux, err, savepath):
	"""Writes separate seasons of hb data to text files                                                                     

    	Parameters 
    	----------                                                                                     
    	filenmes : :class:`array-like`
    	List of names of text files to spit out e.g. 'season1'

    	indices : :class:`array-like` 
    	List of boolean arrays for each season, bounds defined in code by user

    	jd : :class:`array-like`
    	Array of reduced jd (-2457000)

    	flux : :class:`array-like`
    	Array of Hbeta fluxes

    	err : :class:`array-like`
    	Array of Hbeta flux uncertainities
        
        objpath : :class'string'
        Path to object directory

    	Returns
    	-------
    	:class:`.txt file`
    	Outputs are text files named with filename_hb.txt for each
    	season of data
        
        AJF NOTES:
        
        User does not need to change this function.
        
    	"""
	
	for i in range(len(filenames)):
		out = open(savepath+filenames[i]+'/'+filenames[i]+'_hb.txt', 'w')
		[out.write(str(jd[i][j]) + ',' + str(flux[i][j]) + ',' 
				+ str(err[i][j]) +  '\n') for j in range(len(jd[i]))]
		out.close()

def good_seeing(ra, dec, start, end, jd, tol=0.05):
        """Checks ztf data for seeing and selects data with seeing below
        some set limit.

	Parameters
    	----------
    	ra : :class:`float`
    	RA of object in decimal degrees (GET FROM NED)
        
    	dec : :class:`float'
    	Dec of object in decimal degrees (GET FROM NED)

    	start : :class:`string'
    	Date of start of WIRO observations in form yyyy-mm-dd
	
        end : :class:`string'
        Date of end of WIRO observations in form yyyy-mm-dd
        
        jd : :class:`array-like'
    	Array of jd dates for ztf data 
        
        tol : :class:`float'
    	Minimum value of the difference between same dates. 
        Default = 0.05
        
    	Returns                                                                
    	-------                                                                                        
    	:class:`array-like`
    	An array of indices that selects observations with good seeing.
        
        AJF NOTES:
        
        User does not need to change this function.
    	
        """
        
        jdstart = Time(start).jd
        jdend = Time(end).jd
        
        zquery = query.ZTFQuery()
        zquery.load_metadata(kind="sci", radec=[ra, dec], size=0.01, sql_query=f"seeing<1.4 and obsjd BETWEEN {jdstart} AND {jdend}")

        data = zquery.metatable

        jd_goodseeing = data.obsjd
        jd_goodseeing -= 2457000

        ii = [np.argsort(abs(jd-val)< tol) for val in jd_goodseeing]
        index = np.concatenate(ii, axis=None)
        index = np.unique(index)

        return index
        

def create_dir(num):
	for i in range(int(num)):
		if os.path.isdir(f'../pycali/season{int(i)+1}'): 
			break                
		else:
			os.mkdir(f'../brains/season{int(i)+1}')
			os.mkdir(f'../nar_sub/season{int(i)+1}')
			os.mkdir(f'../integrated/season{int(i)+1}')
			os.mkdir(f'../pycali/season{int(i)+1}')
			
	if not os.path.isdir(f'../pycali/all_seasons'):
		os.mkdir(f'../pycali/all_seasons')
                
def main():
    	# TZ initialize parser
	parser = argparse.ArgumentParser(description='Reads in ztf data, converts to flux units, splits into seasons, and writes each season in the form for cali.py')
	parser.add_argument('objpath', metavar='object_path', help='path to separate object directory (probably "original"')
	parser.add_argument('savepath', metavar='save_path', help='path to directory where you want continuum and h-beta files saved (probably "pycali"')
	parser.add_argument('season_num', metavar = 'season_num', help = 'number of total seasons')
	parser.add_argument('ra', metavar='obj_ra', help='object ra in decimal degrees')
	parser.add_argument('dec', metavar='obj_dec', help='object dec in decimal degrees')
	parser.add_argument('arcsearch', metavar='arc_search', help='search in ztf-query around the ra/dec in a radius of this many arcsec (3 is pretty good)')
	arg = parser.parse_args()
        
	create_dir(arg.season_num)
	
	# check to make sure h-beta peak is in limits for all survey filters        		
	param_path = arg.objpath+'parameters'
	parameters = open(param_path).readlines()
	z = [float(i.split()[1]) for i in parameters if i.split()[0] == 'z'][0]
	hb_1 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'hb_1'][0]
	hb_2 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'hb_2'][0]
	hb_peak = (((hb_1+hb_2)/2)*z)+((hb_1+hb_2)/2)
	ag_omit_l = 3390
	ag_omit_h = 6210           
	av_omit_l = 4630
	av_omit_h = 6390
	zg_omit_l = 3970
	zg_omit_h = 5613
	zr_omit_l = 5500
	zr_omit_h = 7394
	if hb_peak > ag_omit_h or hb_peak< ag_omit_l:
		que = input('\nYour h_beta peak may be outside of asas-sn g-band limits. Still proceed?\n')                        
		if que not in ['yes', 'y', 'Yes', 'Y']:
			sys.exit()     
	if hb_peak > av_omit_h or hb_peak< av_omit_l:
		que = input('\nYour h_beta peak may be outside of asas-sn V-band limits. Still proceed?\n')                        
		if que not in ['yes', 'y', 'Yes', 'Y']:
			sys.exit()  
	if hb_peak > zg_omit_h or hb_peak< zg_omit_l:
		que = input('\nYour h_beta peak may be outside of ztf g-band limits. Still proceed?\n')                        
		if que not in ['yes', 'y', 'Yes', 'Y']:
			sys.exit()     
	if hb_peak > zr_omit_h or hb_peak< zr_omit_l:
		que = input('\nYour h_beta peak may be outside of ztf r-band limits. If you are only using g-band, then you are set to proceed. If you are using r-band, consider speaking to Dr. Brotherton, Thea, or Aidan. Still proceed?\n')                        
		if que not in ['yes', 'y', 'Yes', 'Y']:
			sys.exit()
                                                
	asas_hjd, flux_asas, flux_err_asas, asas_cam = get_asas(arg.objpath) 	
	jd, flux, flux_err, filt = get_ztf(arg.ra, arg.dec, arg.objpath, arg.arcsearch)

	#ii = good_seeing(32.394162, 52.442475, '2018-08-23', '2022-03-14', jd) 

	g, asas_cs_dic, wiro, g_all, asas_all_dic, wiro_all_season_cont, cs_dic, g_all_meta, asas_all_dic_meta = split_seasons(jd, flux, flux_err, filt, arg.objpath, arg.season_num, asas_hjd, flux_asas, flux_err_asas, asas_cam)

	names = []
	for i in range(int(arg.season_num)):
        	names.append('season'+str(i+1))	
	write_cont(wiro, g, asas_cs_dic, arg.season_num, names, arg.savepath, g_all, wiro_all_season_cont, asas_all_dic, cs_dic, g_all_meta, asas_all_dic_meta)
	writehb(names,wiro[:,0],wiro[:,3],wiro[:,4], arg.savepath)


if __name__ == '__main__':
	main()


