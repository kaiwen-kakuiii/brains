import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import astropy.units as u
import argparse
import os

# script originally made by Thea Zastrocky, and then subsequently modified by Aidan Ferguson
# finished tidying up on 03/11/24 at 10:43 am in BRAINS/example_brains/software
# edited on 03/06/24 at 8:51 pm, added condition if flux goes negative
# edited 7/1/24 AJF
def readspec(filename1, filename2):
    """Fit and subtract narrow lines from single epoch spectra 

    Parameters
    ----------
    filename1 : :class:'string'
    name of combined.txt_season# file
    filename2 : :class:'string'
    name of nar_sub_combined.txt_season# file
    

    Returns
    -------
    jd: :class:'float':
    JD of epoch

    wave, flux, err: :class:'array-like'
    Arrays of narrow-line-subtracted wavelengths, fluxes, and errors 

    syserr: :class:'float':
    Systematic error of epoch
    
    AJF notes:
    intakes combined.txt_season# and nar_sub_combined.txt_season# lists, which contain all epoch
    names. un-narrow-line-subtracted files contain precaluclated systematic errors, which are accounted for in the lines immediately below.
    nar_sub_combined epoch files contain wavelength, flux, and flux error measurements after narrow-line subtraction, which are needed to 
    create the h-beta line file for both BRAINS and reverberation mapping
       
    """
    l = open(filename1).readlines()
    jd = float(l[0].split()[2])
    if len([float(i) for i in l[0].split('[')[-1].split(']')[0].split()]) > 1:
        #syserr = np.std([float(i) for i in l[0].split('[')[-1].split(']')[0].split()], ddof = 1)
        num_syserr = len([float(i) for i in l[0].split('[')[-1].split(']')[0].split()])
        syserr = np.std([float(i) for i in l[0].split('[')[-1].split(']')[0].split()], ddof = 1) / num_syserr**0.5
    else:
        syserr = 0.0
    
    file2 = np.loadtxt(filename2, dtype = np.float64, delimiter=',')
    wave, flux, err = file2[:,0], file2[:,1], file2[:,2]
    # if want to print with precision, use np.set_printoptions(precision = #)
    return jd, wave, flux, err, syserr


def hb_flux_func(wave, flux, err, syserr, param_path, name1, out3, em):
    """Measure Hb and 5100 fluxes (densities) from single epoch spectra; format brains h-beta file

    Parameters
    ----------
    wave : :class:'numpy-array'
    Array of wavelengths in the rest frame
    
    flux : :class:'numpy-array'
    Array of spectrum fluxes, narrow-line subtracted

    err : :class:'numpy-array'
    Array of errors on flux values

    syserr : :class:'float'
    Systematic error for the epoch

    param_path : :class:'string'
    objectpath1

    Returns
    -------
    wave_hb: class: float
    wavelengths for brains (observed frame)
    
    fhb: class: float
    fluxes for each epoch, used for brains h-beta file
    
    efhb: class: float
    errors of brains h-beta fluxes
    
    wave_hb_int: class: float
    wavelengths of integrated fluxes (rest frame)
    
    fhb_int: :class:'float'
    Hbeta narrow line subtracted integrated flux 

    efhb_int: :class:'float'
    Hbeta narrow line subtracted integrated flux errors 
    
    AJF notes:
    intakes parameter file and wave, flux, err, and syserr from readspec above. uses these to subtract off linear continuum right "under"
    h-beta line, fully isolating the broad line component, and then both writes this data to the brains h-beta line profile file and also
    integrates this data for use with normal reverberation mapping
    """

    parameters = open(param_path+"parameters").readlines()
    z = [float(i.split()[1]) for i in parameters if i.split()[0] == 'z'][0] 
    con_l_1 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'con_l_1'][0]
    con_l_2 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'con_l_2'][0]
    con_r_1 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'con_r_1'][0]
    con_r_2 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'con_r_2'][0]
    hb_1 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'hb_1'][0]
    hb_2 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'hb_2'][0]

    # TZ put wavelengths in rest frame
    # wave /= (1+z)
    # all wavelengths imported from date_combined.txt.nar_sub files are already in rest frame! -AJF 01-31-2024

    # TZ get left continuum window
    index = np.where((wave >= con_l_1) & (wave < con_l_2))
    flux0 = flux[index[0]]
    mediancl = np.median(flux0)
    wavecl = np.mean(wave[index[0]])

    # TZ get right continuum window
    index = np.where((wave >= con_r_1) & (wave < con_r_2))
    flux1 = flux[index[0]]
    mediancr = np.median(flux1)
    emediancr = np.std(flux1) / (float(len(flux1)))**0.5
    emediancr = (emediancr**2 + (syserr * mediancr)**2)**0.5
    wavecr = np.mean(wave[index[0]])

    # TZ fit linear continuum between the right and left continuum regions and subtract from the spectrum
    index = np.where((wave >= hb_1) & (wave < hb_2))
    wave_hb = wave[index[0]]
    fagn = flux[index[0]]
    efagn = err[index[0]]
    f = interpolate.interp1d(np.array([wavecl, wavecr]), np.array([mediancl, mediancr]))
    fcon = f(wave_hb)
     
    # AJF subtract off linear continuum fit
    fhb = fagn - fcon #subtracts off lin. cont.
    wave_hb_int = wave_hb #wavelengths in rest frame for integration    
    fhb_int = np.sum(fhb)*(wave_hb_int[1] - wave_hb_int[0]) # integrates
       
    # if flux goes negative, add offset to make it at least zero
    if len([*filter(lambda x: x<0, fhb)]) > 0:
    	print(f'H-beta flux is negative after linear continuum fit subtraction for file: {name1}, {em}')
    	print(f'Integrated flux before using offset: {fhb_int}')
    	print(f'Min of H-beta flux array is: {min(fhb)}')
    	out3.write('H-beta flux is negative after linear continuum fit subtraction for file: {}, {}\n'.format(name1, em))
    	out3.write('Integrated flux before using offset: {}\n'.format(fhb_int))
    	out3.write('Min of H-beta flux array is: {}\n'.format(min(fhb)))    	
    	offset = min(fhb)*(-1)
    	fhb += offset
    	fhb_int = np.sum(fhb)*(wave_hb_int[1] - wave_hb_int[0])
    	print(f'Integrated flux after using offset: {fhb_int}')    
    	print(f'Min of H-beta flux array is now: {min(fhb)}')
    	print('\n')
    	out3.write('Integrated flux after using offset: {}\n'.format(fhb_int))    
    	out3.write('Min of H-beta flux array is now: {}\n'.format(min(fhb)))
    	out3.write('\n')             
    wave_hb *= (1+z) #puts in observed frame which brains wants
    efhb = efagn
    
    #AJF and TZ integrate fluxes for reverberation mapping
    efhb_int = np.sum(efagn**2)**0.5 * (wave_hb_int[1]-wave_hb_int[0])
    efhb_int = (efhb_int**2 + (syserr*fhb_int)**2)**0.5

    return wave_hb, fhb, efhb, wave_hb_int, fhb_int, efhb_int, out3 # first three outputs used for brains hb file, next three used for integration (RM, CCF)


def main():
    """
    AJF notes:
    main() function runs everything above it; sort of like a large loop; see "if __name__=='__main__'" documentation online 
    """
    parser = argparse.ArgumentParser(description='Calculates narrow line subbed hb fluxes')
    parser.add_argument('objpath1', metavar='objectpath1', help='Path to object directory for narrow_line_subbed and original epoch files (probably path to "test" directory)')
    parser.add_argument('savepath', metavar='savepath', help='path to save brains formatted line profile files to (probably path to "brains" directory)')
    parser.add_argument('savepath2', metavar='savepath2', help='Path to save integrated flux to (probably path to "integrated" directory)')
    parser.add_argument('season', metavar='season_name', help='Season you are working on')
    parser.add_argument('-em', '--ems', metavar='ems', type = str, nargs = '+', help = 'error methods to be used')
    arg = parser.parse_args()
    ems = arg.ems

    # TZ read in epoch names
    os.system(f'cp ../original/* {arg.objpath1}')
    os.system(f'cp ../nar_sub/{arg.season}/* {arg.objpath1}')
    """
    df1 = pd.read_csv(arg.objpath1+'combined.txt_'+arg.season, delim_whitespace=True, header=None)
    names1 = np.array(df1[0])
    df2 = pd.read_csv(arg.objpath1+'nar_sub_combined.txt_'+arg.season, delim_whitespace=True, header=None)
    names2 = np.array(df2[0])
    
    # AJF check nb automatically; check manually by counting number of wavelength points between hb1 and hb2 if needed
    parameters = open(arg.objpath1+"parameters").readlines()
    hb_1 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'hb_1'][0]
    hb_2 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'hb_2'][0]
    check_date = str(names2[0])
    check_nb = np.array((pd.read_csv(arg.objpath1+check_date, delimiter=',', header=None))[0])
    idex = np.where((check_nb >= hb_1) & (check_nb < hb_2))
    nb = len(check_nb[idex[0]]) 
    ne = len(names1)
    
    # TZ and AJF start text files for brains and integrated fluxes
    out = open(arg.savepath+arg.season+'/'+'brains_hb_'+arg.season+'.txt', 'w')
    out2 = open(arg.savepath2+arg.season+'/'+'integrated_'+arg.season+'_hb.txt', 'w')
    out3 = open(arg.savepath2+arg.season+'/cmdoutput_'+arg.season+'.txt', 'w')
    out.write('# {} {}\n'.format(ne,nb))

    for name1 in names1:
    	# TZ read in spectra
        indx = np.where(names1 == name1)
        indx_val = indx[-1][-1]
        name2 = names2[indx_val]
        name1 = arg.objpath1+name1
        name2 = arg.objpath1+name2
        jd, wave, flux, err, syserr = readspec(name1,name2)
        jd = jd- 2457000
        out.write('# {}\n'.format(jd))
        
    	# TZ calculate hb fluxes with narrow lines subtraced
        wave_hb, hb_flux, hb_err, wave_hb_int, fhb_int, efhb_int, out3 = hb_flux_func(np.array(wave), np.array(flux), np.array(err), syserr, 
            arg.objpath1, name1, out3)
        for j in range(len(wave_hb)):
        	out.write('{} {} {}\n'.format(wave_hb[j], hb_flux[j], hb_err[j]))
        out.write('\n')
        
        # AJF write jd, integrated fluxes and flux error to integrated text file
        out2.write('{},{},{}\n'.format(jd, fhb_int, efhb_int))
    out.close()
    out2.close()
    out3.close()
    """
    ###################################### now do for m1, m2, m3, m4 ###########################
    
    # TZ read in epoch names
    for em in ems:
    	df1 = pd.read_csv(arg.objpath1+'combined.txt_'+arg.season, delim_whitespace=True, header=None)
    	names1 = np.array(df1[0])
    	df2 = pd.read_csv(arg.objpath1+'nar_sub_combined.txt_'+arg.season+'_'+em, delim_whitespace=True, header=None)
    	names2 = np.array(df2[0])
    
    	# AJF check nb automatically; check manually by counting number of wavelength points between hb1 and hb2 if needed
    	parameters = open(arg.objpath1+"parameters").readlines()
    	hb_1 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'hb_1'][0]
    	hb_2 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'hb_2'][0]
    	check_date = str(names2[0])
    	check_nb = np.array((pd.read_csv(arg.objpath1+check_date, delimiter=',', header=None))[0])
    	idex = np.where((check_nb >= hb_1) & (check_nb < hb_2))
    	nb = len(check_nb[idex[0]]) 
    	ne = len(names1)
    
    	# TZ and AJF start text files for brains and integrated fluxes
    	out = open(arg.savepath+arg.season+'/'+'brains_hb_'+arg.season+'.txt_'+em, 'w')
    	out2 = open(arg.savepath2+arg.season+'/'+'integrated_'+arg.season+'_hb.txt_'+em, 'w')
    	out3 = open(arg.savepath2+arg.season+'/cmdoutput_brains_and_int_'+arg.season+'.txt_'+em, 'w')
    	out.write('# {} {}\n'.format(ne,nb))

    	for name1 in names1:
    		# TZ read in spectra
        	indx = np.where(names1 == name1)
        	indx_val = indx[-1][-1]
        	name2 = names2[indx_val]
        	name1 = arg.objpath1+name1
        	name2 = arg.objpath1+name2
        	jd, wave, flux, err, syserr = readspec(name1,name2)
        	jd = jd- 2457000
        	out.write('# {}\n'.format(jd))
        
    		# TZ calculate hb fluxes with narrow lines subtraced
        	wave_hb, hb_flux, hb_err, wave_hb_int, fhb_int, efhb_int, out3 = hb_flux_func(np.array(wave), np.array(flux), np.array(err), syserr, arg.objpath1, name1, out3, em)
        	for j in range(len(wave_hb)):
        		out.write('{} {} {}\n'.format(wave_hb[j], hb_flux[j], hb_err[j]))
        	out.write('\n')
        
        	# AJF write jd, integrated fluxes and flux error to integrated text file
        	out2.write('{},{},{}\n'.format(jd, fhb_int, efhb_int))
    	out.close()
    	out2.close()
    	out3.close()
	

if __name__=='__main__':
    main() 


