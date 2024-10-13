# created originally by Thea Zastroxky, modified and edited by Aidan Ferguson
# updated 7/1/24

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lmfit
from lmfit import minimize, Parameters, fit_report
from scipy.optimize import curve_fit
from scipy.stats import skew
import astropy.units as u
import argparse
import os
import sys
from tqdm import tqdm
import corner


def subtract_narrow_lines(wave, flux, err, param_path, savepath, season, names, low, high, out3, q_mcmc):
    """
    Parent function subtract_narrow lines intakes original combined epoch data and fits Gaussians to
    each narrow line component, icluding the OIII lines and the narrow line h-beta component. It also
    fits the AGN power law, subtracting off most of the continuum.
    """
    
    # TZ read in continuum, OIII, hb integration windows and redshift
    
    parameters = open(param_path).readlines()
    z = [float(i.split()[1]) for i in parameters if i.split()[0] == 'z'][0] 
    hb_1 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'hb_1'][0] 
    hb_2 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'hb_2'][0]
    oiii_1 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'oiii_1'][0]
    oiii_2 = [float(i.split()[1]) for i in parameters if i.split()[0] == 'oiii_2'][0]
    oiii_4959l = [float(i.split()[1]) for i in parameters if i.split()[0] == 'oiii4959_l'][0]
    oiii_4959r = [float(i.split()[1]) for i in parameters if i.split()[0] == 'oiii4959_r'][0]
    
    # TZ convert all wavelengths to observed frame and define bounds
    oiii5007_step1_lim1 = oiii_1
    oiii5007_step1_lim2 = oiii_2
    hb_lim1 = hb_1 
    hb_lim2 = hb_2    
    hb_wave_lims = (wave>=(hb_1)) & (wave<=(hb_lim2))
    iii = hb_wave_lims
    # subtract continuum background; define fit function for AGN power law
    ii = (wave >= (low)) & (wave <= (high))

    
    # fit continuum
    np.set_printoptions(formatter={'float':lambda x: '{}\n'.format(x)})           
    def resi_agn_cont(params, wave, flux, err, return_model=False):
        c = params['c']
        a = params['a']
        model = (c*(wave**a))
        if return_model == True:
            return model
        return (flux - model) / err

    params_agn_cont = Parameters()

    params_agn_cont.add('a', value = -0.2)
    params_agn_cont.add('c', value = 1e-13)

    out_agn_cont= minimize(resi_agn_cont, params_agn_cont, args = (wave[ii], flux[ii], err[ii]), scale_covar = False)
    
    agn_cont = resi_agn_cont(out_agn_cont.params, wave, flux, err, True) 
    
    
    # AJF error prop   #############################################################################################################################################
    
    # method 1: no mcmc, just minimize and covariance on the minimize least-sq fit
    #print(f'this is original error, no modifications: {err[iii]}')
    def error_prop_1(covar, wave_1, a_1, c_1):
        df_da = c_1 * np.log(wave_1) * (wave_1**a_1)
        df_dc = wave_1**a_1
        unc_a = covar[0][0]
        unc_c = covar[1][1]
        covar_val = covar[1][0]
        #print(df_da, df_dc, unc_a, unc_c, covar_val)
        under = np.abs( (unc_a * df_da**2) + (unc_c * df_dc**2) + (2 * df_da * df_dc * covar_val) )
        y_err = np.sqrt(under) 
        return y_err
    err_covar_m1 = error_prop_1(out_agn_cont.covar, wave[iii], out_agn_cont.params['a'].value, out_agn_cont.params['c'].value)
    #print(f'this is err on model using method 1, covariances:\n {err_covar_m1}')
    err_m1 = np.sqrt((err[iii]**2+err_covar_m1**2))
    #print(f'this is final error after propogation with original error, using method 1:\n {err_m1}')
    #print(f'this is report fit for method 1: {lmfit.fit_report(out_agn_cont.params)}')
    #print(f'this is covariance for a,c for method 1: {out_agn_cont.covar}\n')
    
    
    # method 2: use mcmc to create sample space of resi_agn_cont parameters a and c; then, use every single set of a and c to produce a corresponding model power law. 
    # then, take standard deviation of column vector which corresponds to stdev of each wavelength for the full set of model fluxes
    if q_mcmc == 'y':
        print(f'\nMethods 2 and 3, mcmc, progress below:')
        burn = 150
        steps = 1000
        thin = 1
        mcmc = minimize(resi_agn_cont, method = 'emcee', nan_policy = 'omit', burn = burn, steps = steps, thin = thin, params = params_agn_cont, progress = True, args = (wave[ii], flux[ii], err[ii]), run_mcmc_kwargs={'skip_initial_state_check':True})
        df_a = mcmc.flatchain['a']
        df_c = mcmc.flatchain['c']
        x1 = np.arange(1, len(df_a)+1, 1)
    
        # plot mcmc results
        fig = plt.figure(figsize = (20, 10))
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        ax1.plot(x1, df_a, 'k-', linewidth = 0.03)
        ax1.set_title(f'MCMC sample of parameter a, with burn {burn}, steps {steps} and thin {thin}')
        ax2.plot(x1, df_c, 'b-', linewidth = 0.03)
        ax2.set_title(f'MCMC sample of parameter c, with burn {burn}, steps {steps} and thin {thin}')
        ax3.plot(mcmc.acceptance_fraction, 'o')
        ax3.set_xlabel('walker')
        ax3.set_ylabel('acceptance fraction')
        plt.tight_layout()
        plt.savefig(savepath+'/'+season+'/'+names[0:8]+'_mcmc_sample.png', format = 'png')
        plt.show()


        # method 2 mcmc here
        mc_arr = np.zeros((len(df_a),len(wave[iii]))) #rows, columns
        check_std = np.zeros((len(df_a)))
    
        def cont_num(a_num, c_num, wave_num, flux_num):
            model_num = (c_num*(wave_num**a_num))
            return model_num
        print(f'\nAssignment to array progress below:')
        for i in tqdm(range(len(df_a))):
            out = cont_num(df_a[i], df_c[i], wave[iii], flux[iii])
            mc_arr[i,:] = out
        
        for i in range(len(df_a)):
            check_std[i] = mc_arr[i, -1]
            
        mc_err = np.zeros(len(wave[iii]))
        mc_err_quantile = np.zeros(len(wave[iii]))
        test_skew = np.zeros(len(wave[iii]))
        mc_err_std = np.zeros(len(wave[iii]))
    
        length = len(wave[iii])
    
        #cols = 5
        #rows = length // cols
        #if length % cols !=0:
	        #rows +=1
        #position = range(1,length+1)

        #fig = plt.figure(figsize = (40,60))
    
        print(f'\nAssignment to array and testing skewness progress below:')
        for i in tqdm(range(length)):
            test_skew[i] = skew(mc_arr[:,i])
            low_err, high_err = np.quantile(mc_arr[:,i], q = ( (1.0-0.683)/2, 1-(0.317/2) ))
            mc_err_quantile[i] = (high_err - low_err)/2.0
            mc_err_std[i] = np.std(mc_arr[:,i])
            if test_skew[i] > 2 or test_skew[i] < -2:
                mc_err[i] = mc_err_quantile[i]
            else:
                mc_err[i] = mc_err_std[i]
        
            #ax = fig.add_subplot(rows, cols, position[i])
            #ax.hist(mc_arr[:,i], density = False, bins = 100, color = 'black')
            #ax.text(0.05, 1.02, wave[iii][i], transform=plt.gca().transAxes, fontsize = 9)
        #plt.show()

        #print(f'this is skew: {test_skew}')    
        #print(f'this is model error using method 2, mcmc create set of fluxes, 1 stdevs:\n {mc_err_std}')
        #print(f'this is model error using method 2, mcmc create set of fluxes, 1 stdevs, quantile:\n {mc_err_quantile}')
        err_m2 = np.sqrt(mc_err**2 + err[iii]**2)
        #print(f'this is final error after propogation with original error, using method 2 (mcmc set of fluxes):\n {err_m2}')
    
        # method 3 mcmc here (covariances)
        print(f'\nthis is fit report for both methods 2 and 3 (mcmc):\n')
        lmfit.report_fit(mcmc.params) # prints results of mcmc
        chq = 1
        while chq == 1:
            try:      
                input2 = float(input(f'\nWhat is the number above in C(a,c)? (should be between -1 and 1):\n'))
                out3.write('\nCorrelation is the number above in C(a,c), should be between -1 and 1:\n')
                out3.write('{}\n'.format(input2))
                chq = 0
                break
            except ValueError:
                print('Try entering the number again.')
                out3.write('\nTry entering a number again\n')    
    
        cov_val = input2 * np.sqrt( mcmc.params['a'].stderr**2 * mcmc.params['c'].stderr**2)
        #print(f'this is covariance for a,c for method 3: {cov_val}')
        #print(f'compared to covariance for a,c for method 1: {out_agn_cont.covar[1][0]}')
    
        mcmc_covar = np.array( [ [mcmc.params['a'].stderr**2,cov_val] , [cov_val,mcmc.params['c'].stderr**2] ] )
        def error_prop_1(covar, wave_1, a_1, c_1):
            y_err = np.sqrt( np.abs( (covar[0][0]*(c_1*np.log(wave_1)**wave_1**a_1)**2) + (covar[1][1]*wave_1**(2*a_1)) + (2*c_1*covar[1][0]*np.log(wave_1)*wave_1**(2*a_1)) ) )
            return y_err
        err_covar_m3 = error_prop_1(mcmc_covar, wave[iii], mcmc.params['a'].value, mcmc.params['c'].value)
        #print(f'this is err on model using method 3, mcmc covariances:\n {err_covar_m3}')
        err_m3 = np.sqrt((err[iii]**2+err_covar_m3**2))
        #print(f'this is final error after propogation with original error, using method 3:\n {err_m3}')
    
    # method 4 bootstrapping here, uses results from method 1
    def boot_method(waves, a, c, a_err, c_err):
        a_ran = np.random.uniform((-1*np.abs(a_err)), np.abs(a_err))
        c_ran = np.random.uniform((-1*np.abs(c_err)), np.abs(c_err))
        new_a = a + a_ran
        new_c = c + c_ran
        mod = new_c * (waves**new_a)
        #print(f'\n\na_ran: {a_ran}, c_ran: {c_ran}, a_err: {a_err}, c_err: {c_err}, a: {a}, c: {c}, new_a: {new_a}, new_c: {new_c}\n')
        return mod
    
    num = 100000
    m4_mod_arr = np.zeros((num, len(wave[iii])))
    m4_std = np.zeros(len(wave[iii]))    
    print(f'\nmethod 4 bootstrap progress below:')
    for i in tqdm(range(num)):
        mod = boot_method(wave[iii], out_agn_cont.params['a'].value, out_agn_cont.params['c'].value, np.sqrt(out_agn_cont.covar[0][0]), np.sqrt(out_agn_cont.covar[1][1]))
        m4_mod_arr[i,:] = mod
           
    for i in range(len(wave[iii])):
        m4_std[i] = np.std(m4_mod_arr[:,i])
          
    #print(f'this is method 4 error on model, m4_std: {m4_std}')
    
    err_m4 = np.sqrt((err[iii]**2 + m4_std**2))
    #print(f'this is final error after prop with og error, model 4: {err_m4}')
    
    # explanations here    
    
    print(f'\n\nSummary: Method 0 is no error propagation. Method 1 uses covariances calculated from minimize least-sq-resid method. Method 2 uses mcmc to calculate a set of model functions, then finds the standard')
    print(f'deviation or quantile measurement of each wavelength bin across all flux samples (i.e. takes column of fluxes corresponding to single wavelength) and uses this as error')
    print(f'Method 3 uses mcmc to explore parameter space, find best value for parameters, and errors on those parameters as well as correlation coeff., which is used to find')
    print(f'covariance between a,c, and error is calculated same way as method 1, just using mcmc to find best parameters/errors rather than least-sq-resid. Method 4 uses the')
    print(f'bootstrapping technique of finding a random value between the limits of error for each a,c pair found in method 1, and adding this random value ot the best parameter')
    print(f'value, creating a new a,c pair, which, in turn, creates a new model perterbed by the random error added on. It does this {num} times, creating a model space from which stdev can be drawn.\n')
    
    
    # AJF set errors:
    err_original = np.zeros(len(err))
    
    m1_err = np.zeros(len(err))
    m2_err = np.zeros(len(err))
    m3_err = np.zeros(len(err))
    m4_err = np.zeros(len(err))
    
    # sets all 4 methods error as just original error
    for i in range(len(err)):
        err_original[i] = err[i]
        m1_err[i] = err[i] 
        m4_err[i] = err[i]
        if q_mcmc == 'y':
            m2_err[i] = err[i]
            m3_err[i] = err[i]
            
    # sets method error in H-beta limits to derived errors from propogation  
    m1_err[iii] = err_m1
    m4_err[iii] = err_m4
    
    if q_mcmc == 'y':
        m2_err[iii] = err_m2
        m3_err[iii] = err_m3

      
    #####################################################################################################################################################################   
   
    # subtract cont   
     
    flux2 = flux - agn_cont   
    
    # check loop   
    out3.write('\nYou have just subtracted {}, date: {} with low and high continuum limits: {},{} with {} mcmc\n'.format(season, names, low, high, q_mcmc))
    print(f'\nYou have just subtracted {season}, date: {names} with low and high continuum limits: {low},{high} with {q_mcmc} mcmc.\n')
    
    if len([*filter(lambda x: x<0, flux2[iii])]) > 0:
    	offset = (min(flux2[iii])*-1)+((min(flux2[iii])*-1)/4)
    	agn_cont = agn_cont - offset
    	print(f'This was minimum of flux inside h-beta limits: {min(flux2[iii])}')
    	out3.write('This was minimum of flux inside h-beta limits: {}\n'.format(min(flux2[iii])))     
    	flux -= agn_cont
    	print(f'Used offset in power law subtraction')
    	out3.write('Used offset in power law subtraction\n')
    	print(f'This is now the minimum of flux inside h-beta limits: {min(flux[iii])}')
    	out3.write('This is now the minimum of flux inside h-beta limits: {}\n'.format(min(flux[iii])))
    else:
    	flux = flux2 	   
       	
       
    # step1: subtract [OIII]4959 to decontaminate 5007
    # A linear background is added below [OIII]4959, 
    # and a small shift is determined as a free parameter.
    # AJF - OIII lines occur at 4959 A and 5007 A
    
    # read in spectra data
    indexoiii = (wave >= oiii5007_step1_lim1) & (wave <= oiii5007_step1_lim2)
    wave_oiii = wave[indexoiii]
    flux_oiii = flux[indexoiii]
    err_oiii = err[indexoiii]

    # define fitting functions
    def resi_narrow(params, wave_narrow, flux_narrow, err_narrow, sub_linear_cont=True):
        wavelength = params['wavelength']
        p0 = params['peak']
        p1 = params['shift'] 
        p2 = params['sigma']
        gauss_model = p0*np.exp(-0.5*(wave_narrow-(wavelength+p1))**2/p2**2)
        if sub_linear_cont == True:
            a = params['a']
            b = params['b']
            gauss_model += (a*wave_narrow + b)
   
        return (flux_narrow - gauss_model) / err_narrow

    def narrow_line(params, wave_narrow, sub_linear_cont=True):
        wavelength = params['wavelength']
        peak = params['peak']
        shift = params['shift']
        sigma = params['sigma']
        gauss_model = peak*np.exp(-0.5*(wave_narrow-(wavelength+shift))**2/sigma**2)
        if sub_linear_cont == True:
            a = params['a']
            b = params['b']
        return gauss_model
    
    # define OIII-5007 line parameters
    params_5007 = Parameters()
    params_5007.add('wavelength', value = 5007, vary=False)
    params_5007.add('a', value = 0)
    params_5007.add('b', value= 0)
    params_5007.add('peak', value = 1e-15, max = 1e-12, min=1e-16)
    params_5007.add('shift', value = 1.0) # TZ center shift
    params_5007.add('sigma', value = 10.0)  # sigma


    # TZ fit 5007 line
    out_5007 = minimize(resi_narrow, params_5007, args = (wave_oiii, flux_oiii, err_oiii))
   
    profile_5007 = np.zeros(len(flux))
    profile_5007[indexoiii] = narrow_line(out_5007.params, wave_oiii, False)
    flux -= profile_5007
    
    # AJF prop. error ################################################
    def narrow_line_err(params, wave_narrow, covar):   
        wavelength = params['wavelength']
        peak = params['peak']
        shift = params['shift']
        sigma = params['sigma']
        
        dpeak = (np.exp(-0.5*(wave_narrow-(wavelength+shift))**2/sigma**2))
        dwl = (1/sigma**2)*peak*(wave_narrow-wavelength-shift)*np.exp((-0.5*(wave_narrow-wavelength-shift)**2 )/sigma**2 )
        dshift = (1/sigma**2)*peak*(wave_narrow-wavelength-shift)*np.exp((-0.5*(wave_narrow-wavelength-shift)**2 )/sigma**2 )
        dsigma = peak*(wave_narrow-wavelength-shift)**2*(1/sigma**3)*np.exp(( (-0.5*(wave_narrow-wavelength-shift)**2))) 
        
        model_unc_cov = np.sqrt( (covar[2][2]*dpeak**2) + (covar[3][3]*dshift**2) + (covar[4][4]*dsigma**2) + 2 * ( (covar[2][3] * dpeak* dshift) + (covar[2][4]*dpeak*dsigma) + (covar[4][3]*dsigma*dshift) )  )   # ADD IN COVARIANCES HERE

        return model_unc_cov   
    
    ##################################################################
    
      
    # TZ integrate the profile to get 5007 flux density
    flux_5007 = np.sum(profile_5007[indexoiii])*(wave_oiii[1] - wave_oiii[0])

    # TZ now fit 4959 using best fit params > need to sub off continuum
    index4959 = (wave >= oiii_4959l) & (wave <= oiii_4959r)
    
    # define OIII-4959 line parameters
    params_4959 = Parameters()
    params_4959.add('wavelength', value = 4959, vary=False)
    params_4959.add('peak', value = out_5007.params['peak'].value/2.98)
    params_4959.add('shift', value = out_5007.params['shift'].value, vary=False) # TZ center shift
    params_4959.add('sigma', value = out_5007.params['sigma'].value, vary=False)  # sigma
    params_4959.add('a', value = 1e-17)
    params_4959.add('b', value= 1e-15)

    # TZ fit 4959 line
    out_4959 = minimize(resi_narrow, params_4959, args = (wave[index4959], flux[index4959], 
        err[index4959]))

    profile_4959 = np.zeros(len(flux))
    profile_4959[index4959] = narrow_line(out_4959.params, wave[index4959], False)
    flux -= profile_4959
       
    # TZ integrate the profile to get 4959 flux density
    flux_4959 = np.sum(profile_4959[index4959])*(wave[index4959][1] - wave[index4959][0])

    # TZ sub off narrow component of Hb
    index_hb = (wave >= hb_lim1) & (wave <= hb_lim2)
    wave_hb = wave[index_hb]
    flux_hb = flux[index_hb]
    err_hb = err[index_hb]
    
    # defne narrow line H-beta line parameters
    params_narrowhb = Parameters()
    params_narrowhb.add('wavelength', value = 4861, vary=False)
    params_narrowhb.add('peak', value = out_5007.params['peak'].value/10, vary=False)
    params_narrowhb.add('shift', value = out_5007.params['shift'].value, vary=False) # TZ center shift
    params_narrowhb.add('sigma', value = out_5007.params['sigma'].value, vary=False)  # sigma
    params_narrowhb.add('a', value = 1e-17)
    params_narrowhb.add('b', value= 1e-16, max=1e-16)
    
    # fit the narrow h-beta line
    out_narrowhb = minimize(resi_narrow, params_narrowhb, args = (wave_hb, flux_hb, 
        err_hb, False))
    
    profile_narrowhb = np.zeros(len(flux))
    profile_narrowhb[index_hb] = narrow_line(out_narrowhb.params, wave_hb, False)
    flux -= profile_narrowhb
     
    """ 
    # AJF prop error for narrow line #############################################################################################################################################################
    # method 1, covariances
         
    unc_hb = np.zeros(len(err_flux_m1))
    unc_hb = narrow_line_err(out_narrowhb.params, wave_hb, out_5007.covar)
    #print(f'this is hb_err narrow line fit using covariances: {unc_hb}')

    
    # method 2: use mcmc to create sample space of resi_agn_cont parameters a and c; then, use every single set of a and c to produce a corresponding model power law, then subtract each power law from the original
    # epoch's flux to get a set of length equal to the number of models (len(df_a)) of power-law-subtracted fluxes. then, take standard deviation of column vector which corresponds to stdev of
    # each wavelength for the full set of subtracted fluxes
    
    burn = 1500
    steps = 5000
    thin = 5
    mcmc = minimize(resi_narrow, method = 'emcee', nan_policy = 'omit', burn = burn, steps = steps, thin = thin, params = params_5007, is_weighted = True, progress = True, args = (wave_oiii, flux_oiii, err_oiii))
    df_p = mcmc.flatchain['peak']
    df_sh = mcmc.flatchain['shift']
    df_si = mcmc.flatchain['sigma']
    x1 = np.arange(1, len(df_p)+1, 1)
    
    fig = plt.figure(figsize = (13, 9))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    ax1.plot(x1, df_p, 'k-')
    ax1.set_title(f'MCMC sample of parameter peak, with burn {burn}, steps {steps} and thin {thin}')
    ax2.plot(x1, df_sh, 'b-')
    ax2.set_title(f'MCMC sample of parameter shift, with burn {burn}, steps {steps} and thin {thin}')
    ax3.plot(x1, df_si, 'g-')
    ax3.set_title(f'MCMC sample of parameter sigma, with burn {burn}, steps {steps} and thin {thin}')
    plt.show()
    
    mc_arr2 = np.zeros((len(df_p),len(wave[iii]))) #rows, columns
    check_std2 = np.zeros((len(df_p)))
    
    def narrow_num(peak, shift, sigma, wave_narrow, flux_num):
        gauss_model = peak*np.exp(-0.5*(wave_narrow-(5007+shift))**2/sigma**2)
        diff = flux_num - gauss_model
        return diff

    for i in tqdm(range(len(df_p))):
        out = narrow_num(df_p[i], df_sh[i], df_si[i], wave[iii], flux[iii])
        if i == 0 or i==3 or i==5:
            print(f'this is df_p, df_sh, df_si: {df_p[i]}, {df_sh[i]}, {df_si[i]}')
        mc_arr2[i,:] = out
    
    print(f'test of flux values for one set of params: {mc_arr2[0,:]}')
    print(f'test of flux values for one set of params: {mc_arr2[3,:]}')
    print(f'test of flux values for one set of params: {mc_arr2[5,:]}')
    print(f'last flux of last one: {mc_arr2[5,-1]}')    
    for i in range(len(df_p)):
        check_std2[i] = mc_arr2[i, -1]
    print(f'last wavelength set of fluxes: {check_std2[0:40]}')        
    mc_err2 = np.zeros(len(wave[iii]))
    for i in range(len(wave[iii])):
        mc_err2[i] = np.std(mc_arr2[:,i])
        
    print(f'this is error using method 2 (MCMC), 3 stdevs:\n {3*mc_err2}')
    #print(f'this is final error using method 2 for HBeta limits {wave[iii][0]} to {wave[iii][-1]}, 3 stdevs: {3*mc_err}')
    
  
    # AJF test errors:
    #print(f'this is og_err, err_flux_m1, 3*mc_err, and unc_hb_m1: {err[iii]}\n\n{err_flux_m1}\n\n{3*mc_err}\n\n{unc_hb}')
    err_og_mcmc_covar =  np.sqrt((3*mc_err)**2 + (unc_hb)**2) # original error is converted to new error through mcmc process, then covar of narrow line added in quad 
    err_og_mcmc_covar_diff =  np.sqrt(err[iii]**2 + (3*mc_err)**2 + (unc_hb)**2) # original error is added to error through mcmc process, then covar of narrow line added in quad 
    err_og_covar_covar = np.sqrt(err_flux_m1**2 + unc_hb**2) # original error is added in quad to covar err from agn_cont, then covar of narrow line added in quad
    err_og_none_covar = np.sqrt(err[iii]**2 + (unc_hb)**2) # original error (no agn_cont error added) is added in quad to covar of narrow line 
    err_og_covar_none = err_flux_m1 # original error (no hb_line error added) is added in quad to covar of agn_cont line
    #print(f'this is err_og_mcmc_covar: {err_og_mcmc_covar}')
    #print(f'this is err_og_mcmc_covar_diff: {err_og_mcmc_covar_diff}')
    #print(f'this is err_og_covar_covar: {err_og_covar_covar}')
    #print(f'this is err_og_none_covar: {err_og_none_covar}') 
    #print(f'this is err_og_covar_none: {err_og_covar_none}')
    
    #print(f'this is unc_hb: {unc_hb[index_hb]}')
    print(f'this is err2 hb after hb: {err2[index_hb]}')
    print(f'this is wave after hb: {wave_hb}')
    
    # AJF set error
    err = err2
        
    #############################################################################################################################################################
    """
    
    # TZ integrate the profile to get narrow line h-beta flux density
    flux_nhb = np.sum(profile_narrowhb[index_hb])*(wave_hb[1] - wave_hb[0])


    # TZ fitting is done. Need to plot nicely to show fits
    # TZ need to assess goodness of fit


    fig = plt.figure(figsize = (13, 9))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
 

    # TZ plot agn power continuum
    ax1.plot(wave,flux+agn_cont+profile_5007+profile_4959+profile_narrowhb, label='data')
    ax1.plot(wave, agn_cont, '--', label='power continuum, fitted to linear continuum model')
    ax1.plot(wave, flux+profile_5007+profile_4959+profile_narrowhb, label='data w/ power cont subtracted')
    #ax1.plot(wave[ii], agn_lin, 'k', label = 'linear model fit to continuum region', alpha = 0.75, linestyle = (0, (5,3)))
    ax1.legend()
    ax1.set_title('AGN power continuum')

    cont_5007 = out_5007.params['a'].value*wave_oiii+out_5007.params['b'].value
    cont_4959 = out_4959.params['a'].value*wave[index4959]+out_4959.params['b'].value
    cont_hb = out_narrowhb.params['a'].value*wave_hb+out_narrowhb.params['b'].value

    # TZ now need to plot everything nicely
    # TZ first plot oiii fit 
    ax2.plot(wave, flux+profile_5007+profile_4959, label='data')
    ax2.plot(wave[index4959], profile_4959[index4959]+cont_4959, '--', color='purple', label='4959 fit')
    ax2.plot(wave[index4959], cont_4959, '--', color='purple')
    ax2.plot(wave[indexoiii], profile_5007[indexoiii]+cont_5007, '--', color= 'green', label='5007 fit')
    ax2.plot(wave_oiii, cont_5007, '--', color='green')
    ax2.set_xlim(wave[index4959][0]-10, wave_oiii[-1]+10)
    ax2.legend()
    ax2.set_title('O III fits')

    # TZ now plot Hbeta narrow line
    ax3.plot(wave, flux, label='data w/ narrow hb subbed')         
    ax3.plot(wave, flux+profile_narrowhb, '--', label='data')
    ax3.plot(wave, profile_narrowhb, '--', label='narrow hb fit')
    ax3.set_xlim(4700, 5050) 
    ax3.legend()
    ax3.set_title(f'Hb narrow line fit')
    
    #plot Date AJF
    ax4.text(0.2, 0.5, r'Date: ' + names[4:6] + '-' + names[6:8] + '-' + names[0:4], size = 20, transform = ax4.transAxes)
    ax4.axis('off')
    
    # save and show fitted figure "date_sub.pdf"
    plt.savefig(savepath+'/'+season+'/'+names[0:8]+'_sub.pdf', format = 'pdf')
    plt.show()
    #plt.close()


    # TZ finally plot subtracted spectrum and residuals, in a separate figure
    # AJF make it zoomed on hb line
    fig2 = plt.figure(figsize = (10, 6))
    axspec = fig2.add_subplot(211)
    axspec.plot(wave, flux)
    axspec.set_ylabel('Flux')
    axspec.axhline(y=0, color = 'r', linestyle = 'dotted', label = 'zero-point reference line')
    axspec.vlines(x=[hb_lim1, hb_lim2], ymin = min(flux), ymax = max(flux), color = 'orange', linestyle = 'dotted', label = 'h_beta line upper/lower limits')
    axspec.set_xlim(hb_lim1-20, hb_lim2+20)
    range1 = (wave>(hb_lim1-20)) & (wave<(hb_lim2+20))
    axspec.set_ylim((min(flux[range1])-np.std(flux[range1])), (max(flux[range1])+np.std(flux[range1])))
    axspec.legend(loc='best')
    axspec.grid(which = 'both', axis = 'y', alpha=0.5)
    
    # add subplot plain
    axspec2 = fig2.add_subplot(212)
    axspec2.plot(wave, flux, label = 'Nar. Line Sub. Spectra')
    axspec2.set_xlabel('Wavelength')
    axspec2.set_ylabel('Flux')
    axspec2.legend(loc='best')
    # AJF add date
    axspec2.text(0.5, 0.9, r'Date: ' + names[4:6] + '-' + names[6:8] + '-' + names[0:4], size=10, ha='center', va='center', transform = axspec2.transAxes)


    # AJF save and show figure as date_spec.pdf
    plt.savefig(savepath+'/'+season+'/'+names[0:8]+'_spec.pdf', format = 'pdf')
    plt.show()
    #plt.close()
    check2 = input('Are you happy with this result?\n')
    out3.write('Are you happy with this result?\n')
    out3.write('{}\n'.format(check2))

    return flux, agn_cont, err, profile_5007, profile_4959, profile_narrowhb, flux_5007, flux_4959, flux_nhb, check2, out3, m1_err, m2_err, m3_err, m4_err


def main():
    
    # TZ initialize parser
    parser = argparse.ArgumentParser(description='Fits a continuum and narrow line components of AGN spectra and subtracts the fit.') 
    parser.add_argument('objpath', metavar='obj_path', help='path the object folder')
    parser.add_argument('savepath', metavar='save_path', help='path to place you want to save the output spectra')
    parser.add_argument('season', metavar='season', help='name of season')
    parser.add_argument('low', metavar = 'low', type = int, help = 'lower limit of the power law fit')
    parser.add_argument('high', metavar = 'high', type = int, help = 'upper limit of the power law fit')
    parser.add_argument('mcmc', metavar = 'mcmc', help = 'would you like to run mcmc sampling (y or n)')
    arg = parser.parse_args()
    np.set_printoptions(threshold=sys.maxsize)
    q_mcmc = arg.mcmc
    
    # TZ read in epoch names
    df = pd.read_csv(arg.objpath+'combined.txt_'+arg.season, delim_whitespace=True, header=None)
    names = np.array(df[0])

    # TZ read in spectra from individual epochs
    spec = [pd.read_csv(arg.objpath+name, delim_whitespace=True, header=None, skiprows=2) for name in names]
    wave = np.array(spec[0][0]) # same wavelengths for all dates (bc convolved) so only need one array
    flux = [np.array(x[1]) for x in spec]
    err = [np.array(x[2]) for x in spec]
    savepath = arg.savepath

    # TZ put wavelengths in rest frame, AJF intake parameter file for z and do de-redshift for all combined files before loop
    param_path = arg.objpath+'parameters'
    parameter_path = param_path
    parameters = open(param_path).readlines()
    z = [float(i.split()[1]) for i in parameters if i.split()[0] == 'z'][0]
    
    # AJF convert to rest frame
    wave /= (1+z)
    
    # begin cmdoutput.txt
    out3 = open(arg.savepath+'/'+arg.season+'/'+'cmdoutput_'+arg.season+'.txt', 'w')
    check = 'y'  
    if q_mcmc == 'y':
        out3.write('\nYou have chosen to include mcmc sampling in propagating the AGN power law.\n')
        print(f'\nYou have chosen to include mcmc sampling in propagating the AGN power law.\n')
    else:
        out3.write('\nYou have chosen to NOT include mcmc sampling in propagating the AGN power law.\n')
        print(f'\nYou have chosen to NOT include mcmc sampling in propagating the AGN power law.\n')
    for i in range(len(flux)):               
        # start subtracting
        low = arg.low
        high = arg.high
        print("="*75)
        out3.write('{}'.format("="*75))   
        while check in ['y']:
            out3.write('\nSubtracting...\n')
            flux_narrowsub, agn_cont, err_final, profile_5007, profile_4959, profile_narrowhb, flux_5007, flux_4959, flux_nhb, check2, out3, m1_err, m2_err, m3_err, m4_err = subtract_narrow_lines(wave, flux[i], err[i], parameter_path, arg.savepath, arg.season, names[i], low, high, out3, q_mcmc)
        
            if check2 in ['y', 'yes', 'Y', 'Yes']: 
        
                # write narrow subtracted fluxes and wavelengths to files
                
                # use original error
                out = open(arg.savepath+'/'+arg.season+'/'+names[i]+'.nar_sub_m0', 'w')
                for j in range(len(wave)):
                    out.write(str(wave[j]) + ',' + str(flux_narrowsub[j]) + ',' + str(err_final[j]) +'\n')
                out.close()
                
                # use err method 1
                out_m1 = open(arg.savepath+'/'+arg.season+'/'+names[i]+'.nar_sub_m1', 'w')
                for j in range(len(wave)):
                    out_m1.write(str(wave[j]) + ',' + str(flux_narrowsub[j]) + ',' + str(m1_err[j]) +'\n')
                out_m1.close()      
                
                if q_mcmc == 'y':          
                    # use err method 2
                    out_m2 = open(arg.savepath+'/'+arg.season+'/'+names[i]+'.nar_sub_m2', 'w')
                    for j in range(len(wave)):
                        out_m2.write(str(wave[j]) + ',' + str(flux_narrowsub[j]) + ',' + str(m2_err[j]) +'\n')
                    out_m2.close()     
                           
                    # use err method 3
                    out_m3 = open(arg.savepath+'/'+arg.season+'/'+names[i]+'.nar_sub_m3', 'w')
                    for j in range(len(wave)):
                        out_m3.write(str(wave[j]) + ',' + str(flux_narrowsub[j]) + ',' + str(m3_err[j]) +'\n')
                    out_m3.close()              
                
                # use err method 4
                out_m4 = open(arg.savepath+'/'+arg.season+'/'+names[i]+'.nar_sub_m4', 'w')
                for j in range(len(wave)):
                    out_m4.write(str(wave[j]) + ',' + str(flux_narrowsub[j]) + ',' + str(m4_err[j]) +'\n')
                out_m4.close()                
        
                # write list of narrow subtracted files
                if i == 0:
                    out2 = open(arg.savepath+'/'+arg.season+'/'+'nar_sub_combined.txt_' + arg.season+'_m0', 'w')
                    for k in range(len(names)):
                        out2.write(str(names[k]+'.nar_sub_m0' + '\n'))
                    out2.close()
                    
                    out2_m1 = open(arg.savepath+'/'+arg.season+'/'+'nar_sub_combined.txt_' + arg.season+'_m1', 'w')
                    for k in range(len(names)):
                        out2_m1.write(str(names[k]+'.nar_sub_m1' + '\n'))
                    out2_m1.close()
                    
                    out2_m4 = open(arg.savepath+'/'+arg.season+'/'+'nar_sub_combined.txt_' + arg.season+'_m4', 'w')
                    for k in range(len(names)):
                        out2_m4.write(str(names[k]+'.nar_sub_m4' + '\n'))
                    out2_m4.close()
                    
                    if q_mcmc == 'y':
                        out2_m2 = open(arg.savepath+'/'+arg.season+'/'+'nar_sub_combined.txt_' + arg.season+'_m2', 'w')
                        for k in range(len(names)):
                            out2_m2.write(str(names[k]+'.nar_sub_m2' + '\n'))
                        out2_m2.close()
                    
                        out2_m3 = open(arg.savepath+'/'+arg.season+'/'+'nar_sub_combined.txt_' + arg.season+'_m3', 'w')
                        for k in range(len(names)):
                            out2_m3.write(str(names[k]+'.nar_sub_m3' + '\n'))
                        out2_m3.close()
                    
                # counter for terminal print
                count = int(i)+1
                count_down = int(len(flux))-count
                print('\n')
                print(f'You have {count_down} days left.')
                print('\n')
                out3.write('You have {} days left.\n'.format(count_down))
                
                break
            else:
                # TZ read in spectra from individual epochs
                wave = np.array(spec[0][0]) # same wavelengths for all dates (bc convolved) so only need one array
                flux = [np.array(x[1]) for x in spec]
                err = [np.array(x[2]) for x in spec]

                # TZ put wavelengths in rest frame, AJF intake parameter file for z and do de-redshift for all combined files before loop
                parameters = open(param_path).readlines()
                z = [float(i.split()[1]) for i in parameters if i.split()[0] == 'z'][0]
    
                # AJF convert to rest frame
                wave /= (1+z)
                
                # start over
                print(f'Choose new low and high continuum fit limits. Min is {min(wave)}, max is {max(wave)}.')
                print('\n')
                chq = 1
                while chq == 1:
                    try:
                        low = int(input('What is the lower limit for the continuum fit you would like (in Angstroms)?\n'))
                        chq = 0
                        break
                    except ValueError:
                        print(f'Not a valid number, enter again.')
                        out3.write('\nNot a valid number, try again.')
                chq = 1
                while chq ==1:
                    try:
                        high = int(input('What is the upper limit for the continuum fit you would like (in Angstroms)?\n'))
                        chq = 0
                        break
                    except ValueError:
                        print(f'Not a valid number, enter again.')
                        out3.write('\nNot a valid number, try again.\n')
                out3.write('Choose new low and high continuum fit limits. Min is {}, max is {}.\n'.format(min(wave), max(wave)))
                out3.write('What is the lower limit for the continuum fit you would like (in Angstroms)?\n')
                out3.write('{}\n'.format(low))
                out3.write('What is the upper limit for the continuum fit you would like (in Angstroms)?\n')
                out3.write('{}\n'.format(high))
   
    out3.close()

    # TZ writes individual narrow line profiles to text files
    '''    
    out = open(savepath+'/narrow_line_profiles.txt', 'w')
    out.write('wave,agn_cont,profile_5007,profile_4959,profile_narrowhb\n')
    for i in range(len(agn_cont)):
        out.write('{},{},{},{},{}\n'.format(wave[i],agn_cont[i],profile_5007[i],profile_4959[i],profile_narrowhb[i]))
    out.close()

    out = open(savepath+'/ratio_narrow_lines.txt', 'w')
    out.write('5007/4959  5007/narrow_hb\n')
    out.write('{:.2f}  {:.2f}'.format(flux_5007/flux_4959, flux_5007/flux_nhb))
    out.close()
    '''
    

if __name__=='__main__':
    main() 


# AJF modified "names" keyword to display date on plots, and axspec to display date, and saved plots. also modified when de-redshift occured so that
# it happened in "main"
# 01/18/24 dated
