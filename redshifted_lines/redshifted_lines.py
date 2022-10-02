
'''
The wavelength_coverage_plot method creates a plot of spectral lines
redshifted into bands/filters of your choice. I made this with inspiration
from work by Taylor Hutchison (@aibhleog). 

This has the JWST/NIRSpec, MIRI, NIRCam WFSS, and NIRISS spectroscopic bands, 
along with Gemini/GMOS preloaded. You can add your own as well.

Credit: Nikko Cleri
		cleri@tamu.edu
		Texas A&M University
'''

_author_ = 'Nikko Cleri'

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PATH = '/Users/alvis/Research/jwst_cycle2_prep/redshifted_lines/'# include the '/' at the end!
PATH_FILTERS = PATH+'filters/' # include the '/' at the end!
PATH_FIGURES = PATH+'figures/' # include the '/' at the end!

# Dictionary of lines
# currently holds the lines I use most frequently in my work
lines = {'lybreak':[911.75,r'Lyman Limit $\lambda$911.75'],
         'lya':[1215.67,r'Ly$\alpha$ $\lambda$1215.67'],
        'nv':[1240,'NV $\lambda$1240'],
        'civ':[1548,'CIV $\lambda$1549'],
        'heii1':[1640.4,'HeII $\lambda$1640'],
        'oiii]':[1664,'OIII] $\lambda$1660,1666'],
        'siiii]':[1883,'SiIII] $\lambda$1883,1892'],
        'ciii]':[1906.8,'CIII] $\lambda$1907\n  & $\lambda$1909'],
        'mgii':[2798,'MgII $\lambda$2796\n  & $\lambda$2803'],
        '[nev]1':[3346,''],
        '[nev]2':[3426,'[NeV] $\lambda$3346\n    & $\lambda$3426'],
        'Hbreak':[3646,r'Balmer Limit $\lambda$3646'],
        '[oii]':[3727,'[OII] $\lambda$3727\n    & $\lambda$3729'],
        '[neiii]':[3869,'[NeIII] $\lambda$3869'],
        'hdelta':[4102,'H$\delta$ $\lambda$4102'],
        'hgamma':[4341,'H$\gamma$ $\lambda$4341'],
        'heii2':[4686,'HeII $\lambda$4686'],
        'hbeta':[4861, r'H$\beta$ $\lambda$4863'],
        '[oiii]1':[4959,''],
        '[oiii]2':[5007,'[OIII] $\lambda$4959\n    & $\lambda$5007'],
        'hei':[5876,'HeI $\lambda$5876'],
        'halpha':[6563, r'H$\alpha$ $\lambda$6563'],
        '[sii]1':[6716,''],
        '[sii]2':[6731,'[SII] $\lambda$6716\n    & $\lambda$6731'],
        'pabreak':[8250,r'Paschen Limit $\lambda$8250'],
        'pabeta':[12818,r'Pa$\beta$ $\lambda$12820'],
        'paalpha':[18751,r'Pa$\alpha$ $\lambda$18750'],
        'brgamma':[21661,r'Br$\gamma$ $\lambda$21661'],
        'brbeta':[26259,r'Br$\beta$ $\lambda$26259'],
        'pfepsilon':[30383,r'Pf$\epsilon$ $\lambda$30383'],
        'pfdelta':[32970,r'Pf$\delta$ $\lambda$32970'],
        'pfgamma':[37406,r'Pf$\gamma$ $\lambda$37406'],
        'bralpha':[40523,r'Br$\alpha$ $\lambda$40523'],
        'pfbeta':[46538,r'Pf$\beta$ $\lambda$46538'],
        'pfalpha':[74599,r'Pf$\alpha$ $\lambda$74599']
        }

linesdf = pd.DataFrame(lines)
linesdf.rename(index={0:'wave', 1:'name'}, inplace=True)
linesdf = linesdf.transpose()
linesdf.to_csv('lines.csv', index=False)

lst = ['-','--','-.',':']
hatching = ['/', '|', '-', '+', 'x', 'o', 'O', '.', '*']

# Reading in a bunch of different spectroscopic bandpasses
# JWST/MIRI, JWST/NIRISS, JWST/NIRCam, JWST/NIRSpec
miri_lrs = pd.DataFrame({'wave':[5e4,12e4],'throughput':[1,1]}) # angstroms
miri_mrs = pd.DataFrame({'wave':[5.9e4,28.1e4],'throughput':[1,1]}) # angstroms
niriss_wfss = pd.DataFrame({'wave':[0.8e4,2.2e4],'throughput':[1,1]}) # angstroms
niriss_soss = pd.DataFrame({'wave':[0.6e4,2.8e4],'throughput':[1,1]}) # angstroms
nircam_wfss = pd.DataFrame({'wave':[2.4e4,5e4],'throughput':[1,1]}) # angstroms
nirspec_prism = pd.DataFrame({'wave':[0.6e4,5.3e4],'throughput':[1,1]}) # angstroms
nirspec_mrs1 = pd.DataFrame({'wave':[0.7e4,1.27e4],'throughput':[1,1]}) # angstroms
nirspec_mrs2 = pd.DataFrame({'wave':[0.97e4,1.89e4],'throughput':[1,1]}) # angstroms
nirspec_mrs3 = pd.DataFrame({'wave':[1.66e4,3.17e4],'throughput':[1,1]}) # angstroms
nirspec_mrs4 = pd.DataFrame({'wave':[2.87e4,5.27e4],'throughput':[1,1]}) # angstroms
nirspec_mrs = pd.DataFrame({'wave':[0.6e4,5.27e4],'throughput':[1,1]})# angstroms
jwst_filts = [miri_lrs,miri_mrs,niriss_wfss,niriss_soss,nircam_wfss,nirspec_prism,nirspec_mrs1,nirspec_mrs2, nirspec_mrs3, nirspec_mrs4]
jwst_filts_names = ['MIRI LRS','MIRI MRS','NIRISS WFSS','NIRISS SOSS','NIRCam WFSS','NIRSpec Prism','NIRSpec MRS1','NIRSpec MRS2','NIRSpec MRS3','NIRSpec MRS4']

# Gemini GMOS
gmos = pd.DataFrame({'wave':[3600, 10300],'throughput':[1,1]}) # angstroms

# HST WFC3 grisms
g102 = pd.read_csv(PATH_FILTERS+'hst_g102_throughput.csv')
g141 = pd.read_csv(PATH_FILTERS+'hst_g141_throughput.csv')

def line_z_range(line, filter):
    return [min(filter['wave'])/lines[line][0]-1, max(filter['wave'])/lines[line][0]-1]

def wavelength_coverage_plot(z_min=0, z_max=20, lambda_min=0, lambda_max=30,
                             filters=[], filter_names=[],
                             show_lines=[],
                             grid=False,
                             lines_cmap='RdBu_r', filts_cmap='plasma',
                             filename='wavelength_coverage_plot.pdf'):
    
    # Selecting the lines of interest
    # this requires linesdf to be read in already
    sliced_lines = linesdf.loc[show_lines] # get a dataframe of just the lines we want
    cmap_lines = plt.get_cmap(lines_cmap)
    colors = [cmap_lines(j) for j in np.linspace(0,1,len(sliced_lines))]
    z = np.linspace(z_min,z_max,num=1000)
    
    # Filters and filter labels
    # requires filters to already be defined
    cmap_filts = plt.get_cmap(filts_cmap)
    filtcolors = [cmap_filts(j) for j in np.linspace(0.1,1,len(filters))]

    # Making the plot
    plt.figure(figsize=(20,11))
    for i in range(len(sliced_lines)):
        plt.plot((1+z)*sliced_lines['wave'][i]/1e4, z, color=colors[i], label=sliced_lines['name'][i],ls=lst[i%(len(lst))],lw=3)
    for i,f in enumerate(filters):
        plt.axvspan(f['wave'][0]/1e4,f['wave'][f['wave'].size-1]/1e4,color=filtcolors[i],alpha=0.2,zorder=0, label=filter_names[i])
    plt.legend(ncol=2, bbox_to_anchor=(1.3,0.5), loc='center')
    plt.xlabel('Observed Wavelength [microns]')
    plt.ylabel('Redshift')
    if grid:
        plt.grid()
    plt.axis([lambda_min,lambda_max,z_min,z_max])
    plt.savefig(PATH_FIGURES+filename)