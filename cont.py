#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 11:26:06 2023

@author: heitorernandes
"""

from astropy.io import fits
import os
import numpy as np
from scipy.io import readsav

import splice_ech

#test

spec='./data/Reduced_red/spec-c.ech'

hdul = fits.open(spec)



#mode = 'blue'
mode = 'red'



#RED------------------
if mode == 'red':
    directory = 'Reduced_red'
    files = [f for f in os.listdir() if f.endswith('.ech')]
    
    
    for filename in files:
        
        hdul = fits.open(spec)
        spec=hdul[1].data['SPEC']
        sig=hdul[1].data['SIG']
        cont=hdul[1].data['CONT']
        
        sig = sig * np.sqrt(cont)
        
        sav_data = readsav('data/'+directory+'/harps_red.ord_default.sav')
        
        blzcoef = sav_data.blzcoef
        col_range=sav_data.col_range
        
        
        #sav_data.keys()
        #Out[41]: dict_keys(['orders', 'or_range', 'ord_err', 'col_range', 'def_xwd', 'def_sxwd', 'blzcoef'])
        
        
        #ee = make_cont(e, blzcoef, colrange=col_range, param=[10, 8e5, 1e7, 1.], SCALING=True, PLOT=True,
                       #yr=[0, 10], wrange=[[5577, 5577.2], [5885, 5902], [6299.9, 6300.2]])
                       
                       
        #output_name = filename[:filename.find('.ech')] + 'c.ech'
        #save...
        #rdech(filename, RAW=True)
        #wdech(output_name, ee.head, ee.spec, sig=ee.sig, cont=ee.cont, wave=e.wave, orders=e.orders, OVERWRITE=True)


#BLUE-----------------
elif mode == 'blue':
    directory = 'Reduced_blue'
    os.chdir(directory)
    files = [f for f in os.listdir() if f.endswith('.ech')]

    for filename in files:
        e = rdech(filename, NOCONT=True)
        e.sig = e.sig * np.sqrt(e.cont)
        restore('harps_blue.ord_default.sav')
        
        ee = make_cont(e, blzcoef, colrange=col_range, param=[10, 1e5, 1e8, 1.], SCALING=True, PLOT=True,
                       yr=[0, 10], order_scales=order_scales, weights=weights)
        output_name = filename[:filename.find('.ech')] + 'c.ech'
        rdech(filename, RAW=True)
        wdech(output_name, ee.head, ee.spec, sig=ee.sig, cont=ee.cont, wave=e.wave, orders=e.orders, OVERWRITE=True)





















#
