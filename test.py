#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 16:17:35 2023

@author: heitorernandes
"""


from astropy.io import fits
import os
import numpy as np
from scipy.io import readsav

from splice_ech import splice



if __name__ == '__main__':
    # 
    spec='./data/Reduced_red/spec-c.ech'

    ech = fits.open(spec)
    #running test just for red arm
    #RED------------------
    directory = 'Reduced_red'

        
    ech = fits.open(spec)
    spec=ech[1].data['SPEC']
    sig=ech[1].data['SIG']
    cont=ech[1].data['CONT']
    wave=ech[1].data['CONT']
        
    sig = sig * np.sqrt(cont)
        
    sav_data = readsav('data/'+directory+'/harps_red.ord_default.sav')
        
    blzcoef = sav_data.blzcoef
    col_range=sav_data.col_range
            
   
    # test splice_ech
    splice(ech[1], wave, spec, blzcoef, index=0, COLRANGE=list(col_range))
    
    print(col_range)
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
#