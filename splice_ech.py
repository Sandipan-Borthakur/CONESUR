#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 11:37:02 2023

@author: heitorernandes
"""

import numpy as np
from scipy.interpolate import CubicSpline


def splice(ech, wave, spec, blaz, index, sig=None, ORDERS=None, COLRANGE=None,
               WEIGHTS=None, SCALING=None, ORDER_SCALES=None, DEBUG=None, YRANGE=None, WRANGE=None):
    
    
    # ech much be the one after the wave solution
    
    # Sanity check for 'SPEC'
    if not 'SPEC' in ech.columns.names :
        print('Ilegitimate ech structure: must include SPEC tag')
        raise SystemExit

    # Sanity check for 'WAVE'
    if not 'WAVE' in ech.columns.names :
        print('splice_ech expects the wavelength solution to be included in the ech structure')
        raise SystemExit

    # Sanity check for 'CONT' (blaze functions)
    if not 'CONT' in ech.columns.names:
        print('splice_ech expects the blaze functions to be included in the ech structure as CONT')
        raise SystemExit

    # Additional optional parameters
    
    if sig is not None:
        # Sanity check for 'SIG'
        if not 'SIG' in ech.columns.names :
            has_sig = 1 
    else:
        has_sig = 0


    npix = ech.data.dtype['SPEC'].shape[1]  # Order length in pixels
    nord = ech.data.dtype['SPEC'].shape[0]  # Number of spectral orders
    
    print("Lenth of orders in the ech :" +str(npix))
    print("Number of orders in the ech :" +str(nord))
    
    
    # creating the weights 
    
    weights = np.zeros((npix, nord)) + 1

    if COLRANGE is not None:
        sz = np.shape(COLRANGE)
        #check the shape of the COLRANGE must be two, the begining and the end of each order pixel
        if sz[1] != 2:
            print('COLRANGE should match the nstructure of a spectral order, infor on begin and end pixel')
            print('Help:\n', COLRANGE)
            raise SystemExit
        #check if the number of orders are correct
        if sz[0] != nord:
            print('COLRANGE should match the number of spectral orders')
            print('Help:\n', COLRANGE)
            raise SystemExit
        colr = COLRANGE
    else:
        colr = np.zeros((2, nord), dtype=int)
        colr[0, :] = 0
        colr[1, :] = npix - 1





    if ORDERS is not None:
        # If a subset of orders was specified
        sp = ech.data['SPEC'][:, ORDERS]
        ww = ech.data['WAVE'][:, ORDERS]
        bb = ech.data['CONT'][:, ORDERS]
    
        if sig is not None:
            unc = ech.data['SIG'][:, ORDERS]

        # DOIT WHY DOING IT AGAIN?
        npix = len(sp[0, 0])  # Order length in pixels
        nord = len(sp[0, :])   # Number of spectral orders

        weights = np.zeros((npix, nord)) + 1

    else:
        sp = ech.data['SPEC']
        ww = ech.data['WAVE']

        if SCALING is not None:
            bb = ech.data['CONT'] > 1.
            

            for iord in range(nord):
                i0 = colr[0, iord]
                i1 = colr[1, iord]

                # Calculate scale
                scale = np.median(ech.data['SPEC'][i0:i1, iord]) / np.median(np.median(ech.data['CONT'][i0:i1, iord], axis=0), axis=0)

                bb[i0:i1, iord] = np.median(ech.data['CONT'][i0:i1, iord], axis=0) * scale
                sp[i0:i1, iord] = ech.data['SPEC'][i0:i1, iord]

        else:
            bb = ech.data['CONT'] > 1.

            for iord in range(nord):
                i0 = colr[0, iord]
                i1 = colr[1, iord]
                bb[i0:i1, iord] = np.median(ech.data['CONT'][i0:i1, iord], axis=0)

        if sig is not None:
            unc = ech.data['SIG']


#sp == spec; bb == cont; unc== sig;

#;==================================================================================

     #OK
    order_scales = np.ones(nord)
    order_overlap = -np.ones((6, nord))

    # Find the order with the largest signal #OK
    # DOIT check if it is right 
    dd=sp / (bb > 0.1)
    dummy = np.max(np.median(dd[0], axis=1))
    iord0 = np.argmax(np.median(dd[0], axis=1))


    beg1 = colr[0, iord0]
    end1 = colr[1, iord0]
    
    w1 = ww[beg1:end1, iord0]
    s1 = sp[beg1:end1, iord0]
    b1 = bb[beg1:end1, iord0]
    
    if has_sig==1:
        sig1 = unc[beg1:end1, iord0]


    if iord0 > 0:
        for iord in range(iord0 - 1, 0, -1):
            beg0 = beg1  # Shift current order to previous
            end0 = end1
            w0 = w1
            s0 = s1
            b0 = b1
            if has_sig:
                sig0 = sig1
            
            beg1 = colr[0, iord]  # New current order
            end1 = colr[1, iord]
            w1 = ww[beg1:end1, iord]
            s1 = sp[beg1:end1, iord]
            b1 = bb[beg1:end1, iord]
            if has_sig:
                sig1 = unc[beg1:end1, iord]

            i0 = np.where((w0 >= np.min(w1)) & (w0 <= np.max(w1)))[0]  # Overlap within the previous order
            ii0 = np.arange(len(i0))
            
            
            #setting the WRANGE
            if WRANGE is not None :  # If exclusion regions are given, check the overlaps
                for iwrange in range(len(WRANGE)//2):
                    iii0 = np.where((w0[i0[ii0]] < WRANGE[0, iwrange]) or (w0[i0[ii0]] > WRANGE[1, iwrange]))[0]
                    if len(iii0) > 0:
                        ii0 = ii0[iii0]
                    
                    else:
                        ni0 = 0
                        
    return 0
            #








if __name__ == '__main__':

    print('Splice_ech test...')
    



























#