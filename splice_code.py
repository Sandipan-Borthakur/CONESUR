#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 11:39:23 2024

@author: heitor
"""


from astropy.io import fits
import os
import numpy as np
from scipy.io import readsav
import matplotlib.pyplot as plt
from splice_ech import splice
from numpy.polynomial.polynomial import Polynomial


import pickle as pkl

import pyreduce

import warnings
warnings.filterwarnings('ignore')


def poly_fit(x, y, degree=1):
    p = Polynomial.fit(x, y, degree)
    return p.convert().coef



#devoloper test, splice not as function







#calling the splice 

#*c-ech with wave solution 
spec='./data/goodech-Blue.fits'

#spec='./data/Reduced_red/spec.ech'


ech = fits.open(spec)
#running test just for red arm
#RED------------------
directory = 'Reduced_blue'

    
ech = fits.open(spec)
print('============================')

print('HEADER')
print('============================')
print(ech[0].header)
print('============================')



spec=ech[1].data['SPEC']
sig=ech[1].data['SIG']

cont=ech[1].data['CONT']
wave=ech[1].data['WAVE']
    
sig = sig * np.sqrt(cont)
    
sav_data = readsav('data/'+directory+'/harps_blue.ord_default.sav')
    
blzcoef = sav_data.blzcoef
col_range=sav_data.col_range



#original funciton and the way it should be called 
#def splice(ech, wave, spec, blaz, index, sig=None, ORDERS=None, COLRANGE=None,
                #  WEIGHTS=None, SCALING=None, ORDER_SCALES=None, DEBUG=None, YRANGE=None, WRANGE=None):
                    
#splice(ech[1], wave, spec, blzcoef, index=0, COLRANGE=list(col_range))





#BUG fixing the the wrong number of col_rang ERROR
aa = np.zeros(shape=(45, 2))
for i in np.linspace(0, len(col_range)-1, num=len(col_range)):
    i=int(i)
    aa[i]=col_range[i]
    aa[i][1]=int(aa[i][1]+1)

    
col_range=aa.astype(int)






#defining the calling parameter 
ech=ech[1]
wave=wave[0]
spec=spec[0]
blaz=blzcoef
index=0

sig=None
ORDERS=None

COLRANGE=list(col_range)

WEIGHTS=None

SCALING=None
ORDER_SCALES=None
DEBUG=None
YRANGE=None
WRANGE=None






#begin of the Splice "function as a code to make it easy for debuging 



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
    print('Applying Sigma...')
    # Sanity check for 'SIG'
    if not 'SIG' in ech.columns.names :
        has_sig = 1 
else:
    has_sig = 0


# defining the number of orders to user and the length in pixel of these orders


npix = ech.data.dtype['SPEC'].shape[1]  # Order length in pixels
nord = ech.data.dtype['SPEC'].shape[0]  # Number of spectral orders

# TEST
#nord=26

print("Lenth of orders in the ech :" +str(npix))
print("Number of orders in the ech :" +str(nord))


# creating the weights 

weights = np.zeros((nord, npix)) + 1 

if COLRANGE is not None:
    print('Applying Col Range...')
    sz = np.shape(COLRANGE)
    #check the shape of the COLRANGE must be two, the begining and the end of each order pixel
    if sz[1] != 2:
        print('COLRANGE should match the structure of a spectral order, infor on begin and end pixel')
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
    colr[1, :] = npix - 1 # CHECK if it need -1 or -0





if ORDERS is not None:
    

    print('Applying Orders...')
    # If a subset of orders was specified
    # the [0] is to fix the format of the ech
    sp = ech.data['SPEC'][0][:, ORDERS]
    ww = ech.data['WAVE'][0][:, ORDERS]
    bb = ech.data['CONT'][0][:, ORDERS]

    if sig is not None:
        unc = ech.data['SIG'][0][:, ORDERS]

    # DOIT WHY DOING IT AGAIN?
    npix = len(sp[0, 0])  # Order length in pixels
    nord = len(sp[0, :])   # Number of spectral orders

    weights = np.zeros((nord, npix)) + 1  # 45 rows (one for each order) and each order has 4096 pixels

    
    
else:
    # the [0] is to fix the format of the ech
    sp = ech.data['SPEC'][0]
    ww = ech.data['WAVE'][0]

    if SCALING is not None:
        print('Applying Scailing...')
        #bb = ech.data['CONT'][0] > 1. #original in IDL
        bb = ech.data['CONT'][0]
        # b = np.where(blaze<1, 1, blaze) maybe use that
        bb[bb<1]=1
        print('Scailing')
        

        for iord in range(nord):
            i0 = colr[0, iord]
            i1 = colr[1, iord]

            # Calculate scale
            scale = np.median(ech.data['SPEC'][0][iord][i0:i1]) / np.median(np.median(ech.data['CONT'][0][iord][i0:i1], axis=0), axis=0)

            bb[i0:i1, iord] = np.median(ech.data['CONT'][0][iord][i0:i1], axis=0) * scale
            sp[i0:i1, iord] = ech.data['SPEC'][0][iord][i0:i1]

    else:
        bb = ech.data['CONT'][0]
        bb[bb<1]=1

        for iord in range(nord):
        
            i0 = colr[iord][0]  # python way to say i0=[0,iord] which means in idl the first element (0) of the iorder (the row)
            i1 = colr[iord][1]
           
            #bb[i0:i1, iord] = np.median(ech.data['CONT'][i0:i1, iord], axis=0)

    if sig is not None:
        print('Applying sig...')
        unc = ech.data['SIG']


#sp == spec; bb == cont; unc== sig;

#================================================================================

print('Beginning the Splice lower orders')
order_scales = np.ones(nord)
order_overlap = -np.ones((6, nord))

# Find the order with the largest signal #OK
# DOIT check if it is right 
bb[bb<0.1]=0.1
dd=sp / (bb)


signal_orders = np.median(dd,axis=1)
signal_orders=list(signal_orders)

iord0=signal_orders.index(max(signal_orders))

#Sanity check the code works until here

beg1 = colr[iord0][0]
end1 = colr[iord0][1]

print('BUGG-------')
print(colr[iord0][1])
print('BUGG-------')


w1 = ww[iord0][beg1:end1]
s1 = sp[iord0][beg1:end1]
b1 = bb[iord0][beg1:end1]

if has_sig==1:
    sig1 = unc[iord0][beg1:end1]


#for plotting all the orders


#ax.set_xlim([ww[-1][0], ww[0][-1]])


#for plotting the overlap
fig1, ax1 = plt.subplots(figsize=(12, 5))



#plt.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])

#=============================================================================

print('Starting backwards loop')
if iord0 > 0:
    for iord in np.flip(np.linspace(0, iord0-1, num=iord0)): #start with the biggest signal order and goes n-1 n-2 ...
        iord=int(iord)    
        beg0 = beg1  # Shift current order to previous
        end0 = end1
        w0 = w1
        s0 = s1
        b0 = b1
        
        
        if has_sig:
            sig0 = sig1
        
        beg1 = colr[iord][0]  # New current order
        end1 = colr[iord][1]
        
        w1 = ww[iord][beg1:end1]
        s1 = sp[iord][beg1:end1]
        b1 = bb[iord][beg1:end1]
        
        
        if has_sig:
            sig1 = unc[iord][beg1:end1]
            
            
        #Defining the overlap region in i0 and ni0 where i0 is an array of the index (pixels) and ni0 its length
        i0_pix = np.where((w0 >= np.min(w1)) & (w0 <= np.max(w1)))[0] # Overlap within the previous order 
        i0 = ((w0 >= np.min(w1)) & (w0 <= np.max(w1)))    #MASK
        
        ni0=len(i0_pix)                     # Overlap within the previous order ()
        ii0 = np.arange(len(i0))
        
        
        
        #setting the WRANGE
        ## WARNING w0[i0] does not work, it should be a mask! not another array. np.take can work 
        #np.take(w0,i0)
       # if WRANGE is not None :  # If exclusion regions are given, check the overlaps
        #    for iwrange in range(len(WRANGE)//2):
         #       iii0 = np.where((w0[i0[ii0]] < WRANGE[0, iwrange]) or (w0[i0[ii0]] > WRANGE[1, iwrange]))[0]
          #      if len(iii0) > 0:
           #         ii0 = ii0[iii0]
                # DOIT check that sctructure =>  WRANGE[0, iwrange]
            #    else:
             #       ni0 = 0
                    
                    
        #Overlap with the current order based in the same logic as above
        i1_pix = np.where((w1 >= np.min(w0)) & (w1 <= np.max(w0)))[0] 
        i1 = ((w1 >= np.min(w0)) & (w1 <= np.max(w0)))    #MASK
        
        ni1=len(i1_pix)                    
        ii1 = np.arange(len(i1))
         
         
        #setting the WRANGE
     #   if WRANGE is not None :  # If exclusion regions are given, check the overlaps
      #      for iwrange in range(len(WRANGE)//2):
       #         iii1 = np.where((w1[i1[ii1]] < WRANGE[0, iwrange]) or (w1[i1[ii1]] > WRANGE[1, iwrange]))[0]
        #        if len(iii1) > 0:
         #           ii1 = ii1[iii1]
                 
       #         else:
        #            ni1 = 0      
        
        
        
        if (ni0 > 0) and (ni1 > 0):
            print('We have overlaping orders')
            
            #IDL original
            #tmpS0=bezier_interp(w1,s1,bezier_init(w1,s1,/DOUBLE),w0[i0],/DOUBLE)
            tempS0=pyreduce.util.bezier_interp(w1,s1,w0[i0])
            
            tempB0=pyreduce.util.bezier_interp(w1,b1,w0[i0])
            
            if has_sig==1:
                tempU0=pyreduce.util.bezier_interp(w1,sig1,w0[i0])

           

            tempS1=pyreduce.util.bezier_interp(w0,s0,w1[i1])
            
            tempB1=pyreduce.util.bezier_interp(w0,b0,w1[i1])
            
            if has_sig==1:
                tempU1=pyreduce.util.bezier_interp(w0,sig0,w1[i1])

        
            # add scaling 
            if SCALING is not None:
                
                scl0 = np.sum(s0[i0[ii0]]/b0[i0[ii0]])/np.sum(tempS0[ii0]/tempB0[ii0])
                scl1 = np.sum(s1[i1[ii1]]/b1[i1[ii1]])/np.sum(tempS1[ii1]/tempB0[ii1])
                
                scl=np.sqrt(scl0/scl1)
                s1=s1*scl
                if has_sig==1:
                    tempS0=tempS0*scl
                    tempU0=tempU0*scl
                    order_scales[iord]=scl
            #end scaling
             
            if (ww[iord][0]) > (ww[iord+1][0]):
                print('Current order (1) is bluer than the previous (0)')
                
                wgt0 = np.linspace(0.0, ni0-1, num=ni0)/ni0-1
                wgt1 = 1.0 - wgt0
                
                #python i0 = colr[iord][0]  # IDL i0=[0,iord]
                weights[iord+1][i0] = weights[iord+1][i0] * s0[i0]/(s0[i0]*wgt0+tempS0*wgt1)
                
                if beg0 > 0:
                    
                    weights[iord+1][0:beg0-1]=0.0
                    s0[i0]=s0[i0]*wgt0+tempS0*wgt1
                    b0[i0]=b0[i0]*wgt0+tempB0*wgt1
                    
                if has_sig==1:
                    sig0[i0] = np.sqrt(sig0[i0]*sig0[i0]*wgt0+tempU0*tempU0*wgt1)
                    wgt1=  np.linspace(0.0, ni1-1, num=ni1)/ni1-1
                    wgt0 = 1.0 - wgt1
                    weights[iord][i1]=weights[iord][i1]*s1[i1]/(s1[i1]*wgt0+tempS1*wgt1)
                    
                if end1 < npix-0:
                    weights[iord][end1+1:npix-1]=0.0
                    s1[i1]=s1[i1]*wgt0+tempS1*wgt1
                    b1[i1]=b1*wgt0+tempB1*wgt1
                    
                if has_sig==1:
                    sig1[i1]=np.sqrt(sig1[i1]*sig1[i1]*wgt0+tempU1+tempU1*wgt1)
                    #plt.plot(w1,s1)
                    
                #IDL do some plot here ...
                    
                #plot the whole thing
                #ax.plot(w1,s1)  
                
                
                ax1.plot(w1,s1)
                ax1.plot(w1,b1)
                
                #ax1.plot(w1,weights[iord]*s1)


                ax1.plot(w0,s0)
                ax1.plot(w0,b0)
                
                #ax1.plot(w0,weights[iord+1]*s0)
            #check the location of this else
            else:
            
                print('in else')
                wgt1 = np.linspace(0.0, ni0-1, num=ni0)/ni0-1
                wgt0 = 1.0 - wgt1
                
                weights[iord+1][i0] = weights[iord+1][i0] * s0[i0]/(s0[i0]*wgt0+tempS0*wgt1)

                
                if end0 < npix-0:
                     
                     weights[iord+1][end0+1:npix-1]=0.0
                     s0[i0]=s0[i0]*wgt0+tempS0*wgt1
                     b0[i0]=b0[i0]*wgt0+tempB0*wgt1
                     
                if has_sig==1:
                     sig[i0] = np.sqrt(sig0[i0]*sig0[i0]*wgt0+tempU0*tempU0*wgt1)
                     wgt1=  np.linspace(0.0, ni1-1, num=ni1)/ni1-1
                     wgt0 = 1.0 - wgt1
                     weights[iord][i1]=weights[iord][i1]*s1[i1]/(s1[i1]*wgt0+tempS1*wgt1)
                     
                if beg1 > 0:
                     weights[iord][0:beg1-1]=0.0
                     s1[i1]=s1[i1]*wgt0+tempS1*wgt1
                     b1[i1]=b1*wgt0+tempB1*wgt1
                     
                if has_sig==1:
                     sig1[i1]=np.sqrt(sig1[i1]*sig1[i1]*wgt0+tempU1+tempU1*wgt1)
                     #plt.plot(w1,s1)
                     
            #Debug option
            if DEBUG==1:
                #some plots
                fig, ax = plt.subplots(figsize=(10, 10))
                ax.plot(ww[iord],sp[iord]/bb[iord])
                ax.plot(w0,s0/b0)
                ax.plot(w1,s1/b1)

                    
            else:
                #Calculate xmid
                xmid = (beg0 - end1) / 2.0

                #Perform element-wise division and find the top values along the columns
                scl0 = pyreduce.util.top(s0/b0, 1, poly=True)
                b0 = b0 * scl0
                
                scl0 = pyreduce.util.top(s0/b0, 1, poly=True) # why twice?

                scl1 = pyreduce.util.top(s1/b1, 1, poly=True)

                #Polynomial fitting (1st degree)

                poly0 = Polynomial.fit(w0, scl0, 1) # the poly itself
                scl0 = poly_fit(w0, scl0, degree=1) # change the meaning of scl0 now it has the coef of poly

                poly1 = Polynomial.fit(w1, scl1, 1) # the poly itself
                scl1 = poly_fit(w1, scl1, degree=1) # change the meaning of scl0 now it has the coef of poly


                #Generate a range of values similar to IDL's dindgen
                xx = np.linspace(0, 1, 101) * (np.min(w0) - np.max(w1)) + np.max(w1)

                if DEBUG == 1:
                    #checking the polynomium 
                    fig, ax = plt.subplots(figsize=(10, 10))
                    ax.plot(ww[iord],sp[iord]/bb[iord])
                    ax.plot(w1,s1/b1)
                    ax.plot(xx,poly1(xx))
                    ax.plot(w0,s0/b0)
                    ax.plot(xx,poly0(xx))
                
                n =  np.sum(scl0[1] * scl1[1] * xx * xx +
                   scl0[1] * scl1[0] * xx +
                   scl1[1] * scl0[0] * xx +
                   scl1[0] * scl0[0])
                
                d = np.sum((scl1[1] * xx + scl1[0]) ** 2)

                
                scl = n/d
                
                s1=s1*scl
                
                order_scales[iord]=scl
                tempS0=tempS0*scl
                
                #end else
            #making spectra
            
            
            # to adapt the IDL   sp[beg0:end0,iord+1]=s0>0.
            # I modified the s0 here
            s0[s0<0] = 0
            sp[iord+1][beg0:end0] = s0
            # I modified the b0 here
            b0[b0<1] = 1
            bb[iord+1][beg0:end0] = b0    
            
            #setting sig0
            if has_sig==1:
                unc[iord+1][beg0:end0] = sig0
             
                
            # I modified the s0 here
            s1[s1<0] = 0
            sp[iord][beg1:end1] = s1
            # I modified the b1 here
            b1[b1<1] = 1
            bb[iord][beg1:end1] = b1

            #setting sig0
            if has_sig==1:
                unc[iord][beg1:end1] = sig1





#=============================================================================

print('Starting forward loop')


beg1 = colr[iord0][0]
end1 = colr[iord0][1]

w1 = ww[iord0][beg1:end1]
s1 = sp[iord0][beg1:end1]
b1 = bb[iord0][beg1:end1]

if has_sig == 1:
    sig1=unc[iord0][beg1:end1]

if iord0 < iord0+1 :
    
  for iord in np.flip(np.linspace(iord0+1, nord-1, num=iord0)): #start with the biggest signal order and goes n-1 n-2 ...
      
      iord=int(iord)    
      beg0 = beg1  # Shift current order to previous
      end0 = end1
      w0 = w1
      s0 = s1
      b0 = b1
      
      if has_sig:
          sig0 = sig1
      
      beg1 = colr[iord][0]  # New current order
      end1 = colr[iord][1]
      
      w1 = ww[iord][beg1:end1]
      s1 = sp[iord][beg1:end1]
      b1 = bb[iord][beg1:end1]
      
      if has_sig:
          sig1 = unc[iord][beg1:end1]
          
          
      #Defining the overlap region in i0 and ni0 where i0 is an array of the index (pixels) and ni0 its length
      i0_pix = np.where((w0 >= np.min(w1)) & (w0 <= np.max(w1)))[0] # Overlap within the previous order 
      i0 = ((w0 >= np.min(w1)) & (w0 <= np.max(w1)))    #MASK
      
      ni0=len(i0_pix)                     # Overlap within the previous order ()
      ii0 = np.arange(len(i0))


      #WRANGE if ()
                    
      #Overlap with the current order based in the same logic as above
      i1_pix = np.where((w1 >= np.min(w0)) & (w1 <= np.max(w0)))[0] 
      i1 = ((w1 >= np.min(w0)) & (w1 <= np.max(w0)))    #MASK
        
      ni1=len(i1_pix)                    
      ii1 = np.arange(len(i1))
         
         
        #setting the WRANGE
     #   if WRANGE is not None :  # If exclusion regions are given, check the overlaps
      #      for iwrange in range(len(WRANGE)//2):
       #         iii1 = np.where((w1[i1[ii1]] < WRANGE[0, iwrange]) or (w1[i1[ii1]] > WRANGE[1, iwrange]))[0]
        #        if len(iii1) > 0:
         #           ii1 = ii1[iii1]
                 
       #         else:
        #            ni1 = 0      
        

#check content
      if (ni0 > 0) and (ni1 > 0):
           print('We have overlaping orders')
           
           #IDL original
           #tmpS0=bezier_interp(w1,s1,bezier_init(w1,s1,/DOUBLE),w0[i0],/DOUBLE)
           tempS0=pyreduce.util.bezier_interp(w1,s1,w0[i0])
           
           tempB0=pyreduce.util.bezier_interp(w1,b1,w0[i0])
           
           if has_sig==1:
               tempU0=pyreduce.util.bezier_interp(w1,sig1,w0[i0])

          

           tempS1=pyreduce.util.bezier_interp(w0,s0,w1[i1])
           
           tempB1=pyreduce.util.bezier_interp(w0,b0,w1[i1])
           
           if has_sig==1:
               tempU1=pyreduce.util.bezier_interp(w0,sig0,w1[i1])

       
           # add scaling 
           if SCALING is not None:
               
               scl0 = np.sum(abs(s0[i0[ii0]]/b0[i0[ii0]]))/np.sum(abs(tempS0[ii0]/tempB0[ii0]))
               scl1 = np.sum(abs(s1[i1[ii1]]/b1[i1[ii1]]))/np.sum(abs(tempS1[ii1]/tempB0[ii1]))
               
               scl=np.sqrt(scl0/scl1)
               s1=s1*scl
               if has_sig==1:
                   tempS0=tempS0*scl
                   tempU0=tempU0*scl
                   order_scales[iord]=scl
           #end scaling
            
           if (ww[iord][0]) > (ww[iord+1][0]):
               print('Current order (1) is redder than the previous (0)')
               
               wgt0 = np.linspace(0.0, ni0-1, num=ni0)/ni0-1
               wgt1 = 1.0 - wgt0
               
               #python i0 = colr[iord][0]  # IDL i0=[0,iord]
               weights[iord-1][i0] = weights[iord-1][i0] * s0[i0]/(s0[i0]*wgt0+tempS0*wgt1)
               
               if end0 < npix-1:
                   
                   weights[iord-1][end0+1:npix-1]=0.0
                   s0[i0]=s0[i0]*wgt0+tempS0*wgt1
                   b0[i0]=b0[i0]*wgt0+tempB0*wgt1
                   
               if has_sig==1:
                   sig0[i0] = np.sqrt(sig0[i0]*sig0[i0]*wgt0+tempU0*tempU0*wgt1)
                   wgt1=  np.linspace(0.0, ni1-1, num=ni1)/ni1-1
                   wgt0 = 1.0 - wgt1
                   weights[iord][i1]=weights[iord][i1]*s1[i1]/(s1[i1]*wgt0+tempS1*wgt1)
                   
               if beg1 > 0 :
                   weights[iord][0:beg1-1]=0.0
                   s1[i1]=s1[i1]*wgt0+tempS1*wgt1
                   b1[i1]=b1*wgt0+tempB1*wgt1
                   
               if has_sig==1:
                   sig1[i1]=np.sqrt(sig1[i1]*sig1[i1]*wgt0+tempU1+tempU1*wgt1)
                   #plt.plot(w1,s1)
                   
               #IDL do some plot here ...
                   
               #plot the whole thing
               #ax.plot(w1,s1)  
               
               
               ax1.plot(w1,s1)
               ax1.plot(w1,b1)
               
               #ax1.plot(w1,weights[iord]*s1)


               ax1.plot(w0,s0)
               ax1.plot(w0,b0)
               
               #ax1.plot(w0,weights[iord+1]*s0)
           #check the location of this else
           else:
           
               print('in else')
               wgt1 = np.linspace(0.0, ni0-1, num=ni0)/ni0-1
               wgt0 = 1.0 - wgt1
               
               weights[iord-1][i0] = weights[iord-1][i0] * s0[i0]/(s0[i0]*wgt0+tempS0*wgt1)

               
               if beg0 > 0:
                    
                    weights[iord-1][0:beg0-1]=0.0
                    s0[i0]=s0[i0]*wgt0+tempS0*wgt1
                    b0[i0]=b0[i0]*wgt0+tempB0*wgt1
                    
               if has_sig==1:
                    sig0[i0] = np.sqrt(sig0[i0]*sig0[i0]*wgt0+tempU0*tempU0*wgt1)
                    wgt1=  np.linspace(0.0, ni1-1, num=ni1)/ni1-1
                    wgt0 = 1.0 - wgt1
                    weights[iord][i1]=weights[iord][i1]*s1[i1]/(s1[i1]*wgt0+tempS1*wgt1)
                    
               if end1 < npix-1 :
                    weights[iord][end1+1:npix-1]=0.0
                    s1[i1]=s1[i1]*wgt0+tempS1*wgt1
                    b1[i1]=b1*wgt0+tempB1*wgt1
                    
               if has_sig==1:
                    sig1[i1]=np.sqrt(sig1[i1]*sig1[i1]*wgt0+tempU1+tempU1*wgt1)
                    #plt.plot(w1,s1)
                    
           #Debug option
           if DEBUG==1:
               #some plots
               fig, ax = plt.subplots(figsize=(10, 10))
               ax.plot(ww[iord],sp[iord]/bb[iord])
               ax.plot(w0,s0/b0)
               ax.plot(w1,s1/b1)

                   
           else:
               #Calculate xmid
               xmid = (beg0 - end1) / 2.0

               #Perform element-wise division and find the top values along the columns
               scl0 = pyreduce.util.top(s0/b0, 1, poly=True)
               
               scl1 = pyreduce.util.top(s1/b1, 1, poly=True)

               b1 = b1 * scl1
               
               scl1 = pyreduce.util.top(s1/b1, 1, poly=True)

               #Polynomial fitting (1st degree)

               poly0 = Polynomial.fit(w0, scl0, 1) # the poly itself
               scl0 = poly_fit(w0, scl0, degree=1) # change the meaning of scl0 now it has the coef of poly

               poly1 = Polynomial.fit(w1, scl1, 1) # the poly itself
               scl1 = poly_fit(w1, scl1, degree=1) # change the meaning of scl0 now it has the coef of poly


               #Generate a range of values similar to IDL's dindgen
               xx = np.linspace(0, 1, 101) * (np.min(w0) - np.max(w1)) + np.max(w1)

               if DEBUG == 1:
                   #checking the polynomium 
                   fig, ax = plt.subplots(figsize=(10, 10))
                   ax.plot(ww[iord],sp[iord]/bb[iord])
                   ax.plot(w1,s1/b1)
                   ax.plot(xx,poly1(xx))
                   ax.plot(w0,s0/b0)
                   ax.plot(xx,poly0(xx))
               
               n =  np.sum(scl0[1] * scl1[1] * xx * xx +
                  scl0[1] * scl1[0] * xx +
                  scl1[1] * scl0[0] * xx +
                  scl1[0] * scl0[0])
               
               d = np.sum((scl1[1] * xx + scl1[0]) ** 2)

               
               scl = n/d
               
               s1=s1*scl
               
               order_scales[iord]=scl
               tempS0=tempS0*scl
               
               #end else
           #making spectra
           
           
           # to adapt the IDL   sp[beg0:end0,iord+1]=s0>0.
           # I modified the s0 here
           s0[s0<0] = 0
           sp[iord-1][beg0:end0] = s0
           # I modified the b0 here
           b0[b0<1] = 1
           bb[iord-1][beg0:end0] = b0    
           
           #setting sig0
           if has_sig==1:
               unc[iord+1][beg0:end0] = sig0
            
               
           # I modified the s0 here
           s1[s1<0] = 0
           sp[iord][beg1:end1] = s1
           # I modified the b1 here
           b1[b1<1] = 1
           bb[iord][beg1:end1] = b1

           #setting sig0
           if has_sig==1:
               unc[iord][beg1:end1] = sig1


#=============================================================================

#Form the output arrays

# INDEX array will track assosciation of pixels with the original CCD pixels


for iord in np.flip(np.linspace(0, nord-1, num=nord)):    
    iord=int(iord)
    
    i1=colr[iord][0]
    i2=colr[iord][1]
    
    
    # if it is the first order, do...
    if (iord==0):
        wave=ww[iord][i1:i2]
        spec=sp[iord][i1:i2]
        blaz=bb[iord][i1:i2]
        
        if has_sig==1:
            
            sig = unc[iord][i1:i2]
            index = np.full(i2 - i1+1, iord)
        
    else:
        # add the orders together in one single spectrum
        wave = np.append(wave, ww[iord][i1:i2])
        spec = np.append(spec, sp[iord][i1:i2])
        blaz = np.append(blaz, bb[iord][i1:i2])

        if has_sig==1:
            
            sig   = np.append(sig, unc[iord][i1:i2])
            index = np.append(index, np.full(i2 - i1+1, iord))

isort = np.argsort(wave)

wave=wave[isort]
spec=spec[isort]
blaz=blaz[isort]


#does not work if has_sig == 0 
try:
    index=index[isort]
except:
    index = 0
    

if has_sig==1:
    sig=sig[isort]
    
  

spliceout = {'wave':wave, 'spec':spec, 'cont':blaz, 'unc':sig, 'index': index}    

# For testing propose save the output in a pickle

with open('./splice-out.p','wb') as sp_save:
    pkl.dump(spliceout, sp_save)
    
  
    

# END    
































































#