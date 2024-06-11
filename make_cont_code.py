#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 11:39:23 2024

@author: heitor
"""

#% Compiled module: $MAIN$.
#% Compiled module: RDECH.
#% Compiled module: MRDFITS.
#% Compiled module: FXPOSIT.
#% Compiled module: MRD_HREAD.
#% Compiled module: FXPAR.
#% Compiled module: GETTOK.
#% Compiled module: VALID_NUM.
#% Compiled module: FXMOVE.
#% Compiled module: MRD_SKIP.
#% Compiled module: MATCH.
#% Compiled module: MRD_STRUCT.
#% Compiled module: MODECH.
#% Compiled module: TAG_EXIST.
#% Compiled module: SXDELPAR.
#% Compiled module: SXPAR.
#% Compiled module: SXADDPAR.
#% Compiled module: MKWAVE.
 # no RADVEL card in header - barycentric wavelengths
#% Compiled module: MAKE_CONT.
#% Compiled module: MIDDLE.
#% Compiled module: SPLICE_ECH.
#% Compiled module: BEZIER_INIT.
#% Compiled module: BEZIER_INTERP.
#% Compiled module: TOP.

from astropy.io import fits
import os
import numpy as np
from scipy.io import readsav
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import Polynomial


import bezier 

import pickle as pkl

import pyreduce

import warnings
warnings.filterwarnings('ignore')


def poly_fit(x, y, degree=1):
    p = Polynomial.fit(x, y, degree)
    return p.convert().coef


def bezier_init(x,y):
    
    #Computes automatic control points for cubic bezier splines
    
    N = len(x)
    isrt = np.argsort(x)
    
    
    # just a checking if X is increasing
    i  = np.where(isrt[1:N-1]-isrt[0:N-2] != 1 )
    ni = len(i)
    

    if ni > 0:
        print('Arrays X and Y in the call to BEZIER_INIT should be sorted so that X is increasing')
        
        
    i  = np.where(x[1:N-1] == x[0:N-2])
    ni = len(i)
    
    
    if ni > 0:
      print('Array X in the call to BEZIER_INIT should not have identical values')
      
    y2=x
    H2=x[1]-x[0]
    DER2=(y[1]-y[0])/H2
    y2[0]=DER2
    
    for I in np.linspace(1, N-2, num=N-2):
        I = int(I)
        H1=H2
        DER1=DER2
        H2=x[I+1]-x[I]
        DER2=(y[I+1]-y[I])/H2
        ALPHA= (1.0+H2/(H1+H2))/3.0
        
        if DER1*DER2 > 0.0:
            y2[I]=DER1*DER2/(ALPHA*DER2+(1.0-ALPHA)*DER1)
        else:
            y2[I]=0.0
            
    
    y2[N-1]=DER2
    
    return y2



# bugs!
def bezier_interp(xa,ya,y2a,x):
    x=np.array(x)
    N=len(xa)
    #M=len(x) #used in case of double
    
    y=x*0
    
    ii  = np.where((x > min(xa)))  and np.where((x < max(xa)))
    
    nii = len(ii)

    if nii == 0:
        return y
    
    indices = np.searchsorted(xa, x)
    
    KLO = np.clip(indices < (N - 2), 0, 1)
    print('KLO')
    print(KLO)
    KHI=KLO+1
    print('KHI')
    print(KHI)
    
    H=xa[KHI]-xa[KLO]
    y1=ya[KLO]
    y2=ya[KHI]
    
    
    A=(xa[KHI]-x)/H
    B=(x-xa[KLO])/H
    
    C0=y1+H/3.0*y2a[KLO]
    C1=y2-H/3.0*y2a[KHI]
    
    y=A*A*A*y1+3.0*A*A*B*C0+3.0*A*B*B*C1+B*B*B*y2
    
    return y 



#This mostly sanitizes the input by removing masked values and d
#uplicate entries Note that in case of duplicate entries (in x_old) 
#the results are not well defined as only one of the entries is used 
#and the other is discarded



#devoloper test, make cont not as function





#calling the make_cont


#*c-ech with wave solution 
spec='./data/goodech-Blue.fits'

#spec='./data/Reduced_red/spec.ech'


ech = fits.open(spec)
#running test just for red arm
#RED------------------
directory = 'Reduced_blue'

ech = fits.open(spec)

eee = fits.open(spec)

print('============================')

print('HEADER')
print('============================')
print(ech.info())
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


#-----------
#reference 


#function make_cont, ech, blaze, TEMPLATE=template, DEBUG=debug $
#                  , COLRANGE=col_range, YRANGE=yrange, PARAM=param $
#                  , WRANGE=wrange, POLY=ipoly, PLOT=plot, SCALING=scaling $
#                  , SPLICE_ONLY=splice_only,WEIGHTS=weights $
#                  , ORDER_SCALES=order_scales, FILENAME=filename



                    
#splice(ech[1], wave, spec, blzcoef, index=0, COLRANGE=list(col_range))



#-----------

#BUG fixing the the wrong number of col_rang ERROR
aa = np.zeros(shape=(45, 2))
for i in np.linspace(0, len(col_range)-1, num=len(col_range)):
    i=int(i)
    aa[i]=col_range[i]
    aa[i][1]=int(aa[i][1]+1)

    
col_range=aa.astype(int)


    
#-----------

#defining the calling parameter for make_cont
ech=ech[1]
wave=wave[0]
spec=spec[0]
blaze=blzcoef

TEMPLATE=None
DEBUG=None

COLRANGE=list(col_range)

YRANGE=None

#PARAM=None

#testing parameters 
PARAM=[10, 8.e5, 1.e7, 1.]


WRANGE=None
POLY=None
PLOT=None
SCALING=None
SPLICE_ONLY=None
WEIGHTS=None
ORDER_SCALES=None
FILENAME=None

#-----------


#begin of the make_cont "function as a code to make it easy for debuging 


#n_params() number of arguments provided in the fucntion calling 



#def make_cont(*args, **kwargs):
#    if len(args) < 2:
#        print('Usage: ech_new = make_cont(echfile, blaze,')
#        print('               [,TEMPLATE=template[,DEBUG=debug')
#        print('               [,PLOT[,YRANGE=yrange[,PARAM=param')
#        print('               [,COLRANGE=colrange[,WRANGE=wrange]]]]]]])')
#        print(' ech       - ech-structure containing the raw spectrum and the wavelength solution.')
#        print(' blaze     - a 2D array containing the blaze function.')
#        print(' template  - a [2,*] array with a template continuum normalized spectrum')
#        print('             with the wavelength in the first column.')
#        print(' inst_mode - a string that is used to select the default set of parameters.')
#        print(' yrange    - sets vertical plot range for the plot.')
#        print(' param     - an array of 4 parameters (see the code).')
#        print(' colrange  - an array in form [2,nord] giving starting, ending columns')
#        print(' wrange    - an array in form [[Wl1,Wr1],[Wl2,Wr2],...[WlK,WrK]] giving')
#        print('             left and right boundaries to exclude from the fit. Useful')
#        print('             when dealing with emission lines or strong cosmics.')
#        print(' filename  - optional label for plots')
#        print(' /POLY     - use polynomial instead of optimal filtering.')
#        print(' /SCALING  - scale spectral orders when splicing (recommended).')
#        print(' /PLOT     - displays interactive adjustments to the continuum fit.')
#        print(' /DEBUG    - allows examining every order and prevents making')
#        print('             modifications in the ECH file. This flag also sets the PLOT flag.')
#        print(' /SPLICE_ONLY - return spliced spectrum, wavelength scale, blaze functions and')
#        print('             uncertainties in a single structure.')
#        print('The output ech structure contains spliced spectrum subdivided to the spectral')
#        print('intervals of the additional orders and the corresponding continuum normalization.')
#        return 0

#  Usage
#if __name__ == "__main__":
    # Simulate the function being called with command line arguments
#    make_cont(*sys.argv[1:])



#check if the ech have all the columns I need. We need 8.
#not the best way to check. But let's keep it that way for now


#skipping all this checking

if len(ech.columns) == 8 :
     e = ech

else:
    print('MAKE_CONT first parameter must be a valid ech structure or a file')


#save the original
e_orig=e


#check all keyword // Args
if FILENAME != None:
    fname= FILENAME
    
if DEBUG != None:
    d=1
    #plot...
    
    
#checking if we have WAVE

#optimize this printing to show all columns what they are 
have_wave=False
for c in ech.columns:
    print('ech structute has, '+ c.name)
    
    if c.name == 'WAVE':
        have_wave=True
   
        
# highlighing the WAVE structure 
if have_wave is False:
    print('MAKE_CONT: the ech structure does not contain the wavelength scale!')
else:
    print('MAKE_CONT found a WAVE in ech structure')



if sig is not None:
    print('Applying Sigma...')
    # Sanity check for 'SIG'
    if not 'SIG' in ech.columns.names :
        has_sig = 1 
else:
    has_sig = 0
       

#numbers based in the first order
npix = ech.data.dtype['SPEC'].shape[1]  # Order length in pixels
nord = ech.data.dtype['SPEC'].shape[0]  # Number of spectral orders



#defining colrange if it was not defined in one of the input parameters, COLRANGE
if COLRANGE != None:
    
    colrange=COLRANGE

else:
    
    colrange = np.zeros((2, nord), dtype=int)
    colrange[0, :] = 0
    colrange[1, :] = npix - 1 # CHECK if it need -1 or -0



#How exactly will be defined the Template ( SKIP FOR NOW )

if TEMPLATE!= None:
    
    template_wave = TEMPLATE[1] # Template wavelength scale
    template_sp   = TEMPLATE[0] # Template spectrum
    #...
    
    
    
#target = e.data['HEAD']['OBJECT'] # fix that becuase it does not work for the fits format I have here
target = 'aa'

w=e.data['WAVE'][0]


#check spec and blaze #DOIT do it better
if npix != len(blaze[0]):
    print('Inconsistent row length in blaze and spectrum arrays')
    
    
if nord != len(blaze):
    print('MAKE_CONT: inconsistent column length in blaze and spectrum arrays')




# min and max of the beging and end of the first and last order
#w[0][colrange[0][0]:colrange[0][1]]
#w[nord-1][colrange[nord-1][0]:colrange[nord-1][1]]


min_w0 = min(w[0][colrange[0][0]:colrange[0][1]]) 
min_wn = min(w[nord-1][colrange[nord-1][0]:colrange[nord-1][1]])

wmin = min([min_w0,min_wn])

max_w0 = max(w[0][colrange[0][0]:colrange[0][1]]) 
max_wn = max(w[nord-1][colrange[nord-1][0]:colrange[nord-1][1]])

wmax = max([min_w0,min_wn])



#set up parameters
if len(PARAM) >= 4:
    print('Reading the parameters provided')
    par_0 = PARAM[0] #Number of iterations
    par_1 = PARAM[1] #Smoothness for initial guess
    par_2 = 1.e-4
    par_3 = PARAM[2] #Smoothness  for the polishing touch
    par_4 = 0.01*(1.e-2<1.e0/np.sqrt(np.median(e.data['SPEC'])))
    par_5 = PARAM[3] #Vertical scaling

#standard values
else:
    print('Adopting standard parameters')
    par_0 = 10
    par_1 = 5.e5
    par_2 = 1.e-4
    par_3 = 5.e6
    par_4 = 0.01*(1.e-2<1.e0/np.sqrt(np.median(e.data['SPEC'])))
    par_5 = 1.e0
    

#avoid 0 blaze function



b = np.where(blaze<1, 1, blaze)


for iord in np.linspace(0, nord-1, num=nord):
    iord=int(iord)
    i0 = colrange[iord][0]
    i1 = colrange[iord][1]
    b[iord][i0:i1]= pyreduce.util.middle(b[iord][i0:i1],1) # whats is middle ?
    
    
#% Compiled module: MIDDLE.
#syntax: middle,f,{filter/order}[,X=xarg,[ITER=iter[,EPS=eps $
#             [,MIN=mn[,[MAX=mx[,/POLY[,LAM2=lam2[,WEIGHT=wgt]]]]]]]]
#where f      is the function to fit.
#      filter is the smoothing parameter for the optimal filter.
#             If POLY is set, it is interpreted as the order
#             of the smoothing polynomial.
#      xarg   optional X argument of data points. Must have the same size as f.
#      iter   is the maximum number of iterations [def: 40].
#      eps    is convergence level [def: 0.001].
#      mn     minimum function values to be considered [def: min(f)].
#      mx     maximum function values to be considered [def: max(f)].
#      lam2   constraint on 2nd derivative.
#      wgt    vector of weights.
#       1




#what is the sB? sB was not created yet created in line 166 after the splice

#YRANGE read or setup
if YRANGE != None:
    yr=YRANGE
else:
    #yr=[0,2*median(sB)] #IDL #FIXIT
    yr = [0, 2*np.median(spec[-1])]
    
    
    
#Create an equispaced wavelength grid covering the whole thing

#Calculate dwave
dwave = np.abs(w[nord // 2][npix // 2] - w[nord // 2][npix // 2 - 1]) * 0.5
# Calculate nwave
nwave = int(np.ceil((wmax - wmin) / dwave) + 1)
# Generate wave array
wave = (wmax - wmin) * np.arange(nwave + 1) / nwave + wmin
   
 
    
for iord in np.linspace(1, nord-1, num=nord): #Mark gaps between non-overlapping orders
    iord=int(iord)    
    
    miss = np.where( ( wave >= w[iord-1][colrange[iord-1][1]-1] ) &
                ( wave <= w[iord][colrange[iord][0]] ))[0]
   
    nmiss = len(miss)

    if nmiss > 0:
       minmax_wave_miss = [np.min(wave[miss]), np.max(wave[miss])] #check if it works
       #defining min and max for the gap
       
       if WRANGE != None:
           WRANGE.append(minmax_wave_miss)  
       else:
          WRANGE = [minmax_wave_miss]
    
    
e.data['CONT']=b    
    

#=========================================================================

#SPLICE

#=========================================================================


#call splice function 

#ADD SPLICE FUNCTION HERE

#what it must return?
    
# Splice will provide the wsort and 

with open('./splice-out.p', 'rb') as f:
    spliceout = pkl.load(f)
    
wsort = spliceout['wave']
ss    = spliceout['spec']
bb    = spliceout['cont']
sig   = spliceout['unc']

index = spliceout['index']

if SPLICE_ONLY != None:
    out={'wave':wsort, 'spec':ss, 'cont':bb, 'unc':sig}
    
    
j = np.unique(wsort, return_index=True)

j_uniq  = j[0]
j_index = j[1]

#wsort_u = wsort[j] original IDL
# picking just the unique ones abd adding to wsort_u also making it a np.array again (not sure if it is needed)
wsort_u = np.array([wsort[x] for x in j_index])


    
 
if TEMPLATE != None:
    temp=1
    #what is bezier_init?
    temp=pyreduce.util.bezier_interp(template_wave-wmin, template_sp,wsort_u-wmin)

else:
    
    temp=1    
    

sB = (np.array([ss[x] for x in j_index]) / np.array([bb[x] for x in j_index])) / temp

   
#cosmic removal function. check this function 
#sB = clean_cosmic(wsort_u-wmin,sB,PLOT=plot,NHIST=n_histogram,SMOOTH1=smooth1,SMOOTH2=smooth2)   
    
weight=pyreduce.util.middle(sB,0.5,x=wsort_u-wmin)    

nwgt=len(weight)
    
#check if it works, I changed index here
we_d = np.concatenate( ([0.], 2. * weight[1:nwgt-1] - weight[0:nwgt-2] - weight[2:nwgt-0], [0.]) )

weight=weight/pyreduce.util.middle(weight,3.*par_1) + we_d
   
 
#weight=weight>0 works
weight[weight<0] = 1


#check if it works
weight=pyreduce.util.bezier_interp(wsort_u-wmin,weight,wave-wmin)


weight = weight/max(weight)

# Incorporate user-rejected regions


if WRANGE != None:    # If exclusion regions are given, set weights to 0
    for i in np.linspace(0, len(WRANGE)//2-1, num=len(WRANGE)//2-1):
      condition = (wave >= WRANGE[0, i]) & (wave <= WRANGE[1, i])
      weight[condition] = 0


spec_B2 = bezier_init(wsort_u-wmin,sB)

ss_B = pyreduce.util.bezier_interp(wsort_u-wmin,sB,wave-wmin)


bbb=pyreduce.util.middle(bb,1.)

niter=0
cont_B=1.0

#iter:

#=========================================================================

#All the same except for the variable names


# if exclusion regions are given, set weights to zero
if WRANGE != None:
    
    for i in np.linspace(0, len(WRANGE)//2-1, num=len(WRANGE)//2-1):
        condition = (wave >= WRANGE[i][0]) & (wave <= WRANGE[i][1])
        ss_B[condition] = 0.0



c=ss_B/cont_B



if POLY!=None:
    if TEMPLATE != None:
        for ii in range(par_0 + 1):
            c = pyreduce.util.middle(c, POLY if POLY > 2 else 15, eps=par_2, weight=weight, POLY=True)
        c = pyreduce.util.middle(c, POLY if POLY > 2 else 15, eps=par_4, weight=weight, POLY=True) * cont_B
    else:
        for ii in range(par_0 + 1):
            c = pyreduce.util.top(c, POLY if POLY > 2 else 15, eps=par_2, weight=weight, POLY=True)
        c = pyreduce.util.top(c, POLY if POLY > 2 else 15, eps=par_4, weight=weight, POLY=True) * cont_B
else:
    if TEMPLATE != None:
        for ii in range(par_0 + 1):
            c = pyreduce.util.middle(c, par_1, eps=par_2, weight=weight, lambda2=par_3)
        c =pyreduce.util. middle(c, par_1, eps=par_4, weight=weight, lambda2=par_3) * cont_B
    else:
        for ii in range(par_0 + 1):
            c = pyreduce.util.top(c, par_1, eps=par_2, weight=weight, lambda2=par_3)
        c = pyreduce.util.top(c, par_1, eps=par_4, weight=weight, lambda2=par_3) * cont_B

print('if Poly defined value of c:', c)   



cont_B = c*par_5
cont_B = pyreduce.util.middle(cont_B,1.0)


#---
# PLOT
#---


#Interpolate the continuum back on every sp. order


#---
# PLOT
#---


#Define the output 

# Perform the operations
#cont_B2 = bezier_init(wave - wmin, cont_B)
cont = pyreduce.util.bezier_interp(wave - wmin, cont_B, wsort - wmin)



for iord in range(nord):
    
    i1 = colrange[iord][0]
    i2 = colrange[iord][1]
    
    i = iord #bug with index
    
    # format issue need [0]
    e.data['WAVE'][0][iord][i1:i2] = wsort[i]
    e.data['SPEC'][0][iord][i1:i2] = ss[i]
    
    has_sig=0 # fix bug, does not indentify has_sig
    if has_sig == 1 :
        e.data['SIG'][0][iord][i1:i2] = sig[i]
        
        
    e.data['CONT'][0][iord][:] = 1.0
    e.data['CONT'][0][iord][i1:i2] = cont[i] * bbb[i]



# Print results (for verification)
print('e.wave:', e.data['WAVE'])
print('e.spec:', e.data['SPEC'])
print('e.sig:', e.data['SIG'])
print('e.cont:', e.data['CONT'])


#Inspect every sp. order


#---
# PLOT
#---























































# END    































































#