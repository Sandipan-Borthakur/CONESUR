#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 09:51:31 2024

@author: heitor
"""

import numpy as np


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



def bezier_interp(xa,ya,y2a,x):
    
    N=len(xa)
    #M=len(x) #used in case of double
    
    y=x*0
    
    ii  = np.where(x >= min(xa) and x <= max(xa))
    nii = len(ii)

    if nii == 0:
        return y
    
    indices = np.searchsorted(xa, x)
    
    KLO = np.clip(indices < (N - 2), 0, 1)
    KHI=KLO+1
    
    H=xa[KHI]-xa[KLO]
    y1=ya[KLO]
    y2=ya[KHI]
    
    
    A=(xa[KHI]-x)/H
    B=(x-xa[KLO])/H
    
    C0=y1+H/3.0*y2a[KLO]
    C1=y2-H/3.0*y2a[KHI]
    
    y=A*A*A*y1+3.0*A*A*B*C0+3.0*A*B*B*C1+B*B*B*y2
    
    return y 
    
    
















































#