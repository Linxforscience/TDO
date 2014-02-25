# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:50:17 2013

@author: X. FABREGES


"""
from scipy import int16, mean, arange, zeros, real, fft, argmax, pi, savetxt, fromfile, r_, ones, convolve, cos, sin, floor
from scipy.fftpack import fftfreq
import sys
import TDO_C as m

import cProfile

# Declaration globale des slices (faster)
slice1 = [slice(None)]
slice2 = [slice(None)]
slice1[-1] = slice(1, None)
slice2[-1] = slice(None, -1)

# Smoothing of the results
def smooth(x,window_len=11):
    ''' Remove fast oscillations from results '''
    s=r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    
    w=ones(window_len,'d')
    y=convolve(w/w.sum(),s,mode='valid')
    return y
  
# Fast trapezoidal rule. No error catching, no constants
def local_trapz(y):
    ''' Trapezoidal integration optimized for speed '''    
    global slice1
    global slice2
    
    return ((y[slice1]+y[slice2])).sum(axis=-1)

def tdo_fft(inputfile, outputfile):
    ''' Perform Fourier transform and return frequency evolution '''
    global approx_sine
    global approx_cosine
    
    # Input parameters
    fourier_le = 4096   # Fourier length
    time_le    = 2048   # timewindow
    dfmin      = 1      # Frequency resolution
    dt         = 2e-8   # timestep of acquisition
    
    # Lecture du fichier
    fid = open(inputfile, 'rb')
    fid.seek(512)   # Skip useless header
    V = fromfile(fid, int16, -1, '')  
    fid.close()
    V = V-mean(V)  # Remove zero frequency contribution to FT
    t = arange(0, fourier_le)*dt
    
    pstart     = 1      # First timewindow
    pend       = int(floor((len(V)-fourier_le)/time_le)) # Last timewindow
    
    # Approximation of main frequency
    Vf         = abs(real(fft(V[0:time_le])))
    tf         = fftfreq(time_le, dt)
    fmax       = zeros((pend+1-pstart, 3))
    fmax[0, 1] = tf[argmax(Vf[0:int(time_le/2)])]
    fmax[0, 0] = time_le*dt/2
    
    # Calculation of constants
    expon   = -2j*pi*t
    deltaf0 =  tf[1]/1000
    
    # Precise determination of oscillating frequency via DFT    
    for i in xrange(pstart, pend):
        # Use last frequency as starting point
        a = fmax[i-1, 1]
        V_temp = V[i*time_le:i*time_le+fourier_le]
        
        # Previous frequency spectral weight    
        deltaf = deltaf0
        essaimax = abs(local_trapz(m.foo(a,V_temp,expon)))
        
        # Calculation of local derivative of Fourier transform
        # If derivative positive, then search for frequency in growing direction
        if abs(local_trapz(m.foo(a+deltaf,V_temp,expon))) > essaimax:
            while abs(deltaf)>dfmin:
                F = abs(local_trapz(m.foo(a+deltaf,V_temp,expon)))
                if F > essaimax:
                    essaimax = F
                    a += deltaf
                else:
                    deltaf = -deltaf/5
                    
            # Store frequency
            fmax[i, 0:2] = [(i*time_le+fourier_le/2)*dt, a-2.5*deltaf]
        # Lower frequency otherwise
        else:        
            while abs(deltaf)>dfmin:
                F = abs(local_trapz(m.foo(a-deltaf,V_temp,expon)))
                if F > essaimax:
                    essaimax = F
                    a -= deltaf
                else:
                    deltaf = -deltaf/5
        
            # Store frequency
            fmax[i, 0:2] = [(i*time_le+fourier_le/2)*dt, a+2.5*deltaf]
        
    # Save calculation in file
    fmax[:,2] = smooth(fmax[:,1])[5:-5]
    savetxt(outputfile, fmax)

i=0
while i<len(sys.argv):
    # read input and output files from arg
    if sys.argv[i]=='-i':
        inputfile=sys.argv[i+1]
    elif sys.argv[i]=='-o':
        outputfile=sys.argv[i+1]
    i += 1

cProfile.runctx('tdo_fft(inputfile,outputfile)', globals(), locals(), 'oupsss.txt')
#tdo_fft(inputfile, outputfile)