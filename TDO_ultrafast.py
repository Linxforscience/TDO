# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:50:17 2013

@author: X. FABREGES


"""
from scipy import int16, mean, arange, zeros, real, fft, argmax, pi, savetxt, fromfile, r_, ones, convolve, exp, floor, cos, sin
from scipy.fftpack import fftfreq
import sys
from numba import autojit
#import cProfile

# Declaration globale des slices (faster)
slice1 = [slice(None)]
slice2 = [slice(None)]
slice1[-1] = slice(1, None)
slice2[-1] = slice(None, -1)


def smooth(x,window_len=11):
    s=r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    
    w=ones(window_len,'d')
    y=convolve(w/w.sum(),s,mode='valid')
    return y
  
# Methode des trapezes modifiee : constante supprimee
def local_trapz(y):
    ''' Trapezoidal integration optimized for speed '''    
    global slice1
    global slice2
    
    return ((y[slice1]+y[slice2])).sum(axis=-1)

def comp_exp(arr):
    return (cos(-1j*arr)+1j*sin(-1j*arr))
    
def tdo_fft(inputfile, outputfile):
    ''' Perform Fourier transform and return frequency evolution '''
    comp_expf=autojit(comp_exp)    
    
    # Input parameters
    fourier_le = 1024   # Fourier length
    time_le    = 1024   # timewindow
    dfmin      = 100      # Frequency resolution
    dt         = 2e-8   # timestep of acquisition
    
    # Lecture du fichier
    fid = open(inputfile, 'rb')
    fid.seek(512)   # Skip useless header
    V = fromfile(fid, int16, -1, '')  
    fid.close()
    
    pstart     = 1      # First timewindow
    pend       = int(floor((len(V)-fourier_le)/time_le)) # Last timewindow
    
    V = V-mean(V)  # Remove zero frequency contribution to FT
    t = arange(0, fourier_le)*dt
    
    # Approximation of main frequency
    Vf         = abs(real(fft(V[0:fourier_le])))
    tf         = fftfreq(fourier_le, dt)
    fmax       = zeros((pend+1-pstart, 3))
    fmax[0, 1] = tf[argmax(Vf[0:int(fourier_le/2)])]
    fmax[0, 0] = time_le*dt/2
    
    # Calculation of constants
    expon   = -2j*pi*t
    deltaf0 =  tf[1]/1000
    
    # Windowing (uncomment desired window function)
    # Rectangular    
    window  =  ones(fourier_le)    
    # Cosinus    
    #window  =  arange(0,fourier_le)
    #window  =  1-cos(window*2*pi/(fourier_le-1))
    # Kaiser-Bessel    
    #window  =  arange(0,fourier_le)
    #window  =  0.402-0.498*cos(window*2*pi/(fourier_le-1))+0.098*cos(window*4*pi/(fourier_le-1))+0.001*cos(window*6*pi/(fourier_le-1))    
    
    
    # Precise determination of oscillating frequency via DFT    
    for i in xrange(pstart, pend):
        # Utilisation de la dernière valeur comme point de depart
        a = fmax[i-1, 1]
        V_temp = window*V[i*time_le:i*time_le+fourier_le]
        
        # Previous frequency spectral weight
        # Complex exponential time consuming ! Need a smarter way to perform this calculations
        deltaf = deltaf0
        essaimax=abs(local_trapz(V_temp*comp_expf(expon*a)))
        
        # Calculation of local derivative of Fourier transform
        # If derivative positive, then search for frequency in growing direction
        if abs(local_trapz(V_temp*comp_expf(expon*(a+deltaf)))) > essaimax:
            while abs(deltaf)>dfmin:
                F = abs(local_trapz(V_temp*comp_expf(expon*(a+deltaf))))
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
                F = abs(local_trapz(V_temp*comp_expf(expon*(a-deltaf))))
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
    # Scans to extract
    if sys.argv[i]=='-i':
        inputfile=sys.argv[i+1]
    elif sys.argv[i]=='-o':
        outputfile=sys.argv[i+1]
    i += 1

#cProfile.runctx('tdo_fft(inputfile,outputfile)', globals(), locals(), 'oupsss.txt')
tdo_fft(inputfile, outputfile)
#testf=autojit(tdo_fft)
#testf(inputfile, outputfile)