#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:50:17 2013

@author: X. FABREGES

"""
from scipy import int16, arange, zeros, real, fft, argmax
from scipy import pi, savetxt, fromfile, floor
from scipy.fftpack import fftfreq
import sys
import pp

#import cProfile

# Methode des trapezes modifiee : constante supprimee
def local_trapz(yarr):
    ''' Trapezoidal integration optimized for speed '''    
    slice1     = [slice(None)]
    slice2     = [slice(None)]
    slice1[-1] = slice(1, None)
    slice2[-1] = slice(None, -1)
    
    return ((yarr[slice1]+yarr[slice2])).sum(axis=-1)

def tdo_fft(inputfile, outputfile):
    ''' Perform Fourier transform and return frequency evolution '''
    # Input parameters
    fourier_le     = 1024   # Fourier length
    time_le        = 1024   # timewindow
    dfmin          = 0.01     # Frequency resolution
    dt             = 2e-8   # timestep of acquisition
    load_balancing = 1
    
    # Lecture du fichier
    fid = open(inputfile, 'rb')
    fid.seek(512)   # Skip useless header
    V = fromfile(fid, int16, -1, '')  
    fid.close()
    
    pstart     = 1      # First timewindow
    pend       = int(floor((len(V)-4*fourier_le)/time_le))+1 # Last timewindow
    
    t          = arange(0, fourier_le)*dt
    
    # Approximation of main frequency
    Vf         = abs(real(fft(V[0:fourier_le])))
    tf         = fftfreq(fourier_le, dt)
    fmax       = zeros((pend+2-pstart, 2))
    fmax[0, 1] = tf[argmax(Vf[0:int(fourier_le/2)])]
    fmax[0, 0] = time_le*dt/2
    
    # Calculation of constants
    expon   = -2j*pi*t
    deltaf0 =  tf[1]/1000
    if deltaf0 < dfmin:
        deltaf0 = 10*dfmin
    
    # Start jobs
    job_server = pp.Server()
    ncpus = int(job_server.get_ncpus())
    
    serv_jobs = []

    # Load-balancing
    # Last processes are faster
    # If nprocess = ncpus, half of the core remains mostly idle    
    if load_balancing == 1:
        nprocess = ncpus*4
    else:
        nprocess = ncpus
        
    pstart_b = pstart        
    for i in range(0, nprocess):
        if nprocess == 1:
            pend_b = pend
        else:
            pend_b = int(pstart_b+floor(pend/nprocess))
            
        print(pstart_b, pend_b)
        
        args_tuple = (pstart_b, pend_b+1, \
        V[pstart_b*time_le:(pend_b+1)*time_le+fourier_le], \
        dt, dfmin, deltaf0, expon, fourier_le, time_le,)
        
        serv_jobs.append(job_server.submit(find_freq, args_tuple, \
        (local_trapz,)))
        
        pstart_b = pend_b

    pstart_b = pstart
    for i in range(0, nprocess):
        if nprocess == 1:
            pend_b = pend
        else:
            pend_b = int(pstart_b+floor(pend/nprocess))
        
        fmax[pstart_b:pend_b, :] = serv_jobs[i]()
        pstart_b = pend_b
        
    # Save calculation in file
    savetxt(outputfile, fmax)
    job_server.print_stats()

def find_freq(ps, pe, V, dt, dfmin, deltaf0, expon, fourier_le, time_le):
    '''Perform DFT of signal and return main frequency'''
    from scipy import zeros, real, fft, argmax, exp, arange , cos, pi, mean
    from scipy.fftpack import fftfreq

    Vf         = abs(real(fft(V[0:fourier_le]-mean(V[0:fourier_le]))))
    tf         = fftfreq(fourier_le, dt)
    fmax       = zeros((pe-ps, 2))
    fmax[0, 1] = tf[argmax(Vf[0:int(fourier_le/2)])]
    fmax[0, 0] = (ps*time_le+fourier_le/2)*dt
    
    # Rectangular
    #window  =  ones(fourier_le)    
    # Cosinus    
    window =  arange(0, fourier_le)
    window =  1-cos(window*2*pi/(fourier_le-1))
    
    for i in xrange(1, pe-ps):
        # Utilisation de la dernière valeur comme point de depart
        a      = fmax[i-1, 1]
        V_temp = window*V[i*time_le:i*time_le+fourier_le]
        
        # Previous frequency spectral weight
        # Complex exponential time consuming
        # Need a smarter way to perform this calculations
        deltaf   = deltaf0
        essaimax = abs(local_trapz(V_temp*exp(expon*a)))
        
        # Calculation of local derivative of Fourier transform
        # If derivative positive, then search for frequency in growing direction
        if abs(local_trapz(V_temp*exp(expon*(a+deltaf)))) > essaimax:
            while abs(deltaf)>dfmin:
                F = abs(local_trapz(V_temp*exp(expon*(a+deltaf))))
                if F > essaimax:
                    essaimax = F
                    a += deltaf
                else:
                    deltaf = -deltaf/5
                    if (abs(deltaf) < dfmin) and (abs(deltaf) > dfmin*4.9):
                        deltaf=deltaf/abs(deltaf)*1.01*dfmin
                    
            # Store frequency
            fmax[i, 0:2] = [((i+ps)*time_le+fourier_le/2)*dt, a-2.5*deltaf]
        # Lower frequency otherwise
        else:        
            while abs(deltaf)>dfmin:
                F = abs(local_trapz(V_temp*exp(expon*(a-deltaf))))
                if F > essaimax:
                    essaimax = F
                    a -= deltaf
                else:
                    deltaf = -deltaf/5
                    if (abs(deltaf) < dfmin) and (abs(deltaf) > dfmin*4.9):
                        deltaf=deltaf/abs(deltaf)*1.01*dfmin
        
            # Store frequency
            fmax[i, 0:2] = [((i+ps)*time_le+fourier_le/2)*dt, a+2.5*deltaf]
    
    return fmax[1:, :]

i = 0
while i < len(sys.argv):
    # Scans to extract
    if sys.argv[i] == '-i':
        inputf = sys.argv[i+1]
    elif sys.argv[i] == '-o':
        outputf = sys.argv[i+1]
    i += 1

#cProfile.runctx('tdo_fft(inputfile,outputfile)', globals(),\
#locals(), 'oupsss.txt')
tdo_fft(inputf, outputf)