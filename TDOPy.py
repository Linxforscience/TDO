#!/usr/bin/env python
''' Affichage d'un periodogramme issu des donnees TDO
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as ml

VALNFFT = 8192
OVERLAP = 256
PADTO = 7168

class TDOpsd:
    def __init__(self,file):
        b = np.fromfile(file,dtype = np.dtype(np.int16))
        # print b.shape
        self.a = b[256:]
        
        # print self.a.shape
        b = []
    def prbuff(self,n):
        for i in xrange (n):
            print self.a[i],
    def pwd(self):
        plt.subplot(111)
        (Pxx, freqs, bins, im) = plt.specgram(self.a, NFFT=VALNFFT,Fs=50000000, noverlap=OVERLAP, scale_by_freq=False, pad_to=PADTO)
        np.savetxt('res',Pxx)

        s = "SPCGRM" + '-N-' + str(VALNFFT) + '-O-'+ str(OVERLAP) + '-P-' + str(PADTO)
        plt.title(s)
        plt.show()



if __name__ == '__main__':
    try :
        t = TDOpsd('input.dat')
        t.pwd()
        exit(0)
    except KeyboardInterrupt :
        sys.exit(3)
