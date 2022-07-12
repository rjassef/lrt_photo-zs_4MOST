#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import SED_Model

def plot(id):
    k = np.argmin(np.abs(id-ids))
    print(id, ids[k])
    print(z[k],k)
    gal = SED_Model.lrt_model()
    gal.jy = jy[k,:]
    gal.ejy = ejy[k,:]
    gal.jyuse = jyuse[k,:]
    gal.zphot = z[k]

    gal.kc_fit()
    gal.plot()
    plt.show()

#Read the photometric zero points.
jyzero = np.loadtxt("bandmag.dat",usecols=[2])
nchan = len(jyzero)

#Read the photometry and redshift.
cat = np.loadtxt("small_sample.txt", usecols=np.concatenate([[0],np.arange(9,nchan*2+9,1,dtype=np.int32)]))
zpfile = np.loadtxt("zphot_all.dat",usecols=[0,2])
z = zpfile[:,1]

#Calculate the magnitudes
ids  = cat[:,0]
mag  = cat[:,1::2]
emag = cat[:,2::2]
jy = np.zeros(mag.shape)
ejy = np.zeros(emag.shape)
jyuse = np.zeros(emag.shape)
for j in range(nchan):

    #Bad magnitudes
    kbad = (np.isnan(mag[:,j])) | (mag[:,j]<0) | (mag[:,j]>90.)

    #Upper magnitudes
    kup = (mag[:,j]>0) & (mag[:,j]<90) & (emag[:,j]>=90.0) & (~kbad)
    jyuse[kup,j] = 2
    ejy[kup,j] = jyzero[j]*10.**(-0.4*mag[kup,j])

    #The rest
    kgood = (mag[:,j]>0) & (~kup) & (~kbad)
    jyuse[kgood,j] = 1
    jy[kgood,j]    = jyzero[j]*10.**(-0.4*mag[kgood,j])
    ejy[kgood,j]   = 0.4*np.log(10.)*emag[kgood,j]*jy[kgood,j]

#Set a noise floor of 10% in the magntide errors.
ejy = np.where((jyuse==1) & (ejy<0.1*jy),0.1*jy,ejy)
