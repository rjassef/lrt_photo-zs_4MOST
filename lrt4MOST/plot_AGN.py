#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import SED_Model
import re

def plot(id, newid, fname=None):
    #k = np.argmin(np.abs(id-ids))
    k = np.argmin(np.abs(newid-newids))
    gal = SED_Model.lrt_model()
    gal.jy = jy[k,:]
    gal.ejy = ejy[k,:]
    gal.jyuse = jyuse[k,:]

    gal.ejy[np.isnan(gal.ejy)] = 0.

    #Zphot fit.
    gal.zphot = z[k]
    gal.kc_fit()
    if fname is None:
        gal.plot()
    else:
        fname_use = re.sub(".png",".zphot.png",fname)
        gal.plot_to_file(fname_use)

#Read the photometric zero points.
jyzero = np.loadtxt("bandmag.dat",usecols=[2])
nchan = len(jyzero)

#Read the photometry and redshift.
data = np.loadtxt("sample.txt")
zpfile = np.loadtxt("zphot_all.dat",usecols=[0,2])
z = zpfile[:,1]

#Set the new id based on the line position. 
newids = np.arange(1,data.shape[0]+1)

#Calculate the magnitudes
ids  = data[:,0]
jy   = data[:,2::2]
ejy  = data[:,3::2]

jyuse = np.ones(jy.shape)

#Bad magnitudes
kbad = (np.isnan(jy)) | (jy<0)
jyuse[kbad] = 0
jy[kbad] = 0.

#Upper magnitudes
kup = (~kbad) & ( (np.isnan(ejy)) | (ejy<0))
jyuse[kup] = 2
ejy[kup] = jy[kup]
jy[kup] = 0.

#Set a noise floor of 10% in the magntide errors.
ejy = np.where((jyuse==1) & (ejy<0.1*jy),0.1*jy,ejy)

#Load the AGN ids.
agn_ids = np.loadtxt("AGN.dat",usecols=[0])
agn_newids = np.loadtxt("AGN.dat", usecols=[-1])

#Get the redshift distribution. 
z_agn = np.zeros(len(agn_ids))
for i,agn_newid in enumerate(agn_newids):
    k = np.argmin(np.abs(agn_newid-newids))
    z_agn[i] = z[k]
plt.hist(z_agn)
plt.savefig("z_agn_hist.png")

#Draw 100 random objects.
np.random.seed(158420)
ran_agn_ids = np.random.choice(agn_ids,100)
for agn_id in ran_agn_ids:
    #plot(agn_id.astype(np.int32))
    #plot(agn_id)
    k = np.argmin(np.abs(agn_id-agn_ids))
    plot(agn_id, agn_newids[k],
       fname="SED_plots/Obj{0:.0f}.png".format(agn_id))
