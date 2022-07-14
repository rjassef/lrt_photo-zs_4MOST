import numpy as np
from astropy.io import fits
import astropy.units as u
import warnings
import pkg_resources
import re

class ProcPhot(object):

    def __init__(self, noise_floor):

        #Process the fluxes to get them ready for the fitting codes.
        self.process_fluxes(noise_floor=noise_floor)

        #Finally, create the bandmag.dat file. 
        self.create_bandmagfile()

        return

    def process_fluxes(self, noise_floor):

        #Figure out the usable magnitudes.
        self.jyuse = np.ones(self.jy.shape)

        with warnings.catch_warnings():
            # Ignore the warning about an invalid value in less. 
            warnings.filterwarnings('ignore')

            #Bad magnitudes
            kbad = (np.isnan(self.jy)) | (self.jy<0)
            self.jyuse[kbad] = 0
            self.jy[kbad] = 0.

            #Upper magnitudes
            kup = (~kbad) & ( (np.isnan(self.ejy)) | (self.ejy<0))
            self.jyuse[kup] = 2
            self.ejy[kup] = self.jy[kup]
            self.jy[kup] = 0.

            #Set a noise floor of 10% in the magntide errors.
            self.ejy = np.where((self.jyuse==1) & (self.ejy<noise_floor*self.jy),noise_floor*self.jy,self.ejy)
        return

    def create_bandmagfile(self):

        #Read the catalog of bands. 
        stream = pkg_resources.resource_stream(__name__, "filters.txt")
        filter_name = dict()
        cal_type = dict()
        jyzero = dict()
        for line in stream:
            x = line.decode("utf-8").split()
            filter_name[x[0]] = x[1]
            cal_type[x[0]] = x[2]
            jyzero[x[0]] = x[3]
        stream.close()

        #Now write the bandmag.dat file.
        cato = open("bandmag.dat","w")
        for col in self.photcols:
            cato.write("{0:25s} {1:5s} {2:s}\n".format(filter_name[col], cal_type[col], jyzero[col]))
        cato.close()

        return
