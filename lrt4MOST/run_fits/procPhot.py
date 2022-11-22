import numpy as np
import astropy.units as u
import warnings
import pkg_resources
from astropy.table import Table
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery

class ProcPhot(object):

    def __init__(self, noise_floor, fg_dust_corr):

        #Apply the foreground dust reddening corrections if requested.
        if fg_dust_corr:
            self.apply_fg_dust_corr()

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

    def apply_fg_dust_corr(self):

        #Query the SFD catalog for the position of all sources.
        coords = SkyCoord(self.ra, self.dec, frame='icrs')
        sfd = SFDQuery()
        ebv = sfd(coords)

        #Should find a classier way to do this for the lrt github, but this should be good enough for now.
        filename = "fg_red/S11_Table6.dat"
        stream = pkg_resources.resource_stream(__name__, filename)
        s11 = Table.read(stream, format='ascii')
        stream.close()

        bandname = "fg_red/S11_band_name_dictionary.txt"
        stream2 = pkg_resources.resource_stream(__name__, bandname)
        s11_band_name = dict()
        for line in stream2:
            line = line.decode("utf-8")
            x = line.split()
            s11_band_name[x[0]] = x[1]
        stream2.close()

        for k, col in enumerate(self.photcols):
            if col not in s11_band_name:
                print("ERROR: band {} not in dictionary for FG dust corrections.".format(col))
                continue
            R_band = s11['R_V_3.1'][s11['Band']==s11_band_name[col]].value
            A_band = R_band*ebv
            self.jy[k] *= 10**(0.4*A_band)

        return


