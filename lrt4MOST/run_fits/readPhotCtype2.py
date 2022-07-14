import numpy as np
from astropy.io import fits
import astropy.units as u
import re
import warnings

from .procPhot import ProcPhot

class ReadPhotCtype2(ProcPhot):

    def __init__(self, fname, n=1, idcol='ls_id', zcol='redshift', input_unit='abmag', output_flux_unit='Jy', noise_floor=0.1, jyzero=None):

        #Open the file. 
        tab = fits.open(fname, format=format)

        #Save the photometry columns. We recognize them as they have the band name preceeded by a err_. 
        self.photcols = list()
        self.magtype = list()
        colnames = tab[n].data.dtype.names
        for k, col in enumerate(colnames):
            m1 = re.match("(.*?)mag_(.*)", colnames[k-1])
            if hasattr(m1,'group') and re.match(m1.group(1)+"magerr_"+m1.group(2), col):
                self.photcols.append(m1.group(2))
                self.magtype.append(m1.group(1))
        self.nchan = len(self.photcols)

        #Now, save the photometry, the IDs and the redshifts.
        self.id    = tab[1].data[idcol]
        self.nobj  = len(self.id)
        if zcol in tab[1].data.dtype.names:
            self.zspec = tab[1].data[zcol]
        else:
            self.zspec = -1. * np.ones(self.id.shape)
        self.jy    = np.zeros((self.nchan, self.nobj))
        self.ejy   = np.zeros(self.jy.shape)

        flux_factor = (u.Jy/u.Unit(output_flux_unit)).to(1)
        if input_unit is 'abmag':
            jyzero = 3631.*flux_factor * np.ones(self.jy.shape[0])
        with warnings.catch_warnings():
            # Ignore the warning about an invalid value in less. 
            warnings.filterwarnings('ignore')
            for k, col in enumerate(self.photcols):
                cname = self.magtype[k]+"mag_"+col
                cerrname = self.magtype[k]+"magerr_"+col
                mag = np.where(tab[n].data[cname]<0, np.nan, tab[n].data[cname])
                magerr = np.where(tab[n].data[cerrname]<0, np.nan, tab[n].data[cerrname])
                self.jy[k]  = jyzero[k] * 10**(-0.4*mag)
                self.ejy[k] = np.where(np.isnan(self.jy[k]), 
                    jyzero[k] * 10**(-0.4*magerr),
                    0.4*np.log(10.) * self.jy[k] * tab[n].data[cerrname]
                )

        #Close the table.
        tab.close()

        #Run the rest of the photometry process initiation. 
        super(ReadPhotCtype2,self).__init__(noise_floor=noise_floor)

        return
