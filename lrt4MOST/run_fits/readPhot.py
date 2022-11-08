import numpy as np
from astropy.io import fits
import astropy.units as u

from .procPhot import ProcPhot

class ReadPhot(ProcPhot):

    def __init__(self, fname, n=1, idcol='id', zcol='redshift', input_unit='mJy', output_flux_unit='Jy', noise_floor=0.1, flux_min=None, flux_min_band=None):

        #Open the file. 
        tab = fits.open(fname)

        #Save the photometry columns. We recognize them as they have a name followed by a _err. 
        self.photcols = list()
        colnames = tab[n].data.dtype.names
        for k, col in enumerate(colnames):
            if col==colnames[k-1]+"_err":
                self.photcols.append(colnames[k-1])
        self.photcols = np.array(self.photcols)
        self.nchan = len(self.photcols)

        if flux_min is None:
            iuse = np.ones(tab[1].data[idcol].shape, dtype=bool)
        else:
            iuse = (~np.isnan(tab[1].data[flux_min_band])) & (tab[1].data[flux_min_band]>=flux_min.to(u.Unit(input_unit)).value)


        #Now, save the photometry, the IDs and the redshifts.
        self.id    = tab[1].data[idcol][iuse]
        self.nobj  = len(self.id)
        self.zspec = tab[1].data[zcol][iuse]
        self.jy    = np.zeros((self.nchan, self.nobj))
        self.ejy   = np.zeros(self.jy.shape)
        flux_factor = (u.Unit(input_unit)/u.Unit(output_flux_unit)).to(1)
        for k, col in enumerate(self.photcols):
            self.jy[k]  = tab[n].data[col][iuse] * flux_factor
            self.ejy[k] = tab[n].data[col+"_err"][iuse] * flux_factor

        #Close the table.
        tab.close()

        #Run the rest of the photometry process initiation. 
        super(ReadPhot,self).__init__(noise_floor=noise_floor)

        return
