import numpy as np
import pathlib

from .agnSelection import AGNSelection

class BIC(AGNSelection):

    def __init__(self, ztype="zspec"):
        self.ztype = ztype
        return

    def selectAGN(self, res, BIC_min=10., chi2_max=100, fname=None, zphot_min_cut=True):

        #Output catalog name.
        if fname is None:
            fname = "AGN_BIC_{}.dat".format(self.ztype)

        #Calculate the BIC index.
        self.BIC = res.chi2_noagn - res.chi2_agn + 2*np.log(res.nb)
        np.savetxt("BIC_{}.dat".format(self.ztype), self.BIC)

        #Run the selection.
        kuse = (self.BIC>BIC_min) 

        #Generate and save the list of objects.
        self.savelist(fname, kuse, res, chi2_max, zphot_min_cut)

        return
        

