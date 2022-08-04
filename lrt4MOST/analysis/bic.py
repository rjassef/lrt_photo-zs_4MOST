import numpy as np
import pathlib

from .agnSelection import AGNSelection

class BIC(AGNSelection):

    def __init__(self, ztype="zspec"):
        self.ztype = ztype
        return

    def selectAGN(self, res, BIC_min=10., chi2_max=100, fname=None, zphot_min_cut=False, res_stars=None, select_stars_function=None):

        #Output catalog name.
        if fname is None:
            fname = "AGN_BIC_{}.dat".format(self.ztype)

        #Calculate the BIC index.
        self.BIC = res.chi2_noagn - res.chi2_agn + 2*np.log(res.nb)
        np.savetxt("BIC_{}.dat".format(self.ztype), self.BIC)

        #Run the selection.
        self.k_agn = (self.BIC>BIC_min) & (self.is_not_star(res, res_stars, select_stars_function))

        #Generate and save the list of objects.
        self.savelist(fname, self.k_agn, res, chi2_max, zphot_min_cut)

        return
        

