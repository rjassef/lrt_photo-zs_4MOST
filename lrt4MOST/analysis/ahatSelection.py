import numpy as np
import pathlib

from .agnSelection import AGNSelection

class ahatSelection(AGNSelection):

    def __init__(self, ztype="zspec"):
        self.ztype = ztype
        return

    def selectAGN(self, res, ahat_min=0.5, chi2_max=100, fname=None, zphot_min_cut=False, res_stars=None, select_stars_function=None):

        #Output catalog name.
        if fname is None:
            fname = "AGN_ahat_{}.dat".format(self.ztype)

        #Run the selection.
        kuse = (res.ahat>ahat_min) & (self.is_not_star(res, res_stars, select_stars_function))

        #Generate and save the list of objects.
        self.savelist(fname, kuse, res, chi2_max, zphot_min_cut)

        return
        

