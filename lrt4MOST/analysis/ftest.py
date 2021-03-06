import numpy as np
from scipy.special import betainc
import pathlib

from .agnSelection import AGNSelection

class Ftest(AGNSelection):

    def __init__(self, ztype="zspec"):
        self.ztype = ztype
        return

    def getF(self, res, outname=None, force=False):

        #Output name
        if outname is None:
            outname = "Ftest_{0:s}.dat".format(self.ztype)
        if not force and pathlib.Path(outname).exists():
            print("Output file already exists. Skipping calculation.")
            return

        #Get the degrees of freedom.
        nu = res.nb - 6.0
        nu_noagn = res.nb - 4.0
        chi2 = res.chi2_agn
        chi2_noagn = res.chi2_noagn

        self.F = np.where((nu>0.) & (chi2_noagn>chi2),
            ((chi2_noagn-chi2)/(2.0)) / ((chi2+1e-32)/(nu+1e-32)),
            0.)
        nnu1 = nu_noagn - nu
        nnu2 = np.where(nu>0., nu, 1.)
        w = nnu1*self.F/(nnu1*self.F+nnu2)
        self.p = 1.-betainc(nnu1/2.,nnu2/2.,w)


        np.savetxt(outname, np.array([self.F,self.p]).T)

        return

    def readF(self, fname=None):
        if fname is None:
            fname = "Ftest_{0:s}.dat".format(self.ztype)
        data = np.loadtxt(fname)
        self.F = data[:,0]
        self.p = data[:,1]
        del(data)
        return


    def selectAGN(self, res, pmax=0.1, chi2_max=100, fname=None, zphot_min_cut=True):

        #Output catalog name.
        if fname is None:
            fname = "AGN_Ftest_{}.dat".format(self.ztype)

        #Run the selection.
        kuse = (self.p<pmax) 

        #Generate and save the list of objects.
        self.savelist(fname, kuse, res, chi2_max, zphot_min_cut)

        return
        

