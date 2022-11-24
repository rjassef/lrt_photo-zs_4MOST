import numpy as np
from scipy.special import betainc
import pathlib

from .agnSelection import AGNSelection

class Ftest(AGNSelection):

    def __init__(self, ztype="zspec"):
        self.ztype = ztype
        return

    def getF(self, res, savefile=True, outname=None, force=False):

        #Output name
        if outname is None:
            outname = "Ftest_{0:s}.dat".format(self.ztype)
        if not force and pathlib.Path(outname).exists():
            print("Output file already exists. Skipping calculation.")
            self.readF(fname=outname)
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

        #To recognize more easily the sources for which we could not calculate F due to the lack of photometric bands, let's change the F values to -1 for those sources. 
        self.F[nu<=0.] = -1.0

        if savefile:
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


    def selectAGN(self, res, pmax=0.1, chi2_max=100, fname=None, zphot_min_cut=False, res_stars=None, select_stars_function=None, save_list=True):

        #Attempt to read or calculate the F values. 
        self.getF(res)
        # if not hasattr(self, 'F'):
        #     self.readF()

        #Output catalog name.
        if fname is None:
            fname = "AGN_Ftest_{}.dat".format(self.ztype)

        #Run the selection.
        self.k_agn = (self.p<pmax) & (self.is_not_star(res, res_stars, select_stars_function))

        #Generate and save the list of objects.
        if save_list:
            self.savelist(fname, self.k_agn, res, chi2_max, zphot_min_cut)

        return
        

