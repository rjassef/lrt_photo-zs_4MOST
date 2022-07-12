import numpy as np
import SED_Model

class ReadResults(object):

    def __init__(self, ztype="zspec"):

        #Save the ztype requested.
        self.ztype = ztype

        #In order to select the AGN we need to save the data from the SED_fit run as well as the chi squared from the SED_fit_noagn run. 
        data = np.loadtxt("SED_fit_{}_all.dat".format(ztype))
        self.id = data[:,0]
        self.nb = data[:,1]
        self.z  = data[:,2]
        self.chi2_agn = data[:,3]
        self.ahat = data[:,4]
        self.ebv  = data[:,5]
        self.igm  = data[:,6]
        self.comp = data[:,7:11]
        self.vec  = data[:,11:15]
        del(data)

        self.chi2_noagn = np.loadtxt("SED_fit_noagn_{}_all.dat".format(ztype), usecols=[3])

        return

    def plot(self, id, phot, hardcopy=False, run_fit=False):

        #Create the lrt SED Model object. 
        gal = SED_Model.lrt_model()

        #Match. Alert if there is more than one match. 
        k = np.where(self.id==id)
        if len(k)>1:
            print("There are {} matches to id {}. Using the first one".format(len(k), id))
        kuse = k[0][0]

        #Populate the model. 
        gal.jy = phot.jy[:,kuse]
        gal.ejy = phot.ejy[:,kuse]
        gal.jyuse = phot.jyuse[:,kuse]
        gal.ejy[np.isnan(gal.ejy)] = 0.
        gal.zspec = self.z[kuse]

        if run_fit:
            gal.kc_fit()
        else:
            gal.comp = self.comp[kuse]
            gal.ebv = self.ebv[kuse]
            gal.igm = self.igm[kuse]
            gal.jymod = phot.jy[:,kuse]
            gal.chi2 = self.chi2_agn[kuse]

        #Plot it. 
        if hardcopy:
            gal.plot_to_file("Obj_{}_{}.png".format(id, self.ztype))
        else:
            gal.plot()

        return

