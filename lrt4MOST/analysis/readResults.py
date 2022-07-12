import numpy as np

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
