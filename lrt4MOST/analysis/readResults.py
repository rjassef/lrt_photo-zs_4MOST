import numpy as np
import SED_Model
from Star_Model import StarModel
import pathlib
import subprocess

class ReadResults(object):

    def __init__(self, ztype="zspec", catname=None, catname_noagn=None):

        #Save the ztype requested.
        self.ztype = ztype

        #Process star fits differently.
        if self.ztype=='star':
            if catname is None:
                catname = "combined_star_fit.dat"

            #Read the header
            self.star_set_name = dict()
            cat = open(catname)
            x = cat.readline().split()
            for i in range(int(len(x)/2)):
                self.star_set_name[x[i*2+1]] = x[i*2]
            cat.close()

            #Read the rest of the file.
            data = np.loadtxt(catname, skiprows=1)
            self.id        = data[:,0]
            self.chi2      = data[:,1]
            self.star_set  = data[:,2]
            self.star_temp = data[:,3]
            self.comp      = data[:,4:6]
            del(data)
        else:
            #In order to select the AGN we need to save the data from the SED_fit run as well as the chi squared from the SED_fit_noagn run. 
            if catname is None:
                catname = "SED_fit_{}_all.dat".format(ztype)
            data = np.loadtxt(catname)
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

            if catname_noagn is None:
                catname_noagn = "SED_fit_noagn_{}_all.dat".format(ztype)
            self.chi2_noagn = np.loadtxt(catname_noagn, usecols=[3])

        return

    def plot(self, id, phot, hardcopy=False, run_fit=False, folder="SED_plots"):

        #Match. Alert if there is more than one match. 
        k = np.where(self.id==id)
        if len(k)>1:
            print("There are {} matches to id {}. Using the first one".format(len(k), id))
        kuse = k[0][0]

        #Create the lrt SED Model object. 
        if self.ztype == 'star':
            st_num = "{0:.0f}".format(self.star_set[kuse])
            obj = StarModel(st=self.star_set_name[st_num])
        else:
            obj = SED_Model.lrt_model()

        #Populate the model. 
        obj.jy = phot.jy[:,kuse]
        obj.ejy = phot.ejy[:,kuse]
        obj.jyuse = phot.jyuse[:,kuse]
        obj.ejy[np.isnan(obj.ejy)] = 0.
        obj.name = phot.id[kuse]
        if self.ztype!='star':
            #obj.zspec = self.z[kuse]
            exec("obj.{} = {}".format(self.ztype,self.z[kuse]))

        if run_fit:
            if self.ztype=='star':
                obj.fit()
            else:
                obj.kc_fit()
        else:
            obj.comp = self.comp[kuse]
            obj.jymod = phot.jy[:,kuse]
            if self.ztype=='star':
                obj.chi2 = self.chi2[kuse]             
                obj.tfit = self.star_temp[kuse]
                obj.ns_best = int(np.floor(self.star_temp[kuse]))
            else:
                obj.ebv = self.ebv[kuse]
                obj.igm = self.igm[kuse]
                obj.chi2 = self.chi2_agn[kuse]

        #Plot it. 
        if hardcopy:
            if not pathlib.Path(folder).exists():
                subprocess.call(["mkdir", folder])
            obj.plot_to_file(folder+"/Obj_{}_{}.png".format(id, self.ztype))
        else:
            obj.plot()

        return

