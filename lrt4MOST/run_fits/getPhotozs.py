import numpy as np
import subprocess

from .runLRT import RunLRT

class GetPhotozs(RunLRT):

    def __init__(self, with_AGN=True, with_prior=False):
        if with_prior:
            self.set_prior()
        else:
            subprocess.call(["rm","-rf","prior.dat"])
        super(GetPhotozs,self).__init__(fit_type="zphot", with_AGN=with_AGN, ztype="zphot")
        return

    def set_prior(self):
        '''Use the Lin et al. (1996) r-band GLF as a prior. This is a relatively light prior, so it does not matter too much whether the band is exactly theh same r as from Lin et al. '''

        #Find the band id. 
        cat = open("bandmag.dat")
        for k,line in enumerate(cat):
            x = line.split()
            if x[0]=='DECAM_r':
                break
        cat.close()

        #Now, create the prior file. 
        cato = open("prior.dat","w")
        cato.write("1\n")
        cato.write("{}\n".format(k+1))
        cato.write("-21.4\n")
        cato.write("-0.7\n")
        cato.close()

        return

