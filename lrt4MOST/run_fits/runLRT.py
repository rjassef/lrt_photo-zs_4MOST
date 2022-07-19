import numpy as np
import multiprocessing
import subprocess
import os 
from time import sleep
import pathlib

class RunLRT(object):

    def __init__(self, fit_type, with_AGN, ztype):

        #Save input flags.
        self.fit_type = fit_type
        self.with_AGN = with_AGN
        self.ztype = ztype

        #Determine the code to use.
        if self.with_AGN:
            self.code = self.fit_type
        else:
            self.code = self.fit_type+"_noagn"
        
        return
    
    def run(self, phot, ncpu=None, tsleep=2, force=False, reset=False, nobj_per_thread=None):

        #Output catalog name. If it already exist, and execution is not forced, do not continue.
        if self.fit_type=="SED_fit":
            catname = "{}_{}_all.dat".format(self.code, self.ztype)
        else:
            catname = "{}_all.dat".format(self.code)
        if not force and pathlib.Path(catname).exists():
            print("Output file already exists. Skipping calculation.")
            return

        #Set the number of CPUs to use if not provided.
        if ncpu is None:
            ncpu = multiprocessing.cpu_count()-2

        #Convenience variables.
        ntot = phot.nobj
        nchan = phot.nchan

        #If needed to, read the photo-zs.
        zphots = None
        if self.fit_type=="SED_fit" and self.ztype=="zphot":
            if self.with_AGN:
                zphots = np.loadtxt("zphot_all.dat",usecols=[2])
            else:
                zphots = np.loadtxt("zphot_noagn_all.dat",usecols=[2])

        #Create the ncpu catalogs and run the code in each of them separately in parallel.
        if nobj_per_thread is None:
            ncat = ncpu
        else:
            ncat = int(ntot/nobj_per_thread)
        if ncat>=1e6:
            print("Cannot run more than a million sub catalogs without loosing track of the order. Please modify code.")
            exit()
        nobjs = np.ones(ncat, dtype=np.int32) * np.int32(ntot/ncat)
        if np.sum(nobjs)<ntot:
            nobjs[-1] += ntot-np.sum(nobjs)

        #Run all of them, making sure there aren't more than ncpu threads running at a time.
        p = list()
        for n in range(ncat):

            #Check that we are not exceeding the max amounts of threads. Wait here until a thread opens. 
            while self.n_thread_active(p)==ncpu:
                sleep(tsleep)

            #Set the file names. 
            finput  = "{0:s}_input_{1:06d}.dat".format(self.code,n)
            foutput = "{0:s}_output_{1:06d}.dat".format(self.code,n)

            #Check if the input file already exists. Since we may be restarting from a stopped process, we may not need to rewrite it. 
            if reset or not pathlib.Path(finput).exists():

                #Set the boundaries of the main array.
                if n==0:
                    k1 = 0
                else:
                    k1 = np.sum(nobjs[:n])
                k2 = np.sum(nobjs[:n+1])

                #Set the output array.
                out = np.zeros((nobjs[n],2+nchan*3))
                out[:,0] = phot.id[k1:k2]
                out[:,1] = phot.zspec[k1:k2]
                out[:,        2:  nchan+2] = phot.jy[:, k1:k2].T
                out[:,  nchan+2:2*nchan+2] = phot.ejy[:, k1:k2].T
                out[:,2*nchan+2:3*nchan+2] = phot.jyuse[:, k1:k2].T

                #If we are using zphots, then add those to the data array instead. 
                if zphots is not None:
                    out[:,1] = zphots[k1:k2]

                #Write the data file.
                cato = open(finput,"w")
                cato.write("{0:d} {1:d}\n".format(nobjs[n],nchan))
                np.savetxt(cato,out,fmt='%15.0f %15.4f'+'%15.3e'*(2*nchan)+'%3.0f'*nchan)
                cato.close()

            #If we are forcing the reset, then erase the outputfile as well.
            if reset:
                subprocess.call(["rm","-f",foutput])

            p.append(subprocess.Popen("{0:s}/lrt4MOST/fortran_codes/{1:s} {2} {3}".format(os.environ.get('LRT4MOST_LOC'),self.code,finput, foutput),shell=True))

        #Wait for all the threads to finish
        while self.n_thread_active(p)>0:
            sleep(tsleep)

        #Finally, combine the output files into one and remove the intermediary file.
        subprocess.call("cat {0:s}_output* > {1:s}".format(self.code, catname), shell=True)
        subprocess.call("rm {0:s}_input_*.dat {0:s}_output_*.dat".format(self.code), shell=True)

        return

    def n_thread_active(self, p):
        i = 0
        for pp in p:
            if pp.poll() is None:
                i+=1
        return i
