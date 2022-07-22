import subprocess
import os
import multiprocessing
import numpy as np

def Runfitzero(phot, catname="fitzero_sample.txt", ncpu=None, nobj_max=None):

    #First, write the fitzero sample file.
    #Write the fitzero file. 
    cato = open(catname,"w")
    if nobj_max is None:
        nobj_max = len(phot.id)
    for k in range(nobj_max):
        cato.write("{0:10.4f} ".format(phot.zspec[k]))
        for j in range(phot.nchan):
            cato.write("{0:10.3e} ".format(phot.jy[j,k]))
        for j in range(phot.nchan):
            cato.write("{0:10.3e} ".format(phot.ejy[j,k]))    
        for j in range(phot.nchan):
            cato.write("{0:5.0f} ".format(phot.jyuse[j,k]))
        cato.write("\n")
    cato.close()

    #Set the number of CPUs to use.
    if ncpu is None:
        ncpu = multiprocessing.cpu_count()-2
    my_env = os.environ.copy()
    my_env['OMP_NUM_THREADS'] = "{}".format(ncpu)

    #Erase the previous channel.zpc file if it exists.
    subprocess.call(["rm","-f","channel.zpc"])

    #Call the fitzero program. 
    subprocess.call(["{0:s}/lrt4MOST/fortran_codes/run_fitzero".format(os.environ.get('LRT4MOST_LOC')),catname],env=my_env)

    #Now, make the corrections file.
    corr = np.loadtxt("channel.zpc")
    cato = open("zero_point_corrections.dat","w")
    for k,col in enumerate(phot.photcols):
        cato.write("{} {}\n".format(col,corr[k]))
    cato.close()

    #Erase the channel.zpc file.
    subprocess.call(["rm","-f","channel.zpc"])

    return