import subprocess
import os

def Runfitzero(phot, catname="fitzero_sample.txt"):

    #First, write the fitzero sample file.
    #Write the fitzero file. 
    cato = open(catname,"w")
    for k in range(len(phot.id)):
        cato.write("{0:10.4f} ".format(phot.zspec[k]))
        for j in range(phot.nchan):
            cato.write("{0:10.3e} ".format(phot.jy[j,k]))
        for j in range(phot.nchan):
            cato.write("{0:10.3e} ".format(phot.ejy[j,k]))    
        for j in range(phot.nchan):
            cato.write("{0:5.0f} ".format(phot.jyuse[j,k]))
        cato.write("\n")
    cato.close()

    #Erase the previous channel.zpc file if it exists.
    subprocess.call(["rm","-f","channel.zpc"])

    #Finallys, call the fitzero program. 
    subprocess.call("{0:s}/lrt4MOST/fortran_codes/run_fitzero {1}".format(os.environ.get('LRT4MOST_LOC'),catname),shell=True)

    return