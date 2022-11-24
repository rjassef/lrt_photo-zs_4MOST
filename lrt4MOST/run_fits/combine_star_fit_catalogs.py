import numpy as np 
import pathlib

def combine_star_fit_catalogs(stypes, output_catalog="combined_star_fit.dat"):

    if pathlib.Path(output_catalog).exists():
        print("Catalog already exists. Skipping step")
        return

    #Write the header of the catalog. 
    cat = open(output_catalog,"w")
    for k, stype in enumerate(stypes):
        cat.write("{} {} ".format(stype,k+1))
    cat.write("\n")

    #Read the input catalogs and combine them by using the type with the lowest chi squared.
    for k, stype in enumerate(stypes):
        data = np.loadtxt("star_fit_{}_all.dat".format(stype))

        if k==0:
            final_data = np.zeros((data.shape[0],6))
            final_data[:,0] = data[:,0] #IDs
            final_data[:,1] = np.inf #Chi2
            continue

        cond = final_data[:,1]>data[:,3]
        final_data[cond,1] = data[cond,3]
        final_data[cond,2] = k+1 #Stellar sequence
        final_data[cond,3] = data[cond,2] + data[cond,5]/(data[cond,4]+data[cond,5]) #Template number.
        final_data[cond,4] = data[cond,4]
        final_data[cond,5] = data[cond,5]

    #Flag unfittable sources.
    cond = np.isnan(final_data[:,3])
    final_data[cond,2:6] = -1.

    #Save the file.
    np.savetxt(cat,final_data,fmt="%20d %10.3e %5d %10.2f %10.3e %10.3e")
    cat.close()

    return
