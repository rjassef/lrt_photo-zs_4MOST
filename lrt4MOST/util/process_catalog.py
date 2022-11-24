import numpy as np
import astropy.units as u 
from astropy.table import Table
import pathlib

from .. import ReadPhot, GetPhotozs, GetSEDFits, GetStarFits, combine_star_fit_catalogs
from .. import ReadResults, Ftest, BIC

def process_catalog(catalog_file_name, nobj_per_thread = 50000, ncpu=None, mag_save_bname='r_prime', fout_name=None):

    #Final output catalog. Skip if it exists already.
    if fout_name is None:
        fout_name = "fullproc_"+catalog_file_name
    if pathlib.Path(fout_name).exists():
        print("Catalog fully processed already. Skipping all operationds.")
        return

    #Read the data
    print('Reading photometry for AGN/Gal SED fits...')
    phot = ReadPhot(catalog_file_name)

    #Get the photo-zs with and without the AGN template.
    print('Getting AGN/Gal photo-zs...')
    zphot_calc = GetPhotozs(with_AGN=True)
    zphot_calc.run(phot, ncpu=ncpu, nobj_per_thread=nobj_per_thread)

    print('Getting Gal only photo-zs...')
    zphot_noagn_calc = GetPhotozs(with_AGN=False, zmax=2.0)
    zphot_noagn_calc.run(phot, ncpu=ncpu, nobj_per_thread=nobj_per_thread)

    #Fit the SEDs with z_phot
    print('Running full AGN/Gal SED fit...')
    SED_fit_zphot = GetSEDFits(with_AGN=True, ztype="zphot")
    SED_fit_zphot.run(phot, ncpu=ncpu, nobj_per_thread=nobj_per_thread)

    print('Running full Gal only SED fit...')
    SED_fit_zphot_noagn = GetSEDFits(with_AGN=False, ztype="zphot")
    SED_fit_zphot_noagn.run(phot, ncpu=ncpu, nobj_per_thread=nobj_per_thread)

    #Remove the photometry file from memory, and re-read the photometry by without foreground reddening corrections.
    print('Reading photometry for Stellar fits...')
    del phot
    phot2 = ReadPhot(catalog_file_name, fg_dust_corr=False)

    #Fit as stars. First MS, then GS, then SGS and finally BDs.
    stypes = ["MS","GS","SGS","BDs"]
    for stype in stypes:
        print('Running {} stellar SED fits...'.format(stype))
        stars = GetStarFits(stype)
        stars.run(phot2, ncpu=ncpu, nobj_per_thread=nobj_per_thread)

    #Combine into final stellar catalog.
    print('Generating combined stellar fits catalog...')
    combine_star_fit_catalogs(stypes)

    #Save the requested magnitudes. 
    k = np.argwhere(phot2.photcols==mag_save_bname)
    zp = (3631.*u.Jy).to(phot2.output_flux_unit).value
    mag = -2.5*np.log10(phot2.jy[k]/zp)
    mag[np.isnan(mag)] = 99.

    #Load the fit results.
    print("Reading AGN/Gal fit results...")
    res_agngal = ReadResults(ztype='zphot', catname="SED_fit_zphot_all.dat", catname_noagn="SED_fit_noagn_zphot_all.dat")

    print("Reading Stellar fit results...")
    res_stars  = ReadResults(ztype='star', catname="combined_star_fit.dat")

    #Calculate the BIC index.
    BIC = BIC(ztype="zphot")
    BIC.getBIC(res_agngal)

    #Calculate the Fprob
    Fp = Ftest(ztype="zphot")
    Fp.getF(res_agngal)

    #Save the processed catalog.
    proc_table = Table()
    proc_table['id'] = phot2.id
    proc_table['ra'] = phot2.ra
    proc_table['dec'] = phot2.dec
    proc_table['mag'] = mag
    proc_table['nb'] = res_agngal.nb
    proc_table['chi2_agngal'] = res_agngal.chi2_agn
    proc_table['chi2_gal'] = res_agngal.chi2_noagn
    proc_table['chi2_star'] = res_stars.chi2
    proc_table['BIC'] = BIC.BIC
    proc_table['Fp'] = Fp.p
    proc_table.write(fout_name, format='fits')

    return
