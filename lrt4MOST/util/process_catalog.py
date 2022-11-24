from .. import ReadPhot, GetPhotozs, GetSEDFits, GetStarFits, combine_star_fit_catalogs

def process_catalog(catalog_file_name, nobj_per_thread = 50000, ncpu=None):

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

    #Combine into final catalog.
    print('Generating combined stellar fits catalog...')
    combine_star_fit_catalogs(stypes)

    return
