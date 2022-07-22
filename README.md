## Instructions and Description

This is a package that implements the AGN SED selection for ChANGES, the S16 group of 4MOST, based on the LRT templates of [Assef et al. (2010)](https://ui.adsabs.harvard.edu/abs/2010ApJ...713..970A/abstract).

### Installation

After cloning this repository, you will need to install the LRT software. The selection can be run with the latest stable release of the soft which you can find [here](https://astro.udp.cl/~rjassef/sed_templates.html). Additional capabilities, like plotting the SEDs, need the Python wrappers for the LRT libraries, which require the latest unstable release of the code. You can find it in this [repo](https://github.com/rjassef/LRT). I recommend using the repository version. 

Once you have installed the LRT libraries, you need to compile the Fortran codes used by this package. 

    cd lrt4MOST/fortran_codes
    make all

Finally, you either need to add this package to your PYTHONPATH. You can do it, for example in you .bashrc or .tcshrc files. Alternatively, if you use either bash or zsh, you can go to the root folder of this code and execute

    source setup_bash

### Example to Select AGN

We assume that we start either from the pre-processed catalog produced by Alejandra for CIGALE. 

First, we read the photometric catalog using the [ReadPhot](lrt4MOST/run_fits/readPhot.py) function: 

```python
from lrt4MOST import ReadPhot
phot = ReadPhot("catalog.fits")
```

If starting from Victoria's catalog, use [ReadPhotCtype2](lrt4MOST/run_fits/readPhotCtype2.py) instead of [ReadPhot](lrt4MOST/run_fits/readPhot.py).

Then, we calculate the photometric redshifts using the [GetPhotozs](lrt4MOST/run_fits/getPhotozs.py) class, allowing first for all templates, and the only for the galaxy one. 

```python
from lrt4MOST import GetPhotozs
zphot = GetPhotozs(with_AGN=True)
zphot.run(phot)
zphot_noagn = GetPhotozs(with_AGN=False)
zphot_noagn.run(phot)
```

Next, we get the full SED modeling using the [GetSEDFits](lrt4MOST/run_fits/getSEDFits.py) class (done for technical reasons in a different step than the photo-z fitting). Assuming we do not have spec-zs, 

```python
from lrt4MOST import GetSEDFits
SED_fits = GetSEDFits(with_AGN=True, ztype="zphot")
SED_fits.run(phot)
SED_fits_noagn = GetSEDFits(with_AGN=False, ztype="zphot")
SED_fits_noagn.run(phot)
```

To select the AGN from the sample, we can use an F-test, implemented in the [Ftest](lrt4MOST/analysis/ftest.py) class. 

```python 
from lrt4MOST import Ftest
ftest_zphot = Ftest(ztype="zphot")
ftest_zphot.getF(res_zphot)
ftest_zphot.selectAGN(res_zphot)
```

We can also use the BIC selection implemented in the [BIC](lrt4MOST/analysis/bic.py) class:

```python
from lrt4MOST import BIC
bic_zphot = BIC(ztype="zphot")
bic_zphot.selectAGN(res_zphot)
```

Or we can use a selection based on $\hat{a}$, the fraction of luminosity of the best-fit SED coming from the AGN template. This is implemented in the [ahatSelection](lrt4MOST/analysis/ahatSelection.py) class.

```python
from lrt4MOST import ahatSelection
ahat_sel_zphot = ahatSelection(ztype="zphot")
ahat_sel_zphot.selectAGN(res_zphot)
```
