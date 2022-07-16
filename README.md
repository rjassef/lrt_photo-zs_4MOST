## Instructions and Description

This is a package that implements the AGN SED selection for ChANGES, the S16 group of 4MOST, based on the LRT templates of [Assef et al. (2010)](https://ui.adsabs.harvard.edu/abs/2010ApJ...713..970A/abstract).

### Installation

After cloning this repository, you will need to install the LRT software. The selection can be run with the latest stable release of the soft which you can find [here](https://astro.udp.cl/~rjassef/sed_templates.html). Additional capabilities, like plotting the SEDs, need the Python wrappers for the LRT libraries, which require the latest unstable release of the code. You can find it in this [repo](https://github.com/rjassef/LRT). I recommend using the repository version. 

Once you have installed the LRT libraries, you need to compile the Fortran codes used by this package. 

    cd lrt4MOST/fortran_codes
    make all

Finally, you either need to add this package to your PYTHONPATH. You can do it, for example in you .bashrc or .tcshrc files. Alternatively, if you use either bash or zsh, you can go to the root folder of this code and execute

    source setup_bash
