# List of Improvements Needed

- [ ] Implement foreground reddening corrections
   - *Importance*: Critical
   - *Difficulty*: Medium
   - *Status*: Mostly implemented. The remaining issue is that we are currently using the coefficients from Table 6 of Schlafly et al. (2006)

- [ ] Fix bug in script to combine stellar fits. Because of round-off, it is not always possible to plot the correct fit using the effective stellar template number. 
   - *Importance*: Medium
   - *Difficulty*: Low
   - *Required changes*: Add an extra column with the reference template number for SED plotting. For example, an object that has 0.1% of template 6 and 99.9% of template 7 will be written with an effective template number of 7.00, but we need the reference template number of 6 to do SED = T6*0.001 + T7*0.999

- [ ] Create documentation
   - *Importance*: High
   - *Difficulty*: Low
   - *Required changes*: Add documentation to all python classes and subclasses. 

- [ ] Only use fits files to save space.
   - *Importance*: Medium
   - *Difficulty*: Medium
   - *Required changes*: Modify fortran codes and data reading structures to only read fits files. Do not use the astropy.table module as the memory usage is too large.

- [ ] Avoid loading data files to conserve RAM. 
   - *Importance*: Unclear until full sample to process is agreed on.
   - *Difficulty*: Hard
   - *Required changes*: All codes would need work. Some may need to be re-written in C or Fortran to allow for this without affecting execution speed.
        
- [ ] Add photo-zs to final catalogs. Add also reddening and Lagn/Lhost.
