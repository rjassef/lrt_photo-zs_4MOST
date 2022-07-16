# List of Improvements Needed

- [ ] Only use fits files to save space.
        * *Importance*: Medium
        * *Difficulty*: Medium
        * *Required changes*: Modify fortran codes and data reading structures to only read fits files. Do not use the astropy.table module as the memory usage is too large.

- [ ] Avoid loading data files to conserve RAM. 
        * *Importance*: Unclear until full sample to process is agreed on.
        * *Difficulty*: Hard
        * *Required changes*: All codes would need work. Some may need to be re-written in C or Fortran to allow for this without affecting execution speed.
        