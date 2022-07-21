#Fortran 77 Compiler
FC = gfortran

#Compiler Flags
FCFLAGS = -ffixed-line-length-none -std=legacy

#Removal Tool
RM = /bin/rm -f

#Source files
SRC = zphot.f zphot_noagn.f SED_fit.f SED_fit_noagn.f

#Executables.
EXEC :=$(SRC:.f=)

#Libraries
LIBS = restart.f -L$(LRTPATH) -llrt -llapack

all: $(EXEC) run_fitzero

$(EXEC): $(SRC)
	$(FC) -o $@ $@.f $(FCFLAGS) $(LIBS)

run_fitzero: run_fitzero.f fitzero_omp.f
	$(FC) -o run_fitzero run_fitzero.f fitzero_omp.f $(FCFLAGS) $(LIBS) -fopenmp

clean: 
	$(RM) $(EXEC)

