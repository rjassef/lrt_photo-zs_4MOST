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
LIBS = -L$(LRTPATH) -llrt -llapack

all: $(EXEC)

$(EXEC): $(SRC)
	$(FC) -o $@ $@.f $(FCFLAGS) $(LIBS)

clean: 
	$(RM) $(EXEC)

