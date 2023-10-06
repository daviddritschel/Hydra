 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
dipole: $(objects) $(sourcedir)/init/dipole.f90
	$(f90) parameters.o constants.o contours.o generic.o $(sourcedir)/init/dipole.f90 -o dipole $(flags)

# Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


