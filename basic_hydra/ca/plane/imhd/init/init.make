 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
random: $(objects) $(fft_lib) $(sourcedir)/init/random.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/random.f90 -o random $(flags)

vstrip: $(objects) $(sourcedir)/init/vstrip.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/vstrip.f90 -o vstrip $(flags)

ellipse: $(objects) $(sourcedir)/init/ellipse.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/ellipse.f90 -o ellipse $(flags)

gauss: $(objects) $(sourcedir)/init/gauss.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/gauss.f90 -o gauss $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


