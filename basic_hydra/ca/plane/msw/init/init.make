 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
ranpv: $(objects) $(fft_lib) $(sourcedir)/init/ranpv.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/ranpv.f90 -o ranpv $(flags)

gauss: $(objects) $(fft_lib) $(sourcedir)/init/gauss.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/init/gauss.f90 -o gauss $(flags)

vstrip: $(objects) $(sourcedir)/init/vstrip.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/vstrip.f90 -o vstrip $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


