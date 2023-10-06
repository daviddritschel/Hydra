 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:

patch: $(objects) $(fft_lib) $(sourcedir)/init/patch.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/patch.f90 -o patch $(flags)

twovort: $(objects) $(fft_lib) $(sourcedir)/init/twovort.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/twovort.f90 -o twovort $(flags)

ranvor: $(objects) $(fft_lib) $(sourcedir)/init/ranvor.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/ranvor.f90 -o ranvor $(flags)

bickley: $(objects) $(fft_lib) $(sourcedir)/init/bickley.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/bickley.f90 -o bickley $(flags)

rest-state: $(objects) $(sourcedir)/init/rest-state.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/rest-state.f90 -o rest-state $(flags)

staircase: $(objects) $(sourcedir)/init/staircase.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/staircase.f90 -o staircase $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)
