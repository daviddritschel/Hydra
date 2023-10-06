 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
modon: $(objects) $(sourcedir)/init/modon.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/modon.f90 -o modon $(flags)

straka: $(objects) $(sourcedir)/init/straka.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/straka.f90 -o straka $(flags)

slug: $(objects) $(sourcedir)/init/slug.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/slug.f90 -o slug $(flags)

tide: $(objects) $(sourcedir)/init/tide.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/tide.f90 -o tide $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


