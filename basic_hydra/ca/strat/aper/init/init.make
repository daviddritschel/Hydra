 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
dam-break: $(fft_lib) $(objects) $(sourcedir)/init/dam-break.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/dam-break.f90 -o dam-break $(flags)

slug: $(fft_lib) $(objects) $(sourcedir)/init/slug.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/slug.f90 -o slug $(flags)

vort: $(fft_lib) $(objects) $(sourcedir)/init/vort.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/vort.f90 -o vort $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


