 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:

inigamma: $(objects) $(fft_lib) $(sourcedir)/init/inigamma.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/inigamma.f90 -o inigamma $(flags)

randompv: $(objects) $(fft_lib) $(sourcedir)/init/randompv.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/randompv.f90 -o randompv $(flags)

adjust: $(objects) $(fft_lib) $(sourcedir)/init/adjust.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/adjust.f90 -o adjust $(flags)

rest: $(objects) $(fft_lib) $(sourcedir)/init/rest.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/rest.f90 -o rest $(flags)

rossby: $(sourcedir)/init/rossby.f90
	$(f90) $(sourcedir)/init/rossby.f90 -o rossby $(flags)

dambreak: $(objects) $(sourcedir)/init/dambreak.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/dambreak.f90 -o dambreak $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


