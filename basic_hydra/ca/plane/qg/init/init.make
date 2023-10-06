 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
ranpv: $(objects) $(fft_lib) $(sourcedir)/init/ranpv.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/ranpv.f90 -o ranpv $(flags)

narrow: $(objects) $(fft_lib) $(sourcedir)/init/narrow.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/narrow.f90 -o narrow $(flags)

vstrip: $(objects) $(sourcedir)/init/vstrip.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/vstrip.f90 -o vstrip $(flags)

dipole: $(objects) $(sourcedir)/init/dipole.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/dipole.f90 -o dipole $(flags)

staircase: $(objects) $(sourcedir)/init/staircase.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/staircase.f90 -o staircase $(flags)

ellipse: $(objects) $(sourcedir)/init/ellipse.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/ellipse.f90 -o ellipse $(flags)

rossby: $(objects) $(sourcedir)/init/rossby.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/rossby.f90 -o rossby $(flags)

rest: $(objects) $(sourcedir)/init/rest.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/rest.f90 -o rest $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


