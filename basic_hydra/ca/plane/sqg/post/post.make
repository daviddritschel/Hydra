 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
genfg: $(objects) $(sourcedir)/congen.f90 $(sourcedir)/post/genfg.f90
	$(f90) parameters.o constants.o generic.o contours.o $(sourcedir)/congen.f90 $(sourcedir)/post/genfg.f90 -o genfg $(flags) 

energy: $(objects) $(fft_lib) $(sourcedir)/post/energy.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/energy.f90 -o energy $(flags)

zonal: $(objects) $(fft_lib) $(sourcedir)/post/zonal.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/zonal.f90 -o zonal $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

