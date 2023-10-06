 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
#dgbal: $(objects) $(fft_lib) $(sourcedir)/post/dgbal.f90
#	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/dgbal.f90 -o dgbal $(flags)

#genfg: $(objects) $(fft_lib) $(sourcedir)/congen.f90 $(sourcedir)/post/genfg.f90
#	$(f90) $(fft_lib) $(objects) $(sourcedir)/congen.f90 $(sourcedir)/post/genfg.f90 -o genfg $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

