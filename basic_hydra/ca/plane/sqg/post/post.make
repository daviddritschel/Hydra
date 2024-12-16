 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
conprint: $(objects) $(sourcedir)/post/conprint.f90
	$(f90) parameters.o constants.o variables.o contours.o $(sourcedir)/post/conprint.f90 -o conprint $(flags)

area: $(objects) $(sourcedir)/post/area.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/area.f90 -o area $(flags)

genfg: $(objects) $(sourcedir)/congen.f90 $(sourcedir)/post/genfg.f90
	$(f90) parameters.o constants.o contours.o $(sourcedir)/congen.f90 $(sourcedir)/post/genfg.f90 -o genfg $(flags) 

froude: $(objects) $(fft_lib) $(sourcedir)/post/froude.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/froude.f90 -o froude $(flags)

cross: $(objects) $(fft_lib) $(sourcedir)/post/cross.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/cross.f90 -o cross $(flags)

new_cross: $(objects) $(fft_lib) $(sourcedir)/post/new_cross.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/new_cross.f90 -o new_cross $(flags)

zvel: $(objects) $(fft_lib) $(sourcedir)/post/zvel.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/zvel.f90 -o zvel $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

