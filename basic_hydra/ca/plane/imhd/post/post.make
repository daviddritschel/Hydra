 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
conprint: $(objects) $(sourcedir)/post/conprint.f90
	$(f90) parameters.o constants.o contours.o $(sourcedir)/post/conprint.f90 -o conprint $(flags)

tracercont: $(objects) $(sourcedir)/post/tracercont.f90
	$(f90) parameters.o constants.o contours.o $(sourcedir)/post/tracercont.f90 -o tracercont $(flags)

source: $(objects) $(fft_lib) $(sourcedir)/post/source.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/source.f90 -o source $(flags)

polar: $(objects) $(fft_lib) $(sourcedir)/post/polar.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/polar.f90 -o polar $(flags)

energy: $(objects) $(fft_lib) $(sourcedir)/post/energy.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/energy.f90 -o energy $(flags)

potentials: $(objects) $(fft_lib) $(sourcedir)/post/potentials.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/potentials.f90 -o potentials $(flags)

genfg: $(objects) $(fft_lib) $(sourcedir)/congen.f90 $(sourcedir)/post/genfg.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/congen.f90 $(sourcedir)/post/genfg.f90 -o genfg $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

