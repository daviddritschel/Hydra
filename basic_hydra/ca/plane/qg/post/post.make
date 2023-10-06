 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
genfg: $(objects) $(sourcedir)/congen.f90 $(sourcedir)/post/genfg.f90
	$(f90) parameters.o constants.o generic.o contours.o $(sourcedir)/congen.f90 $(sourcedir)/post/genfg.f90 -o genfg $(flags) 

decay: $(objects) $(sourcedir)/post/decay.f90
	$(f90) parameters.o constants.o $(sourcedir)/post/decay.f90 -o decay $(flags) 

fspec: $(objects) $(fft_lib) $(sourcedir)/post/fspec.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/fspec.f90 -o fspec $(flags)

psi-q: $(objects) $(fft_lib) $(sourcedir)/post/psi-q.f90
	$(f90) $(fft_lib) parameters.o constants.o $(sourcedir)/post/psi-q.f90 -o psi-q $(flags) 

extract: $(objects) $(fft_lib) $(sourcedir)/post/extract.f90
	$(f90) $(fft_lib) parameters.o constants.o $(sourcedir)/post/extract.f90 -o extract $(flags) 

zonal: $(objects) $(fft_lib) $(sourcedir)/post/zonal.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/zonal.f90 -o zonal $(flags)

average: $(objects) $(fft_lib) $(sourcedir)/post/average.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/average.f90 -o average $(flags)

energy: $(objects) $(fft_lib) $(sourcedir)/post/energy.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/energy.f90 -o energy $(flags)

kinetic: $(objects) $(fft_lib) $(sourcedir)/post/kinetic.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/kinetic.f90 -o kinetic $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

