 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
pressure: $(objects) $(fft_lib) $(sourcedir)/post/pressure.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/pressure.f90 -o pressure $(flags)

gamma-tilde: $(objects) $(fft_lib) $(sourcedir)/post/gamma-tilde.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/gamma-tilde.f90 -o gamma-tilde $(flags)

fopvbal: $(objects) $(fft_lib) $(sourcedir)/post/fopvbal.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/fopvbal.f90 -o fopvbal $(flags)

dgbal: $(objects) $(fft_lib) $(sourcedir)/post/dgbal.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/dgbal.f90 -o dgbal $(flags)

ispectra: $(objects) $(fft_lib) $(sourcedir)/post/ispectra.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/ispectra.f90 -o ispectra $(flags)

genfg: $(objects) $(fft_lib) $(sourcedir)/congen.f90 $(sourcedir)/post/genfg.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/congen.f90 $(sourcedir)/post/genfg.f90 -o genfg $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

