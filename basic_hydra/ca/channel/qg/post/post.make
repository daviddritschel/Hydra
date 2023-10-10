 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------
 #Rules:
#conprint: $(objects) post/conprint.f90
#	$(f90) parameters.o constants.o post/conprint.f90 -o conprint $(flags)

#genfg: $(objects) post/genfg.f90
#	$(f90) parameters.o constants.o post/genfg.f90 -o genfg $(flags)

pvsep: $(objects) post/pvsep.f90
	$(f90) parameters.o constants.o post/pvsep.f90 -o pvsep $(flags)

zonal: $(objects) $(fft_lib) post/zonal.f90
	$(f90) $(fft_lib) $(objects) post/zonal.f90 -o zonal $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
#post_all: $(present_post_files)
post_all: zonal

