 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
ortho: $(objects) post/ortho.f90
	$(f90) parameters.o constants.o post/ortho.f90 -o ortho $(flags)
measure: $(objects) post/measure.f90
	$(f90) parameters.o constants.o post/measure.f90 -o measure $(flags)
enstrophy: $(objects) post/enstrophy.f90
	$(f90) parameters.o constants.o post/enstrophy.f90 -o enstrophy $(flags)
genfg: $(objects) $(fft_lib) post/genfg.f90
	$(f90) $(fft_lib) $(objects) post/genfg.f90 -o genfg $(flags)
zonal: $(objects) $(fft_lib) post/zonal.f90
	$(f90) $(fft_lib) $(objects) post/zonal.f90 -o zonal $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

