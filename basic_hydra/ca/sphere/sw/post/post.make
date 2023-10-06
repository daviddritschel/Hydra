 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
dgbal: $(objects) $(fft_lib) $(sourcedir)/post/dgbal.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/dgbal.f90 -o dgbal $(flags)
ispectra: $(objects) $(fft_lib) $(sourcedir)/post/ispectra.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/ispectra.f90 -o ispectra $(flags)
profile: $(objects) $(fft_lib) $(sourcedir)/post/profile.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/profile.f90 -o profile $(flags)
zonal: $(objects) $(fft_lib) $(sourcedir)/post/zonal.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/zonal.f90 -o zonal $(flags)
genfg: $(objects) $(fft_lib) post/genfg.f90
	$(f90) $(fft_lib) $(objects) post/genfg.f90 -o genfg $(flags)
ortho: $(objects) post/ortho.f90
	$(f90) parameters.o constants.o post/ortho.f90 -o ortho $(flags)
measure: $(objects) post/measure.f90
	$(f90) parameters.o constants.o post/measure.f90 -o measure $(flags)
hov: $(objects) post/hov.f90
	$(f90) parameters.o post/hov.f90 -o hov $(flags)
jets: $(objects) post/jets.f90
	$(f90) parameters.o post/jets.f90 -o jets $(flags)
angm: $(objects) post/angm.f90
	$(f90) parameters.o post/angm.f90 -o angm $(flags)
avgdiag: $(objects) post/avgdiag.f90
	$(f90) parameters.o post/avgdiag.f90 -o avgdiag $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

