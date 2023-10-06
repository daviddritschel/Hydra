 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------
 #Rules:
conprint: $(objects) $(fft_lib) post/conprint.f90
	$(f90) $(fft_lib) $(objects) post/conprint.f90 -o conprint $(flags)
decomp: $(objects) $(fft_lib) post/decomp.f90
	$(f90) $(fft_lib) $(objects) post/decomp.f90 -o decomp $(flags)
energy: $(objects) $(fft_lib) post/energy.f90
	$(f90) $(fft_lib) $(objects) post/energy.f90 -o energy $(flags)
genfg: $(objects) $(fft_lib) post/genfg.f90
	$(f90) $(fft_lib) $(objects) post/genfg.f90 -o genfg $(flags)
modes: $(objects) $(fft_lib) post/modes.f90
	$(f90) $(fft_lib) $(objects) post/modes.f90 -o modes $(flags)
vormod: $(objects) $(fft_lib) post/vormod.f90
	$(f90) $(fft_lib) $(objects) post/vormod.f90 -o vormod $(flags)
zonal: $(objects) $(fft_lib) post/zonal.f90
	$(f90) $(fft_lib) $(objects) post/zonal.f90 -o zonal $(flags)
zspec: $(objects) $(fft_lib) post/zspec.f90
	$(f90) $(fft_lib) $(objects) post/zspec.f90 -o zspec $(flags)
displace: $(objects) post/displace.f90
	$(f90) parameters.o constants.o post/displace.f90 -o displace $(flags)
stair: $(objects) post/stair.f90
	$(f90) parameters.o constants.o post/stair.f90 -o stair $(flags)
slope: $(objects) post/slope.f90
	$(f90) parameters.o constants.o post/slope.f90 -o slope $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

