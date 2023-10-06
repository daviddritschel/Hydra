 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
r4toc2: $(fft_lib) $(objects) $(sourcedir)/post/r4toc2.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/r4toc2.f90 -o r4toc2 $(flags) 

genfg: $(objects) $(fft_lib) $(sourcedir)/post/genfg.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/genfg.f90 -o genfg $(flags)

conprint: $(objects) $(fft_lib) $(sourcedir)/post/conprint.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/conprint.f90 -o conprint $(flags)

split: $(fft_lib) $(objects) $(sourcedir)/post/split.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/split.f90 -o split $(flags)

cofm: $(fft_lib) $(objects) $(sourcedir)/post/cofm.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/cofm.f90 -o cofm $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

