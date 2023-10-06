 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
vertical-shear: $(objects) $(sourcedir)/init/vertical-shear.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/vertical-shear.f90 -o vertical-shear $(flags)
eq-displace: $(objects) $(fft_lib) $(sourcedir)/init/eq-displace.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/eq-displace.f90 -o eq-displace $(flags)

jet: $(objects) $(fft_lib) $(sourcedir)/init/jet.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/jet.f90 -o jet $(flags)
rest: $(objects) $(sourcedir)/init/rest.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/rest.f90 -o rest $(flags)

ellipse_topo: $(objects) $(sourcedir)/init/ellipse_topo.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/ellipse_topo.f90 -o ellipse_topo $(flags)
sine_topo: $(objects) $(sourcedir)/init/sine_topo.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/sine_topo.f90 -o sine_topo $(flags)
random_topo: $(objects) $(fft_lib) $(sourcedir)/init/random_topo.f90
	$(f90) $(fft_lib) parameters.o constants.o $(sourcedir)/init/random_topo.f90 -o random_topo $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


