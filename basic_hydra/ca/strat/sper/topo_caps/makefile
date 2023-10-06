 # Set up source and install directories,
 # these can be edited if the source/install structure needed is different to the standard.
sourcedir = .
installdir = ..

 # Set main source and object files:
sources = $(sourcedir)/parameters.f90 $(sourcedir)/constants.f90 $(sourcedir)/variables.f90 $(sourcedir)/spectral.f90 \
          $(sourcedir)/contours.f90 $(sourcedir)/common.f90 $(sourcedir)/congen.f90 $(sourcedir)/evolution.f90 #Main f90 sources
objects = parameters.o constants.o variables.o spectral.o \
          contours.o common.o congen.o evolution.o #Main .o files 

 # Set location of fft library:
fft_lib = ~/hydra/lib/stafft/stafft.f90 ~/hydra/lib/stafft/sta2dfft.f90 ~/hydra/lib/stafft/deriv1d.f90 

 # Set flag options: 
norm_flags = -O3 -mcmodel=medium #Compiler flags
test_flags = -O3 -mcmodel=medium -fbounds-check -fbacktrace  #Compiler flags

 # Set compiler and flag used:
f90=gfortran
#flags=$(test_flags)
flags=$(norm_flags)

#----------------------------------------------------------------------
 # Rules:
all: topo_caps setup proxy_init_all proxy_post_all install

topo_caps: $(objects) $(fft_lib) topo_caps.f90
	$(f90) $(fft_lib) $(objects) topo_caps.f90 -o topo_caps $(flags)

genXYfg: $(objects) genXYfg.f90
	$(f90) parameters.o constants.o contours.o genXYfg.f90 -o genXYfg $(flags)

mix: $(objects) $(fft_lib) mix.f90
	$(f90) $(fft_lib) $(objects) mix.f90 -o mix $(flags)

$(objects): $(fft_lib) $(sources) 
	$(f90) $(fft_lib) $(sources) -c $(flags)

.PHONY: clean setup install
clean:
	rm -f *.o *.mod
setup:
	test -d $(installdir)/cont || mkdir $(installdir)/cont
	test -d $(installdir)/fine || mkdir $(installdir)/fine

install:
	find -type f -perm /a+x -exec /bin/mv {} $(installdir)/ \;

 # Include intialisation and post-processing rules 
 # if these exist:
-include init/init.make
-include post/post.make

 # Proxy rules to ensure whatever exists in the sub-directories
 # also gets built:
ifdef init_exists
proxy_init_all: init_all
else 
proxy_init_all:
	 @echo No init directory found 
endif

ifdef post_exists
proxy_post_all: post_all
else 
proxy_post_all:
	 @echo No post directory found 
endif

