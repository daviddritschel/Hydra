 # Set up source and install directories,
 # these can be edited if the source/install structure needed is different to the standard.
sourcedir = .
installdir = ..

 # Set main source and object files:
sources = $(sourcedir)/parameters.f90 $(sourcedir)/constants.f90 $(sourcedir)/variables.f90 $(sourcedir)/spectral.f90 $(sourcedir)/contours.f90
objects = parameters.o constants.o variables.o spectral.o contours.o

 # Set location of fft library:
fft_lib = ~/hydra/lib/stafft/stafft.f90 ~/hydra/lib/stafft/deriv1d.f90 

 # Set flag options: 
norm_flags = -O3 -mcmodel=medium #Compiler flags
test_flags = -O3 -mcmodel=medium -fbounds-check -fbacktrace  #Compiler flags

 # Set compiler and flag used:
f90=gfortran
flags=$(norm_flags)

#----------------------------------------------------------------------
 # Rules:
all: casl setup proxy_init_all proxy_post_all install

casl: $(objects) $(fft_lib) $(sourcedir)/casl.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/casl.f90 -o casl $(flags)

grid: $(objects) $(sourcedir)/grid.f90
	$(f90) parameters.o constants.o variables.o contours.o $(sourcedir)/grid.f90 -o grid $(flags)

enstro: $(objects) $(sourcedir)/enstro.f90
	$(f90) parameters.o $(sourcedir)/enstro.f90 -o enstro $(flags)

area: $(objects) $(sourcedir)/area.f90
	$(f90) parameters.o $(sourcedir)/area.f90 -o area $(flags)

$(objects): $(fft_lib) $(sources) 
	$(f90) $(fft_lib) $(sources) -c $(flags)

 # Phony rules for setup etc.: 
.PHONY: all clean setup install
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
