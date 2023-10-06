 # Set up source and install directories,
 # these can be edited if the source/install structure needed is different to the standard.
sourcedir = .
installdir = ..

 # Set main source and object files:
sources = $(sourcedir)/parameters.f90 $(sourcedir)/constants.f90 $(sourcedir)/spectral.f90 \
          $(sourcedir)/generic.f90 $(sourcedir)/contours.f90 $(sourcedir)/common.f90 $(sourcedir)/congen.f90 $(sourcedir)/evolution.f90 #Main f90 sources
objects = parameters.o constants.o spectral.o \
          generic.o contours.o common.o congen.o evolution.o #Main .o files 

 # Set location of fft library:
fft_lib = ~/hydra/lib/stafft/stafft.f90 ~/hydra/lib/stafft/sta2dfft.f90 

 # Set flag options: 
norm_flags = -O3 -mcmodel=medium #Compiler flags
test_flags = -O3 -mcmodel=medium -fbounds-check -fbacktrace  #Compiler flags

 # Set compiler and flag used:
f90=gfortran
flags=$(norm_flags)

#----------------------------------------------------------------------
 # Rules:
all: caps variance setup proxy_init_all install

caps: $(objects) $(fft_lib) caps.f90
	$(f90) $(fft_lib) $(objects) caps.f90 -o caps $(flags)

variance: $(objects) $(sourcedir)/variance.f90
	$(f90) parameters.o constants.o $(sourcedir)/variance.f90 -o variance $(flags)

$(objects): $(fft_lib) $(sources) 
	$(f90) $(fft_lib) $(sources) -c $(flags)

.PHONY: clean setup install
clean:
	rm -f *.o *.mod
setup:
	test -d $(installdir)/cont || mkdir $(installdir)/cont
	test -d $(installdir)/fine || mkdir $(installdir)/fine
	test -d $(installdir)/evolution || mkdir $(installdir)/evolution
	test -d $(installdir)/spectra || mkdir $(installdir)/spectra

install:
	find -type f -perm /a+x -exec /bin/mv {} $(installdir)/ \;

 # Include intialisation rules if these exist:
-include init/init.make

 # Proxy rules to ensure whatever exists in the sub-directories
 # also gets built:
ifdef init_exists
proxy_init_all: init_all
else 
proxy_init_all:
	 @echo No init directory found 
endif

