sources = parameters.f90 constants.f90 variables.f90 spectral.f90 \
          generic.f90 contours.f90 common.f90 congen.f90 evolution.f90 #Main f90 sources

objects = parameters.o constants.o variables.o spectral.o \
          generic.o contours.o common.o congen.o evolution.o #Main .o files 

flags = -O3 -mcmodel=medium #Compiler flags
# *** Uncomment below and comment above if testing ***
#flags = -O3 -mcmodel=medium -fbounds-check -fbacktrace  #Compiler flags

all: casl fgrid image dam-break slug split setup

casl: $(objects) ~/Libdir/f02agf.f ~/Libdir/double_fft.f casl.f90
	gfortran ~/Libdir/f02agf.f ~/Libdir/double_fft.f $(objects) \
                 casl.f90 -o casl $(flags)

image: parameters.f90 constants.f90 image.f90
	gfortran parameters.f90 constants.f90 image.f90 -o image $(flags) 

fgrid: $(objects) ~/Libdir/f02agf.f ~/Libdir/double_fft.f fgrid.f90
	gfortran ~/Libdir/f02agf.f ~/Libdir/double_fft.f $(objects) \
                 fgrid.f90 -o fgrid $(flags)

dam-break: parameters.o constants.o dam-break.f90
	gfortran parameters.o constants.o dam-break.f90 -o dam-break $(flags)

slug: parameters.o constants.o slug.f90
	gfortran parameters.o constants.o slug.f90 -o slug $(flags)

split: $(objects) ~/Libdir/f02agf.f ~/Libdir/double_fft.f split.f90
	gfortran ~/Libdir/f02agf.f ~/Libdir/double_fft.f $(objects) \
                 split.f90 -o split $(flags)

extract: parameters.o constants.o extract.f90
	gfortran parameters.o constants.o extract.f90 -o extract $(flags)

$(objects): $(sources) 
	gfortran $(sources) -c $(flags)

.PHONY: clean setup
clean:
	rm -f *.o *.mod
setup:
	test -d contours || mkdir contours

