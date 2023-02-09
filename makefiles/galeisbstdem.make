
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

executable  = $(exedir)/galeisbstdem.exe
includes    = -I$(srcdir)
includes   += -I$(cpputilssrc)
includes   += -I$(csv_include)
includes   += -I$(eigen_include)

libs       = $(if $(FFTW_DIR),-L$(FFTW_DIR) -lfftw3,-lfftw3) 
libs       += -fopenmp
cxxflags   += -fopenmp
cxxflags   += -D_MPI_ENABLED
#cxxflags   += -DUSEGLOBALSTACKTRACE

ifeq ($(HAVE_NETCDF),1)
    cxxflags += -DHAVE_NETCDF
    includes +=  -I$(geophysics_netcdf_include)
    includes +=  -I$(marray_include)
    libs     +=  -lnetcdf -lnetcdf_c++4
endif

ifeq ($(HAVE_GDAL),1)
    cxxflags += -DHAVE_GDAL
    includes += -I$(gdal_include)
    libs     += -lgdal
    objects  += $(cpputilssrc)/gdal_utils.o
endif

all: compile link
allclean: clean compile link

objects += $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
objects += $(srcdir)/galeisbstdem.o

%.o : %.cpp
	@echo ' '
	@echo 'Compiling ' $<
	$(mpicxx) -c $(includes) $(cxxflags) $< -o $@

compile: $(objects)

link: $(objects)
	mkdir -p $(exedir)
	@echo ' '
	@echo Linking
	$(mpicxx) $(objects) $(libs) -o $(executable)

clean:
	@echo ' '
	@echo Cleaning
	rm -f $(objects)
	rm -f $(executable)

