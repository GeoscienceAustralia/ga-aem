
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

executable  = $(exedir)/garjmcmctdem.exe
testexe     = $(exedir)/runtests_garjmcmc.exe

includes    = -I$(srcdir)
includes   += -I$(cpputilssrc)
libs        = $(if $(FFTW_DIR),-L$(FFTW_DIR) -lfftw3,-lfftw3)
libs       += $(ldflags)
testlibs    = $(libs) -lgtest -lgmock -lgtest_main -lpthread

cxxflags   += -D_MPI_ENABLED
#cxxflags   += -DUSEGLOBALSTACKTRACE

ifeq ($(HAVE_NETCDF),1)
    cxxflags += -DHAVE_NETCDF
    includes +=  $(if $(geophysics_netcdf_root),-I$(geophysics_netcdf_root)/src,)
    includes +=  $(if $(geophysics_netcdf_root),-I$(geophysics_netcdf_root)/submodules/marray/include/andres,)
    libs     +=  -lnetcdf -lnetcdf_c++4
endif

ifeq ($(HAVE_GDAL),1)
    cxxflags += -DHAVE_GDAL
    libs     += -lgdal
    objects  += $(cpputilssrc)/gdal_utils.o
endif

ifeq ($(HAVE_CGAL),1)
    cxxflags += -DHAVE_CGAL
    libs     += -lCGAL_Core
    objects  += $(cpputilssrc)/cgal_utils.o
endif


all: compile link
allclean: clean compile link

test: tcompile tlink
testclean: clean tcompile tlink

objects += $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
objects += $(srcdir)/garjmcmctdem.o

testobjects += $(tstdir)/test_rjmcmc1d.o
testobjects += $(tstdir)/test_rjmcmc1dtdeminverter.o
testobjects += $(tstdir)/test_ptrvec.o
testobjects += $(cpputilssrc)/general_utils.o
testobjects += $(cpputilssrc)/file_utils.o

%.o : %.cpp
	@echo ' '
	@echo 'Compiling ' $<
	$(mpicxx) -c $(includes) $(cxxflags) $< -o $@

compile: $(objects)

tcompile: $(testobjects)

link: $(objects)
	mkdir -p $(exedir)
	@echo ' '
	@echo Linking
	$(mpicxx) $(objects) $(libs) -o $(executable)

tlink: $(testobjects)
	mkdir -p $(exedir)
	@echo ' '
	@echo 'Linking test executable'
	$(mpicxx) $(testobjects) $(testlibs) -o $(testexe)

clean:
	@echo ' '
	@echo Cleaning
	rm -f $(objects)
	rm -f $(testobjects)
	rm -f $(executable)
	rm -f $(testexe)
