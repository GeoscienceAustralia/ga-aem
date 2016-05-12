
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

#GNU compiler on raijin.nci.org.au
cxx        = mpiCC
cxxflags   = -D_MPI_ENABLED -std=c++11 -O3 -Wall
libs       = -L$(FFTW_DIR) -lfftw3
exedir     = ../bin/raijin/gnu

#Intel compiler on raijin.nci.org.au
#cxx        = mpiCC
#cxxflags   = -D_MPI_ENABLED -std=c++11 -O3 -Wall -diag-disable remark
#libs       = -L$(FFTW_DIR) -lfftw3
#exedir     = ../bin/raijin/intel

srcdir     = ../src
tntdir     = ../third_party/tnt
objdir     = ./obj
includes   = -I$(srcdir) -I$(tntdir)
executable = $(exedir)/garjmcmctdem.exe

all: compile link
allclean: clean compile link

objects = \
	$(objdir)/blocklanguage.o \
	$(objdir)/fielddefinition.o \
	$(objdir)/geometry3d.o \
	$(objdir)/le.o \
	$(objdir)/matrix_ops.o \
	$(objdir)/general_utils.o \
	$(objdir)/file_utils.o \
	$(objdir)/tdemsystem.o \
	$(objdir)/random.o \
	$(objdir)/rjmcmc1d.o \
	$(objdir)/rjmcmc1dtdeminverter.o \
	$(objdir)/garjmcmctdem.o

$(objects): $(objdir)/%.o: $(srcdir)/%.cpp
	mkdir -p $(objdir)
	@echo ' '
	@echo Compiling $<
	$(cxx) -c $(includes) $(cxxflags) $< -o $@

compile: $(objects)

link: $(objects)
	mkdir -p $(exedir)
	@echo
	@echo Linking
	$(cxx) $(objects) $(libs) -o $(executable)

clean:
	@echo
	@echo Cleaning
	rm -f $(objects)
	rm -f $(executable)

