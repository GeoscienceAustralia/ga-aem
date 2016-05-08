
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

#GNU compiler on raijin.nci.org.au
cxx        = g++
cxxflags   = -std=c++11 -O3 -Wall -fPIC
libs       = -L$(FFTW_DIR) -lfftw3
ldflags    = -shared
bindir     = ../matlab/bin/raijin-gnu

#Intel compiler on raijin.nci.org.au
#cxx        = icpc
#cxxflags   = -std=c++11 -O3 -Wall -fPIC -diag-disable remark
#libs       = -L$(FFTW_DIR) -lfftw3
ldflags     = -shared
#bindir     = ../matlab/bin/raijin-intel

srcdir     = ../src
tntdir     = ../third_party/tnt
objdir     = ./obj
includes   = -I$(srcdir) -I$(tntdir)
library    = $(bindir)/gatdaem1d.so

all: compile link
allclean: clean compile link

objects = \
	$(objdir)/blocklanguage.o \
	$(objdir)/geometry3d.o \
	$(objdir)/le.o \
	$(objdir)/general_utils.o \
	$(objdir)/file_utils.o \
	$(objdir)/tdemsystem.o \
	$(objdir)/gatdaem1d.o 

$(objects): $(objdir)/%.o: $(srcdir)/%.cpp
	mkdir -p $(objdir)
	@echo ' '
	@echo Compiling $<
	$(cxx) -c $(includes) $(cxxflags) $< -o $@

compile: $(objects)

link: $(objects)
	mkdir -p $(bindir)
	@echo ' '
	@echo Creating library
	$(cxx) $(ldflags) $(objects) $(libs) -o $(library)

clean: 
	@echo ' '
	@echo Cleaning
	rm -f $(objects)
	rm -f $(library)
