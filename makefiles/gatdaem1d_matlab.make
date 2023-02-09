
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

cxxflags   += -fPIC
ldflags    = -shared
bindir     = ../matlab/bin/linux
includes   = -I$(srcdir) -I$(cpputilssrc)
libs       = -L$(FFTW_DIR) -lfftw3
library    = $(bindir)/gatdaem1d.so

all: compile link
allclean: clean compile link

objects += $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
objects += $(srcdir)/gatdaem1d.o

%.o : %.cpp
	@echo ' '
	@echo 'Compiling ' $<
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

