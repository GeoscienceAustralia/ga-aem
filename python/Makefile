SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

cxxflags   += -fPIC -O3
ldflags    += -shared
bindir     = gatdaem1d

srcdir     = ../src
objdir     = ./obj
includes   = -I$(srcdir) -I$(FFTW_ROOT)/include -I../submodules/cpp-utils/include
libs       = -L$(FFTW_ROOT)/lib -lfftw3
library    = $(bindir)/gatdaem1d.so

all: compile link
allclean: clean compile link

objects = \
	$(objdir)/gatdaem1d.o \

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