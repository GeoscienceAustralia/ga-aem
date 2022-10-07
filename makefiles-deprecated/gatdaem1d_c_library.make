
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

libdir   = ../lib

includes   = -I$(srcdir) -I$(cpputilssrc)
library    = $(libdir)/gatdaem1d_c_library.a

all: compile bind
allclean: clean compile bind

objects += $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
objects += $(srcdir)/gatdaem1d.o

%.o : %.cpp
	@echo ' '
	@echo 'Compiling ' $<
	$(cxx) -c $(includes) $(cxxflags) $< -o $@

compile: $(objects)

bind: $(objects)
	mkdir -p $(libdir)
	@echo ' '
	@echo Binding library
	ar rvs $(library) $(objects)

clean: 
	@echo ' '
	@echo Cleaning
	rm -f $(objects)
	rm -f $(library)


