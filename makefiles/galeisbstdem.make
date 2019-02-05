
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

cxxflags   += -D_MPI_ENABLED -fopenmp
includes   = -I/usr/include/openmpi-x86_64 -I$(srcdir) -I$(cpputilssrc) -I$(tntdir)
libs       = -L$(FFTW_DIR) -lfftw3 -fopenmp
executable = $(exedir)/galeisbstdem.exe

all: compile link
allclean: clean compile link

objects += $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
objects += $(cpputilssrc)/geometry3d.o
objects += $(cpputilssrc)/matrix_ops.o
objects += $(srcdir)/tdemsystem.o
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

