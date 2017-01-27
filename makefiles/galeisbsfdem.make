
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

cxxflags   += -D_MPI_ENABLED -fopenmp
includes   = -I$(srcdir) -I$(cpputilssrc) -I$(tntdir)
libs       = -fopenmp
executable = $(exedir)/galeisbsfdem.exe

all: compile link
allclean: clean compile link

objects += $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
objects += $(cpputilssrc)/blocklanguage.o
objects += $(cpputilssrc)/geometry3d.o
objects += $(cpputilssrc)/fielddefinition.o
objects += $(cpputilssrc)/matrix_ops.o
objects += $(cpputilssrc)/rollpitchyaw.o
objects += $(srcdir)/layeredearthmodeller.o
objects += $(srcdir)/fdemsystem.o
objects += $(srcdir)/galeisbsfdem.o

%.o : %.cpp
	@echo ' '
	@echo 'Compiling ' $<
	$(cxx) -c $(includes) $(cxxflags) $< -o $@

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

