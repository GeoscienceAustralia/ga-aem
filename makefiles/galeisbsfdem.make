
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

executable  = $(exedir)/galeisbsfdem.exe
includes    = -I$(srcdir)
includes   += -I$(cpputilssrc)
libs        = -fopenmp
cxxflags    = -fopenmp
cxxflags   += -D_MPI_ENABLED
#cxxflags   += -DUSEGLOBALSTACKTRACE


all: compile link
allclean: clean compile link

objects += $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
objects += $(srcdir)/galeisbsfdem.o

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

