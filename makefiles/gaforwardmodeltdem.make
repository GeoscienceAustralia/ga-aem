
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

includes   = -I$(srcdir) -I$(cpputilssrc) -I$(tntdir)
libs       = -L$(FFTW_DIR) -lfftw3
executable = $(exedir)/gaforwardmodeltdem.exe

all: compile link
allclean: clean compile link

objects += $(cpputilssrc)/blocklanguage.o
objects += $(cpputilssrc)/geometry3d.o
objects += $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
objects += $(srcdir)/tdemsystem.o
objects += $(srcdir)/gaforwardmodeltdem.o

%.o : %.cpp
	@echo ' '
	@echo 'Compiling ' $<
	$(cxx) -c $(includes) $(cxxflags) $< -o $@

compile: $(objects)

link: $(objects)
	mkdir -p $(exedir)
	@echo ' '
	@echo Linking
	$(cxx) $(objects) $(libs) -o $(executable)

clean: 
	@echo ' '
	@echo Cleaning
	rm -f $(objects)
	rm -f $(executable)


