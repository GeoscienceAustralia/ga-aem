
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

includes   = -I$(srcdir) -I$(cpputilssrc)
executable = $(exedir)/ctlinedata2sgrid.exe

all: compile link
allclean: clean compile link

objects += $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
objects += $(cpputilssrc)/blocklanguage.o
objects += $(cpputilssrc)/geometry3d.o
objects += $(srcdir)/ctlinedata2sgrid.o

%.o : %.cpp
	@echo ' '
	@echo 'Compiling ' $<
	$(cxx) -c $(includes) $(cxxflags) $< -o $@

compile: $(objects)

link: $(objects)
	mkdir -p $(exedir)
	@echo ' '
	@echo Linking
	$(cxx) $(objects) -o $(executable)

clean: 
	@echo ' '
	@echo Cleaning
	rm -f $(objects)
	rm -f $(executable)


