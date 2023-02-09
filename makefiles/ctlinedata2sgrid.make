
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

includes   = -I$(srcdir) -I$(cpputilssrc) -I$(ticpp_dir)
cxxflags   += -DTIXML_USE_TICPP
executable = $(exedir)/ctlinedata2sgrid.exe

all: compile link
allclean: clean compile link

objects += $(ticpp_dir)/tinyxmlparser.o
objects += $(ticpp_dir)/tinyxmlerror.o
objects += $(ticpp_dir)/tinyxml.o
objects += $(ticpp_dir)/tinystr.o
objects += $(ticpp_dir)/ticpp.o
objects += $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
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


