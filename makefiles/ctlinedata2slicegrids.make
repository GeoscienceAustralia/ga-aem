
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

includes   = -I$(srcdir) -I$(cpputilssrc) -I$(csv_include) -I$(gdal_include)
libs       = -lgdal
executable = $(exedir)/ctlinedata2slicegrids.exe

all: compile link
allclean: clean compile link

objects += $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
objects += $(srcdir)/ctlinedata2slicegrids.o

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


