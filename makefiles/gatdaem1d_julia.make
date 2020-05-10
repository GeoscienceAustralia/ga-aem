
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

cxxflags   += -fPIC
ldflags    += -shared
juliadir   = ../julia

includes   = -I$(juliadir) -I$(srcdir) -I$(cpputilssrc)
libs       = -L$(FFTW_DIR) -lfftw3
library    = $(juliadir)/gatdaem1d_julia.so

all: compile link
allclean: clean compile link

objects += $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
objects += $(juliadir)/gatdaem1d_julia.o

%.o : %.cpp
	@echo ' '
	@echo 'Compiling ' $<
	$(cxx) -c $(includes) $(cxxflags) $< -o $@

compile: $(objects)

link: $(objects)
	mkdir -p $(juliadir)
	@echo ' '
	@echo Creating shared library for Julia
	$(cxx) $(ldflags) $(objects) $(libs) -o $(library)

clean: 
	@echo ' '
	@echo Cleaning
	rm -f $(objects)
	rm -f $(library)


#GNU compiler on raijin.nci.org.au
#export cxxflags='-std=c++11 -O3 -Wall -fdiagnostics-color=always -D_GLIBCXX_USE_CXX11_ABI=1 -Wno-unknown-pragmas'
#module load gcc/5.2.0
#module load fftw3/3.3.7-gcc #Use the same module as Julia just to be sure
#module list

#rm gatdaem1d_julia.so
#Compile as shared lib
#g++ -fPIC -shared -O3 $cxxflags \
#-I$srcdir \
#-I$cpputilssrc \
#-L$FFTW_DIR -lfftw3 \
#$cpputilssrc/general_utils.cpp \
#$cpputilssrc/file_utils.cpp \
#gatdaem1d_julia.cpp \
#-o gatdaem1d_julia.so


