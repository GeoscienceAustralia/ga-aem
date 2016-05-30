
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

cxxflags   += -D_MPI_ENABLED
srcdir     = ../src
tntdir     = ../third_party/tnt
objdir     = ./obj
includes   = -I$(srcdir) -I$(tntdir)
libs       = -L$(FFTW_DIR) -lfftw3
executable = $(exedir)/garjmcmctdem.exe

all: compile link
allclean: clean compile link

objects = \
	$(objdir)/blocklanguage.o \
	$(objdir)/fielddefinition.o \
	$(objdir)/geometry3d.o \
	$(objdir)/le.o \
	$(objdir)/matrix_ops.o \
	$(objdir)/general_utils.o \
	$(objdir)/file_utils.o \
	$(objdir)/tdemsystem.o \
	$(objdir)/random.o \
	$(objdir)/rjmcmc1d.o \
	$(objdir)/rjmcmc1dtdeminverter.o \
	$(objdir)/garjmcmctdem.o

$(objects): $(objdir)/%.o: $(srcdir)/%.cpp
	mkdir -p $(objdir)
	@echo ' '
	@echo Compiling $<
	$(mpicxx) -c $(includes) $(cxxflags) $< -o $@

compile: $(objects)

link: $(objects)
	mkdir -p $(exedir)
	@echo
	@echo Linking
	$(mpicxx) $(objects) $(libs) -o $(executable)

clean:
	@echo
	@echo Cleaning
	rm -f $(objects)
	rm -f $(executable)

