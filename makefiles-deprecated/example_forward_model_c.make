
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .c .o

includes   = -I$(srcdir) 
libs       = ../lib/gatdaem1d_c_library.a -L$(FFTW_DIR) -lfftw3 -lstdc++ -lm
executable = $(exedir)/example_forward_model_c.exe

all: compile link
allclean: clean compile link

objects += $(srcdir)/example_forward_model_c.o

%.o : %.c
	@echo ' '
	@echo 'Compiling ' $<
	$(cc) -c $(includes) $(ccflags) $< -o $@

compile: $(objects)

link: $(objects)
	mkdir -p $(exedir)
	@echo ' '
	@echo Linking
	$(cc) $(objects) $(libs) -o $(executable)

clean: 
	@echo ' '
	@echo Cleaning
	rm -f $(objects)
	rm -f $(executable)


