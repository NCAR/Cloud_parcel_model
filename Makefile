#see forum https://stackoverflow.com/questions/35234003/how-to-create-a-makefile-for-a-fortran-program-using-modules
# details of makefile tutorial http://makepp.sourceforge.net/1.19/makepp_tutorial.html
# compiler
FC=gfortran-mp-8

# compile flag

# source file and objects
SRCS = $(wildcard *.F90)
OBJS = $(patsubst %.F90,%.o,$(wildcard *.F90 ) )

# Ditto for mods
MODS=$(wildcard mod_*.F90)
MOD_OBJS=$(patsubst %.F90,%.o,$(MODS))

# program name
PROGRAM = cloud_parcel
PRG_OBJ = $(PROGRAM)

#clean the suffixes
.SUFFIXES:

# Set the suffixes we are interested in
.SUFFIXES: .F90 .o

#make without parameters will make first target found.
default: $(MOD_OBJS) $(OBJS) $(PROGRAM) clean
module: $(MOD_OBJS)

# compiler steps for all objs
# compile module objects first
$(MOD_OBJS) : mod_%.o : mod_%.F90 

# compile all objects

$(OBJS) : %.o : %.F90
	$(FC) $(FCFLAGS) -c -o $@ $<;

# Linker
$(PROGRAM) : $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $^ $(LIBS)

debug:
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJS)"
	@echo "MODS = $(MODS)"
	@echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "PRG_OBJ = $(PRG_OBJ)"


clean:
	rm -f *.o *.mod IUGG.* stdout
#.PHONY: debug default clean

$(PRG_OBJ) : $(MOD_OBJS)
#$(OBJS) : $(MOD_OBJS)
