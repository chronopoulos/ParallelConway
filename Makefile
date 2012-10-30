# this Makefile is for monolith

#Compilers
MPIF90C=mpif90
MPICC=mpicc
GFC=/usr/bin/gfortran

# Flags
FFLAGS = -ffree-line-length-none

# Object Files
#F_FILES = pgm_mod.f90 read_pgm.f90
#C_FILES = util.c
O_FILES = pgm.o life.o util.o

# PreRequisites

%.o: %.f90
	$(MPIF90C) $(FFLAGS) -c $< -o $@

%.o: %.c
	$(MPICC) $(CFLAGS) -c $< -o $@

# Main Targets

all: clean $(O_FILES)
	$(MPIF90C) -o life $(O_FILES) $(LDFLAGS)

# Housekeeping

.PHONY: clean
clean:		
	/bin/rm -f core life $(O_FILES) *.mod read_pgm *~
