####
TYPE = intel
CC = icpc
CPP = icpc
F90 = ifort
F77 = ifort
OPTFLAGS = -O2 -I../prolib -L../prolib
CFLAGS =
CPPFLAGS =
F90FLAGS =
FFLAGS =
####
ifeq ($(MAKECMDGOALS),gdb)
	TYPE = gdb
	CC = gcc
	CPP = g++
	F77 = g77
	F90 = gfortran
	CFLAGS =
	CPPFLAGS =
	FFLAGS =
	F90FLAGS =
	OPTFLAGS = -g -p
endif
####
ifeq ($(MAKECMDGOALS),gnu)
	TYPE = gnu
	CC = gcc
	CPP = g++
	F77 = g77
	F90 = gfortran
	CFLAGS =
	CPPFLAGS =
	FFLAGS =
	F90FLAGS =
	OPTFLAGS = -O2
endif
####
OBJ_1 = $(patsubst %.cc, .$(TYPE).%.o, $(SRC))
OBJ_2 = $(patsubst %.c, .$(TYPE).%.o, $(OBJ_1))
OBJ_3 = $(patsubst %.f, .$(TYPE).%.o, $(OBJ_2))
OBJ = $(patsubst %.f90, .$(TYPE).%.o, $(OBJ_3))
#
.$(TYPE).%.o: %.c
	$(CC) $(CFLAGS) $(OPTFLAGS) -c $< -o $@
.$(TYPE).%.o: %.cc
	$(CPP) $(CPPFLAGS) $(OPTFLAGS) -c $< -o $@
.$(TYPE).%.o: %.f
	$(F90) $(FFLAGS) $(OPTFLAGS) -c $< -o $@
.$(TYPE).%.o: %.f90
	$(F90) $(F90FLAGS) $(OPTFLAGS) -c $< -o $@
