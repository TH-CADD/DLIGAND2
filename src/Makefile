include Makefile.in.intel

SRC = dfire_misc.cc molecule.cc
#
all:  score2.$(TYPE)  #stats.$(TYPE)
gdb: all
gnu: all
#
stats.$(TYPE): $(OBJ) .$(TYPE).stats.o
	$(CPP) $(CPPFLAGS) $(OPTFLAGS) $^ -o $@ -lpro_gnu
score.$(TYPE): $(OBJ) .$(TYPE).score.o
	$(CPP) $(CPPFLAGS) $(OPTFLAGS) $^ -o $@ -lpro_gnu
score2.$(TYPE): $(OBJ) .$(TYPE).score2.o
	$(CPP) $(CPPFLAGS) $(OPTFLAGS) $^ -o $@ -lpro_gnu
clean:
	rm -f .*.o score2.intel
