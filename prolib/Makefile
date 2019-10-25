include Makefile.in

# quatfit* not necessary
SRC = restype.cc protein.cc xyzatm.cc restpl.cc rotamer.cc pbb_std.cc obsolete.cc protein2.cc protein3.cc quatfit.cc quatfit2.cc
#
all: libpro_$(TYPE).a 
all: main.$(TYPE)
gdb: all
gnu: all
pg: all
mp: all
intel: all
CPPFLAGS2 =
#
libpro_$(TYPE).a: $(OBJ)
	ar cr $@ $^
$(OBJ): *.h

%.$(TYPE): .$(TYPE).%.o libpro_$(TYPE).a
	$(CPP) $(CPPFLAGS) $< -o $@ -L. -lpro_$(TYPE)

clean:
	rm -f .*.o
