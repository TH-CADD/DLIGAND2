#include "protein.h"
#include "dfire.h"
#include "array1.h"

class Molecule:public Protein{
public:
	Molecule(string fn):Protein(fn){};
	void rdmol2(string fn);
	void stats();
	double score();
};
