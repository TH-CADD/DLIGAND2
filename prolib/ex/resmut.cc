#include "protein.h"

class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
};
int main(int argc, char *argv[]){
	if(argc != 2) die("usage: %s pdb\n", argv[0]);
	Molecule *mol=new Molecule(argv[1]);
	mol->resMutation(0, "GLU");
	mol->wrpdb("1.pdb");
	mol->mutRotamer(0, 1);
	mol->wrpdb("2.pdb");
	delete mol;
}
