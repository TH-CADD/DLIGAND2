#include "protein.h"

class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
};
int main(int argc, char *argv[]){
	if(argc != 3) die("usage: %s pdb out\n", argv[0]);
	Molecule *mol=new Molecule(argv[1]);
	mol->refill();
	mol->mutRotamer(0, 0);
	mol->wrpdb(argv[2]);
	delete mol;
}
