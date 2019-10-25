#include "protein.h"
#include "aastd.h"

class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
	void addCB();
	void calHSEa();
};
void Molecule::addCB(){
	initneib();
	for(int ir=0; ir<nres; ir++){
		int ib = rstart[ir]+4;
		if(rnam[ir] == "GLY") continue;
		initAtomint(ib); initAtomzcrd0(ib);
		if(! atoms.at(ib).filled()) makexyz1(ib); //atoms.at(ib).xyzatm(atoms);
	}
	for(int ir=0; ir<nres; ir++){
		if(rnam[ir] != "GLY") continue;
		addCB_GLY(ir);
		int ib = rstart[ir]+4;
		if(! atoms.at(ib).filled()) makexyz1(ib); //atoms.at(ib).xyzatm(atoms);
	}
}
double cutoff = 13.;
bool balpha = true;
int main(int argc, char *argv[]){
	if(argc < 2){
		fprintf(stderr, "Usage: %s pdb\n", argv[0]);
		fprintf(stderr, "Attn: set env DATADIR=[/data1/yueyang/source/lib/]\n");
		exit(0);
	}
	Molecule *mol=new Molecule(argv[1]);
	mol -> addCB();
	mol -> wrpdb(stdout);
	delete mol;
}
