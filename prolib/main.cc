#include "protein.h"

class Molecule:public Protein{
public:
	Molecule():Protein(){};
	Molecule(const string str):Protein(str){};
	void test();
};
void Molecule::test(){
	prtdim(rstart.size(), rstart);
	chnseq.resize(nres, 0);
	initneib();
}
int main(int argc, char *argv[]){
	Molecule *mol=new Molecule();
	initRestypes(getdatadir() + "aminodna.dat");
	for(int i=0; i<4; i++) mol ->addRes("ALA");
	mol -> test();
}
