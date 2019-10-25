#include "protein.h"

class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
	void calchi();
};
void Molecule::calchi(){
	initneib();
	for(int ir=0; ir<nres; ir++){
		cout<<rnam[ir]<<" ";
		int nrot = 6;
		if(rnam[ir] == "PHE") nrot = 2;
		if(rnam[ir] == "TRP") nrot = 2;
		if(rnam[ir] == "TYR") nrot = 2;
		for(int ii=5; ii<rstart[ir+1]-rstart[ir]; ii++){
			if(ii > nrot + 5) continue;
			int i = rstart[ir] + ii;
			initAtomint(i);
			int *iz = atoms[i].getint();
			if(iz[3] != 0) continue;
//			cout<<atoms[i].getname()<<' ';
			double tor = 360.;
			if(atoms[i].filled() && atoms[iz[0]].filled() && atoms[iz[1]].filled()
						&& atoms[iz[2]].filled()){
				tor = torsion(i, iz[0], iz[1], iz[2]);
			}
			cout<<atoms[i].getname()<<' '<<tor<<' ';
		}
		cout<<endl;
	}
}
int main(int argc, char *argv[]){
	for(int i=1; i<argc; i++){
		Molecule *mol=new Molecule(argv[i]);
		cout<<argv[i]<<endl;
		mol -> calchi();
		delete mol;
	}
}
