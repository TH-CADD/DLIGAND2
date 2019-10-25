#include "protein.h"

class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
	void caldis_cb();
};
void Molecule::caldis_cb(){
	int icb[nres];
	for(int ir=0; ir<nres; ir++){
		icb[ir] = rstart[ir] + 4;
		if(rnam[ir]=="GLY") icb[ir] = rstart[ir] + 1;
	}
	for(int ir=0; ir<nres; ir++)
	for(int jr=ir+3; jr<nres; jr++){
		double r2 = distance2(icb[ir], icb[jr]);
//		if(r2 < 9.5*9.5)
		printf("%d %d %c%c %f %d %d\n", ir+1, jr+1, rnam1_std[resid[ir]], rnam1_std[resid[jr]], sqrt(r2), seq0[ir], seq0[jr]);
	}
}
int main(int argc, char *argv[]){
	for(int i=1; i<argc; i++){
		Molecule *mol=new Molecule(argv[i]);
		mol->caldis_cb();
		delete mol;
	}
}
