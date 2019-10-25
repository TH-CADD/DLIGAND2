#include "protein.h"

class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
	void calzcrd();
};
void Molecule::calzcrd(){
	printf("# rn Phi Psi Omega\n");
	for(int ir=0; ir<nres; ir++){
		double ag[3];
		for(int m=0; m<3; m++) ag[m] = 360;
		int i1 = rstart[ir];
		if(ir > 0) {
			int i0 = rstart[ir-1];
			ag[0] = torsion(i0+2, i1, i1+1, i1+2);
			ag[2] = torsion(i0+1, i0+2, i1, i1+1);
		}
		if(ir < nres-1) {
			int i2 = rstart[ir+1];
			ag[1] = torsion(i1, i1+1, i1+2, i2);
		}
		int id = resid[ir];
		printf("%d %c %.2f %.2f %.2f\n", seq0[ir], rnam1_std[id], ag[0], ag[1], ag[2]);
	}
}
int main(int argc, char *argv[]){
	for(int i=1; i<argc; i++){
		Molecule *mol=new Molecule(argv[i]);
		printf("# %s\n", argv[i]);
		mol -> calzcrd();
		delete mol;
	}
}
