#include "protein.h"

class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
	void caltheta();
};
void Molecule::caltheta(){
	bool bgap[nres];
	for(int i=0; i<nres; i++) bgap[i] = 1;
	for(int i=0; i<nres-1; i++){
		int i1=rstart[i]+1, i2=rstart[i+1]+1;
		if(! (isfilled(i1) && isfilled(i2)) ) continue;
		double r2 = distance2(i1, i2);
		double r = sqrt(r2);
		if(fabs(r-3.8) < 0.5) bgap[i] = 0;
//		else fprintf(stderr, "%d %f\n", i, r);
	}
//
	double theta[nres], tors[nres];
	for(int i=0; i<nres; i++) theta[i] = tors[i] = 360;
	for(int i=0; i<nres-2; i++){
		int i1=rstart[i]+1;
		int i2=rstart[i+1]+1;
		int i3=rstart[i+2]+1;
		theta[i+1] = angle(i1, i2, i3);
		if(i >= nres-3) continue;
		int i4=rstart[i+3]+1;
		tors[i+2] = torsion(i1, i2, i3, i4);
	}
	for(int i=0; i<nres; i++){
		int id = resid[i];
		int c = 4;
		if(i<nres-3) c = bgap[i] + bgap[i+1] + bgap[i+2] + bgap[i+3];
		printf("%d %c %.1f %.1f %d\n", seq0[i], rnam1_std[id], theta[i], tors[i], c);
	}
}
int main(int argc, char *argv[]){
	for(int i=1; i<argc; i++){
		Molecule *mol=new Molecule(argv[i]);
		printf("#%s\n", argv[i]);
		mol -> caltheta();
		delete mol;
	}
}
