#include "protein.h"

class Molecule:public Protein{
public:
	Molecule(const string fn):Protein(fn){};
	double gyrate(void);
};
double Molecule::gyrate(void){
	double x0[3], s1[3], s2[3], rg=0;
	int nr=0;
	for(int ir=0; ir<nres; ir++){
		int i0 = rstart[ir];
		if(! atoms.at(i0+1).filled()) continue;
		double *xp = atoms.at(i0+1).getx();
		if(nr == 0) {
			for(int m=0; m<3; m++) x0[m] = xp[m];
		}
		for(int m=0; m<3; m++){
			double dt = xp[m] - x0[m];
			s1[m] += dt; s2[m] += dt*dt;
		}
		nr ++;
	}
	for(int m=0; m<3; m++){
		s1[m] /= nr; s2[m] /= nr;
		rg += s2[m] - s1[m]*s1[m];
	}
	rg = sqrt(rg);
	double rg0 = 3.*pow(double(nr), 1./3);   // 3*rg0 + 2 in ROSETTA
	printf("%s %d %.1f # %d %.2f\n", pdbnm.c_str(), nres, rg, nr, rg/rg0);
	return rg;
}
int main(int argc, char *argv[]){
	for(int i=1; i<argc; i++){
		Molecule *mol = new Molecule(argv[i]);
		mol -> gyrate();
		delete mol;
	}
}
