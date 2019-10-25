#include "protein.h"
#include "surface.cc"

class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
	void calASA();
};

string ANAM_std = "CNOSPH";
double RADIUS_std[] = {1.91, 1.82, 1.70, 2., 2.1, 0.};
void Molecule::calASA(){
	map<char, double> ANAM2RAD;
	for(int i=0; i<ANAM_std.size(); i++) ANAM2RAD[ANAM_std[i]] = RADIUS_std[i];
	double xa[natm][3], radii[natm], w[natm], area[natm];
	for(int i=0; i<natm; i++){
		double *xp = atoms[i].getx();
		for(int m=0; m<3; m++) xa[i][m] = xp[m];
		w[i] = 1.;
		char c1 = atoms[i].getname()[0];
		radii[i] = ANAM2RAD[c1];
	}
	surface(natm, xa[0], radii, w, area);
	double ASA0[nres], ASA1[nres];
	bzero(ASA0, sizeof(ASA0));
	for(int i=0; i<natm; i++) ASA0[rseq[i]] += area[i];
	double ASA0_tot = sum(nres, ASA0);

	for(int ir=0; ir<nres; ir++){
		double asa0 = 0., asa1=0;;
		for(int i=rstart[ir]; i<rstart[ir+1]; i++){
			asa0 += area[i];
			if(i>=rstart[ir]+5) asa1 += area[i];
		}
		printf("%d %c %d %.3f %.3f\n", seq0[ir], rnam1_std[resid[ir]], rstart[ir+1]-(rstart[ir]+5), asa0, asa1);
	}
	exit(1);


	double radii2[natm];
	for(int ir=0; ir<nres; ir++){
		int n1 = 0;
		for(int i=0; i<natm; i++){
			if(i>=rstart[ir]+5 && i<rstart[ir+1]) continue;
			double *xp = atoms[i].getx();
			for(int m=0; m<3; m++) xa[n1][m] = xp[m];
			radii2[n1] = radii[i];
			n1 ++;
		}
		bzero(area, sizeof(area));
		surface(n1, xa[0], radii, w, area);
		double a1 = sum(n1, area);
		ASA1[ir] = ASA0_tot - a1;
		printf("%d %c %d %.3f %.3f\n", seq0[ir], rnam1_std[resid[ir]], natm-n1, ASA0[ir], ASA1[ir]);
	}
}

int main(int argc, char *argv[]){
	Molecule *mol = new Molecule(argv[1]);
	mol -> calASA();
}
