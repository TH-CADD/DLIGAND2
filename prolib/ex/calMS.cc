#include "protein.h"

class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
	void calCM();
};
void Molecule::calCM(){
	bool lcomp[nres];
	for(int ir=0; ir<nres; ir++){
		lcomp[ir] = true;
		for(int i=rstart[ir]; i<rstart[ir+1]; i++){
			if(! atoms[i].filled()) {lcomp[ir]=false; break;}
		}
	}
//
	for(int ir=0; ir<nres; ir++){
		if(! lcomp[ir]) continue;
		int nsc = rstart[ir+1] - rstart[ir] - 4;
		if(nsc < 1) continue;
		int icb = rstart[ir]+4;
		double xm[3], xcb[3], *xp;
		xp = atoms[icb].getx();
		for(int m=0; m<3; m++) xcb[m] = xp[m];
		bzero(xm, sizeof(xm));
//
		for(int i=rstart[ir]+4; i<rstart[ir+1]; i++){
			xp = atoms[i].getx();
			for(int m=0; m<3; m++) xm[m] += xp[m];
		}
		for(int m=0; m<3; m++) xm[m] /= nsc;
		atoms[icb].setx(xm);
		int i0 = rstart[ir];
		double ag[10];
		ag[0] = sqrt( distance2(icb, i0+1) );
		ag[1] = angle(icb, i0+1, i0);
		ag[2] = angle(icb, i0+1, i0+2);
		ag[3] = torsion(icb, i0+1, i0, i0+2);
		ag[4] = angle(i0, i0+1, i0+2);
		printf("%d %s ", ir, rnam[ir].c_str());
		prtdim(5, ag);
	}
}
int main(int argc, char *argv[]){
	for(int i=1; i<argc; i++){
		Molecule *mol=new Molecule(argv[i]);
		cout<<argv[i]<<endl;
		mol -> calCM();
		delete mol;
	}
}
