#include "protein.h"
#include "array1.h"

class Molecule:public Protein{
public:
	Molecule(const string str):Protein(str){};
	void calrmin(char);
	double calrmin_res(int,int);
};
double Molecule::calrmin_res(int ir, int jr){
	double rmin2 = 1.e6;
	for(int i=rstart[ir]; i<rstart[ir+1]; i++)
	for(int j=rstart[jr]; j<rstart[jr+1]; j++){
		if(! isfilled(i) || ! isfilled(j)) continue; 
		double r2 = distance2(i, j);
		if(rmin2 > r2) rmin2 = r2;
	}
	return sqrt(rmin2);
}
void Molecule::calrmin(char ch1){
	if(ch1 == '*') ch1 = chainID[0];
	for(int ir=0; ir<nres; ir++){
		if(chainID[chnseq[ir]] != ch1) continue;
		double rmin=1.e3; int imin=-1;
		for(int jr=0; jr<nres; jr++){
			if(chnseq[jr] == chnseq[ir]) continue;
			double r = calrmin_res(ir, jr);
			if(rmin > r) {rmin = r; imin = jr;}
		}
		if(imin>=0) printf("%d%c %.1f %s\n", seq0[ir], rnam1_std[resid[ir]], rmin, resinfo[imin].c_str());
	}
//	printf("%c %s", ch1, s1.c_str());
}
int main(int argc, char *argv[]){
	if(argc < 2) die("Usage: %s pdb [ch]\n", argv[0]);
	char ch1 = '*';
	if(argc >= 3) ch1 = argv[2][0];
	Molecule *mol=new Molecule(argv[1]);
	mol -> calrmin(ch1);
	delete mol;
}
