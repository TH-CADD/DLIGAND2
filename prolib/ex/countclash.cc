#include "protein.h"

class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
	int countclash();
};
int Molecule::countclash(){
	int nclash = 0;
	double rmin2_res[nres][nres];
	for(int i=0; i<nres; i++)
	for(int j=i+1; j<nres; j++) rmin2_res[i][j] = 10000.;
//
	for(int i=0; i<natm; i++)
	for(int j=i+1; j<natm; j++){
//		if(ir>0 && (i>=rstart[ir] || j<=rstart[ir-1])) continue;
		if(rseq[j]-rseq[i] <= 1) continue;
		if(getatype(i) < 0) continue;
		if(getatype(j) < 0) continue;
		if(! atoms[i].filled() || ! atoms[j].filled()) continue;
		double r2 = distance2(i, j);
		rmin2_res[rseq[i]][rseq[j]] = min(rmin2_res[rseq[i]][rseq[j]], r2);
		if(r2 > 2.*2.) continue;
		printf("%d %s %s -- %d %s %s %.2f\n", rseq[i], resinfo[rseq[i]].c_str(), atoms[i].getname().c_str(), rseq[j], resinfo[rseq[j]].c_str(), atoms[j].getname().c_str(), sqrt(r2));
		nclash ++;
	}
	printf("number of clashes: %d\n\n", nclash);
	for(int i=0; i<nres; i++)
	for(int j=i+2; j<nres; j++){
		double r2 = rmin2_res[i][j];
		if(r2 > 2.*2.) continue;
		printf("%s -- %s: %f\n", (resinfo[i]).c_str(), (resinfo[j]).c_str(), sqrt(r2));
	}
	return nclash;
}
int main(int argc, char *argv[]){
	for(int i=1; i<argc; i++){
		Molecule *mol=new Molecule(argv[i]);
		mol->countclash();
		delete mol;
	}
}
