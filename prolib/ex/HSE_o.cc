#include "protein.h"

class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
	void addCB();
	void calHSEa(FILE*);
};
void Molecule::addCB(){
	initneib();
	for(int ir=0; ir<nres; ir++){
		int ib = rstart[ir]+4;
		if(rnam[ir] == "GLY") continue;
		initAtomint(ib); initAtomzcrd0(ib);
		if(! atoms.at(ib).filled()) makexyz1(ib); //atoms.at(ib).xyzatm(atoms);
	}
	for(int ir=0; ir<nres; ir++){
		if(rnam[ir] != "GLY") continue;
		addCB_GLY(ir);
		int ib = rstart[ir]+4;
		if(! atoms.at(ib).filled()) makexyz1(ib);
	}
}
double cutoff = 13.;
bool balpha = true;
void Molecule::calHSEa(FILE *fp){
	int hse[3][nres]; bzero(hse, sizeof(hse));
	double vsc[nres][3], xd[3];
	for(int ir=0; ir<nres; ir++){
		xdiff(rstart[ir]+1, rstart[ir]+4, vsc[ir]);
	}
	for(int ir=0; ir<nres; ir++)
	for(int jr=ir+1; jr<nres; jr++){
		if(balpha) xdiff(rstart[ir]+1, rstart[jr]+1, xd);
		else xdiff(rstart[ir]+4, rstart[jr]+4, xd);

		if(dot_product(xd, xd) > cutoff*cutoff) continue;
		hse[0][ir] ++; hse[0][jr] ++;

		double dt = dot_product(xd, vsc[ir]);
		if(dt > 0) hse[1][ir] ++;
		else hse[2][ir] ++;
		dt = -dot_product(xd, vsc[jr]);
		if(dt > 0) hse[1][jr] ++;
		else hse[2][jr] ++;
	}
	cout<<pdbnm.c_str()<<endl;
	for(int ir=0; ir<nres; ir++){
		int id = resid[ir];
		fprintf(fp, "%s %c %d %d %d\n", seq0_str[ir].c_str(), rnam1_std[id], hse[0][ir], 
				hse[1][ir], hse[2][ir]);
	}
}
int main(int argc, char *argv[]){
	if(argc < 2){
		fprintf(stderr, "Usage: %s [-Beta] [-cutoff cutoff] PDBs\n", argv[0]);
		fprintf(stderr, "OUTPUT: CNT, HSEAU, HSEAD\n"
				"Attn: set env DATADIR=[/data1/yueyang/source/lib/]\n");
		exit(0);
	}
	int iflag= 0, it;
	it = findargs(argc, argv, "-Beta");
	if(it > 0){
		balpha = false; iflag = max(iflag, it);
	}
	it = findargs(argc, argv, "-cutoff");
	if(it > 0){
		cutoff = strtod(argv[it+1], NULL);
		iflag = max(iflag, it + 1);
	}
	for(int i=iflag+1; i<argc; i++){
		string fno = argv[i] + string(".hse");
		FILE *fpo = fopen(fno.c_str(), "r");
		if(fpo != NULL) {fclose(fpo); continue;}
		fpo = fopen(fno.c_str(), "w");
		Molecule *mol=new Molecule(argv[i]);
		mol -> addCB();
		mol -> calHSEa(fpo);
		delete mol;
	}
}
