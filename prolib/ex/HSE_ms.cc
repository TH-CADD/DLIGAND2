#include "protein.h"

class MS_cent{
	double xint[20][3];
public:
	MS_cent(string fn);
	double *getxint(int id){return xint[id];}
};
MS_cent::MS_cent(string fn){
	FILE *fp = openfile(getenv("DATADIR") + fn, "r");
	char str[121], ss[10][15];
	bzero(xint, sizeof(xint));
	while(fgets(str, 120, fp) != NULL){
		if(str[0] == '#') continue;
		int nt = str2dat(str, ss);
		if(nt <= 2) continue;
		int id = resDefine(ss[0]);
		for(int m=0; m<3; m++) xint[id][m] = strtod(ss[m+1], NULL);
	}
}
class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
	void calms_center();
	void calHSEa();
};
void Molecule::calms_center(){
	static MS_cent *ms_aa=new MS_cent("ms_cent.dat");
	initneib();
/*	for(int ir=0; ir<nres; ir++){
		int ib = rstart[ir]+4;
		if(rnam[ir] == "GLY") continue;
		initAtomint(ib); initAtomintcrd0(ib);
		if(! atoms.at(ib).filled()) atoms.at(ib).xyzatm(atoms);
	}*/
	for(int ir=0; ir<nres; ir++){
		int i0 = rstart[ir];
		int ib = i0 + 4; 
		if(rnam[ir] == "GLY"){
			addCB_GLY(ir);
		} else {
			string an = atoms[i0+2].getname();
			if(an != "C") die("%d: %s != C\n", an.c_str());
			int iz[4] = {i0+1, i0, i0+2, 1};
			atoms[ib].setint( iz );
			double *xp = ms_aa -> getxint( resid[ir]  );
			atoms[ib].setzcrd(xp);
			bool lms=true;
			for(int m=0; m<3; m++)
				if(! atoms.at(iz[0]).filled()) lms=false;
			if(! lms && atoms.at(ib).filled()) {
				fprintf(stderr, "Warning, cb replace ms in %d\n", ir);
			}
		}
		atoms.at(ib).xyzatm1(atoms);
		if(! atoms.at(ib).filled()) cout<<"not filled:"<<ir<<endl;
	}
}
double cutoff = 13.;
bool balpha = false;
void Molecule::calHSEa(){
	int hse[3][nres]; bzero(hse, sizeof(hse));
	double vsc[nres][3], xd[3];
	int icb[nres];
	for(int ir=0; ir<nres; ir++){
		xdiff(rstart[ir]+1, rstart[ir]+4, vsc[ir]);

		icb[ir] = rstart[ir]+1;
		if(!balpha && rnam[ir]!="GLY") icb[ir] = rstart[ir]+4;
	}
	for(int ir=0; ir<nres; ir++)
	for(int jr=ir+1; jr<nres; jr++){
		xdiff(icb[ir], icb[jr], xd);
/*		if(dot_product(xd, xd) < 5.*5.) {
			if(atoms.at(rstart[ir]+4).filled() && atoms.at(rstart[jr]+4).filled())
					printf("dis %d %d %f\n", seq0[ir], seq0[jr], sqrt(dot_product(xd, xd)));
		}*/

		if(dot_product(xd, xd) > cutoff*cutoff) continue;
		hse[0][ir] ++; hse[0][jr] ++;

		double dt = dot_product(xd, vsc[ir]);
		if(dt > 0) hse[1][ir] ++;
		else hse[2][ir] ++;
		dt = -dot_product(xd, vsc[jr]);
		if(dt > 0) hse[1][jr] ++;
		else hse[2][jr] ++;
	}
//
	cout<<pdbnm.c_str()<<endl;
	for(int ir=0; ir<nres; ir++){
		int id = resid[ir];
		if(! atoms.at(icb[ir]).filled()) continue;
		printf("%s %c %d %d %d\n", seq0_str[ir].c_str(), rnam1_std[id], hse[0][ir], 
				hse[1][ir], hse[2][ir]);
	}
}
int main(int argc, char *argv[]){
	if(argc < 2){
		fprintf(stderr, "Usage: %s [-alpha] [-cutoff cutoff] PDBs\n", argv[0]);
		fprintf(stderr, "OUTPUT: CNT, HSEAU, HSEAD\n"
				"Attn: set env DATADIR=[/data1/yueyang/source/lib/]\n");
		exit(0);
	}
	int iflag= 0, it;
	it = findargs(argc, argv, "-alpha");
	if(it > 0){
		balpha = true; iflag = max(iflag, it);
	}
	it = findargs(argc, argv, "-cutoff");
	if(it > 0){
		cutoff = strtod(argv[it+1], NULL);
		iflag = max(iflag, it + 1);
	}
	for(int i=iflag+1; i<argc; i++){
		Molecule *mol=new Molecule(argv[i]);
		mol -> calms_center();
		mol -> calHSEa();
		delete mol;
	}
}
