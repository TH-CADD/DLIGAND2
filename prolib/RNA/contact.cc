#include "protein.h"
// print distances between residue pairs < 2*d0

class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
	void calculate();
};
double cutoff = 8.; int seqcut = 5;
bool balpha = true, bseq0 = 0, bpair=false;
void Molecule::calculate(){
//	string seq1 = getfasta();
	int ncnt=0, total=0, nr=0;
	int idx1[nres];
	for(int ir=0; ir<nres; ir++){
		idx1[ir] = findatom(ir, "C1'");
		if(idx1[ir]>=0 && ! isfilled(idx1[ir])) idx1[ir] = -1;
		if(idx1[ir] >= 0) nr ++;
	}
	fprintf(stderr, "effective/ nres: %d %d\n", nres, nr);
	if(! bpair){
		printf("# ");
		for(int ir=0; ir<nres; ir++){
			if(idx1[ir] < 0) continue;
			char rn1 = rnam[ir][rnam[ir].size()-1];
			printf("%c%d_%d ", rn1, ir, seq0[ir]);
		}
		printf("\n");
//
		for(int ir=0; ir<nres; ir++){
			if(idx1[ir] < 0) continue;
			for(int jr=0; jr<ir; jr++){
				if(idx1[jr] < 0) continue;
				double r2 = distance2(idx1[ir], idx1[jr]);
				printf(" %.1f", sqrt(r2));
			}
			printf("\n");
		}
		return;
	}
//
	for(int ir=0; ir<nres; ir++)
	for(int jr=0; jr<ir; jr++){
		if(fabs(jr-ir) < seqcut) continue;
		if(! isfilled(idx1[ir]) || ! isfilled(idx1[jr])) continue; 
		double r2 = distance2(idx1[ir], idx1[jr]);
		if(r2 > 4*cutoff*cutoff) continue;
		ncnt ++;
		total += (ir-jr);
		char rn1, rn2;
		if(resid[ir] >= 20) {
			rn1 = rnam[ir][rnam[ir].size()-1];
			rn2 = rnam[jr][rnam[jr].size()-1];
		} else {
			rn1 = rnam1_std[resid[ir]];
			rn2 = rnam1_std[resid[jr]];
		}
		if(bseq0) {
			printf("%c%d %c%d %.1f\n", rn1, seq0[ir], rn2, seq0[jr], sqrt(r2));
		} else {
			printf("%c%d %c%d %.1f\n", rn1, ir, rn2, jr, sqrt(r2));
		}
	}
//	printf("#%s %d %d %d\n", pdbnm.c_str(), nr, ncnt, total);
}
int main(int argc, char *argv[]){
	if(argc < 2){
		fprintf(stderr, "Usage: %s [-pair [-seqcut sc] [-cutoff cutoff]] [-seq0] PDBs\n", argv[0]);
		exit(0);
	}
	bool bopts[argc]; bzero(bopts, sizeof(bopts));
	int it = findargs(argc, argv, "-cutoff");
	if(it > 0){
		cutoff = strtod(argv[it+1], NULL);
		bopts[it+1] = 1;
	}
	DEBUG = findargs(argc, argv, "-v");
	bseq0 = findargs(argc, argv, "-seq0") > 0;
	it = findargs(argc, argv, "-seqcut");
	if(it > 0){seqcut = atoi(argv[it+1]); bopts[it+1] = 1;}
	bpair = findargs(argc, argv, "-pair") > 0;
//
	initRestypes(getdatadir() + "amino+na.dat");
	for(int i=1; i<argc; i++){
		if(argv[i][0] == '-' || bopts[i]==1) continue;
		Molecule *mol=new Molecule(argv[i]);
		mol -> calculate();
		delete mol;
	}
}
