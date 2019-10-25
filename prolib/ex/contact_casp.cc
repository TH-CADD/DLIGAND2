#include "protein.h"
// print distances between residue pairs < 2*d0

class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
	void calculate();
};
double cutoff = 8.; int seqcut = 5;
bool bseq0 = 0;
void Molecule::calculate(){
	string seq1 = getfasta();
	printf("#%d %s\n", seq1.size(), seq1.c_str());
	int ncnt=0, total=0, nr=0;
	int iCB[nres];
	for(int ir=0; ir<nres; ir++){
		if(rnam1_std[resid[ir]] == 'G') iCB[ir] = findatom(ir, "CA");
		else iCB[ir] = findatom(ir, "CB");
	}
/*	for(int ir=0; ir<nres; ir++)
		if(isfilled(iCB[ir])) nr ++;
*/
	if(bseq0){
		printf("#");
		for(int ir=0; ir<nres; ir++) {
			if(! isfilled(iCB[ir])) continue;
			printf(" %c%d", rnam1_std[resid[ir]], seq0[ir]);
		}
		printf("\n");
	}
	for(int ir=0; ir<nres; ir++)
	for(int jr=ir+1; jr<nres; jr++){
		if(jr-ir < seqcut) continue;
		if(! isfilled(iCB[ir]) || ! isfilled(iCB[jr])) continue; 
		double r2 = distance2(iCB[ir], iCB[jr]);
		if(r2 > 4*cutoff*cutoff) continue;
		ncnt ++;
		total += (jr-ir);
		if(bseq0) {
			printf("%c%d %c%d %.1f\n", rnam1_std[resid[ir]], seq0[ir], rnam1_std[resid[jr]], seq0[jr], sqrt(r2));
		} else {
			printf("%c%d %c%d %.1f\n", rnam1_std[resid[ir]], ir, rnam1_std[resid[jr]], jr, sqrt(r2));
		}
	}
//	printf("#%s %d %d %d\n", pdbnm.c_str(), nr, ncnt, total);
}
int main(int argc, char *argv[]){
	if(argc < 2){
		fprintf(stderr, "Usage: %s [-seqcut sc] [-cutoff cutoff] [-seq0] PDBs\n", argv[0]);
		exit(0);
	}
	int iflag= 0, it;
	it = findargs(argc, argv, "-cutoff");
	if(it > 0){
		cutoff = strtod(argv[it+1], NULL);
		iflag = max(iflag, it + 1);
	}
	it = findargs(argc, argv, "-seq0");
	if(it > 0){iflag = max(iflag,it); bseq0 = 1;}
	it = findargs(argc, argv, "-seqcut");
	if(it > 0){iflag = max(iflag,it+1); seqcut = atoi(argv[it+1]);}
//
	for(int i=iflag+1; i<argc; i++){
		Molecule *mol=new Molecule(argv[i]);
		mol -> calculate();
		delete mol;
	}
}
