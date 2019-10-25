#include "molecule.h"
#include "misc.h"
#include "dfire.h"

int main(int argc, char *argv[]){
	vector<string> pdblist;
	for(int i=1; i<argc; i++){
		string opt1 = argv[i];
		if(opt1 == "-v") DEBUG ++;
		else if(opt1 == "-nosym") bsym = false;
		else if(opt1 == "-noag") bag = false;
		else if(argv[i][0] == '-') die("not known opt: %s\n", argv[i]);
		pdblist.push_back(argv[i]);
	}
	if(pdblist.size() < 1) die("no PDB input\n");
//
	bzero(edfire, sizeof(edfire));
	void rdpolar(string); rdpolar("");
	for(int i=1; i<argc; i++){
		if(argv[i] == string("-v")) DEBUG ++;
		if(argv[i][0] == '-') continue;
		fprintf(stderr, "\r%d/%d %s", i, argc, argv[i]);
		string bn = argv[i];
		string fn = bn + '/' + bn + "_protein.pdb";
		if(! file_existed(fn)){
			fprintf(stderr, "%s not exists", fn.c_str()); continue;
		}
		Molecule *mol = new Molecule(fn);
		fn =  bn + '/' + bn + "_ligand.mol2";
		mol -> rdmol2(fn);
//		mol->initneib();
		mol -> stats();
		delete mol;
	}
	void prtobs(string); prtobs("Nobs1");
}
void prtobs(string fn){
	FILE *fp = openfile(fn, "w");
	for(int i=0; i<mbin; i++){
		if(i==0) fprintf(fp,"%d: %d\n", i, nkind);
		else fprintf(fp,"%d:\n", i);
		for(int j=0; j<nkind; j++){
			for(int k=0; k<nkind; k++) fprintf(fp, " %d", int(edfire[j][k][i]+0.5));
			fprintf(fp,"\n");
		}
	}
	fprintf(fp, "Nobs0::\n");
	for(int i=0; i<mbin; i++){
		for(int j=0; j<nkind; j++) fprintf(fp, " %d", int(Nobs0[j][i]));
		fprintf(fp, "\n");
	}
	fclose(fp);
}
