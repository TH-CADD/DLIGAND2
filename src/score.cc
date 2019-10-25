#include "molecule.h"
#include "misc.h"

void rdlist(string fn, vector<string> &list){
	char str[201];
	FILE *fp = openfile(fn, "r");
	while(fgets(str, 200, fp) != NULL){
		sscanf(str, "%s", str); list.push_back(str);
	}
	fclose(fp);
}
int main(int argc, char *argv[]){
	if(argc <= 1) {
		fprintf(stderr, "Usage: %s [-v] [-s dfire.2] [-d dir] [-etype 1,2] [-l list|pdbs]\n", argv[0]);
		exit(0);
	}
	vector<string> pdblist;
	string sdir = "", slib = "";
	for(int i=1; i<argc; i++){
		string opt1 = argv[i];
		if(opt1 == "-v") DEBUG ++;
		else if(opt1 == "-d") {
			sdir = argv[++i];
		}else if(opt1 == "-s"){
			slib = argv[++i];
		}else if(opt1 == "-etype"){
			etype = strtol(argv[++i], NULL, 10);
		}else if(opt1 == "-l"){
			i++; rdlist(argv[i], pdblist);
		} else if(argv[i][0] == '-') die("not known opt: %s\n", argv[i]);
		else pdblist.push_back(argv[i]);
	}
	if(pdblist.size() < 1) die("no PDB input\n");
//
	void initDFIRE(string fn); initDFIRE(slib);
//
	for(int i=0; i<pdblist.size(); i++){
/*		if(argv[i] == string("-V")) DEBUG ++;
		if(argv[i][0] == '-') continue;*/
		string bn = pdblist[i];
		string fn = sdir + bn + "_protein.pdb";
		if(! file_existed(fn)){
			printf("%s not exists", fn.c_str()); continue;
		}
//		fprintf(stderr, "\r%d/%d %s", i, pdblist.size(), bn.c_str());
		Molecule *mol = new Molecule(fn);
		fn =  sdir + bn + "_ligand.mol2";
		mol -> rdmol2(fn);
//		mol->initneib();
		cout<<bn<<' ';
		double e0 = mol->score();
//		cout<<bn<<' '<<e0<<endl;
		delete mol;
	}
}
