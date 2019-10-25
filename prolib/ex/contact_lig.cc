#include "protein3.h"

double RCUT = 0.5, RCUT2=-1;
bool bprint_atom = 0;
map<string, double> vdw_dat;
set<string> vdw_unk;
void rdvdw(const string fn){
	FILE *fp = openfile(fn);
	string line, ss[4];
	while(getline1(fp, line)){
		if(line[0] == '#') continue;
		str2dat(line, 2, ss);
		vdw_dat[ss[0]] = atof(ss[1].c_str());
	}
}
double getvdw(string an0){
	if(vdw_dat.size() < 1) {
		if(getenv("DATADIR") != NULL) rdvdw(string(getenv("DATADIR")) + "vdw.dat");
	}
	string an = an0.substr(0,2);
	if(an[0] == 'H') return 0.;
	char a1 = an[0];
	if(vdw_dat.count(an) > 0) return vdw_dat[an];
	else if(a1 == 'P') return 1.80;
	else if(a1 == 'F') return 1.47;
	else if(a1 == 'S') return 1.80;
	else if(a1 == 'C') return 1.70;
	else if(a1 == 'N') return 1.55;
	else if(a1 == 'O') return 1.52;
	else if(a1 == 'I') return 1.98;
	if(vdw_unk.count(an) == 0){
		fprintf(stderr, "unknown atom type for vdw: %s\n", an.c_str());
		vdw_unk.insert(an);
	}
	return 1.0;
	
/*	if(an=="CL") return 1.75;
	else if(an=="ZN") return 1.39;
	else if(an=="MG") return 1.73;
	else if(a1 == 'P') return 1.80;
	else if(a1 == 'F') return 1.47;
	else if(a1 == 'S') return 1.80;
	else if(a1 == 'C') return 1.70;
	else if(a1 == 'N') return 1.55;
	else if(a1 == 'O') return 1.52;
	if(vdw_unk.count(an) == 0){
		fprintf(stderr, "unknown atom type for vdw: %s\n", an.c_str());
		vdw_unk.insert(an);
	}*/
	return 0.;
}
class Molecule:public Protein3{
public:
	Molecule(const string str):Protein3(str, 1){};
	void test(const string);
	bool contact_res(int ir, int jr);
};
bool Molecule::contact_res(int ir, int jr){
	vector<Atom3> &pa = residues[ir].getatoms();
	vector<Atom3> &pb = residues[jr].getatoms();
	for(int i=0; i<pa.size(); i++)
	for(int j=0; j<pb.size(); j++){
		double r1 = getvdw(pa[i].getname()), r2 = getvdw(pb[j].getname());
		if(! (pa[i].isfilled() && pb[j].isfilled())) continue;
		double r = sqrt(distance2(pa[i].getx(), pb[j].getx()));
		if(RCUT2 > 0) {
			if(r < RCUT2) return 1;
		} else {
			if(r < r1+r2+RCUT) return 1;
		}
	}
	return 0;
}
void Molecule::test(const string fn2){
	int nres0 = nres;
	rdpdb(fn2);
	assert(nres > nres0);

	string sinfo = "", sinfo0 = ""; char s1[20];
	for(int ir=0; ir<nres0; ir++){
		bool bcont = 0;
		for(int jr=nres0; jr<nres; jr++){
			if(contact_res(ir, jr)) {bcont = 1; break;}
		}
		if(! bcont) continue;
		sprintf(s1, "%d ", ir);
		sinfo += rnam1_std[residues[ir].gettype()] + string(s1);
		sprintf(s1, "%d ", residues[ir].getrseq0());
		sinfo0 += rnam1_std[residues[ir].gettype()] + string(s1);
//		cout<<rnam1_std[residues[ir].gettype()]<<ir<<' ';

		if(! bprint_atom) continue;
		vector<Atom3> &pa = residues[ir].getatoms();
		for(int i=0; i<pa.size(); i++){
			printf("%s", pa[i].pdbline.c_str());
		}
	}

	cout<<">\n";
	for(int ir=0; ir<nres0; ir++) {
		cout<<rnam1_std[residues[ir].gettype()];
	}
	cout<<endl;

	cout<<"#bnd_idx: "<<sinfo<<endl;
	cout<<"#bnd_idx0: "<<sinfo0<<endl;
}
int main(int argc, char *argv[]){
	if(argc < 3) die("usage: RUN pdb lig [-c CUT(v1+v2+0.5) | -c2 CUT2(3.5) ] [--print_atom]");

	DEBUG = findargs(argc, argv, "-v") > 0;
	int it;
	if((it=findargs(argc, argv, "-c")) > 0) {RCUT = atof(argv[it+1]);}
	if((it=findargs(argc, argv, "-c2")) > 0) {RCUT2 = atof(argv[it+1]);}
	if((it=findargs(argc, argv, "--print_atom")) > 0) bprint_atom = 1;

	Molecule *mol = new Molecule(argv[1]);
	mol -> test(argv[2]);
}
