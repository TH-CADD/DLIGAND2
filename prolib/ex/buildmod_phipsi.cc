#include "protein.h"
#include "array1.h"

class Molecule:public Protein{
	vector<double> dat_pred[9], phipsi1[2], *phipsi0;  // phi/psi, theta/tau
	vector<Xvec> list_xca0, list_xca1;
	string seq_spd3;
public:
	Molecule(const string str):Protein(str){};
	Molecule(){};
	void rdspd3(string str);
	void init(string bn);
	void run_res(int ir);
	void runMC();
	void minimize1();

	void build_phi(vector<double>*phipsi,int);
	void build_theta(vector<double>*phipsi);
	void caltheta(vector<double>*tor);
	friend int checkali(Molecule *mol1, Molecule *mol2, vector<int> &iali, string rali[]);
	void getfrag(vector<int> &);
	void getxca(vector<Xvec> &xc);
};
void Molecule::getxca(vector<Xvec> &xc){
	xc.clear();
	for(int i=0; i<nres; i++) xc.push_back( Xvec(atoms[rstart[i]+1].getx()) );
}
void Molecule::build_theta(vector<double>*phipsi){	//theta
	int nr = phipsi[0].size();
	for(int i=0; i<nr; i++) addRes("ALA");
	if(DEBUG) cerr<<"nres: "<<nres<<" natm: "<<natm<<endl;
	chnseq.resize(nr, 0);
/*	initneib();
	for(int i=0; i<natm; i++) initAtomint(i);
	for(int i=0; i<natm; i++) initAtomzcrd0(i);*/
	int i0, i1, i2, i3;
	for(int ir=0; ir<nres; ir++){
		i1 = i2 = i3 = -1;
		if(ir>=3) i3 = rstart[ir-3]+1;
		if(ir>=2) i2 = rstart[ir-2]+1;
		if(ir>=1) i1 = rstart[ir-1]+1;
		i0 = rstart[ir]+1;
		int iz[] = {i1, i2, i3, 0};
		atoms[i0].setint(iz);
		double d1=360, d2=360.;
		if(ir > 0){
			d1 = phipsi[0][ir-1], d2 = phipsi[1][ir-1];
		}
		if(d1 > 359) d1 = 90.;
		if(d2 > 359) d2 = 180.;
		double dz[] = {3.8, d1, d2};
		atoms[i0].setzcrd(dz);
		makexyz1(i0);
	}
//	wrpdb("1.pdb");
}
void Molecule::caltheta(vector<double> *torsions){
	bool bgap[nres];
	for(int i=0; i<nres; i++) bgap[i] = 1;
	for(int i=0; i<nres-1; i++){
		int i1=rstart[i]+1, i2=rstart[i+1]+1;
		if(! (isfilled(i1) && isfilled(i2)) ) continue;
		double r2 = distance2(i1, i2);
		double r = sqrt(r2);
		if(fabs(r-3.8) < 0.5) bgap[i] = 0;
//		else fprintf(stderr, "%d %f\n", i, r);
	}
//
//	double theta[nres], tors[nres];
//	for(int i=0; i<nres; i++) theta[i] = tors[i] = 360;
	for(int m=0; m<2; m++){
		torsions[m].resize(nres, 360.);
		for(int i=0; i<nres; i++) torsions[m][i] = 360.;
	}
	for(int i=0; i<nres-2; i++){
		int i1=rstart[i]+1;
		int i2=rstart[i+1]+1;
		int i3=rstart[i+2]+1;
		torsions[0][i+1] = angle(i1, i2, i3);
		if(i >= nres-3) continue;
		int i4=rstart[i+3]+1;
		torsions[1][i+2] = torsion(i1, i2, i3, i4);
	}
/*	for(int i=0; i<nres; i++){
		int id = resid[i];
		int c = 4;
		if(i<nres-3) c = bgap[i] + bgap[i+1] + bgap[i+2] + bgap[i+3];
		printf("%d %c %.1f %.1f %d\n", i+1, rnam1_std[id], theta[i], tors[i], c);
	}*/
}
void Molecule::rdspd3(string fn){
	string line;
	vector<string> ss;
	FILE *fp = openfile(fn);
	seq_spd3 = "";
	while(getline1(fp, line)){
		if(line[0]=='#') continue;
		int nt = str2dat(line, ss);
		for(int m=0; m<4; m++){
			double dt = atof(ss[4+m].c_str());
			dat_pred[m].push_back(dt);
		}
		seq_spd3 += ss[1];
	}
}
void Molecule::build_phi(vector<double>*phipsi, int ir){
	if(ir < 0){
		if(nres >= 0){
			reset();
			int nr = phipsi[0].size();
			for(int i=0; i<nr; i++) addRes(rnam3_std[resid[i]]);
			if(DEBUG) cerr<<"nres: "<<nres<<" natm: "<<natm<<endl;
			chnseq.resize(nr, 0);
			initneib();
			for(int i=0; i<natm; i++) initAtomint(i);
			for(int i=0; i<natm; i++) initAtomzcrd0(i);
		}
		for(int ir=0; ir<nres; ir++){
			int i0 = rstart[ir];
			atoms[i0+2].settorsion(phipsi[0][ir]);
			atoms[i0+3].settorsion(phipsi[1][ir] + 180.);
		}
		makexyz_part(0, natm);
	} else {
		int i0 = rstart[ir], i1 = rstart[ir+1];
		atoms[i0+2].settorsion(phipsi[0][ir]);
		atoms[i0+3].settorsion(phipsi[1][ir] + 180.);
		makexyz_part(i0+2, i1+2);
	}
}
void Molecule::init(string fn){
	string line;
	vector<string> ss;
	FILE *fp = openfile(fn);
	resid.clear();
	while(getline1(fp, line)){
		if(line[0]=='#') continue;
		if(line.substr(0,3) == "phi"){
			phipsi1[0].clear();
			int nt = str2dat(line, ss);
			for(int i=1; i<nt; i++) phipsi1[0].push_back(atof(ss[i].c_str()));
		} else if(line.substr(0,3) == "psi"){
			phipsi1[1].clear();
			int nt = str2dat(line, ss);
			for(int i=1; i<nt; i++) phipsi1[1].push_back(atof(ss[i].c_str()));
		}
	}
	assert(phipsi1[0].size() == phipsi1[1].size());
	resid.resize(phipsi1[0].size(), 0);
	build_phi(phipsi1, -1);
	fclose(fp);
}
int main(int argc, char *argv[]){
	if(argc < 2) die("usage: RUN phipsi");
	initRestypes_pro();
	Molecule *mol1 = new Molecule();
	mol1 -> init(argv[1]);
	mol1 -> wrpdb("mod1.pdb");
	delete mol1;
}
