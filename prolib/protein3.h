#ifndef _PROTEIN3
#define _PROTEIN3
#include "aastd.h"
#include "restype.h"
#include "restpl.h"
#include "rotamer.h"
using namespace std;

class Residue;
class Protein3;
class Atom3{
protected:
	bool bfill;
	double x[3];
	string name, info, pdbline;
	int type;
	vector<Residue> *res_all;
	vector<int> neib[2]; // res & atom idx
//	
	vector<int> iz[2]; int chiral;
	double ag[3];
public:
	Atom3(int it, const string &an){init(); type=it; name=an;}
	Atom3(){init();}
	void init(){bfill=false; res_all= NULL;}
	bool makexyz1();
//	void setresidues(vector<Residue> *res){res_all = res;}
	void setx(const string &scrd){ str2dat(scrd, 3, x); bfill = true; };
	void setx(const double *xt){for(int m=0; m<3; m++) x[m]=xt[m]; bfill=true;};
	void setint(int ic, int p[][3]){chiral = ic; for(int m=0; m<3; m++) {iz[0][m]=p[0][m]; iz[1][m] = p[1][m];}}
	void delete_neib1(int);
	void clear_neib(){neib[0].clear(); neib[1].clear();}
	void initint(int ir0, int aseq);
	void addneib(int ir, int ia);
	void setchiral(int it){chiral = it;}
	void setzcrd(double *ag0){for(int m=0; m<3; m++) ag[m] = ag0[m];}
	void setztor(double a0){ag[2] = a0;}
	void setinfo(const string &s){info = s;}
	void setline(const string &s){pdbline = s;}
	void setname(const string &an, int it){name=an; type=it;}
// no change to the class
	double calztor();
	void calzcrd();
	Atom3 &getatom(int ir, int ia);
	void prtneib(int);
	void prtint();
//
	const char *getname_str(){return name.c_str();}
	string &getname(){return name;}
	bool isfilled(){return bfill;}
	double *getx(){return x;}
	double *getzcrd(){return ag;}
	double getztor(){return ag[2];}
	int getchiral(){return chiral;}
	int gettype(){return type;}
	vector<int> *getneib(){return neib;}
	vector<int> *getiz(int &ch){ch = chiral; return iz;}
	vector<Residue>* &getresidues(){return res_all;}
	friend class Protein3;
	friend class Molecule;
};
class Residue{
protected:
	string name, info;
	int type, rseq, rseq0; // residue type; res. index in chain
	vector<Atom3> atoms;
//
	vector<double> phipsi;
	vector<int> idx_chi;
public:
	Residue(){}
	Residue(int id, int idx){init(id, idx);}
	Residue(string rn0, int idx){init(rn0, idx);}
	void init(int id, int idx);
	void init(string id, int idx);
	void init_het(const string rn, int idx){name=rn; type=-1; rseq=idx;}
	bool initAtom_pdb(const string&);
	void initneib_sc();
	void initAtomint(int);
	void initAtomzcrd0();
	void mutation(int);
	int mutRotamer(int irot);
	int mutRotamer(int irot, double&);
	void mutRotamer(const vector<double> &chi);
	Rotamers &determine_Rotamers();
	void setresidues(vector<Residue> *res){for(int i=0; i<atoms.size(); i++) atoms[i].getresidues() = res;}
// 
	vector<Atom3> &getatoms(){return atoms;}
	int findatom(string an, bool);
	void calphipsi();
	void calidx_chi();
	void calchi();
	void getchi(vector<double>& chi);
	double calProt_nat();
	void prtneib(){cout<<'#'<<name<<endl;for(int m=0; m<atoms.size(); m++) atoms[m].prtneib(rseq);}
//
	const char *getname_str(){return name.c_str();}
	string &getname(){return name;}
	string &getinfo(){return info;}
	int gettype(){return type;}
	int getrseq0(){return rseq0;}
	vector<double> &getphipsi(){return phipsi;}
	friend class Protein3;
	friend class Molecule;
};
class Protein3{
protected:
	bool bhet;
	string name, chainID;
	int nres, natm, nchain;
	vector<Residue> residues;
	vector<int> chnseq, chnstart; // temporarily for multi-chains
public:
	Protein3(){reset();}
	Protein3(const string fn, bool bh=0);
	Protein3(FILE *fp, bool bh=0){reset(); bhet=bh; rdpdb(fp);}
	void reset();
	void rdpdb(FILE*);
	void rdpdb(const string);
	bool addRes(const string rn0);
	bool addRes_het(const string rn0);
	void initneib();
	void initAtomint();
	void initAtomzcrd0();
	void resMutation(int ir, const string rn);
//
	void wrpdb(FILE*);
	void wrpdb(const string fn);
	void prtneib(){for(int i=0; i<residues.size(); i++) residues[i].prtneib();}
//
	friend class Molecule;
};
#endif
