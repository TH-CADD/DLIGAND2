#ifndef _PROTEIN
#define _PROTEIN
#include "aastd.h"
#include "restype.h"
#include "restpl.h"
using namespace std;

class Atom{
protected:
	bool bfill;
	double x[3];
	string name;
	int type, aseq;
	vector <int> ineib;
//
	int iz[4]; double ag[3];
/*	int ia, ib, ic, chiral;
	double bond, angle1, angle2;*/
public:
	Atom(){bfill=false; type=-1;}
	Atom(int is, int it, string an);
	void setname(string an, int it){name=an; type=it;}
	void settype(int it){type = it;}
	void setint(int *it){for(int i=0; i<4; i++) iz[i] = it[i];}
	void setchiral(int it){iz[3] = it;}
	void setzcrd(double *at){for(int i=0; i<3; i++) ag[i] = at[i];}
	void setfill(bool bf){bfill = bf;}
	void setx(string scrd){ str2dat(scrd, 3, x); bfill = true; };
	void setx(const double *xt){for(int m=0; m<3; m++) x[m]=xt[m]; bfill=true;};
	void addneib(int it);
	void clearNeib(){ineib.clear();}
//	bool xyzatm1(vector <Atom> &ap);
	void settorsion(double at){ag[2] = at;}
//
	double *getx(){return x;};
	string &getname(){return name;};
	bool filled(){return bfill;};
	int gettype(){return type;};
	vector<int>& getineib(){return ineib;};
	int *getint(){return iz;}
	double *getzcrd(){return ag;}
	double gettorsion(){return ag[2];}
//
	friend class Molecule;
};
class Protein{
protected:
	string pdbnm;
	int nres, natm, nchain;
// chains
	vector <char> chainID;
// residues-
	vector <string> rnam, resinfo;
	vector <int> resid, seq0, chnseq;
// res & atoms
	vector <int> rstart;
// atoms
	vector <int> pdbidx;	//pos of pdb, for decoys trj
	vector <int> rseq;		// natm
	vector <Atom> atoms;		//natm
public:
	Protein(){reset();}
	Protein(FILE *fp){reset(); rdpdb(fp);}
	Protein(const string fn);
	vector<int> &getseq0(){return seq0;}
	vector<int> &getresid(){return resid;}
	string getfasta();
	void rdpdb(FILE*);
	void reset();
	bool addRes(int);
	bool addRes(string);
	bool initAtom(string);
	void addAtom(string);
	void initneib();
	void initneib1_sc(int ir);		// for RES mutation only, Be careful!
	void initAtomint(int ia);
	void initAtomzcrd0(int i0);
	void initint(){for(int i=0; i<natm; i++) initAtomint(i);}
	void initzcrd0(){for(int i=0; i<natm; i++) initAtomzcrd0(i);}
	void calAtomzcrd(int ia);
	void insertAtoms(int idx, int na);
	void eraseAtoms(int idx, int na);
	void resMutation(int ir, string rn);
	void popResidues(int nr);
	void mutAA(int ir, int id);
	int mutRotamer(int ir, int irot);
	void addCB_GLY(int ir);
	bool makexyz1(int ia);
	void makexyz_part(int i1, int i2){for(int i=i1; i<i2; i++) makexyz1(i);}
//
	void wrpdb(FILE*);
	void wrpdb(string fn);
	void xdiff(int ia, int ib, double *xd);
	double distance2(int ia,int ib);
	double angle(int,int,int);
	double torsion(int,int,int,int);
	bool iscomplete();
	bool refill();
	double calResclash(int ir);
	int getReschi(int ir, int *idx);
	int findatom(int ires, string an);
	int getnres(){return nres;};
	void calphipsi(int ir, double*ag);
	bool isfilled(int ia){return atoms[ia].filled();}
//
	double *getx(int ia){return atoms[ia].getx();}
	void setx(int ia, double *x0){atoms[ia].setx(x0);}
	int getatype(int ia){return atoms[ia].gettype();}
//
	friend class Molecule;
};
#endif
