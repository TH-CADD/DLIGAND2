#ifndef _RESTD
#define _RESTD
#include "misc.h"

extern int DEBUG;
extern string datadir;
static const double PI=acos(-1.), RADIAN=180.0/PI;
//
static const int nres_std = 20;
static const char rnam3_std[][4] = 
		{"ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU",
		"MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"};
static const char rnam1_std0[] = "XACDEFGHIKLMNPQRSTVWY";
static const char *rnam1_std = rnam1_std0 + 1;
inline int aaDefine(string rn0, int DEBUG=0){
	string rn = rn0;
	if(rn=="HSD" || rn=="HSE") rn = "HIS";
	for(int i=0; i<nres_std; i++){
		if(rn==rnam3_std[i]) return i;
	}
	if(DEBUG > 0){
		fprintf(stderr, "Unrecognized residue name: %s\n", rn.c_str());
	}
	return -1;
}
//
inline int aaDefine1(char rn, int DEBUG=0){
/*	static map<char, int> aa1_tab;
	static bool bfirst = 0;
	if(! bfirst) {
		bfirst = 1;
		for(int i=0; i<nres_std; i++) aa1_tab[rnam1_std[i]] = i;
	}
	if(aa1_tab.count(rn) > 0) return aa1_tab[rn];*/
	for(int i=0; i<nres_std; i++){
		if(rn==rnam1_std[i]) return i;
	}
	if(DEBUG > 0){
		fprintf(stderr, "Unrecognized residue name: %c\n", rn);
	}
	return -1;
}

inline int ssDefine(char cs){
	if(cs=='C' || cs=='S' || cs=='T' || cs==' ') return 0;
	if(cs=='H' || cs=='G' || cs=='I') return 1;
	if(cs=='E' || cs=='B' || cs=='b') return 2;
	die("not known ss: %c\n", cs);
	return -1;
}
inline string const convert_hetres(const string &rn3){
	if(rn3 == "MSE") return "MET";
	else if(rn3.substr(0,2)=="CS" || rn3=="CCS") return "CYS";
	return rn3;
}
//
void xdiff( double *xa,  double *xb, double *xab);
double distance2( double *xa,  double *xb);
double angle( double *xa,  double *xb,  double *xc);
double torsion( double *xa,  double *xb,  double *xc,  double *xd);
bool xyzatm(double *xi, double *xa, double *xb, double *xc, double *ag, int chiral);
void getphipsi_bin(int id, int *iphipsi);
#endif
