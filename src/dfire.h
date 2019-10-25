#ifndef _DFIRE_H
#define _DFIRE_H

// constant
const int mbin=30, mabin=6, mpolar=60, nbin=mbin;
const int matype_mol2 = 11, matype_pro=167, matype=matype_pro + matype_mol2;
const int mkind = matype;
const double Rcut0=mbin*0.5, Rcut=nbin*0.5, dRcut=0.5;
const double ert=0.01, epen1=10.*ert, epen2=2.*ert;
//
extern int etype, Atype;
extern double ALPHA, rcut;
extern bool bsym, bag;
extern int nkind, iskind[matype];
extern float edfire[mkind][mkind][mbin], Nobs0[mkind][mbin];
//
inline int r2bin(double r){
	if(r >= Rcut0) return -1;
	else return int(r/0.5);
}
inline double bin2r(int i){
	return (i+0.5) * 0.5;
}
inline int angle2bin(double ag){
	const double PI = acos(-1.), radian=180./PI;
   static double agbin[mabin] = {-2./3,-1./3, 0., 1/3., 2/3., 1.};
	if(fabs(ag)-1.>1.0e-5){
		fprintf(stderr, "Warning, error angle: %f\n", ag);
	}
	for(int i=0; i<mabin; i++)
		if(ag<=agbin[i]) return i;
	return mabin-1;
}
//
int mol2Define(string s1);

#endif
