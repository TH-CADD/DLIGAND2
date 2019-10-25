#include "restype.h"
#include "dfire.h"
#include "misc.h"

// global variables for DFIRE
int etype = 1;  // 1;167*14   2; 14*14
int Atype = 0;	//0,all; 1,CA; 12,CB; 2,CA+CB; 5,main-chain
double rcut = 15.;
double ALPHA=1.61;
bool bsym=true, bag=false;
int nkind, iskind[matype], nkind2, ismol2[matype];
float edfire[mkind][mkind][mbin], Nobs0[mkind][mbin];
//
void initDFIRE(string fn){
	void rdpolar(string); rdpolar("");
	void initEobs(string); initEobs(fn);

}
int mol2Define(string s1){
/*	char mol2_type[matype_mol2][6] = {"C.2", "C.3", "C.ar", "C.cat", "N.4",
		"N.am", "N.pl3", "O.2", "O.3", "O.co2", "S.3", "P.3", "N.2", "N.ar"
	};*/
	char mol2_type[matype_mol2][6] = {"C.2", "C.3", "C.ar", "C.cat",
		"O.2", "O.3", "O.co2", "N.4", "N.am", "N.pl3", "S.3" //, "P.3", "N.2", "N.ar"
	};
	if(s1.substr(0,2) == "S.") s1 = "S.3";
	else if(s1=="P.3" || s1=="Cl" || s1=="Br" || s1=="Met") s1 = "S.3";

	else if(s1 == "F") s1 = "O.co2";
	else if(s1=="C.1") s1 = "C.3";

	else if(s1=="N.3" || s1=="N.1") s1 = "N.2";
	if(s1=="N.2" || s1=="N.ar") s1 = "N.pl3";

	for(int i=0; i<matype_mol2; i++){
		if(s1 == mol2_type[i]) return i;
	}
	if(DEBUG > 0) fprintf(stderr, "not known atom: %s\n", s1.c_str());
	return -1;
}
void rdpolar(string fn){
	void initRestypes_pro(); initRestypes_pro();
	nkind = matype;
	for(int i=0; i<matype; i++) iskind[i] = i;
	if(Atype == 0) return;
//
	for(int i=0; i<matype_pro; i++){
		string &an = atomtypes[i]->getname();
		if(Atype == 1){
			if(an == "CA") continue;
		} else if (Atype == 12){
			if(an == "CB") continue;
		} else if (Atype == 2){
			if(an=="CA" || an=="CB") continue;
		} else if (Atype == 5){
			if(an=="CA" || an=="CB" || an=="N" || an=="O" || an=="C") continue;
		}
		iskind[i] = -1;
	}
}
void rdmol2type(string fn){
	if(fn == "") fn = datadir + "amino.mol2";
	FILE *fp = openfile(fn, "r");
	char str[121], ss[9][10];
	fgets(str, 120, fp);
	while(fgets(str, 120, fp) != NULL){
		if(strstr(str, "END") == str) break;
		if(strstr(str, "ATOM ") != str) continue;
		str2dat(str+4, 3, ss);
		int ia = atomDefine(ss[0], ss[1]);
		if(ia < 0) continue;
		int imol2 = mol2Define(ss[2]);
		ismol2[ia] = imol2;
	}
	nkind2 = matype_mol2;
	for(int i=0; i<matype_mol2; i++){
		ismol2[matype_pro+i] = i;// + matype_mol2;
	}
//	nkind2 += matype_mol2;
	if(DEBUG > 0){
		fprintf(stderr, "nkind2: %d\n", nkind2);
	}
}
void initEobs(string fn){
//	void initpobs1(string fn); initpobs1("gyr.dat");
	if(fn == ""){
		if(datadir=="" && getenv("DATADIR") != NULL){
			datadir = getenv("DATADIR");
			if(datadir != "") datadir += string("/");
		}
		fn = datadir + "dfire.2";
	}
	FILE *fp=openfile(fn, "r");
	char str[6001], *strs, *stre;
	if(DEBUG) cerr<<fn<<endl;
// read the Nobs libary
	int Nobs[nkind][nkind][mbin], Nobs0[nkind][mbin], dat[nkind];
	nkind = matype_pro;
	for(int i=0; i<mbin; i++){
		fgets(str,120,fp);
		if(strchr(str, ':') == NULL) die("Error lib: %s\n", str);
		for(int j=0; j<nkind; j++){
			fgets(str,6000,fp); int it = str2dat(str, dat);
			if(it != nkind) die("wrong lib: %d -- %d\n%s\n", it, nkind, str);
			for(int k=0; k<nkind; k++) Nobs[j][k][i] = dat[k];
		}
	}
	if(fgets(str,120,fp) != NULL){		// Nobs0
		for(int i=0; i<mbin; i++){
			fgets(str,6000,fp); int it = str2dat(str, dat);
			if(it != nkind) die("wrong Nobs0: %s", str);
			for(int k=0; k<nkind; k++) Nobs0[k][i] = dat[k];
		}
	}
//
	nkind = matype;
/* symtry
if(bsym){
	for(int i=0; i<mbin; i++){
		for(int j=0; j<nkind; j++)
		for(int k=j+1; k<nkind; k++){
			Nobs[k][j][i] = Nobs[j][k][i] + Nobs[k][j][i];
			Nobs[j][k][i] = Nobs[j][k][i];
		}
	}
}*/
// mol2
	rdmol2type("");
	int Nobs2[nkind2][nkind2][mbin]; bzero(Nobs2, sizeof(Nobs2));
	int Nobs3[nkind][nkind2][mbin]; bzero(Nobs3, sizeof(Nobs3));
	double pvol[3][mbin]; bzero(pvol, sizeof(pvol));
	for(int b=0; b<mbin; b++){
		for(int i=0; i<matype; i++)
		for(int j=i; j<matype; j++){
			int ik = ismol2[i], jk = ismol2[j];
			int dt = Nobs[i][j][b] + Nobs[j][i][b];
			if(j == i) dt /= 2;
			Nobs2[jk][ik][b] += dt; Nobs3[i][jk][b] += dt;
			Nobs2[ik][jk][b] += dt; Nobs3[j][ik][b] += dt;
		}
	}
/*	for(int i=0; i<nkind; i++)
	for(int j=0; j<nkind; j++){
		int ik = ismol2[i], jk = ismol2[j];
		for(int b=0; b<mbin; b++) Nobs[i][j][b] = Nobs2[ik][jk][b];
	}*/
// pobs0
	double pobs[mbin], pobs0[mbin];
	for(int m=0; m<nbin; m++) pobs0[m] = pow(bin2r(m), ALPHA);
	double ds = sum(nbin, pobs0);
	for(int m=0; m<nbin; m++) pobs0[m] /= ds;
// calculte Eobs
	double dat1[mbin];
	for(int i=0; i<matype_pro; i++)
	for(int j=matype_pro; j<matype; j++){
		int ik=ismol2[i], jk=ismol2[j]; int n0 = 75;
		if(etype == 2){
			ds = sum(nbin, Nobs2[ik][jk]) + n0;
			for(int m=0; m<nbin; m++) pobs[m] = (Nobs2[ik][jk][m]+n0*pobs0[m]) / ds;
		}else{
//		ds = sum(nbin, Nobs[i][j]) + n0;
//		for(int m=0; m<nbin; m++) pobs[m] = (Nobs[i][j][m]+n0*pobs[m]) / ds;
			ds = sum(nbin, Nobs3[i][jk]) + n0;
			for(int m=0; m<nbin; m++) pobs[m] = (Nobs3[i][jk][m]+n0*pobs0[m]) / ds;
		}
//
		for(int m=0; m<nbin; m++){
			double e, r = bin2r(m), rc=bin2r(nbin-1);
			e = -ert*log(pobs[m] / pobs0[m]);
//			e = -ert*log(pobs[m] / pow(r/rc, ALPHA) /  pobs[nbin-1]);
//			e = min(e, epen1*2);
			dat1[m] = e;
		}
		for(int m=0; m<nbin; m++) {
			edfire[j][i][m] = edfire[i][j][m] = min(epen1, dat1[m] - dat1[nbin-1]);
		}
	}
	if(DEBUG > 2)
	{
	for(int m=0; m<nbin; m++){
		printf("%d %f %f\n", m, edfire[2][167+1][m], edfire[20][167+7][m]);
	}
	exit(0);
	}
	return;
}
