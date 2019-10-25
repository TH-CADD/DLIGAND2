#include "protein.h"

Restpl *restpl = NULL;
Restpl::Restpl(string fn){
	if(fn == "") fn = getdatadir() + "template.pdb";
	FILE *fp = openfile(fn, "r");
	char str[121], an[5]; string line, s0, s1;
	nres = natm = 0;
	while(fgets(str, 120, fp) != NULL){
		if(strstr(str, "ATOM ") != str) continue;
		line = str;
		s1 = line.substr(17,10);
		if(s0 != s1) {
			sscanf(str+17, "%3s", an);
			an[3]='\0'; rnam[nres] = an;
			rstart[nres] = natm;
			nres ++;  s0 = s1;
			if(nres >= mres) die("Increase mres in protpl!\n");
		}
		sscanf(str+12, "%s", an); sscanf(str+30, "%lf%lf%lf", x[natm], &x[natm][1], &x[natm][2]);
		anam[natm] = an; natm ++;
		if(natm >= matm) die("Increase matm in protpl!\n");
	}
	rstart[nres] = natm;
}
int Restpl::findatom_tpl(string rn, string an){
	int ires=-1;
	for(int i=0; i<nres; i++){
		if(rn != rnam[i]) continue;
		ires = i; break;
	}
	if(ires == -1) {
		fprintf(stderr, "Fail to find residue in template: %s %d\n", rn.c_str(), nres);
		exit(1);
	}
//
	for(int i=rstart[ires]; i<rstart[ires+1]; i++){
		if(an == anam[i]) return i;
	}
	fprintf(stderr, "Fail to find atom in template: %s %s\n", rn.c_str(), an.c_str());
	exit(1);
}
void Restpl::calint_tpl(int *idx, int &chiral, double *ag){
	int i0=idx[0], i1=idx[1], i2=idx[2], i3=idx[3];
	if(i1 < 0) return;
	ag[0] = sqrt( distance2(x[i0], x[i1]) );
	if(i2 < 0) return;
	ag[1] = angle(x[i0], x[i1], x[i2]);
	if(i3 < 0) return;
	ag[2] = torsion(x[i0], x[i1], x[i2], x[i3]);
	if(chiral != 0) {
		chiral = 1;
		if( ag[2] > 0) chiral = -1;
		ag[2] = angle(x[i0], x[i1], x[i3]);
	}
}
void initRestpl(string fn){
	if(restpl != NULL) return;
	restpl = new Restpl(fn);
}
