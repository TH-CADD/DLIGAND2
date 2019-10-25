#include "molecule.h"

double Molecule::score(){
	double etot=0;
	int Nb[mbin]; bzero(Nb, sizeof(Nb));

        for (int ka=0;ka<167;ka++) {
            for (int kb=167;kb<178;kb++) {
            printf("|%d |%d |",ka,kb);
            prtdim(30,edfire[ka][kb]);
           // printf("\n");
            }
        }
        exit(0);

	for(int i=0; i<natm; i++){
		double e0 = 0., r2min=1000;
		int ta = atoms[i].gettype(); if(ta < 0) continue;
		if(ta >= matype_pro) continue;
		int ka = iskind[ta]; if(ka < 0) continue;
		for(int j=0; j<natm; j++){
			if(rseq[i] == rseq[j]) continue;
			int tb = atoms[j].gettype(); if(tb < 0) continue;
			if(tb < matype_pro) continue;
			int kb = iskind[tb]; if(kb < 0) continue;
			double r2 = distance2(i, j);
			r2min = min(r2, r2min);
//			if(r2 > rcut*rcut) continue;
			if(r2 > Rcut*Rcut) continue;
			//printf("Rcut is: %f",rcut);
			//if(r2 > 11.75*11.75) continue;
			int b = r2bin(sqrt(r2));
//			printf("%d %d %d %f\n", ka, kb, b, edfire[ka][kb][b]);
			Nb[b] ++;
			e0 += double(edfire[ka][kb][b]);
		}
		if(r2min > rcut*rcut) continue;
//		e0 = min(e0, epen1);
		etot += e0;
	}
//	prtdim(mbin, Nb);
	return etot;
}
void Molecule::stats(){
	for(int i=0; i<natm; i++){
		int ta = atoms[i].gettype(); if(ta < 0) continue;
		int ka = iskind[ta]; if(ka < 0) continue;
		for(int j=i+1; j<natm; j++){
			int tb = atoms[j].gettype(); if(tb < 0) continue;
			int kb = iskind[tb]; if(kb < 0) continue;
			double r2 = distance2(i, j);
			if(r2 > Rcut*Rcut) continue;
			int b = r2bin(sqrt(r2));
			if(rseq[i] == rseq[j]) {Nobs0[ka][b] ++; Nobs0[kb][b] ++; continue;}
			edfire[ka][kb][b] ++;
		}
	}
}
void Molecule::rdmol2(string fn){
	string ligname = fn;
	int mol2Define(string s1);
	FILE *fp = openfile(fn, "r");
	char str[201], s1[20];
	while(fgets(str, 120, fp) != NULL){
		if(strstr(str, "@<TRIPOS>ATOM") != NULL) break;
	}
	string s0;
	while(fgets(str, 200, fp) != NULL){
		if(strstr(str, "@<TRIPOS>BOND") != NULL) break;
		sscanf(str+59, "%s", s1);
		strcpy(s1, "LIG");
		if(s0 != s1){
			s0 = s1;
			rnam.push_back(s1);
			resid.push_back(-1);
			chnseq[nres] = nchain;
			if(rstart.size() == nres) rstart.push_back(natm);
			nres ++;
		}
		sscanf(str+47, "%s", s1);
		string an = s1;
		if(an[0]=='H' || an=="LP") continue;
//		if(an=="F" || an=="Cl" || an=="Met") continue;
		int it = mol2Define(s1);
		atoms.push_back( Atom(natm, matype_pro+it, s1) );
		Atom *ap = &atoms.back();
		ap -> setx(str+15);
		rseq.push_back(nres-1);
		natm ++;
	}
	fclose(fp);
	rstart.push_back(natm);
	nchain ++;
}
