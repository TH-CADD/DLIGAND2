#include "rotamer.h"
#include "restype.h"

Rotamers rotamers_aa[20][36][36];
void Rotamers::addRotamer(char *line){
	double dat[10];
	int nt = str2dat(line+31, dat);
	if(nt != 9) die("Error rotamer line: %d\n%s\n", nt, line);
	addRotamer(*dat, dat+1);
}
void Rotamers::addRotamer(double p, double *dim){
	int nc = 0;
	for(int m=0; m<4; m++){
		if(dim[m+4] > 1.e-3) nc++;
		else break;
	}
	chi.push_back(vector<int>(2*nc));
	int nchi = chi[0].size() / 2;
	if(nchi != nc) die("wrong rotamer: %d %d %f %f\n", nchi, nc, p, *dim);
	vector<int> &ip = chi.back();
	for(int i=0; i<nchi; i++){
		ip[i] = int(dim[i]*10.);
		ip[i+nchi] = int(dim[i+4]*10.);
	}
	ps.push_back( p );
}
void initRotamers(){
	static bool binit = false;
	if(binit) return;
	binit = true;
	string fn = getdatadir() + "bbdep02.cbin";
	if(file_existed(fn)){
		bool initRotamers_bin(string fn);
		initRotamers_bin(fn); return;
	}
//
	fn = getdatadir() + "bbdep02.May.sortlib";
	fprintf(stderr, "reading %s\n", fn.c_str());
	FILE *fp = openfile(fn, "r");
	char str[121], rn[4];
	while( fgets(str, 120, fp) != NULL){
		if(str[0] == '#') continue;
		int phi, psi;
		sscanf(str, "%s%d%d", rn, &phi, &psi);
		int id = resDefine(rn);
		if(id < 0) die("Error residue: %s\n", rn);
		if(phi==-180 || psi==-180) continue;
		if(phi < 0) phi +=  360;
		if(psi < 0) psi +=  360;
		phi /= 10; psi /= 10;
		rotamers_aa[id][phi][psi].addRotamer(str);
	}
	fclose(fp);
}
bool initRotamers_bin(string fn){
	FILE *fp = openfile(fn, "rb");
   char str[201], rn[9];
   short bsame, idat[9];
	fprintf(stderr, "reading %s\n", fn.c_str());
   fread(str, sizeof(char), 50, fp);
   int id = -1, phi=-1, psi=-1;
   double p, chi[8];
  	while(1){
		fread(&bsame, sizeof(short), 1, fp);
		if(bsame == 2) break;
		if(! bsame){
			fread(rn, sizeof(char), 3, fp); rn[3] = '\0';
			id = resDefine(rn);
			fread(idat, sizeof(short), 3, fp);
			phi = idat[0]; psi = idat[1];
			if(phi < 0) phi +=  360;
			if(psi < 0) psi +=  360;
			phi /= 10; psi /= 10;
		}
   	float p0;
		fread(&p0, sizeof(float), 1, fp);
		fread(idat, sizeof(short), 8, fp);
		p = p0;
		for(int i=0; i<8; i++) chi[i] = idat[i] / 10.;
		rotamers_aa[id][phi][psi].addRotamer(p, chi);
	}
	fread(rn, sizeof(char), 3, fp); rn[3] = '\0';
	if(string(rn) != "END") die("ERROR ROTAMER END: %s\n", rn);
	fclose(fp);
	return true;
}
