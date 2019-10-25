#include "misc.h"
#include "aastd.h"
#include "array1.h"

class Protein2{
protected:
	string pdbnm;
	int nres, natm;
	vector<string> rnam, resinfo, anam;	// 1,2: [nres]; 3: [natm]
	vector<int> rseq, rstart;		// [natm], [nres]
	vector<int> rstart_ch;
	vector<Xvec> x;		// [natm]
public:
	Protein2(string fn);
//	double distance2(int i, int j){return distance2(x[i], x[j]);}
	void rdpdb(string);
	void rdpdb(FILE*);
	void wrpdb(string fn);
	void wrpdb(FILE*);
	int getnres(){return nres;}
	int getnatm(){return natm;}
	double torsion(int ia, int ib, int ic, int id);
	void calzcrd();
//
	friend class Molecule;
};

Protein2::Protein2(const string fn){
	pdbnm = fn;
	natm = nres = 0;
	rstart_ch.push_back(0);
	rdpdb(fn);
}
void Protein2::rdpdb(string fn){
	FILE *fp=openfile(fn, "r");
	rdpdb(fp); fclose(fp);
}
void Protein2::rdpdb(FILE *fp){
	char str[121], an[5];
	string line, resinfo0="";
	double xt[3];
	while(fgets(str,120,fp) != NULL){
		line = str;
		if(line.substr(0,3) == "END") break;
		if(line.substr(0,4) != "ATOM") continue;
		string info = line.substr(17, 10);
		sscanf(str+12, "%4s", an); an[3] = '\0';
		if(an[0] == 'H') continue;
		if(nres == 0 || info != resinfo0) {
			resinfo.push_back(info);
			resinfo0 = info;
			rnam.push_back(info.substr(0,3));
			if(rstart.size() <= nres) rstart.push_back(natm); // false when reading 2nd chain
			nres ++;
		}
		anam.push_back(an);
		rseq.push_back(nres-1);
		sscanf(str+30, "%lf%lf%lf", xt, xt+1, xt+2);
		x.push_back( Xvec(xt) );
		natm ++;
	}
	rstart.push_back(natm);
	rstart_ch.push_back(natm);
}
void Protein2::wrpdb(string fn){
	FILE *fp = openfile(fn, "w");
	wrpdb(fp); fclose(fp);
}
void Protein2::wrpdb(FILE *fp){
	char fmt1[] = "ATOM%7d  %-4s%3s %c%4d    %8.3f%8.3f%8.3f\n";
	for(int ii=0; ii<rstart_ch.size(); ii++){
		char chn = 'A' + ii;
		for(int i=rstart_ch[ii]; i<rstart_ch[ii+1]; i++){
			int ir = rseq[i];
			double *xt = x[i].getx();
			fprintf(fp, fmt1, i+1, anam.at(i).c_str(), rnam.at(ir).c_str(), chn,
					ir+1, *xt, xt[1], xt[2]);
		}
		fprintf(fp, "TER\n");
	}
	fprintf(fp, "END\n");
}
void xdiff(double *xa, double *xb, double *xba){
	for(int m=0; m<3; m++) xba[m] = xa[m] - xb[m];
}
double Protein2::torsion(int ia, int ib, int ic, int id){
	double xab[3], xbc[3], xcd[3], x1[3], x2[3], val;
	double *xa=x[ia].getx(), *xb=x[ib].getx(), *xc=x[ic].getx(), *xd=x[id].getx();
	xdiff(xb, xa, xab);
	xdiff(xc, xb, xbc);
	xdiff(xd, xc, xcd);
	cross_product(xab, xbc, x1);
	cross_product(xbc, xcd, x2);
	val = dot_product(x1, x1)*dot_product(x2, x2);
	if(val < 1.0e-99) return 0.;
	val = acos(dot_product(x1, x2)/sqrt(val));
	if(dot_product(xab, x2)<0) val = -val;
	return val*RADIAN;
}
void Protein2::calzcrd(){
	printf("# rn Phi Psi Omega\n");
	for(int ir=0; ir<nres; ir++){
		double ag[3];
		for(int m=0; m<3; m++) ag[m] = 360;
		int i1 = rstart[ir];
		if(ir > 0) {
			int i0 = rstart[ir-1];
			ag[0] = torsion(i0+2, i1, i1+1, i1+2);
			ag[2] = torsion(i0+1, i0+2, i1, i1+1);
		}
		if(ir < nres-1) {
			int i2 = rstart[ir+1];
			ag[1] = torsion(i1, i1+1, i1+2, i2);
		}
		printf("%d %s %f %f %f\n", ir+1, rnam[ir].c_str(), ag[0], ag[1], ag[2]);
	}
}
int main(int argc, char *argv[]){
	for(int i=1; i<argc; i++){
		Protein2 *mol=new Protein2(argv[i]);
		cout<<argv[i]<<endl;
		mol -> calzcrd();
		delete mol;
	}
}
