#include "protein.h"

class Molecule:public Protein{
public:
	Molecule(const char *str):Protein(str){};
	void test1(int,double);
};
void calRmat_axis(double *xv, double ag, double R[][3]);
void calXrot(double R[][3], double *xa, double *x0, double *xn);
void Molecule::test1(int ir0, double ag){
	double xd[3], R[3][3], x2[3];
	int i0 = rstart[ir0];
	xdiff(i0, i0+1, xd);
	normalize(xd);
	calRmat_axis(xd, ag, R);
//	calXrot(R, atoms[i0].getx(), atoms[natm-1].getx(), x2);
/*	double distance2(double*, double*);
	cout<<ir<<' '<<sqrt(distance2(atoms[natm-1].getx(), x2))<<endl;

	continue;*/
	for(int i=i0+2; i<natm; i++){
		calXrot(R, atoms[i0].getx(), atoms[i].getx(), atoms[i].getx());
	}
}
int main(int argc, char *argv[]){
	if(argc < 2) die("Usage: %s PDBs ir ag\n", argv[0]);
	Molecule *mol = new Molecule(argv[1]);
	int ir0 = atoi(argv[2]);
	double ag = atof(argv[3]);
	mol -> test1(ir0, ag/RADIAN);
	mol -> wrpdb(stdout);
}
