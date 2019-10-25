#include "protein.h"
#include "array1.h"

class Molecule:public Protein{
public:
	Molecule(const string str):Protein(str){};
	void getxca(vector<Xvec> &xc);
};
void Molecule::getxca(vector<Xvec> &xc){
	xc.clear();
	for(int i=0; i<nres; i++) xc.push_back( Xvec(atoms[rstart[i]+1].getx()) );
}
double calrms1(vector<Xvec> &xa, vector<Xvec> &xb);
int main(int argc, char *argv[]){
	if(argc < 2) die("usage: RUN pdb1 pdb2"); 
	vector<Xvec> list_xca1, list_xca2;
	Molecule *mol1 = new Molecule(argv[1]);
	mol1 -> getxca(list_xca1);

	for(int m=2; m<argc; m++){
		Molecule *mol2 = new Molecule(argv[m]);
		mol2 -> getxca(list_xca2);
		assert(list_xca1.size() == list_xca2.size());
		cout<<argv[m]<<": "<<calrms1(list_xca1, list_xca2)<<endl;
	}
}
