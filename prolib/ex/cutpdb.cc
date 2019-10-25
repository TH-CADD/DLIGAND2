#include "misc.h"
#include "array1.h"

vector<string> list_atoms;
vector<Xvec> xref;
void readpdb(string fn){
	FILE *fp = openfile(fn);
	char str[201];
	while(fgets(str, 200, fp) != NULL){
		if(strstr(str, "ATOM")==NULL && strstr(str, "HETATM")==NULL) continue;
		str[54] = '\0'; list_atoms.push_back(str);
	}
	fclose(fp);
}
void getcoords(char chn){
	double x[3];
	for(int i=0; i<list_atoms.size(); i++){
		string& line = list_atoms[i];
		if(line[21] != chn) continue;
		str2dat(line.substr(30, 24), 3, x);
		xref.push_back(x);
	}
}
void print_all(char chn, double rcut){
	Xvec x;
	getcoords(chn);
	for(int i=0; i<list_atoms.size(); i++){
		bool bclose = 0;
		string& line = list_atoms[i];
		if(line[21] == chn) {printf("%s\n", line.c_str()); continue;}
		str2dat(line.substr(30, 24), 3, x.getdat());
		for(int k=0; k<xref.size(); k++){
			if(xref[k].distance2(x) < rcut*rcut) {bclose = 1; break;}
		}
		if(bclose) printf("%s\n", line.c_str());
	}
}
int main(int argc, char *argv[]){
	if(argc < 2) die("uage: RUN pdb chain cut");
	list_atoms.clear();
	xref.clear();
	readpdb(argv[1]);
	print_all(argv[2][0], atof(argv[3]));
}
