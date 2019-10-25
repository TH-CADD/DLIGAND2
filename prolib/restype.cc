#include "restype.h"

vector <Restype> restypes;
vector<Atomtype*> atomtypes;
int natype = 0;
int atomDefine(string rn, string an){
	int rid = resDefine(rn), pos;
	if(rid < 0) return -1;
	int type = restypes.at(rid).define(an, pos);
	if(DEBUG>0 && type<0) fprintf(stderr, "Unkown atom: %s %s\n", rn.c_str(), an.c_str());
	return type;
}
void collect_atomtypes(){
	atomtypes.clear();
	for(unsigned i=0; i<restypes.size(); i++){
		vector <Atomtype>& ap = restypes[i].getatypes();
		for(unsigned m=0; m<ap.size(); m++){
			atomtypes.push_back( &ap.at(m) );
		}
	}
}
void initRestypes_pro(){
	int flag = 0;
	restypes.clear(); natype = 0;
	for(int i=0; i<nres_std; i++){
		restypes.push_back(Restype(i, rnam3_std[i]));
		string rn(rnam3_std[i]);

		Restype *res = &restypes[i];
		res->addatom("N"); res->addatom("CA"); res->addatom("C"); res->addatom("O");
//		if(flag==2){res->addatom("CB"); continue;}				// only main+Cb
//		else if(flag==1 && rn!="PRO") res->addatom("H");		// mainchain H included
//
		if(rn != "GLY") res->addatom("CB");
		if(flag!=0 && rn!="GLY" && rn!="ALA") res->addatom("HB1");
//
		if(rn == "VAL"){
			res->addatom("CG1"); res->addatom("CG2");
		} else if(rn == "LEU"){
			res->addatom("CG"); res->addatom("CD1"); res->addatom("CD2"); 
		} else if(rn == "ILE"){
			res->addatom("CG1"); res->addatom("CG2"); res->addatom("CD1");
		} else if(rn == "SER"){
			res->addatom("OG");
		} else if(rn == "THR"){
			res->addatom("OG1"); res->addatom("CG2");
		} else if(rn == "CYS"){
			res->addatom("SG");
		}else if(rn == "PRO"){
			res->addatom("CG"); res->addatom("CD");
			res->addneib_match("CD", "N");
		}else if(rn == "PHE"){
			res->addatom("CG"); res->addatom("CD1"); res->addatom("CD2"); res->addatom("CE1");
			res->addatom("CE2"); res->addatom("CZ");
		}else if(rn == "TYR"){
			res->addatom("CG"); res->addatom("CD1"); res->addatom("CD2"); res->addatom("CE1");
			res->addatom("CE2"); res->addatom("CZ"); res->addatom("OH"); 
		}else if(rn == "TRP"){
			res->addatom("CG"); res->addatom("CD1"); res->addatom("CD2"); res->addatom("NE1");
			res->addatom("CE2"); res->addatom("CE3"); res->addatom("CZ2"); res->addatom("CZ3");
			res->addatom("CH2");
			res->addneib_match("CE3", "CD2");
		}else if(rn == "HIS"){
			res->addatom("CG"); res->addatom("ND1"); res->addatom("CD2"); res->addatom("CE1");
			res->addatom("NE2");
			res->addneib_match("NE2", "CE1");
		}else if(rn == "ASP"){
			res->addatom("CG"); res->addatom("OD1"); res->addatom("OD2");
		}else if(rn == "ASN"){
			res->addatom("CG"); res->addatom("OD1"); res->addatom("ND2");
		}else if(rn == "GLU"){
			res->addatom("CG"); res->addatom("CD"); res->addatom("OE1"); res->addatom("OE2");
		}else if(rn == "GLN"){
			res->addatom("CG"); res->addatom("CD"); res->addatom("OE1"); res->addatom("NE2");
		}else if(rn == "MET"){
			res->addatom("CG"); res->addatom("SD"); res->addatom("CE");
		}else if(rn == "LYS"){
			res->addatom("CG"); res->addatom("CD"); res->addatom("CE"); res->addatom("NZ");
		}else if(rn == "ARG"){
			res->addatom("CG"); res->addatom("CD"); res->addatom("NE"); res->addatom("CZ");
			res->addatom("NH1"); res->addatom("NH2");
		}
	}
	collect_atomtypes();
}
void Restype::addatom(const string an){
	atypes.push_back(Atomtype(natype++, an));
	if(an=="CA") addneib_match(an,"N");
	else if(an=="C") addneib_match(an, "CA");
	else if(an=="O") addneib_match(an, "C");
	else if(an=="H") addneib_match(an, "N");
	else if(an=="CB") addneib_match(an, "CA");
	else if(an=="HB1") addneib_match(an, "CB");
	else if(an.size()>1){
//
		string str;
		if(an.at(1) == 'G') str="*B";
		else if(an.at(1) == 'D') str="*G";
		else if(an.at(1) == 'E') str="*D";
		else if(an.at(1) == 'Z') str="*E";
		else if(an.at(1) == 'H') str="*Z";
		if(an.size()>=3) str += an.substr(2);
		addneib_match(an, str);
	}
}
void Restype::addneib(const string an1, const string an2){
	int i1, i2;
	i1 = i2 = -1;
	for(unsigned i=0; i<atypes.size(); i++){
		if(atypes.at(i).getname() == an1) i1 = i;
		else if(atypes.at(i).getname() == an2) i2 = i;
	}
	if(i1<0 || i2<0) {
		if(DEBUG > 0){
			fprintf(stderr, "not found pair in res %s: %s %s %d %d\n", name.c_str(), an1.c_str(), an2.c_str(), i1, i2);
		}
		return;
	}
	neibs[0].push_back(an1);
	neibs[1].push_back(an2);
}
void Restype::addneib_match(const string an1, const string an2){
	int i1, i2;
	vector <int> ilist;
	i1 = i2 = -1;
	for(unsigned i=0; i<atypes.size(); i++){
		if(atypes.at(i).isname(an2)) ilist.push_back(i);
		else if(atypes.at(i).isname(an1)) {i1=i; break;}
	}
	i2 = ilist.size();
	if(i1<0 || i2<=0) {
		if(DEBUG > 0){
			fprintf(stderr, "not found pair in res %s: %d %d %s %s\n", name.c_str(), i1, i2, an1.c_str(), an2.c_str());
		}
/*		for(int i=0; i<atypes.size(); i++) cerr<<*atypes.at(i).getname()<<' ';
		cerr<<"Error: "<<name<<' '<<an1<<' '<<an2<<endl;*/
		return;
	}
	for(unsigned i=0; i<ilist.size(); i++){
		i2 = ilist.at(i);
		neibs[0].push_back(atypes.at(i1).getname());
		neibs[1].push_back(atypes.at(i2).getname());
	}
}
int Restype::define(string an, int &pos){
	for(int i=0; i<atypes.size(); i++){
		if(atypes.at(i).getname()==an) {
			pos = i;
			return atypes.at(i).gettype();
		}
	}
	return -1;
}
bool Atomtype::isname(const string an){
	int is, ie, flag;
	if(an[0]=='*'){
		is=1; ie = min(name.size(), an.size()) - is;
		ie = max(ie, 1);
		flag = name.compare(is, ie, an, is, ie);
	} else flag = name.compare(an);
	return (flag == 0);
}
int resDefine(string rn){
	if(restypes.size() < 1) return aaDefine(rn);
	for(int i=0; i<restypes.size(); i++){
		if(restypes[i].getname() == rn) return i;
	}
	if(DEBUG > 0){
		fprintf(stderr, "unknown residue name: %s\n", rn.c_str());
	}
	return -1;
}
void initRestypes(string fn){
	if(fn == "") fn = getdatadir() + "amino.dat";
	if(DEBUG > 0) fprintf(stderr, "reading file: %s\n", fn.c_str());
	FILE *fp = openfile(fn, "r");
	char str[201], ss[15][20]; 
	string rn0="", an0="";
	restypes.clear();
	int nr = 0, na = 0; Restype *rp = NULL;
	while(fgets(str, 200, fp) != NULL){
		if(strstr(str, "END") == str) break;
		if(strstr(str, "ATOM") != str) continue;
		int nt = str2dat(str, ss);
		if(rn0 != ss[1]){
			restypes.push_back(Restype(nr, ss[1]));
			rp = &restypes[nr];
			nr ++; rn0 = ss[1];
		}
		rp -> addatom(na, ss[2]);
		int nn = strtol(ss[5], NULL, 10);
		if(nn == 1) rp -> addneib(ss[2], an0);
		else if(nn < 0){
			for(int m=0; m<-nn; m++) rp -> addneib(ss[2], ss[6+m]);
		}
		na ++; an0 = ss[2];
	}
	natype = na;
	fclose(fp);
	if(DEBUG){
		fprintf(stderr, "Number of Restype: %d; Atomtype: %d\n", nr, na);
	}
	collect_atomtypes();
}
string getdatadir(){
	static bool binit = false;
	if(binit) return datadir;
	binit = true;
	if(getenv("DATADIR") != NULL){
		datadir = getenv("DATADIR");
		if(datadir!="" && datadir.at(datadir.size()-1)!='/') datadir += '/';
	}
	return datadir;
}
