#ifndef _RESTYPE
#define _RESTYPE
#include "aastd.h"
#include "misc.h"
using namespace std;

class Atomtype{
	int type;
	string name;
public:
	Atomtype(int ia, const string an){type=ia; name=an;};
	bool isname(const string an);
	string &getname(){return name;};
	int gettype(){return type;};
};
//
class Restype{
	int type;
	string name;
	vector <string> neibs[2];
	vector <Atomtype> atypes;
public:
	Restype(){};
	Restype(int t, const string s){type=t; name=s;}; 
	void setRestype(const string s){name=s;}; 
	void addatom(int type, const string s){atypes.push_back(Atomtype(type,s));};
	void addneib(const string s, const string s2);
	void addatom(const string s);
	void addneib_match(const string s, const string s2);
	string &getname(){return name;};
	vector <Atomtype>& getatypes(){return atypes;};
	vector <string> *getneibs(){return neibs;};
	void prtneib();
	int getnatm(){return atypes.size();};
	int define(string anam, int &pos);
};
extern vector <Restype> restypes;
extern vector<Atomtype*> atomtypes;
extern int natype;
//
void initRestypes_pro();
void initRestypes(string);
int atomDefine(string rn, string an);
int resDefine(string);
string getdatadir();

#endif
