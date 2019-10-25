#ifndef _ROTAMER
#define _ROTAMER

#include "protein.h"

class Rotamers{
	vector<float> ps;
//	int nchi;
	vector< vector<int> > chi;			// chi[0] [ag0,ag1.., dag0,dag1..]
public:
	Rotamers(){}
	~Rotamers(){}
	void addRotamer(char*);
	void addRotamer(double, double*);
//
	int getnrot(){return chi.size();}
	int getnchi(){if(chi.size()==0) return 0; return chi[0].size()/2;}
	vector<int> &getchi(int k, double &p){p = ps.at(k); return chi.at(k);}
	vector<int> &getchi(int k){return chi.at(k);}
	float &getp(int i){return ps.at(i);}
};
extern Rotamers rotamers_aa[20][36][36];
void initRotamers();

#endif
