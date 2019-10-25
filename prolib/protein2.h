#ifndef _PROTEIN2
#define _PROTEIN2
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
//
	friend class Molecule;
};
#endif
