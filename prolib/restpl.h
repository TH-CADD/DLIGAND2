#ifndef _PRO_TPL
#define _PRO_TPL
#include "misc.h"
class Restpl{
	static const int mres=80, matm=900;
	int nres, natm, rstart[mres];
	double x[matm][3];
	string anam[matm], rnam[mres];
public:
	Restpl(string);
	int findatom_tpl(string rn, string an);
	void calint_tpl(int *idx, int &chiral, double *ag);
};
extern Restpl *restpl;
void initRestpl(string);
#endif
