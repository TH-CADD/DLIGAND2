#include "aastd.h"

namespace PBB_dat{
	bool binit_PBB = false;
	double pbb0[20][36][36], pbb_sum[20][36*36];
};
void initpbb_aa(){	// read the probability of amini acid on backbone torsion
	using namespace PBB_dat;
	extern string getdatadir();
	FILE *fp = openfile(getdatadir()+"phipsi_ra1.lib", "r");
	int id = 0; char str[201];
	while(fgets(str, 200, fp) != NULL){
		if(*str == '#') continue;
		str[3] = '\0';
		id = aaDefine(str);
		if(id < 0) die("wrong resi: %s\n", str);
		for(int m=0; m<36; m++){
			fgets(str, 200, fp);
			int na = str2dat(str, pbb0[id][m]);
			if(na != 36) die("wrong nele: %d != 36\n%s\n", na, str);
		}
	}
	fclose(fp);
	const int mtot = 36*36;
	for(int i=0; i<20; i++){
//		double ds = sum(mtot, &pbb0[i][0][0]);
		int it = 0; double ds = 0.;
		for(int j=0; j<36; j++)
		for(int k=0; k<36; k++) {
			ds += pbb0[i][j][k];
			pbb_sum[i][it++] = ds;
		}
		assert(ds > 0);
		it = 0;
		for(int j=0; j<36; j++)
		for(int k=0; k<36; k++) {
			pbb0[i][j][k] /= ds;
			pbb_sum[i][it++] /= ds;
		}
		assert( abs(pbb_sum[i][mtot-1]-1.) < 1.0e-3 );
		pbb_sum[i][mtot-1] = 1.;
	}
	binit_PBB = true;
}
double getpbb_aa(int id, int iphi, int ipsi){
	using namespace PBB_dat;
	if(! binit_PBB) initpbb_aa();
	return pbb0[id][iphi][ipsi];
}
double getpbb_aa(int id, double phi, double psi){
	int iphi = int(fmod(phi+365., 360) / 10.);
	int ipsi = int(fmod(psi+365., 360) / 10.);
	return getpbb_aa(id, iphi, ipsi);
}
void getphipsi_bin(int id, int *iphipsi){
	using namespace PBB_dat;
	if(! binit_PBB) initpbb_aa();
	int it = findElement(36*36, pbb_sum[id], drand48(), "phipsi");
	iphipsi[0] = it / 36;
	iphipsi[1] = it - iphipsi[0]*36;
}
