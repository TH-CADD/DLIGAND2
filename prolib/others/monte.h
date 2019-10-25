#ifndef _MONTE
#define _MONTE
#include "misc.h"

extern long istep;
class Monte{
	long nstep, nprint;
	int bSA;		// 0: fixed T0; 1: SA of T0 --> T1; 2: accept ratio 0.2~0.6
	double T0, T1;
	double (*calEnergy)();
	void (*do_mutation)();
	void (*save_para)(double);
	void (*retrieve_para)();
public:
	Monte(long n, double (*s1)(), void (*mut)(), void (*save)(double), void (*ret)()){
		nstep = n; calEnergy = s1; do_mutation = mut;
		save_para = save; retrieve_para = ret;
		bSA = 0; T0 = T1 = 1.;
		nprint = 10000;
	}
	void run1(){
		cpu_time();
		double T = T0, A0 = 1, eold, enew, emin;
		if(bSA == 1) A0 = exp(log(T1/T0) / nstep);
		emin = eold = (*calEnergy)();
		save_para(emin);
		printf("Starting energy: %g\n", eold);
//
		long naccept=0; double pacc0=-1.0;
//
		for(istep=1; istep<=nstep; istep++){
			(*do_mutation)();
			enew = (*calEnergy)();
			bool bsucc = 1;
			if(enew > eold){
				double dt = exp(-(enew-eold) / T);
				if(dt < drand48()) bsucc = 0;
			}
			if(bsucc) {
				eold = enew;
				(*save_para)(enew);
				naccept ++;
				if(emin > eold) emin = eold;
			} else {
				(*retrieve_para)();
			}
			if(bSA == 1) T *= A0;
			int nrec = 1e3;
			if(istep%nrec == 0){
				double p1 = naccept / double(nrec);
				if(pacc0 < 0) pacc0 = p1;
				pacc0 += 0.05 * (p1 - pacc0);
				if(bSA == 2){
					if(pacc0 > 0.6) T /= 1.1;
					else if(pacc0 < 0.2) T *= 1.1;
				}
				naccept = 0;
			}
			if(istep%nprint==0 || istep==nstep){
				printf("%ld Time %.1fm: %.3f %g %g %g %g\n", istep, cpu_time(), pacc0, enew, eold, emin, T);
			}
		}
		printf("Lowest energy: %g\n", emin);
	}
	friend void runMC();
};
#endif
