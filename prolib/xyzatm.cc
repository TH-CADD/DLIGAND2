#include "protein.h"

// att! in this subroutine, xab is vector b-->a
bool xyzatm(double *xi, double *xa, double *xb, double *xc, double *ag, int chiral){
	double bond=ag[0], angle1=ag[1], angle2=ag[2];
	const double eps=1.0e-99; // RADIAN=180.0/acos(-1.);
	double sin1, cos1, sin2, cos2, rab, rbc, rac, xab[3], xbc[3], xac[3], xt[3], xu[3];
	double sinb, cosb, sing, cosg, cosine, sine, sine2, xtmp,ztmp,a,b,c;
	int m, i3;
//
	sin1 = sin(angle1/RADIAN); cos1=cos(angle1/RADIAN);
	sin2 = sin(angle2/RADIAN); cos2=cos(angle2/RADIAN);
	if(xa == NULL){
		for(m=0; m<3; m++) xi[m]=0.;
	} else if(xb == NULL){
		for(m=0; m<3; m++) xi[m]=xa[m];
		xi[2] += bond;
//	  if no fourth site given, place the atom in the x,z-plane
	} else if(xc == NULL){
		xdiff(xa, xb, xab);
		rab = sqrt(dot_product(xab, xab));
		for(m=0; m<3; m++) xab[m] /= rab;
		cosb = xab[2];
		sinb = sqrt(xab[0]*xab[0] + xab[1]*xab[1]);
		if(sinb < eps){
			cosg=1.0; sing=0.0; 
		} else{
			cosg = xab[1] / sinb; sing = xab[0] / sinb;
		}
		xtmp = bond*sin1;
		ztmp = rab - bond*cos1;
		xi[0] = xb[0] + xtmp*cosg + ztmp*sing*sinb;
		xi[1] = xb[1] - xtmp*sing + ztmp*cosg*sinb;
		xi[2] = xb[2] + ztmp*cosb;
//	  general case where the second angle is a dihedral angle
	} else if(chiral == 0){
		xdiff(xa, xb, xab);
		xdiff(xb, xc, xbc);
		rab = sqrt(dot_product(xab, xab));
		rbc = sqrt(dot_product(xbc, xbc));
		for(m=0; m<3; m++){
			xab[m] /= rab; xbc[m] /= rbc;
		}
		cross_product(xbc, xab, xt);
		cosine = dot_product(xab, xbc);
		sine = sqrt(max(1.0-cosine*cosine, eps));
		if (fabs(cosine) > 1.0){
			fprintf(stderr, "XYZATM -- Undefined Dihedral, Angle at Atom\n");
			return false;
		}
		for(m=0; m<3; m++) xt[m] /= sine;
		cross_product(xt, xab, xu);

		for(m=0; m<3; m++) xi[m] = xa[m] + bond*(xu[m]*sin1*cos2 + xt[m]*sin1*sin2 - xab[m]*cos1);
//	  general case where the second angle is a bond angle
	}else if((int)fabs(chiral) == 1){
		xdiff(xa, xb, xab);
		xdiff(xa, xc, xac);
		rab = sqrt(dot_product(xab, xab));
		rac = sqrt(dot_product(xac, xac));
		for(m=0; m<3; m++){
			xab[m] /= -rab; xac[m] /= rac;		// here xab is xba
		}
		cross_product(xac, xab, xt);
		cosine = dot_product(xab, xac);
		sine2 = max(1.0-cosine*cosine, eps);
		if (fabs(cosine) > 1.0){
			fprintf(stderr, "XYZATM --Defining Atoms Colinear at Atom\n");
			return false;
		}
		a = (-cos2 - cosine*cos1) / sine2;
		b = (cos1 + cosine*cos2) / sine2;
		c = (1.0 + a*cos2 - b*cos1) / sine2;
		if (c > eps) c = chiral * sqrt(c);
		else if (c < -eps){
			c = .0;
			for(m=0; m<3; m++){
				xtmp = (a*xac[m] + b* xab[m]);
				c += xtmp*xtmp;
			}
			c = sqrt(c); a = a / c; b = b / c; c = 0.0;
		} else c = 0.0;
		for(m=0; m<3; m++) xi[m] = xa[m] + bond*(a*xac[m] + b*xab[m] + c*xt[m]);
	}
	return true;
}
void xdiff(double *xa, double *xb, double *xba){
	for(int m=0; m<3; m++) xba[m] = xa[m] - xb[m];
}
double distance2(double *xa, double *xb){
	double xab[3];
	xdiff(xb, xa, xab);
	return dot_product(xab, xab);
}
double angle(double *xa, double *xb, double *xc){
	double xba[3], xbc[3], val;
	xdiff(xa, xb, xba);
	xdiff(xc, xb, xbc);
	val = dot_product(xba, xba)*dot_product(xbc, xbc);
	if(val < 1.0e-99) return 0.;
	val = dot_product(xba, xbc)/sqrt(val);
	return acos(val)*RADIAN;
}
double torsion(double *xa, double *xb, double *xc, double *xd){
	double xab[3], xbc[3], xcd[3], x1[3], x2[3], val;
	xdiff(xb, xa, xab);
	xdiff(xc, xb, xbc);
	xdiff(xd, xc, xcd);
	cross_product(xab, xbc, x1);
	cross_product(xbc, xcd, x2);
	val = dot_product(x1, x1)*dot_product(x2, x2);
	if(val < 1.0e-99) return 0.;
	val = acos(dot_product(x1, x2)/sqrt(val));
	if(dot_product(xab, x2)<0) val = -val;
	return val*RADIAN;
}
