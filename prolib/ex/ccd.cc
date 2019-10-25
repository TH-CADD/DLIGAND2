#include "misc.h"

static const double radian = 180./acos(-1);
void calRmat_fit(int nfit, double *w, double xref[][3], double xmv[][3], double U[][3]){
	double xm1[3], xm2[3];
	bzero(xm1, sizeof(xm1)); bzero(xm2, sizeof(xm2));
	double x1[nfit][3], x2[nfit][3];
	for(int i=0; i<nfit; i++)
	for(int m=0; m<3; m++){
		xm1[m] += xref[i][m]; xm2[m] += xmv[i][m];
	}
	for(int m=0; m<3; m++) {xm1[m] /= nfit; xm2[m] /= nfit;}
	for(int i=0; i<nfit; i++)
	for(int m=0; m<3; m++){
		x1[i][m] = xref[i][m] - xm1[m];
		x2[i][m] = xmv[i][m] - xm2[m];
	}
//
	void quatfit(int n, double *w, double x1[][3], double x2[][3], double u[][3]);
	quatfit(nfit, w, x1, x2, U);
	for(int m=0; m<3; m++) U[3][m] = xm1[m] - dot_product(U[m], xm2);
}
void calRmat_axis(double *xv, double ag, double R[][3]){
	double v1, v2, v3;
	double t1, t2, t3, t6, t7, t8, t9, t11, t12, t15, t19, t20, t24;
	v1 = xv[0]; v2 = xv[1]; v3 = xv[2];
	t1 =  cos(ag); t2 =  1 - t1; t3 =  v1*v1;
	t6 =  t2*v1; t7 =  t6*v2; t8 =  sin(ag); t9 =  t8*v3;
	t11 = t6*v3; t12 = t8*v2; t15 = v2*v2;
	t19 = t2*v2*v3; t20 = t8*v1; t24 = v3*v3;
	R[0][0] = t1 + t2*t3;
	R[1][0] = t7 - t9;
	R[2][0] = t11 + t12;
	R[0][1] = t7 + t9;
	R[1][1] = t1 + t2*t15;
	R[2][1] = t19 - t20;
	R[0][2] = t11 - t12;
	R[1][2] = t19 + t20;
	R[2][2] = t1 + t2*t24;
}
void calXrot(double R[][3], double *xa, double *x0, double *xn){
	double xd[3];
	for(int m=0; m<3; m++) xd[m] = x0[m] - xa[m];
	for(int m=0; m<3; m++){
		xn[m] = dot_product(R[m], xd) + xa[m];
	}
}
int ccd1(double *axis, double *x0, int nfix, double xfix[][3], double xmv[][3], double *val){
	double vf[3], vm[3], vx[3], vy[3];	// vx, r; vy, s; vz, theta
	for(int m=0; m<3; m++) val[m] = 0.;
	for(int i=0; i<nfix; i++){
		for(int m=0; m<3; m++){
			vm[m] = xmv[i][m] - x0[m];
			vf[m] = xfix[i][m] - x0[m];
		}
		double rz = dot_product(axis, vm);
		for(int m=0; m<3; m++){
			vx[m] = vm[m] - rz * axis[m];
			vf[m] -= rz * axis[m];
		}
		double rx = normalize(vx);
		cross_product(axis, vx, vy);
		val[0] += rx*rx + dot_product(vf, vf);
		val[1] += 2.*rx * dot_product(vf, vx);
		val[2] += 2.*rx * dot_product(vf, vy);
	}
	double rx = sqrt(val[1]*val[1] + val[2]*val[2]);
//	cout<<val[1]<<' '<<rx<<endl;
	double ag = 0.;
	if(rx > 1.0e-6) ag = acos(val[1] / rx); // * radian;
	if(val[2] > 0) ag = -ag;
	val[1] = rx; val[2] = ag;
// result = val[0] - val[1] * cos(A0-val[2])
	return 0;
}
/*int main(){
	double x[9][3], axis[3], x0[3], val[3], R[3][3];
	bzero(x, sizeof(x));
	bzero(x0, sizeof(x0));
	bzero(axis, sizeof(axis));
	axis[2] = 1.;
	x[0][1] = 1.;
	x[1][0] = 2.; x[1][1] = 2; x[1][2] = 2.;
	ccd1(axis, x0, 1, x, &x[1], val);
	cout<<val[2] * 180. / acos(-1)<<endl;
	prtdim(3, val);
	getRmat(axis, val[2], R);
	calXrot(R, x[0], x[0]);
	prtdim(3, x[0]);
	for(int m=0; m<3; m++) x[0][m] -= x[1][m];
	cout<<(dot_product(x[0], x[0]))<<endl;
}*/
