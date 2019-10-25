#include "array1.h"
#include "misc.h"

void jacobi(double a[][4], double *d, double v[][4]);
void quatfit(int n, vector<double> &w, vector<Xvec> &x1, vector<Xvec> &x2, double u[][3]){
	double q[4];
	double xxyx, xxyy, xxyz;
	double xyyx, xyyy, xyyz;
	double xzyx, xzyy, xzyz;
	double c[4][4], v[4][4];
	double d[4];
// generate the upper triangle of the quadratic form matrix
	xxyx = xxyy = xxyz = 0.;
	xyyx = xyyy = xyyz = 0.;
	xzyx = xzyy = xzyz = 0.;
	
	for(int i=0; i<n; i++){
		xxyx += x2[i][0] * x1[i][0] * w[i];
		xxyy += x2[i][0] * x1[i][1] * w[i];
		xxyz += x2[i][0] * x1[i][2] * w[i];
		xyyx += x2[i][1] * x1[i][0] * w[i];
		xyyy += x2[i][1] * x1[i][1] * w[i];
		xyyz += x2[i][1] * x1[i][2] * w[i];
		xzyx += x2[i][2] * x1[i][0] * w[i];
		xzyy += x2[i][2] * x1[i][1] * w[i];
		xzyz += x2[i][2] * x1[i][2] * w[i];
	}

	c[0][0] = xxyx + xyyy + xzyz;
	c[0][1] = xzyy - xyyz;
	c[1][1] = xxyx - xyyy - xzyz;
	c[0][2] = xxyz - xzyx;
	c[1][2] = xxyy + xyyx;
	c[2][2] = xyyy - xzyz - xxyx;
	c[0][3] = xyyx - xxyy;
	c[1][3] = xzyx + xxyz;
	c[2][3] = xyyz + xzyy;
	c[3][3] = xzyz - xxyx - xyyy;

// diagonalize c
	jacobi (c, d, v);


// extract the desired quaternion
	q[0] = v[0][3];
	q[1] = v[1][3];
	q[2] = v[2][3];
	q[3] = v[3][3];

// generate the rotation matrix
	u[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
	u[1][0] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
	u[2][0] = 2.0 * (q[1] * q[3] + q[0] * q[2]);

	u[0][1] = 2.0 * (q[2] * q[1] + q[0] * q[3]);
	u[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
	u[2][1] = 2.0 * (q[2] * q[3] - q[0] * q[1]);

	u[0][2] = 2.0 *(q[3] * q[1] - q[0] * q[2]);
	u[1][2] = 2.0 * (q[3] * q[2] + q[0] * q[1]);
	u[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
}

void translate1(int na, vector<Xvec> &x0, vector<Xvec> &xn, double *xc){
	for(int m=0; m<3; m++) xc[m] = 0.;
	for(int i=0; i<na; i++){
		for(int m=0; m<3; m++) xc[m] += x0[i][m];
	}
	for(int m=0; m<3; m++) xc[m] /= na;
	for(int i=0; i<na; i++){
		for(int m=0; m<3; m++) xn[i][m] = x0[i][m] - xc[m];
	}
}
void translate2(int na, vector<double> wfit, vector<Xvec> &x0, vector<Xvec> &xn, double *xc){
	double ds = 0.;
	for(int m=0; m<3; m++) xc[m] = 0.;
	for(int i=0; i<na; i++){
		for(int m=0; m<3; m++) xc[m] += x0[i][m] * wfit[i];
		ds += wfit[i];
	}
	for(int m=0; m<3; m++) xc[m] /= ds;
	for(int i=0; i<na; i++){
		for(int m=0; m<3; m++) xn[i][m] = x0[i][m] - xc[m];
	}
}
void rotmol(int na, double u[][3], vector<Xvec> &x0, vector<Xvec> &xn){
	double xt[3];
	for(int i=0; i<na; i++){
		for(int m=0; m<3; m++) xt[m] = dot_product(u[m], x0[i].getx());
		for(int m=0; m<3; m++) xn[i][m] = xt[m] + u[3][m];
	}
}
double calrms1(vector<Xvec> &xa, vector<Xvec> &xb){
	int na=xa.size();
	double xc[3], u[4][3]; //bzero(u, sizeof(u));
	vector<Xvec> x1(na), x2(na);
	vector<double> w(na, 1.);
	translate1(na, xa, x1, xc);
	translate1(na, xb, x2, xc);
	quatfit(na, w, x1, x2, u);
	for(int m=0; m<3; m++) u[3][m] = 0.;		// very important
	rotmol(na, u, x2, x2);
	double rms1 = 0.;
	for(int i=0; i<na; i++){
		rms1 += x1[i].distance2(x2[i]);
	}
	return sqrt(rms1/na);
}
