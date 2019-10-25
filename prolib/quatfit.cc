#include "misc.h"
/*=======================================================
	a - input: matrix to diagonalize
	v - output: eigenvectors
	d - output: eigenvalues
=========================================================*/

void jacobi(double a[][4], double *d, double v[][4]){
	const int nrot = 30;
	double onorm, dnorm;
	double b, dma, q, t, c, s;
	double atemp, vtemp, dtemp;
	int i, j, k, l;

	for (j = 0; j <= 3; j++) {
		for (i = 0; i <= 3; i++) v[i][j] = 0.0;
		v[j][j] = 1.0;
		d[j] = a[j][j];
	}

	for (l = 1; l <= nrot; l++) {
		dnorm = 0.0;
		onorm = 0.0;
		for (j = 0; j <= 3; j++) {
			dnorm = dnorm + fabs(d[j]);
			for (i = 0; i <= j - 1; i++) {
				onorm = onorm + fabs(a[i][j]);
		  }
		}
		if((onorm/dnorm) <= 1.0e-12) goto Exit_now;
		for (j = 1; j <= 3; j++) {
			for (i = 0; i <= j - 1; i++) {
				b = a[i][j];
				if(fabs(b) > 0.0) {
					dma = d[j] - d[i];
					if((fabs(dma) + fabs(b)) <=  fabs(dma)) {
						t = b / dma;
				}
				else {
					q = 0.5 * dma / b;
					t = 1.0/(fabs(q) + sqrt(1.0+q*q));
					if(q < 0.0) {
						t = -t;
					}
				}
				c = 1.0/sqrt(t * t + 1.0);
				s = t * c;
				a[i][j] = 0.0;
				for (k = 0; k <= i-1; k++) {
					atemp = c * a[k][i] - s * a[k][j];
					a[k][j] = s * a[k][i] + c * a[k][j];
					a[k][i] = atemp;
					}
				for (k = i+1; k <= j-1; k++) {
					atemp = c * a[i][k] - s * a[k][j];
					a[k][j] = s * a[i][k] + c * a[k][j];
					a[i][k] = atemp;
					}
				for (k = j+1; k <= 3; k++) {
					atemp = c * a[i][k] - s * a[j][k];
					a[j][k] = s * a[i][k] + c * a[j][k];
					a[i][k] = atemp;
					}
				for (k = 0; k <= 3; k++) {
					vtemp = c * v[k][i] - s * v[k][j];
					v[k][j] = s * v[k][i] + c * v[k][j];
					v[k][i] = vtemp;
				}
				dtemp = c*c*d[i] + s*s*d[j] - 2.0*c*s*b;
				d[j] = s*s*d[i] + c*c*d[j] +  2.0*c*s*b;
				d[i] = dtemp;
		    }  /* end if */
		  } /* end for i */
		} /* end for j */
	} /* end for l */
	
Exit_now:

//	nrot = l;

	for (j = 0; j <= 2; j++) {
		k = j; dtemp = d[k];
		for (i = j+1; i <= 3; i++) {
			if(d[i] < dtemp) {
				k = i; dtemp = d[k];
			}
		}

		if(k > j) {
			d[k] = d[j];
			d[j] = dtemp;
			for (i = 0; i <= 3; i++) {
				dtemp = v[i][k];
				v[i][k] = v[i][j];
				v[i][j] = dtemp;
			}
		}
	}
}
void q2mat (double *q, double u[][3]) {
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
void quatfit(int n, double *w, double x1[][3], double x2[][3], double u[][3]){
	double q[4];
	double xxyx, xxyy, xxyz;
	double xyyx, xyyy, xyyz;
	double xzyx, xzyy, xzyz;
	double c[4][4], v[4][4];
	double d[4];
// generate the upper triangle of the quadratic form matrix
	xxyx = xxyy = xxyz = xyyx = xyyy = 0.;
	xyyz = xzyx = xzyy = xzyz = 0.0;
	
	for(int i = 1; i <= n; i++) {
		xxyx = xxyx + x2[i][0] * x1[i][0] * w[i];
		xxyy = xxyy + x2[i][0] * x1[i][1] * w[i];
		xxyz = xxyz + x2[i][0] * x1[i][2] * w[i];
		xyyx = xyyx + x2[i][1] * x1[i][0] * w[i];
		xyyy = xyyy + x2[i][1] * x1[i][1] * w[i];
		xyyz = xyyz + x2[i][1] * x1[i][2] * w[i];
		xzyx = xzyx + x2[i][2] * x1[i][0] * w[i];
		xzyy = xzyy + x2[i][2] * x1[i][1] * w[i];
		xzyz = xzyz + x2[i][2] * x1[i][2] * w[i];
	}
	for(int i=0; i<=3; i++)
	for(int j=0; j<=3; j++){
		c[i][j] = 0.0;
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
	q2mat(q, u);
}
void translate1(int na, double *w, double x[][3], double xn[][3], double *xc){
	double ws=0.;
	for(int m=0; m<3; m++) xc[m] = 0.;
	for(int i=0; i<na; i++){
		for(int m=0; m<3; m++) xc[m] += x[i][m] * w[i];
		ws += w[i];
	}
	if(ws < 1.0e-6) die("too small w: %g\n", ws);
	for(int m=0; m<3; m++) xc[m] /= ws;
	for(int i=0; i<na; i++){
		for(int m=0; m<3; m++) xn[i][m] = x[i][m] - xc[m];
	}
}
void rotmol(int na, double u[][3], double x[][3], double xn[][3]){
	double xt[3];
	for(int i=0; i<na; i++){
		for(int m=0; m<3; m++) xt[m] = dot_product(u[m], x[i]);
		for(int m=0; m<3; m++) xn[i][m] = xt[m];
	}
}
/*int main(){
	double x0[][3] = {{-14.152,0.961,4.712},{-13.296,0.028,3.924},{-11.822,0.338,4.193}};
	double y0[][3] = {{-11.121,-0.642,4.703},{-9.669,-0.447,4.998},{-8.861,-1.586,4.373}};
	double w[] = {1., 1., 1.}; double u[4][3];
	double x[3][3], y[3][3], xc1[3], xc2[3], xt[3];
	translate(3, w, x0, x, xc1);
	translate(3, w, y0, y, xc2);
	for(int m=0; m<3; m++) prtdim(3, x[m]);
	for(int m=0; m<3; m++) prtdim(3, y[m]);
	quatfit(3, w, x, y, u);
	for(int m=0; m<3; m++) prtdim(3, u[m]);
	cout<<endl;
	for(int m=0; m<3; m++) u[3][m] = xc1[m] - dot_product(u[m], xc2);
//	rotmol(3, u, y, y);
	for(int i=0; i<3; i++){
		for(int m=0; m<3; m++) xt[m] = dot_product(u[m], y0[i]) + u[3][m];
		prtdim(3, xt);
	}
}*/
