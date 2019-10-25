#ifndef _XVEC
#define _XVEC

class Xvec{
	double xv, yv, zv;
	friend inline ostream & operator <<( ostream & os, Xvec const & v );
	friend inline void cross_product(const Xvec &v1, const Xvec &v2, Xvec &v3);
 public:

	Xvec(double xi=0.0, double yi=0.0, double zi=0.0) : xv(xi), yv(yi), zv(zi){}

	Xvec(const Xvec &from){xv=from.xv; yv=from.yv; zv=from.zv;}
	Xvec(const double *x){xv=*x; yv=x[1]; zv=x[2];}

	inline double operator[](int inx) const{ return (inx==0?xv:inx==1?yv:zv); }
	inline Xvec operator-() const{ return Xvec(-xv, -yv, -zv);}
	inline Xvec operator-(const Xvec &v) const{ return Xvec(xv-v.xv, yv-v.yv, zv-v.zv); }
	inline Xvec operator+(const Xvec &v) const{ return Xvec(xv+v.xv, yv+v.yv, zv+v.zv); }
	inline double operator*(Xvec &v) const{ return xv*v.xv + yv*v.yv + zv*v.zv; }
	inline Xvec operator/(double v) const{ return Xvec(xv/v, yv/v, zv/v); }
	inline Xvec operator*(double v) const{ return Xvec(xv*v, yv*v, zv*v); }

	inline Xvec &operator=(const Xvec& from){
		xv=from.xv; yv=from.yv; zv=from.zv; return(*this);
	}
	inline Xvec &operator+=(const Xvec& from){
		xv+=from.xv; yv+=from.yv; zv+=from.zv; return(*this);
	}
	inline Xvec &operator-=(const Xvec& from){
		xv-=from.xv; yv-=from.yv; zv-=from.zv; return(*this);
	}

	inline double normalize(){
		double r=sqrt(xv*xv + yv*yv + zv*zv), len=r;
		if(len < 1.0e-8) len=1.0e-8;
		xv/=len; yv/=len; zv/=len;
		return r;
	}
};
inline ostream & operator <<( ostream & os, Xvec const & v ){
	os<<v.xv<<' '<<v.yv<<' '<<v.zv; return os;
}
inline void cross_product(const Xvec &v1, const Xvec &v2, Xvec &v3){
	v3.xv = v1.yv*v2.zv - v1.zv*v2.yv;
	v3.yv = v1.zv*v2.xv - v1.xv*v2.zv;
	v3.zv = v1.xv*v2.yv - v1.yv*v2.xv;
}

#endif
