#ifndef _ARRAY1
#define _ARRAY1
#include "misc.h"

template <class T, size_t m_dat>
class Array1{
	T dat[m_dat];
 public:
	Array1(){}
	Array1(T v){for(int m=0; m<m_dat; m++) dat[m] = v;}
	Array1(const Array1 &from){for(int m=0; m<m_dat; m++) dat[m] = from.dat[m];}
	Array1(const T *xt){for(int m=0; m<m_dat; m++) dat[m] = xt[m];}
	~Array1(){}
//
	void setx(T vs[]){for(int m=0; m<m_dat; m++) dat[m] = vs[m];}
	void fill(T v){for(int m=0; m<m_dat; m++) dat[m] = v;}
	operator T*(){return dat;}
	T *getdat(){return dat;}
	T *getx(){return dat;}
	T &operator[](size_t it){return dat[it];}
	Array1 operator+(const Array1 &v){
		T xt[m_dat];
		for(int m=0; m<m_dat; m++) xt[m] = dat[m] + v.dat[m];
		return Array1(xt);
	}
	Array1 operator-(const Array1 &v){
		T xt[m_dat];
		for(int m=0; m<m_dat; m++) xt[m] = dat[m] - v.dat[m];
		return Array1(xt);
	}
	Array1 operator*(const T v){
		T xt[m_dat];
		for(int m=0; m<m_dat; m++) xt[m] = dat[m] * v;
		return Array1(xt);
	}
	Array1 operator/(const T v){
		T xt[m_dat];
		for(int m=0; m<m_dat; m++) xt[m] = dat[m] / v;
		return Array1(xt);
	}
	T normalize(){
		double s = 0;
		for(int m=0; m<m_dat; m++) s += dat[m] * dat[m];
		s = max(1e-99, sqrt(s));
		for(int m=0; m<m_dat; m++) dat[m] /= s;
		return s;
	}
	Array1 &operator=(const Array1 &from){
		for(int m=0; m<m_dat; m++) dat[m] = from.dat[m]; return(*this);}
	Array1 &operator+=(const Array1 &from){
		for(int m=0; m<m_dat; m++) dat[m] += from.dat[m]; return(*this);}
	Array1 &operator-=(const Array1 &from) {
		for(int m=0; m<m_dat; m++) dat[m] -= from.dat[m]; return(*this);}
	Array1 &operator*=(const T v) {
		for(int m=0; m<m_dat; m++) dat[m] *= v; return (*this);}
	Array1 &operator/=(const T v) {
		for(int m=0; m<m_dat; m++) dat[m] /= v; return (*this);}
	T dot_product(const Array1 &v){
		T s = 0;
		for(int m=0; m<m_dat; m++) s += dat[m] * v.dat[m];
		return s;
	}
	T distance2(Array1 const &v){
		T x_, s_ = 0;
		for(int m=0; m<m_dat; m++){
			x_ = dat[m] - v.dat[m]; s_ += x_ * x_;
		}
		return s_;
	}
//
	friend inline ostream &operator<<(ostream &os, Array1 const &v){
		for(int m=0; m<m_dat; m++) os<<' '<<v.dat[m];
		os<<endl;
		return os;
	}
	friend inline int operator<(Array1 const &v1, Array1 const &v2){
		for(int m=0; m<m_dat; m++){
			if(v1.dat[m] < v2.dat[m]) return 1;
			else if(v1.dat[m] > v2.dat[m]) return 0;
		}
		return 0;
	}
};

typedef Array1<double, 3> Xvec;

#endif
