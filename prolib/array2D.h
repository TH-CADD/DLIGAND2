#ifndef _ARRAY
#define _ARRAY
#include "misc.h"

// The order is x(0,0..n-1), x(1,..)...(m-1,n-1)
template<class T>
class Array1D{
protected:
	long length, size;
	T *dat;
public:
	Array1D(long n0):length(n0){dat = new T[n0]; size = n0;}
	Array1D():length(0){dat = NULL; size=0;}
	~Array1D(){clean();}
	void clean(){if(dat != NULL) delete [] dat; dat = NULL;}
	void init(T v){for(long i=0; i<length; i++) dat[i] = v;}
	void init(T *v){for(long i=0; i<length; i++) dat[i] = v[i];}
	void init(){init(T(0));}
	T *getx(){return dat;}
	T& operator[](long i0){
# ifdef GDB
		if(i0>=length || i0<0) die("overflow: %ld %ld", i0, length);
#endif
		return dat[i0];
	}
	void operator/=(T t0){ for(int i=0; i<length; i++) dat[i] /= t0; return;}
	void operator*=(T t0){ for(int i=0; i<length; i++) dat[i] *= t0; return;}
	void operator+=(T t0){ for(int i=0; i<length; i++) dat[i] += t0; return;}
	void operator-=(T t0){ for(int i=0; i<length; i++) dat[i] -= t0; return;}
	void resize0(const long l, bool bcopy=1){
		if(l<=size) {length = l; return;}
		else if(size == 0) {
			dat = new T[l]; length = size = l;
		} else {
			T *dat2 = new T[l];
			if(bcopy) memcpy(dat2, dat, length*sizeof(T));
			delete [] dat; dat = dat2;
			length = size = l;
		}
	}
};
//
template<class T>
class Array2D:public Array1D<T>::Array1D{
private:
	long n0, n1;
public:
	Array2D(Array2D const &v1):Array1D<T>::Array1D(0){resize(v1.n0, v1.n1, 0); init(v1.dat); }
	Array2D(){n0=n1=0;}
	Array2D(long i0, long i1):n0(i0),n1(i1),Array1D<T>::Array1D(i0*i1){}
	T& operator()(long i0, long i1){return (*this)[i1+n1*i0];}
	void resize(const long i0, const long i1, bool bcopy=1){n0=i0; n1=i1; Array1D<T>::resize0(n0*n1, bcopy);}
	void reshape(const long i0, const long i1){
		assert(i0*i1 <= Array1D<T>::length); n0 = i0; n1 = i1;
	}
	Array2D<T>& operator=(const Array2D<T>& v1 ){
		resize(v1.n0, v1.n1, 0); init(v1.dat);
		return *this;
	}
};
//
template<class T>
class Array3D:public Array1D<T>::Array1D{
private:
	long n0, n1, n2;
public:
	Array3D(){n0=n1=n2=0;}
	Array3D(const long i0, const long i1, const long i2):n0(i0),n1(i1),n2(i2),Array1D<T>::Array1D(i0*i1*i2){}
	T& operator()(const long i0, const long i1, const long i2){return (*this)[i2+n2*(i1+n1*i0)];}
	void reshape(const long i0, const long i1, const long i2){
		assert(i0*i1*i2 <= Array1D<T>::length);
		n0 = i0; n1 = i1; n2 = i2;
	}
	void resize(const long i0, const long i1, const long i2){n0=i0; n1=i1; n2=i2; Array1D<T>::resize0(n0*n1*n2);}
	Array3D<T>& operator=(const Array3D<T>& v1 ){
		resize(v1.n0, v1.n1, v1.n2, 0); init(v1.dat);
		return *this;
	}
};
//
template<class T>
class Array4D:public Array1D<T>::Array1D{
private:
	long n0, n1, n2, n3;
public:
	Array4D(){n0=n1=n2=n3=0;}
	Array4D(const long i0, const long i1, const long i2, const long i3):n0(i0),n1(i1),n2(i2),n3(i3),Array1D<T>::Array1D(i0*i1*i2*i3){}
	T& operator()(long i0, long i1, long i2, long i3){return (*this)[i3+n3*(i2+n2*(i1+i0*n1))];}
	void reshape(const long i0, const long i1, const long i2, const long i3){
		assert(i0*i1*i2*i3 <= Array1D<T>::length);
		n0 = i0; n1 = i1; n2 = i2; n3 = i3;
	}
	void resize(const long i0, const long i1, const long i2, const long i3){n0=i0; n1=i1; n2=i2; n3=i3; Array1D<T>::resize0(n0*n1*n2*n3);}
	Array4D<T>& operator=(const Array4D<T>& v1 ){
		resize(v1.n0, v1.n1, v1.n2, v1.n3, 0); init(v1.dat);
		return *this;
	}
};

#endif
