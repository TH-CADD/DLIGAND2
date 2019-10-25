#ifndef _MISC
#define _MISC
#include <algorithm>
#include <cassert>
#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <map>
#include <ctime>
using namespace std;

extern int DEBUG;
static const double SMALL_NUM=1.e-99, BIG_NUM=1.e99;
#define die0() do {fprintf(stderr, "\nERROR ON LINE %d OF FILE %s!\n\n", __LINE__, __FILE__); exit(1);} while(0)
#define die1(...) do {warn(__VA_ARGS__); die0();} while(0)
#define die(...) do {warn(__VA_ARGS__); exit(0);} while(0)
//
// string
template <class T>
inline int str2dat(const string &line, T *dat){
   istringstream iss(line);
	int n=0;
   while(! (iss>>dat[n]).fail() ) n++;
	return n;
}
template <class T>
inline bool str2dat(const string &line, const int n, T* dat){
   istringstream iss(line);
	for(int i=0; i<n; i++){
   	if( (iss>>dat[i]).fail() )  return false;
	}
	return true;
}
template <class T>
inline int str2dat(const string &line, vector<T> &dat){
   istringstream iss(line);
	T dt; int n=0, m=dat.size();
   while(! (iss>>dt).fail() ){
		if(n < m) dat[n] = dt;
		else dat.push_back(dt);
		n ++;
	}
	dat.resize(n);
	return n;
}
inline void split_str(const string &s, char delim, vector<string> &elems){
	elems.clear();
	stringstream ss(s+delim);
	string item;
	while(getline(ss, item, delim)){
		elems.push_back(item);
	}
}
inline string trim(const string line){
	int i, is, ie, size=line.size();
	if(size==0) return line;
	for(i=0; i<size; i++)
		if(line.at(i) != ' ') break;
	is = i;
	for(i=size-1; i>=0; i--)
		if(line.at(i) != ' ') break;
	ie = i;
	return line.substr(is,ie-is+1);
}
inline bool endswith(const string s1, const string s2){
	if(s1.size()>=s2.size() && s1.substr(s1.size()-s2.size(), s2.size())==s2) return true;
	return false;
}
inline string rdelete(const string s1, const string s2){
	if(s1.size()>=s2.size() && s1.substr(s1.size()-s2.size(), s2.size())==s2)
		return s1.substr(0, s1.size()-s2.size());
	return s1;
}
//
template <class T>
inline double sum(int n, T x){
	double s=0;
	for(int i=0; i<n; i++) s += x[i];
	return s;
}
inline void cross_product(double *r1,double *r2,double *val){
    val[0] = r1[1]*r2[2]-r1[2]*r2[1];
    val[1] = -r1[0]*r2[2]+r1[2]*r2[0];
    val[2] = r1[0]*r2[1]-r1[1]*r2[0];
}
template <class T>
inline double dot_product(int n, T r1, T r2){
    double value=0;
    for(int i=0;i<n;i++) value += r1[i]*r2[i];
    return value;
}
template <class T>
inline double dot_product(T r1, T r2){
	return dot_product(3, r1, r2);
}
inline double normalize(int n, double *x){
	double r = sqrt(dot_product(n, x, x));
	for(int m=0; m<n; m++) x[m] /= max(r, SMALL_NUM);
	return r;
}
inline double normalize(double *x){
	return normalize(3, x);
}
//
//
inline int getline1(FILE *fp, string &line){
	char *str_line = NULL; size_t str_size;
	int status = getline(&str_line, &str_size, fp);
	if(status < 0) return 0;
	line = string(str_line);
	delete str_line;
	return status;
}
inline FILE *openfile(const string file, const string action="r"){
	FILE *fp = NULL;
	if(file == "--"){
		if(action[0] == 'r') return stdin;
		else if(action[0] == 'w') return stdout;
	}
	if(file=="STDIN" && action[0]=='r') return stdin;
	if(file=="STDOUT" && action[0]=='w') return stdout;
	if(action == "r"){
		if(endswith(file, ".gz")) {
			string cmd1 = "gunzip -c " + file;
			return popen(cmd1.c_str(), "r");
		}
	}
//
	fp = fopen(file.c_str(), action.c_str());
	if(fp == NULL) fprintf(stderr, "Warning, Fail to open file: %s\n", file.c_str());
	return fp;
}
inline bool file_existed(string fn){
	FILE *fp = fopen(fn.c_str(), "r");
	if(fp == NULL) return false;
	fclose(fp); return true;
}
inline bool file_backup(string fn){
	if(! file_existed(fn)) return false;
	rename(fn.c_str(), (fn+'~').c_str()); return true;
}
inline int findargs(int argc, char *argv[], const string tstr){
	for(int i=0; i<argc; i++){
		if(tstr == argv[i]) return i;
	}
	return -1;
}
template <class T>
inline void prtdim(int n, const string fmt1, T dat){
	for(int i=0; i<n; i++) printf(fmt1.c_str(), dat[i]);
	printf("\n");
}
template <class T>
inline void prtdim(int n, T dat, string send="\n"){
	for(int i=0; i<n; i++) cout<<' '<<dat[i];
	cout<<send;
}
template <class T>
inline int minloc(int n, T dat){
	int imin = 0;
	for(int i=1; i<n; i++){
		if(dat[i] < dat[imin]) imin = i;
	}
	return imin;
}
template <class T>
inline int maxloc(int n, T dat){
	int imax = 0;
	for(int i=1; i<n; i++){
		if(dat[i] > dat[imax]) imax = i;
	}
	return imax;
}
template <class T>
inline void copy_dim(int n, T src, T targ){
	for(int i=0; i<n; i++) targ[i] = src[i];
}
template <class T>
inline T sign(T v, T s){
	if(s < 0) return -v;
	return v;
}
inline void warn(const char *fmt, ...){
	va_list ptr; va_start(ptr, fmt);
	vfprintf(stderr, fmt, ptr);
	va_end(ptr);
	fprintf(stderr, "\n");
}
// sort in increasing order
template <class T>
inline void sorti(int n, T *dat){
	typeof(dat[0]) dt;
	for(int i=1; i<n; i++){
		dt = dat[i]; int j;
		for(j=i-1; j>=0; j--){
			if(dat[j] <= dt) break;
			dat[j+1] = dat[j];
		}
		dat[j+1] = dt;
	}
}
//
// T: int,float;  T2: int* or vector<int>   --> increasing order
template <class T, class T2>
inline void sortd(int n, T dat, T2 &seq){
	for(int i=0; i<n; i++) seq[i] = i;
	typeof(dat[0]) dt(0);
	for(int i=1; i<n; i++){
		int it=seq[i], j; dt = dat[it];
		for(j=i-1; j>=0; j--){
			if(dat[seq[j]] <= dt) break;
			seq[j+1] = seq[j];
		}
		seq[j+1] = it;
	}
}
template <class T>
inline T powi(T x, int n){
	if(n == 0) return 1.;
	if(n < 0) return 1. / powi(x, -n);
	T w = x, y = 1;
	if (n & 1) y = x;
	n >>= 1;
	while(n) {
		w = w * w;
		if( n & 1 ) y *= w;
		n >>= 1;
	}
	return y;
}
inline char* chomp(char *str){
	int nt = strlen(str);
	if(str[nt-1] == '\n') str[nt-1] = '\0';
	return str;
}
template<class T, class T2>
inline int findElement(int n, T edim, T2 e, const char *info){
	if(edim[n-1] < e) {
		fprintf(stderr, "%s: findEle bigger than last one: %f %f\n", info, e, edim[n-1]); return n-1;
	} else if(edim[0] > e) {
		return 0;
	}
	int is=0, ie=n-1;
	while(is < ie-1){
		int it = (is+ie) / 2;
		if(edim[it] >= e) ie = it;
		else is = it;
	}
	return ie;
}
inline double cpu_time(){
	static long c0 = time(NULL);
	return (time(NULL) - c0) / 60.;
}
#endif
