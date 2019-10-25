#include "protein.h"

/*template <class T>
inline bool str2dat(const string &line, T &dat, ios_base & (*fmt)(ios_base&)) {
   istringstream iss(line);
   return !(iss>>fmt>>dat).fail();
}*/
/*inline void die(const char *fmt, ...){
	va_list ptr; va_start(ptr, fmt);
	warn(fmt, ptr);
	exit(1);
}*/
void initPRtype(int flag){
	die("initPRtype replaced by initRestypes_pro or initRestypes\n");
}
int defineAtype(string rn, string an){
	die("defineAtype replaced by atomDefine\n");
	return 0;
}
