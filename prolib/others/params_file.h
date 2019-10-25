#ifndef __PARAMS_FILE
#define __PARAMS_FILE
#include "misc.h"

class Params_input{
	map<string, string> keywords;
public:
	Params_input(){}
	Params_input(string fn){readfile(fn);}
	void readfile(string fn){
		char str[401], ss[2][201];
		FILE *fp = openfile(fparams.c_str(), "r");
		while(fgets(str, 400, fp) != NULL){
			if(*str == '#') continue;
			if(! str2dat(str, 1, ss)) continue;
			keywords[ss[0]] = str;
		}
		fclose(fp);
	}
	bool has_key(string key1){
		return (keywords.count(key1) > 0);
	}
	string get(string key1){
		if(keywords.count(key1) < 1) return "";
		char ss[2][201];
		str2dat(keywords[key1], 1, ss);
		return ss[1];
	}
};
#endif
