#include "protein2.h"

Protein2::Protein2(const string fn){
	pdbnm = fn;
	natm = nres = 0;
	rstart_ch.push_back(0);
	rdpdb(fn);
}
void Protein2::rdpdb(string fn){
	FILE *fp=openfile(fn, "r");
	rdpdb(fp); fclose(fp);
}
void Protein2::rdpdb(FILE *fp){
	char str[121], an[5];
	string line, resinfo0="";
	double xt[3];
	while(fgets(str,120,fp) != NULL){
		line = str;
		if(line.substr(0,3) == "END") break;
		if(line.substr(0,4) != "ATOM") continue;
		string info = line.substr(17, 10);
		sscanf(str+12, "%4s", an); an[3] = '\0';
		if(an[0] == 'H') continue;
		if(nres == 0 || info != resinfo0) {
			resinfo.push_back(info);
			resinfo0 = info;
			rnam.push_back(info.substr(0,3));
			if(rstart.size() <= nres) rstart.push_back(natm); // false when reading 2nd chain
			nres ++;
		}
		anam.push_back(an);
		rseq.push_back(nres-1);
		sscanf(str+30, "%lf%lf%lf", xt, xt+1, xt+2);
		x.push_back( Xvec(xt) );
		natm ++;
	}
	rstart.push_back(natm);
	rstart_ch.push_back(natm);
}
void Protein2::wrpdb(string fn){
	FILE *fp = openfile(fn, "w");
	wrpdb(fp); fclose(fp);
}
void Protein2::wrpdb(FILE *fp){
	char fmt1[] = "ATOM%7d  %-4s%3s %c%4d    %8.3f%8.3f%8.3f\n";
	for(int ii=0; ii<rstart_ch.size(); ii++){
		char chn = 'A' + ii;
		for(int i=rstart_ch[ii]; i<rstart_ch[ii+1]; i++){
			int ir = rseq[i];
			double *xt = x[i].getx();
			fprintf(fp, fmt1, i+1, anam.at(i).c_str(), rnam.at(ir).c_str(), chn,
					ir+1, *xt, xt[1], xt[2]);
		}
		fprintf(fp, "TER\n");
	}
	fprintf(fp, "END\n");
}
