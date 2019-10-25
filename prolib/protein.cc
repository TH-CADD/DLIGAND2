#include "protein.h"
#include "rotamer.h"

int DEBUG = 0;
string datadir = "";
Protein::Protein(const string fn){
	reset();
	if(fn == "--"){
		pdbnm = "STDIN";
		rdpdb(stdin); return;
	}
	pdbnm = fn;
	FILE *fp=openfile(fn, "r");
	rdpdb(fp); fclose(fp);
}
void Protein::reset(){
	pdbnm = "";
	nchain = natm = nres = 0;
	chainID.clear();
	rnam.clear(); resinfo.clear();
	resid.clear(); seq0.clear(); chnseq.clear();
	rstart.clear(); pdbidx.clear();
	rseq.clear(); atoms.clear();
	rstart.push_back(0);
}
void Protein::popResidues(int nr){
	if(nr <= 0) return;
	assert(nres >= nr);
	nres -= nr;
	natm = rstart[nres];
	rnam.resize(nres);
	resinfo.resize(nres);
	resid.resize(nres);
	seq0.resize(nres);
	chnseq.resize(nres);
	rstart.resize(nres+1);
//
	pdbidx.resize(natm);
	rseq.resize(natm);
	atoms.resize(natm);
}
void Protein::rdpdb(FILE *fp){
	char str[121], an[5], chn0 = '\n';
	string line, s0, s1, status0, rn3;
	int ires, oires;
//
	if(restypes.size()==0) initRestypes_pro();
	if(rstart.size() == 0) rstart.push_back(0);
//
	oires=-99999;
//	nchain = natm = nres = 0;
	bool bres = false;
	while(fgets(str,120,fp) != NULL){
		line = str;
		if(line.substr(0,3) == "END") break;
		if(line.substr(0,3) == "TER") {status0 = "TER"; continue;}
		if(line.size() < 20) continue;
		rn3 = line.substr(17,3);
		if(line.substr(0,6) == "HETATM"){
			if(convert_hetres(rn3) == rn3) continue;
			rn3 = convert_hetres(rn3);
		} else if(line.substr(0,4) != "ATOM") continue;
//		if(line.at(12)=='H' || line.at(13)=='H') {pdbidx.push_back(-1); continue;}
		if(str[21] != chn0 || status0 == "TER") {nchain ++; chn0 = str[21]; chainID.push_back(chn0);}
//
		s1 = line.substr(17,10); ires=strtol(str+22,NULL,10);
		if(s0!=s1 || oires!=ires){
			if(DEBUG && ires-oires!=1 && oires!=-9999 && !chnseq.empty() && nchain==chnseq.back()){
				fprintf(stderr, "Warning, Residue seqens not adjacent: %d %d\n", oires, ires);
			}
			bres = addRes(rn3);
			if(bres){
				seq0.push_back(ires); resinfo.push_back(s1);
				chnseq.push_back(nchain - 1);
			}
			s0 = s1; oires = ires;
		}
		if(! bres) {
			pdbidx.push_back(-1);
		} else {
			if(! initAtom(line)) pdbidx.push_back(-1);
		}
		status0 = "ATOM";
	}
	if(DEBUG > 0){
		fprintf(stderr, "nres/natm = %d/%d\n", nres, natm);
	}
//	fclose(fp);
/*	vector<Amino>::iterator pos;
	for(pos=AAs.begin(); pos!=AAs.end(); pos++) pos->print();*/
}
void Protein::addAtom(string line){
}
// require to call "initRestypes" first
bool Protein::addRes(string rn0){
	char rn[4]; sscanf(rn0.c_str(), "%s", rn);
	int rid = resDefine(rn);
	return addRes(rid);
}
bool Protein::addRes(int rid){
	if(rid<0) return false;
	if(rid >= restypes.size()) die("Resid too big: %d\n", rid);
	resid.push_back(rid);
	rnam.push_back(restypes[rid].getname());
	vector <Atomtype> ap = restypes.at(rid).getatypes();
	for(int i=0; i<ap.size(); i++){
		string an = ap[i].getname();
		int it = ap[i].gettype();
		atoms.push_back(Atom(natm, it, an));
		rseq.push_back(nres);
		natm ++;
	}
	rstart.push_back(natm);
	nres ++;
	return true;
}
bool Protein::initAtom(string line){
	char an1 = line[12], an2 = line[16], an0[6];
	sscanf(line.substr(12,4).c_str(), "%s", an0);
	string an = an0;
	if(an1=='1' || an1=='2' || an1=='3'){
		an = string(an0+1) + an1;
	} else if(an2=='1' || an2=='2' || an2=='3'){
		an += an2;
	}
//
	string rn = rnam.at(nres-1);
	if(an == "OXT") return false;
	else if(an[0] == 'H') return false;
	else if(an=="CD" && rn=="ILE") an="CD1";
	else if(rn.size() < 3){		// DNA or RNA
		if(an == "OP1") an = "O1P";
		else if(an == "OP2") an = "O2P";
		else if(an == "OP3") return false;
		else if(rn=="DT" && an=="C7") an="C5M";
	}
	for(int i=rstart[nres-1]; i<rstart[nres]; i++){
		if(an == atoms[i].getname()){
			atoms[i].setx(line.substr(30, 24));
			pdbidx.push_back(i);
			return true;
		}
	}
	if(DEBUG>0) {
		fprintf(stderr, "Warning, inconsistent atoms: %s %s\n", rnam[nres-1].c_str(), an.c_str());
	}
	return false;
}
int Protein::findatom(int ir, string an){
	if(ir>=nres || ir<0) return -1;
	for(int i=rstart.at(ir); i<rstart.at(ir+1); i++){
		if(atoms.at(i).getname() == an) return i;
	}
	fprintf(stderr, "Fail to find atom: %d, %s\n", ir, an.c_str());
	return -1;
}
void Protein::wrpdb(string fn){
	FILE *fp = openfile(fn, "w");
	wrpdb(fp); fclose(fp);
}
void Protein::wrpdb(FILE *fp){
	char fmt1[] = "ATOM%7d  %-4s%3s %c%4d%12.3f%8.3f%8.3f\n";
	char chn0 = '\n';
	for(int i=0; i<atoms.size(); i++){
		Atom *ap = &atoms.at(i);
		if(! ap->filled()) continue;
		int ir = rseq.at(i); double *x=ap->getx();
		char chn = 'A';
		if(chnseq.size()>ir && chainID.size()>chnseq[ir]){
			chn = chainID.at(chnseq[ir]);
			if(chn0!='\n' && chn!=chn0) fprintf(fp, "TER\n");
			chn0 = chn;
		} else if(DEBUG){
			fprintf(stderr, "not defined chain name for RES%d: %d %d\n", ir, chnseq.size(), chainID.size());
		}
//		if(chainID.size() > ir) chn = chainID[ir];
		int rid = ir + 1;
		if(seq0.size() > ir) rid = seq0[ir];
		fprintf(fp, fmt1, i+1, ap->getname().c_str(), rnam.at(ir).c_str(), chn,
				rid, *x, x[1], x[2]);
	}
	fprintf(fp, "TER\nEND\n");
}
//
//
string Protein::getfasta(){
	string s1;
	for(int i=0; i<nres; i++) s1 += rnam1_std[resid[i]];
	return s1;
}
//
void Protein::xdiff(int ia, int ib, double *xd){
	void xdiff(double *xa, double *xb, double *xab);
	xdiff(atoms[ia].getx(), atoms[ib].getx(), xd);
}
double Protein::distance2(int ia, int ib){
	double distance2(double *xa, double *xb);
	return distance2(atoms[ia].getx(), atoms[ib].getx());
}
double Protein::angle(int ia, int ib, int ic){
	double angle(double *xa, double *xb, double *xc);
	return angle(atoms[ia].getx(), atoms[ib].getx(), atoms[ic].getx());
}
double Protein::torsion(int ia, int ib, int ic, int id){
	double torsion(double *xa, double *xb, double *xc, double *xd);
	return torsion(atoms[ia].getx(), atoms[ib].getx(), atoms[ic].getx(), atoms[id].getx());
}

void Protein::calphipsi(int ir, double *ag){
	int i1 = rstart[ir];
	if(ir > 0){
		int i0 = rstart[ir-1];
		ag[0] = torsion(i1+2, i1+1, i1, i0+2);
	} else ag[0] = 360.;
	ag[1] = torsion(i1+3, i1+2, i1+1, i1);
	ag[1] = fmod(ag[1]+360., 360.) - 180.;
}
bool Protein::iscomplete(){
	bool bcomp=true;
	for(int ir=0; ir<nres; ir++){
		int id = resid.at(ir);
		if(id < 0){
			if(ir!=0 && ir!=nres-1) bcomp = false;
			fprintf(stderr, "Error residue type: %s %d", rnam[ir].c_str(), ir);
			continue;
		}
		Restype *rp = &restypes.at(id);
		vector <Atomtype> ap=rp->getatypes();
		for(int i=0; i<ap.size(); i++){
			bool finish=false;
			string at = ap[i].getname();
			for(int ia=rstart[ir]; ia<rstart[ir+1]; ia++){
				if(atoms.at(ia).getname() == at) {finish=true; break;}
			}
			if(! finish){
				if(ir!=0 && ir!=nres-1) bcomp = false;
				fprintf(stderr, "Lost atoms: %s %d %s", rnam[ir].c_str(), seq0[ir], at.c_str());
			}
		}
	}
	return bcomp;
}
void Protein::initneib(){
	for(int i=0; i<natm; i++) atoms[i].clearNeib();
	for(int ir=1; ir<nres; ir++){
		int id = resid.at(ir);
		if(id < 0) continue;
		if(chnseq[ir] != chnseq[ir-1]) continue;
		if(rnam[ir-1].size() != rnam[ir].size()) continue;		// to prevent between pro/DNA
		if(rnam[ir].size() == 3){		// amino acid
			int i1=rstart[ir-1]+2, i2=rstart[ir];
			atoms.at(i2).addneib(i1);		//N --- C
			atoms.at(i1).addneib(i2);
			if(atoms[i2].getname()!="N" || atoms[i1].getname()!="C"){
				die("wrong atom name for AA neib: %d %d\n", i1, i2);
			}
		} else {
			int i1=rstart[ir-1]+8, i2=rstart[ir];
			atoms.at(i2).addneib(i1);		//P --- O3'
			atoms.at(i1).addneib(i2);
			if(atoms[i2].getname()!="P" || atoms[i1].getname()!="O3'"){
				die("wrong atom name for AA neib: %d %d\n", i1, i2);
			}
		}
	}
//
	for(int ir=0; ir<nres; ir++){
		int id = resid.at(ir);
		if(id < 0) continue;
		Restype *rp = &restypes.at(id);
		vector <string> *neib = rp->getneibs();
		for(int i=0; i<neib->size(); i++){
			int i1=-1, i2=-1;
			for(int ia=rstart[ir]; ia<rstart[ir+1]; ia++){
				if(atoms.at(ia).getname() == neib[0].at(i)) i1=ia;
				else if(atoms.at(ia).getname() == neib[1].at(i)) i2=ia;
				if(i1>=0 && i2>=0) break;
			}
			atoms.at(i1).addneib(i2);
			atoms.at(i2).addneib(i1);
		}
	}
}
void Protein::initneib1_sc(int ir){		// specially for mutation
	int id = resid.at(ir);
	for(int i=rstart[ir]; i<rstart[ir+1]; i++) atoms[i].clearNeib();
	Restype *rp = &restypes.at(id);
	vector <string> *neib = rp->getneibs();
	for(int i=0; i<neib->size(); i++){
		int i0 = findatom(ir, neib[0][i]);
		int i1 = findatom(ir, neib[1][i]);
		if(i0<0 || i1<0) continue;
		atoms.at(i0).addneib(i1);
		atoms.at(i1).addneib(i0);
	}
}
// Att! aseq afterwards not changed
void Protein::eraseAtoms(int idx, int na){
	if(idx<=0 || na<=0) die("Error erase: %d %d\n", idx, na);
	for(int ir=rseq[idx]+1; ir<=nres; ir++) rstart[ir] -= na;
	rseq.erase(rseq.begin()+idx, rseq.begin()+idx+na);
	atoms.erase(atoms.begin()+idx, atoms.begin()+idx+na);
	natm -= na;
}
void Protein::insertAtoms(int idx, int na){
	if(idx<=0 || na<=0) die("Error insert: %d %d\n", idx, na);	// in front
	for(int ir=rseq[idx-1]+1; ir<=nres; ir++) rstart[ir] += na;
	rseq.insert(rseq.begin()+idx, na, rseq[idx-1]);
	atoms.insert(atoms.begin()+idx, na, Atom());
	natm += na;
}
void Protein::initAtomint(int i0){
	int ir = rseq[i0];
	vector<int> np0 = atoms[i0].getineib();
	int ia, ib, ic, chiral;
	ia = ib = ic = -1; chiral=0;
	string rn = rnam[ir];
	if(np0.size() == 0) {
		fprintf(stderr, "No neib: %d %s\n", ir, atoms[i0].getname().c_str()); exit(1);
	}
	if(np0[0] < i0) {
		ia = np0[0];
// in PRO, to replace N
		if(rn=="PRO" && atoms[i0].getname()=="CD") ia = findatom(ir, "CG");

		vector<int> np1 = atoms[ia].getineib();
		if(np1.at(0) < i0) {
			ib = np1[0];
			if(np1.at(1) < i0) {ic = np1[1]; chiral = 1;}
			else {
				vector<int> np2 = atoms[ib].getineib();
				ic = np2.at(0);
				if(ic == ia) {
					if(np2.size()>1 && np2[1]<i0) ic = np2[1];
					else ic = -1;
				}
			}
		}
	}
	int iz[] = {ia, ib, ic, chiral};
	atoms[i0].setint(iz);
}
void Protein::initAtomzcrd0(int i0){
	int *iz = atoms[i0].getint();
	int ia=iz[0], ib=iz[1], ic=iz[2], chiral=iz[3];
	if(restpl == NULL) restpl = new Restpl("");
//
	string rn = rnam[ rseq[i0] ];
	int idx[4]; double ag[3];
	for(int i=0; i<4; i++) idx[i] = -1;
	if(ic<0 || rseq[ic]==rseq[i0]) {
		idx[0] = restpl->findatom_tpl(rn, atoms[i0].getname());
		if(ia >= 0) idx[1] = restpl->findatom_tpl(rn, atoms[ia].getname());
		if(ib >= 0) idx[2] = restpl->findatom_tpl(rn, atoms[ib].getname());
		if(ic >= 0) idx[3] = restpl->findatom_tpl(rn, atoms[ic].getname());
//
	} else if (rseq[ic] != rseq[i0]) {		// main-chain atoms
		string rn1 = "ALA", rn2 = "CYS";		// adjacent residues in template
		if(rnam[rseq[i0]].size() == 2) {rn1="DA"; rn2="DT";} // DNA
		else if(rnam[rseq[i0]].size() == 1) {rn1="A"; rn2="U";} // RNA
		idx[0] = restpl->findatom_tpl(rn2, atoms[i0].getname());
		if(rseq[ia] == rseq[i0]) idx[1] = restpl->findatom_tpl(rn2, atoms[ia].getname());
		else idx[1] = restpl->findatom_tpl(rn1, atoms[ia].getname());
		if(rseq[ib] == rseq[i0]) idx[2] = restpl->findatom_tpl(rn2, atoms[ib].getname());
		else idx[2] = restpl->findatom_tpl(rn1, atoms[ib].getname());
		idx[3] = restpl->findatom_tpl(rn1, atoms[ic].getname());
	}
	restpl->calint_tpl(idx, chiral, ag);
//
	if(ia < 0) chiral = 10;
	else if(ib < 0) chiral = 11;
	else if(ic < 0) chiral = 12;
	atoms[i0].setchiral(chiral);
	atoms[i0].setzcrd(ag);
}
void Protein::calAtomzcrd(int ia){
	int *iz = atoms[ia].getint();
	double *ag = atoms[ia].getzcrd();
	if(iz[0] < 0) return;
	ag[0] = sqrt( distance2(ia, iz[0]) );
	if(iz[1] < 0) return;
	ag[1] = angle(ia, iz[0], iz[1]);
	if(iz[2] < 0) return;
	ag[2] = torsion(ia, iz[0], iz[1], iz[2]);
	if(iz[3] != 0) {
		iz[3] = 1;
		if( ag[2] > 0) iz[3] = -1;
		ag[2] = angle(ia, iz[0], iz[2]);
	}
}
void Protein::addCB_GLY(int ir){
	if(rnam[ir] != "GLY") die("not GLY: %d %s\n", ir, rnam[ir].c_str()); 
	int ib = findatom(ir, "O") + 1;
	insertAtoms(ib, 1);
	atoms[ib].setname("CB", -1);
//
	if(restpl == NULL) restpl = new Restpl("");
	string rn="ALA"; int idx[4], chiral=1; double ag[4];
	idx[0] = restpl->findatom_tpl(rn, "CB"); idx[1] = restpl->findatom_tpl(rn, "CA");
	idx[2] = restpl->findatom_tpl(rn, "N"); idx[3] = restpl->findatom_tpl(rn, "C");
	restpl->calint_tpl(idx, chiral, ag);
	atoms[ib].setzcrd(ag);
//
	int it[4];
	it[0] = findatom(ir, "CA"); it[1] = findatom(ir, "N");
	it[2] = findatom(ir, "C"); it[3] = chiral;
	atoms[ib].setint(it); makexyz1(ib);
}
int Protein::getReschi(int ir, int *idx){
	string rn = rnam[ir];
	int nchi=0, mchi=4;
	if(rn=="PHE" || rn=="TYR" || rn=="TRP" || rn=="HIS") mchi = 2;
	for(int i=rstart[ir]+5; i<rstart[ir+1]; i++){
		int *iz = atoms[i].getint();
		if(iz[3] != 0) continue;
		idx[nchi++] = i;
		if(nchi >= mchi) break;
	}
	return nchi;
}
void Protein::resMutation(int ir, string rn){
	int id = resDefine(rn);
	mutAA(ir, id);
}
void Protein::mutAA(int ir, int id){
	if(ir<0 || ir>=nres) die("Error residue No.: %d %d\n", ir, nres);
	if(id < 0) die("Error resid: %d\n", id);
	int id0 = resid[ir];
	if(id == id0) return;
	int na0 = rstart[ir+1] - rstart[ir]; //restypes[id0].getnatm();
// save old CHI
	initneib1_sc(ir);
	for(int i=4; i<na0; i++){
		initAtomint(rstart[ir]+i); calAtomzcrd(rstart[ir]+i);
	}
	int nchi0, idx[4]; double chi0[4];
	nchi0 = getReschi(ir, idx);
	for(int i=0; i<nchi0; i++) chi0[i] = atoms.at(idx[i]).gettorsion();
//
	resid[ir] = id; rnam[ir] = restypes.at(id).getname();
	int na = restypes[id].getnatm();
	int nd = na0 - na;
	if(nd > 0) eraseAtoms(rstart[ir+1]-nd, nd);
	else if(nd < 0) insertAtoms(rstart[ir+1], -nd);
//
	vector <Atomtype> ap = restypes[id].getatypes();
	for(int i=0; i<na; i++){
		int ia = rstart[ir] + i;
		string an = ap.at(i).getname();
		int it = ap.at(i).gettype();
		atoms[ia].setname(an, it);
	}
	initneib1_sc(ir);
//
	for(int i=4; i<na; i++){
		initAtomint(rstart[ir]+i);
		initAtomzcrd0(rstart[ir]+i);
	}
	if(rstart[ir]+4 < rstart[ir+1]) {
		if(! atoms.at(rstart[ir]+4).filled()) makexyz1(rstart[ir]+4); // CB atom
	}
	int nchi = getReschi(ir, idx);
	if(nchi <= 0) return;
// set the rotamer
	if(nchi <= nchi0) {		// for residues Big->small
		for(int i=0; i<nchi; i++) atoms.at(idx[i]).settorsion( chi0[i] );
		makexyz_part(rstart[ir]+5, rstart[ir+1]);
//		for(int i=5; i<na; i++) atoms.at(rstart[ir]+i).xyzatm1(atoms);
		return;
	}
	mutRotamer(ir, 0);
}
int Protein::mutRotamer(int ir, int irot){
// from standard rotamers
	double phi = atoms.at(rstart[ir]+2).gettorsion();
	double psi = atoms.at(rstart[ir]+3).gettorsion() + 180.;
	int iphi = int(fmod(phi+365., 360) / 10.);
	int ipsi = int(fmod(psi+365., 360) / 10.);
// set the chi angle
	initRotamers();
	int id = resid[ir];
	Rotamers *rp = &rotamers_aa[id][iphi][ipsi];
	if(irot >= rp->getnrot()) return -1;
	int idx[4];
	int nchi = getReschi(ir, idx);
	vector<int> &ichi0 = rp -> getchi(irot);
	if(ichi0.size()/2 != nchi) die("Diff chi: %s %d %d %d\n", rnam[ir].c_str(), ir, nchi, ichi0.size()/2);
	for(int i=0; i<nchi; i++) atoms.at(idx[i]).settorsion( ichi0[i]/10. );
//	for(int i=rstart[ir]+5; i<rstart[ir+1]; i++) atoms.at(i).xyzatm1(atoms);
	makexyz_part(rstart[ir]+5, rstart[ir+1]);
	return 0;
}
bool Protein::refill(){
	initneib();
	for(int i=0; i<natm; i++){
		if( atoms.at(i).filled() ) continue;
		initAtomint(i);
		initAtomzcrd0(i);
		if(! makexyz1(i)) return false;
	}
	return true;
}
double Protein::calResclash(int ir){
	double nclash=0.;
	for(int i=0; i<natm; i++){
		if(rseq[i] != ir) continue;
		for(int j=0; j<natm; j++){
			if(ir<0 && i<j) continue;
			if(rseq[i] == rseq[j]) continue;
			double r = sqrt( distance2(i, j) );
			if(r < 2.) nclash ++;
		}
	}
	return nclash;
}
//
Atom::Atom(int is, int it, string an){
	aseq=is; type=it; name=an;
	bfill=false;
}
void Atom::addneib(int ia){
	for(int i=0; i<ineib.size(); i++){
		if(ia == ineib[i]) {die("Duplicate neib: %d", ia);}
		else if(ia > ineib[i]) continue;
		ineib.insert(ineib.begin()+i, ia);
		return;
	}
	ineib.push_back(ia);
}
bool Protein::makexyz1(int ia){
	double *xa, *xb, *xc;
	xa = xb = xc = NULL; atoms[ia].setfill(false); 
	int *iz = atoms[ia].getint();
	int i0=iz[0], i1=iz[1], i2=iz[2];
	if(i0 >= 0){
		if(! atoms[i0].filled()) return false;
		xa = atoms[i0].getx();
	}
	if(i1 >= 0){
		if(! atoms[i1].filled()) return false;
		xb = atoms[i1].getx();
	}
	if(iz[2] >= 0){
		if(! atoms[i2].filled()) return false;
		xc = atoms[i2].getx();
	}
	double *x0 = atoms[ia].getx();
	double *ag = atoms[ia].getzcrd();
	bool bfill = xyzatm(x0, xa, xb, xc, ag, iz[3]);
	atoms[ia].setfill(bfill); 
	return bfill;
}
