#include "protein3.h"
#include "rotamer.h"

bool Atom3::makexyz1(){
	double *xp[3];
	if(iz[0].size() < 3) die("initint need before makexyz");
	bfill = false;
	for(int m=0; m<3; m++) xp[m] = NULL;
	for(int m=0; m<3; m++){
		if(iz[0][m] < 0) break;
		int ir = iz[0][m], ia = iz[1][m];
		Atom3 *ap = &getatom(iz[0][m], iz[1][m]);
		if(! ap->isfilled()) return bfill;
		xp[m] = ap->getx();
	}
	bfill = xyzatm(x, xp[0], xp[1], xp[2], ag, chiral);
	return bfill;
}
inline int ATOMSEQ_CMP(int ir1, int ia1, int ir2, int ia2){
	if(ir1 > ir2) return 1;
	if(ir1 < ir2) return -1;
	return ia1 - ia2;
}
void Atom3::addneib(int ir, int ia){
	for(int i=0; i<neib[0].size(); i++){
		int icmp = ATOMSEQ_CMP(ir, ia, neib[0][i], neib[1][i]);
		if(icmp == 0){
//			cout<<neib[0].at(-1)<<endl;
			die("Duplicate neib: %d %d\n", ir,ia);
		}
		else if(icmp > 0) continue;
		neib[0].insert(neib[0].begin()+i, ir);
		neib[1].insert(neib[1].begin()+i, ia);
		return;
	}
	neib[0].push_back(ir);
	neib[1].push_back(ia);
}
void Atom3::delete_neib1(int ir){	// delete neib within resi.
	vector<int> idx;
	for(int i=0; i<neib[0].size(); i++){
		if(neib[0][i] != ir) idx.push_back(i);
	}
	for(int i=0; i<idx.size(); i++){
		neib[0][i] = neib[0][idx[i]];
		neib[1][i] = neib[1][idx[i]];
	}
	neib[0].resize(idx.size());
	neib[1].resize(idx.size());
}
void Atom3::prtneib(int ir0){
	cout<<name<<": ";
	for(int m=0; m<neib[0].size(); m++){
		if(neib[0][m] > ir0) cout<<'+';
		else if(neib[0][m] < ir0) cout<<'-';
		cout<<getatom(neib[0][m], neib[1][m]).getname()<<' ';
	}
	cout<<endl;
}
void Atom3::prtint(){
	cout<<name<<": ";
	for(int m=0; m<iz[0].size(); m++){
		if(iz[0][m] < 0) break;
		cout<<iz[0][m]<<getatom(iz[0][m], iz[1][m]).getname()<<' ';
	}
	cout<<chiral<<endl;
}
double Atom3::calztor(){
	if(iz[0][0]<0 || iz[0][1]<0 || iz[0][2]<0) return 360.;
	double *x1 = getatom(iz[0][0], iz[1][0]).getx();
	double *x2 = getatom(iz[0][1], iz[1][1]).getx();
	double *x3 = getatom(iz[0][2], iz[1][2]).getx();
	return torsion(x, x1, x2, x3);
}
void Atom3::calzcrd(){
	if(iz[0][0] < 0) return;
	Atom3 *ap1 = &getatom(iz[0][0], iz[1][0]);
	ag[0] = sqrt( distance2(x, ap1->getx()) );
	if(iz[0][1] < 0) return;
	Atom3 *ap2 = &getatom(iz[0][1], iz[1][1]);
	ag[1] = angle(x, ap1->getx(), ap2->getx());
	if(iz[0][2] < 0) return;
	Atom3 *ap3 = &getatom(iz[0][2], iz[1][2]);
	ag[2] = torsion(x, ap1->getx(), ap2->getx(), ap3->getx());
	if(chiral != 0) {
		chiral = 1;
		if(ag[2] > 0) chiral = -1;
		ag[2] = angle(x, ap1->getx(), ap3->getx());
	}
}
Atom3 &Atom3::getatom(int ir, int ia){
	return res_all->at(ir).getatoms().at(ia);
}
void Atom3::initint(int ir0, int aseq){
	iz[0].resize(3); iz[1].resize(3);
	for(int m=0; m<3; m++) iz[0][m] = iz[1][m] = -1;
	chiral = -100;
	if(neib[0].size() <= 0) return;
	int icmp = ATOMSEQ_CMP(neib[0][0], neib[1][0], ir0, aseq);
	if(icmp >= 0) return;
	iz[0][0] = neib[0][0]; iz[1][0] = neib[1][0];
	Atom3 *ap1 = &getatom(iz[0][0], iz[1][0]);
// in PRO, to replace N
	if(res_all->at(iz[0][0]).getname()=="PRO" && name=="CD"){
		ap1 = &getatom(neib[0][1], neib[1][1]);
		assert(ap1->getname() == "CG");
		iz[0][0] = neib[0][1]; iz[1][0] = neib[1][1];
	}
//
	vector<int> *np1 = ap1->getneib();
	icmp = ATOMSEQ_CMP(np1[0][0], np1[1][0], ir0, aseq);
	if(icmp >= 0) return;
	iz[0][1] = np1[0][0]; iz[1][1] = np1[1][0];
	icmp = ATOMSEQ_CMP(np1[0][1], np1[1][1], ir0, aseq);
	if(icmp < 0){
		iz[0][2] = np1[0][1]; iz[1][2] = np1[1][1];
		chiral = 1;
	} else {
		Atom3 *ap2 = &getatom(iz[0][1], iz[1][1]);
		vector<int> *np2 = ap2->getneib();
		icmp = ATOMSEQ_CMP(np2[0][0], np2[1][0], iz[0][0], iz[1][0]);
		if(icmp < 0){
			iz[0][2] = np2[0][0]; iz[1][2] = np2[1][0];
			chiral = 0;
		}else{
			assert(icmp == 0);
			if(np2[0].size() <= 1) return;
			icmp = ATOMSEQ_CMP(np2[0][1], np2[1][1], iz[0][0], iz[1][0]);
			if(icmp >= 0) return;
			iz[0][2] = np2[0][1]; iz[1][2] = np2[1][1];
			chiral = 0;
		}
	}
}
//
//
void Residue::initAtomzcrd0(){
	if(restpl == NULL) restpl = new Restpl("");
	int idx[4]; double ag[3];
	for(int i=0; i<atoms.size(); i++){
		idx[0] = restpl->findatom_tpl(name, atoms[i].getname());
		int chiral;
		vector<int> *iz = atoms[i].getiz(chiral);
		for(int m=0; m<3; m++){
			idx[m+1] = -1;
			if(iz[0][m] < 0) break; 
			Atom3 *ap = &atoms[i].getatom(iz[0][m], iz[1][m]);
			idx[m+1] = restpl->findatom_tpl(name, ap->getname());
		}
		restpl->calint_tpl(idx, chiral, ag);
		atoms[i].setchiral(chiral);
		atoms[i].setzcrd(ag);
	}
}
void Residue::init(string rn0, int idx){
	int rid = resDefine(rn0);
	init(rid, idx);
	name = rn0;
}
// require to call "initRestypes_pro" first
void Residue::init(int rid, int idx){
	type = rid; rseq = idx;
	atoms.clear();
	if(type < 0) return;
	if(type >= restypes.size()) die("Resid too big: %d\n", rid);
	name = restypes[type].getname();
	vector <Atomtype> &ap = restypes[type].getatypes();
	for(int i=0; i<ap.size(); i++){
		string an = ap[i].getname();
		int it = ap[i].gettype();
		atoms.push_back(Atom3(it, an));
	}
}
Rotamers &Residue::determine_Rotamers(){
	initRotamers();
	calphipsi();
	int iphi = int(fmod(phipsi[0]+365., 360) / 10.);
	int ipsi = int(fmod(phipsi[1]+365., 360) / 10.);
	return rotamers_aa[type][iphi][ipsi];
}
int Residue::mutRotamer(int irot, double &prot){
	Rotamers *rp = &determine_Rotamers();
	if(irot >= rp->getnrot()) return -1;
	if(idx_chi.size() == 0) calidx_chi();
	vector<int> &ichi0 = rp->getchi(irot, prot);
	vector<double> chi0;
	int nchi0 = ichi0.size() / 2;
	for(int i=0; i<nchi0; i++) chi0.push_back(ichi0[i]/10.);
	mutRotamer(chi0);
	return 0;
}
void Residue::mutRotamer(const vector<double> &chi0){
	calidx_chi();
	if(chi0.size() != idx_chi.size()){
		die("Diff chi: %s %d %d\n", name.c_str(), idx_chi.size(), chi0.size());
	}
	for(int i=0; i<idx_chi.size(); i++) atoms.at(idx_chi[i]).setztor(chi0[i]);
	for(int i=4; i<atoms.size(); i++) atoms[i].makexyz1();
}
int Residue::mutRotamer(int irot){
	double p;
	return mutRotamer(irot, p);
}
void Residue::mutation(int id){
	if(id < 0) die("Error resid: %d\n", id);
	if(type == id) return;
	int type0 = type;
	type = id; name = restypes[id].getname();
	idx_chi.clear();
//
	vector <Atomtype> ap = restypes[id].getatypes();
	atoms.resize(ap.size());
	for(int i=0; i<ap.size(); i++){
		string an = ap.at(i).getname();
		int it = ap[i].gettype();
		atoms[i].setname(an, it);
		atoms[i].getresidues() = atoms[0].getresidues();
//		atoms[i].setresidues(atoms[0].getresidues());
	}
//
	initneib_sc();
	initAtomint(4);
	initAtomzcrd0();
//
	calidx_chi();
	if(idx_chi.size() > 0){
	/* set the rotamer
		if(nchi <= nchi0) {		// for residues Big->small
			for(int i=0; i<nchi; i++) atoms.at(idx[i]).settorsion( chi0[i] );
			makexyz_part(rstart[ir]+5, rstart[ir+1]);
	//		for(int i=5; i<na; i++) atoms.at(rstart[ir]+i).xyzatm1(atoms);
			return;
		}*/
		mutRotamer(0);
	} else {
		for(int i=4; i<atoms.size(); i++) atoms[i].makexyz1();
	}
}
void Residue::calphipsi(){
	if(phipsi.size() == 2) return;
	phipsi.clear(); phipsi.resize(2, 360.);
	if(rseq > 0) phipsi[0] = atoms.at(2).calztor();
	phipsi[1] = atoms.at(3).calztor() + 180.;
}
double Residue::calProt_nat(){
	calphipsi();
	int iphi = int(fmod(phipsi[0]+365., 360) / 10.);
	int ipsi = int(fmod(phipsi[1]+365., 360) / 10.);
	Rotamers &rp = rotamers_aa[type][iphi][ipsi];
	assert(rp.getchi(0).size() == idx_chi.size()*2);
//
	double p0, prot=0., dag_min;
	int nchi = idx_chi.size();
	for(int irot=0; irot<rp.getnrot(); irot++){
		vector<int> &ichi0 = rp.getchi(irot, p0);
		double dag = 0.;
		for(int i=0; i<idx_chi.size(); i++){
			double dt = atoms.at(idx_chi[i]).getztor() - ichi0[i]/10.;
			dt = fmod(fabs(dt), 360.);
			if(dt > 180.) dt = 360.- dt;
			dag += dt;
		}
		if(irot==0 || dag_min > dag) {
			prot = p0; dag_min = dag;
		}
	}
	return prot;
}
void Residue::getchi(vector<double> &chi){
	chi.clear();
	for(int i=0; i<idx_chi.size(); i++){
		double ag = atoms[idx_chi[i]].getztor();
		chi.push_back(ag);
	}
}
void Residue::calchi(){
	if(idx_chi.size() == 0) calidx_chi();
	for(int i=0; i<idx_chi.size(); i++){
		double ag = atoms.at(idx_chi[i]).calztor();
		atoms[idx_chi[i]].setztor(ag);
	}
}
void Residue::calidx_chi(){
	idx_chi.clear();
	int mchi=4;
	if(name=="PHE" || name=="TYR" || name=="TRP" || name=="HIS") mchi = 2;
	for(int i=5; i<atoms.size(); i++){
		int it = atoms[i].getchiral();
		if(it != 0) continue;
		idx_chi.push_back(i);
		if(idx_chi.size() >= mchi) break;
	}
}
int Residue::findatom(string an, bool bprint=0){
	for(int i=0; i<atoms.size(); i++){
		if(atoms.at(i).getname() == an) return i;
	}
	if(DEBUG>0 || bprint){
		fprintf(stderr, "findatom error: %s-%s\n", name.c_str(), an.c_str());
	}
	return -1;
}
inline string fix_pdbanam(const string rn, const string &an){
	if(an == "OXT") return "";
	else if(an[0] == 'H') return "";
	else if(an=="CD" && rn=="ILE") return "CD1";
	else if(rn.size() < 3){		// DNA or RNA
		if(an == "OP1") return "O1P";
		else if(an == "OP2") return "O2P";
		else if(an == "OP3") return false;
		else if(rn=="DT" && an=="C7") return "C5M";
	}
	return an;
}
bool Residue::initAtom_pdb(const string &line){
	char an1 = line[12], an2 = line[16], an0[6];
	sscanf(line.substr(12,4).c_str(), "%s", an0);
	string an = an0;
	if(an1=='1' || an1=='2' || an1=='3'){
		an = string(an0+1) + an1;
	} else if(an2=='1' || an2=='2' || an2=='3'){
		an += an2;
	}
	string ainf = line.substr(6, 10);
//
	if(type < 0){
		atoms.push_back(Atom3(-1,an));
		atoms.back().setx( line.substr(30, 24) );
		atoms.back().setinfo( ainf );
		atoms.back().setline( line );
		return true;
	}
//
	string an_fix = fix_pdbanam(name, an);
	if(an_fix == "") return false;
	int it = findatom(an_fix);
	if(it >= 0) {
		atoms[it].setx(line.substr(30, 24));
		atoms[it].setinfo(ainf);
		atoms[it].setline( line );
		return true;
	}
	return false;
}
void Residue::initneib_sc(){
	for(int i=0; i<atoms.size(); i++) atoms[i].delete_neib1(rseq);
	Restype *rp = &restypes.at(type);
	vector <string> *neib0 = rp->getneibs();
	for(int i=0; i<neib0->size(); i++){
		int i0 = findatom(neib0[0][i]);
		int i1 = findatom(neib0[1][i]);
		if(DEBUG > 0) printf("%d %d\n", i0, i1);
		if(i0<0 || i1<0) continue;
		atoms.at(i0).addneib(rseq, i1);
		atoms.at(i1).addneib(rseq, i0);
	}
}
void Residue::initAtomint(int i0=0){
	for(int i=i0; i<atoms.size(); i++) atoms[i].initint(rseq, i);
}
//
//
Protein3::Protein3(const string fn, bool bh){
	reset();
	bhet = bh;
	if(fn == "--"){
		name = "STDIN";
		rdpdb(stdin); return;
	}
	name = fn;
	FILE *fp=openfile(fn, "r");
	rdpdb(fp); fclose(fp);
}
void Protein3::reset(){
	name = chainID = "";
	chnseq.clear(); chnstart.clear();
	natm = nres = nchain = 0;
	residues.clear();
}
void Protein3::rdpdb(const string fn){
	FILE *fp = openfile(fn);
	rdpdb(fp);
	fclose(fp);
}
void Protein3::rdpdb(FILE *fp){
	char str[121], an[5], chn0 = '\n';
	string line, s0, s1, status0;
// init standard amino acids
	if(restypes.size()==0) initRestypes_pro();
	if(chnstart.size()>0 && chnstart.back() >= nres) chnstart.pop_back();
//
	bool bres = false;
	while(fgets(str,120,fp) != NULL){
		bool bh = 0;
		line = str;
		if(line.substr(0,3) == "END") break;
		if(line.substr(0,3) == "TER") {status0 = "TER"; continue;}
//		if(line.substr(0,3) == "TER") nchain ++;
		if(bhet && line.substr(0,6) == "HETATM") bh = 1;
		else if(line.substr(0,4) != "ATOM") continue;
//
		if(line.substr(17,3) == "HOH") continue;
// chain
		if(str[21] != chn0 || status0 == "TER") {
			nchain ++; chn0 = str[21];
			chainID += chn0;
			chnstart.push_back(nres);
		}
// residue
		s1 = line.substr(17,10); int ires=strtol(str+22,NULL,10);
		if(s0 != s1){
			if(DEBUG>0 && !chnseq.empty() && nchain==chnseq.back()){
				fprintf(stderr, "Warning, Residue not adjacent: %s --- %s\n", s0.c_str(), s1.c_str());
			}
			s0 = s1;
			string rn3 = line.substr(17,3);
			if(bh){
				string rn30 = rn3; rn3 = convert_hetres(rn30);
				if(rn3 != rn30) bh=0;
			}
			if(bh) bres = addRes_het(rn3);
			else bres = addRes(rn3);
			if(! bres) continue;
			residues.back().rseq0 = ires;
			residues.back().setresidues(&residues);
			residues.back().getinfo() = s1;
			chnseq.push_back(nchain - 1);
			nres ++;
		}
		if(bres){
			bool batom = residues.back().initAtom_pdb(line);
			if(batom) natm ++;
		}
		/*if(! bres) pdbidx.push_back(-1);
		else {
			if(! initAtom(line)) pdbidx.push_back(-1);
		}*/
		status0 = "ATOM";
	}
	if(DEBUG > 0){
		fprintf(stderr, "nres/natm = %d/%d\n", nres, natm);
	}
//	fclose(fp);
//	vector<Amino>::iterator pos;
//	for(pos=AAs.begin(); pos!=AAs.end(); pos++) pos->print();
}
bool Protein3::addRes_het(const string rn0){
	char rn[4]; sscanf(rn0.c_str(), "%s", rn);
	Residue *rp = new Residue();
	rp -> init_het(rn, nres);
	residues.push_back(*rp);
	return true;
}
bool Protein3::addRes(const string rn0){
	char rn[4]; sscanf(rn0.c_str(), "%s", rn);
	int rid = resDefine(rn);
	if(rid < 0) return false;
	residues.push_back(Residue(rid, nres));
	return true;
}
void Protein3::wrpdb(const string fn){
	file_backup(fn);
	FILE *fp = openfile(fn, "w");
	wrpdb(fp); fclose(fp);
}
void Protein3::wrpdb(FILE *fp){
	char fmt1[] = "ATOM%7d  %-4s%3s %c%4d    %8.3f%8.3f%8.3f\n";
	char chn0 = '\n';
	int na = 0;
	for(int ir=0; ir<residues.size(); ir++){
		vector<Atom3> &ap = residues[ir].getatoms();
		string &rnam = residues[ir].getname();
		for(int i=0; i<ap.size(); i++){
			if(! ap[i].isfilled()) continue;
			double *x = ap[i].getx();
			char chn = 'A';
			if(chnseq.size()>ir && chainID.size()>chnseq[ir]){
				chn = chainID.at(chnseq[ir]);
				if(chn0!='\n' && chn!=chn0) fprintf(fp, "TER\n");
				chn0 = chn;
			} else if(DEBUG){
				fprintf(stderr, "not defined chain name for RES%d: %d %d\n", ir, chnseq.size(), chainID.size());
			}
			int rid = ir + 1;
//			if(seq0.size() > ir) rid = seq0[ir];
			fprintf(fp, fmt1, ++na, ap[i].getname().c_str(), rnam.c_str(),
					chn, rid, *x, x[1], x[2]);
		}
	}
	fprintf(fp, "TER\nEND\n");
}
void Protein3::initneib(){
	for(int ir=0; ir<residues.size(); ir++){
		vector<Atom3> &ap = residues[ir].getatoms();
		for(int i=0; i<ap.size(); i++) ap[i].clear_neib();
	}
//
	for(int ir=0; ir<residues.size(); ir++) residues[ir].initneib_sc();
	for(int ir=1; ir<residues.size(); ir++){
		Residue &rp = residues[ir];
		int id = rp.gettype();
		if(id < 0) continue;
		if(chnseq[ir] != chnseq[ir-1]) continue;
		if(rp.getname().size() == 3){		// amino acid
			Atom3 &ap1 = residues[ir-1].getatoms()[2];  // C of res. ir-1
			Atom3 &ap2 = rp.getatoms()[0];		// N of resi. ir
			ap1.addneib(ir, 0);
			ap2.addneib(ir-1, 2);
			if(ap2.getname()!="N" || ap1.getname()!="C"){
				die("wrong AA neib: %s %s -- %s %s\n", residues[ir-1].getname_str(), ap1.getname_str(), rp.getname_str(), ap2.getname_str());
			}
		} else {
			Atom3 &ap1 = residues[ir-1].getatoms()[8];  // O3' for base ir-1
			Atom3 &ap2 = rp.getatoms()[0];		// P for base
			ap1.addneib(ir, 0);
			ap2.addneib(ir-1, 8);
			if(ap2.getname()!="P" || ap1.getname()!="O3'"){
				die("wrong Base neib: %s %s -- %s %s\n", residues[ir-1].getname_str(), ap1.getname_str(), rp.getname_str(), ap2.getname_str());
			}
		}
	}
}
void Protein3::initAtomint(){
	for(int ir=0; ir<residues.size(); ir++) residues[ir].initAtomint(0);
}
void Protein3::initAtomzcrd0(){
	for(int ir=0; ir<residues.size(); ir++) residues[ir].initAtomzcrd0();
}
void Protein3::resMutation(int ir, string rn){
	int id = resDefine(rn);
	if(id < 0) die("unknown residue for mut: %s\n", rn.c_str());
	residues[ir].mutation(id);
}
/*
bool Protein3::iscomplete(){
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
// Att! aseq afterwards not changed
void Protein3::addCB_GLY(int ir){
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
bool Protein3::refill(){
	initneib();
	for(int i=0; i<natm; i++){
		if( atoms.at(i).filled() ) continue;
		initAtomint(i);
		initAtomzcrd0(i);
		if(! makexyz1(i)) return false;
	}
	return true;
}
*/
