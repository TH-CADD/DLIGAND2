#include "misc.h"

double probe=1.5; 
void surface(int na, vector<double> &xv, vector<double> &radius, vector<double> &weight, vector<double> &area){
	void sort2(int n, double *list, int *key);
	const int maxarc=300;
	const double pi = 3.141592653589793238;
	double force[na], total;
	double deal, decl, arcf[maxarc], arci[maxarc];
	int narc, ider[maxarc];
	double thec, ccsq, bsqk, bsql, kent[maxarc], dsql, ther[maxarc];
	int skip[na], omit[maxarc];
	double wght, area_tot, risq[maxarc];
	int iout;
	double kout[maxarc], rrsq, xysq, b[maxarc],p;
	int i, j, k, l, m, i3;
	double r[na], s, t, v, x[na], y[na], z[na], delta, 
		dtkal, dtlal, dtkcl;
	int intag[maxarc];
	double exang, dtlcl;
	int moved;
	double therk, b1[maxarc];
	int komit;
	double rmove, risqk, risql;
	int sign_yder[maxarc];
	double t1, t2, rplus, delta2, wxlsq;
	int intag1[maxarc];
	double cc, bg[maxarc];
	int ib, jb;
	double bk, dk, gi;
	int ii;
	double gl, gk;
	int mi, in, io, ni, ir;
	double td, tb, ti, tf, xc[maxarc], yc[maxarc], zc[maxarc], arclen, rr, 
		tr, ex[maxarc], tt, lt[maxarc], xr, tx, ty, tz, yr, cosine, zr, ux[maxarc],
		 uy[maxarc], uz[maxarc], ri[maxarc], gr[maxarc], arcsum, xc1[maxarc], rminus, 
		yc1[maxarc], zc1[maxarc], tk1, tk2, tr2, bgl, dax, day, daz, the, rcn, 
		rik, bsq[maxarc], eps;
	int key[maxarc];
	double dsq[maxarc], txb, tyb, axx, axy, axz, ayx, ayy, azx, azy, 
		azz, wxl, uxl, uyl, uzl, txr, tyr, txk, tyk, tzk, txl, tyl, tzl;
	int top, lforce;
	double pid2, bsq1[maxarc], dsq1[maxarc], pix2, pix4, faca, facb, 
		facc, rrx2, gaca, gacb;
	double *darea;

	area.resize(na);
	lforce = false;
	if(lforce) darea=new double[3*na];
	else darea=NULL;
// initialization
	iout=6; area_tot=0;
	double eff_x=1.;
	for(i=0;i<na;i++){
		i3=i*3;
		x[i]=xv[i3] * eff_x;
		y[i]=xv[i3+1] * eff_x;
		z[i]=xv[i3+2] * eff_x;
	}
// zero the area and derivatives, and set the sphere radii
	total=0;
	for(i=0;i<na;i++){
		i3=i*3;
		area[i]=0;
		if(lforce) darea[i3]=darea[i3+1]=darea[i3+2]=0;
		r[i]=radius[i];
		if(r[i]>=0) r[i]+=probe;
	}
// set pi multiples, overlap criterion and tolerances
	pix2=2*pi; pix4=4*pi; pid2=pi/2.0;
	delta=1.0E-8; delta2=delta*delta;
	eps=rmove=1.0E-8;
	for(i=0;i<maxarc;i++){
		ider[i]=0;
		sign_yder[i]=0;
	}
// set the "skip" array to exclude all inactive atoms
// that do not overlap any of the current active atoms
	for(i=0;i<na;i++) skip[i]=1;
	for(i=0;i<na;i++){
		xr=x[i]; yr=y[i]; zr=z[i]; rr=r[i];
		for(k=0;k<na;k++){
			rplus=rr+r[k];
			ccsq=(xr-x[k])*(xr-x[k])+(yr-y[k])*(yr-y[k])+(zr-z[k])*(zr-z[k]);
			if(ccsq <= rplus) skip[k]=0;
		}
	}
// compute the area and derivatives of current "ir" sphere
	for(ir=0;ir<na;ir++){
		if(skip[ir]) continue;
		xr=x[ir]; yr=y[ir]; zr=z[ir];
		rr=r[ir];
		rrx2=rr*2.0;
		rrsq=rr * rr;
		wght=weight[ir];
		moved=0;
// initialize some counters and sums for the "ir" sphere
l10: ;
		io=ib=jb=0;
		arclen=exang=0.0;
// test each sphere to see if it overlaps the "ir" sphere
		for(i=0;i<na;i++){
			if(i==ir) continue;
			rplus=rr+r[i];
			tx=x[i]-xr;
			if(fabs(tx)>=rplus) continue;
			ty=y[i]-yr;
			if(fabs(ty)>=rplus) continue;
			tz=z[i]-zr;
			if(fabs(tz)>=rplus) continue;
// check for overlap of spheres by testing center to
// center distance against sum and difference of radii
			xysq=tx*tx+ty*ty;
			if(xysq<delta2){
				tx=delta; ty=0; xysq=delta2;
			}
			ccsq=xysq+tz*tz;
			cc=sqrt(ccsq);
			if(rplus-cc<=delta) continue;
			rminus=rr-r[i];
// check for a completely buried "ir" sphere
			if(cc-fabs(rminus)<=delta){
				if(rminus<0.0) goto l180;
				continue;
			}
// calculate overlap parameters between "i" and "ir" sphere
			if(io>=maxarc){
				printf("SURFACE  --  Increase the Value of maxarc\n");
			}
			xc1[io]=tx;
			yc1[io]=ty;
			zc1[io]=tz;
			dsq1[io]=xysq;
			bsq1[io]=ccsq;
			b1[io]=cc;
			gr[io]=(ccsq+rplus*rminus)/(rrx2*b1[io]);
			intag1[io]=i;
			io++;
		}
// case where no other spheres overlap the current sphere
		if(io==0){
			area[ir]=pix4;
			goto l160;
		}
// case where only one sphere overlaps the current sphere
		if(io==1){
			k=0;
			txk=xc1[0];
			tyk=yc1[0];
			tzk=zc1[0];
			bsqk=bsq1[0];
			bk=b1[0];
			intag[0]=intag1[0];
			arcsum=pix2;
			ib++;
			arclen+=gr[k]*arcsum;
			if(lforce){
			if(moved==0){
				in=intag[k];
				t1=arcsum*rrsq*(bsqk-rrsq+r[in]*r[in])/(rrx2*bsqk*bk);
				darea[3*ir] -= txk*t1*wght;
				darea[3*ir+1] -= tyk*t1*wght;
				darea[3*ir+2] -= tzk*t1*wght;
				darea[3*in] += txk*t1*wght;
				darea[3*in+1] += tyk*t1*wght;
				darea[3*in+2] += tzk*t1*wght;
			}}
			goto l150;
		}
// general case where more than one sphere intersects the
// current sphere; sort intersecting spheres by their degree
// of overlap with the current main sphere
		sort2(io,gr,key);
		for(i=0;i<io;i++){
			k=key[i];
			intag[i]=intag1[k];
			xc[i]=xc1[k];
			yc[i]=yc1[k];
			zc[i]=zc1[k];
			dsq[i]=dsq1[k];
			b[i]=b1[k];
			bsq[i]=bsq1[k];
			omit[i]=0;
		}
// radius of the each circle on the surface of the "ir" sphere
		for(i=0;i<io;i++){
			gi=gr[i]*rr;
			bg[i]=b[i]*gi;
			risq[i]=rrsq-gi*gi;
			ri[i]=sqrt(risq[i]);
			ther[i]=pid2-asin(min(1.0,max(-1.0,gr[i])));
		}
// find boundary of inaccessible area on "ir" sphere
		for(k=0;k<io-1;k++){
			if(omit[k]==0){
				txk=xc[k];
				tyk=yc[k];
				tzk=zc[k];
				bk = b[k];
				therk=ther[k];
				for(j=k+1;j<io;j++){
					if(omit[j]) goto l60;
// check to see if J circle is intersecting K circle
// get distance between circle centers and sum of radii
					cc=(txk*xc[j]+tyk*yc[j]+tzk*zc[j])/(bk*b[j]);
					cc = acos(min(1.0,max(-1.0,cc)));
					td = therk + ther[j];
// check to see if circles enclose separate regions
					if(cc>=td) goto l60;
// check for circle J completely inside circle K
					if(cc+ther[j]<therk) goto l40;
// check for circles essentially parallel
					if(cc>delta) goto l50;
l40: ;
					omit[j]=1;
					goto l60;
// check for "ir" sphere completely buried
l50: ;
					if(pix2-cc<=td) goto l180;
l60: ;
				}
			}
		}
// find T value of circle intersections
		for(k=0;k<io;k++){
			if(omit[k]) goto l110;
			komit=omit[k];
			omit[k]=1;
			narc=0;
			top=0;
			txk=xc[k];
			tyk=yc[k];
			tzk=zc[k];
			dk=sqrt(dsq[k]);
			bsqk=bsq[k];
			bk=b[k];
			gk=gr[k]*rr;
			risqk = risq[k];
			rik = ri[k];
			therk = ther[k];
// rotation matrix elements
			t1 = tzk / (bk*dk);
			axx = txk * t1;
			axy = tyk * t1;
			axz = dk / bk;
			ayx = tyk / dk;
			ayy = txk / dk;
			azx = txk / bk;
			azy = tyk / bk;
			azz = tzk / bk;
			for(l=0;l<io;l++){
				if(omit[l]==0){
					txl=xc[l];
					tyl=yc[l];
					tzl=zc[l];
// rotate spheres so K vector colinear with z-axis
					uxl = txl*axx + tyl*axy - tzl*axz;
					uyl = tyl*ayy - txl*ayx;
					uzl = txl*azx + tyl*azy + tzl*azz;
					cosine = min(1.0,max(-1.0,uzl/b[l]));
					if(acos(cosine)<therk+ther[l]){
						dsql=uxl*uxl+uyl*uyl;
						tb = uzl*gk - bg[l];
						txb = uxl * tb;
						tyb = uyl * tb;
						td = rik * dsql;
						tr2 = risqk*dsql - tb*tb;
						tr2 = max(eps,tr2);
						tr = sqrt(tr2);
						txr = uxl * tr;
						tyr = uyl * tr;
// get T values of intersection for K circle
						tb = (txb+tyr) / td;
						tb = min(1.0,max(-1.0,tb));
						tk1 = acos(tb);
						if (tyb-txr < 0.0)  tk1 = pix2 - tk1;
						tb = (txb-tyr) / td;
						tb = min(1.0,max(-1.,tb));
						tk2 = acos(tb);
						if (tyb+txr < 0.0)  tk2 = pix2 - tk2;
						thec = (rrsq*uzl-gk*bg[l]) / (rik*ri[l]*b[l]);
						if (fabs(thec) < 1.0){
							the = -acos(thec);
						} else if(thec >=1.0){
							the=0.0;
						}else{
							the=-pi;
						}
// see if "tk1" is entry or exit point; check t=0 point;
// "ti" is exit point, "tf" is entry point
						cosine=min(1., max(-1.,(uzl*gk-uxl*rik)/(b[l]*rr)));
						if((acos(cosine)-ther[l])*(tk2-tk1)<=0.0){
							ti=tk2; tf=tk1;
						} else{
							ti=tk2; tf=tk1;
						}
						if(tf<=ti){
							arcf[narc]=tf;
							arci[narc]=0.0;
							tf=pix2;
							lt[narc]=l;
							ex[narc]=the;
							top=1;
							narc++;
						}
						arcf[narc]=tf;
						arci[narc]=ti;
						lt[narc]=l;
						ex[narc]=the;
						ux[l]=uxl;
						uy[l]=uyl;
						uz[l]=uzl;
						narc++;
						if(narc>=maxarc){
							printf("SURFACE  --  Increase the Value of maxarc");
						}
					}
				}
			}
			omit[k]=komit;
// special case; K circle without intersections
			if(narc<=0) goto l90;
// general case; sum up arclength and set connectivity code
			sort2(narc,arci,key);
			arcsum=arci[0];
			mi=key[0];
			t=arcf[mi];
			ni=mi;
			if(narc>1){
				for(j=1;j<narc;j++){
					m=key[j];
					if(t<arci[j]){
						arcsum+=arci[j]-t;
						exang+=ex[ni];
						l=(int)lt[ni];
						ider[l]++;
						sign_yder[l]++;
						kent[jb]=l*maxarc+k;
						l=(int)lt[m];
						ider[l]++;
						sign_yder[l]--;
						kout[jb] = maxarc*k + l;
						jb++;
						if(jb>=maxarc){
							printf("SURFACE  --  Increase the Value of maxarc\n");
						}
					}
					tt=arcf[m];
					if(tt>=t){t=tt;ni=m;}
				}
			}
			arcsum+=pix2-t;
			if(top==0){
				exang+=ex[ni];
				l=(int)lt[ni];
				ider[l]++;
				sign_yder[l]++;
				kent[jb]=l*maxarc + k;
				l=(int)lt[mi];
				ider[l]++;
				sign_yder[l]--;
				kout[jb]=k*maxarc+l;
				jb++;
			}
// calculate the surface area derivatives
			if(lforce){
			for(l=0;l<io;l++){
				if(ider[l]!=0){
					rcn=ider[l]*rrsq;
					ider[l]=0;
					uzl=uz[l];
					gl=gr[l]*rr;
					bgl=bg[l];
					bsql=bsq[l];
					risql=risq[l];
					wxlsq=bsql-uzl*uzl;
					wxl=sqrt(wxlsq);
					p=bgl-gk*uzl;
					v=risqk*wxlsq-p*p;
					v=max(eps,v);
					v=sqrt(v);
					t1=rr*(gk*(bgl-bsql)+uzl*(bgl-rrsq))/ (v*risql*bsql);
					deal = -wxl*t1;
					decl = -uzl*t1 - rr/v;
					dtkal = (wxlsq-p) / (wxl*v);
					dtkcl = (uzl-gk) / v;
					s = gk*b[l] - gl*uzl;
					t1 = 2.0*gk - uzl;
					t2 = rrsq - bgl;
					dtlal = -(risql*wxlsq*b[l]*t1-s*(wxlsq*t2+risql*bsql))/
						(risql*wxl*bsql*v);
					dtlcl = -(risql*b[l]*(uzl*t1-bgl)-uzl*t2*s)/ (risql*bsql*v);
					gaca = rcn * (deal-(gk*dtkal-gl*dtlal)/rr) / wxl;
					gacb = (gk-uzl*gl/b[l]) * sign_yder[l] * rr / wxlsq;
					sign_yder[l] = 0;
					if(moved==0){
						faca=ux[l]*gaca-uy[l]*gacb;
						facb=uy[l]*gaca+ux[l]*gacb;
						facc = rcn * (decl-(gk*dtkcl-gl*dtlcl)/rr);
						dax = axx*faca - ayx*facb + azx*facc;
						day = axy*faca + ayy*facb + azy*facc;
						daz = azz*facc - axz*faca;
						in=intag[l];
						darea[3*ir]+=dax*wght;
						darea[3*ir+1]+=day*wght;
						darea[3*ir+2]+=daz*wght;
						darea[3*in]-=dax*wght;
						darea[3*in+1]-=day*wght;
						darea[3*in+2]-=daz*wght;
					}
				}
			}}
			goto l100;
l90: ;
			arcsum=pix2;
			ib++;
l100: ;
			arclen+=gr[k]*arcsum;
			if(lforce){
			if(moved==0){
				in=intag[k];
				t1=arcsum*rrsq*(bsqk-rrsq+r[in]*r[in]) / (rrx2*bsqk*bk);
				darea[3*ir]-=txk*t1*wght;
				darea[3*ir+1]-=tyk*t1*wght;
				darea[3*ir+2]-=tzk*t1*wght;
				darea[3*in]+=txk*t1*wght;
				darea[3*in+1]+=tyk*t1*wght;
				darea[3*in+2]+=tzk*t1*wght;
			}}
l110: ;
		}
		if(arclen==0.0) goto l180;
		if(jb==0) goto l150;
// find number of independent boundaries and check connectivity
		j=0;
		for(k=0;k<jb;k++){
			if(kout[k]!=-1){
				i=k;
l120: ;
				m=(int)kout[i];
				kout[i]=-1;
				j++;
				for(ii=0;ii<jb;ii++){
					if(m==kent[ii]){
						if(ii==k){
							ib++;
							if(j==jb) goto l150;
							goto l130;
						}
						i=ii;
						goto l120;
					}
				}
l130: ;
			}
		}
		ib++;
// attempt to fix connectivity error by moving atom slightly
		if(moved){
			printf("SURFACE  --  Connectivity Error at atom %d\n",ir+1);
		}else{
			moved=1;
			xr+=rmove;
			yr+=rmove;
			zr+=rmove;
			goto l10;
		}
// form the accessible area for the current atom
l150: ;
		area[ir]=ib*pix2+exang+arclen;
		area[ir]=fmod(area[ir],pix4);
l160: ;
		area[ir]*=rrsq;
// attempt to fix negative area by moving atom slightly
		if(area[ir]<0.0){
			if(moved){
				printf("SURFACE  --  Negative Area at atom %d\n",ir);
			}else{
				moved=1;
				xr+=rmove;
				yr+=rmove;
				zr+=rmove;
				goto l10;
			}
		}
// weight the accessible area by the scale factor
		area_tot+=area[ir];
		area[ir]*=wght;
		total+=area[ir];
//		printf("%d %f\n",ir+1,area[ir]);
l180: ;
	}
	if(lforce){
		for(i=0;i<na;i++){
			force[i*3]=darea[3*i];
			force[i*3+1]=darea[3*i+1];
			force[i*3+2]=darea[3*i+2];
		}
	}
}

/*	 ############################################################## 
	 ##                                                          ##
	 ##  subroutine sort2  --  heapsort of real array with keys  ##
	 ##                                                          ##
	 ##############################################################


	 "sort2" takes an input list of reals and sorts it
	 into ascending order using the Heapsort algorithm;
	 it also returns a key into the original ordering
*/
void sort2(int n, double *list, int *key) {

	int keys, i, j, k, index;
	double lists;
// initialize index into the original ordering 
	for(i=0;i<n;i++){
		key[i]=i;
	}
// perform the heapsort of the input list 
	k = n / 2;
	index = n - 1;
	while(n > 1) {
		if (k > 0) {
			k--;
			lists = list[k];
			keys = key[k];
		} else {
			lists = list[index];
			keys = key[index];
			list[index] = list[0];
			key[index] = key[0];
			index--;
			if (index <= 0) {
				list[0] = lists;
				key[0] = keys;
				return;
			}
		}
		i = k;
		j = k*2+1;
		while(j <= index) {
			if (j<index && list[j]<list[j+1]) j++;
			if (lists < list[j]) {
				list[i] = list[j];
				key[i] = key[j];
				i = j;
				j = 2*j + 1;
			} else {
				j = index + 1;
			}
		}
		list[i] = lists;
		key[i] = keys;
	}
	return;
}
//
//
int main(int argc, char *argv[]){
	if(argc < 2) die("usage: RUN1 xyzrn [-atom|-total]");
	bool btotal = findargs(argc, argv, "-total") > 0;
	bool batom = findargs(argc, argv, "-atom") > 0;
	
	string fpdb="";
	for(int i=1; i<argc; i++){
		if(argv[i][0] != '-') {fpdb = argv[i]; break;}
	}
	FILE *fp = openfile(fpdb);
	string line, ss[9];
	vector<double> xlist, rlist, wlist, area;
	vector<string> list_ch, list_anam;
	while(getline1(fp, line)){
		int nt = str2dat(line, ss);
		if (nt < 1) continue;
		for(int m=0; m<3; m++) xlist.push_back(strtod(ss[m].c_str(), NULL));
		rlist.push_back(strtod(ss[3].c_str(), NULL));
		list_ch.push_back(ss[4]);
		list_anam.push_back(ss[5]);
		wlist.push_back(1.);
	}
	int natm = rlist.size();
	surface(natm, xlist, rlist, wlist, area);

	cout<<"#total area: "<<sum(area.size(), area)<<endl;
	if(btotal) exit(0);
	if(batom){
		for(int i=0; i<natm; i++) printf("%s %s %g\n", list_anam[i].c_str(), list_ch[i].c_str(), area[i]);
	} else {
		int nr=0; vector<string> list_rnam, ss2;
		vector<double> list_rasa;
		vector<int> list_natm;
		for(int i=0; i<natm; i++){
			split_str(list_anam[i], '_', ss2);
			assert(ss2.size() == 3);
			string rn1 = ss2[1]+"_"+ss2[2] + "_" + list_ch[i];
			if(i==0 || rn1!=list_rnam[nr-1]) {
				list_rnam.push_back(rn1);
				list_rasa.push_back(area[i]);
				list_natm.push_back(1);
				nr ++;
			} else {
				list_rasa[nr-1] += area[i];
				list_natm[nr-1] ++;
			}
		}
		for(int i=0; i<list_rnam.size(); i++) printf("%s %g %d\n", list_rnam[i].c_str(), list_rasa[i], list_natm[i]);
	}
}
