#ifndef _MACRO_Y
#define _MACRO_Y

// the dieX() series of macros
/*#define  die0() { fprintf(stderr, \
		"\nERROR ON LINE %d OF FILE %s, EXITING!\n", \
		__LINE__, __FILE__); exit(1); }*/
#define die1(a) { fprintf(stderr, a); die0(); }
#define die2(a, b) { fprintf(stderr, a, b); die0(); }
// memory allocation macros
#define build_1D_array(a, t, x) { \
		if((a=new t[x]) == NULL) \
		die2("\nallocation of %ldB failed\n", sizeof(t)*((long int)(x))); }
//           printf("# calloc_1D_array on line %d of file %s: %lx from %p\n",
//		  __LINE__, __FILE__, sizeof(t)*((long int)(x)), a); }
#define build_2D_array(a,t,x,y) { \
		build_1D_array(a,t*,x); build_1D_array(a[0],t,(x)*(y)); \
		for(int _i=1; _i<(x); _i++) a[_i]=a[0]+_i*(y); }
#define build_3D_array(a,t,x,y,z) { \
		build_2D_array(a,t*,x,y); \
		build_1D_array(a[0][0],t,(x)*(y)*(z)); \
		for(int _x=0;_x<(x);_x++) \
		for(int _y=0;_y<(y);_y++) a[_x][_y] = a[0][0] + (z)*(_y + (y)*_x); }
#define  delete_1D_array(a) { \
		if(a!=NULL){ delete [] a; a=NULL; } }
#define  delete_2D_array(a) { \
		if(a!=NULL){ delete [] a[0]; delete [] a; a = NULL;} }
#define  delete_3D_array(a) { \
		if(a!=NULL){ delete [] a[0][0]; delete [] a[0]; delete [] a; a=NULL;} }

#endif
