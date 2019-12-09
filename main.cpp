/*

	TSApredictor

	Marc Meidlinger, 2019

	Fast routine to predict a refinement level at which interior cells can be detected using
	the Figueiredo et. al cell mapping/interval arithmetics algorithm: 

	"Images of Julia sets that you can trust"
	by Luiz-Henrique de Figueiredo, Diego Nehab, Jofge Stolfi, Joao Batista Oliveira-
	from 2013

	"Rigorous bounds for polynomial Julia sets"
	by Luiz-Henrique de Figueiredo, Diego Nehab, Jofge Stolfi, Joao Batista Oliveira-

*/

#include "math.h"
#include "stdio.h"
#include "stdint.h"
#include "string.h"
#include "quadmath.h"
#include "time.h"

typedef uint8_t BYTE;
typedef uint32_t DDBYTE;
typedef DDBYTE *PDDBYTE;

// chunk size dependend on operating system
// win32 => 512 MB
//const uint64_t CHUNKSIZE=( (uint64_t)1 << 29 );
// win64 => 1 GB
const uint64_t CHUNKSIZE=( (uint64_t)1 << 30 );

// used floating type
// comment out or in what is needed
#define _DOUBLE
//#define _LONGDOUBLE
//#define _QUADMATH

const int64_t DENOM225=( (int64_t)1 << 25);
const int32_t MAXZEROS=1024;
const int32_t SHIFTPERDDBYTE=5; // 32 bit-wide

const BYTE SQUARE_GRAY=0;
const DDBYTE ALL32GRAY=0;

const BYTE SQUARE_POTW=1;
const DDBYTE DDBYTEMAX=0b11111111111111111111111111111111;
const DDBYTE ALL32POTW=DDBYTEMAX;

#ifdef _QUADMATH
typedef __float128 NTYP;
const char NNTYPSTR[]="qd";
#endif

#ifdef  _LONGDOUBLE
typedef long double NTYP;
const char NNTYPSTR[]="ld";
#endif

#ifdef _DOUBLE
typedef double NTYP;
const char NNTYPSTR[]="d";
#endif

// two orbit points are identical if:
const NTYP ZEROEPSILON=1E-15;
// the number 0 is
const NTYP COEFFZEROLIMIT=1E-40;
// maximal degree for class Polynom
const int32_t MAXDEGREE=32;

// array manager
const int32_t MAXPTR=2048;

enum { 
	FUNC_Z2C=0,FUNC_Z2AZC=1,FUNC_Z3AZC=2,
	FUNC_Z4AZC=3,FUNC_Z5AZC=4,FUNC_Z6AZC=5,
	
	FUNCANZ
};

const char funcname[][32] = {
	"Z2C","Z2AZC","Z3AZC","Z4AZC","Z5AZC","Z6AZC"
};


// structs

struct PlaneRect {
	NTYP x0,x1,y0,y1;
};

struct ScreenRect {
	int32_t x0,x1,y0,y1;
};

struct ArrayDDByteManager {
	DDBYTE* current;
	int32_t allocatedIdx,freeFromIdx,allocatePerBlock;
	PDDBYTE ptr[MAXPTR];
	int32_t anzptr;
	
	ArrayDDByteManager ();
	virtual ~ArrayDDByteManager ();
	void FreeAll(void);
	PDDBYTE getMemory(const int32_t);
};

struct Complex {
	NTYP re,im;
	
	Complex(const NTYP r_=0, const NTYP i_=0);
		
	friend Complex operator+(const Complex&,const Complex&);
	friend Complex operator-(const Complex&,const Complex&);
	Complex operator*(const Complex);
	Complex& operator=(const Complex&);
	friend Complex operator/(const Complex&,const Complex&);
	
	NTYP norm(void);
	NTYP normQ(void);
	void output(FILE*);
};

struct Polynom {
	int32_t grad;
	Complex coeff[MAXDEGREE];
	int8_t coeffnull[MAXDEGREE];
	
	Polynom();
	void setCoeff(const int32_t,const Complex);
	void setCoeff(const int32_t,const NTYP,const NTYP);
	void setCoeff(const int32_t,const NTYP);
	virtual void eval_arg_f(const Complex,Complex&); // über Horner
	void clearCoeff(void);
	void output(FILE*);
};

struct PeriodicPoint {
	Complex pp;
	int32_t mem0,mem1;
	int32_t y0,y1;
};

struct Root {
	Complex attractor;
	PeriodicPoint* cycle;
	PlaneRect ps_basinrect;
	int interiorfound;
	int cyclelen;
	int cyclenumber;
	double multiplier;
	
	void clear(void);
};


// globals

FILE *flog=NULL;
DDBYTE _STARTWITH=ALL32GRAY;
int MAXIT=25000;
int LEVEL0=8,LEVEL1=24;
Root zero[MAXZEROS];
char COMPUTECOMMANDLINE[4096];
int fctr=1;
int _ENCLOSEMENTWIDTH=128;
int nbr_of_cp;
int PERIODICLEN0=-1,PERIODICLEN1=-1;
PlaneRect local;
void (*getBoundingBoxfA)(PlaneRect&,PlaneRect&) = NULL;
Polynom fkt;
int _FUNC;
NTYP seedC0re,seedC1re,seedC0im,seedC1im; 
NTYP FAKTORAre,FAKTORAim;
NTYP scaleRangePerPixel,scalePixelPerRange;
NTYP COMPLETE0,COMPLETE1;
Complex cplxA,cplxC;

// forward declarations

inline int32_t scrcoord_as_lowerleft(const NTYP);
inline NTYP minimumD(const NTYP,const NTYP);
inline NTYP maximumD(const NTYP,const NTYP);
inline NTYP minimumD(const NTYP,const NTYP,const NTYP,const NTYP);
inline NTYP maximumD(const NTYP,const NTYP,const NTYP,const NTYP);


#define LOGMSG(TT) \
{\
	fprintf(flog,TT); fflush(flog);\
	printf(TT);\
}

#define LOGMSG2(TT,AA) \
{\
	fprintf(flog,TT,AA); fflush(flog);\
	printf(TT,AA);\
}

#define LOGMSG3(TT,AA,BB) \
{\
	fprintf(flog,TT,AA,BB); fflush(flog);\
	printf(TT,AA,BB);\
}

#define LOGMSG4(TT,AA,BB,CC) \
{\
	fprintf(flog,TT,AA,BB,CC); fflush(flog);\
	printf(TT,AA,BB,CC);\
}

#define LOGMSG5(TT,AA,BB,CC,DD) \
{\
	fprintf(flog,TT,AA,BB,CC,DD); fflush(flog);\
	printf(TT,AA,BB,CC,DD);\
}

#define SQUARE_LIES_ENTIRELY_IN_LOCAL(BBX) \
	(\
		(local.x0 <= BBX.x0) &&\
		(BBX.x1 <= local.x1) &&\
		(local.y0 <= BBX.y0) &&\
		(BBX.y1 <= local.y1)\
	)

#define SQUARE_LIES_ENTIRELY_IN_COMPLETE(BBX) \
	(\
		(COMPLETE0 <= BBX.x0) &&\
		(BBX.x1 <= COMPLETE1) &&\
		(COMPLETE0 <= BBX.y0) &&\
		(BBX.y1 <= COMPLETE1)\
	)


// non-struct function

// functions

inline int32_t scrcoord_as_lowerleft(const NTYP a) {
	// calculating the screen coordinte of the pixel that contains the coordinate
	// if the coordinate lies on an edge/corner (and belongs to more than one pixel)
	// the pixel where it lies on the left,bottom edge/corner is returned
	#ifdef _QUADMATH
	return (int)floorq( (a - COMPLETE0) * scalePixelPerRange );
	#else
	return (int)floor( (a - COMPLETE0) * scalePixelPerRange );
	#endif
}

inline NTYP maximumD(const NTYP a,const NTYP b,const NTYP c,const NTYP d) {
	NTYP m=a;
	if (b > m) m=b;
	if (c > m) m=c;
	if (d > m) m=d;
	return m;
}

inline NTYP minimumD(const NTYP a,const NTYP b,const NTYP c,const NTYP d) {
	NTYP m=a;
	if (b < m) m=b;
	if (c < m) m=c;
	if (d < m) m=d;
	return m;
}

inline NTYP minimumD(const NTYP a,const NTYP b) {
	if (a < b) return a;
	return b;
}

inline NTYP maximumD(const NTYP a,const NTYP b) {
	if (a > b) return a;
	return b;
}

// z^2+c
void getBoundingBoxfA_z2c(PlaneRect& A,PlaneRect& fA) {
	fA.x0=minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)+seedC0re;
	fA.x1=maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)+seedC1re;
	fA.y0=2*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)+seedC0im;
	fA.y1=2*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)+seedC1im;
}

// z^2+A*z+c
void getBoundingBoxfA_z2azc(PlaneRect& A,PlaneRect& fA) {
	fA.x0=seedC0re+minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1)+minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1)-maximumD(A.y0*A.y0,A.y1*A.y1);
	fA.x1=seedC1re+maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1)+maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1)-minimumD(A.y0*A.y0,A.y1*A.y1);
	fA.y0=seedC0im+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+2*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1);
	fA.y1=seedC1im+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+2*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1);
}

// z^3+A*z+c
void getBoundingBoxfA_z3azc(PlaneRect& A,PlaneRect& fA) {
	fA.x0=minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+A.x0*A.x0*A.x0-(3*maximumD(A.x0*minimumD(A.y0*A.y0,A.y1*A.y1),A.x0*maximumD(A.y0*A.y0,A.y1*A.y1),A.x1*minimumD(A.y0*A.y0,A.y1*A.y1),A.x1*maximumD(A.y0*A.y0,A.y1*A.y1)))+seedC0re;
	fA.x1=maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+A.x1*A.x1*A.x1-(3*minimumD(A.x0*minimumD(A.y0*A.y0,A.y1*A.y1),A.x0*maximumD(A.y0*A.y0,A.y1*A.y1),A.x1*minimumD(A.y0*A.y0,A.y1*A.y1),A.x1*maximumD(A.y0*A.y0,A.y1*A.y1)))+seedC1re;
	fA.y0=minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+3*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0,A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y1)-(A.y1*A.y1*A.y1)+seedC0im;
	fA.y1=maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+3*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0,A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y1)-(A.y0*A.y0*A.y0)+seedC1im;
}

// z^4+A*z+c
void getBoundingBoxfA_z4azc(PlaneRect& A,PlaneRect& fA) {
	fA.x0=minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1,FAKTORAre*A.x0,FAKTORAre*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1,FAKTORAim*A.y0,FAKTORAim*A.y1)+minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-(6*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1)))+minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+seedC0re;
	fA.x1=maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1,FAKTORAre*A.x0,FAKTORAre*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1,FAKTORAim*A.y0,FAKTORAim*A.y1)+maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-(6*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1)))+maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+seedC1re;
	fA.y0=minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1,FAKTORAre*A.y0,FAKTORAre*A.y1)+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1,FAKTORAim*A.x0,FAKTORAim*A.x1)+4*minimumD((A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1)*A.y1)-(4*maximumD(A.x0*(A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1)))+seedC0im;
	fA.y1=maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1,FAKTORAre*A.y0,FAKTORAre*A.y1)+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1,FAKTORAim*A.x0,FAKTORAim*A.x1)+4*maximumD((A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1)*A.y1)-(4*minimumD(A.x0*(A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1)))+seedC1im;
}

// z^5+A*z+c
void getBoundingBoxfA_z5azc(PlaneRect& A,PlaneRect& fA) {
	fA.x0=minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1,FAKTORAre*A.x0,FAKTORAre*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1,FAKTORAim*A.y0,FAKTORAim*A.y1)+A.x0*A.x0*A.x0*A.x0*A.x0-(2*(5*maximumD((A.x0*A.x0*A.x0)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x0*A.x0*A.x0)*maximumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+5*minimumD(A.x0*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x0*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1))+seedC0re;
	fA.x1=maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1,FAKTORAre*A.x0,FAKTORAre*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1,FAKTORAim*A.y0,FAKTORAim*A.y1)+A.x1*A.x1*A.x1*A.x1*A.x1-(2*(5*minimumD((A.x0*A.x0*A.x0)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x0*A.x0*A.x0)*maximumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+5*maximumD(A.x0*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x0*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1))+seedC1re;
	fA.y0=minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1,FAKTORAre*A.y0,FAKTORAre*A.y1)+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1,FAKTORAim*A.x0,FAKTORAim*A.x1)+5*minimumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1)-(2*(5*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1))))+A.y0*A.y0*A.y0*A.y0*A.y0+seedC0im;
	fA.y1=maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1,FAKTORAre*A.y0,FAKTORAre*A.y1)+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1,FAKTORAim*A.x0,FAKTORAim*A.x1)+5*maximumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1)-(2*(5*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1))))+A.y1*A.y1*A.y1*A.y1*A.y1+seedC1im;
}

// z^5+c*z+A
// c kann IA sein, A ist fix
void getBoundingBoxfA_z5cza(PlaneRect& A,PlaneRect& fA) {
	fA.x0=minimumD(seedC0re*A.x0,seedC0re*A.x1,seedC1re*A.x0,seedC1re*A.x1)-maximumD(seedC0im*A.y0,seedC0im*A.y1,seedC1im*A.y0,seedC1im*A.y1)+A.x0*A.x0*A.x0*A.x0*A.x0-(2*(5*maximumD((A.x0*A.x0*A.x0)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x0*A.x0*A.x0)*maximumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+5*minimumD(A.x0*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x0*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1))+FAKTORAre;
	fA.x1=maximumD(seedC0re*A.x0,seedC0re*A.x1,seedC1re*A.x0,seedC1re*A.x1)-minimumD(seedC0im*A.y0,seedC0im*A.y1,seedC1im*A.y0,seedC1im*A.y1)+A.x1*A.x1*A.x1*A.x1*A.x1-(2*(5*minimumD((A.x0*A.x0*A.x0)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x0*A.x0*A.x0)*maximumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+5*maximumD(A.x0*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x0*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1))+FAKTORAre;
	fA.y0=minimumD(seedC0re*A.y0,seedC0re*A.y1,seedC1re*A.y0,seedC1re*A.y1)+minimumD(seedC0im*A.x0,seedC0im*A.x1,seedC1im*A.x0,seedC1im*A.x1)+5*minimumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1)-(2*(5*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1))))+A.y0*A.y0*A.y0*A.y0*A.y0+FAKTORAim;
	fA.y1=maximumD(seedC0re*A.y0,seedC0re*A.y1,seedC1re*A.y0,seedC1re*A.y1)+maximumD(seedC0im*A.x0,seedC0im*A.x1,seedC1im*A.x0,seedC1im*A.x1)+5*maximumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1)-(2*(5*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1))))+A.y1*A.y1*A.y1*A.y1*A.y1+FAKTORAim;
}

// z^6+A*z+c
void getBoundingBoxfA_z6azc(PlaneRect& A,PlaneRect& fA) {
	fA.x0=seedC0re+minimumD(FAKTORAre*A.x0,FAKTORAre*A.x1)-maximumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+minimumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(3*(5*maximumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+3*(5*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)))-maximumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1);
	fA.x1=seedC1re+maximumD(FAKTORAre*A.x0,FAKTORAre*A.x1)-minimumD(FAKTORAim*A.y0,FAKTORAim*A.y1)+maximumD(A.x0*A.x0*A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1*A.x1*A.x1)-(3*(5*minimumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+3*(5*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)))-minimumD(A.y0*A.y0*A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1*A.y1*A.y1);
	fA.y0=minimumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+minimumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+6*minimumD((A.x0*A.x0*A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1*A.x1*A.x1)*A.y1)-(4*(5*maximumD((A.x0*A.x0*A.x0)*(A.y0*A.y0*A.y0),(A.x0*A.x0*A.x0)*(A.y1*A.y1*A.y1),(A.x1*A.x1*A.x1)*(A.y0*A.y0*A.y0),(A.x1*A.x1*A.x1)*(A.y1*A.y1*A.y1))))+6*minimumD(A.x0*(A.y0*A.y0*A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1*A.y1*A.y1))+seedC0im;
	fA.y1=maximumD(FAKTORAre*A.y0,FAKTORAre*A.y1)+maximumD(FAKTORAim*A.x0,FAKTORAim*A.x1)+6*maximumD((A.x0*A.x0*A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1*A.x1*A.x1)*A.y1)-(4*(5*minimumD((A.x0*A.x0*A.x0)*(A.y0*A.y0*A.y0),(A.x0*A.x0*A.x0)*(A.y1*A.y1*A.y1),(A.x1*A.x1*A.x1)*(A.y0*A.y0*A.y0),(A.x1*A.x1*A.x1)*(A.y1*A.y1*A.y1))))+6*maximumD(A.x0*(A.y0*A.y0*A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1*A.y1*A.y1))+seedC1im;
}

char* seedCstr225(char* erg) {
	sprintf(erg,"c_ia_%I64d_%I64d_x_%I64d_%I64d",
		(int64_t)floor(DENOM225*seedC0re),
		(int64_t)floor(DENOM225*seedC1re),
		(int64_t)floor(DENOM225*seedC0im),
		(int64_t)floor(DENOM225*seedC1im)
	);
	return erg;
}

char* FAKTORAstr225(char* erg) {
	sprintf(erg,"A_%I64d_%I64d",
		(int64_t)floor(DENOM225*FAKTORAre),
		(int64_t)floor(DENOM225*FAKTORAim)
	);
		
	return erg;
}

int getfuncidx(const char* s) {
	for(int32_t i=0;i<FUNCANZ;i++) {
		if (!strcmp(s,funcname[i])) return i;
	}
	
	return -1;
}

void setfunc(const int32_t afunc) {
	fkt.clearCoeff();
	COMPUTECOMMANDLINE[0]=0;

	switch (afunc) {
		case FUNC_Z3AZC: {
			getBoundingBoxfA=getBoundingBoxfA_z3azc;
			fkt.setCoeff(3,1);
			fkt.setCoeff(1,cplxA);
			fkt.setCoeff(0,cplxC);
			sprintf(COMPUTECOMMANDLINE,"func=z3azc c=%.20lg,%.20lg A=%.20lg,%.20lg cmd=period,-1",
				(double)cplxC.re,(double)cplxC.im,
				(double)cplxA.re,(double)cplxA.im
			);
			break;
		} 
		case FUNC_Z4AZC: {
			getBoundingBoxfA=getBoundingBoxfA_z4azc;
			fkt.setCoeff(4,1);
			fkt.setCoeff(1,cplxA);
			fkt.setCoeff(0,cplxC);
			sprintf(COMPUTECOMMANDLINE,"func=z4azc c=%.20lg,%.20lg A=%.20lg,%.20lg cmd=period,-1",
				(double)cplxC.re,(double)cplxC.im,
				(double)cplxA.re,(double)cplxA.im
			);
			break;
		}
		case FUNC_Z5AZC: {
			getBoundingBoxfA=getBoundingBoxfA_z5azc;
			fkt.setCoeff(5,1);
			fkt.setCoeff(1,cplxA);
			fkt.setCoeff(0,cplxC);
			sprintf(COMPUTECOMMANDLINE,"func=z5azc c=%.20lg,%.20lg A=%.20lg,%.20lg cmd=period,-1",
				(double)cplxC.re,(double)cplxC.im,
				(double)cplxA.re,(double)cplxA.im
			);
			break;
		}
		case FUNC_Z6AZC: {
			getBoundingBoxfA=getBoundingBoxfA_z6azc;
			fkt.setCoeff(6,1);
			fkt.setCoeff(1,cplxA);
			fkt.setCoeff(0,cplxC);
			sprintf(COMPUTECOMMANDLINE,"func=z6azc c=%.20lg,%.20lg A=%.20lg,%.20lg cmd=period,-1",
				(double)cplxC.re,(double)cplxC.im,
				(double)cplxA.re,(double)cplxA.im
			);
			break;
		}
		case FUNC_Z2AZC: {
			getBoundingBoxfA=getBoundingBoxfA_z2azc;
			fkt.setCoeff(2,1);
			fkt.setCoeff(1,cplxA);
			fkt.setCoeff(0,cplxC);
			sprintf(COMPUTECOMMANDLINE,"func=z2azc c=%.20lg,%.20lg A=%.20lg,%.20lg cmd=period,-1",
				(double)cplxC.re,(double)cplxC.im,
				(double)cplxA.re,(double)cplxA.im
			);
			break;
		}
		default: {
			getBoundingBoxfA=getBoundingBoxfA_z2c;
			fkt.setCoeff(2,1);
			fkt.setCoeff(0,cplxC);
			sprintf(COMPUTECOMMANDLINE,"func=z2c c=%.20lg,%.20lg cmd=period,-1",
				(double)cplxC.re,(double)cplxC.im
			);
			break;
		}
	} // switch
}

// struct ArrayByteManager

ArrayDDByteManager::ArrayDDByteManager() {
	current=NULL;
	allocatedIdx=0;
	freeFromIdx=-1;
	anzptr=0;
	double d=CHUNKSIZE; d /= sizeof(DDBYTE);
	allocatePerBlock=(int)floor(d);
}

void ArrayDDByteManager::FreeAll(void) {
	for(int32_t i=0;i<anzptr;i++) {
		delete[] ptr[i];
	}
	current=NULL;
	anzptr=0;
}

ArrayDDByteManager::~ArrayDDByteManager() {
	FreeAll();
}

PDDBYTE ArrayDDByteManager::getMemory(const int32_t aanz) {
	if (anzptr >= (MAXPTR-8)) {
		LOGMSG("ArrayDDByteManager:: Zu wenig Speicher\n");
		exit(99);
	}
	if (
		(!current) ||
		((freeFromIdx + aanz + 2) >= allocatedIdx)
	) {
		printf("x");
		ptr[anzptr]=current=new DDBYTE[allocatePerBlock];
		anzptr++;
		if (!current) {
			LOGMSG("Memory-Fehler. ArrayByteManager.\n");
			exit(99);
		}
		freeFromIdx=0;
		allocatedIdx=allocatePerBlock;
	}
	
	PDDBYTE p=&current[freeFromIdx];
	freeFromIdx += aanz;
	return p;
}

char* upper(char* s) {
	if (!s) return NULL;
	for(int32_t i=(strlen(s)-1);i>=0;i--) {
		if ((s[i]>='a')&&(s[i]<='z')) s[i]=s[i]-'a'+'A';
	}

	return s;
}

int newton(
	Polynom& polyf,
	Polynom& polyabl,
	const Complex astart,
	Complex& erg) 
{
	Complex z=astart,f,fabl,zletzt,d;
	
	for(int32_t i=1;i<MAXIT;i++) {
		zletzt=z;
		// z(n+1)=zn-f(zn) / f'(zn);
		polyf.eval_arg_f(z,f);
		polyabl.eval_arg_f(z,fabl);
		z = z - f / fabl;
		d = z-zletzt;
		if (d.normQ() < ZEROEPSILON) {
			erg=z;
			return i;
		}
	}
	
	return 0;
}

int getNullstellenIdx(const Complex aw,const int32_t ait) {
	for(int32_t i=0;i<nbr_of_cp;i++) {
		Complex d=zero[i].attractor- aw;
		if (d.normQ() < ZEROEPSILON) {
			return i;
		}
	}
	if (nbr_of_cp > (MAXZEROS-8)) {
		LOGMSG("Error. Too many roots.\n");
		exit(99);
	}
	
	zero[nbr_of_cp].clear();
	zero[nbr_of_cp].attractor=aw;
	nbr_of_cp++;

	return (nbr_of_cp-1);
}

int getLagrange(Polynom& f) {
	// i.e. Douady
	double res=1.0;
	for(int32_t i=0;i<=f.grad;i++) {
		res += fkt.coeff[i].norm();
	}
	res /= fkt.coeff[fkt.grad].norm();
	int expo=(int)ceil(log(ceil(res))/log(2.0));
	return (1 << expo);
}

void ableitenFA(Polynom& infkt,Polynom& ergabl) {
	ergabl.clearCoeff();
	
	for(int32_t i=1;i<=infkt.grad;i++) {
		if (infkt.coeffnull[i]==0) {
			ergabl.setCoeff(
				i-1,
				Complex(i,0)*infkt.coeff[i]
			);
		}
	}
}

int ps_construct_critical_orbits(void) {
	double escapeQ=COMPLETE1*COMPLETE1;
	int cyclenumber=1;
	int returnvalue=0;

	Complex *orbit=new Complex[MAXIT];
	int orbitlen=0;
	
	Polynom polyabl;
	ableitenFA(fkt,polyabl);
	
	for(int32_t cp=0;cp<nbr_of_cp;cp++) {
		Complex z0=zero[cp].attractor;
		
		Complex zn=z0;
		int esc=0;
		orbitlen=0;
		for(int32_t i=0;i<MAXIT;i++) {
			orbit[orbitlen]=zn;
			orbitlen++;
			if (zn.normQ() > escapeQ) {
				esc=1;
				break;
			}
			Complex tmp;
			fkt.eval_arg_f(zn,tmp);
			zn=tmp;
		} // i
		
		if (esc>0) {
			zero[cp].clear();
			continue;
		}
		
		// bounded critical orbit
		// is it periodic ?
		int cyclestart=-1,cycleend=orbitlen-1;
		for(int32_t i=(orbitlen-2);i>=0;i--) {
			Complex d=orbit[i]-orbit[orbitlen-1];
			if (d.normQ() < ZEROEPSILON) {
				cyclestart=i;
				break;
			}
		}
		
		// not periodic
		if (cyclestart < 0) {
			// not periodic
			zero[cp].clear();
			continue;
		}
		
		// periodic orbit
		// start: cyclestart+1 .. cycleend
		// (orbit point at [cyclestart] is identical to [cycleend]
		
		// has an earlier analyzed critical point already
		// found that orbit ?
		// then discard this one here
		int found=0;
		for(int32_t cpprev=0;cpprev<cp;cpprev++) {
			for(int32_t k=0;k<zero[cpprev].cyclelen;k++) {
				Complex d=zero[cpprev].cycle[k].pp - orbit[cycleend];
				if (d.normQ() < ZEROEPSILON) {
					found=1;
					break;
				}
			}
			if (found>0) break;
		} // cpprev
		
		if (found>0) {
			// previously found cycle
			// this cp is marked as uninteresting
			zero[cp].clear();
			continue;
		}
		
		zero[cp].cyclelen=cycleend - (cyclestart+1) + 1;
		zero[cp].cycle=new PeriodicPoint[zero[cp].cyclelen + 8];
		Complex multiplier=Complex(1,0);
		for(int32_t i=(cyclestart+1);i<=cycleend;i++) {
			zero[cp].cycle[i-(cyclestart+1)].pp=orbit[i];
			Complex der;
			polyabl.eval_arg_f(orbit[i],der);
			multiplier = multiplier*der;
		}
		zero[cp].multiplier=multiplier.norm();
		zero[cp].cyclenumber=cyclenumber;
		cyclenumber++;
		if (zero[cp].multiplier > 1.00001) {
			// a little buffer
			// repelling
			zero[cp].cyclelen=0;
			delete[] zero[cp].cycle;
			zero[cp].cycle=NULL;
			zero[cp].interiorfound=0;
		}
		
		returnvalue++;
		
	} // cp
	
	delete[] orbit;
	
	return returnvalue;
}

void ps_find_critical_points(void) {
	// fkt: iterated function
	// cp: zeros of 1st order derivative
	Polynom fktforcp,ablforcp;
	ableitenFA(fkt,fktforcp);
	ableitenFA(fktforcp,ablforcp);
	
	double ESCAPEQ=COMPLETE1*COMPLETE1;

	int32_t LEN=1024;
	// 3 times ESCAPER: far away from the roots, so dynamics
	// here are tame: cvhannels to infinity as in paper byy Schleicher
	// about the universal set of Newton starting points
	int32_t SCRRE0=-3*ESCAPEQ;
	int32_t SCRIM0=-3*ESCAPEQ;
	int32_t SCRIM1= 3*ESCAPEQ;
	NTYP sk=SCRIM1-SCRIM0; sk /= LEN;
	
	Complex start;
	nbr_of_cp=0;
	#define SUCHE(X0,Y0,X1,Y1)\
	{\
	int yd=1,xd=1;\
	if ((Y0) > (Y1)) yd=-1; else if ((Y0) < (Y1)) yd=1; \
	if ((X0) > (X1)) xd=-1; else if ((X0) < (X1)) xd=1; \
	for(int32_t y=(Y0);y<=(Y1);y+=yd) {\
		start.im=y*sk + SCRIM0;\
		for(int32_t x=(X0);x<=(X1);x+=xd) {\
			start.re=x*sk + SCRRE0; \
			\
			Complex nulls;\
			int32_t it=newton(fktforcp,ablforcp,start,nulls);\
			if (it > 0) {\
				getNullstellenIdx(nulls,it);\
				if (nbr_of_cp >= fktforcp.grad) break;\
			} \
		} \
		if (nbr_of_cp >= fktforcp.grad) break;\
	} \
	}

	// searches all border pixels of the underlying image in the 3 times Lagrange estimate
	// (modified after Hubbard, Schleicher, Sutherland. How to find all roots of complex polynomials, 2001)
	
	SUCHE(0,0,0,LEN-1)
	if (nbr_of_cp < fktforcp.grad) {
		SUCHE(0,LEN-1,LEN-1,LEN-1)
	}
	if (nbr_of_cp < fktforcp.grad) {
		SUCHE(LEN-1,LEN-1,LEN-1,0)
	}
	if (nbr_of_cp < fktforcp.grad) {
		SUCHE(LEN-1,0,0,0)
	}
}

// Complex

void Complex::output(FILE* f) {
	fprintf(f,"%.20lg%+.20lgi",(double)re,(double)im);
}

Complex Complex::operator*(const Complex a) {
    Complex ret(a.re*re-a.im*im, a.im*re+a.re*im);

    return ret;
}

Complex::Complex(const NTYP r, const NTYP i) {
	re=r;
	im=i;
}

Complex operator/(const Complex& a,const Complex& b) {
    const NTYP n2=b.re*b.re+b.im*b.im;
    Complex ret(
		(a.re*b.re+a.im*b.im)/n2,
		(a.im*b.re-a.re*b.im) / n2
    );

    return ret;
}

NTYP Complex::norm(void) {
	return sqrt(re*re+im*im);
}

NTYP Complex::normQ(void) {
	return (re*re+im*im);
}

Complex& Complex::operator=(const Complex& c) {
	if (this != &c) {
		re=c.re;
		im=c.im;
	}
	
	return *this;
}

bool operator!=(const Complex a, const Complex b) {
    return ( (a.re != b.re) || (a.im != b.im) );
}

Complex operator+(const Complex& a, const Complex& b) {
    Complex ret(a.re+b.re, a.im+b.im);
    return ret;
}

Complex operator-(const Complex& a, const Complex& b) {
    Complex ret(a.re-b.re, a.im-b.im);
    return ret;
}

// Polynom

void Polynom::setCoeff(const int32_t aidx,const NTYP ar) {
	setCoeff(aidx,Complex(ar,0.0));
}

void Polynom::setCoeff(const int32_t aidx,const Complex ac) {
	coeff[aidx]=ac;
	if (coeff[aidx].normQ() < COEFFZEROLIMIT) coeffnull[aidx]=1;
	else {
		coeffnull[aidx]=0;
		if (aidx > grad) grad=aidx;
	}
}

void Polynom::clearCoeff(void) {
	for(int32_t i=0;i<MAXDEGREE;i++) {
		coeffnull[i]=1;
		coeff[i]=Complex(0.0,0.0);
	}
	grad=0;
}

void Polynom::output(FILE* f) {
	fprintf(f,"p(z)=");
	int32_t erster=1;
	for(int32_t i=grad;i>=0;i--) {
		if (coeffnull[i]==0) {
			if (erster) erster=0; else fprintf(f,"+");
			fprintf(f,"(");
			coeff[i].output(f);
			if (i>1) fprintf(f,")*z^%i",i);
			else if (i==1) fprintf(f,")*z");
			else fprintf(f,")");
		}
	}
	fprintf(f,"\n");
}

void Polynom::eval_arg_f(const Complex az,Complex& erg) {
	erg=coeff[grad];
	
	for(int32_t i=grad;i>0;i--) {
		erg = erg*az;
		erg = erg + coeff[i-1];
	}
}

Polynom::Polynom() {
	grad=0;
	clearCoeff();
}

void Polynom::setCoeff(const int32_t aidx,const NTYP ar,const NTYP ai) {
	setCoeff(aidx,Complex(ar,ai));
}

int cm_local(Root& onecycle,const DDBYTE startwith) {
	// a bounding box for ALL cyclic points
	// small rectangles around every cyclic point
	ScreenRect enclosementall;
	
	ArrayDDByteManager mgr;
	PDDBYTE *ispotwY=NULL;
	
	#define SET32_MY(MM,YY,FF32) \
	{\
		if (\
			( (MM) >= mem0 ) &&\
			( (MM) <= mem1 ) &&\
			( (YY) >= enclosementall.y0 ) &&\
			( (YY) <= enclosementall.y1 ) \
		) {\
			if (ispotwY[ (YY)-enclosementall.y0 ]) {\
				ispotwY[ (YY)-enclosementall.y0 ][MM-mem0] = FF32;\
			} else {\
				LOGMSG4("Error/mil. Set32 %i,%i,%i\n",MM,YY,FF32);\
				exit(99);\
			}\
		} else {\
			LOGMSG4("Error. Set32 %i,%i,%i\n",MM,YY,FF32);\
			exit(99);\
		}\
	}
	
	#define GET32_MY(MM,YY,ERG32) \
	{\
		ERG32=ALL32POTW;\
		if (\
			( (YY) >= enclosementall.y0) &&\
			( (YY) <= enclosementall.y1)\
		) {\
			if ( (ispotwY[ (YY)-enclosementall.y0 ]) &&\
				( (MM) >= mem0) &&\
				( (MM) <= mem1) \
			) {\
				ERG32=ispotwY[ (YY)-enclosementall.y0][ (MM)-mem0 ];\
			}\
		}\
	}

	// YY is ABSOLUTE screen coordinates
	// MM is ABSOLUTE mem position, i.e. *32 is pixel coordinate
	#define CELLCOLOR_XY(XX,YY,ERG) \
	{\
		ERG=SQUARE_POTW;\
		if (\
			( (XX) >= enclosementall.x0) &&\
			( (XX) <= enclosementall.x1) &&\
			( (YY) >= enclosementall.y0) &&\
			( (YY) <= enclosementall.y1) \
		) {\
			int bmem=(XX) >> SHIFTPERDDBYTE;\
			DDBYTE f32;\
			GET32_MY(bmem,YY,f32);\
			int bbit=(XX) % (1 << SHIFTPERDDBYTE);\
			ERG=(f32 >> bbit) & 0b1;\
		}\
	}
	
	int32_t interiorpresentat=0;
	onecycle.interiorfound=0;
	int8_t* ywithgray=NULL;
	for(int32_t REFINEMENT=LEVEL0;REFINEMENT<=LEVEL1;REFINEMENT++) {
		printf("\nchecking level %i ",REFINEMENT);
		int32_t SCREENWIDTH=(1 << REFINEMENT);
		int32_t MAXMEM=SCREENWIDTH >> SHIFTPERDDBYTE;
		scaleRangePerPixel=(COMPLETE1-COMPLETE0)/(double)SCREENWIDTH;
		scalePixelPerRange=(double)SCREENWIDTH/(COMPLETE1-COMPLETE0);
		enclosementall.x0=enclosementall.y0=SCREENWIDTH;
		enclosementall.x1=enclosementall.y1=0;
		for(int32_t k=0;k<onecycle.cyclelen;k++) {
			int32_t xx=scrcoord_as_lowerleft(onecycle.cycle[k].pp.re);
			int32_t yy=scrcoord_as_lowerleft(onecycle.cycle[k].pp.im);
			ScreenRect scr;
			scr.x0=xx-_ENCLOSEMENTWIDTH; 
			scr.x1=xx+_ENCLOSEMENTWIDTH;
			scr.y0=yy-_ENCLOSEMENTWIDTH;
			scr.y1=yy+_ENCLOSEMENTWIDTH;

			#define TRIM(WW) \
			{\
				if ( WW < 0) WW=0;\
				else if (WW >= SCREENWIDTH) WW=SCREENWIDTH-1;\
			}

			TRIM(scr.x0)
			TRIM(scr.x1)
			TRIM(scr.y0)
			TRIM(scr.y1)
			
			if (scr.x0 < enclosementall.x0) enclosementall.x0=scr.x0;
			if (scr.x1 > enclosementall.x1) enclosementall.x1=scr.x1;
			if (scr.y0 < enclosementall.y0) enclosementall.y0=scr.y0;
			if (scr.y1 > enclosementall.y1) enclosementall.y1=scr.y1;
			
			onecycle.cycle[k].mem0=scr.x0 >> SHIFTPERDDBYTE;
			onecycle.cycle[k].mem1=scr.x1 >> SHIFTPERDDBYTE;
			if (onecycle.cycle[k].mem1 >= MAXMEM) {
				LOGMSG("Implementation error. Maxmem reached\n");
				exit(99);
			}
			onecycle.cycle[k].y0=scr.y0;
			onecycle.cycle[k].y1=scr.y1;
		}
		
		int32_t mem0=enclosementall.x0 >> SHIFTPERDDBYTE;
		int32_t mem1=enclosementall.x1 >> SHIFTPERDDBYTE;
		if (mem1 >= MAXMEM) {
			LOGMSG("Implementation error. Maxmem/2 reached\n");
			exit(99);
		}
		
		// translate enclosementall into complex coordinates
		local.x0=onecycle.ps_basinrect.x0=enclosementall.x0*scaleRangePerPixel + COMPLETE0;
		local.x1=onecycle.ps_basinrect.x1=(enclosementall.x1+1)*scaleRangePerPixel + COMPLETE0;
		local.y0=onecycle.ps_basinrect.y0=enclosementall.y0*scaleRangePerPixel + COMPLETE0;
		local.y1=onecycle.ps_basinrect.y1=(enclosementall.y1+1)*scaleRangePerPixel + COMPLETE0;
		
		if (REFINEMENT==LEVEL0) printf("allocating ");
		
		int32_t LOCALLENY=enclosementall.y1-enclosementall.y0 + 1;
		int32_t LOCALLENX=mem1-mem0+1;

		if (ywithgray) delete[] ywithgray;
		ywithgray=new int8_t[LOCALLENY];
		for(int32_t y=0;y<LOCALLENY;y++) {
			ywithgray[y]=0;
		}

		// now go over the enclosements again and
		// set rows intersecting an enclosement to 1
		// so memory gets allocated and the row checked
		for(int32_t k=0;k<onecycle.cyclelen;k++) {
			for(int32_t y=onecycle.cycle[k].y0;y<=onecycle.cycle[k].y1;y++) {
				ywithgray[y-enclosementall.y0]=1;
			}
		}

		// allocate enough memory
		mgr.FreeAll();
		// now pointers in zeilenY are invalid
		// now delete that array
		if (ispotwY) delete[] ispotwY;
		ispotwY=new PDDBYTE[LOCALLENY];
		if (!ispotwY) {
			LOGMSG("Memory error. ispotwY\n");
			exit(99);
		}

		for(int32_t y=0;y<LOCALLENY;y++) {
			if (ywithgray[y]>0) {
				ispotwY[y]=mgr.getMemory(LOCALLENX);
				if (!ispotwY[y]) {
					LOGMSG("Memory error. ispotwY/2\n");
					exit(99);
				}
				// set ALL to startvalue
				for(int32_t m=0;m<LOCALLENX;m++) {
					ispotwY[y][m]=startwith;
				}
			} else {
				ispotwY[y]=NULL;
			}
		}
		
		// setting all enclosements of periodic points
		// to GRAY
		for(int32_t k=0;k<onecycle.cyclelen;k++) {
			for(int32_t y=onecycle.cycle[k].y0;y<=onecycle.cycle[k].y1;y++) {
				for(int32_t m=onecycle.cycle[k].mem0;m<=onecycle.cycle[k].mem1;m++) {
					SET32_MY(m,y,ALL32GRAY);
				}
			}
		}
		
		// going over the interesting (GRAY) regions
		// and try to propagate POTENTIALLY_WHITE until
		// no change occurs
		// anything still GRAY is then bounded
		// and means, BLACK will emerge at level REFINEMENT
		// in the fullTSA
		
		if (REFINEMENT==LEVEL0) printf(" analyzing ");
		else printf(" ");
		
		int8_t changed=1;
		PlaneRect A,bbxfA;
		int32_t noch0=256*(24-REFINEMENT);
		if (noch0<1) noch0=1;
		int32_t noch=1;
		
		while (changed>0) {
			changed=0;
			if ((--noch)<=0) {
				printf(".");
				noch=noch0;
			}
			
			// go over enclosementall, so points in overlapping
			// enclosement[k]'s will not be analyzed twice
			// in that round of the while-loop
			for(int32_t y=enclosementall.y0;y<=enclosementall.y1;y++) {
				if (ywithgray[y-enclosementall.y0]<=0) continue;
				
				int8_t graythere=0;
				
				A.y0=y*scaleRangePerPixel + COMPLETE0;
				A.y1=A.y0+scaleRangePerPixel;
				for(int32_t m=mem0;m<=mem1;m++) {
					DDBYTE ff;
					GET32_MY(m,y,ff);
					if (ff == ALL32POTW) continue;
					BYTE fchanged=0;
					DDBYTE fneu=ff;
					
					uint32_t xcoord0=m << SHIFTPERDDBYTE;
					
					for(int32_t bit=0;bit<32;bit++) {
						BYTE tmpf=ff & 0b1;
						ff >>= 1;
						if (tmpf == SQUARE_POTW) continue;
						
						graythere=1;
						
						int xc=xcoord0+bit;
						A.x0=xc*scaleRangePerPixel + COMPLETE0;
						A.x1=A.x0+scaleRangePerPixel;
						
						getBoundingBoxfA(A,bbxfA);
						
						// bbxfA overlaps with outside of cycle enclosement (local)
						if (
							(SQUARE_LIES_ENTIRELY_IN_LOCAL(bbxfA) <= 0) ||
							(SQUARE_LIES_ENTIRELY_IN_COMPLETE(bbxfA) <= 0)
						) {
							fchanged=1;
							fneu |= (1 << bit);
							continue;
						}
						
						ScreenRect scr;
						scr.x0=scrcoord_as_lowerleft(bbxfA.x0);
						scr.x1=scrcoord_as_lowerleft(bbxfA.x1);
						scr.y0=scrcoord_as_lowerleft(bbxfA.y0);
						scr.y1=scrcoord_as_lowerleft(bbxfA.y1);
						// scr is in "screen", i.e. >= 0 and < SCREENWIDTH in all coordinates
						
						// check the intersected with pixels
						int8_t hitspotentiallywhite=0;
						for(int32_t by=scr.y0;by<=scr.y1;by++) {
							for(int32_t bx=scr.x0;bx<=scr.x1;bx++) {
								BYTE bf;
								CELLCOLOR_XY(bx,by,bf);
								if (bf == SQUARE_POTW) {
									hitspotentiallywhite=1;
									break;
								}
							}
							if (hitspotentiallywhite>0) break;
						}
						
						if (hitspotentiallywhite>0) {
							fchanged=1;
							fneu |= (1 << bit);
						}
						
					} // bit
					
					if (fchanged>0) {
						changed=1;
						SET32_MY(m,y,fneu);
					}
				} // m
				
				if (graythere<=0) {
					ywithgray[y-enclosementall.y0]=0;
				}
			} // y
		} // main-while loop as long as new information
		// is being created
		
		// if GRAY cells are present = black emerges
		interiorpresentat=0;
		// ywithgray: not indicative any more of pointer present
		for(int32_t y=0;y<LOCALLENY;y++) {
			if (!ispotwY[y]) continue;
			
			for(int32_t m=0;m<LOCALLENX;m++) {
				// is cell still GRAY
				if (ispotwY[y][m] != ALL32POTW) {
					interiorpresentat=REFINEMENT;
					break;
				} // m
				if (interiorpresentat>0) break;
			}
			if (interiorpresentat>0) break;
		} // k

		if (interiorpresentat>0) break;
	} // REFINEMENT
	
	onecycle.interiorfound=interiorpresentat;
	delete[] ywithgray;
	mgr.FreeAll();

	return interiorpresentat;
}

// struct Root

void Root::clear(void) {
	cyclelen=0;
	cycle=NULL;
	interiorfound=0;
	multiplier=0.0;
}

int32_t main(int32_t argc,char** argv) {
	int32_t c0=clock();
	
	flog=fopen("tsapredictor.log","at");
	fprintf(flog,"\n-----------------\n");
	
	printf("  FUNC=string / c=re,im / A=re,im / ENCW=n / LEVEL=n,m / PERIODS=n,m\n");
	
	// standard
	getBoundingBoxfA=getBoundingBoxfA_z2c;
	_FUNC=FUNC_Z2C;
	COMPLETE0=-2;
	COMPLETE1=2;
	seedC0re=seedC1re=floor(-1.0*DENOM225) / DENOM225; 
	seedC0im=seedC1im=floor(0.0*DENOM225) / DENOM225; 
	FAKTORAre=FAKTORAim=0.0;
	_ENCLOSEMENTWIDTH=128;
	LEVEL0=10;
	LEVEL1=24;
	_STARTWITH=ALL32POTW;
		
	for(int32_t i=1;i<argc;i++) {
		upper(argv[i]);
		if (strstr(argv[i],"FUNC=")==argv[i]) {
			_FUNC=getfuncidx(&argv[i][5]);
		} else if (strstr(argv[i],"C=")==argv[i]) {
			double r0,i0; // not NTYP
			// command line parameters are always considered double no matter the datatype used
			if (sscanf(&argv[i][2],"%lf,%lf",&r0,&i0) == 2) {
				seedC0re=seedC1re=floor(r0*DENOM225)/DENOM225;
				seedC0im=seedC1im=floor(i0*DENOM225)/DENOM225;
			}
		} else if (strstr(argv[i],"A=")==argv[i]) {
			double r0,i0;
			if (sscanf(&argv[i][2],"%lf,%lf",&r0,&i0) == 2) {
				FAKTORAre=floor(r0*DENOM225)/DENOM225;
				FAKTORAim=floor(i0*DENOM225)/DENOM225;
			}
		} else if (strstr(argv[i],"ENCW=")==argv[i]) {
			int32_t a;
			if (sscanf(&argv[i][5],"%i",&a) == 1) {
				if (a < 0) {
					a=-a;
					_STARTWITH=ALL32GRAY; // all gray, i.e. to analyze
				} else {
					_STARTWITH=ALL32POTW; 
				}
				_ENCLOSEMENTWIDTH=a;
			}
		} else if (strstr(argv[i],"LEVEL=")==argv[i]) {
			int32_t a,b;
			if (sscanf(&argv[i][6],"%i,%i",&a,&b) == 2) {
				LEVEL0=a;
				LEVEL1=b;
			}
		} else if (strstr(argv[i],"PERIODS=")==argv[i]) {
			int32_t a,b;
			if (sscanf(&argv[i][8],"%i,%i",&a,&b) == 2) {
				PERIODICLEN0=a;
				PERIODICLEN1=b;
			}
		} 
	} // i
	
	cplxC=Complex(seedC0re,seedC0im);
	cplxA=Complex(FAKTORAre,FAKTORAim);
	
	if (LEVEL0<8) LEVEL0=8;
	if (LEVEL1>31) LEVEL1=31;

	// setting function pointers
	setfunc(_FUNC);
	fkt.output(stdout);
	fkt.output(flog);
	fprintf(flog,"ENCW=%i pixels\n",_ENCLOSEMENTWIDTH);
	if (_STARTWITH == ALL32GRAY) {
		LOGMSG("  per cycle: analyzing whole rectangle around all periodic points\n");
	} else {
		LOGMSG("  per cycle: analyzing small neighbourhoods around periodic point\n");
	}
	
	// must be AFTER setfunc
	// enclosement for filled-in Julia set is computed
	COMPLETE1=getLagrange(fkt);
	COMPLETE0=-COMPLETE1;
	LOGMSG2("Filled-in set is contained in %.0lg-square\n",(double)COMPLETE1);
	LOGMSG2("numerical type: %s\n",NNTYPSTR);
	
	// searching for zeros
	ps_find_critical_points();
	
	if (nbr_of_cp<=0) {
		LOGMSG("No critical points found.\n");
		exit(99);
	}
	
	for(int32_t i=0;i<nbr_of_cp;i++) {
		LOGMSG("critical point: ");
		zero[i].attractor.output(stdout);
		zero[i].attractor.output(flog);
		LOGMSG("\n");
	}
	LOGMSG("\n");
	
	// construct orbits of critical points of bounded
	if (ps_construct_critical_orbits() <= 0) {
		LOGMSG("No critical orbit found.\n(Does an attractor exist at all?)");
		exit(99);
	}
	
	// there may be critical points that go into the
	// same cycle. But only one of them is marked in
	// the array zero
	for(int32_t cp=0;cp<nbr_of_cp;cp++) {
		if (zero[cp].cyclelen>0) {
			// a valid cycle
			LOGMSG4("cycle #%i |multiplier|=%.5lg len=%i: ",
				zero[cp].cyclenumber,
				zero[cp].multiplier,
				zero[cp].cyclelen);
			for(int32_t k=0;k<zero[cp].cyclelen;k++) {
				zero[cp].cycle[k].pp.output(stdout);
				zero[cp].cycle[k].pp.output(flog);
				LOGMSG(" -> ");
			}
			LOGMSG("(reentering ");
			Complex reenter;
			fkt.eval_arg_f(zero[cp].cycle[zero[cp].cyclelen-1].pp,reenter);
			reenter.output(stdout);
			reenter.output(flog);
			LOGMSG(")\n");
		}
	}
	
	// all zeros with cyclelen>0 => analyze per
	// cell mapping
	for(int32_t cp=0;cp<nbr_of_cp;cp++) {
		if (zero[cp].cyclelen<=0) continue;
		
		if (
			(PERIODICLEN0>0)
			&& ( !(
				(PERIODICLEN0 <= zero[cp].cyclelen) &&
				(zero[cp].cyclelen <= PERIODICLEN1)
			) )
		) {
			continue;
		}
		
		LOGMSG3("\nanalyzing cycle #%i (period %i) ...\n",zero[cp].cyclenumber,zero[cp].cyclelen);
			
		int32_t interiorpresent=cm_local(zero[cp],_STARTWITH);
		zero[cp].interiorfound=interiorpresent;
			
		if (interiorpresent>0) {
			LOGMSG2("\n  black present at refinement level %i\n",interiorpresent);
			LOGMSG("  computing this and at latest here emerging cycles from scratch in command-line:\n");
		    LOGMSG5("    juliatsacore_%s range=%.0lg len=%i %s\n",NNTYPSTR,ceil(COMPLETE1),interiorpresent,COMPUTECOMMANDLINE);
		    if (interiorpresent > 12) {
				LOGMSG("  (but level-by-level computation using already calculated data is recommended for speed reasons)\n");
			}
		} else {
			LOGMSG3("\n  NO black found in levels %i..%i at current parameters\n",LEVEL0,LEVEL1);
		}
	} // cp
	
	// do the enclosements of different cycles overlap
	// if so: detected black for a given cycle might
	// have actually detected another (earlier emerging) one
	// only valid if all cycles are actually analyzed (PERIODS command-line)
	int8_t overlapping=0;
	for(int32_t i=0;i<nbr_of_cp;i++) {
		if (
			(zero[i].cyclelen <= 0) ||
			(zero[i].interiorfound <= 0) 
		) continue;
		
		for(int32_t k=0;k<nbr_of_cp;k++) {
			if (
				(i==k) ||
				(zero[k].cyclelen <= 0) ||
				(zero[k].interiorfound <= 0) 
			) continue;
			
			// does enclosement i overlap with k
			if (
				(zero[i].ps_basinrect.x1 < zero[k].ps_basinrect.x0) ||
				(zero[i].ps_basinrect.x0 > zero[k].ps_basinrect.x1) ||
				(zero[i].ps_basinrect.y1 < zero[k].ps_basinrect.y0) ||
				(zero[i].ps_basinrect.y0 > zero[k].ps_basinrect.y1)
			) {
			} else {
				overlapping=1;
				break;
			}
		} // x
		if (overlapping>0) break;
	}
	
	if (overlapping>0) {
		LOGMSG("\n\n!!!!! CAVE !!!!!\n  Enclosements of periodic points of different cycles overlap.\n");
		LOGMSG("  Black when detected for a specific cycle might have actually detected a different one.\n");
	}
	
	int32_t c1=clock();
	LOGMSG2("%.0lf sec duration\n",
		(double)(c1-c0)/CLOCKS_PER_SEC);
		
    return 0;
}

