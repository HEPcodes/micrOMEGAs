#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<unistd.h>

/*======================================
  Higgs, LEP limits and deltarho
  ======================================*/
extern double deltarho_(void);
extern int    masslimits_(void);
extern int    HBblocks(char * fname);
extern int    LiLithF(char*fname);

#define deltarho   deltarho_
#define masslimits masslimits_

/*==========
  UMSSMTools
  ==========*/
extern int    assignValFunc(char * name, double val);
extern int    umssmtools_(void);
extern int    read_prmU(int n);
extern double UparC(int r);
extern double uparctof_(int *r);

extern double bsgnlo_(double *M, double*P);
extern double deltamd_(double *M, double*P);
extern double deltams_(double *M, double*P);
extern double bsmumu_(double *M, double*P);
extern double btaunu_(double *M, double*P);
extern double gmuon_(double *M, double*P);
extern double bxismulow_(double *M, double*P);
extern double bxismuhigh_(double *M, double*P);


#define umssmtools umssmtools_
#define bsgnlo     bsgnlo_
#define deltamd    deltamd_
#define deltams    deltams_
#define bsmumu     bsmumu_
#define btaunu     btaunu_
#define gmuon      gmuon_
#define bxismulow  bxismulow_
#define bxismuhigh bxismuhigh_

/*=====================
  Les Houches interface
  =====================*/
extern int  readLesH(char*fname, int SM );
extern int  lesHinput(char * fname);
extern void FillVal(int mode);

/*=====
  Other
  =====*/
extern void   o1Contents(FILE * f);
extern double randpar(char * xx, double start, double amp, int pwr, int sign, int logmod, int expval);
extern double SignParam (double x);
