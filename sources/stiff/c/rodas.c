#include<stdio.h>

/* rodas.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer mle, mue, mbjac, mbb, mdiag, mdiff, mbdiag;
} linal_;

#define linal_1 linal_

union {
    struct {
	doublereal xold, hout;
	integer nn;
    } _1;
    struct {
	doublereal xold, h__;
	integer n;
    } _2;
} conros_;

#define conros_1 (conros_._1)
#define conros_2 (conros_._2)

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static doublereal c_b61 = 1.;
static logical c_false = FALSE_;
static logical c_true = TRUE_;
static doublereal c_b74 = 0.;
static doublereal c_b78 = .25;

/*<        >*/
/* Subroutine */ int rodas_(integer *n, U_fp fcn, integer *ifcn, doublereal *
	x, doublereal *y, doublereal *xend, doublereal *h__, doublereal *rtol,
	 doublereal *atol, integer *itol, U_fp jac, integer *ijac, integer *
	mljac, integer *mujac, U_fp dfx, integer *idfx, U_fp mas, integer *
	imas, integer *mlmas, integer *mumas, U_fp solout, integer *iout, 
	doublereal *work, integer *lwork, integer *iwork, integer *liwork, 
	doublereal *rpar, integer *ipar, integer *idid)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
//    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
//	    e_wsle(void);

    /* Local variables */
    static integer i__, m1, m2, nm1, iee, lde;
    static doublereal fac1, fac2;
    static integer ndec, njac;
    static doublereal safe;
    static integer ijob, nfcn, ieip;
    static logical pred;
    static integer meth;
    static doublereal hmax;
    static integer iedy, iefx, nmax, nsol, ieak1, ieak2, ieak3, ieak4, ieak5, 
	    ieak6, iedy1, iejac, ldjac;
    static logical jband;
    static integer iecon, iemas, ldmas;
    static logical arret;
    static integer nstep, ldmas2, naccpt, nrejct;
    static logical implct;
    static integer ieynew, istore;
    static logical autnms;
    extern /* Subroutine */ int roscor_(integer *, U_fp, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, U_fp, integer *, integer *,
	     integer *, U_fp, integer *, U_fp, integer *, integer *, U_fp, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, logical *, logical *,
	     logical *, logical *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *);
    static doublereal uround;

    /* Fortran I/O blocks */
//    static cilist io___10 = { 0, 6, 0, 0, 0 };
//    static cilist io___12 = { 0, 6, 0, 0, 0 };
//    static cilist io___17 = { 0, 6, 0, 0, 0 };
//    static cilist io___19 = { 0, 6, 0, 0, 0 };
//    static cilist io___23 = { 0, 6, 0, 0, 0 };
//    static cilist io___25 = { 0, 6, 0, 0, 0 };
//    static cilist io___26 = { 0, 6, 0, 0, 0 };
//    static cilist io___28 = { 0, 6, 0, 0, 0 };
//    static cilist io___36 = { 0, 6, 0, 0, 0 };
//    static cilist io___53 = { 0, 6, 0, 0, 0 };
//    static cilist io___55 = { 0, 6, 0, 0, 0 };


/* ---------------------------------------------------------- */
/*     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC) */
/*     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS  MY'=F(X,Y). */
/*     THIS IS AN EMBEDDED ROSENBROCK METHOD OF ORDER (3)4 */
/*     (WITH STEP SIZE CONTROL). */
/*     C.F. SECTIONS IV.7  AND VI.3 */

/*     AUTHORS: E. HAIRER AND G. WANNER */
/*              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES */
/*              CH-1211 GENEVE 24, SWITZERLAND */
/*              E-MAIL:  Ernst.Hairer@math.unige.ch */
/*                       Gerhard.Wanner@math.unige.ch */

/*     THIS CODE IS PART OF THE BOOK: */
/*         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL */
/*         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS. */
/*         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14, */
/*         SPRINGER-VERLAG 1991, SECOND EDITION 1996. */

/*     VERSION OF OCTOBER 28, 1996 */

/*     INPUT PARAMETERS */
/*     ---------------- */
/*     N           DIMENSION OF THE SYSTEM */

/*     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE */
/*                 VALUE OF F(X,Y): */
/*                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR) */
/*                    DOUBLE PRECISION X,Y(N),F(N) */
/*                    F(1)=...   ETC. */
/*                 RPAR, IPAR (SEE BELOW) */

/*     IFCN        GIVES INFORMATION ON FCN: */
/*                    IFCN=0: F(X,Y) INDEPENDENT OF X (AUTONOMOUS) */
/*                    IFCN=1: F(X,Y) MAY DEPEND ON X (NON-AUTONOMOUS) */

/*     X           INITIAL X-VALUE */

/*     Y(N)        INITIAL VALUES FOR Y */

/*     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE) */

/*     H           INITIAL STEP SIZE GUESS; */
/*                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT, */
/*                 H=1.D0/(NORM OF F'), USUALLY 1.D-2 OR 1.D-3, IS GOOD. */
/*                 THIS CHOICE IS NOT VERY IMPORTANT, THE CODE QUICKLY */
/*                 ADAPTS ITS STEP SIZE (IF H=0.D0, THE CODE PUTS H=1.D-6). */

/*     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY */
/*                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N. */

/*     ITOL        SWITCH FOR RTOL AND ATOL: */
/*                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS. */
/*                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF */
/*                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL */
/*                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS. */
/*                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW */
/*                     RTOL(I)*ABS(Y(I))+ATOL(I). */

/*     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES */
/*                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y */
/*                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY */
/*                 A DUMMY SUBROUTINE IN THE CASE IJAC=0). */
/*                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM */
/*                    SUBROUTINE JAC(N,X,Y,DFY,LDFY,RPAR,IPAR) */
/*                    DOUBLE PRECISION X,Y(N),DFY(LDFY,N) */
/*                    DFY(1,1)= ... */
/*                 LDFY, THE COLOMN-LENGTH OF THE ARRAY, IS */
/*                 FURNISHED BY THE CALLING PROGRAM. */
/*                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO */
/*                    BE FULL AND THE PARTIAL DERIVATIVES ARE */
/*                    STORED IN DFY AS */
/*                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J) */
/*                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND */
/*                    THE PARTIAL DERIVATIVES ARE STORED */
/*                    DIAGONAL-WISE AS */
/*                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J). */

/*     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN: */
/*                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE */
/*                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED. */
/*                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC. */

/*     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN: */
/*                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR */
/*                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION. */
/*                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN */
/*                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW */
/*                       THE MAIN DIAGONAL). */

/*     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON- */
/*                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL). */
/*                 NEED NOT BE DEFINED IF MLJAC=N. */

/*     DFX         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES */
/*                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO X */
/*                 (THIS ROUTINE IS ONLY CALLED IF IDFX=1 AND IFCN=1; */
/*                 SUPPLY A DUMMY SUBROUTINE IN THE CASE IDFX=0 OR IFCN=0). */
/*                 OTHERWISE, THIS SUBROUTINE MUST HAVE THE FORM */
/*                    SUBROUTINE DFX(N,X,Y,FX,RPAR,IPAR) */
/*                    DOUBLE PRECISION X,Y(N),FX(N) */
/*                    FX(1)= ... */

/*     IDFX        SWITCH FOR THE COMPUTATION OF THE DF/DX: */
/*                    IDFX=0: DF/DX IS COMPUTED INTERNALLY BY FINITE */
/*                       DIFFERENCES, SUBROUTINE "DFX" IS NEVER CALLED. */
/*                    IDFX=1: DF/DX IS SUPPLIED BY SUBROUTINE DFX. */

/*     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      ----- */
/*     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): - */

/*     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS- */
/*                 MATRIX M. */
/*                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY */
/*                 MATRIX AND NEEDS NOT TO BE DEFINED; */
/*                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE. */
/*                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM */
/*                    SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR) */
/*                    DOUBLE PRECISION AM(LMAS,N) */
/*                    AM(1,1)= .... */
/*                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED */
/*                    AS FULL MATRIX LIKE */
/*                         AM(I,J) = M(I,J) */
/*                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED */
/*                    DIAGONAL-WISE AS */
/*                         AM(I-J+MUMAS+1,J) = M(I,J). */

/*     IMAS       GIVES INFORMATION ON THE MASS-MATRIX: */
/*                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY */
/*                       MATRIX, MAS IS NEVER CALLED. */
/*                    IMAS=1: MASS-MATRIX  IS SUPPLIED. */

/*     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX: */
/*                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR */
/*                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION. */
/*                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE */
/*                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW */
/*                       THE MAIN DIAGONAL). */
/*                 MLMAS IS SUPPOSED TO BE .LE. MLJAC. */

/*     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON- */
/*                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL). */
/*                 NEED NOT BE DEFINED IF MLMAS=N. */
/*                 MUMAS IS SUPPOSED TO BE .LE. MUJAC. */

/*     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE */
/*                 NUMERICAL SOLUTION DURING INTEGRATION. */
/*                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP. */
/*                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. */
/*                 IT MUST HAVE THE FORM */
/*                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N, */
/*                                       RPAR,IPAR,IRTRN) */
/*                    DOUBLE PRECISION X,Y(N),CONT(LRC) */
/*                    .... */
/*                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH */
/*                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS */
/*                    THE FIRST GRID-POINT). */
/*                 "XOLD" IS THE PRECEEDING GRID-POINT. */
/*                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN */
/*                    IS SET <0, RODAS RETURNS TO THE CALLING PROGRAM. */

/*          -----  CONTINUOUS OUTPUT: ----- */
/*                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION */
/*                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH */
/*                 THE FUNCTION */
/*                        >>>   CONTRO(I,S,CONT,LRC)   <<< */
/*                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH */
/*                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE */
/*                 S SHOULD LIE IN THE INTERVAL [XOLD,X]. */

/*     IOUT        GIVES INFORMATION ON THE SUBROUTINE SOLOUT: */
/*                    IOUT=0: SUBROUTINE IS NEVER CALLED */
/*                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT */

/*     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK". */
/*                 SERVES AS WORKING SPACE FOR ALL VECTORS AND MATRICES. */
/*                 "LWORK" MUST BE AT LEAST */
/*                             N*(LJAC+LMAS+LE1+14)+20 */
/*                 WHERE */
/*                    LJAC=N              IF MLJAC=N (FULL JACOBIAN) */
/*                    LJAC=MLJAC+MUJAC+1  IF MLJAC<N (BANDED JAC.) */
/*                 AND */
/*                    LMAS=0              IF IMAS=0 */
/*                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL) */
/*                    LMAS=MLMAS+MUMAS+1  IF MLMAS<N (BANDED MASS-M.) */
/*                 AND */
/*                    LE1=N               IF MLJAC=N (FULL JACOBIAN) */
/*                    LE1=2*MLJAC+MUJAC+1 IF MLJAC<N (BANDED JAC.). */
/*                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE */
/*                 MASS-MATRIX IS THE INDENTITY (IMAS=0), THE MINIMUM */
/*                 STORAGE REQUIREMENT IS */
/*                             LWORK = 2*N*N+14*N+20. */
/*                 IF IWORK(9)=M1>0 THEN "LWORK" MUST BE AT LEAST */
/*                          N*(LJAC+14)+(N-M1)*(LMAS+LE1)+20 */
/*                 WHERE IN THE DEFINITIONS OF LJAC, LMAS AND LE1 THE */
/*                 NUMBER N CAN BE REPLACED BY N-M1. */

/*     LWORK       DECLARED LENGTH OF ARRAY "WORK". */

/*     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK". */
/*                 "LIWORK" MUST BE AT LEAST N+20. */

/*     LIWORK      DECLARED LENGTH OF ARRAY "IWORK". */

/*     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH */
/*                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING */
/*                 PROGRAM AND THE FCN, DFX, JAC, MAS, SOLOUT SUBROUTINES. */

/* ---------------------------------------------------------------------- */

/*     SOPHISTICATED SETTING OF PARAMETERS */
/*     ----------------------------------- */
/*              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK */
/*              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),..,WORK(4) */
/*              AS WELL AS IWORK(1),IWORK(2) DIFFERENT FROM ZERO. */
/*              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES: */

/*    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS. */
/*              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 100000. */

/*    IWORK(2)  SWITCH FOR THE CHOICE OF THE COEFFICIENTS */
/*              IF IWORK(2).EQ.1  METHOD (SEE BOOK, PAGE 452) */
/*              IF IWORK(2).EQ.2  SAME METHOD WITH DIFFERENT PARAMETERS */
/*              IF IWORK(2).EQ.3  METHOD WITH COEFF. OF GERD STEINEBACH */
/*              THE DEFAULT VALUE (FOR IWORK(2)=0) IS IWORK(2)=1. */

/*    IWORK(3)  SWITCH FOR STEP SIZE STRATEGY */
/*              IF IWORK(3).EQ.1  MOD. PREDICTIVE CONTROLLER (GUSTAFSSON) */
/*              IF IWORK(3).EQ.2  CLASSICAL APPROACH */
/*              THE DEFAULT VALUE (FOR IWORK(3)=0) IS IWORK(3)=1. */

/*       IF THE DIFFERENTIAL SYSTEM HAS THE SPECIAL STRUCTURE THAT */
/*            Y(I)' = Y(I+M2)   FOR  I=1,...,M1, */
/*       WITH M1 A MULTIPLE OF M2, A SUBSTANTIAL GAIN IN COMPUTERTIME */
/*       CAN BE ACHIEVED BY SETTING THE PARAMETERS IWORK(9) AND IWORK(10). */
/*       E.G., FOR SECOND ORDER SYSTEMS P'=V, V'=G(P,V), WHERE P AND V ARE */
/*       VECTORS OF DIMENSION N/2, ONE HAS TO PUT M1=M2=N/2. */
/*       FOR M1>0 SOME OF THE INPUT PARAMETERS HAVE DIFFERENT MEANINGS: */
/*       - JAC: ONLY THE ELEMENTS OF THE NON-TRIVIAL PART OF THE */
/*              JACOBIAN HAVE TO BE STORED */
/*              IF (MLJAC.EQ.N-M1) THE JACOBIAN IS SUPPOSED TO BE FULL */
/*                 DFY(I,J) = PARTIAL F(I+M1) / PARTIAL Y(J) */
/*                FOR I=1,N-M1 AND J=1,N. */
/*              ELSE, THE JACOBIAN IS BANDED ( M1 = M2 * MM ) */
/*                 DFY(I-J+MUJAC+1,J+K*M2) = PARTIAL F(I+M1) / PARTIAL Y(J+K*M2) */
/*                FOR I=1,MLJAC+MUJAC+1 AND J=1,M2 AND K=0,MM. */
/*       - MLJAC: MLJAC=N-M1: IF THE NON-TRIVIAL PART OF THE JACOBIAN IS FULL */
/*                0<=MLJAC<N-M1: IF THE (MM+1) SUBMATRICES (FOR K=0,MM) */
/*                     PARTIAL F(I+M1) / PARTIAL Y(J+K*M2),  I,J=1,M2 */
/*                    ARE BANDED, MLJAC IS THE MAXIMAL LOWER BANDWIDTH */
/*                    OF THESE MM+1 SUBMATRICES */
/*       - MUJAC: MAXIMAL UPPER BANDWIDTH OF THESE MM+1 SUBMATRICES */
/*                NEED NOT BE DEFINED IF MLJAC=N-M1 */
/*       - MAS: IF IMAS=0 THIS MATRIX IS ASSUMED TO BE THE IDENTITY AND */
/*              NEED NOT BE DEFINED. SUPPLY A DUMMY SUBROUTINE IN THIS CASE. */
/*              IT IS ASSUMED THAT ONLY THE ELEMENTS OF RIGHT LOWER BLOCK OF */
/*              DIMENSION N-M1 DIFFER FROM THAT OF THE IDENTITY MATRIX. */
/*              IF (MLMAS.EQ.N-M1) THIS SUBMATRIX IS SUPPOSED TO BE FULL */
/*                 AM(I,J) = M(I+M1,J+M1)     FOR I=1,N-M1 AND J=1,N-M1. */
/*              ELSE, THE MASS MATRIX IS BANDED */
/*                 AM(I-J+MUMAS+1,J) = M(I+M1,J+M1) */
/*       - MLMAS: MLMAS=N-M1: IF THE NON-TRIVIAL PART OF M IS FULL */
/*                0<=MLMAS<N-M1: LOWER BANDWIDTH OF THE MASS MATRIX */
/*       - MUMAS: UPPER BANDWIDTH OF THE MASS MATRIX */
/*                NEED NOT BE DEFINED IF MLMAS=N-M1 */

/*    IWORK(9)  THE VALUE OF M1.  DEFAULT M1=0. */

/*    IWORK(10) THE VALUE OF M2.  DEFAULT M2=M1. */

/*    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16. */

/*    WORK(2)   MAXIMAL STEP SIZE, DEFAULT XEND-X. */

/*    WORK(3), WORK(4)   PARAMETERS FOR STEP SIZE SELECTION */
/*              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION */
/*                 WORK(3) <= HNEW/HOLD <= WORK(4) */
/*              DEFAULT VALUES: WORK(3)=0.2D0, WORK(4)=6.D0 */

/*    WORK(5)   THE SAFETY FACTOR IN STEP SIZE PREDICTION, */
/*              DEFAULT 0.9D0. */

/* ----------------------------------------------------------------------- */

/*     OUTPUT PARAMETERS */
/*     ----------------- */
/*     X           X-VALUE WHERE THE SOLUTION IS COMPUTED */
/*                 (AFTER SUCCESSFUL RETURN X=XEND) */

/*     Y(N)        SOLUTION AT X */

/*     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP */

/*     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN: */
/*                   IDID= 1  COMPUTATION SUCCESSFUL, */
/*                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT) */
/*                   IDID=-1  INPUT IS NOT CONSISTENT, */
/*                   IDID=-2  LARGER NMAX IS NEEDED, */
/*                   IDID=-3  STEP SIZE BECOMES TOO SMALL, */
/*                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR. */

/*   IWORK(14)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL */
/*                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED) */
/*   IWORK(15)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY */
/*                      OR NUMERICALLY) */
/*   IWORK(16)  NSTEP   NUMBER OF COMPUTED STEPS */
/*   IWORK(17)  NACCPT  NUMBER OF ACCEPTED STEPS */
/*   IWORK(18)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST), */
/*                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED) */
/*   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS (N-DIMENSIONAL MATRIX) */
/*   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS */
/* --------------------------------------------------------- */
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/*          DECLARATIONS */
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       DIMENSION Y(N),ATOL(*),RTOL(*),WORK(LWORK),IWORK(LIWORK) >*/
/*<       DIMENSION RPAR(*),IPAR(*) >*/
/*<       LOGICAL AUTNMS,IMPLCT,JBAND,ARRET,PRED >*/
/*<       EXTERNAL FCN,JAC,DFX,MAS,SOLOUT >*/
/* *** *** *** *** *** *** *** */
/*        SETTING THE PARAMETERS */
/* *** *** *** *** *** *** *** */
/*<       NFCN=0 >*/
    /* Parameter adjustments */
    --y;
    --rtol;
    --atol;
    --work;
    --iwork;
    --rpar;
    --ipar;

    /* Function Body */
    nfcn = 0;
/*<       NACCPT=0 >*/
    naccpt = 0;
/*<       NREJCT=0 >*/
    nrejct = 0;
/*<       NSTEP=0 >*/
    nstep = 0;
/*<       NJAC=0 >*/
    njac = 0;
/*<       NDEC=0 >*/
    ndec = 0;
/*<       NSOL=0 >*/
    nsol = 0;
/*<       ARRET=.FALSE. >*/
    arret = FALSE_;
/* -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- */
/*<       IF(IWORK(1).EQ.0)THEN >*/
    if (iwork[1] == 0) {
/*<          NMAX=10000000 >*/
	nmax = 10000000;
/*<       ELSE >*/
    } else {
/*<          NMAX=IWORK(1) >*/
	nmax = iwork[1];
/*<          IF(NMAX.LE.0)THEN >*/
	if (nmax <= 0) {
/*<             WRITE(6,*)' WRONG INPUT IWORK(1)=',IWORK(1) >*/
//	    s_wsle(&io___10);
//	    do_lio(&c__9, &c__1, " WRONG INPUT IWORK(1)=", (ftnlen)22);
//	    do_lio(&c__3, &c__1, (char *)&iwork[1], (ftnlen)sizeof(integer));
//	    e_wsle();
            printf(" WRONG INPUT IWORK(1)=%d\n", iwork[1]);
/*<             ARRET=.TRUE. >*/
	    arret = TRUE_;
/*<          END IF >*/
	}
/*<       END IF >*/
    }
/* -------- METH   COEFFICIENTS OF THE METHOD */
/*<       IF(IWORK(2).EQ.0)THEN >*/
    if (iwork[2] == 0) {
/*<          METH=1 >*/
	meth = 1;
/*<       ELSE >*/
    } else {
/*<          METH=IWORK(2) >*/
	meth = iwork[2];
/*<          IF(METH.LE.0.OR.METH.GE.4)THEN >*/
	if (meth <= 0 || meth >= 4) {
/*<             WRITE(6,*)' CURIOUS INPUT IWORK(2)=',IWORK(2) >*/
//	    s_wsle(&io___12);
//	    do_lio(&c__9, &c__1, " CURIOUS INPUT IWORK(2)=", (ftnlen)24);
//	    do_lio(&c__3, &c__1, (char *)&iwork[2], (ftnlen)sizeof(integer));
//	    e_wsle();
           printf("CURIOUS INPUT IWORK(2)=%d\n", iwork[2]);
/*<             ARRET=.TRUE. >*/
	    arret = TRUE_;
/*<          END IF >*/
	}
/*<       END IF >*/
    }
/* -------- PRED   STEP SIZE CONTROL */
/*<       IF(IWORK(3).LE.1)THEN >*/
    if (iwork[3] <= 1) {
/*<          PRED=.TRUE. >*/
	pred = TRUE_;
/*<       ELSE >*/
    } else {
/*<          PRED=.FALSE. >*/
	pred = FALSE_;
/*<       END IF >*/
    }
/* -------- PARAMETER FOR SECOND ORDER EQUATIONS */
/*<       M1=IWORK(9) >*/
    m1 = iwork[9];
/*<       M2=IWORK(10) >*/
    m2 = iwork[10];
/*<       NM1=N-M1 >*/
    nm1 = *n - m1;
/*<       IF (M1.EQ.0) M2=N >*/
    if (m1 == 0) {
	m2 = *n;
    }
/*<       IF (M2.EQ.0) M2=M1 >*/
    if (m2 == 0) {
	m2 = m1;
    }
/*<       IF (M1.LT.0.OR.M2.LT.0.OR.M1+M2.GT.N) THEN >*/
    if (m1 < 0 || m2 < 0 || m1 + m2 > *n) {
/*<        WRITE(6,*)' CURIOUS INPUT FOR IWORK(9,10)=',M1,M2 >*/
//	s_wsle(&io___17);
//	do_lio(&c__9, &c__1, " CURIOUS INPUT FOR IWORK(9,10)=", (ftnlen)31);
//	do_lio(&c__3, &c__1, (char *)&m1, (ftnlen)sizeof(integer));
//	do_lio(&c__3, &c__1, (char *)&m2, (ftnlen)sizeof(integer));
//	e_wsle();
        printf("CURIOUS INPUT FOR IWORK(9,10)=%d,%d\n",m1,m2);
/*<        ARRET=.TRUE. >*/
	arret = TRUE_;
/*<       END IF >*/
    }
/* -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0 */
/*<       IF(WORK(1).EQ.0.D0)THEN >*/
    if (work[1] == 0.) {
/*<          UROUND=1.D-16 >*/
	uround = 1e-16;
/*<       ELSE >*/
    } else {
/*<          UROUND=WORK(1) >*/
	uround = work[1];
/*<          IF(UROUND.LT.1.D-16.OR.UROUND.GE.1.D0)THEN >*/
	if (uround < 1e-16 || uround >= 1.) {
/*<             WRITE(6,*)' COEFFICIENTS HAVE 16 DIGITS, UROUND=',WORK(1) >*/
//	    s_wsle(&io___19);
//	    do_lio(&c__9, &c__1, " COEFFICIENTS HAVE 16 DIGITS, UROUND=", (
//		    ftnlen)37);
//	    do_lio(&c__5, &c__1, (char *)&work[1], (ftnlen)sizeof(doublereal))
		    ;
//	    e_wsle();
        printf("COEFFICIENTS HAVE 16 DIGITS, UROUND=%d\n", work[1]);
/*<             ARRET=.TRUE. >*/
	    arret = TRUE_;
/*<          END IF >*/
	}
/*<       END IF >*/
    }
/* -------- MAXIMAL STEP SIZE */
/*<       IF(WORK(2).EQ.0.D0)THEN >*/
    if (work[2] == 0.) {
/*<          HMAX=XEND-X >*/
	hmax = *xend - *x;
/*<       ELSE >*/
    } else {
/*<          HMAX=WORK(2) >*/
	hmax = work[2];
/*<       END IF >*/
    }
/* -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION */
/*<       IF(WORK(3).EQ.0.D0)THEN >*/
    if (work[3] == 0.) {
/*<          FAC1=5.D0 >*/
	fac1 = 5.;
/*<       ELSE >*/
    } else {
/*<          FAC1=1.D0/WORK(3) >*/
	fac1 = 1. / work[3];
/*<       END IF >*/
    }
/*<       IF(WORK(4).EQ.0.D0)THEN >*/
    if (work[4] == 0.) {
/*<          FAC2=1.D0/6.0D0 >*/
	fac2 = .16666666666666666;
/*<       ELSE >*/
    } else {
/*<          FAC2=1.D0/WORK(4) >*/
	fac2 = 1. / work[4];
/*<       END IF >*/
    }
/*<       IF (FAC1.LT.1.0D0.OR.FAC2.GT.1.0D0) THEN >*/
    if (fac1 < 1. || fac2 > 1.) {
/*<             WRITE(6,*)' CURIOUS INPUT WORK(3,4)=',WORK(3),WORK(4) >*/
//	s_wsle(&io___23);
//	do_lio(&c__9, &c__1, " CURIOUS INPUT WORK(3,4)=", (ftnlen)25);
//	do_lio(&c__5, &c__1, (char *)&work[3], (ftnlen)sizeof(doublereal));
//	do_lio(&c__5, &c__1, (char *)&work[4], (ftnlen)sizeof(doublereal));
//	e_wsle();
        printf("' CURIOUS INPUT WORK(3,4)= %d %d\n", work[3],work[4]);
/*<             ARRET=.TRUE. >*/
	arret = TRUE_;
/*<          END IF >*/
    }
/* --------- SAFE     SAFETY FACTOR IN STEP SIZE PREDICTION */
/*<       IF (WORK(5).EQ.0.0D0) THEN >*/
    if (work[5] == 0.) {
/*<          SAFE=0.9D0 >*/
	safe = .9;
/*<       ELSE >*/
    } else {
/*<          SAFE=WORK(5) >*/
	safe = work[5];
/*<          IF (SAFE.LE.0.001D0.OR.SAFE.GE.1.0D0) THEN >*/
	if (safe <= .001 || safe >= 1.) {
/*<             WRITE(6,*)' CURIOUS INPUT FOR WORK(5)=',WORK(5) >*/
//	    s_wsle(&io___25);
//	    do_lio(&c__9, &c__1, " CURIOUS INPUT FOR WORK(5)=", (ftnlen)27);
//	    do_lio(&c__5, &c__1, (char *)&work[5], (ftnlen)sizeof(doublereal))
//		    ;
//	    e_wsle();
           printf("CURIOUS INPUT FOR WORK(5)=%d\n", work[5]);
/*<             ARRET=.TRUE. >*/
	    arret = TRUE_;
/*<          END IF >*/
	}
/*<       END IF >*/
    }
/* --------- CHECK IF TOLERANCES ARE O.K. */
/*<       IF (ITOL.EQ.0) THEN >*/
    if (*itol == 0) {
/*<           IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE.10.D0*UROUND) THEN >*/
	if (atol[1] <= 0. || rtol[1] <= uround * 10.) {
/*<               WRITE (6,*) ' TOLERANCES ARE TOO SMALL' >*/
//	    s_wsle(&io___26);
//	    do_lio(&c__9, &c__1, " TOLERANCES ARE TOO SMALL", (ftnlen)25);
//	    e_wsle();
            printf(" TOLERANCES ARE TOO SMALL\n");
/*<               ARRET=.TRUE. >*/
	    arret = TRUE_;
/*<           END IF >*/
	}
/*<       ELSE >*/
    } else {
/*<           DO I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<              IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE.10.D0*UROUND) THEN >*/
	    if (atol[i__] <= 0. || rtol[i__] <= uround * 10.) {
/*<                 WRITE (6,*) ' TOLERANCES(',I,') ARE TOO SMALL' >*/
//		s_wsle(&io___28);
//		do_lio(&c__9, &c__1, " TOLERANCES(", (ftnlen)12);
//		do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
//		do_lio(&c__9, &c__1, ") ARE TOO SMALL", (ftnlen)15);
//		e_wsle();
           printf("TOLERANCES(%d) ARE TOO SMALL\n",i__);
/*<                 ARRET=.TRUE. >*/
		arret = TRUE_;
/*<              END IF >*/
	    }
/*<           END DO >*/
	}
/*<       END IF >*/
    }
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/*         COMPUTATION OF ARRAY ENTRIES */
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/* ---- AUTONOMOUS, IMPLICIT, BANDED OR NOT ? */
/*<       AUTNMS=IFCN.EQ.0 >*/
    autnms = *ifcn == 0;
/*<       IMPLCT=IMAS.NE.0 >*/
    implct = *imas != 0;
/*<       JBAND=MLJAC.LT.NM1 >*/
    jband = *mljac < nm1;
/* -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS --- */
/* -- JACOBIAN AND MATRIX E */
/*<       IF(JBAND)THEN >*/
    if (jband) {
/*<          LDJAC=MLJAC+MUJAC+1 >*/
	ldjac = *mljac + *mujac + 1;
/*<          LDE=MLJAC+LDJAC >*/
	lde = *mljac + ldjac;
/*<       ELSE >*/
    } else {
/*<          MLJAC=NM1 >*/
	*mljac = nm1;
/*<          MUJAC=NM1 >*/
	*mujac = nm1;
/*<          LDJAC=NM1 >*/
	ldjac = nm1;
/*<          LDE=NM1 >*/
	lde = nm1;
/*<       END IF >*/
    }
/* -- MASS MATRIX */
/*<       IF (IMPLCT) THEN >*/
    if (implct) {
/*<           IF (MLMAS.NE.NM1) THEN >*/
	if (*mlmas != nm1) {
/*<               LDMAS=MLMAS+MUMAS+1 >*/
	    ldmas = *mlmas + *mumas + 1;
/*<               IF (JBAND) THEN >*/
	    if (jband) {
/*<                  IJOB=4 >*/
		ijob = 4;
/*<               ELSE >*/
	    } else {
/*<                  IJOB=3 >*/
		ijob = 3;
/*<               END IF >*/
	    }
/*<           ELSE >*/
	} else {
/*<               LDMAS=NM1 >*/
	    ldmas = nm1;
/*<               IJOB=5 >*/
	    ijob = 5;
/*<           END IF >*/
	}
/* ------ BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF "JAC" */
/*<           IF (MLMAS.GT.MLJAC.OR.MUMAS.GT.MUJAC) THEN >*/
	if (*mlmas > *mljac || *mumas > *mujac) {
/*<        >*/
//	    s_wsle(&io___36);
//	    do_lio(&c__9, &c__1, "BANDWITH OF \"MAS\" NOT LARGER THAN BANDWI"
//		    "TH OF \"JAC\"", (ftnlen)51);
//	    e_wsle();
           printf("BANDWITH OF 'MAS' NOT LARGER THAN BANDWITH OF 'JAC'\n"); 
/*<             ARRET=.TRUE. >*/
	    arret = TRUE_;
/*<           END IF >*/
	}
/*<       ELSE >*/
    } else {
/*<           LDMAS=0 >*/
	ldmas = 0;
/*<           IF (JBAND) THEN >*/
	if (jband) {
/*<              IJOB=2 >*/
	    ijob = 2;
/*<           ELSE >*/
	} else {
/*<              IJOB=1 >*/
	    ijob = 1;
/*<           END IF >*/
	}
/*<       END IF >*/
    }
/*<       LDMAS2=MAX(1,LDMAS) >*/
    ldmas2 = max(1,ldmas);
/* ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK ----- */
/*<       IEYNEW=21 >*/
    ieynew = 21;
/*<       IEDY1=IEYNEW+N >*/
    iedy1 = ieynew + *n;
/*<       IEDY=IEDY1+N >*/
    iedy = iedy1 + *n;
/*<       IEAK1=IEDY+N >*/
    ieak1 = iedy + *n;
/*<       IEAK2=IEAK1+N >*/
    ieak2 = ieak1 + *n;
/*<       IEAK3=IEAK2+N >*/
    ieak3 = ieak2 + *n;
/*<       IEAK4=IEAK3+N >*/
    ieak4 = ieak3 + *n;
/*<       IEAK5=IEAK4+N >*/
    ieak5 = ieak4 + *n;
/*<       IEAK6=IEAK5+N >*/
    ieak6 = ieak5 + *n;
/*<       IEFX =IEAK6+N >*/
    iefx = ieak6 + *n;
/*<       IECON=IEFX+N >*/
    iecon = iefx + *n;
/*<       IEJAC=IECON+4*N >*/
    iejac = iecon + (*n << 2);
/*<       IEMAS=IEJAC+N*LDJAC >*/
    iemas = iejac + *n * ldjac;
/*<       IEE  =IEMAS+NM1*LDMAS >*/
    iee = iemas + nm1 * ldmas;
/* ------ TOTAL STORAGE REQUIREMENT ----------- */
/*<       ISTORE=IEE+NM1*LDE-1 >*/
    istore = iee + nm1 * lde - 1;
/*<       IF(ISTORE.GT.LWORK)THEN >*/
    if (istore > *lwork) {
/*<          WRITE(6,*)' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',ISTORE >*/
//	s_wsle(&io___53);
//	do_lio(&c__9, &c__1, " INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=", (
//		ftnlen)43);
//	do_lio(&c__3, &c__1, (char *)&istore, (ftnlen)sizeof(integer));
//	e_wsle();
         printf("INSUFFICIENT STORAGE FOR WORK, MIN. LWORK= %d\n", istore);
/*<          ARRET=.TRUE. >*/
	arret = TRUE_;
/*<       END IF >*/
    }
/* ------- ENTRY POINTS FOR INTEGER WORKSPACE ----- */
/*<       IEIP=21 >*/
    ieip = 21;
/*<       ISTORE=IEIP+NM1-1 >*/
    istore = ieip + nm1 - 1;
/*<       IF(ISTORE.GT.LIWORK)THEN >*/
    if (istore > *liwork) {
/*<          WRITE(6,*)' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=',ISTORE >*/
//	s_wsle(&io___55);
//	do_lio(&c__9, &c__1, " INSUFF. STORAGE FOR IWORK, MIN. LIWORK=", (
//		ftnlen)40);
//	do_lio(&c__3, &c__1, (char *)&istore, (ftnlen)sizeof(integer));
//	e_wsle();
          printf("' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=%d\n",istore);
/*<          ARRET=.TRUE. >*/
	arret = TRUE_;
/*<       END IF >*/
    }
/* ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1 */
/*<       IF (ARRET) THEN >*/
    if (arret) {
/*<          IDID=-1 >*/
	*idid = -1;
/*<          RETURN >*/
	return 0;
/*<       END IF >*/
    }
/* -------- CALL TO CORE INTEGRATOR ------------ */
/*<        >*/
    roscor_(n, (U_fp)fcn, x, &y[1], xend, &hmax, h__, &rtol[1], &atol[1], 
	    itol, (U_fp)jac, ijac, mljac, mujac, (U_fp)dfx, idfx, (U_fp)mas, 
	    mlmas, mumas, (U_fp)solout, iout, idid, &nmax, &uround, &meth, &
	    ijob, &fac1, &fac2, &safe, &autnms, &implct, &jband, &pred, &
	    ldjac, &lde, &ldmas2, &work[ieynew], &work[iedy1], &work[iedy], &
	    work[ieak1], &work[ieak2], &work[ieak3], &work[ieak4], &work[
	    ieak5], &work[ieak6], &work[iefx], &work[iejac], &work[iee], &
	    work[iemas], &iwork[ieip], &work[iecon], &m1, &m2, &nm1, &nfcn, &
	    njac, &nstep, &naccpt, &nrejct, &ndec, &nsol, &rpar[1], &ipar[1]);
/*<       IWORK(14)=NFCN >*/
    iwork[14] = nfcn;
/*<       IWORK(15)=NJAC >*/
    iwork[15] = njac;
/*<       IWORK(16)=NSTEP >*/
    iwork[16] = nstep;
/*<       IWORK(17)=NACCPT >*/
    iwork[17] = naccpt;
/*<       IWORK(18)=NREJCT >*/
    iwork[18] = nrejct;
/*<       IWORK(19)=NDEC >*/
    iwork[19] = ndec;
/*<       IWORK(20)=NSOL >*/
    iwork[20] = nsol;
/* ----------- RETURN ----------- */
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* rodas_ */




/*  ----- ... AND HERE IS THE CORE INTEGRATOR  ---------- */

/*<        >*/
/* Subroutine */ int roscor_(integer *n, S_fp fcn, doublereal *x, doublereal *
	y, doublereal *xend, doublereal *hmax, doublereal *h__, doublereal *
	rtol, doublereal *atol, integer *itol, S_fp jac, integer *ijac, 
	integer *mljac, integer *mujac, S_fp dfx, integer *idfx, S_fp mas, 
	integer *mlmas, integer *mumas, S_fp solout, integer *iout, integer *
	idid, integer *nmax, doublereal *uround, integer *meth, integer *ijob,
	 doublereal *fac1, doublereal *fac2, doublereal *safe, logical *
	autnms, logical *implct, logical *banded, logical *pred, integer *
	ldjac, integer *lde, integer *ldmas, doublereal *ynew, doublereal *
	dy1, doublereal *dy, doublereal *ak1, doublereal *ak2, doublereal *
	ak3, doublereal *ak4, doublereal *ak5, doublereal *ak6, doublereal *
	fx, doublereal *fjac, doublereal *e, doublereal *fmas, integer *ip, 
	doublereal *cont, integer *m1, integer *m2, integer *nm1, integer *
	nfcn, integer *njac, integer *nstep, integer *naccpt, integer *nrejct,
	 integer *ndec, integer *nsol, doublereal *rpar, integer *ipar)
{
    /* Format strings */
    static char fmt_979[] = "(\002 EXIT OF RODAS AT X=\002,e18.4)";

    /* System generated locals */
    integer fjac_dim1, fjac_offset, e_dim1, e_offset, fmas_dim1, fmas_offset, 
	    i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal), pow_dd(
	    doublereal *, doublereal *);
//    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
//	     s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
//	    e_wsle(void);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal c2, c3, c4, d1, d2, d3, d4;
    static integer j1;
    static doublereal a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, c21, 
	    c31, c32, c41, c42, c43, c51, c52, c53, c54, c61, c62, c63, c64, 
	    c65, d21, d22, d23, d24, d25, d31, d32, d33, d34, d35;
    static integer md, mm;
    static doublereal sk, hd1, hd2, hd3, hd4;
    static integer nn2, nn3;
    static doublereal fac, hc21, hc31, hc32, hc41, hc42, hc43, hc51, hc52, 
	    hc53, hc54, hc61, hc62;
    static integer ier, lrc;
    static doublereal hc63, hc64, hc65, err, hacc;
    static integer lbeg, lend;
    static doublereal delt, hnew;
    static logical last;
    static doublereal hopt, gamma;
    extern /* Subroutine */ int rocoe_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal ysafe, hmaxn;
    static integer nsing;
    static doublereal xdelt;
    static integer irtrn;
    static doublereal erracc;
    static integer mujacj;
    extern /* Subroutine */ int decomr_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, logical *, integer *);
    static doublereal facgus;
    static logical reject;
    static integer mujacp;
    static doublereal posneg;
    extern /* Subroutine */ int slvrod_(integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, logical *);

    /* Fortran I/O blocks */
//    static cilist io___150 = { 0, 6, 0, fmt_979, 0 };
//    static cilist io___151 = { 0, 6, 0, 0, 0 };
//    static cilist io___152 = { 0, 6, 0, fmt_979, 0 };
//    static cilist io___153 = { 0, 6, 0, 0, 0 };
//    static cilist io___154 = { 0, 6, 0, fmt_979, 0 };
//    static cilist io___155 = { 0, 6, 0, 0, 0 };
//    static cilist io___156 = { 0, 6, 0, fmt_979, 0 };


/* ---------------------------------------------------------- */
/*     CORE INTEGRATOR FOR RODAS */
/*     PARAMETERS SAME AS IN RODAS WITH WORKSPACE ADDED */
/* ---------------------------------------------------------- */
/*         DECLARATIONS */
/* ---------------------------------------------------------- */
/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<        >*/
/*<       DIMENSION CONT(4*N) >*/
/*<       INTEGER IP(NM1) >*/
/*<       LOGICAL REJECT,AUTNMS,IMPLCT,BANDED,LAST,PRED >*/
/*<       COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG >*/
/*<       COMMON /CONROS/XOLD,HOUT,NN >*/
/* *** *** *** *** *** *** *** */
/*  INITIALISATIONS */
/* *** *** *** *** *** *** *** */
/*<       NN=N  >*/
    /* Parameter adjustments */
    --cont;
    --fx;
    --ak6;
    --ak5;
    --ak4;
    --ak3;
    --ak2;
    --ak1;
    --dy;
    --dy1;
    --ynew;
    --y;
    --rtol;
    --atol;
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    --rpar;
    --ipar;

    /* Function Body */
    conros_1.nn = *n;
/*<       NN2=2*N >*/
    nn2 = *n << 1;
/*<       NN3=3*N >*/
    nn3 = *n * 3;
/*<       LRC=4*N >*/
    lrc = *n << 2;
/* ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ---------- */
/*<       IF (IMPLCT) CALL MAS (NM1,FMAS,LDMAS,RPAR,IPAR) >*/
    if (*implct) {
	(*mas)(nm1, &fmas[fmas_offset], ldmas, &rpar[1], &ipar[1]);
    }
/* ------ SET THE PARAMETERS OF THE METHOD ----- */
/*<        >*/
    rocoe_(meth, &a21, &a31, &a32, &a41, &a42, &a43, &a51, &a52, &a53, &a54, &
	    c21, &c31, &c32, &c41, &c42, &c43, &c51, &c52, &c53, &c54, &c61, &
	    c62, &c63, &c64, &c65, &gamma, &c2, &c3, &c4, &d1, &d2, &d3, &d4, 
	    &d21, &d22, &d23, &d24, &d25, &d31, &d32, &d33, &d34, &d35);
/* --- INITIAL PREPARATIONS */
/*<       IF (M1.GT.0) IJOB=IJOB+10 >*/
    if (*m1 > 0) {
	*ijob += 10;
    }
/*<       POSNEG=SIGN(1.D0,XEND-X) >*/
    d__1 = *xend - *x;
    posneg = d_sign(&c_b61, &d__1);
/*<       HMAXN=MIN(ABS(HMAX),ABS(XEND-X)) >*/
/* Computing MIN */
    d__2 = abs(*hmax), d__3 = (d__1 = *xend - *x, abs(d__1));
    hmaxn = min(d__2,d__3);
/*<       IF (ABS(H).LE.10.D0*UROUND) H=1.0D-6 >*/
    if (abs(*h__) <= *uround * 10.) {
	*h__ = 1e-6;
    }
/*<       H=MIN(ABS(H),HMAXN)  >*/
/* Computing MIN */
    d__1 = abs(*h__);
    *h__ = min(d__1,hmaxn);
/*<       H=SIGN(H,POSNEG)  >*/
    *h__ = d_sign(h__, &posneg);
/*<       REJECT=.FALSE. >*/
    reject = FALSE_;
/*<       LAST=.FALSE. >*/
    last = FALSE_;
/*<       NSING=0 >*/
    nsing = 0;
/*<       IRTRN=1 >*/
    irtrn = 1;
/*<       IF (AUTNMS) THEN >*/
    if (*autnms) {
/*<          HD1=0.0D0 >*/
	hd1 = 0.;
/*<          HD2=0.0D0 >*/
	hd2 = 0.;
/*<          HD3=0.0D0 >*/
	hd3 = 0.;
/*<          HD4=0.0D0 >*/
	hd4 = 0.;
/*<       END IF >*/
    }
/* -------- PREPARE BAND-WIDTHS -------- */
/*<       MBDIAG=MUMAS+1 >*/
    linal_1.mbdiag = *mumas + 1;
/*<       IF (BANDED) THEN >*/
    if (*banded) {
/*<           MLE=MLJAC >*/
	linal_1.mle = *mljac;
/*<           MUE=MUJAC >*/
	linal_1.mue = *mujac;
/*<           MBJAC=MLJAC+MUJAC+1 >*/
	linal_1.mbjac = *mljac + *mujac + 1;
/*<           MBB=MLMAS+MUMAS+1 >*/
	linal_1.mbb = *mlmas + *mumas + 1;
/*<           MDIAG=MLE+MUE+1 >*/
	linal_1.mdiag = linal_1.mle + linal_1.mue + 1;
/*<           MDIFF=MLE+MUE-MUMAS >*/
	linal_1.mdiff = linal_1.mle + linal_1.mue - *mumas;
/*<       END IF >*/
    }
/*<       IF (IOUT.NE.0) THEN  >*/
    if (*iout != 0) {
/*<           XOLD=X >*/
	conros_1.xold = *x;
/*<           IRTRN=1 >*/
	irtrn = 1;
/*<           HOUT=H >*/
	conros_1.hout = *h__;
/*<        >*/
	i__1 = *naccpt + 1;
	(*solout)(&i__1, &conros_1.xold, x, &y[1], &cont[1], &lrc, n, &rpar[1]
		, &ipar[1], &irtrn);
/*<           IF (IRTRN.LT.0) GOTO 179 >*/
	if (irtrn < 0) {
	    goto L179;
	}
/*<       END IF >*/
    }
/* --- BASIC INTEGRATION STEP */
/*<  1    IF (NSTEP.GT.NMAX) GOTO 178 >*/
L1:
    if (*nstep > *nmax) {
	goto L178;
    }
/*<       IF (0.1D0*ABS(H).LE.ABS(X)*UROUND) GOTO 177 >*/
    if (abs(*h__) * .1 <= abs(*x) * *uround) {
	goto L177;
    }
/*<       IF (LAST) THEN >*/
    if (last) {
/*<           H=HOPT >*/
	*h__ = hopt;
/*<           IDID=1 >*/
	*idid = 1;
/*<           RETURN >*/
	return 0;
/*<       END IF >*/
    }
/*<       HOPT=H >*/
    hopt = *h__;
/*<       IF ((X+H*1.0001D0-XEND)*POSNEG.GE.0.D0) THEN >*/
    if ((*x + *h__ * 1.0001 - *xend) * posneg >= 0.) {
/*<          H=XEND-X >*/
	*h__ = *xend - *x;
/*<          LAST=.TRUE. >*/
	last = TRUE_;
/*<       END IF >*/
    }
/* *** *** *** *** *** *** *** */
/*  COMPUTATION OF THE JACOBIAN */
/* *** *** *** *** *** *** *** */
/*<       CALL FCN(N,X,Y,DY1,RPAR,IPAR) >*/
    (*fcn)(n, x, &y[1], &dy1[1], &rpar[1], &ipar[1]);
/*<       NFCN=NFCN+1 >*/
    ++(*nfcn);
/*<       NJAC=NJAC+1 >*/
    ++(*njac);
/*<       IF (IJAC.EQ.0) THEN >*/
    if (*ijac == 0) {
/* --- COMPUTE JACOBIAN MATRIX NUMERICALLY */
/*<          IF (BANDED) THEN >*/
	if (*banded) {
/* --- JACOBIAN IS BANDED */
/*<             MUJACP=MUJAC+1 >*/
	    mujacp = *mujac + 1;
/*<             MD=MIN(MBJAC,N) >*/
	    md = min(linal_1.mbjac,*n);
/*<             DO MM=1,M1/M2+1 >*/
	    i__1 = *m1 / *m2 + 1;
	    for (mm = 1; mm <= i__1; ++mm) {
/*<                DO K=1,MD >*/
		i__2 = md;
		for (k = 1; k <= i__2; ++k) {
/*<                   J=K+(MM-1)*M2 >*/
		    j = k + (mm - 1) * *m2;
/*<  12               AK2(J)=Y(J) >*/
L12:
		    ak2[j] = y[j];
/*<                   AK3(J)=DSQRT(UROUND*MAX(1.D-5,ABS(Y(J)))) >*/
/* Computing MAX */
		    d__2 = 1e-5, d__3 = (d__1 = y[j], abs(d__1));
		    ak3[j] = sqrt(*uround * max(d__2,d__3));
/*<                   Y(J)=Y(J)+AK3(J) >*/
		    y[j] += ak3[j];
/*<                   J=J+MD >*/
		    j += md;
/*<                   IF (J.LE.MM*M2) GOTO 12  >*/
		    if (j <= mm * *m2) {
			goto L12;
		    }
/*<                   CALL FCN(N,X,Y,AK1,RPAR,IPAR) >*/
		    (*fcn)(n, x, &y[1], &ak1[1], &rpar[1], &ipar[1]);
/*<                   J=K+(MM-1)*M2 >*/
		    j = k + (mm - 1) * *m2;
/*<                   J1=K >*/
		    j1 = k;
/*<                   LBEG=MAX(1,J1-MUJAC)+M1 >*/
/* Computing MAX */
		    i__3 = 1, i__4 = j1 - *mujac;
		    lbeg = max(i__3,i__4) + *m1;
/*<  14               LEND=MIN(M2,J1+MLJAC)+M1 >*/
L14:
/* Computing MIN */
		    i__3 = *m2, i__4 = j1 + *mljac;
		    lend = min(i__3,i__4) + *m1;
/*<                   Y(J)=AK2(J) >*/
		    y[j] = ak2[j];
/*<                   MUJACJ=MUJACP-J1-M1 >*/
		    mujacj = mujacp - j1 - *m1;
/*<                   DO L=LBEG,LEND >*/
		    i__3 = lend;
		    for (l = lbeg; l <= i__3; ++l) {
/*<                      FJAC(L+MUJACJ,J)=(AK1(L)-DY1(L))/AK3(J)  >*/
			fjac[l + mujacj + j * fjac_dim1] = (ak1[l] - dy1[l]) /
				 ak3[j];
/*<                   END DO >*/
		    }
/*<                   J=J+MD >*/
		    j += md;
/*<                   J1=J1+MD >*/
		    j1 += md;
/*<                   LBEG=LEND+1 >*/
		    lbeg = lend + 1;
/*<                   IF (J.LE.MM*M2) GOTO 14 >*/
		    if (j <= mm * *m2) {
			goto L14;
		    }
/*<                END DO >*/
		}
/*<             END DO >*/
	    }
/*<          ELSE >*/
	} else {
/* --- JACOBIAN IS FULL */
/*<             DO I=1,N >*/
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*<                YSAFE=Y(I) >*/
		ysafe = y[i__];
/*<                DELT=DSQRT(UROUND*MAX(1.D-5,ABS(YSAFE))) >*/
/* Computing MAX */
		d__1 = 1e-5, d__2 = abs(ysafe);
		delt = sqrt(*uround * max(d__1,d__2));
/*<                Y(I)=YSAFE+DELT >*/
		y[i__] = ysafe + delt;
/*<                CALL FCN(N,X,Y,AK1,RPAR,IPAR) >*/
		(*fcn)(n, x, &y[1], &ak1[1], &rpar[1], &ipar[1]);
/*<                DO J=M1+1,N >*/
		i__2 = *n;
		for (j = *m1 + 1; j <= i__2; ++j) {
/*<                  FJAC(J-M1,I)=(AK1(J)-DY1(J))/DELT >*/
		    fjac[j - *m1 + i__ * fjac_dim1] = (ak1[j] - dy1[j]) / 
			    delt;
/*<                END DO >*/
		}
/*<                Y(I)=YSAFE >*/
		y[i__] = ysafe;
/*<             END DO >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE >*/
    } else {
/* --- COMPUTE JACOBIAN MATRIX ANALYTICALLY */
/*<          CALL JAC(N,X,Y,FJAC,LDJAC,RPAR,IPAR) >*/
	(*jac)(n, x, &y[1], &fjac[fjac_offset], ldjac, &rpar[1], &ipar[1]);
/*<       END IF >*/
    }
/*<       IF (.NOT.AUTNMS) THEN >*/
    if (! (*autnms)) {
/*<          IF (IDFX.EQ.0) THEN >*/
	if (*idfx == 0) {
/* --- COMPUTE NUMERICALLY THE DERIVATIVE WITH RESPECT TO X */
/*<             DELT=SQRT(UROUND*MAX(1.D-5,ABS(X))) >*/
/* Computing MAX */
	    d__1 = 1e-5, d__2 = abs(*x);
	    delt = sqrt(*uround * max(d__1,d__2));
/*<             XDELT=X+DELT >*/
	    xdelt = *x + delt;
/*<             CALL FCN(N,XDELT,Y,AK1,RPAR,IPAR) >*/
	    (*fcn)(n, &xdelt, &y[1], &ak1[1], &rpar[1], &ipar[1]);
/*<             DO J=1,N >*/
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/*<                FX(J)=(AK1(J)-DY1(J))/DELT >*/
		fx[j] = (ak1[j] - dy1[j]) / delt;
/*<             END DO >*/
	    }
/*<          ELSE >*/
	} else {
/* --- COMPUTE ANALYTICALLY THE DERIVATIVE WITH RESPECT TO X */
/*<             CALL DFX(N,X,Y,FX,RPAR,IPAR) >*/
	    (*dfx)(n, x, &y[1], &fx[1], &rpar[1], &ipar[1]);
/*<          END IF >*/
	}
/*<       END IF >*/
    }
/*<    2  CONTINUE >*/
L2:
/* *** *** *** *** *** *** *** */
/*  COMPUTE THE STAGES */
/* *** *** *** *** *** *** *** */
/*<       FAC=1.D0/(H*GAMMA) >*/
    fac = 1. / (*h__ * gamma);
/*<        >*/
    decomr_(n, &fjac[fjac_offset], ldjac, &fmas[fmas_offset], ldmas, mlmas, 
	    mumas, m1, m2, nm1, &fac, &e[e_offset], lde, &ip[1], &ier, ijob, 
	    implct, &ip[1]);
/*<       IF (IER.NE.0) GOTO 80 >*/
    if (ier != 0) {
	goto L80;
    }
/*<       NDEC=NDEC+1 >*/
    ++(*ndec);
/* --- PREPARE FOR THE COMPUTATION OF THE 6 STAGES */
/*<       HC21=C21/H >*/
    hc21 = c21 / *h__;
/*<       HC31=C31/H >*/
    hc31 = c31 / *h__;
/*<       HC32=C32/H >*/
    hc32 = c32 / *h__;
/*<       HC41=C41/H >*/
    hc41 = c41 / *h__;
/*<       HC42=C42/H >*/
    hc42 = c42 / *h__;
/*<       HC43=C43/H >*/
    hc43 = c43 / *h__;
/*<       HC51=C51/H >*/
    hc51 = c51 / *h__;
/*<       HC52=C52/H >*/
    hc52 = c52 / *h__;
/*<       HC53=C53/H >*/
    hc53 = c53 / *h__;
/*<       HC54=C54/H >*/
    hc54 = c54 / *h__;
/*<       HC61=C61/H >*/
    hc61 = c61 / *h__;
/*<       HC62=C62/H >*/
    hc62 = c62 / *h__;
/*<       HC63=C63/H >*/
    hc63 = c63 / *h__;
/*<       HC64=C64/H >*/
    hc64 = c64 / *h__;
/*<       HC65=C65/H >*/
    hc65 = c65 / *h__;
/*<       IF (.NOT.AUTNMS) THEN >*/
    if (! (*autnms)) {
/*<          HD1=H*D1 >*/
	hd1 = *h__ * d1;
/*<          HD2=H*D2 >*/
	hd2 = *h__ * d2;
/*<          HD3=H*D3 >*/
	hd3 = *h__ * d3;
/*<          HD4=H*D4 >*/
	hd4 = *h__ * d4;
/*<       END IF >*/
    }
/* --- THE STAGES */
/*<        >*/
    slvrod_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset], 
	    ldmas, mlmas, mumas, m1, m2, nm1, &fac, &e[e_offset], lde, &ip[1],
	     &dy1[1], &ak1[1], &fx[1], &ynew[1], &hd1, ijob, &c_false);
/*<       DO  I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          YNEW(I)=Y(I)+A21*AK1(I) >*/
	ynew[i__] = y[i__] + a21 * ak1[i__];
/*<       END DO >*/
    }
/*<       CALL FCN(N,X+C2*H,YNEW,DY,RPAR,IPAR) >*/
    d__1 = *x + c2 * *h__;
    (*fcn)(n, &d__1, &ynew[1], &dy[1], &rpar[1], &ipar[1]);
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          YNEW(I)=HC21*AK1(I) >*/
	ynew[i__] = hc21 * ak1[i__];
/*<       END DO >*/
    }
/*<        >*/
    slvrod_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset], 
	    ldmas, mlmas, mumas, m1, m2, nm1, &fac, &e[e_offset], lde, &ip[1],
	     &dy[1], &ak2[1], &fx[1], &ynew[1], &hd2, ijob, &c_true);
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)  >*/
	ynew[i__] = y[i__] + a31 * ak1[i__] + a32 * ak2[i__];
/*<       END DO >*/
    }
/*<       CALL FCN(N,X+C3*H,YNEW,DY,RPAR,IPAR) >*/
    d__1 = *x + c3 * *h__;
    (*fcn)(n, &d__1, &ynew[1], &dy[1], &rpar[1], &ipar[1]);
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          YNEW(I)=HC31*AK1(I)+HC32*AK2(I) >*/
	ynew[i__] = hc31 * ak1[i__] + hc32 * ak2[i__];
/*<       END DO >*/
    }
/*<        >*/
    slvrod_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset], 
	    ldmas, mlmas, mumas, m1, m2, nm1, &fac, &e[e_offset], lde, &ip[1],
	     &dy[1], &ak3[1], &fx[1], &ynew[1], &hd3, ijob, &c_true);
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          YNEW(I)=Y(I)+A41*AK1(I)+A42*AK2(I)+A43*AK3(I) >*/
	ynew[i__] = y[i__] + a41 * ak1[i__] + a42 * ak2[i__] + a43 * ak3[i__];
/*<       END DO >*/
    }
/*<       CALL FCN(N,X+C4*H,YNEW,DY,RPAR,IPAR) >*/
    d__1 = *x + c4 * *h__;
    (*fcn)(n, &d__1, &ynew[1], &dy[1], &rpar[1], &ipar[1]);
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          YNEW(I)=HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I)  >*/
	ynew[i__] = hc41 * ak1[i__] + hc42 * ak2[i__] + hc43 * ak3[i__];
/*<       END DO >*/
    }
/*<        >*/
    slvrod_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset], 
	    ldmas, mlmas, mumas, m1, m2, nm1, &fac, &e[e_offset], lde, &ip[1],
	     &dy[1], &ak4[1], &fx[1], &ynew[1], &hd4, ijob, &c_true);
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          YNEW(I)=Y(I)+A51*AK1(I)+A52*AK2(I)+A53*AK3(I)+A54*AK4(I) >*/
	ynew[i__] = y[i__] + a51 * ak1[i__] + a52 * ak2[i__] + a53 * ak3[i__] 
		+ a54 * ak4[i__];
/*<       END DO >*/
    }
/*<       CALL FCN(N,X+H,YNEW,DY,RPAR,IPAR) >*/
    d__1 = *x + *h__;
    (*fcn)(n, &d__1, &ynew[1], &dy[1], &rpar[1], &ipar[1]);
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          AK6(I)=HC52*AK2(I)+HC54*AK4(I)+HC51*AK1(I)+HC53*AK3(I)  >*/
	ak6[i__] = hc52 * ak2[i__] + hc54 * ak4[i__] + hc51 * ak1[i__] + hc53 
		* ak3[i__];
/*<       END DO >*/
    }
/*<        >*/
    slvrod_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset], 
	    ldmas, mlmas, mumas, m1, m2, nm1, &fac, &e[e_offset], lde, &ip[1],
	     &dy[1], &ak5[1], &fx[1], &ak6[1], &c_b74, ijob, &c_true);
/* ------------ EMBEDDED SOLUTION --------------- */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          YNEW(I)=YNEW(I)+AK5(I)   >*/
	ynew[i__] += ak5[i__];
/*<       END DO >*/
    }
/*<       CALL FCN(N,X+H,YNEW,DY,RPAR,IPAR) >*/
    d__1 = *x + *h__;
    (*fcn)(n, &d__1, &ynew[1], &dy[1], &rpar[1], &ipar[1]);
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<        >*/
	cont[i__] = hc61 * ak1[i__] + hc62 * ak2[i__] + hc65 * ak5[i__] + 
		hc64 * ak4[i__] + hc63 * ak3[i__];
/*<       END DO >*/
    }
/*<        >*/
    slvrod_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset], 
	    ldmas, mlmas, mumas, m1, m2, nm1, &fac, &e[e_offset], lde, &ip[1],
	     &dy[1], &ak6[1], &fx[1], &cont[1], &c_b74, ijob, &c_true);
/* ------------ NEW SOLUTION --------------- */
/*<       DO  I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          YNEW(I)=YNEW(I)+AK6(I)   >*/
	ynew[i__] += ak6[i__];
/*<       END DO >*/
    }
/*<       NSOL=NSOL+6 >*/
    *nsol += 6;
/*<       NFCN=NFCN+5  >*/
    *nfcn += 5;
/* ------------ DENSE OUTPUT ---------- */
/*<       IF (IOUT.NE.0) THEN >*/
    if (*iout != 0) {
/*<          DO I=1,N  >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             CONT(I)=Y(I) >*/
	    cont[i__] = y[i__];
/*<        >*/
	    cont[i__ + nn2] = d21 * ak1[i__] + d22 * ak2[i__] + d23 * ak3[i__]
		     + d24 * ak4[i__] + d25 * ak5[i__];
/*<        >*/
	    cont[i__ + nn3] = d31 * ak1[i__] + d32 * ak2[i__] + d33 * ak3[i__]
		     + d34 * ak4[i__] + d35 * ak5[i__];
/*<          END DO >*/
	}
/*<       END IF >*/
    }
/* *** *** *** *** *** *** *** */
/*  ERROR ESTIMATION */
/* *** *** *** *** *** *** *** */
/*<       NSTEP=NSTEP+1 >*/
    ++(*nstep);
/* ------------ COMPUTE ERROR ESTIMATION ---------------- */
/*<       ERR=0.D0 >*/
    err = 0.;
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IF (ITOL.EQ.0) THEN >*/
	if (*itol == 0) {
/*<             SK=ATOL(1)+RTOL(1)*MAX(ABS(Y(I)),ABS(YNEW(I))) >*/
/* Computing MAX */
	    d__3 = (d__1 = y[i__], abs(d__1)), d__4 = (d__2 = ynew[i__], abs(
		    d__2));
	    sk = atol[1] + rtol[1] * max(d__3,d__4);
/*<          ELSE >*/
	} else {
/*<             SK=ATOL(I)+RTOL(I)*MAX(ABS(Y(I)),ABS(YNEW(I))) >*/
/* Computing MAX */
	    d__3 = (d__1 = y[i__], abs(d__1)), d__4 = (d__2 = ynew[i__], abs(
		    d__2));
	    sk = atol[i__] + rtol[i__] * max(d__3,d__4);
/*<          END IF >*/
	}
/*<          ERR=ERR+(AK6(I)/SK)**2 >*/
/* Computing 2nd power */
	d__1 = ak6[i__] / sk;
	err += d__1 * d__1;
/*<       END DO >*/
    }
/*<       ERR=SQRT(ERR/N) >*/
    err = sqrt(err / *n);
/* --- COMPUTATION OF HNEW */
/* --- WE REQUIRE .2<=HNEW/H<=6. */
/*<       FAC=MAX(FAC2,MIN(FAC1,(ERR)**0.25D0/SAFE)) >*/
/* Computing MAX */
/* Computing MIN */
    d__3 = *fac1, d__4 = pow_dd(&err, &c_b78) / *safe;
    d__1 = *fac2, d__2 = min(d__3,d__4);
    fac = max(d__1,d__2);
/*<       HNEW=H/FAC   >*/
    hnew = *h__ / fac;
/* *** *** *** *** *** *** *** */
/*  IS THE ERROR SMALL ENOUGH ? */
/* *** *** *** *** *** *** *** */
/*<       IF (ERR.LE.1.D0) THEN >*/
    if (err <= 1.) {
/* --- STEP IS ACCEPTED */
/*<          NACCPT=NACCPT+1 >*/
	++(*naccpt);
/*<          IF (PRED) THEN >*/
	if (*pred) {
/*       --- PREDICTIVE CONTROLLER OF GUSTAFSSON */
/*<             IF (NACCPT.GT.1) THEN >*/
	    if (*naccpt > 1) {
/*<                FACGUS=(HACC/H)*(ERR**2/ERRACC)**0.25D0/SAFE >*/
/* Computing 2nd power */
		d__2 = err;
		d__1 = d__2 * d__2 / erracc;
		facgus = hacc / *h__ * pow_dd(&d__1, &c_b78) / *safe;
/*<                FACGUS=MAX(FAC2,MIN(FAC1,FACGUS)) >*/
/* Computing MAX */
		d__1 = *fac2, d__2 = min(*fac1,facgus);
		facgus = max(d__1,d__2);
/*<                FAC=MAX(FAC,FACGUS) >*/
		fac = max(fac,facgus);
/*<                HNEW=H/FAC >*/
		hnew = *h__ / fac;
/*<             END IF >*/
	    }
/*<             HACC=H >*/
	    hacc = *h__;
/*<             ERRACC=MAX(1.0D-2,ERR) >*/
	    erracc = max(.01,err);
/*<          END IF >*/
	}
/*<          DO I=1,N  >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             Y(I)=YNEW(I) >*/
	    y[i__] = ynew[i__];
/*<          END DO >*/
	}
/*<          XOLD=X  >*/
	conros_1.xold = *x;
/*<          X=X+H >*/
	*x += *h__;
/*<          IF (IOUT.NE.0) THEN  >*/
	if (*iout != 0) {
/*<             DO I=1,N  >*/
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*<                CONT(NN+I)=Y(I) >*/
		cont[conros_1.nn + i__] = y[i__];
/*<             END DO >*/
	    }
/*<             IRTRN=1 >*/
	    irtrn = 1;
/*<             HOUT=H >*/
	    conros_1.hout = *h__;
/*<        >*/
	    i__1 = *naccpt + 1;
	    (*solout)(&i__1, &conros_1.xold, x, &y[1], &cont[1], &lrc, n, &
		    rpar[1], &ipar[1], &irtrn);
/*<             IF (IRTRN.LT.0) GOTO 179 >*/
	    if (irtrn < 0) {
		goto L179;
	    }
/*<          END IF >*/
	}
/*<          IF (ABS(HNEW).GT.HMAXN) HNEW=POSNEG*HMAXN >*/
	if (abs(hnew) > hmaxn) {
	    hnew = posneg * hmaxn;
	}
/*<          IF (REJECT) HNEW=POSNEG*MIN(ABS(HNEW),ABS(H))  >*/
	if (reject) {
/* Computing MIN */
	    d__1 = abs(hnew), d__2 = abs(*h__);
	    hnew = posneg * min(d__1,d__2);
	}
/*<          REJECT=.FALSE. >*/
	reject = FALSE_;
/*<          H=HNEW >*/
	*h__ = hnew;
/*<          GOTO 1 >*/
	goto L1;
/*<       ELSE >*/
    } else {
/* --- STEP IS REJECTED */
/*<          REJECT=.TRUE. >*/
	reject = TRUE_;
/*<          LAST=.FALSE. >*/
	last = FALSE_;
/*<          H=HNEW >*/
	*h__ = hnew;
/*<          IF (NACCPT.GE.1) NREJCT=NREJCT+1 >*/
	if (*naccpt >= 1) {
	    ++(*nrejct);
	}
/*<          GOTO 2 >*/
	goto L2;
/*<       END IF >*/
    }
/* --- SINGULAR MATRIX */
/*<   80  NSING=NSING+1 >*/
L80:
    ++nsing;
/*<       IF (NSING.GE.5) GOTO 176 >*/
    if (nsing >= 5) {
	goto L176;
    }
/*<       H=H*0.5D0 >*/
    *h__ *= .5;
/*<       REJECT=.TRUE. >*/
    reject = TRUE_;
/*<       LAST=.FALSE. >*/
    last = FALSE_;
/*<       GOTO 2 >*/
    goto L2;
/* --- FAIL EXIT */
/*<  176  CONTINUE >*/
L176:
/*<       WRITE(6,979)X    >*/
//    s_wsfe(&io___150);
//    do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
//    e_wsfe();
/*<       WRITE(6,*) ' MATRIX IS REPEATEDLY SINGULAR, IER=',IER >*/
//    s_wsle(&io___151);
//    do_lio(&c__9, &c__1, " MATRIX IS REPEATEDLY SINGULAR, IER=", (ftnlen)36);
//    do_lio(&c__3, &c__1, (char *)&ier, (ftnlen)sizeof(integer));
//    e_wsle();
     printf(" MATRIX IS REPEATEDLY SINGULAR, IER=%d\n",ier);
/*<       IDID=-4 >*/
    *idid = -4;
/*<       RETURN >*/
    return 0;
/*<  177  CONTINUE >*/
L177:
/*<       WRITE(6,979)X    >*/
//    s_wsfe(&io___152);
//    do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
//    e_wsfe();
/*<       WRITE(6,*) ' STEP SIZE T0O SMALL, H=',H >*/
//    s_wsle(&io___153);
//    do_lio(&c__9, &c__1, " STEP SIZE T0O SMALL, H=", (ftnlen)24);
//    do_lio(&c__5, &c__1, (char *)&(*h__), (ftnlen)sizeof(doublereal));
//    e_wsle();
     printf("STEP SIZE T0O SMALL, H=%e\n",h__);
/*<       IDID=-3 >*/
    *idid = -3;
/*<       RETURN >*/
    return 0;
/*<  178  CONTINUE >*/
L178:
/*<       WRITE(6,979)X    >*/
//    s_wsfe(&io___154);
//    do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
//    e_wsfe();
/*<       WRITE(6,*) ' MORE THAN NMAX =',NMAX,'STEPS ARE NEEDED'  >*/
//    s_wsle(&io___155);
//    do_lio(&c__9, &c__1, " MORE THAN NMAX =", (ftnlen)17);
//    do_lio(&c__3, &c__1, (char *)&(*nmax), (ftnlen)sizeof(integer));
//    do_lio(&c__9, &c__1, "STEPS ARE NEEDED", (ftnlen)16);
//    e_wsle();
       printf("MORE THAN NMAX =%d STEPS ARE NEEDED\n",nmax); 
/*<       IDID=-2 >*/
    *idid = -2;
/*<       RETURN >*/
    return 0;
/* --- EXIT CAUSED BY SOLOUT */
/*<  179  CONTINUE >*/
L179:
/*<       WRITE(6,979)X >*/
//    s_wsfe(&io___156);
//    do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
//    e_wsfe();
/*<  979  FORMAT(' EXIT OF RODAS AT X=',E18.4)  >*/
/*<       IDID=2 >*/
    *idid = 2;
/*<       RETURN >*/
    return 0;
/*<       END  >*/
} /* roscor_ */


/*<       FUNCTION CONTRO(I,X,CONT,LRC) >*/
doublereal contro_(integer *i__, doublereal *x, doublereal *cont, integer *
	lrc)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal s;

/* ---------------------------------------------------------- */
/*     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT IN CONNECTION */
/*     WITH THE OUTPUT-SUBROUTINE FOR RODAS. IT PROVIDES AN */
/*     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X. */
/* ---------------------------------------------------------- */
/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       DIMENSION CONT(LRC) >*/
/*<       COMMON /CONROS/XOLD,H,N >*/
/*<       S=(X-XOLD)/H  >*/
    /* Parameter adjustments */
    --cont;

    /* Function Body */
    s = (*x - conros_2.xold) / conros_2.h__;
/*<        >*/
    ret_val = cont[*i__] * (1 - s) + s * (cont[*i__ + conros_2.n] + (1 - s) * 
	    (cont[*i__ + (conros_2.n << 1)] + s * cont[*i__ + conros_2.n * 3])
	    );
/*<       RETURN >*/
    return ret_val;
/*<       END >*/
} /* contro_ */


/*<        >*/
/* Subroutine */ int rocoe_(integer *meth, doublereal *a21, doublereal *a31, 
	doublereal *a32, doublereal *a41, doublereal *a42, doublereal *a43, 
	doublereal *a51, doublereal *a52, doublereal *a53, doublereal *a54, 
	doublereal *c21, doublereal *c31, doublereal *c32, doublereal *c41, 
	doublereal *c42, doublereal *c43, doublereal *c51, doublereal *c52, 
	doublereal *c53, doublereal *c54, doublereal *c61, doublereal *c62, 
	doublereal *c63, doublereal *c64, doublereal *c65, doublereal *gamma, 
	doublereal *c2, doublereal *c3, doublereal *c4, doublereal *d1, 
	doublereal *d2, doublereal *d3, doublereal *d4, doublereal *d21, 
	doublereal *d22, doublereal *d23, doublereal *d24, doublereal *d25, 
	doublereal *d31, doublereal *d32, doublereal *d33, doublereal *d34, 
	doublereal *d35)
{
    static doublereal bet2p, bet3p, bet4p;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       GOTO (1,2,3), METH >*/
    switch (*meth) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
    }
/*<  1      C2=0.386D0 >*/
L1:
    *c2 = .386;
/*<         C3=0.21D0  >*/
    *c3 = .21;
/*<         C4=0.63D0 >*/
    *c4 = .63;
/*<         BET2P=0.0317D0 >*/
    bet2p = .0317;
/*<         BET3P=0.0635D0 >*/
    bet3p = .0635;
/*<         BET4P=0.3438D0  >*/
    bet4p = .3438;
/*<        D1= 0.2500000000000000D+00 >*/
    *d1 = .25;
/*<        D2=-0.1043000000000000D+00 >*/
    *d2 = -.1043;
/*<        D3= 0.1035000000000000D+00 >*/
    *d3 = .1035;
/*<        D4=-0.3620000000000023D-01 >*/
    *d4 = -.03620000000000023;
/*<        A21= 0.1544000000000000D+01 >*/
    *a21 = 1.544;
/*<        A31= 0.9466785280815826D+00 >*/
    *a31 = .9466785280815826;
/*<        A32= 0.2557011698983284D+00 >*/
    *a32 = .2557011698983284;
/*<        A41= 0.3314825187068521D+01 >*/
    *a41 = 3.314825187068521;
/*<        A42= 0.2896124015972201D+01 >*/
    *a42 = 2.896124015972201;
/*<        A43= 0.9986419139977817D+00 >*/
    *a43 = .9986419139977817;
/*<        A51= 0.1221224509226641D+01 >*/
    *a51 = 1.221224509226641;
/*<        A52= 0.6019134481288629D+01 >*/
    *a52 = 6.019134481288629;
/*<        A53= 0.1253708332932087D+02 >*/
    *a53 = 12.53708332932087;
/*<        A54=-0.6878860361058950D+00 >*/
    *a54 = -.687886036105895;
/*<        C21=-0.5668800000000000D+01 >*/
    *c21 = -5.6688;
/*<        C31=-0.2430093356833875D+01 >*/
    *c31 = -2.430093356833875;
/*<        C32=-0.2063599157091915D+00 >*/
    *c32 = -.2063599157091915;
/*<        C41=-0.1073529058151375D+00 >*/
    *c41 = -.1073529058151375;
/*<        C42=-0.9594562251023355D+01 >*/
    *c42 = -9.594562251023355;
/*<        C43=-0.2047028614809616D+02 >*/
    *c43 = -20.47028614809616;
/*<        C51= 0.7496443313967647D+01 >*/
    *c51 = 7.496443313967647;
/*<        C52=-0.1024680431464352D+02 >*/
    *c52 = -10.24680431464352;
/*<        C53=-0.3399990352819905D+02 >*/
    *c53 = -33.99990352819905;
/*<        C54= 0.1170890893206160D+02 >*/
    *c54 = 11.7089089320616;
/*<        C61= 0.8083246795921522D+01 >*/
    *c61 = 8.083246795921522;
/*<        C62=-0.7981132988064893D+01 >*/
    *c62 = -7.981132988064893;
/*<        C63=-0.3152159432874371D+02 >*/
    *c63 = -31.52159432874371;
/*<        C64= 0.1631930543123136D+02 >*/
    *c64 = 16.31930543123136;
/*<        C65=-0.6058818238834054D+01 >*/
    *c65 = -6.058818238834054;
/*<        GAMMA= 0.2500000000000000D+00   >*/
    *gamma = .25;
/*<        D21= 0.1012623508344586D+02 >*/
    *d21 = 10.12623508344586;
/*<        D22=-0.7487995877610167D+01 >*/
    *d22 = -7.487995877610167;
/*<        D23=-0.3480091861555747D+02 >*/
    *d23 = -34.80091861555747;
/*<        D24=-0.7992771707568823D+01 >*/
    *d24 = -7.992771707568823;
/*<        D25= 0.1025137723295662D+01 >*/
    *d25 = 1.025137723295662;
/*<        D31=-0.6762803392801253D+00 >*/
    *d31 = -.6762803392801253;
/*<        D32= 0.6087714651680015D+01 >*/
    *d32 = 6.087714651680015;
/*<        D33= 0.1643084320892478D+02 >*/
    *d33 = 16.43084320892478;
/*<        D34= 0.2476722511418386D+02 >*/
    *d34 = 24.76722511418386;
/*<        D35=-0.6594389125716872D+01 >*/
    *d35 = -6.594389125716872;
/*<       RETURN >*/
    return 0;

/*<  2      C2=0.3507221D0 >*/
L2:
    *c2 = .3507221;
/*<         C3=0.2557041D0  >*/
    *c3 = .2557041;
/*<         C4=0.6817790D0 >*/
    *c4 = .681779;
/*<         BET2P=0.0317D0 >*/
    bet2p = .0317;
/*<         BET3P=0.0047369D0 >*/
    bet3p = .0047369;
/*<         BET4P=0.3438D0  >*/
    bet4p = .3438;
/*<        D1= 0.2500000000000000D+00 >*/
    *d1 = .25;
/*<        D2=-0.6902209999999998D-01 >*/
    *d2 = -.06902209999999998;
/*<        D3=-0.9671999999999459D-03 >*/
    *d3 = -9.671999999999459e-4;
/*<        D4=-0.8797900000000025D-01 >*/
    *d4 = -.08797900000000025;
/*<        A21= 0.1402888400000000D+01 >*/
    *a21 = 1.4028884;
/*<        A31= 0.6581212688557198D+00 >*/
    *a31 = .6581212688557198;
/*<        A32=-0.1320936088384301D+01 >*/
    *a32 = -1.320936088384301;
/*<        A41= 0.7131197445744498D+01 >*/
    *a41 = 7.131197445744498;
/*<        A42= 0.1602964143958207D+02 >*/
    *a42 = 16.02964143958207;
/*<        A43=-0.5561572550509766D+01 >*/
    *a43 = -5.561572550509766;
/*<        A51= 0.2273885722420363D+02 >*/
    *a51 = 22.73885722420363;
/*<        A52= 0.6738147284535289D+02 >*/
    *a52 = 67.38147284535289;
/*<        A53=-0.3121877493038560D+02 >*/
    *a53 = -31.2187749303856;
/*<        A54= 0.7285641833203814D+00 >*/
    *a54 = .7285641833203814;
/*<        C21=-0.5104353600000000D+01 >*/
    *c21 = -5.1043536;
/*<        C31=-0.2899967805418783D+01 >*/
    *c31 = -2.899967805418783;
/*<        C32= 0.4040399359702244D+01 >*/
    *c32 = 4.040399359702244;
/*<        C41=-0.3264449927841361D+02 >*/
    *c41 = -32.64449927841361;
/*<        C42=-0.9935311008728094D+02 >*/
    *c42 = -99.35311008728094;
/*<        C43= 0.4999119122405989D+02 >*/
    *c43 = 49.99119122405989;
/*<        C51=-0.7646023087151691D+02 >*/
    *c51 = -76.46023087151691;
/*<        C52=-0.2785942120829058D+03 >*/
    *c52 = -278.5942120829058;
/*<        C53= 0.1539294840910643D+03 >*/
    *c53 = 153.9294840910643;
/*<        C54= 0.1097101866258358D+02 >*/
    *c54 = 10.97101866258358;
/*<        C61=-0.7629701586804983D+02 >*/
    *c61 = -76.29701586804983;
/*<        C62=-0.2942795630511232D+03 >*/
    *c62 = -294.2795630511232;
/*<        C63= 0.1620029695867566D+03 >*/
    *c63 = 162.0029695867566;
/*<        C64= 0.2365166903095270D+02 >*/
    *c64 = 23.6516690309527;
/*<        C65=-0.7652977706771382D+01 >*/
    *c65 = -7.652977706771382;
/*<        GAMMA= 0.2500000000000000D+00   >*/
    *gamma = .25;
/*<        D21=-0.3871940424117216D+02 >*/
    *d21 = -38.71940424117216;
/*<        D22=-0.1358025833007622D+03 >*/
    *d22 = -135.8025833007622;
/*<        D23= 0.6451068857505875D+02 >*/
    *d23 = 64.51068857505875;
/*<        D24=-0.4192663174613162D+01 >*/
    *d24 = -4.192663174613162;
/*<        D25=-0.2531932050335060D+01 >*/
    *d25 = -2.53193205033506;
/*<        D31=-0.1499268484949843D+02 >*/
    *d31 = -14.99268484949843;
/*<        D32=-0.7630242396627033D+02 >*/
    *d32 = -76.30242396627033;
/*<        D33= 0.5865928432851416D+02 >*/
    *d33 = 58.65928432851416;
/*<        D34= 0.1661359034616402D+02 >*/
    *d34 = 16.61359034616402;
/*<        D35=-0.6758691794084156D+00 >*/
    *d35 = -.6758691794084156;
/*<       RETURN >*/
    return 0;

/* Coefficients for RODAS with order 4 for linear parabolic problems */
/* Gerd Steinebach (1993) */
/*<   3     GAMMA = 0.25D0 >*/
L3:
    *gamma = .25;
/*< 	C2=3.d0*GAMMA >*/
    *c2 = *gamma * 3.;
/*<         C3=0.21D0  >*/
    *c3 = .21;
/*<         C4=0.63D0 >*/
    *c4 = .63;
/*< 	BET2P=0.D0 >*/
    bet2p = 0.;
/*< 	BET3P=c3*c3*(c3/6.d0-GAMMA/2.d0)/(GAMMA*GAMMA) >*/
    bet3p = *c3 * *c3 * (*c3 / 6. - *gamma / 2.) / (*gamma * *gamma);
/*<         BET4P=0.3438D0  >*/
    bet4p = .3438;
/*<        D1= 0.2500000000000000D+00 >*/
    *d1 = .25;
/*<        D2=-0.5000000000000000D+00 >*/
    *d2 = -.5;
/*<        D3=-0.2350400000000000D-01 >*/
    *d3 = -.023504;
/*<        D4=-0.3620000000000000D-01 >*/
    *d4 = -.0362;
/*<        A21= 0.3000000000000000D+01 >*/
    *a21 = 3.;
/*<        A31= 0.1831036793486759D+01 >*/
    *a31 = 1.831036793486759;
/*<        A32= 0.4955183967433795D+00 >*/
    *a32 = .4955183967433795;
/*<        A41= 0.2304376582692669D+01 >*/
    *a41 = 2.304376582692669;
/*<        A42=-0.5249275245743001D-01 >*/
    *a42 = -.05249275245743001;
/*<        A43=-0.1176798761832782D+01 >*/
    *a43 = -1.176798761832782;
/*<        A51=-0.7170454962423024D+01 >*/
    *a51 = -7.170454962423024;
/*<        A52=-0.4741636671481785D+01 >*/
    *a52 = -4.741636671481785;
/*<        A53=-0.1631002631330971D+02 >*/
    *a53 = -16.31002631330971;
/*<        A54=-0.1062004044111401D+01 >*/
    *a54 = -1.062004044111401;
/*<        C21=-0.1200000000000000D+02 >*/
    *c21 = -12.;
/*<        C31=-0.8791795173947035D+01 >*/
    *c31 = -8.791795173947035;
/*<        C32=-0.2207865586973518D+01 >*/
    *c32 = -2.207865586973518;
/*<        C41= 0.1081793056857153D+02 >*/
    *c41 = 10.81793056857153;
/*<        C42= 0.6780270611428266D+01 >*/
    *c42 = 6.780270611428266;
/*<        C43= 0.1953485944642410D+02 >*/
    *c43 = 19.5348594464241;
/*<        C51= 0.3419095006749676D+02 >*/
    *c51 = 34.19095006749676;
/*<        C52= 0.1549671153725963D+02 >*/
    *c52 = 15.49671153725963;
/*<        C53= 0.5474760875964130D+02 >*/
    *c53 = 54.7476087596413;
/*<        C54= 0.1416005392148534D+02 >*/
    *c54 = 14.16005392148534;
/*<        C61= 0.3462605830930532D+02 >*/
    *c61 = 34.62605830930532;
/*<        C62= 0.1530084976114473D+02 >*/
    *c62 = 15.30084976114473;
/*<        C63= 0.5699955578662667D+02 >*/
    *c63 = 56.99955578662667;
/*<        C64= 0.1840807009793095D+02 >*/
    *c64 = 18.40807009793095;
/*<        C65=-0.5714285714285717D+01 >*/
    *c65 = -5.714285714285717;

/*<        D21= 0.2509876703708589D+02 >*/
    *d21 = 25.09876703708589;
/*<        D22= 0.1162013104361867D+02 >*/
    *d22 = 11.62013104361867;
/*<        D23= 0.2849148307714626D+02 >*/
    *d23 = 28.49148307714626;
/*<        D24=-0.5664021568594133D+01 >*/
    *d24 = -5.664021568594133;
/*<        D25= 0.0000000000000000D+00 >*/
    *d25 = 0.;
/*<        D31= 0.1638054557396973D+01 >*/
    *d31 = 1.638054557396973;
/*<        D32=-0.7373619806678748D+00 >*/
    *d32 = -.7373619806678748;
/*<        D33= 0.8477918219238990D+01 >*/
    *d33 = 8.47791821923899;
/*<        D34= 0.1599253148779520D+02 >*/
    *d34 = 15.9925314877952;
/*<        D35=-0.1882352941176471D+01 >*/
    *d35 = -1.882352941176471;
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* rocoe_ */

