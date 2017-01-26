#include<stdio.h>
/* main.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__0 = 0;

/*     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC) */
/*     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS */
/*                     M*Y'=F(X,Y). */
/*     THE SYSTEM CAN BE (LINEARLY) IMPLICIT (MASS-MATRIX M .NE. I) */
/*     OR EXPLICIT (M=I). */
/*     THE METHOD USED IS AN IMPLICIT RUNGE-KUTTA METHOD (RADAU IIA) */
/*     OF ORDER 5 WITH STEP SIZE CONTROL AND CONTINUOUS OUTPUT. */
/*     CF. SECTION IV.8 */

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

/*     VERSION OF JULY 9, 1996 */
/*     (latest small correction: January 18, 2002) */

/*     INPUT PARAMETERS */
/*     ---------------- */
/*     N           DIMENSION OF THE SYSTEM */

/*     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE */
/*                 VALUE OF F(X,Y): */
/*                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR) */
/*                    DOUBLE PRECISION X,Y(N),F(N) */
/*                    F(1)=...   ETC. */
/*                 RPAR, IPAR (SEE BELOW) */

/*     X           INITIAL X-VALUE */

/*     Y(N)        INITIAL VALUES FOR Y */

/*     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE) */

/*     H           INITIAL STEP SIZE GUESS; */
/*                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT, */
/*                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD. */
/*                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS */
/*                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6). */

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
/*                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS */
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
/*                    IS SET <0, RADAU5 RETURNS TO THE CALLING PROGRAM. */

/*          -----  CONTINUOUS OUTPUT: ----- */
/*                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION */
/*                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH */
/*                 THE FUNCTION */
/*                        >>>   CONTR5(I,S,CONT,LRC)   <<< */
/*                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH */
/*                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE */
/*                 S SHOULD LIE IN THE INTERVAL [XOLD,X]. */
/*                 DO NOT CHANGE THE ENTRIES OF CONT(LRC), IF THE */
/*                 DENSE OUTPUT FUNCTION IS USED. */

/*     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT: */
/*                    IOUT=0: SUBROUTINE IS NEVER CALLED */
/*                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT. */

/*     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK". */
/*                 WORK(1), WORK(2),.., WORK(20) SERVE AS PARAMETERS */
/*                 FOR THE CODE. FOR STANDARD USE OF THE CODE */
/*                 WORK(1),..,WORK(20) MUST BE SET TO ZERO BEFORE */
/*                 CALLING. SEE BELOW FOR A MORE SOPHISTICATED USE. */
/*                 WORK(21),..,WORK(LWORK) SERVE AS WORKING SPACE */
/*                 FOR ALL VECTORS AND MATRICES. */
/*                 "LWORK" MUST BE AT LEAST */
/*                             N*(LJAC+LMAS+3*LE+12)+20 */
/*                 WHERE */
/*                    LJAC=N              IF MLJAC=N (FULL JACOBIAN) */
/*                    LJAC=MLJAC+MUJAC+1  IF MLJAC<N (BANDED JAC.) */
/*                 AND */
/*                    LMAS=0              IF IMAS=0 */
/*                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL) */
/*                    LMAS=MLMAS+MUMAS+1  IF MLMAS<N (BANDED MASS-M.) */
/*                 AND */
/*                    LE=N               IF MLJAC=N (FULL JACOBIAN) */
/*                    LE=2*MLJAC+MUJAC+1 IF MLJAC<N (BANDED JAC.) */

/*                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE */
/*                 MASS-MATRIX IS THE INDENTITY (IMAS=0), THE MINIMUM */
/*                 STORAGE REQUIREMENT IS */
/*                             LWORK = 4*N*N+12*N+20. */
/*                 IF IWORK(9)=M1>0 THEN "LWORK" MUST BE AT LEAST */
/*                          N*(LJAC+12)+(N-M1)*(LMAS+3*LE)+20 */
/*                 WHERE IN THE DEFINITIONS OF LJAC, LMAS AND LE THE */
/*                 NUMBER N CAN BE REPLACED BY N-M1. */

/*     LWORK       DECLARED LENGTH OF ARRAY "WORK". */

/*     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK". */
/*                 IWORK(1),IWORK(2),...,IWORK(20) SERVE AS PARAMETERS */
/*                 FOR THE CODE. FOR STANDARD USE, SET IWORK(1),.., */
/*                 IWORK(20) TO ZERO BEFORE CALLING. */
/*                 IWORK(21),...,IWORK(LIWORK) SERVE AS WORKING AREA. */
/*                 "LIWORK" MUST BE AT LEAST 3*N+20. */

/*     LIWORK      DECLARED LENGTH OF ARRAY "IWORK". */

/*     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH */
/*                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING */
/*                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES. */

/* ---------------------------------------------------------------------- */

/*     SOPHISTICATED SETTING OF PARAMETERS */
/*     ----------------------------------- */
/*              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK */
/*              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),... */
/*              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO. */
/*              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES: */

/*    IWORK(1)  IF IWORK(1).NE.0, THE CODE TRANSFORMS THE JACOBIAN */
/*              MATRIX TO HESSENBERG FORM. THIS IS PARTICULARLY */
/*              ADVANTAGEOUS FOR LARGE SYSTEMS WITH FULL JACOBIAN. */
/*              IT DOES NOT WORK FOR BANDED JACOBIAN (MLJAC<N) */
/*              AND NOT FOR IMPLICIT SYSTEMS (IMAS=1). */

/*    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS. */
/*              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000. */

/*    IWORK(3)  THE MAXIMUM NUMBER OF NEWTON ITERATIONS FOR THE */
/*              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP. */
/*              THE DEFAULT VALUE (FOR IWORK(3)=0) IS 7. */

/*    IWORK(4)  IF IWORK(4).EQ.0 THE EXTRAPOLATED COLLOCATION SOLUTION */
/*              IS TAKEN AS STARTING VALUE FOR NEWTON'S METHOD. */
/*              IF IWORK(4).NE.0 ZERO STARTING VALUES ARE USED. */
/*              THE LATTER IS RECOMMENDED IF NEWTON'S METHOD HAS */
/*              DIFFICULTIES WITH CONVERGENCE (THIS IS THE CASE WHEN */
/*              NSTEP IS LARGER THAN NACCPT + NREJCT; SEE OUTPUT PARAM.). */
/*              DEFAULT IS IWORK(4)=0. */

/*       THE FOLLOWING 3 PARAMETERS ARE IMPORTANT FOR */
/*       DIFFERENTIAL-ALGEBRAIC SYSTEMS OF INDEX > 1. */
/*       THE FUNCTION-SUBROUTINE SHOULD BE WRITTEN SUCH THAT */
/*       THE INDEX 1,2,3 VARIABLES APPEAR IN THIS ORDER. */
/*       IN ESTIMATING THE ERROR THE INDEX 2 VARIABLES ARE */
/*       MULTIPLIED BY H, THE INDEX 3 VARIABLES BY H**2. */

/*    IWORK(5)  DIMENSION OF THE INDEX 1 VARIABLES (MUST BE > 0). FOR */
/*              ODE'S THIS EQUALS THE DIMENSION OF THE SYSTEM. */
/*              DEFAULT IWORK(5)=N. */

/*    IWORK(6)  DIMENSION OF THE INDEX 2 VARIABLES. DEFAULT IWORK(6)=0. */

/*    IWORK(7)  DIMENSION OF THE INDEX 3 VARIABLES. DEFAULT IWORK(7)=0. */

/*    IWORK(8)  SWITCH FOR STEP SIZE STRATEGY */
/*              IF IWORK(8).EQ.1  MOD. PREDICTIVE CONTROLLER (GUSTAFSSON) */
/*              IF IWORK(8).EQ.2  CLASSICAL STEP SIZE CONTROL */
/*              THE DEFAULT VALUE (FOR IWORK(8)=0) IS IWORK(8)=1. */
/*              THE CHOICE IWORK(8).EQ.1 SEEMS TO PRODUCE SAFER RESULTS; */
/*              FOR SIMPLE PROBLEMS, THE CHOICE IWORK(8).EQ.2 PRODUCES */
/*              OFTEN SLIGHTLY FASTER RUNS */

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

/* ---------- */

/*    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16. */

/*    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION, */
/*              DEFAULT 0.9D0. */

/*    WORK(3)   DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; */
/*              INCREASE WORK(3), TO 0.1 SAY, WHEN JACOBIAN EVALUATIONS */
/*              ARE COSTLY. FOR SMALL SYSTEMS WORK(3) SHOULD BE SMALLER */
/*              (0.001D0, SAY). NEGATIV WORK(3) FORCES THE CODE TO */
/*              COMPUTE THE JACOBIAN AFTER EVERY ACCEPTED STEP. */
/*              DEFAULT 0.001D0. */

/*    WORK(4)   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1. */
/*              SMALLER VALUES OF WORK(4) MAKE THE CODE SLOWER, BUT SAFER. */
/*              DEFAULT MIN(0.03D0,RTOL(1)**0.5D0) */

/*    WORK(5) AND WORK(6) : IF WORK(5) < HNEW/HOLD < WORK(6), THEN THE */
/*              STEP SIZE IS NOT CHANGED. THIS SAVES, TOGETHER WITH A */
/*              LARGE WORK(3), LU-DECOMPOSITIONS AND COMPUTING TIME FOR */
/*              LARGE SYSTEMS. FOR SMALL SYSTEMS ONE MAY HAVE */
/*              WORK(5)=1.D0, WORK(6)=1.2D0, FOR LARGE FULL SYSTEMS */
/*              WORK(5)=0.99D0, WORK(6)=2.D0 MIGHT BE GOOD. */
/*              DEFAULTS WORK(5)=1.D0, WORK(6)=1.2D0 . */

/*    WORK(7)   MAXIMAL STEP SIZE, DEFAULT XEND-X. */

/*    WORK(8), WORK(9)   PARAMETERS FOR STEP SIZE SELECTION */
/*              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION */
/*                 WORK(8) <= HNEW/HOLD <= WORK(9) */
/*              DEFAULT VALUES: WORK(8)=0.2D0, WORK(9)=8.D0 */

/* ----------------------------------------------------------------------- */

/*     OUTPUT PARAMETERS */
/*     ----------------- */
/*     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED */
/*                 (AFTER SUCCESSFUL RETURN X=XEND). */

/*     Y(N)        NUMERICAL SOLUTION AT X */

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
/*   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS OF BOTH MATRICES */
/*   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS, OF BOTH */
/*                      SYSTEMS; THE NSTEP FORWARD-BACKWARD SUBSTITUTIONS, */
/*                      NEEDED FOR STEP SIZE SELECTION, ARE NOT COUNTED */
/* ----------------------------------------------------------------------- */
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/*          DECLARATIONS */
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/*<       SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR) >*/
/* Subroutine */ int mas_(integer *n, doublereal *am, integer *lmas, real *
	rpar, integer *ipar)
{
    /* System generated locals */
    integer am_dim1, am_offset;

    /* Builtin functions */
//    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
//	    e_wsle(void);

    /* Fortran I/O blocks */
//    static cilist io___1 = { 0, 6, 0, 0, 0 };


/*<       DOUBLE PRECISION AM(LMAS,N) >*/
/*<        write(*,*) 'Call of MASS is a mistake' >*/
    /* Parameter adjustments */
    am_dim1 = *lmas;
    am_offset = 1 + am_dim1;
    am -= am_offset;

    /* Function Body */
//    s_wsle(&io___1);
//    do_lio(&c__9, &c__1, "Call of MASS is a mistake", (ftnlen)25);
//    e_wsle();
    printf("Call of MASS is a mistake\n");
/*<       END  >*/
    return 0;
} /* mas_ */

/*<       SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR) >*/
/* Subroutine */ int fcn_(integer *n, doublereal *x, doublereal *y, 
	doublereal *f, real *rpar, integer *ipar)
{
    /* Builtin functions */
    double exp(doublereal);

/*<       DOUBLE PRECISION X,Y(N),F(N) >*/
/*<       F(1)=-Y(1) >*/
    /* Parameter adjustments */
    --f;
    --y;

    /* Function Body */
    f[1] = -y[1];
/*<       F(2)=exp(-Y(2)) >*/
    f[2] = exp(-y[2]);
/*<       END >*/
    return 0;
} /* fcn_ */

/*<       SUBROUTINE JAC(N,X,Y,DFY,LDFY,RPAR,IPAR) >*/
/* Subroutine */ int jac_(integer *n, doublereal *x, doublereal *y, 
	doublereal *dfy, integer *ldfy, real *rpar, integer *ipar)
{
    /* System generated locals */
    integer dfy_dim1, dfy_offset;

    /* Builtin functions */
//    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
//	    e_wsle(void);

    /* Fortran I/O blocks */
//    static cilist io___2 = { 0, 6, 0, 0, 0 };


/*<       DOUBLE PRECISION X,Y(N),DFY(LDFY,N) >*/
/*<       write(*,*) 'Call of JAC is a mistake' >*/
    /* Parameter adjustments */
    --y;
    dfy_dim1 = *ldfy;
    dfy_offset = 1 + dfy_dim1;
    dfy -= dfy_offset;

    /* Function Body */
//    s_wsle(&io___2);
//    do_lio(&c__9, &c__1, "Call of JAC is a mistake", (ftnlen)24);
//    e_wsle();
      printf("Call of JAC is a mistake\n");
/*<       END  >*/
    return 0;
} /* jac_ */

/*<        >*/
/* Subroutine */ int solout_(integer *nr, real *xold, doublereal *x, 
	doublereal *y, doublereal *cont, integer *lrc, integer *n, real *rpar,
	 integer *ipar, integer *irtrn)
{
    /* Builtin functions */
//    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
//	    e_wsle(void);

    /* Fortran I/O blocks */
//    static cilist io___3 = { 0, 6, 0, 0, 0 };


/*<       DOUBLE PRECISION X,Y(N),CONT(LRC) >*/
/*<       write(*,*) 'Call of SOLOUT is a mistake'             >*/
    /* Parameter adjustments */
    --cont;
    --y;

    /* Function Body */
//    s_wsle(&io___3);
//    do_lio(&c__9, &c__1, "Call of SOLOUT is a mistake", (ftnlen)27);
//    e_wsle();
       printf("Call of SOLOUT is a mistake\n");
/*<       END >*/
    return 0;
} /* solout_ */

/*<       integer function  smallRadu2(N,FCN,X,Y,XEND,H,eps) >*/
integer smallradu2_(integer *n, U_fp fcn, doublereal *x, doublereal *y, 
	doublereal *xend, doublereal *h__, doublereal *eps)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int jac_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, real *, integer *), mas_(integer *, 
	    doublereal *, integer *, real *, integer *);
    static integer ijac, idid, imas, ipar[5];
    static doublereal atol[2], rpar[5];
    static integer itol, iout;
    static doublereal rtol[2], work[100];
    static integer mljac, mujac;
    extern /* Subroutine */ int rodas_(integer *, U_fp, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, S_fp, integer *, integer *, integer *, 
	    U_fp, integer *, S_fp, integer *, integer *, integer *, S_fp, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *);
    static integer mlmas, mumas, iwork[100], lwork, liwork;
    extern /* Subroutine */ int solout_(integer *, real *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, real *, integer 
	    *, integer *);

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER IDID,N,LWORK,IWORK,IMAS,MLMAS,MUMAS,ITOL,IOUT >*/
/*<       DIMENSION  Y(2),ATOL(2),RTOL(2),WORK(100),IWORK(100) >*/
/*<       DIMENSION RPAR(5),IPAR(5) >*/
/*<       LOGICAL   IMPLCT,JBAND,ARRET,STARTN,PRED >*/
/*<       INTEGER I >*/
/*<       EXTERNAL FCN,JAC,MAS,SOLOUT >*/
/* DIMENSION */
/* PRECISION */
/*<       ITOL=0 >*/
    /* Parameter adjustments */
    --y;

    /* Function Body */
    itol = 0;
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; (int)i__<=(int)i__1; ++i__) {
    
/*<         ATOL(i)=Y(i)*eps*0.001 >*/
	atol[i__ - 1] = y[i__] * *eps * .001f;
/*<         RTOL(i)=eps >*/
	rtol[i__ - 1] = *eps;
/*<       ENDDO >*/
    }
/* MASS MATRIX */
/*<       IMAS=0 >*/
    imas = 0;
/*<       MLMAS=0 >*/
    mlmas = 0;
/*<       MUMAS=0 >*/
    mumas = 0;
/* JACOBIAN */
/*<       IJAC=0 >*/
    ijac = 0;
/*<       MLJAC=0 >*/
    mljac = 0;
/*<       MUJAC=0 >*/
    mujac = 0;
/*<       LWORK=100 >*/
    lwork = 100;
/*<       LIWORK=100  >*/
    liwork = 100;
/*<       DO I=1,LIWORK >*/
    i__1 = liwork;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         IWORK(I)=0  >*/
	iwork[i__ - 1] = 0;
/*<       ENDDO >*/
    }
/*<       DO I=1,LWORK >*/
    i__1 = lwork;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         WORK(I)=0 >*/
	work[i__ - 1] = 0.;
/*<       ENDDO >*/
    }
/*      Y(1)=1 */
/*      Y(2)=1 */
/*      X=0 */
/*      XEND=1 */
/*      H=0.1 */
/*<       IOUT=0 >*/
    iout = 0;
/*      CALL    RADAU5(N,FCN,X,Y,XEND,H, */
/*     &                  RTOL,ATOL,ITOL, */
/*     &                  JAC ,IJAC,MLJAC,MUJAC, */
/*     &                  MAS ,IMAS,MLMAS,MUMAS, */
/*     &                  SOLOUT,IOUT, */
/*     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID) */
/*<        >*/
    rodas_(n, (U_fp)fcn, &c__1, x, &y[1], xend, h__, rtol, atol, &itol, (S_fp)
	    jac_, &ijac, &mljac, &mujac, (U_fp)fcn, &c__0, (S_fp)mas_, &imas, 
	    &mlmas, &mumas, (S_fp)solout_, &iout, work, &lwork, iwork, &
	    liwork, rpar, ipar, &idid);
/*       write(*,*) 'IDID=',IDID */
/*      write(*,*) 'Y(1)=',Y(1),' Y(2)=',Y(2) */
/*<        smallRadu2=IDID >*/
    ret_val = idid;
/*<       END  >*/
    return ret_val;
} /* smallradu2_ */

