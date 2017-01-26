/* dc_decsol.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static integer c__1 = 1;

/* ****************************************** */
/*     VERSION OF SEPTEMBER 18, 1995 */
/* ****************************************** */

/*<        >*/
/* Subroutine */ int decomr_(integer *n, doublereal *fjac, integer *ldjac, 
	doublereal *fmas, integer *ldmas, integer *mlmas, integer *mumas, 
	integer *m1, integer *m2, integer *nm1, doublereal *fac1, doublereal *
	e1, integer *lde1, integer *ip1, integer *ier, integer *ijob, logical 
	*calhes, integer *iphes)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset, i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, j, k, j1, ib, mm, jm1;
    extern /* Subroutine */ int dec_(integer *, integer *, doublereal *, 
	    integer *, integer *);
    static doublereal sum;
    extern /* Subroutine */ int decb_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *), dech_(integer *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    elmhes_(integer *, integer *, integer *, integer *, doublereal *, 
	    integer *);

/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<        >*/
/*<       LOGICAL CALHES >*/
/*<       COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG >*/

/*<       GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB >*/
    /* Parameter adjustments */
    --iphes;
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip1;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    e1_dim1 = *lde1;
    e1_offset = 1 + e1_dim1;
    e1 -= e1_offset;

    /* Function Body */
    switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L14;
	case 15:  goto L15;
    }

/* ----------------------------------------------------------- */

/*<    1  CONTINUE >*/
L1:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
/*<       DO J=1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          DO  I=1,N >*/
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             E1(I,J)=-FJAC(I,J) >*/
	    e1[i__ + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
/*<          END DO >*/
	}
/*<          E1(J,J)=E1(J,J)+FAC1 >*/
	e1[j + j * e1_dim1] += *fac1;
/*<       END DO >*/
    }
/*<       CALL DEC (N,LDE1,E1,IP1,IER) >*/
    dec_(n, lde1, &e1[e1_offset], &ip1[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   11  CONTINUE >*/
L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO J=1,NM1 >*/
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
/*<          JM1=J+M1 >*/
	jm1 = j + *m1;
/*<          DO I=1,NM1 >*/
	i__2 = *nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             E1(I,J)=-FJAC(I,JM1) >*/
	    e1[i__ + j * e1_dim1] = -fjac[i__ + jm1 * fjac_dim1];
/*<          END DO >*/
	}
/*<          E1(J,J)=E1(J,J)+FAC1 >*/
	e1[j + j * e1_dim1] += *fac1;
/*<       END DO >*/
    }
/*<  45   MM=M1/M2 >*/
L45:
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=1,NM1 >*/
	i__2 = *nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             SUM=0.D0 >*/
	    sum = 0.;
/*<             DO K=0,MM-1 >*/
	    i__3 = mm - 1;
	    for (k = 0; k <= i__3; ++k) {
/*<                SUM=(SUM+FJAC(I,J+K*M2))/FAC1 >*/
		sum = (sum + fjac[i__ + (j + k * *m2) * fjac_dim1]) / *fac1;
/*<             END DO >*/
	    }
/*<             E1(I,J)=E1(I,J)-SUM >*/
	    e1[i__ + j * e1_dim1] -= sum;
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL DEC (NM1,LDE1,E1,IP1,IER) >*/
    dec_(nm1, lde1, &e1[e1_offset], &ip1[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    2  CONTINUE >*/
L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
/*<       DO J=1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=1,MBJAC >*/
	i__2 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             E1(I+MLE,J)=-FJAC(I,J) >*/
	    e1[i__ + linal_1.mle + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
/*<          END DO >*/
	}
/*<          E1(MDIAG,J)=E1(MDIAG,J)+FAC1 >*/
	e1[linal_1.mdiag + j * e1_dim1] += *fac1;
/*<       END DO >*/
    }
/*<       CALL DECB (N,LDE1,E1,MLE,MUE,IP1,IER) >*/
    decb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   12  CONTINUE >*/
L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
/*<       DO J=1,NM1 >*/
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
/*<          JM1=J+M1 >*/
	jm1 = j + *m1;
/*<          DO I=1,MBJAC >*/
	i__2 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             E1(I+MLE,J)=-FJAC(I,JM1) >*/
	    e1[i__ + linal_1.mle + j * e1_dim1] = -fjac[i__ + jm1 * fjac_dim1]
		    ;
/*<          END DO >*/
	}
/*<          E1(MDIAG,J)=E1(MDIAG,J)+FAC1 >*/
	e1[linal_1.mdiag + j * e1_dim1] += *fac1;
/*<       END DO >*/
    }
/*<   46  MM=M1/M2 >*/
L46:
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=1,MBJAC >*/
	i__2 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             SUM=0.D0 >*/
	    sum = 0.;
/*<             DO K=0,MM-1 >*/
	    i__3 = mm - 1;
	    for (k = 0; k <= i__3; ++k) {
/*<                SUM=(SUM+FJAC(I,J+K*M2))/FAC1 >*/
		sum = (sum + fjac[i__ + (j + k * *m2) * fjac_dim1]) / *fac1;
/*<             END DO >*/
	    }
/*<             E1(I+MLE,J)=E1(I+MLE,J)-SUM >*/
	    e1[i__ + linal_1.mle + j * e1_dim1] -= sum;
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL DECB (NM1,LDE1,E1,MLE,MUE,IP1,IER) >*/
    decb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier)
	    ;
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    3  CONTINUE >*/
L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
/*<       DO J=1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=1,N >*/
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             E1(I,J)=-FJAC(I,J) >*/
	    e1[i__ + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
/*<          END DO >*/
	}
/*<          DO I=MAX(1,J-MUMAS),MIN(N,J+MLMAS) >*/
/* Computing MAX */
	i__2 = 1, i__3 = j - *mumas;
/* Computing MIN */
	i__5 = *n, i__6 = j + *mlmas;
	i__4 = min(i__5,i__6);
	for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
/*<             E1(I,J)=E1(I,J)+FAC1*FMAS(I-J+MBDIAG,J) >*/
	    e1[i__ + j * e1_dim1] += *fac1 * fmas[i__ - j + linal_1.mbdiag + 
		    j * fmas_dim1];
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL DEC (N,LDE1,E1,IP1,IER) >*/
    dec_(n, lde1, &e1[e1_offset], &ip1[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   13  CONTINUE >*/
L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO J=1,NM1 >*/
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
/*<          JM1=J+M1 >*/
	jm1 = j + *m1;
/*<          DO I=1,NM1 >*/
	i__4 = *nm1;
	for (i__ = 1; i__ <= i__4; ++i__) {
/*<             E1(I,J)=-FJAC(I,JM1) >*/
	    e1[i__ + j * e1_dim1] = -fjac[i__ + jm1 * fjac_dim1];
/*<          END DO >*/
	}
/*<          DO I=MAX(1,J-MUMAS),MIN(NM1,J+MLMAS) >*/
/* Computing MAX */
	i__4 = 1, i__2 = j - *mumas;
/* Computing MIN */
	i__5 = *nm1, i__6 = j + *mlmas;
	i__3 = min(i__5,i__6);
	for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
/*<             E1(I,J)=E1(I,J)+FAC1*FMAS(I-J+MBDIAG,J) >*/
	    e1[i__ + j * e1_dim1] += *fac1 * fmas[i__ - j + linal_1.mbdiag + 
		    j * fmas_dim1];
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       GOTO 45 >*/
    goto L45;

/* ----------------------------------------------------------- */

/*<    4  CONTINUE >*/
L4:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
/*<       DO J=1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=1,MBJAC >*/
	i__3 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<             E1(I+MLE,J)=-FJAC(I,J) >*/
	    e1[i__ + linal_1.mle + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
/*<          END DO >*/
	}
/*<          DO I=1,MBB >*/
	i__3 = linal_1.mbb;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<             IB=I+MDIFF >*/
	    ib = i__ + linal_1.mdiff;
/*<             E1(IB,J)=E1(IB,J)+FAC1*FMAS(I,J) >*/
	    e1[ib + j * e1_dim1] += *fac1 * fmas[i__ + j * fmas_dim1];
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL DECB (N,LDE1,E1,MLE,MUE,IP1,IER) >*/
    decb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   14  CONTINUE >*/
L14:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
/*<       DO J=1,NM1 >*/
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
/*<          JM1=J+M1 >*/
	jm1 = j + *m1;
/*<          DO I=1,MBJAC >*/
	i__3 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<             E1(I+MLE,J)=-FJAC(I,JM1) >*/
	    e1[i__ + linal_1.mle + j * e1_dim1] = -fjac[i__ + jm1 * fjac_dim1]
		    ;
/*<          END DO >*/
	}
/*<          DO I=1,MBB >*/
	i__3 = linal_1.mbb;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<             IB=I+MDIFF >*/
	    ib = i__ + linal_1.mdiff;
/*<             E1(IB,J)=E1(IB,J)+FAC1*FMAS(I,J) >*/
	    e1[ib + j * e1_dim1] += *fac1 * fmas[i__ + j * fmas_dim1];
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       GOTO 46 >*/
    goto L46;

/* ----------------------------------------------------------- */

/*<    5  CONTINUE >*/
L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
/*<       DO J=1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=1,N >*/
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<             E1(I,J)=FMAS(I,J)*FAC1-FJAC(I,J) >*/
	    e1[i__ + j * e1_dim1] = fmas[i__ + j * fmas_dim1] * *fac1 - fjac[
		    i__ + j * fjac_dim1];
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL DEC (N,LDE1,E1,IP1,IER) >*/
    dec_(n, lde1, &e1[e1_offset], &ip1[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   15  CONTINUE >*/
L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO J=1,NM1 >*/
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
/*<          JM1=J+M1 >*/
	jm1 = j + *m1;
/*<          DO I=1,NM1 >*/
	i__3 = *nm1;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<             E1(I,J)=FMAS(I,J)*FAC1-FJAC(I,JM1) >*/
	    e1[i__ + j * e1_dim1] = fmas[i__ + j * fmas_dim1] * *fac1 - fjac[
		    i__ + jm1 * fjac_dim1];
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       GOTO 45 >*/
    goto L45;

/* ----------------------------------------------------------- */

/*<    6  CONTINUE >*/
L6:
/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
/* ---  THIS OPTION IS NOT PROVIDED */
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    7  CONTINUE >*/
L7:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
/*<       IF (CALHES) CALL ELMHES (LDJAC,N,1,N,FJAC,IPHES)  >*/
    if (*calhes) {
	elmhes_(ldjac, n, &c__1, n, &fjac[fjac_offset], &iphes[1]);
    }
/*<       CALHES=.FALSE. >*/
    *calhes = FALSE_;
/*<       DO J=1,N-1 >*/
    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j) {
/*<          J1=J+1 >*/
	j1 = j + 1;
/*<          E1(J1,J)=-FJAC(J1,J) >*/
	e1[j1 + j * e1_dim1] = -fjac[j1 + j * fjac_dim1];
/*<       END DO >*/
    }
/*<       DO J=1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=1,J >*/
	i__3 = j;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<             E1(I,J)=-FJAC(I,J) >*/
	    e1[i__ + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
/*<          END DO >*/
	}
/*<          E1(J,J)=E1(J,J)+FAC1 >*/
	e1[j + j * e1_dim1] += *fac1;
/*<       END DO >*/
    }
/*<       CALL DECH(N,LDE1,E1,1,IP1,IER) >*/
    dech_(n, lde1, &e1[e1_offset], &c__1, &ip1[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   55  CONTINUE >*/
L55:
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* decomr_ */


/*     END OF SUBROUTINE DECOMR */

/* *********************************************************** */

/*<        >*/
/* Subroutine */ int decomc_(integer *n, doublereal *fjac, integer *ldjac, 
	doublereal *fmas, integer *ldmas, integer *mlmas, integer *mumas, 
	integer *m1, integer *m2, integer *nm1, doublereal *alphn, doublereal 
	*betan, doublereal *e2r, doublereal *e2i, integer *lde1, integer *ip2,
	 integer *ier, integer *ijob)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e2r_dim1, 
	    e2r_offset, e2i_dim1, e2i_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, j1;
    static doublereal bb;
    static integer ib, mm, jm1;
    static doublereal bet, alp;
    extern /* Subroutine */ int decc_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *);
    static doublereal ffma, abno;
    static integer imle;
    static doublereal sumi, sumr, sums;
    extern /* Subroutine */ int decbc_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *), dechc_(
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *);

/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<        >*/
/*<       COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG >*/

/*<       GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB >*/
    /* Parameter adjustments */
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip2;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    e2i_dim1 = *lde1;
    e2i_offset = 1 + e2i_dim1;
    e2i -= e2i_offset;
    e2r_dim1 = *lde1;
    e2r_offset = 1 + e2r_dim1;
    e2r -= e2r_offset;

    /* Function Body */
    switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L14;
	case 15:  goto L15;
    }

/* ----------------------------------------------------------- */

/*<    1  CONTINUE >*/
L1:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
/*<       DO J=1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=1,N >*/
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             E2R(I,J)=-FJAC(I,J) >*/
	    e2r[i__ + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
/*<             E2I(I,J)=0.D0 >*/
	    e2i[i__ + j * e2i_dim1] = 0.;
/*<          END DO >*/
	}
/*<          E2R(J,J)=E2R(J,J)+ALPHN >*/
	e2r[j + j * e2r_dim1] += *alphn;
/*<          E2I(J,J)=BETAN >*/
	e2i[j + j * e2i_dim1] = *betan;
/*<       END DO >*/
    }
/*<       CALL DECC (N,LDE1,E2R,E2I,IP2,IER) >*/
    decc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   11  CONTINUE >*/
L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO J=1,NM1 >*/
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
/*<          JM1=J+M1 >*/
	jm1 = j + *m1;
/*<          DO I=1,NM1 >*/
	i__2 = *nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             E2R(I,J)=-FJAC(I,JM1) >*/
	    e2r[i__ + j * e2r_dim1] = -fjac[i__ + jm1 * fjac_dim1];
/*<             E2I(I,J)=0.D0 >*/
	    e2i[i__ + j * e2i_dim1] = 0.;
/*<          END DO >*/
	}
/*<          E2R(J,J)=E2R(J,J)+ALPHN >*/
	e2r[j + j * e2r_dim1] += *alphn;
/*<          E2I(J,J)=BETAN >*/
	e2i[j + j * e2i_dim1] = *betan;
/*<       END DO >*/
    }
/*<   45  MM=M1/M2 >*/
L45:
    mm = *m1 / *m2;
/*<       ABNO=ALPHN**2+BETAN**2 >*/
/* Computing 2nd power */
    d__1 = *alphn;
/* Computing 2nd power */
    d__2 = *betan;
    abno = d__1 * d__1 + d__2 * d__2;
/*<       ALP=ALPHN/ABNO >*/
    alp = *alphn / abno;
/*<       BET=BETAN/ABNO >*/
    bet = *betan / abno;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=1,NM1 >*/
	i__2 = *nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             SUMR=0.D0 >*/
	    sumr = 0.;
/*<             SUMI=0.D0 >*/
	    sumi = 0.;
/*<             DO K=0,MM-1 >*/
	    i__3 = mm - 1;
	    for (k = 0; k <= i__3; ++k) {
/*<                SUMS=SUMR+FJAC(I,J+K*M2) >*/
		sums = sumr + fjac[i__ + (j + k * *m2) * fjac_dim1];
/*<                SUMR=SUMS*ALP+SUMI*BET >*/
		sumr = sums * alp + sumi * bet;
/*<                SUMI=SUMI*ALP-SUMS*BET >*/
		sumi = sumi * alp - sums * bet;
/*<             END DO >*/
	    }
/*<             E2R(I,J)=E2R(I,J)-SUMR >*/
	    e2r[i__ + j * e2r_dim1] -= sumr;
/*<             E2I(I,J)=E2I(I,J)-SUMI >*/
	    e2i[i__ + j * e2i_dim1] -= sumi;
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL DECC (NM1,LDE1,E2R,E2I,IP2,IER) >*/
    decc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    2  CONTINUE >*/
L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
/*<       DO J=1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=1,MBJAC >*/
	i__2 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             IMLE=I+MLE >*/
	    imle = i__ + linal_1.mle;
/*<             E2R(IMLE,J)=-FJAC(I,J) >*/
	    e2r[imle + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
/*<             E2I(IMLE,J)=0.D0 >*/
	    e2i[imle + j * e2i_dim1] = 0.;
/*<          END DO >*/
	}
/*<          E2R(MDIAG,J)=E2R(MDIAG,J)+ALPHN >*/
	e2r[linal_1.mdiag + j * e2r_dim1] += *alphn;
/*<          E2I(MDIAG,J)=BETAN >*/
	e2i[linal_1.mdiag + j * e2i_dim1] = *betan;
/*<       END DO >*/
    }
/*<       CALL DECBC (N,LDE1,E2R,E2I,MLE,MUE,IP2,IER) >*/
    decbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &ip2[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   12  CONTINUE >*/
L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
/*<       DO J=1,NM1 >*/
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
/*<          JM1=J+M1 >*/
	jm1 = j + *m1;
/*<          DO I=1,MBJAC >*/
	i__2 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             E2R(I+MLE,J)=-FJAC(I,JM1) >*/
	    e2r[i__ + linal_1.mle + j * e2r_dim1] = -fjac[i__ + jm1 * 
		    fjac_dim1];
/*<             E2I(I+MLE,J)=0.D0 >*/
	    e2i[i__ + linal_1.mle + j * e2i_dim1] = 0.;
/*<          END DO >*/
	}
/*<          E2R(MDIAG,J)=E2R(MDIAG,J)+ALPHN >*/
	e2r[linal_1.mdiag + j * e2r_dim1] += *alphn;
/*<          E2I(MDIAG,J)=E2I(MDIAG,J)+BETAN >*/
	e2i[linal_1.mdiag + j * e2i_dim1] += *betan;
/*<       END DO >*/
    }
/*<   46  MM=M1/M2 >*/
L46:
    mm = *m1 / *m2;
/*<       ABNO=ALPHN**2+BETAN**2 >*/
/* Computing 2nd power */
    d__1 = *alphn;
/* Computing 2nd power */
    d__2 = *betan;
    abno = d__1 * d__1 + d__2 * d__2;
/*<       ALP=ALPHN/ABNO >*/
    alp = *alphn / abno;
/*<       BET=BETAN/ABNO >*/
    bet = *betan / abno;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=1,MBJAC >*/
	i__2 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             SUMR=0.D0 >*/
	    sumr = 0.;
/*<             SUMI=0.D0 >*/
	    sumi = 0.;
/*<             DO K=0,MM-1 >*/
	    i__3 = mm - 1;
	    for (k = 0; k <= i__3; ++k) {
/*<                SUMS=SUMR+FJAC(I,J+K*M2) >*/
		sums = sumr + fjac[i__ + (j + k * *m2) * fjac_dim1];
/*<                SUMR=SUMS*ALP+SUMI*BET >*/
		sumr = sums * alp + sumi * bet;
/*<                SUMI=SUMI*ALP-SUMS*BET >*/
		sumi = sumi * alp - sums * bet;
/*<             END DO >*/
	    }
/*<             IMLE=I+MLE >*/
	    imle = i__ + linal_1.mle;
/*<             E2R(IMLE,J)=E2R(IMLE,J)-SUMR >*/
	    e2r[imle + j * e2r_dim1] -= sumr;
/*<             E2I(IMLE,J)=E2I(IMLE,J)-SUMI >*/
	    e2i[imle + j * e2i_dim1] -= sumi;
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL DECBC (NM1,LDE1,E2R,E2I,MLE,MUE,IP2,IER) >*/
    decbc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &ip2[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    3  CONTINUE >*/
L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
/*<       DO  J=1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          DO  I=1,N >*/
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             E2R(I,J)=-FJAC(I,J) >*/
	    e2r[i__ + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
/*<             E2I(I,J)=0.D0 >*/
	    e2i[i__ + j * e2i_dim1] = 0.;
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       DO J=1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=MAX(1,J-MUMAS),MIN(N,J+MLMAS) >*/
/* Computing MAX */
	i__2 = 1, i__3 = j - *mumas;
/* Computing MIN */
	i__5 = *n, i__6 = j + *mlmas;
	i__4 = min(i__5,i__6);
	for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
/*<             BB=FMAS(I-J+MBDIAG,J) >*/
	    bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
/*<             E2R(I,J)=E2R(I,J)+ALPHN*BB >*/
	    e2r[i__ + j * e2r_dim1] += *alphn * bb;
/*<             E2I(I,J)=BETAN*BB >*/
	    e2i[i__ + j * e2i_dim1] = *betan * bb;
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL DECC(N,LDE1,E2R,E2I,IP2,IER) >*/
    decc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   13  CONTINUE >*/
L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO J=1,NM1 >*/
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
/*<          JM1=J+M1 >*/
	jm1 = j + *m1;
/*<          DO I=1,NM1 >*/
	i__4 = *nm1;
	for (i__ = 1; i__ <= i__4; ++i__) {
/*<             E2R(I,J)=-FJAC(I,JM1) >*/
	    e2r[i__ + j * e2r_dim1] = -fjac[i__ + jm1 * fjac_dim1];
/*<             E2I(I,J)=0.D0 >*/
	    e2i[i__ + j * e2i_dim1] = 0.;
/*<          END DO >*/
	}
/*<          DO I=MAX(1,J-MUMAS),MIN(NM1,J+MLMAS) >*/
/* Computing MAX */
	i__4 = 1, i__2 = j - *mumas;
/* Computing MIN */
	i__5 = *nm1, i__6 = j + *mlmas;
	i__3 = min(i__5,i__6);
	for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
/*<             FFMA=FMAS(I-J+MBDIAG,J) >*/
	    ffma = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
/*<             E2R(I,J)=E2R(I,J)+ALPHN*FFMA >*/
	    e2r[i__ + j * e2r_dim1] += *alphn * ffma;
/*<             E2I(I,J)=E2I(I,J)+BETAN*FFMA >*/
	    e2i[i__ + j * e2i_dim1] += *betan * ffma;
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       GOTO 45 >*/
    goto L45;

/* ----------------------------------------------------------- */

/*<    4  CONTINUE >*/
L4:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
/*<       DO J=1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=1,MBJAC >*/
	i__3 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<             IMLE=I+MLE >*/
	    imle = i__ + linal_1.mle;
/*<             E2R(IMLE,J)=-FJAC(I,J) >*/
	    e2r[imle + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
/*<             E2I(IMLE,J)=0.D0 >*/
	    e2i[imle + j * e2i_dim1] = 0.;
/*<          END DO >*/
	}
/*<          DO I=MAX(1,MUMAS+2-J),MIN(MBB,MUMAS+1-J+N) >*/
/* Computing MAX */
	i__3 = 1, i__4 = *mumas + 2 - j;
/* Computing MIN */
	i__5 = linal_1.mbb, i__6 = *mumas + 1 - j + *n;
	i__2 = min(i__5,i__6);
	for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
/*<             IB=I+MDIFF >*/
	    ib = i__ + linal_1.mdiff;
/*<             BB=FMAS(I,J) >*/
	    bb = fmas[i__ + j * fmas_dim1];
/*<             E2R(IB,J)=E2R(IB,J)+ALPHN*BB >*/
	    e2r[ib + j * e2r_dim1] += *alphn * bb;
/*<             E2I(IB,J)=BETAN*BB >*/
	    e2i[ib + j * e2i_dim1] = *betan * bb;
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL DECBC (N,LDE1,E2R,E2I,MLE,MUE,IP2,IER) >*/
    decbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &ip2[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   14  CONTINUE >*/
L14:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
/*<       DO J=1,NM1 >*/
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
/*<          JM1=J+M1 >*/
	jm1 = j + *m1;
/*<          DO I=1,MBJAC >*/
	i__2 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             E2R(I+MLE,J)=-FJAC(I,JM1) >*/
	    e2r[i__ + linal_1.mle + j * e2r_dim1] = -fjac[i__ + jm1 * 
		    fjac_dim1];
/*<             E2I(I+MLE,J)=0.D0 >*/
	    e2i[i__ + linal_1.mle + j * e2i_dim1] = 0.;
/*<          END DO >*/
	}
/*<          DO I=1,MBB >*/
	i__2 = linal_1.mbb;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             IB=I+MDIFF >*/
	    ib = i__ + linal_1.mdiff;
/*<             FFMA=FMAS(I,J) >*/
	    ffma = fmas[i__ + j * fmas_dim1];
/*<             E2R(IB,J)=E2R(IB,J)+ALPHN*FFMA >*/
	    e2r[ib + j * e2r_dim1] += *alphn * ffma;
/*<             E2I(IB,J)=E2I(IB,J)+BETAN*FFMA >*/
	    e2i[ib + j * e2i_dim1] += *betan * ffma;
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       GOTO 46 >*/
    goto L46;

/* ----------------------------------------------------------- */

/*<    5  CONTINUE >*/
L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
/*<       DO J=1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=1,N >*/
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             BB=FMAS(I,J) >*/
	    bb = fmas[i__ + j * fmas_dim1];
/*<             E2R(I,J)=BB*ALPHN-FJAC(I,J) >*/
	    e2r[i__ + j * e2r_dim1] = bb * *alphn - fjac[i__ + j * fjac_dim1];
/*<             E2I(I,J)=BB*BETAN >*/
	    e2i[i__ + j * e2i_dim1] = bb * *betan;
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL DECC(N,LDE1,E2R,E2I,IP2,IER) >*/
    decc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   15  CONTINUE >*/
L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO J=1,NM1 >*/
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
/*<          JM1=J+M1 >*/
	jm1 = j + *m1;
/*<          DO I=1,NM1 >*/
	i__2 = *nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             E2R(I,J)=ALPHN*FMAS(I,J)-FJAC(I,JM1) >*/
	    e2r[i__ + j * e2r_dim1] = *alphn * fmas[i__ + j * fmas_dim1] - 
		    fjac[i__ + jm1 * fjac_dim1];
/*<             E2I(I,J)=BETAN*FMAS(I,J) >*/
	    e2i[i__ + j * e2i_dim1] = *betan * fmas[i__ + j * fmas_dim1];
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       GOTO 45 >*/
    goto L45;

/* ----------------------------------------------------------- */

/*<    6  CONTINUE >*/
L6:
/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
/* ---  THIS OPTION IS NOT PROVIDED */
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    7  CONTINUE >*/
L7:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
/*<       DO J=1,N-1 >*/
    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j) {
/*<          J1=J+1 >*/
	j1 = j + 1;
/*<          E2R(J1,J)=-FJAC(J1,J) >*/
	e2r[j1 + j * e2r_dim1] = -fjac[j1 + j * fjac_dim1];
/*<          E2I(J1,J)=0.D0 >*/
	e2i[j1 + j * e2i_dim1] = 0.;
/*<       END DO >*/
    }
/*<       DO J=1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          DO I=1,J >*/
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             E2I(I,J)=0.D0 >*/
	    e2i[i__ + j * e2i_dim1] = 0.;
/*<             E2R(I,J)=-FJAC(I,J) >*/
	    e2r[i__ + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
/*<          END DO >*/
	}
/*<          E2R(J,J)=E2R(J,J)+ALPHN >*/
	e2r[j + j * e2r_dim1] += *alphn;
/*<          E2I(J,J)=BETAN >*/
	e2i[j + j * e2i_dim1] = *betan;
/*<       END DO >*/
    }
/*<       CALL DECHC(N,LDE1,E2R,E2I,1,IP2,IER) >*/
    dechc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &c__1, &ip2[1], ier);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   55  CONTINUE >*/
L55:
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* decomc_ */


/*     END OF SUBROUTINE DECOMC */

/* *********************************************************** */

/*<        >*/
/* Subroutine */ int slvrar_(integer *n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, integer *m1, integer *m2, integer *
	nm1, doublereal *fac1, doublereal *e1, integer *lde1, doublereal *z1, 
	doublereal *f1, integer *ip1, integer *iphes, integer *ier, integer *
	ijob)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset, i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, j, k;
    static doublereal s1;
    static integer mm, mp, im1, mp1, jkm;
    extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal sum1;
    extern /* Subroutine */ int solb_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), solh_(integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *);
    static doublereal zsafe;

/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<        >*/
/*<       COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG >*/

/*<       GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,13,15), IJOB >*/
    /* Parameter adjustments */
    --iphes;
    --f1;
    --z1;
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip1;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    e1_dim1 = *lde1;
    e1_offset = 1 + e1_dim1;
    e1 -= e1_offset;

    /* Function Body */
    switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L13;
	case 15:  goto L15;
    }

/* ----------------------------------------------------------- */

/*<    1  CONTINUE >*/
L1:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Z1(I)=Z1(I)-F1(I)*FAC1 >*/
	z1[i__] -= f1[i__] * *fac1;
/*<       END DO >*/
    }
/*<       CALL SOL (N,LDE1,E1,Z1,IP1) >*/
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   11  CONTINUE >*/
L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Z1(I)=Z1(I)-F1(I)*FAC1 >*/
	z1[i__] -= f1[i__] * *fac1;
/*<       END DO >*/
    }
/*<  48   CONTINUE >*/
L48:
/*<       MM=M1/M2 >*/
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM1=0.D0 >*/
	sum1 = 0.;
/*<          DO K=MM-1,0,-1 >*/
	for (k = mm - 1; k >= 0; --k) {
/*<             JKM=J+K*M2 >*/
	    jkm = j + k * *m2;
/*<             SUM1=(Z1(JKM)+SUM1)/FAC1 >*/
	    sum1 = (z1[jkm] + sum1) / *fac1;
/*<             DO I=1,NM1 >*/
	    i__2 = *nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<                IM1=I+M1 >*/
		im1 = i__ + *m1;
/*<                Z1(IM1)=Z1(IM1)+FJAC(I,JKM)*SUM1 >*/
		z1[im1] += fjac[i__ + jkm * fjac_dim1] * sum1;
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOL (NM1,LDE1,E1,Z1(M1+1),IP1) >*/
    sol_(nm1, lde1, &e1[e1_offset], &z1[*m1 + 1], &ip1[1]);
/*<  49   CONTINUE >*/
L49:
/*<       DO I=M1,1,-1 >*/
    for (i__ = *m1; i__ >= 1; --i__) {
/*<          Z1(I)=(Z1(I)+Z1(M2+I))/FAC1 >*/
	z1[i__] = (z1[i__] + z1[*m2 + i__]) / *fac1;
/*<       END DO >*/
    }
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    2  CONTINUE >*/
L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Z1(I)=Z1(I)-F1(I)*FAC1 >*/
	z1[i__] -= f1[i__] * *fac1;
/*<       END DO >*/
    }
/*<       CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1) >*/
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]
	    );
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   12  CONTINUE >*/
L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Z1(I)=Z1(I)-F1(I)*FAC1 >*/
	z1[i__] -= f1[i__] * *fac1;
/*<       END DO >*/
    }
/*<   45  CONTINUE >*/
L45:
/*<       MM=M1/M2 >*/
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM1=0.D0 >*/
	sum1 = 0.;
/*<          DO K=MM-1,0,-1 >*/
	for (k = mm - 1; k >= 0; --k) {
/*<             JKM=J+K*M2 >*/
	    jkm = j + k * *m2;
/*<             SUM1=(Z1(JKM)+SUM1)/FAC1 >*/
	    sum1 = (z1[jkm] + sum1) / *fac1;
/*<             DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC) >*/
/* Computing MAX */
	    i__2 = 1, i__3 = j - *mujac;
/* Computing MIN */
	    i__5 = *nm1, i__6 = j + *mljac;
	    i__4 = min(i__5,i__6);
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
/*<                IM1=I+M1 >*/
		im1 = i__ + *m1;
/*<                Z1(IM1)=Z1(IM1)+FJAC(I+MUJAC+1-J,JKM)*SUM1 >*/
		z1[im1] += fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1] * 
			sum1;
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOLB (NM1,LDE1,E1,MLE,MUE,Z1(M1+1),IP1) >*/
    solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[*m1 + 1],
	     &ip1[1]);
/*<       GOTO 49 >*/
    goto L49;

/* ----------------------------------------------------------- */

/*<    3  CONTINUE >*/
L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S1=0.0D0 >*/
	s1 = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS) >*/
/* Computing MAX */
	i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__3 = min(i__5,i__6);
	for (j = max(i__4,i__2); j <= i__3; ++j) {
/*<             S1=S1-FMAS(I-J+MBDIAG,J)*F1(J) >*/
	    s1 -= fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
/*<          END DO >*/
	}
/*<          Z1(I)=Z1(I)+S1*FAC1 >*/
	z1[i__] += s1 * *fac1;
/*<       END DO >*/
    }
/*<       CALL SOL (N,LDE1,E1,Z1,IP1) >*/
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   13  CONTINUE >*/
L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO I=1,M1 >*/
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Z1(I)=Z1(I)-F1(I)*FAC1 >*/
	z1[i__] -= f1[i__] * *fac1;
/*<       END DO >*/
    }
/*<       DO I=1,NM1 >*/
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IM1=I+M1 >*/
	im1 = i__ + *m1;
/*<          S1=0.0D0 >*/
	s1 = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS) >*/
/* Computing MAX */
	i__3 = 1, i__4 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *nm1, i__6 = i__ + *mumas;
	i__2 = min(i__5,i__6);
	for (j = max(i__3,i__4); j <= i__2; ++j) {
/*<             S1=S1-FMAS(I-J+MBDIAG,J)*F1(J+M1) >*/
	    s1 -= fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j + *m1]
		    ;
/*<          END DO >*/
	}
/*<          Z1(IM1)=Z1(IM1)+S1*FAC1 >*/
	z1[im1] += s1 * *fac1;
/*<       END DO >*/
    }
/*<       IF (IJOB.EQ.14) GOTO 45 >*/
    if (*ijob == 14) {
	goto L45;
    }
/*<       GOTO 48 >*/
    goto L48;

/* ----------------------------------------------------------- */

/*<    4  CONTINUE >*/
L4:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S1=0.0D0 >*/
	s1 = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS) >*/
/* Computing MAX */
	i__2 = 1, i__3 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__4 = min(i__5,i__6);
	for (j = max(i__2,i__3); j <= i__4; ++j) {
/*<             S1=S1-FMAS(I-J+MBDIAG,J)*F1(J) >*/
	    s1 -= fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
/*<          END DO >*/
	}
/*<          Z1(I)=Z1(I)+S1*FAC1 >*/
	z1[i__] += s1 * *fac1;
/*<       END DO >*/
    }
/*<       CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1) >*/
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]
	    );
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    5  CONTINUE >*/
L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S1=0.0D0 >*/
	s1 = 0.;
/*<          DO J=1,N >*/
	i__4 = *n;
	for (j = 1; j <= i__4; ++j) {
/*<             S1=S1-FMAS(I,J)*F1(J) >*/
	    s1 -= fmas[i__ + j * fmas_dim1] * f1[j];
/*<          END DO >*/
	}
/*<          Z1(I)=Z1(I)+S1*FAC1 >*/
	z1[i__] += s1 * *fac1;
/*<       END DO >*/
    }
/*<       CALL SOL (N,LDE1,E1,Z1,IP1) >*/
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   15  CONTINUE >*/
L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO I=1,M1 >*/
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Z1(I)=Z1(I)-F1(I)*FAC1 >*/
	z1[i__] -= f1[i__] * *fac1;
/*<       END DO >*/
    }
/*<       DO I=1,NM1 >*/
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IM1=I+M1 >*/
	im1 = i__ + *m1;
/*<          S1=0.0D0 >*/
	s1 = 0.;
/*<          DO J=1,NM1 >*/
	i__4 = *nm1;
	for (j = 1; j <= i__4; ++j) {
/*<             S1=S1-FMAS(I,J)*F1(J+M1) >*/
	    s1 -= fmas[i__ + j * fmas_dim1] * f1[j + *m1];
/*<          END DO >*/
	}
/*<          Z1(IM1)=Z1(IM1)+S1*FAC1 >*/
	z1[im1] += s1 * *fac1;
/*<       END DO >*/
    }
/*<       GOTO 48 >*/
    goto L48;

/* ----------------------------------------------------------- */

/*<    6  CONTINUE >*/
L6:
/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
/* ---  THIS OPTION IS NOT PROVIDED */
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    7  CONTINUE >*/
L7:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Z1(I)=Z1(I)-F1(I)*FAC1 >*/
	z1[i__] -= f1[i__] * *fac1;
/*<       END DO >*/
    }
/*<       DO MM=N-2,1,-1 >*/
    for (mm = *n - 2; mm >= 1; --mm) {
/*<           MP=N-MM >*/
	mp = *n - mm;
/*<           MP1=MP-1 >*/
	mp1 = mp - 1;
/*<           I=IPHES(MP) >*/
	i__ = iphes[mp];
/*<           IF (I.EQ.MP) GOTO 746 >*/
	if (i__ == mp) {
	    goto L746;
	}
/*<           ZSAFE=Z1(MP) >*/
	zsafe = z1[mp];
/*<           Z1(MP)=Z1(I) >*/
	z1[mp] = z1[i__];
/*<           Z1(I)=ZSAFE >*/
	z1[i__] = zsafe;
/*<  746      CONTINUE >*/
L746:
/*<           DO I=MP+1,N  >*/
	i__1 = *n;
	for (i__ = mp + 1; i__ <= i__1; ++i__) {
/*<              Z1(I)=Z1(I)-FJAC(I,MP1)*Z1(MP) >*/
	    z1[i__] -= fjac[i__ + mp1 * fjac_dim1] * z1[mp];
/*<           END DO >*/
	}
/*<        END DO >*/
    }
/*<        CALL SOLH(N,LDE1,E1,1,Z1,IP1) >*/
    solh_(n, lde1, &e1[e1_offset], &c__1, &z1[1], &ip1[1]);
/*<        DO MM=1,N-2 >*/
    i__1 = *n - 2;
    for (mm = 1; mm <= i__1; ++mm) {
/*<           MP=N-MM >*/
	mp = *n - mm;
/*<           MP1=MP-1 >*/
	mp1 = mp - 1;
/*<           DO I=MP+1,N  >*/
	i__4 = *n;
	for (i__ = mp + 1; i__ <= i__4; ++i__) {
/*<              Z1(I)=Z1(I)+FJAC(I,MP1)*Z1(MP) >*/
	    z1[i__] += fjac[i__ + mp1 * fjac_dim1] * z1[mp];
/*<           END DO >*/
	}
/*<           I=IPHES(MP) >*/
	i__ = iphes[mp];
/*<           IF (I.EQ.MP) GOTO 750 >*/
	if (i__ == mp) {
	    goto L750;
	}
/*<           ZSAFE=Z1(MP) >*/
	zsafe = z1[mp];
/*<           Z1(MP)=Z1(I) >*/
	z1[mp] = z1[i__];
/*<           Z1(I)=ZSAFE >*/
	z1[i__] = zsafe;
/*<  750      CONTINUE >*/
L750:
/*<       END DO >*/
	;
    }
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   55  CONTINUE >*/
L55:
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* slvrar_ */


/*     END OF SUBROUTINE SLVRAR */

/* *********************************************************** */

/*<        >*/
/* Subroutine */ int slvrai_(integer *n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, integer *m1, integer *m2, integer *
	nm1, doublereal *alphn, doublereal *betan, doublereal *e2r, 
	doublereal *e2i, integer *lde1, doublereal *z2, doublereal *z3, 
	doublereal *f2, doublereal *f3, doublereal *cont, integer *ip2, 
	integer *iphes, integer *ier, integer *ijob)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e2r_dim1, 
	    e2r_offset, e2i_dim1, e2i_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal s2, s3, bb;
    static integer mm, mp, im1, jm1, mp1;
    static doublereal z2i, z3i;
    static integer jkm, mpi;
    static doublereal sum2, sum3, abno;
    extern /* Subroutine */ int solc_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static integer iimu;
    static doublereal sumh, e1imp;
    extern /* Subroutine */ int solbc_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal zsafe;
    extern /* Subroutine */ int solhc_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);

/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<        >*/
/*<       DIMENSION E2R(LDE1,NM1),E2I(LDE1,NM1) >*/
/*<       COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG >*/

/*<       GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,13,15), IJOB >*/
    /* Parameter adjustments */
    --iphes;
    --f3;
    --f2;
    --z3;
    --z2;
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip2;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    e2i_dim1 = *lde1;
    e2i_offset = 1 + e2i_dim1;
    e2i -= e2i_offset;
    e2r_dim1 = *lde1;
    e2r_offset = 1 + e2r_dim1;
    e2r -= e2r_offset;

    /* Function Body */
    switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L13;
	case 15:  goto L15;
    }

/* ----------------------------------------------------------- */

/*<    1  CONTINUE >*/
L1:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=-F2(I) >*/
	s2 = -f2[i__];
/*<          S3=-F3(I) >*/
	s3 = -f3[i__];
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       CALL SOLC (N,LDE1,E2R,E2I,Z2,Z3,IP2) >*/
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	    );
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   11  CONTINUE >*/
L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=-F2(I) >*/
	s2 = -f2[i__];
/*<          S3=-F3(I) >*/
	s3 = -f3[i__];
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<  48   ABNO=ALPHN**2+BETAN**2 >*/
L48:
/* Computing 2nd power */
    d__1 = *alphn;
/* Computing 2nd power */
    d__2 = *betan;
    abno = d__1 * d__1 + d__2 * d__2;
/*<       MM=M1/M2 >*/
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM2=0.D0 >*/
	sum2 = 0.;
/*<          SUM3=0.D0 >*/
	sum3 = 0.;
/*<          DO K=MM-1,0,-1 >*/
	for (k = mm - 1; k >= 0; --k) {
/*<             JKM=J+K*M2 >*/
	    jkm = j + k * *m2;
/*<             SUMH=(Z2(JKM)+SUM2)/ABNO >*/
	    sumh = (z2[jkm] + sum2) / abno;
/*<             SUM3=(Z3(JKM)+SUM3)/ABNO >*/
	    sum3 = (z3[jkm] + sum3) / abno;
/*<             SUM2=SUMH*ALPHN+SUM3*BETAN >*/
	    sum2 = sumh * *alphn + sum3 * *betan;
/*<             SUM3=SUM3*ALPHN-SUMH*BETAN >*/
	    sum3 = sum3 * *alphn - sumh * *betan;
/*<             DO I=1,NM1 >*/
	    i__2 = *nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<                IM1=I+M1 >*/
		im1 = i__ + *m1;
/*<                Z2(IM1)=Z2(IM1)+FJAC(I,JKM)*SUM2 >*/
		z2[im1] += fjac[i__ + jkm * fjac_dim1] * sum2;
/*<                Z3(IM1)=Z3(IM1)+FJAC(I,JKM)*SUM3 >*/
		z3[im1] += fjac[i__ + jkm * fjac_dim1] * sum3;
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOLC (NM1,LDE1,E2R,E2I,Z2(M1+1),Z3(M1+1),IP2) >*/
    solc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[*m1 + 1], &z3[*
	    m1 + 1], &ip2[1]);
/*<  49   CONTINUE >*/
L49:
/*<       DO I=M1,1,-1 >*/
    for (i__ = *m1; i__ >= 1; --i__) {
/*<          MPI=M2+I >*/
	mpi = *m2 + i__;
/*<          Z2I=Z2(I)+Z2(MPI) >*/
	z2i = z2[i__] + z2[mpi];
/*<          Z3I=Z3(I)+Z3(MPI) >*/
	z3i = z3[i__] + z3[mpi];
/*<          Z3(I)=(Z3I*ALPHN-Z2I*BETAN)/ABNO >*/
	z3[i__] = (z3i * *alphn - z2i * *betan) / abno;
/*<          Z2(I)=(Z2I*ALPHN+Z3I*BETAN)/ABNO >*/
	z2[i__] = (z2i * *alphn + z3i * *betan) / abno;
/*<       END DO >*/
    }
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    2  CONTINUE >*/
L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=-F2(I) >*/
	s2 = -f2[i__];
/*<          S3=-F3(I) >*/
	s3 = -f3[i__];
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       CALL SOLBC (N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2) >*/
    solbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &z2[1], &z3[1], &ip2[1]);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   12  CONTINUE >*/
L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=-F2(I) >*/
	s2 = -f2[i__];
/*<          S3=-F3(I) >*/
	s3 = -f3[i__];
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<   45  ABNO=ALPHN**2+BETAN**2 >*/
L45:
/* Computing 2nd power */
    d__1 = *alphn;
/* Computing 2nd power */
    d__2 = *betan;
    abno = d__1 * d__1 + d__2 * d__2;
/*<       MM=M1/M2 >*/
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM2=0.D0 >*/
	sum2 = 0.;
/*<          SUM3=0.D0 >*/
	sum3 = 0.;
/*<          DO K=MM-1,0,-1 >*/
	for (k = mm - 1; k >= 0; --k) {
/*<             JKM=J+K*M2 >*/
	    jkm = j + k * *m2;
/*<             SUMH=(Z2(JKM)+SUM2)/ABNO >*/
	    sumh = (z2[jkm] + sum2) / abno;
/*<             SUM3=(Z3(JKM)+SUM3)/ABNO >*/
	    sum3 = (z3[jkm] + sum3) / abno;
/*<             SUM2=SUMH*ALPHN+SUM3*BETAN >*/
	    sum2 = sumh * *alphn + sum3 * *betan;
/*<             SUM3=SUM3*ALPHN-SUMH*BETAN >*/
	    sum3 = sum3 * *alphn - sumh * *betan;
/*<             DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC) >*/
/* Computing MAX */
	    i__2 = 1, i__3 = j - *mujac;
/* Computing MIN */
	    i__5 = *nm1, i__6 = j + *mljac;
	    i__4 = min(i__5,i__6);
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
/*<                IM1=I+M1 >*/
		im1 = i__ + *m1;
/*<                IIMU=I+MUJAC+1-J >*/
		iimu = i__ + *mujac + 1 - j;
/*<                Z2(IM1)=Z2(IM1)+FJAC(IIMU,JKM)*SUM2 >*/
		z2[im1] += fjac[iimu + jkm * fjac_dim1] * sum2;
/*<                Z3(IM1)=Z3(IM1)+FJAC(IIMU,JKM)*SUM3 >*/
		z3[im1] += fjac[iimu + jkm * fjac_dim1] * sum3;
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOLBC (NM1,LDE1,E2R,E2I,MLE,MUE,Z2(M1+1),Z3(M1+1),IP2) >*/
    solbc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &z2[*m1 + 1], &z3[*m1 + 1], &ip2[1]);
/*<       GOTO 49 >*/
    goto L49;

/* ----------------------------------------------------------- */

/*<    3  CONTINUE >*/
L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=0.0D0 >*/
	s2 = 0.;
/*<          S3=0.0D0 >*/
	s3 = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS) >*/
/* Computing MAX */
	i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__3 = min(i__5,i__6);
	for (j = max(i__4,i__2); j <= i__3; ++j) {
/*<             BB=FMAS(I-J+MBDIAG,J) >*/
	    bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
/*<             S2=S2-BB*F2(J) >*/
	    s2 -= bb * f2[j];
/*<             S3=S3-BB*F3(J) >*/
	    s3 -= bb * f3[j];
/*<          END DO >*/
	}
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2) >*/
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	    );
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   13  CONTINUE >*/
L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO I=1,M1 >*/
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=-F2(I) >*/
	s2 = -f2[i__];
/*<          S3=-F3(I) >*/
	s3 = -f3[i__];
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       DO I=1,NM1 >*/
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IM1=I+M1 >*/
	im1 = i__ + *m1;
/*<          S2=0.0D0 >*/
	s2 = 0.;
/*<          S3=0.0D0 >*/
	s3 = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS) >*/
/* Computing MAX */
	i__3 = 1, i__4 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *nm1, i__6 = i__ + *mumas;
	i__2 = min(i__5,i__6);
	for (j = max(i__3,i__4); j <= i__2; ++j) {
/*<             JM1=J+M1 >*/
	    jm1 = j + *m1;
/*<             BB=FMAS(I-J+MBDIAG,J) >*/
	    bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
/*<             S2=S2-BB*F2(JM1) >*/
	    s2 -= bb * f2[jm1];
/*<             S3=S3-BB*F3(JM1) >*/
	    s3 -= bb * f3[jm1];
/*<          END DO >*/
	}
/*<          Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN >*/
	z2[im1] = z2[im1] + s2 * *alphn - s3 * *betan;
/*<          Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN >*/
	z3[im1] = z3[im1] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       IF (IJOB.EQ.14) GOTO 45 >*/
    if (*ijob == 14) {
	goto L45;
    }
/*<       GOTO 48 >*/
    goto L48;

/* ----------------------------------------------------------- */

/*<    4  CONTINUE >*/
L4:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=0.0D0 >*/
	s2 = 0.;
/*<          S3=0.0D0 >*/
	s3 = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS) >*/
/* Computing MAX */
	i__2 = 1, i__3 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__4 = min(i__5,i__6);
	for (j = max(i__2,i__3); j <= i__4; ++j) {
/*<             BB=FMAS(I-J+MBDIAG,J) >*/
	    bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
/*<             S2=S2-BB*F2(J) >*/
	    s2 -= bb * f2[j];
/*<             S3=S3-BB*F3(J) >*/
	    s3 -= bb * f3[j];
/*<          END DO >*/
	}
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       CALL SOLBC(N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2) >*/
    solbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &z2[1], &z3[1], &ip2[1]);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    5  CONTINUE >*/
L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=0.0D0 >*/
	s2 = 0.;
/*<          S3=0.0D0 >*/
	s3 = 0.;
/*<          DO J=1,N >*/
	i__4 = *n;
	for (j = 1; j <= i__4; ++j) {
/*<             BB=FMAS(I,J) >*/
	    bb = fmas[i__ + j * fmas_dim1];
/*<             S2=S2-BB*F2(J) >*/
	    s2 -= bb * f2[j];
/*<             S3=S3-BB*F3(J) >*/
	    s3 -= bb * f3[j];
/*<          END DO >*/
	}
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2) >*/
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	    );
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   15  CONTINUE >*/
L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO I=1,M1 >*/
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=-F2(I) >*/
	s2 = -f2[i__];
/*<          S3=-F3(I) >*/
	s3 = -f3[i__];
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       DO I=1,NM1 >*/
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IM1=I+M1 >*/
	im1 = i__ + *m1;
/*<          S2=0.0D0 >*/
	s2 = 0.;
/*<          S3=0.0D0 >*/
	s3 = 0.;
/*<          DO J=1,NM1 >*/
	i__4 = *nm1;
	for (j = 1; j <= i__4; ++j) {
/*<             JM1=J+M1 >*/
	    jm1 = j + *m1;
/*<             BB=FMAS(I,J) >*/
	    bb = fmas[i__ + j * fmas_dim1];
/*<             S2=S2-BB*F2(JM1) >*/
	    s2 -= bb * f2[jm1];
/*<             S3=S3-BB*F3(JM1) >*/
	    s3 -= bb * f3[jm1];
/*<          END DO >*/
	}
/*<          Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN >*/
	z2[im1] = z2[im1] + s2 * *alphn - s3 * *betan;
/*<          Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN >*/
	z3[im1] = z3[im1] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       GOTO 48 >*/
    goto L48;

/* ----------------------------------------------------------- */

/*<    6  CONTINUE >*/
L6:
/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
/* ---  THIS OPTION IS NOT PROVIDED */
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    7  CONTINUE >*/
L7:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=-F2(I) >*/
	s2 = -f2[i__];
/*<          S3=-F3(I) >*/
	s3 = -f3[i__];
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       DO MM=N-2,1,-1 >*/
    for (mm = *n - 2; mm >= 1; --mm) {
/*<           MP=N-MM >*/
	mp = *n - mm;
/*<           MP1=MP-1 >*/
	mp1 = mp - 1;
/*<           I=IPHES(MP) >*/
	i__ = iphes[mp];
/*<           IF (I.EQ.MP) GOTO 746 >*/
	if (i__ == mp) {
	    goto L746;
	}
/*<           ZSAFE=Z2(MP) >*/
	zsafe = z2[mp];
/*<           Z2(MP)=Z2(I) >*/
	z2[mp] = z2[i__];
/*<           Z2(I)=ZSAFE  >*/
	z2[i__] = zsafe;
/*<           ZSAFE=Z3(MP) >*/
	zsafe = z3[mp];
/*<           Z3(MP)=Z3(I) >*/
	z3[mp] = z3[i__];
/*<           Z3(I)=ZSAFE >*/
	z3[i__] = zsafe;
/*<  746      CONTINUE >*/
L746:
/*<           DO I=MP+1,N  >*/
	i__1 = *n;
	for (i__ = mp + 1; i__ <= i__1; ++i__) {
/*<              E1IMP=FJAC(I,MP1) >*/
	    e1imp = fjac[i__ + mp1 * fjac_dim1];
/*<              Z2(I)=Z2(I)-E1IMP*Z2(MP) >*/
	    z2[i__] -= e1imp * z2[mp];
/*<              Z3(I)=Z3(I)-E1IMP*Z3(MP) >*/
	    z3[i__] -= e1imp * z3[mp];
/*<           END DO >*/
	}
/*<        END DO >*/
    }
/*<        CALL SOLHC(N,LDE1,E2R,E2I,1,Z2,Z3,IP2) >*/
    solhc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &c__1, &z2[1], &z3[1],
	     &ip2[1]);
/*<        DO MM=1,N-2 >*/
    i__1 = *n - 2;
    for (mm = 1; mm <= i__1; ++mm) {
/*<           MP=N-MM >*/
	mp = *n - mm;
/*<           MP1=MP-1 >*/
	mp1 = mp - 1;
/*<           DO I=MP+1,N  >*/
	i__4 = *n;
	for (i__ = mp + 1; i__ <= i__4; ++i__) {
/*<              E1IMP=FJAC(I,MP1) >*/
	    e1imp = fjac[i__ + mp1 * fjac_dim1];
/*<              Z2(I)=Z2(I)+E1IMP*Z2(MP) >*/
	    z2[i__] += e1imp * z2[mp];
/*<              Z3(I)=Z3(I)+E1IMP*Z3(MP) >*/
	    z3[i__] += e1imp * z3[mp];
/*<           END DO >*/
	}
/*<           I=IPHES(MP) >*/
	i__ = iphes[mp];
/*<           IF (I.EQ.MP) GOTO 750 >*/
	if (i__ == mp) {
	    goto L750;
	}
/*<           ZSAFE=Z2(MP) >*/
	zsafe = z2[mp];
/*<           Z2(MP)=Z2(I) >*/
	z2[mp] = z2[i__];
/*<           Z2(I)=ZSAFE  >*/
	z2[i__] = zsafe;
/*<           ZSAFE=Z3(MP) >*/
	zsafe = z3[mp];
/*<           Z3(MP)=Z3(I) >*/
	z3[mp] = z3[i__];
/*<           Z3(I)=ZSAFE >*/
	z3[i__] = zsafe;
/*<  750      CONTINUE >*/
L750:
/*<       END DO >*/
	;
    }
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   55  CONTINUE >*/
L55:
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* slvrai_ */


/*     END OF SUBROUTINE SLVRAI */

/* *********************************************************** */

/*<        >*/
/* Subroutine */ int slvrad_(integer *n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, integer *m1, integer *m2, integer *
	nm1, doublereal *fac1, doublereal *alphn, doublereal *betan, 
	doublereal *e1, doublereal *e2r, doublereal *e2i, integer *lde1, 
	doublereal *z1, doublereal *z2, doublereal *z3, doublereal *f1, 
	doublereal *f2, doublereal *f3, doublereal *cont, integer *ip1, 
	integer *ip2, integer *iphes, integer *ier, integer *ijob)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset, e2r_dim1, e2r_offset, e2i_dim1, e2i_offset, i__1, i__2,
	     i__3, i__4, i__5, i__6;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal s1, s2, s3, bb;
    static integer mm, mp, j1b, j2b, im1, jm1, mp1;
    static doublereal z2i, z3i;
    static integer jkm, mpi;
    extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal sum1, sum2, sum3, ffja, abno;
    extern /* Subroutine */ int solb_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), solc_(integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *), solh_(integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal sumh, e1imp;
    extern /* Subroutine */ int solbc_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal zsafe;
    extern /* Subroutine */ int solhc_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);

/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<        >*/
/*<       COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG >*/

/*<       GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,13,15), IJOB >*/
    /* Parameter adjustments */
    --iphes;
    --f3;
    --f2;
    --f1;
    --z3;
    --z2;
    --z1;
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip2;
    --ip1;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    e2i_dim1 = *lde1;
    e2i_offset = 1 + e2i_dim1;
    e2i -= e2i_offset;
    e2r_dim1 = *lde1;
    e2r_offset = 1 + e2r_dim1;
    e2r -= e2r_offset;
    e1_dim1 = *lde1;
    e1_offset = 1 + e1_dim1;
    e1 -= e1_offset;

    /* Function Body */
    switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L13;
	case 15:  goto L15;
    }

/* ----------------------------------------------------------- */

/*<    1  CONTINUE >*/
L1:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=-F2(I) >*/
	s2 = -f2[i__];
/*<          S3=-F3(I) >*/
	s3 = -f3[i__];
/*<          Z1(I)=Z1(I)-F1(I)*FAC1 >*/
	z1[i__] -= f1[i__] * *fac1;
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       CALL SOL (N,LDE1,E1,Z1,IP1) >*/
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
/*<       CALL SOLC (N,LDE1,E2R,E2I,Z2,Z3,IP2) >*/
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	    );
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   11  CONTINUE >*/
L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=-F2(I) >*/
	s2 = -f2[i__];
/*<          S3=-F3(I) >*/
	s3 = -f3[i__];
/*<          Z1(I)=Z1(I)-F1(I)*FAC1 >*/
	z1[i__] -= f1[i__] * *fac1;
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<  48   ABNO=ALPHN**2+BETAN**2 >*/
L48:
/* Computing 2nd power */
    d__1 = *alphn;
/* Computing 2nd power */
    d__2 = *betan;
    abno = d__1 * d__1 + d__2 * d__2;
/*<       MM=M1/M2 >*/
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM1=0.D0 >*/
	sum1 = 0.;
/*<          SUM2=0.D0 >*/
	sum2 = 0.;
/*<          SUM3=0.D0 >*/
	sum3 = 0.;
/*<          DO K=MM-1,0,-1 >*/
	for (k = mm - 1; k >= 0; --k) {
/*<             JKM=J+K*M2 >*/
	    jkm = j + k * *m2;
/*<             SUM1=(Z1(JKM)+SUM1)/FAC1 >*/
	    sum1 = (z1[jkm] + sum1) / *fac1;
/*<             SUMH=(Z2(JKM)+SUM2)/ABNO >*/
	    sumh = (z2[jkm] + sum2) / abno;
/*<             SUM3=(Z3(JKM)+SUM3)/ABNO >*/
	    sum3 = (z3[jkm] + sum3) / abno;
/*<             SUM2=SUMH*ALPHN+SUM3*BETAN >*/
	    sum2 = sumh * *alphn + sum3 * *betan;
/*<             SUM3=SUM3*ALPHN-SUMH*BETAN >*/
	    sum3 = sum3 * *alphn - sumh * *betan;
/*<             DO I=1,NM1 >*/
	    i__2 = *nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<                IM1=I+M1 >*/
		im1 = i__ + *m1;
/*<                Z1(IM1)=Z1(IM1)+FJAC(I,JKM)*SUM1 >*/
		z1[im1] += fjac[i__ + jkm * fjac_dim1] * sum1;
/*<                Z2(IM1)=Z2(IM1)+FJAC(I,JKM)*SUM2 >*/
		z2[im1] += fjac[i__ + jkm * fjac_dim1] * sum2;
/*<                Z3(IM1)=Z3(IM1)+FJAC(I,JKM)*SUM3 >*/
		z3[im1] += fjac[i__ + jkm * fjac_dim1] * sum3;
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOL (NM1,LDE1,E1,Z1(M1+1),IP1) >*/
    sol_(nm1, lde1, &e1[e1_offset], &z1[*m1 + 1], &ip1[1]);
/*<       CALL SOLC (NM1,LDE1,E2R,E2I,Z2(M1+1),Z3(M1+1),IP2) >*/
    solc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[*m1 + 1], &z3[*
	    m1 + 1], &ip2[1]);
/*<  49   CONTINUE >*/
L49:
/*<       DO I=M1,1,-1 >*/
    for (i__ = *m1; i__ >= 1; --i__) {
/*<          MPI=M2+I >*/
	mpi = *m2 + i__;
/*<          Z1(I)=(Z1(I)+Z1(MPI))/FAC1 >*/
	z1[i__] = (z1[i__] + z1[mpi]) / *fac1;
/*<          Z2I=Z2(I)+Z2(MPI) >*/
	z2i = z2[i__] + z2[mpi];
/*<          Z3I=Z3(I)+Z3(MPI) >*/
	z3i = z3[i__] + z3[mpi];
/*<          Z3(I)=(Z3I*ALPHN-Z2I*BETAN)/ABNO >*/
	z3[i__] = (z3i * *alphn - z2i * *betan) / abno;
/*<          Z2(I)=(Z2I*ALPHN+Z3I*BETAN)/ABNO >*/
	z2[i__] = (z2i * *alphn + z3i * *betan) / abno;
/*<       END DO >*/
    }
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    2  CONTINUE >*/
L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=-F2(I) >*/
	s2 = -f2[i__];
/*<          S3=-F3(I) >*/
	s3 = -f3[i__];
/*<          Z1(I)=Z1(I)-F1(I)*FAC1 >*/
	z1[i__] -= f1[i__] * *fac1;
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1) >*/
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]
	    );
/*<       CALL SOLBC (N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2) >*/
    solbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &z2[1], &z3[1], &ip2[1]);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   12  CONTINUE >*/
L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=-F2(I) >*/
	s2 = -f2[i__];
/*<          S3=-F3(I) >*/
	s3 = -f3[i__];
/*<          Z1(I)=Z1(I)-F1(I)*FAC1 >*/
	z1[i__] -= f1[i__] * *fac1;
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<   45  ABNO=ALPHN**2+BETAN**2 >*/
L45:
/* Computing 2nd power */
    d__1 = *alphn;
/* Computing 2nd power */
    d__2 = *betan;
    abno = d__1 * d__1 + d__2 * d__2;
/*<       MM=M1/M2 >*/
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM1=0.D0 >*/
	sum1 = 0.;
/*<          SUM2=0.D0 >*/
	sum2 = 0.;
/*<          SUM3=0.D0 >*/
	sum3 = 0.;
/*<          DO K=MM-1,0,-1 >*/
	for (k = mm - 1; k >= 0; --k) {
/*<             JKM=J+K*M2 >*/
	    jkm = j + k * *m2;
/*<             SUM1=(Z1(JKM)+SUM1)/FAC1 >*/
	    sum1 = (z1[jkm] + sum1) / *fac1;
/*<             SUMH=(Z2(JKM)+SUM2)/ABNO >*/
	    sumh = (z2[jkm] + sum2) / abno;
/*<             SUM3=(Z3(JKM)+SUM3)/ABNO >*/
	    sum3 = (z3[jkm] + sum3) / abno;
/*<             SUM2=SUMH*ALPHN+SUM3*BETAN >*/
	    sum2 = sumh * *alphn + sum3 * *betan;
/*<             SUM3=SUM3*ALPHN-SUMH*BETAN >*/
	    sum3 = sum3 * *alphn - sumh * *betan;
/*<             DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC) >*/
/* Computing MAX */
	    i__2 = 1, i__3 = j - *mujac;
/* Computing MIN */
	    i__5 = *nm1, i__6 = j + *mljac;
	    i__4 = min(i__5,i__6);
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
/*<                IM1=I+M1 >*/
		im1 = i__ + *m1;
/*<                FFJA=FJAC(I+MUJAC+1-J,JKM) >*/
		ffja = fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1];
/*<                Z1(IM1)=Z1(IM1)+FFJA*SUM1 >*/
		z1[im1] += ffja * sum1;
/*<                Z2(IM1)=Z2(IM1)+FFJA*SUM2 >*/
		z2[im1] += ffja * sum2;
/*<                Z3(IM1)=Z3(IM1)+FFJA*SUM3 >*/
		z3[im1] += ffja * sum3;
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOLB (NM1,LDE1,E1,MLE,MUE,Z1(M1+1),IP1) >*/
    solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[*m1 + 1],
	     &ip1[1]);
/*<       CALL SOLBC (NM1,LDE1,E2R,E2I,MLE,MUE,Z2(M1+1),Z3(M1+1),IP2) >*/
    solbc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &z2[*m1 + 1], &z3[*m1 + 1], &ip2[1]);
/*<       GOTO 49 >*/
    goto L49;

/* ----------------------------------------------------------- */

/*<    3  CONTINUE >*/
L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S1=0.0D0 >*/
	s1 = 0.;
/*<          S2=0.0D0 >*/
	s2 = 0.;
/*<          S3=0.0D0 >*/
	s3 = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS) >*/
/* Computing MAX */
	i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__3 = min(i__5,i__6);
	for (j = max(i__4,i__2); j <= i__3; ++j) {
/*<             BB=FMAS(I-J+MBDIAG,J) >*/
	    bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
/*<             S1=S1-BB*F1(J) >*/
	    s1 -= bb * f1[j];
/*<             S2=S2-BB*F2(J) >*/
	    s2 -= bb * f2[j];
/*<             S3=S3-BB*F3(J) >*/
	    s3 -= bb * f3[j];
/*<          END DO >*/
	}
/*<          Z1(I)=Z1(I)+S1*FAC1 >*/
	z1[i__] += s1 * *fac1;
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       CALL SOL (N,LDE1,E1,Z1,IP1) >*/
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
/*<       CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2) >*/
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	    );
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   13  CONTINUE >*/
L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO I=1,M1 >*/
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=-F2(I) >*/
	s2 = -f2[i__];
/*<          S3=-F3(I) >*/
	s3 = -f3[i__];
/*<          Z1(I)=Z1(I)-F1(I)*FAC1 >*/
	z1[i__] -= f1[i__] * *fac1;
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       DO I=1,NM1 >*/
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IM1=I+M1 >*/
	im1 = i__ + *m1;
/*<          S1=0.0D0 >*/
	s1 = 0.;
/*<          S2=0.0D0 >*/
	s2 = 0.;
/*<          S3=0.0D0 >*/
	s3 = 0.;
/*<          J1B=MAX(1,I-MLMAS) >*/
/* Computing MAX */
	i__3 = 1, i__4 = i__ - *mlmas;
	j1b = max(i__3,i__4);
/*<          J2B=MIN(NM1,I+MUMAS) >*/
/* Computing MIN */
	i__3 = *nm1, i__4 = i__ + *mumas;
	j2b = min(i__3,i__4);
/*<          DO J=J1B,J2B >*/
	i__3 = j2b;
	for (j = j1b; j <= i__3; ++j) {
/*<             JM1=J+M1 >*/
	    jm1 = j + *m1;
/*<             BB=FMAS(I-J+MBDIAG,J) >*/
	    bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
/*<             S1=S1-BB*F1(JM1) >*/
	    s1 -= bb * f1[jm1];
/*<             S2=S2-BB*F2(JM1) >*/
	    s2 -= bb * f2[jm1];
/*<             S3=S3-BB*F3(JM1) >*/
	    s3 -= bb * f3[jm1];
/*<          END DO >*/
	}
/*<          Z1(IM1)=Z1(IM1)+S1*FAC1 >*/
	z1[im1] += s1 * *fac1;
/*<          Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN >*/
	z2[im1] = z2[im1] + s2 * *alphn - s3 * *betan;
/*<          Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN >*/
	z3[im1] = z3[im1] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       IF (IJOB.EQ.14) GOTO 45 >*/
    if (*ijob == 14) {
	goto L45;
    }
/*<       GOTO 48 >*/
    goto L48;

/* ----------------------------------------------------------- */

/*<    4  CONTINUE >*/
L4:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S1=0.0D0 >*/
	s1 = 0.;
/*<          S2=0.0D0 >*/
	s2 = 0.;
/*<          S3=0.0D0 >*/
	s3 = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS) >*/
/* Computing MAX */
	i__3 = 1, i__4 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__2 = min(i__5,i__6);
	for (j = max(i__3,i__4); j <= i__2; ++j) {
/*<             BB=FMAS(I-J+MBDIAG,J) >*/
	    bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
/*<             S1=S1-BB*F1(J) >*/
	    s1 -= bb * f1[j];
/*<             S2=S2-BB*F2(J) >*/
	    s2 -= bb * f2[j];
/*<             S3=S3-BB*F3(J) >*/
	    s3 -= bb * f3[j];
/*<          END DO >*/
	}
/*<          Z1(I)=Z1(I)+S1*FAC1 >*/
	z1[i__] += s1 * *fac1;
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1) >*/
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]
	    );
/*<       CALL SOLBC(N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2) >*/
    solbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &z2[1], &z3[1], &ip2[1]);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    5  CONTINUE >*/
L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S1=0.0D0 >*/
	s1 = 0.;
/*<          S2=0.0D0 >*/
	s2 = 0.;
/*<          S3=0.0D0 >*/
	s3 = 0.;
/*<          DO J=1,N >*/
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/*<             BB=FMAS(I,J) >*/
	    bb = fmas[i__ + j * fmas_dim1];
/*<             S1=S1-BB*F1(J) >*/
	    s1 -= bb * f1[j];
/*<             S2=S2-BB*F2(J) >*/
	    s2 -= bb * f2[j];
/*<             S3=S3-BB*F3(J) >*/
	    s3 -= bb * f3[j];
/*<          END DO >*/
	}
/*<          Z1(I)=Z1(I)+S1*FAC1 >*/
	z1[i__] += s1 * *fac1;
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       CALL SOL (N,LDE1,E1,Z1,IP1) >*/
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
/*<       CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2) >*/
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	    );
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   15  CONTINUE >*/
L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO I=1,M1 >*/
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=-F2(I) >*/
	s2 = -f2[i__];
/*<          S3=-F3(I) >*/
	s3 = -f3[i__];
/*<          Z1(I)=Z1(I)-F1(I)*FAC1 >*/
	z1[i__] -= f1[i__] * *fac1;
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       DO I=1,NM1 >*/
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IM1=I+M1 >*/
	im1 = i__ + *m1;
/*<          S1=0.0D0 >*/
	s1 = 0.;
/*<          S2=0.0D0 >*/
	s2 = 0.;
/*<          S3=0.0D0 >*/
	s3 = 0.;
/*<          DO J=1,NM1 >*/
	i__2 = *nm1;
	for (j = 1; j <= i__2; ++j) {
/*<             JM1=J+M1 >*/
	    jm1 = j + *m1;
/*<             BB=FMAS(I,J) >*/
	    bb = fmas[i__ + j * fmas_dim1];
/*<             S1=S1-BB*F1(JM1) >*/
	    s1 -= bb * f1[jm1];
/*<             S2=S2-BB*F2(JM1) >*/
	    s2 -= bb * f2[jm1];
/*<             S3=S3-BB*F3(JM1) >*/
	    s3 -= bb * f3[jm1];
/*<          END DO >*/
	}
/*<          Z1(IM1)=Z1(IM1)+S1*FAC1 >*/
	z1[im1] += s1 * *fac1;
/*<          Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN >*/
	z2[im1] = z2[im1] + s2 * *alphn - s3 * *betan;
/*<          Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN >*/
	z3[im1] = z3[im1] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       GOTO 48 >*/
    goto L48;

/* ----------------------------------------------------------- */

/*<    6  CONTINUE >*/
L6:
/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
/* ---  THIS OPTION IS NOT PROVIDED */
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    7  CONTINUE >*/
L7:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          S2=-F2(I) >*/
	s2 = -f2[i__];
/*<          S3=-F3(I) >*/
	s3 = -f3[i__];
/*<          Z1(I)=Z1(I)-F1(I)*FAC1 >*/
	z1[i__] -= f1[i__] * *fac1;
/*<          Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN >*/
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
/*<          Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN >*/
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
/*<       END DO >*/
    }
/*<       DO MM=N-2,1,-1 >*/
    for (mm = *n - 2; mm >= 1; --mm) {
/*<           MP=N-MM >*/
	mp = *n - mm;
/*<           MP1=MP-1 >*/
	mp1 = mp - 1;
/*<           I=IPHES(MP) >*/
	i__ = iphes[mp];
/*<           IF (I.EQ.MP) GOTO 746 >*/
	if (i__ == mp) {
	    goto L746;
	}
/*<           ZSAFE=Z1(MP) >*/
	zsafe = z1[mp];
/*<           Z1(MP)=Z1(I) >*/
	z1[mp] = z1[i__];
/*<           Z1(I)=ZSAFE >*/
	z1[i__] = zsafe;
/*<           ZSAFE=Z2(MP) >*/
	zsafe = z2[mp];
/*<           Z2(MP)=Z2(I) >*/
	z2[mp] = z2[i__];
/*<           Z2(I)=ZSAFE  >*/
	z2[i__] = zsafe;
/*<           ZSAFE=Z3(MP) >*/
	zsafe = z3[mp];
/*<           Z3(MP)=Z3(I) >*/
	z3[mp] = z3[i__];
/*<           Z3(I)=ZSAFE >*/
	z3[i__] = zsafe;
/*<  746      CONTINUE >*/
L746:
/*<           DO I=MP+1,N  >*/
	i__1 = *n;
	for (i__ = mp + 1; i__ <= i__1; ++i__) {
/*<              E1IMP=FJAC(I,MP1) >*/
	    e1imp = fjac[i__ + mp1 * fjac_dim1];
/*<              Z1(I)=Z1(I)-E1IMP*Z1(MP) >*/
	    z1[i__] -= e1imp * z1[mp];
/*<              Z2(I)=Z2(I)-E1IMP*Z2(MP) >*/
	    z2[i__] -= e1imp * z2[mp];
/*<              Z3(I)=Z3(I)-E1IMP*Z3(MP) >*/
	    z3[i__] -= e1imp * z3[mp];
/*<           END DO >*/
	}
/*<        END DO >*/
    }
/*<        CALL SOLH(N,LDE1,E1,1,Z1,IP1) >*/
    solh_(n, lde1, &e1[e1_offset], &c__1, &z1[1], &ip1[1]);
/*<        CALL SOLHC(N,LDE1,E2R,E2I,1,Z2,Z3,IP2) >*/
    solhc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &c__1, &z2[1], &z3[1],
	     &ip2[1]);
/*<        DO MM=1,N-2 >*/
    i__1 = *n - 2;
    for (mm = 1; mm <= i__1; ++mm) {
/*<           MP=N-MM >*/
	mp = *n - mm;
/*<           MP1=MP-1 >*/
	mp1 = mp - 1;
/*<           DO I=MP+1,N  >*/
	i__2 = *n;
	for (i__ = mp + 1; i__ <= i__2; ++i__) {
/*<              E1IMP=FJAC(I,MP1) >*/
	    e1imp = fjac[i__ + mp1 * fjac_dim1];
/*<              Z1(I)=Z1(I)+E1IMP*Z1(MP) >*/
	    z1[i__] += e1imp * z1[mp];
/*<              Z2(I)=Z2(I)+E1IMP*Z2(MP) >*/
	    z2[i__] += e1imp * z2[mp];
/*<              Z3(I)=Z3(I)+E1IMP*Z3(MP) >*/
	    z3[i__] += e1imp * z3[mp];
/*<           END DO >*/
	}
/*<           I=IPHES(MP) >*/
	i__ = iphes[mp];
/*<           IF (I.EQ.MP) GOTO 750 >*/
	if (i__ == mp) {
	    goto L750;
	}
/*<           ZSAFE=Z1(MP) >*/
	zsafe = z1[mp];
/*<           Z1(MP)=Z1(I) >*/
	z1[mp] = z1[i__];
/*<           Z1(I)=ZSAFE >*/
	z1[i__] = zsafe;
/*<           ZSAFE=Z2(MP) >*/
	zsafe = z2[mp];
/*<           Z2(MP)=Z2(I) >*/
	z2[mp] = z2[i__];
/*<           Z2(I)=ZSAFE  >*/
	z2[i__] = zsafe;
/*<           ZSAFE=Z3(MP) >*/
	zsafe = z3[mp];
/*<           Z3(MP)=Z3(I) >*/
	z3[mp] = z3[i__];
/*<           Z3(I)=ZSAFE >*/
	z3[i__] = zsafe;
/*<  750      CONTINUE >*/
L750:
/*<       END DO >*/
	;
    }
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   55  CONTINUE >*/
L55:
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* slvrad_ */


/*     END OF SUBROUTINE SLVRAD */

/* *********************************************************** */

/*<        >*/
/* Subroutine */ int estrad_(integer *n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, doublereal *h__, doublereal *dd1, 
	doublereal *dd2, doublereal *dd3, S_fp fcn, integer *nfcn, doublereal 
	*y0, doublereal *y, integer *ijob, doublereal *x, integer *m1, 
	integer *m2, integer *nm1, doublereal *e1, integer *lde1, doublereal *
	z1, doublereal *z2, doublereal *z3, doublereal *cont, doublereal *f1, 
	doublereal *f2, integer *ip1, integer *iphes, doublereal *scal, 
	doublereal *err, logical *first, logical *reject, doublereal *fac1, 
	doublereal *rpar, integer *ipar)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, mm, mp, im1;
    extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal sum, hee1, hee2, hee3, sum1;
    extern /* Subroutine */ int solb_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), solh_(integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *);
    static doublereal zsafe;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<        >*/
/*<       DIMENSION CONT(N),RPAR(1),IPAR(1) >*/
/*<       LOGICAL FIRST,REJECT >*/
/*<       COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG >*/
/*<       HEE1=DD1/H >*/
    /* Parameter adjustments */
    --scal;
    --iphes;
    --f2;
    --f1;
    --cont;
    --z3;
    --z2;
    --z1;
    --y;
    --y0;
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip1;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    e1_dim1 = *lde1;
    e1_offset = 1 + e1_dim1;
    e1 -= e1_offset;
    --rpar;
    --ipar;

    /* Function Body */
    hee1 = *dd1 / *h__;
/*<       HEE2=DD2/H >*/
    hee2 = *dd2 / *h__;
/*<       HEE3=DD3/H >*/
    hee3 = *dd3 / *h__;
/*<       GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB >*/
    switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L14;
	case 15:  goto L15;
    }

/*<    1  CONTINUE >*/
L1:
/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX */
/*<       DO  I=1,N  >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I) >*/
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
/*<          CONT(I)=F2(I)+Y0(I) >*/
	cont[i__] = f2[i__] + y0[i__];
/*<       END DO >*/
    }
/*<       CALL SOL (N,LDE1,E1,CONT,IP1)  >*/
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
/*<       GOTO 77 >*/
    goto L77;

/*<   11  CONTINUE >*/
L11:
/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO I=1,N  >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I) >*/
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
/*<          CONT(I)=F2(I)+Y0(I) >*/
	cont[i__] = f2[i__] + y0[i__];
/*<       END DO >*/
    }
/*<   48  MM=M1/M2 >*/
L48:
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM1=0.D0 >*/
	sum1 = 0.;
/*<          DO K=MM-1,0,-1 >*/
	for (k = mm - 1; k >= 0; --k) {
/*<             SUM1=(CONT(J+K*M2)+SUM1)/FAC1 >*/
	    sum1 = (cont[j + k * *m2] + sum1) / *fac1;
/*<             DO I=1,NM1 >*/
	    i__2 = *nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<                IM1=I+M1 >*/
		im1 = i__ + *m1;
/*<                CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1 >*/
		cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOL (NM1,LDE1,E1,CONT(M1+1),IP1)  >*/
    sol_(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
/*<       DO I=M1,1,-1 >*/
    for (i__ = *m1; i__ >= 1; --i__) {
/*<          CONT(I)=(CONT(I)+CONT(M2+I))/FAC1 >*/
	cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
/*<       END DO >*/
    }
/*<       GOTO 77 >*/
    goto L77;

/*<    2  CONTINUE >*/
L2:
/* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX */
/*<       DO I=1,N  >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I) >*/
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
/*<          CONT(I)=F2(I)+Y0(I) >*/
	cont[i__] = f2[i__] + y0[i__];
/*<       END DO >*/
    }
/*<       CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1) >*/
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[
	    1]);
/*<       GOTO 77 >*/
    goto L77;

/*<   12  CONTINUE >*/
L12:
/* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
/*<       DO I=1,N  >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I) >*/
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
/*<          CONT(I)=F2(I)+Y0(I) >*/
	cont[i__] = f2[i__] + y0[i__];
/*<       END DO >*/
    }
/*<   45  MM=M1/M2 >*/
L45:
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM1=0.D0 >*/
	sum1 = 0.;
/*<          DO K=MM-1,0,-1 >*/
	for (k = mm - 1; k >= 0; --k) {
/*<             SUM1=(CONT(J+K*M2)+SUM1)/FAC1 >*/
	    sum1 = (cont[j + k * *m2] + sum1) / *fac1;
/*<             DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC) >*/
/* Computing MAX */
	    i__2 = 1, i__3 = j - *mujac;
/* Computing MIN */
	    i__5 = *nm1, i__6 = j + *mljac;
	    i__4 = min(i__5,i__6);
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
/*<                IM1=I+M1 >*/
		im1 = i__ + *m1;
/*<                CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1 >*/
		cont[im1] += fjac[i__ + *mujac + 1 - j + (j + k * *m2) * 
			fjac_dim1] * sum1;
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1) >*/
    solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*m1 + 
	    1], &ip1[1]);
/*<       DO I=M1,1,-1 >*/
    for (i__ = *m1; i__ >= 1; --i__) {
/*<          CONT(I)=(CONT(I)+CONT(M2+I))/FAC1 >*/
	cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
/*<       END DO >*/
    }
/*<       GOTO 77 >*/
    goto L77;

/*<    3  CONTINUE >*/
L3:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I) >*/
	f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
/*<       END DO >*/
    }
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS) >*/
/* Computing MAX */
	i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__3 = min(i__5,i__6);
	for (j = max(i__4,i__2); j <= i__3; ++j) {
/*<             SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J) >*/
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
/*<          END DO >*/
	}
/*<          F2(I)=SUM >*/
	f2[i__] = sum;
/*<          CONT(I)=SUM+Y0(I) >*/
	cont[i__] = sum + y0[i__];
/*<       END DO >*/
    }
/*<       CALL SOL (N,LDE1,E1,CONT,IP1)  >*/
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
/*<       GOTO 77 >*/
    goto L77;

/*<   13  CONTINUE >*/
L13:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO I=1,M1 >*/
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I) >*/
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
/*<          CONT(I)=F2(I)+Y0(I) >*/
	cont[i__] = f2[i__] + y0[i__];
/*<       END DO >*/
    }
/*<       DO I=M1+1,N >*/
    i__1 = *n;
    for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
/*<          F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I) >*/
	f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
/*<       END DO >*/
    }
/*<       DO I=1,NM1 >*/
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS) >*/
/* Computing MAX */
	i__3 = 1, i__4 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *nm1, i__6 = i__ + *mumas;
	i__2 = min(i__5,i__6);
	for (j = max(i__3,i__4); j <= i__2; ++j) {
/*<             SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J+M1) >*/
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j + *
		    m1];
/*<          END DO >*/
	}
/*<          IM1=I+M1 >*/
	im1 = i__ + *m1;
/*<          F2(IM1)=SUM >*/
	f2[im1] = sum;
/*<          CONT(IM1)=SUM+Y0(IM1) >*/
	cont[im1] = sum + y0[im1];
/*<       END DO >*/
    }
/*<       GOTO 48 >*/
    goto L48;

/*<    4  CONTINUE >*/
L4:
/* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I) >*/
	f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
/*<       END DO >*/
    }
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS) >*/
/* Computing MAX */
	i__2 = 1, i__3 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__4 = min(i__5,i__6);
	for (j = max(i__2,i__3); j <= i__4; ++j) {
/*<             SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J) >*/
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
/*<          END DO >*/
	}
/*<          F2(I)=SUM >*/
	f2[i__] = sum;
/*<          CONT(I)=SUM+Y0(I) >*/
	cont[i__] = sum + y0[i__];
/*<       END DO >*/
    }
/*<       CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1) >*/
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[
	    1]);
/*<       GOTO 77 >*/
    goto L77;

/*<   14  CONTINUE >*/
L14:
/* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
/*<       DO I=1,M1 >*/
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I) >*/
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
/*<          CONT(I)=F2(I)+Y0(I) >*/
	cont[i__] = f2[i__] + y0[i__];
/*<       END DO >*/
    }
/*<       DO I=M1+1,N >*/
    i__1 = *n;
    for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
/*<          F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I) >*/
	f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
/*<       END DO >*/
    }
/*<       DO I=1,NM1 >*/
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS) >*/
/* Computing MAX */
	i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *nm1, i__6 = i__ + *mumas;
	i__3 = min(i__5,i__6);
	for (j = max(i__4,i__2); j <= i__3; ++j) {
/*<             SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J+M1) >*/
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j + *
		    m1];
/*<          END DO >*/
	}
/*<          IM1=I+M1 >*/
	im1 = i__ + *m1;
/*<          F2(IM1)=SUM >*/
	f2[im1] = sum;
/*<          CONT(IM1)=SUM+Y0(IM1) >*/
	cont[im1] = sum + y0[im1];
/*<       END DO >*/
    }
/*<       GOTO 45 >*/
    goto L45;

/*<    5  CONTINUE >*/
L5:
/* ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I) >*/
	f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
/*<       END DO >*/
    }
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO J=1,N >*/
	i__3 = *n;
	for (j = 1; j <= i__3; ++j) {
/*<             SUM=SUM+FMAS(I,J)*F1(J) >*/
	    sum += fmas[i__ + j * fmas_dim1] * f1[j];
/*<          END DO >*/
	}
/*<          F2(I)=SUM >*/
	f2[i__] = sum;
/*<          CONT(I)=SUM+Y0(I) >*/
	cont[i__] = sum + y0[i__];
/*<       END DO >*/
    }
/*<       CALL SOL (N,LDE1,E1,CONT,IP1)  >*/
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
/*<       GOTO 77 >*/
    goto L77;

/*<   15  CONTINUE >*/
L15:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO I=1,M1 >*/
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I) >*/
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
/*<          CONT(I)=F2(I)+Y0(I) >*/
	cont[i__] = f2[i__] + y0[i__];
/*<       END DO >*/
    }
/*<       DO I=M1+1,N >*/
    i__1 = *n;
    for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
/*<          F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I) >*/
	f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
/*<       END DO >*/
    }
/*<       DO I=1,NM1 >*/
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO J=1,NM1 >*/
	i__3 = *nm1;
	for (j = 1; j <= i__3; ++j) {
/*<             SUM=SUM+FMAS(I,J)*F1(J+M1) >*/
	    sum += fmas[i__ + j * fmas_dim1] * f1[j + *m1];
/*<          END DO >*/
	}
/*<          IM1=I+M1 >*/
	im1 = i__ + *m1;
/*<          F2(IM1)=SUM >*/
	f2[im1] = sum;
/*<          CONT(IM1)=SUM+Y0(IM1) >*/
	cont[im1] = sum + y0[im1];
/*<       END DO >*/
    }
/*<       GOTO 48 >*/
    goto L48;

/*<    6  CONTINUE >*/
L6:
/* ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
/* ------  THIS OPTION IS NOT PROVIDED */
/*<       RETURN >*/
    return 0;

/*<    7  CONTINUE >*/
L7:
/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
/*<       DO I=1,N  >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I) >*/
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
/*<          CONT(I)=F2(I)+Y0(I) >*/
	cont[i__] = f2[i__] + y0[i__];
/*<       END DO >*/
    }
/*<       DO MM=N-2,1,-1 >*/
    for (mm = *n - 2; mm >= 1; --mm) {
/*<          MP=N-MM >*/
	mp = *n - mm;
/*<          I=IPHES(MP) >*/
	i__ = iphes[mp];
/*<          IF (I.EQ.MP) GOTO 310 >*/
	if (i__ == mp) {
	    goto L310;
	}
/*<          ZSAFE=CONT(MP) >*/
	zsafe = cont[mp];
/*<          CONT(MP)=CONT(I) >*/
	cont[mp] = cont[i__];
/*<          CONT(I)=ZSAFE >*/
	cont[i__] = zsafe;
/*<  310     CONTINUE >*/
L310:
/*<          DO I=MP+1,N  >*/
	i__1 = *n;
	for (i__ = mp + 1; i__ <= i__1; ++i__) {
/*<             CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP) >*/
	    cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOLH(N,LDE1,E1,1,CONT,IP1) >*/
    solh_(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
/*<       DO MM=1,N-2 >*/
    i__1 = *n - 2;
    for (mm = 1; mm <= i__1; ++mm) {
/*<          MP=N-MM >*/
	mp = *n - mm;
/*<          DO I=MP+1,N  >*/
	i__3 = *n;
	for (i__ = mp + 1; i__ <= i__3; ++i__) {
/*<             CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP) >*/
	    cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
/*<          END DO >*/
	}
/*<          I=IPHES(MP) >*/
	i__ = iphes[mp];
/*<          IF (I.EQ.MP) GOTO 440 >*/
	if (i__ == mp) {
	    goto L440;
	}
/*<          ZSAFE=CONT(MP) >*/
	zsafe = cont[mp];
/*<          CONT(MP)=CONT(I) >*/
	cont[mp] = cont[i__];
/*<          CONT(I)=ZSAFE >*/
	cont[i__] = zsafe;
/*<  440     CONTINUE >*/
L440:
/*<       END DO >*/
	;
    }

/* -------------------------------------- */

/*<   77  CONTINUE >*/
L77:
/*<       ERR=0.D0 >*/
    *err = 0.;
/*<       DO  I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          ERR=ERR+(CONT(I)/SCAL(I))**2 >*/
/* Computing 2nd power */
	d__1 = cont[i__] / scal[i__];
	*err += d__1 * d__1;
/*<       END DO >*/
    }
/*<       ERR=MAX(SQRT(ERR/N),1.D-10) >*/
/* Computing MAX */
    d__1 = sqrt(*err / *n);
    *err = max(d__1,1e-10);

/*<       IF (ERR.LT.1.D0) RETURN >*/
    if (*err < 1.) {
	return 0;
    }
/*<       IF (FIRST.OR.REJECT) THEN >*/
    if (*first || *reject) {
/*<           DO I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<              CONT(I)=Y(I)+CONT(I) >*/
	    cont[i__] = y[i__] + cont[i__];
/*<           END DO >*/
	}
/*<           CALL FCN(N,X,CONT,F1,RPAR,IPAR) >*/
	(*fcn)(n, x, &cont[1], &f1[1], &rpar[1], &ipar[1]);
/*<           NFCN=NFCN+1 >*/
	++(*nfcn);
/*<           DO I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<              CONT(I)=F1(I)+F2(I) >*/
	    cont[i__] = f1[i__] + f2[i__];
/*<           END DO >*/
	}
/*<           GOTO (31,32,31,32,31,32,33,55,55,55,41,42,41,42,41), IJOB >*/
	switch (*ijob) {
	    case 1:  goto L31;
	    case 2:  goto L32;
	    case 3:  goto L31;
	    case 4:  goto L32;
	    case 5:  goto L31;
	    case 6:  goto L32;
	    case 7:  goto L33;
	    case 8:  goto L55;
	    case 9:  goto L55;
	    case 10:  goto L55;
	    case 11:  goto L41;
	    case 12:  goto L42;
	    case 13:  goto L41;
	    case 14:  goto L42;
	    case 15:  goto L41;
	}
/* ------ FULL MATRIX OPTION */
/*<   31      CONTINUE >*/
L31:
/*<           CALL SOL(N,LDE1,E1,CONT,IP1)  >*/
	sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
/*<           GOTO 88 >*/
	goto L88;
/* ------ FULL MATRIX OPTION, SECOND ORDER */
/*<  41      CONTINUE >*/
L41:
/*<          DO J=1,M2 >*/
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
/*<             SUM1=0.D0 >*/
	    sum1 = 0.;
/*<             DO K=MM-1,0,-1 >*/
	    for (k = mm - 1; k >= 0; --k) {
/*<                SUM1=(CONT(J+K*M2)+SUM1)/FAC1 >*/
		sum1 = (cont[j + k * *m2] + sum1) / *fac1;
/*<                DO I=1,NM1 >*/
		i__3 = *nm1;
		for (i__ = 1; i__ <= i__3; ++i__) {
/*<                   IM1=I+M1 >*/
		    im1 = i__ + *m1;
/*<                   CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1 >*/
		    cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
/*<                END DO >*/
		}
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<          CALL SOL(NM1,LDE1,E1,CONT(M1+1),IP1)  >*/
	sol_(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
/*<          DO I=M1,1,-1 >*/
	for (i__ = *m1; i__ >= 1; --i__) {
/*<             CONT(I)=(CONT(I)+CONT(M2+I))/FAC1 >*/
	    cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
/*<          END DO >*/
	}
/*<          GOTO 88 >*/
	goto L88;
/* ------ BANDED MATRIX OPTION */
/*<  32      CONTINUE >*/
L32:
/*<          CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1) >*/
	solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &
		ip1[1]);
/*<          GOTO 88 >*/
	goto L88;
/* ------ BANDED MATRIX OPTION, SECOND ORDER */
/*<  42      CONTINUE >*/
L42:
/*<          DO J=1,M2 >*/
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
/*<             SUM1=0.D0 >*/
	    sum1 = 0.;
/*<             DO K=MM-1,0,-1 >*/
	    for (k = mm - 1; k >= 0; --k) {
/*<                SUM1=(CONT(J+K*M2)+SUM1)/FAC1 >*/
		sum1 = (cont[j + k * *m2] + sum1) / *fac1;
/*<                DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC) >*/
/* Computing MAX */
		i__3 = 1, i__4 = j - *mujac;
/* Computing MIN */
		i__5 = *nm1, i__6 = j + *mljac;
		i__2 = min(i__5,i__6);
		for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
/*<                   IM1=I+M1 >*/
		    im1 = i__ + *m1;
/*<                   CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1 >*/
		    cont[im1] += fjac[i__ + *mujac + 1 - j + (j + k * *m2) * 
			    fjac_dim1] * sum1;
/*<                END DO >*/
		}
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<          CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1) >*/
	solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*
		m1 + 1], &ip1[1]);
/*<          DO I=M1,1,-1 >*/
	for (i__ = *m1; i__ >= 1; --i__) {
/*<             CONT(I)=(CONT(I)+CONT(M2+I))/FAC1 >*/
	    cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
/*<          END DO >*/
	}
/*<           GOTO 88 >*/
	goto L88;
/* ------ HESSENBERG MATRIX OPTION */
/*<   33      CONTINUE >*/
L33:
/*<           DO MM=N-2,1,-1 >*/
	for (mm = *n - 2; mm >= 1; --mm) {
/*<              MP=N-MM >*/
	    mp = *n - mm;
/*<              I=IPHES(MP) >*/
	    i__ = iphes[mp];
/*<              IF (I.EQ.MP) GOTO 510 >*/
	    if (i__ == mp) {
		goto L510;
	    }
/*<              ZSAFE=CONT(MP) >*/
	    zsafe = cont[mp];
/*<              CONT(MP)=CONT(I) >*/
	    cont[mp] = cont[i__];
/*<              CONT(I)=ZSAFE >*/
	    cont[i__] = zsafe;
/*<  510         CONTINUE >*/
L510:
/*<              DO I=MP+1,N  >*/
	    i__1 = *n;
	    for (i__ = mp + 1; i__ <= i__1; ++i__) {
/*<                 CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP) >*/
		cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
/*<              END DO >*/
	    }
/*<           END DO >*/
	}
/*<           CALL SOLH(N,LDE1,E1,1,CONT,IP1) >*/
	solh_(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
/*<           DO MM=1,N-2 >*/
	i__1 = *n - 2;
	for (mm = 1; mm <= i__1; ++mm) {
/*<              MP=N-MM >*/
	    mp = *n - mm;
/*<              DO I=MP+1,N  >*/
	    i__2 = *n;
	    for (i__ = mp + 1; i__ <= i__2; ++i__) {
/*<                 CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP) >*/
		cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
/*<              END DO >*/
	    }
/*<              I=IPHES(MP) >*/
	    i__ = iphes[mp];
/*<              IF (I.EQ.MP) GOTO 640 >*/
	    if (i__ == mp) {
		goto L640;
	    }
/*<              ZSAFE=CONT(MP) >*/
	    zsafe = cont[mp];
/*<              CONT(MP)=CONT(I) >*/
	    cont[mp] = cont[i__];
/*<              CONT(I)=ZSAFE >*/
	    cont[i__] = zsafe;
/*<  640         CONTINUE >*/
L640:
/*<           END DO >*/
	    ;
	}
/* ----------------------------------- */
/*<    88     CONTINUE >*/
L88:
/*<           ERR=0.D0  >*/
	*err = 0.;
/*<           DO I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<              ERR=ERR+(CONT(I)/SCAL(I))**2 >*/
/* Computing 2nd power */
	    d__1 = cont[i__] / scal[i__];
	    *err += d__1 * d__1;
/*<           END DO >*/
	}
/*<           ERR=MAX(SQRT(ERR/N),1.D-10) >*/
/* Computing MAX */
	d__1 = sqrt(*err / *n);
	*err = max(d__1,1e-10);
/*<        END IF >*/
    }
/*<        RETURN >*/
    return 0;
/* ----------------------------------------------------------- */
/*<   55   CONTINUE >*/
L55:
/*<        RETURN >*/
    return 0;
/*<        END >*/
} /* estrad_ */


/*     END OF SUBROUTINE ESTRAD */

/* *********************************************************** */

/*<        >*/
/* Subroutine */ int estrav_(integer *n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, doublereal *h__, doublereal *dd, S_fp 
	fcn, integer *nfcn, doublereal *y0, doublereal *y, integer *ijob, 
	doublereal *x, integer *m1, integer *m2, integer *nm1, integer *ns, 
	integer *nns, doublereal *e1, integer *lde1, doublereal *zz, 
	doublereal *cont, doublereal *ff, integer *ip1, integer *iphes, 
	doublereal *scal, doublereal *err, logical *first, logical *reject, 
	doublereal *fac1, doublereal *rpar, integer *ipar)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, mm, mp, im1;
    extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal sum, sum1;
    extern /* Subroutine */ int solb_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), solh_(integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *);
    static doublereal zsafe;

/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<        >*/
/*<       DIMENSION DD(NS),CONT(N),RPAR(1),IPAR(1) >*/
/*<       LOGICAL FIRST,REJECT >*/
/*<       COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG >*/
/*<       GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB >*/
    /* Parameter adjustments */
    --scal;
    --iphes;
    --cont;
    --y;
    --y0;
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip1;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    --dd;
    --ff;
    --zz;
    e1_dim1 = *lde1;
    e1_offset = 1 + e1_dim1;
    e1 -= e1_offset;
    --rpar;
    --ipar;

    /* Function Body */
    switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L14;
	case 15:  goto L15;
    }

/*<    1  CONTINUE >*/
L1:
/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX */
/*<       DO  I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=1,NS >*/
	i__2 = *ns;
	for (k = 1; k <= i__2; ++k) {
/*<             SUM=SUM+DD(K)*ZZ(I+(K-1)*N) >*/
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
/*<          END DO >*/
	}
/*<          FF(I+N)=SUM/H >*/
	ff[i__ + *n] = sum / *h__;
/*<          CONT(I)=FF(I+N)+Y0(I) >*/
	cont[i__] = ff[i__ + *n] + y0[i__];
/*<       END DO >*/
    }
/*<       CALL SOL (N,LDE1,E1,CONT,IP1)  >*/
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
/*<       GOTO 77 >*/
    goto L77;

/*<   11  CONTINUE >*/
L11:
/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO  I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=1,NS >*/
	i__2 = *ns;
	for (k = 1; k <= i__2; ++k) {
/*<             SUM=SUM+DD(K)*ZZ(I+(K-1)*N) >*/
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
/*<          END DO >*/
	}
/*<          FF(I+N)=SUM/H >*/
	ff[i__ + *n] = sum / *h__;
/*<          CONT(I)=FF(I+N)+Y0(I) >*/
	cont[i__] = ff[i__ + *n] + y0[i__];
/*<       END DO >*/
    }
/*<   48  MM=M1/M2 >*/
L48:
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM1=0.D0 >*/
	sum1 = 0.;
/*<          DO K=MM-1,0,-1 >*/
	for (k = mm - 1; k >= 0; --k) {
/*<             SUM1=(CONT(J+K*M2)+SUM1)/FAC1 >*/
	    sum1 = (cont[j + k * *m2] + sum1) / *fac1;
/*<             DO I=1,NM1 >*/
	    i__2 = *nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<                IM1=I+M1 >*/
		im1 = i__ + *m1;
/*<                CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1 >*/
		cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOL (NM1,LDE1,E1,CONT(M1+1),IP1)  >*/
    sol_(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
/*<       DO I=M1,1,-1 >*/
    for (i__ = *m1; i__ >= 1; --i__) {
/*<          CONT(I)=(CONT(I)+CONT(M2+I))/FAC1 >*/
	cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
/*<       END DO >*/
    }
/*<       GOTO 77 >*/
    goto L77;

/*<    2  CONTINUE >*/
L2:
/* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX */
/*<       DO  I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=1,NS >*/
	i__2 = *ns;
	for (k = 1; k <= i__2; ++k) {
/*<             SUM=SUM+DD(K)*ZZ(I+(K-1)*N) >*/
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
/*<          END DO >*/
	}
/*<          FF(I+N)=SUM/H >*/
	ff[i__ + *n] = sum / *h__;
/*<          CONT(I)=FF(I+N)+Y0(I) >*/
	cont[i__] = ff[i__ + *n] + y0[i__];
/*<       END DO >*/
    }
/*<       CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1) >*/
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[
	    1]);
/*<       GOTO 77 >*/
    goto L77;

/*<   12  CONTINUE >*/
L12:
/* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
/*<       DO  I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=1,NS >*/
	i__2 = *ns;
	for (k = 1; k <= i__2; ++k) {
/*<             SUM=SUM+DD(K)*ZZ(I+(K-1)*N) >*/
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
/*<          END DO >*/
	}
/*<          FF(I+N)=SUM/H >*/
	ff[i__ + *n] = sum / *h__;
/*<          CONT(I)=FF(I+N)+Y0(I) >*/
	cont[i__] = ff[i__ + *n] + y0[i__];
/*<       END DO >*/
    }
/*<   45  MM=M1/M2 >*/
L45:
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM1=0.D0 >*/
	sum1 = 0.;
/*<          DO K=MM-1,0,-1 >*/
	for (k = mm - 1; k >= 0; --k) {
/*<             SUM1=(CONT(J+K*M2)+SUM1)/FAC1 >*/
	    sum1 = (cont[j + k * *m2] + sum1) / *fac1;
/*<             DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC) >*/
/* Computing MAX */
	    i__2 = 1, i__3 = j - *mujac;
/* Computing MIN */
	    i__5 = *nm1, i__6 = j + *mljac;
	    i__4 = min(i__5,i__6);
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
/*<                IM1=I+M1 >*/
		im1 = i__ + *m1;
/*<                CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1 >*/
		cont[im1] += fjac[i__ + *mujac + 1 - j + (j + k * *m2) * 
			fjac_dim1] * sum1;
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1) >*/
    solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*m1 + 
	    1], &ip1[1]);
/*<       DO I=M1,1,-1 >*/
    for (i__ = *m1; i__ >= 1; --i__) {
/*<          CONT(I)=(CONT(I)+CONT(M2+I))/FAC1 >*/
	cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
/*<       END DO >*/
    }
/*<       GOTO 77 >*/
    goto L77;

/*<    3  CONTINUE >*/
L3:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
/*<       DO  I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=1,NS >*/
	i__4 = *ns;
	for (k = 1; k <= i__4; ++k) {
/*<             SUM=SUM+DD(K)*ZZ(I+(K-1)*N) >*/
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
/*<          END DO >*/
	}
/*<          FF(I)=SUM/H >*/
	ff[i__] = sum / *h__;
/*<       END DO >*/
    }
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS) >*/
/* Computing MAX */
	i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__3 = min(i__5,i__6);
	for (j = max(i__4,i__2); j <= i__3; ++j) {
/*<             SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J) >*/
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j];
/*<          END DO >*/
	}
/*<          FF(I+N)=SUM >*/
	ff[i__ + *n] = sum;
/*<          CONT(I)=SUM+Y0(I) >*/
	cont[i__] = sum + y0[i__];
/*<       END DO >*/
    }
/*<       CALL SOL (N,LDE1,E1,CONT,IP1)  >*/
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
/*<       GOTO 77 >*/
    goto L77;

/*<   13  CONTINUE >*/
L13:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO  I=1,M1 >*/
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=1,NS >*/
	i__3 = *ns;
	for (k = 1; k <= i__3; ++k) {
/*<             SUM=SUM+DD(K)*ZZ(I+(K-1)*N) >*/
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
/*<          END DO >*/
	}
/*<          FF(I+N)=SUM/H >*/
	ff[i__ + *n] = sum / *h__;
/*<          CONT(I)=FF(I+N)+Y0(I) >*/
	cont[i__] = ff[i__ + *n] + y0[i__];
/*<       END DO >*/
    }
/*<       DO I=M1+1,N >*/
    i__1 = *n;
    for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=1,NS >*/
	i__3 = *ns;
	for (k = 1; k <= i__3; ++k) {
/*<             SUM=SUM+DD(K)*ZZ(I+(K-1)*N) >*/
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
/*<          END DO >*/
	}
/*<          FF(I)=SUM/H >*/
	ff[i__] = sum / *h__;
/*<       END DO >*/
    }
/*<       DO I=1,NM1 >*/
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS) >*/
/* Computing MAX */
	i__3 = 1, i__4 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *nm1, i__6 = i__ + *mumas;
	i__2 = min(i__5,i__6);
	for (j = max(i__3,i__4); j <= i__2; ++j) {
/*<             SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J+M1) >*/
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j + *
		    m1];
/*<          END DO >*/
	}
/*<          IM1=I+M1 >*/
	im1 = i__ + *m1;
/*<          FF(IM1+N)=SUM >*/
	ff[im1 + *n] = sum;
/*<          CONT(IM1)=SUM+Y0(IM1) >*/
	cont[im1] = sum + y0[im1];
/*<       END DO >*/
    }
/*<       GOTO 48 >*/
    goto L48;

/*<    4  CONTINUE >*/
L4:
/* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
/*<       DO  I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=1,NS >*/
	i__2 = *ns;
	for (k = 1; k <= i__2; ++k) {
/*<             SUM=SUM+DD(K)*ZZ(I+(K-1)*N) >*/
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
/*<          END DO >*/
	}
/*<          FF(I)=SUM/H >*/
	ff[i__] = sum / *h__;
/*<       END DO >*/
    }
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS) >*/
/* Computing MAX */
	i__2 = 1, i__3 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__4 = min(i__5,i__6);
	for (j = max(i__2,i__3); j <= i__4; ++j) {
/*<             SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J) >*/
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j];
/*<          END DO >*/
	}
/*<          FF(I+N)=SUM >*/
	ff[i__ + *n] = sum;
/*<          CONT(I)=SUM+Y0(I) >*/
	cont[i__] = sum + y0[i__];
/*<       END DO >*/
    }
/*<       CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1) >*/
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[
	    1]);
/*<       GOTO 77 >*/
    goto L77;

/*<   14  CONTINUE >*/
L14:
/* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
/*<       DO  I=1,M1 >*/
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=1,NS >*/
	i__4 = *ns;
	for (k = 1; k <= i__4; ++k) {
/*<             SUM=SUM+DD(K)*ZZ(I+(K-1)*N) >*/
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
/*<          END DO >*/
	}
/*<          FF(I+N)=SUM/H >*/
	ff[i__ + *n] = sum / *h__;
/*<          CONT(I)=FF(I+N)+Y0(I) >*/
	cont[i__] = ff[i__ + *n] + y0[i__];
/*<       END DO >*/
    }
/*<       DO I=M1+1,N >*/
    i__1 = *n;
    for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=1,NS >*/
	i__4 = *ns;
	for (k = 1; k <= i__4; ++k) {
/*<             SUM=SUM+DD(K)*ZZ(I+(K-1)*N) >*/
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
/*<          END DO >*/
	}
/*<          FF(I)=SUM/H >*/
	ff[i__] = sum / *h__;
/*<       END DO >*/
    }
/*<       DO I=1,NM1 >*/
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS) >*/
/* Computing MAX */
	i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *nm1, i__6 = i__ + *mumas;
	i__3 = min(i__5,i__6);
	for (j = max(i__4,i__2); j <= i__3; ++j) {
/*<             SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J+M1) >*/
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j + *
		    m1];
/*<          END DO >*/
	}
/*<          IM1=I+M1 >*/
	im1 = i__ + *m1;
/*<          FF(IM1+N)=SUM >*/
	ff[im1 + *n] = sum;
/*<          CONT(IM1)=SUM+Y0(IM1) >*/
	cont[im1] = sum + y0[im1];
/*<       END DO >*/
    }
/*<       GOTO 45 >*/
    goto L45;

/*<    5  CONTINUE >*/
L5:
/* ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=1,NS >*/
	i__3 = *ns;
	for (k = 1; k <= i__3; ++k) {
/*<             SUM=SUM+DD(K)*ZZ(I+(K-1)*N) >*/
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
/*<          END DO >*/
	}
/*<          FF(I)=SUM/H >*/
	ff[i__] = sum / *h__;
/*<       END DO >*/
    }
/*<       DO I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO J=1,N >*/
	i__3 = *n;
	for (j = 1; j <= i__3; ++j) {
/*<             SUM=SUM+FMAS(I,J)*FF(J) >*/
	    sum += fmas[i__ + j * fmas_dim1] * ff[j];
/*<          END DO >*/
	}
/*<          FF(I+N)=SUM >*/
	ff[i__ + *n] = sum;
/*<          CONT(I)=SUM+Y0(I) >*/
	cont[i__] = sum + y0[i__];
/*<       END DO >*/
    }
/*<       CALL SOL (N,LDE1,E1,CONT,IP1)  >*/
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
/*<       GOTO 77 >*/
    goto L77;

/*<   15  CONTINUE >*/
L15:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       DO  I=1,M1 >*/
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=1,NS >*/
	i__3 = *ns;
	for (k = 1; k <= i__3; ++k) {
/*<             SUM=SUM+DD(K)*ZZ(I+(K-1)*N) >*/
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
/*<          END DO >*/
	}
/*<          FF(I+N)=SUM/H >*/
	ff[i__ + *n] = sum / *h__;
/*<          CONT(I)=FF(I+N)+Y0(I) >*/
	cont[i__] = ff[i__ + *n] + y0[i__];
/*<       END DO >*/
    }
/*<       DO I=M1+1,N >*/
    i__1 = *n;
    for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=1,NS >*/
	i__3 = *ns;
	for (k = 1; k <= i__3; ++k) {
/*<             SUM=SUM+DD(K)*ZZ(I+(K-1)*N) >*/
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
/*<          END DO >*/
	}
/*<          FF(I)=SUM/H >*/
	ff[i__] = sum / *h__;
/*<       END DO >*/
    }
/*<       DO I=1,NM1 >*/
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO J=1,NM1 >*/
	i__3 = *nm1;
	for (j = 1; j <= i__3; ++j) {
/*<             SUM=SUM+FMAS(I,J)*FF(J+M1) >*/
	    sum += fmas[i__ + j * fmas_dim1] * ff[j + *m1];
/*<          END DO >*/
	}
/*<          IM1=I+M1 >*/
	im1 = i__ + *m1;
/*<          FF(IM1+N)=SUM >*/
	ff[im1 + *n] = sum;
/*<          CONT(IM1)=SUM+Y0(IM1) >*/
	cont[im1] = sum + y0[im1];
/*<       END DO >*/
    }
/*<       GOTO 48 >*/
    goto L48;

/*<    6  CONTINUE >*/
L6:
/* ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
/* ------  THIS OPTION IS NOT PROVIDED */
/*<       RETURN >*/
    return 0;

/*<    7  CONTINUE >*/
L7:
/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
/*<       DO  I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=1,NS >*/
	i__3 = *ns;
	for (k = 1; k <= i__3; ++k) {
/*<             SUM=SUM+DD(K)*ZZ(I+(K-1)*N) >*/
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
/*<          END DO >*/
	}
/*<          FF(I+N)=SUM/H >*/
	ff[i__ + *n] = sum / *h__;
/*<          CONT(I)=FF(I+N)+Y0(I) >*/
	cont[i__] = ff[i__ + *n] + y0[i__];
/*<       END DO >*/
    }
/*<       DO MM=N-2,1,-1 >*/
    for (mm = *n - 2; mm >= 1; --mm) {
/*<          MP=N-MM >*/
	mp = *n - mm;
/*<          I=IPHES(MP) >*/
	i__ = iphes[mp];
/*<          IF (I.EQ.MP) GOTO 310 >*/
	if (i__ == mp) {
	    goto L310;
	}
/*<          ZSAFE=CONT(MP) >*/
	zsafe = cont[mp];
/*<          CONT(MP)=CONT(I) >*/
	cont[mp] = cont[i__];
/*<          CONT(I)=ZSAFE >*/
	cont[i__] = zsafe;
/*<  310     CONTINUE >*/
L310:
/*<          DO I=MP+1,N  >*/
	i__1 = *n;
	for (i__ = mp + 1; i__ <= i__1; ++i__) {
/*<             CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP) >*/
	    cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOLH(N,LDE1,E1,1,CONT,IP1) >*/
    solh_(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
/*<       DO MM=1,N-2 >*/
    i__1 = *n - 2;
    for (mm = 1; mm <= i__1; ++mm) {
/*<          MP=N-MM >*/
	mp = *n - mm;
/*<          DO I=MP+1,N  >*/
	i__3 = *n;
	for (i__ = mp + 1; i__ <= i__3; ++i__) {
/*<             CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP) >*/
	    cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
/*<          END DO >*/
	}
/*<          I=IPHES(MP) >*/
	i__ = iphes[mp];
/*<          IF (I.EQ.MP) GOTO 440 >*/
	if (i__ == mp) {
	    goto L440;
	}
/*<          ZSAFE=CONT(MP) >*/
	zsafe = cont[mp];
/*<          CONT(MP)=CONT(I) >*/
	cont[mp] = cont[i__];
/*<          CONT(I)=ZSAFE >*/
	cont[i__] = zsafe;
/*<  440     CONTINUE >*/
L440:
/*<       END DO >*/
	;
    }

/* -------------------------------------- */

/*<   77  CONTINUE >*/
L77:
/*<       ERR=0.D0 >*/
    *err = 0.;
/*<       DO  I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          ERR=ERR+(CONT(I)/SCAL(I))**2 >*/
/* Computing 2nd power */
	d__1 = cont[i__] / scal[i__];
	*err += d__1 * d__1;
/*<       END DO >*/
    }
/*<       ERR=MAX(SQRT(ERR/N),1.D-10) >*/
/* Computing MAX */
    d__1 = sqrt(*err / *n);
    *err = max(d__1,1e-10);

/*<       IF (ERR.LT.1.D0) RETURN >*/
    if (*err < 1.) {
	return 0;
    }
/*<       IF (FIRST.OR.REJECT) THEN >*/
    if (*first || *reject) {
/*<           DO I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<              CONT(I)=Y(I)+CONT(I) >*/
	    cont[i__] = y[i__] + cont[i__];
/*<           END DO >*/
	}
/*<           CALL FCN(N,X,CONT,FF,RPAR,IPAR) >*/
	(*fcn)(n, x, &cont[1], &ff[1], &rpar[1], &ipar[1]);
/*<           NFCN=NFCN+1 >*/
	++(*nfcn);
/*<           DO I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<              CONT(I)=FF(I)+FF(I+N) >*/
	    cont[i__] = ff[i__] + ff[i__ + *n];
/*<           END DO >*/
	}
/*<           GOTO (31,32,31,32,31,32,33,55,55,55,41,42,41,42,41), IJOB >*/
	switch (*ijob) {
	    case 1:  goto L31;
	    case 2:  goto L32;
	    case 3:  goto L31;
	    case 4:  goto L32;
	    case 5:  goto L31;
	    case 6:  goto L32;
	    case 7:  goto L33;
	    case 8:  goto L55;
	    case 9:  goto L55;
	    case 10:  goto L55;
	    case 11:  goto L41;
	    case 12:  goto L42;
	    case 13:  goto L41;
	    case 14:  goto L42;
	    case 15:  goto L41;
	}
/* ------ FULL MATRIX OPTION */
/*<  31      CONTINUE >*/
L31:
/*<          CALL SOL (N,LDE1,E1,CONT,IP1)  >*/
	sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
/*<           GOTO 88 >*/
	goto L88;
/* ------ FULL MATRIX OPTION, SECOND ORDER */
/*<  41      CONTINUE >*/
L41:
/*<          DO J=1,M2 >*/
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
/*<             SUM1=0.D0 >*/
	    sum1 = 0.;
/*<             DO K=MM-1,0,-1 >*/
	    for (k = mm - 1; k >= 0; --k) {
/*<                SUM1=(CONT(J+K*M2)+SUM1)/FAC1 >*/
		sum1 = (cont[j + k * *m2] + sum1) / *fac1;
/*<                DO I=1,NM1 >*/
		i__3 = *nm1;
		for (i__ = 1; i__ <= i__3; ++i__) {
/*<                   IM1=I+M1 >*/
		    im1 = i__ + *m1;
/*<                   CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1 >*/
		    cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
/*<                END DO >*/
		}
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<          CALL SOL (NM1,LDE1,E1,CONT(M1+1),IP1)  >*/
	sol_(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
/*<          DO I=M1,1,-1 >*/
	for (i__ = *m1; i__ >= 1; --i__) {
/*<             CONT(I)=(CONT(I)+CONT(M2+I))/FAC1 >*/
	    cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
/*<          END DO >*/
	}
/*<           GOTO 88 >*/
	goto L88;
/* ------ BANDED MATRIX OPTION */
/*<  32      CONTINUE >*/
L32:
/*<          CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1) >*/
	solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &
		ip1[1]);
/*<           GOTO 88 >*/
	goto L88;
/* ------ BANDED MATRIX OPTION, SECOND ORDER */
/*<  42      CONTINUE >*/
L42:
/*<          DO J=1,M2 >*/
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
/*<             SUM1=0.D0 >*/
	    sum1 = 0.;
/*<             DO K=MM-1,0,-1 >*/
	    for (k = mm - 1; k >= 0; --k) {
/*<                SUM1=(CONT(J+K*M2)+SUM1)/FAC1 >*/
		sum1 = (cont[j + k * *m2] + sum1) / *fac1;
/*<                DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC) >*/
/* Computing MAX */
		i__3 = 1, i__4 = j - *mujac;
/* Computing MIN */
		i__5 = *nm1, i__6 = j + *mljac;
		i__2 = min(i__5,i__6);
		for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
/*<                   IM1=I+M1 >*/
		    im1 = i__ + *m1;
/*<                   CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1 >*/
		    cont[im1] += fjac[i__ + *mujac + 1 - j + (j + k * *m2) * 
			    fjac_dim1] * sum1;
/*<                END DO >*/
		}
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<          CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1) >*/
	solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*
		m1 + 1], &ip1[1]);
/*<          DO I=M1,1,-1 >*/
	for (i__ = *m1; i__ >= 1; --i__) {
/*<             CONT(I)=(CONT(I)+CONT(M2+I))/FAC1 >*/
	    cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
/*<          END DO >*/
	}
/*<           GOTO 88 >*/
	goto L88;
/* ------ HESSENBERG MATRIX OPTION */
/*<   33      CONTINUE >*/
L33:
/*<           DO MM=N-2,1,-1 >*/
	for (mm = *n - 2; mm >= 1; --mm) {
/*<              MP=N-MM >*/
	    mp = *n - mm;
/*<              I=IPHES(MP) >*/
	    i__ = iphes[mp];
/*<              IF (I.EQ.MP) GOTO 510 >*/
	    if (i__ == mp) {
		goto L510;
	    }
/*<              ZSAFE=CONT(MP) >*/
	    zsafe = cont[mp];
/*<              CONT(MP)=CONT(I) >*/
	    cont[mp] = cont[i__];
/*<              CONT(I)=ZSAFE >*/
	    cont[i__] = zsafe;
/*<  510         CONTINUE >*/
L510:
/*<              DO I=MP+1,N  >*/
	    i__1 = *n;
	    for (i__ = mp + 1; i__ <= i__1; ++i__) {
/*<                 CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP) >*/
		cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
/*<              END DO >*/
	    }
/*<           END DO >*/
	}
/*<           CALL SOLH(N,LDE1,E1,1,CONT,IP1) >*/
	solh_(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
/*<           DO MM=1,N-2 >*/
	i__1 = *n - 2;
	for (mm = 1; mm <= i__1; ++mm) {
/*<              MP=N-MM >*/
	    mp = *n - mm;
/*<              DO I=MP+1,N  >*/
	    i__2 = *n;
	    for (i__ = mp + 1; i__ <= i__2; ++i__) {
/*<                 CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP) >*/
		cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
/*<              END DO >*/
	    }
/*<              I=IPHES(MP) >*/
	    i__ = iphes[mp];
/*<              IF (I.EQ.MP) GOTO 640 >*/
	    if (i__ == mp) {
		goto L640;
	    }
/*<              ZSAFE=CONT(MP) >*/
	    zsafe = cont[mp];
/*<              CONT(MP)=CONT(I) >*/
	    cont[mp] = cont[i__];
/*<              CONT(I)=ZSAFE >*/
	    cont[i__] = zsafe;
/*<  640         CONTINUE >*/
L640:
/*<           END DO >*/
	    ;
	}
/* ----------------------------------- */
/*<   88      CONTINUE >*/
L88:
/*<           ERR=0.D0  >*/
	*err = 0.;
/*<           DO I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<              ERR=ERR+(CONT(I)/SCAL(I))**2 >*/
/* Computing 2nd power */
	    d__1 = cont[i__] / scal[i__];
	    *err += d__1 * d__1;
/*<           END DO >*/
	}
/*<           ERR=MAX(SQRT(ERR/N),1.D-10) >*/
/* Computing MAX */
	d__1 = sqrt(*err / *n);
	*err = max(d__1,1e-10);
/*<        END IF >*/
    }
/*<        RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   55  CONTINUE >*/
L55:
/*<       RETURN >*/
    return 0;
/*<        END >*/
} /* estrav_ */


/*     END OF SUBROUTINE ESTRAV */

/* *********************************************************** */

/*<        >*/
/* Subroutine */ int slvrod_(integer *n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, integer *m1, integer *m2, integer *
	nm1, doublereal *fac1, doublereal *e, integer *lde, integer *ip, 
	doublereal *dy, doublereal *ak, doublereal *fx, doublereal *ynew, 
	doublereal *hd, integer *ijob, logical *stage1)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e_dim1, e_offset, 
	    i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, j, k, mm, im1, jkm;
    extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal sum;
    extern /* Subroutine */ int solb_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *);

/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<        >*/
/*<       LOGICAL STAGE1 >*/
/*<       COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG >*/

/*<       IF (HD.EQ.0.D0) THEN >*/
    /* Parameter adjustments */
    --ynew;
    --fx;
    --ak;
    --dy;
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

    /* Function Body */
    if (*hd == 0.) {
/*<          DO  I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<            AK(I)=DY(I) >*/
	    ak[i__] = dy[i__];
/*<          END DO >*/
	}
/*<       ELSE >*/
    } else {
/*<          DO I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             AK(I)=DY(I)+HD*FX(I) >*/
	    ak[i__] = dy[i__] + *hd * fx[i__];
/*<          END DO >*/
	}
/*<       END IF >*/
    }

/*<       GOTO (1,2,3,4,5,6,55,55,55,55,11,12,13,13,15), IJOB >*/
    switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L55;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L13;
	case 15:  goto L15;
    }

/* ----------------------------------------------------------- */

/*<    1  CONTINUE >*/
L1:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
/*<       IF (STAGE1) THEN >*/
    if (*stage1) {
/*<          DO I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             AK(I)=AK(I)+YNEW(I) >*/
	    ak[i__] += ynew[i__];
/*<          END DO >*/
	}
/*<       END IF >*/
    }
/*<       CALL SOL (N,LDE,E,AK,IP) >*/
    sol_(n, lde, &e[e_offset], &ak[1], &ip[1]);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   11  CONTINUE >*/
L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       IF (STAGE1) THEN >*/
    if (*stage1) {
/*<          DO I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             AK(I)=AK(I)+YNEW(I) >*/
	    ak[i__] += ynew[i__];
/*<          END DO >*/
	}
/*<       END IF >*/
    }
/*<  48   MM=M1/M2 >*/
L48:
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=MM-1,0,-1 >*/
	for (k = mm - 1; k >= 0; --k) {
/*<             JKM=J+K*M2 >*/
	    jkm = j + k * *m2;
/*<             SUM=(AK(JKM)+SUM)/FAC1 >*/
	    sum = (ak[jkm] + sum) / *fac1;
/*<             DO I=1,NM1 >*/
	    i__2 = *nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<                IM1=I+M1 >*/
		im1 = i__ + *m1;
/*<                AK(IM1)=AK(IM1)+FJAC(I,JKM)*SUM >*/
		ak[im1] += fjac[i__ + jkm * fjac_dim1] * sum;
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOL (NM1,LDE,E,AK(M1+1),IP) >*/
    sol_(nm1, lde, &e[e_offset], &ak[*m1 + 1], &ip[1]);
/*<       DO I=M1,1,-1 >*/
    for (i__ = *m1; i__ >= 1; --i__) {
/*<          AK(I)=(AK(I)+AK(M2+I))/FAC1 >*/
	ak[i__] = (ak[i__] + ak[*m2 + i__]) / *fac1;
/*<       END DO >*/
    }
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    2  CONTINUE >*/
L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
/*<       IF (STAGE1) THEN >*/
    if (*stage1) {
/*<          DO I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             AK(I)=AK(I)+YNEW(I) >*/
	    ak[i__] += ynew[i__];
/*<          END DO >*/
	}
/*<       END IF >*/
    }
/*<       CALL SOLB (N,LDE,E,MLE,MUE,AK,IP) >*/
    solb_(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[1], &ip[1]);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   12  CONTINUE >*/
L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
/*<       IF (STAGE1) THEN >*/
    if (*stage1) {
/*<          DO I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             AK(I)=AK(I)+YNEW(I) >*/
	    ak[i__] += ynew[i__];
/*<          END DO >*/
	}
/*<       END IF >*/
    }
/*<   45  MM=M1/M2 >*/
L45:
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=MM-1,0,-1 >*/
	for (k = mm - 1; k >= 0; --k) {
/*<             JKM=J+K*M2 >*/
	    jkm = j + k * *m2;
/*<             SUM=(AK(JKM)+SUM)/FAC1 >*/
	    sum = (ak[jkm] + sum) / *fac1;
/*<             DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC) >*/
/* Computing MAX */
	    i__2 = 1, i__3 = j - *mujac;
/* Computing MIN */
	    i__5 = *nm1, i__6 = j + *mljac;
	    i__4 = min(i__5,i__6);
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
/*<                IM1=I+M1 >*/
		im1 = i__ + *m1;
/*<                AK(IM1)=AK(IM1)+FJAC(I+MUJAC+1-J,JKM)*SUM >*/
		ak[im1] += fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1] * sum;
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOLB (NM1,LDE,E,MLE,MUE,AK(M1+1),IP) >*/
    solb_(nm1, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[*m1 + 1], &
	    ip[1]);
/*<       DO I=M1,1,-1 >*/
    for (i__ = *m1; i__ >= 1; --i__) {
/*<          AK(I)=(AK(I)+AK(M2+I))/FAC1 >*/
	ak[i__] = (ak[i__] + ak[*m2 + i__]) / *fac1;
/*<       END DO >*/
    }
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    3  CONTINUE >*/
L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
/*<       IF (STAGE1) THEN >*/
    if (*stage1) {
/*<       DO  I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	    sum = 0.;
/*<          DO  J=MAX(1,I-MLMAS),MIN(N,I+MUMAS) >*/
/* Computing MAX */
	    i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	    i__5 = *n, i__6 = i__ + *mumas;
	    i__3 = min(i__5,i__6);
	    for (j = max(i__4,i__2); j <= i__3; ++j) {
/*<             SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J) >*/
		sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ynew[
			j];
/*<          END DO >*/
	    }
/*<          AK(I)=AK(I)+SUM >*/
	    ak[i__] += sum;
/*<       END DO >*/
	}
/*<       END IF >*/
    }
/*<       CALL SOL (N,LDE,E,AK,IP) >*/
    sol_(n, lde, &e[e_offset], &ak[1], &ip[1]);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   13  CONTINUE >*/
L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       IF (STAGE1) THEN >*/
    if (*stage1) {
/*<          DO I=1,M1 >*/
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             AK(I)=AK(I)+YNEW(I) >*/
	    ak[i__] += ynew[i__];
/*<          END DO >*/
	}
/*<          DO I=1,NM1 >*/
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             SUM=0.D0 >*/
	    sum = 0.;
/*<             DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS) >*/
/* Computing MAX */
	    i__3 = 1, i__4 = i__ - *mlmas;
/* Computing MIN */
	    i__5 = *nm1, i__6 = i__ + *mumas;
	    i__2 = min(i__5,i__6);
	    for (j = max(i__3,i__4); j <= i__2; ++j) {
/*<                 SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J+M1) >*/
		sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ynew[
			j + *m1];
/*<             END DO >*/
	    }
/*<             IM1=I+M1 >*/
	    im1 = i__ + *m1;
/*<             AK(IM1)=AK(IM1)+SUM >*/
	    ak[im1] += sum;
/*<          END DO >*/
	}
/*<       END IF >*/
    }
/*<       IF (IJOB.EQ.14) GOTO 45 >*/
    if (*ijob == 14) {
	goto L45;
    }
/*<       GOTO 48 >*/
    goto L48;

/* ----------------------------------------------------------- */

/*<    4  CONTINUE >*/
L4:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
/*<       IF (STAGE1) THEN >*/
    if (*stage1) {
/*<       DO I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	    sum = 0.;
/*<          DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS) >*/
/* Computing MAX */
	    i__2 = 1, i__3 = i__ - *mlmas;
/* Computing MIN */
	    i__5 = *n, i__6 = i__ + *mumas;
	    i__4 = min(i__5,i__6);
	    for (j = max(i__2,i__3); j <= i__4; ++j) {
/*<             SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J) >*/
		sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ynew[
			j];
/*<          END DO >*/
	    }
/*<          AK(I)=AK(I)+SUM >*/
	    ak[i__] += sum;
/*<       END DO >*/
	}
/*<       END IF >*/
    }
/*<       CALL SOLB (N,LDE,E,MLE,MUE,AK,IP) >*/
    solb_(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[1], &ip[1]);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    5  CONTINUE >*/
L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
/*<       IF (STAGE1) THEN >*/
    if (*stage1) {
/*<       DO I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	    sum = 0.;
/*<          DO J=1,N >*/
	    i__4 = *n;
	    for (j = 1; j <= i__4; ++j) {
/*<             SUM=SUM+FMAS(I,J)*YNEW(J) >*/
		sum += fmas[i__ + j * fmas_dim1] * ynew[j];
/*<          END DO >*/
	    }
/*<          AK(I)=AK(I)+SUM >*/
	    ak[i__] += sum;
/*<       END DO >*/
	}
/*<       END IF >*/
    }
/*<       CALL SOL (N,LDE,E,AK,IP) >*/
    sol_(n, lde, &e[e_offset], &ak[1], &ip[1]);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   15  CONTINUE >*/
L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       IF (STAGE1) THEN >*/
    if (*stage1) {
/*<          DO I=1,M1 >*/
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             AK(I)=AK(I)+YNEW(I) >*/
	    ak[i__] += ynew[i__];
/*<          END DO >*/
	}
/*<          DO I=1,NM1 >*/
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             SUM=0.D0 >*/
	    sum = 0.;
/*<             DO J=1,NM1 >*/
	    i__4 = *nm1;
	    for (j = 1; j <= i__4; ++j) {
/*<                SUM=SUM+FMAS(I,J)*YNEW(J+M1) >*/
		sum += fmas[i__ + j * fmas_dim1] * ynew[j + *m1];
/*<             END DO >*/
	    }
/*<             IM1=I+M1 >*/
	    im1 = i__ + *m1;
/*<             AK(IM1)=AK(IM1)+SUM >*/
	    ak[im1] += sum;
/*<          END DO >*/
	}
/*<       END IF >*/
    }
/*<       GOTO 48 >*/
    goto L48;

/* ----------------------------------------------------------- */

/*<    6  CONTINUE >*/
L6:
/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
/* ---  THIS OPTION IS NOT PROVIDED */
/*<       IF (STAGE1) THEN >*/
    if (*stage1) {
/*<       DO 624 I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SUM=0.D0 >*/
	    sum = 0.;
/*<          DO 623 J=1,N >*/
	    i__4 = *n;
	    for (j = 1; j <= i__4; ++j) {
/*<   623       SUM=SUM+FMAS(I,J)*YNEW(J) >*/
/* L623: */
		sum += fmas[i__ + j * fmas_dim1] * ynew[j];
	    }
/*<   624    AK(I)=AK(I)+SUM >*/
/* L624: */
	    ak[i__] += sum;
	}
/*<       CALL SOLB (N,LDE,E,MLE,MUE,AK,IP) >*/
	solb_(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[1], &ip[1]
		);
/*<       END IF >*/
    }
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   55  CONTINUE >*/
L55:
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* slvrod_ */


/*     END OF SUBROUTINE SLVROD */


/* *********************************************************** */

/*<        >*/
/* Subroutine */ int slvseu_(integer *n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, integer *m1, integer *m2, integer *
	nm1, doublereal *fac1, doublereal *e, integer *lde, integer *ip, 
	integer *iphes, doublereal *del, integer *ijob)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e_dim1, e_offset, 
	    i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, j, k, mm, mp, im1, mp1, jkm, mmm;
    extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal sum;
    extern /* Subroutine */ int solb_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), solh_(integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *);
    static doublereal zsafe;

/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<       DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E(LDE,NM1),DEL(N) >*/
/*<       DIMENSION IP(NM1),IPHES(N) >*/
/*<       COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG >*/

/*<       GOTO (1,2,1,2,1,55,7,55,55,55,11,12,11,12,11), IJOB >*/
    /* Parameter adjustments */
    --del;
    --iphes;
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

    /* Function Body */
    switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L1;
	case 4:  goto L2;
	case 5:  goto L1;
	case 6:  goto L55;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L11;
	case 14:  goto L12;
	case 15:  goto L11;
    }

/* ----------------------------------------------------------- */

/*<    1  CONTINUE >*/
L1:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
/*<       CALL SOL (N,LDE,E,DEL,IP) >*/
    sol_(n, lde, &e[e_offset], &del[1], &ip[1]);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   11  CONTINUE >*/
L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
/*<       MM=M1/M2 >*/
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=MM-1,0,-1 >*/
	for (k = mm - 1; k >= 0; --k) {
/*<             JKM=J+K*M2 >*/
	    jkm = j + k * *m2;
/*<             SUM=(DEL(JKM)+SUM)/FAC1 >*/
	    sum = (del[jkm] + sum) / *fac1;
/*<             DO I=1,NM1 >*/
	    i__2 = *nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<                IM1=I+M1 >*/
		im1 = i__ + *m1;
/*<                DEL(IM1)=DEL(IM1)+FJAC(I,JKM)*SUM >*/
		del[im1] += fjac[i__ + jkm * fjac_dim1] * sum;
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOL (NM1,LDE,E,DEL(M1+1),IP) >*/
    sol_(nm1, lde, &e[e_offset], &del[*m1 + 1], &ip[1]);
/*<       DO I=M1,1,-1 >*/
    for (i__ = *m1; i__ >= 1; --i__) {
/*<          DEL(I)=(DEL(I)+DEL(M2+I))/FAC1 >*/
	del[i__] = (del[i__] + del[*m2 + i__]) / *fac1;
/*<       END DO >*/
    }
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    2  CONTINUE >*/
L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
/*<       CALL SOLB (N,LDE,E,MLE,MUE,DEL,IP) >*/
    solb_(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &del[1], &ip[1]);
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   12  CONTINUE >*/
L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
/*<       MM=M1/M2 >*/
    mm = *m1 / *m2;
/*<       DO J=1,M2 >*/
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM=0.D0 >*/
	sum = 0.;
/*<          DO K=MM-1,0,-1 >*/
	for (k = mm - 1; k >= 0; --k) {
/*<             JKM=J+K*M2 >*/
	    jkm = j + k * *m2;
/*<             SUM=(DEL(JKM)+SUM)/FAC1 >*/
	    sum = (del[jkm] + sum) / *fac1;
/*<             DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC) >*/
/* Computing MAX */
	    i__2 = 1, i__3 = j - *mujac;
/* Computing MIN */
	    i__5 = *nm1, i__6 = j + *mljac;
	    i__4 = min(i__5,i__6);
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
/*<                IM1=I+M1 >*/
		im1 = i__ + *m1;
/*<                DEL(IM1)=DEL(IM1)+FJAC(I+MUJAC+1-J,JKM)*SUM >*/
		del[im1] += fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1] * 
			sum;
/*<             END DO >*/
	    }
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOLB (NM1,LDE,E,MLE,MUE,DEL(M1+1),IP) >*/
    solb_(nm1, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &del[*m1 + 1], &
	    ip[1]);
/*<       DO I=M1,1,-1 >*/
    for (i__ = *m1; i__ >= 1; --i__) {
/*<          DEL(I)=(DEL(I)+DEL(M2+I))/FAC1 >*/
	del[i__] = (del[i__] + del[*m2 + i__]) / *fac1;
/*<       END DO >*/
    }
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<    7  CONTINUE >*/
L7:
/* ---  HESSENBERG OPTION */
/*<       DO MMM=N-2,1,-1 >*/
    for (mmm = *n - 2; mmm >= 1; --mmm) {
/*<          MP=N-MMM >*/
	mp = *n - mmm;
/*<          MP1=MP-1 >*/
	mp1 = mp - 1;
/*<          I=IPHES(MP) >*/
	i__ = iphes[mp];
/*<          IF (I.EQ.MP) GOTO 110 >*/
	if (i__ == mp) {
	    goto L110;
	}
/*<          ZSAFE=DEL(MP) >*/
	zsafe = del[mp];
/*<          DEL(MP)=DEL(I) >*/
	del[mp] = del[i__];
/*<          DEL(I)=ZSAFE >*/
	del[i__] = zsafe;
/*<  110     CONTINUE >*/
L110:
/*<          DO I=MP+1,N  >*/
	i__1 = *n;
	for (i__ = mp + 1; i__ <= i__1; ++i__) {
/*<             DEL(I)=DEL(I)-FJAC(I,MP1)*DEL(MP) >*/
	    del[i__] -= fjac[i__ + mp1 * fjac_dim1] * del[mp];
/*<          END DO >*/
	}
/*<       END DO >*/
    }
/*<       CALL SOLH(N,LDE,E,1,DEL,IP) >*/
    solh_(n, lde, &e[e_offset], &c__1, &del[1], &ip[1]);
/*<       DO MMM=1,N-2 >*/
    i__1 = *n - 2;
    for (mmm = 1; mmm <= i__1; ++mmm) {
/*<          MP=N-MMM >*/
	mp = *n - mmm;
/*<          MP1=MP-1 >*/
	mp1 = mp - 1;
/*<          DO I=MP+1,N  >*/
	i__4 = *n;
	for (i__ = mp + 1; i__ <= i__4; ++i__) {
/*<             DEL(I)=DEL(I)+FJAC(I,MP1)*DEL(MP) >*/
	    del[i__] += fjac[i__ + mp1 * fjac_dim1] * del[mp];
/*<          END DO >*/
	}
/*<          I=IPHES(MP) >*/
	i__ = iphes[mp];
/*<          IF (I.EQ.MP) GOTO 240 >*/
	if (i__ == mp) {
	    goto L240;
	}
/*<          ZSAFE=DEL(MP) >*/
	zsafe = del[mp];
/*<          DEL(MP)=DEL(I) >*/
	del[mp] = del[i__];
/*<          DEL(I)=ZSAFE >*/
	del[i__] = zsafe;
/*<  240     CONTINUE >*/
L240:
/*<       END DO >*/
	;
    }
/*<       RETURN >*/
    return 0;

/* ----------------------------------------------------------- */

/*<   55  CONTINUE >*/
L55:
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* slvseu_ */

