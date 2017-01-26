/* decsol.f -- translated by f2c (version 20100827).
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

/*<       SUBROUTINE DEC (N, NDIM, A, IP, IER) >*/
/* Subroutine */ int dec_(integer *n, integer *ndim, doublereal *a, integer *
	ip, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal t;
    static integer nm1, kp1;

/* VERSION REAL DOUBLE PRECISION */
/*<       INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J >*/
/*<       DOUBLE PRECISION A,T >*/
/*<       DIMENSION A(NDIM,N), IP(N) >*/
/* ----------------------------------------------------------------------- */
/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION. */
/*  INPUT.. */
/*     N = ORDER OF MATRIX. */
/*     NDIM = DECLARED DIMENSION OF ARRAY  A . */
/*     A = MATRIX TO BE TRIANGULARIZED. */
/*  OUTPUT.. */
/*     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U . */
/*     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
/*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
/*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
/*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
/*           SINGULAR AT STAGE K. */
/*  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
/*  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N). */
/*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

/*  REFERENCE.. */
/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
/*     C.A.C.M. 15 (1972), P. 274. */
/* ----------------------------------------------------------------------- */
/*<       IER = 0 >*/
    /* Parameter adjustments */
    --ip;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *ier = 0;
/*<       IP(N) = 1 >*/
    ip[*n] = 1;
/*<       IF (N .EQ. 1) GO TO 70 >*/
    if (*n == 1) {
	goto L70;
    }
/*<       NM1 = N - 1 >*/
    nm1 = *n - 1;
/*<       DO 60 K = 1,NM1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/*<         KP1 = K + 1 >*/
	kp1 = k + 1;
/*<         M = K >*/
	m = k;
/*<         DO 10 I = KP1,N >*/
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<           IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I   >*/
	    if ((d__1 = a[i__ + k * a_dim1], abs(d__1)) > (d__2 = a[m + k * 
		    a_dim1], abs(d__2))) {
		m = i__;
	    }
/*<  10     CONTINUE >*/
/* L10: */
	}
/*<         IP(K) = M >*/
	ip[k] = m;
/*<         T = A(M,K) >*/
	t = a[m + k * a_dim1];
/*<         IF (M .EQ. K) GO TO 20 >*/
	if (m == k) {
	    goto L20;
	}
/*<         IP(N) = -IP(N) >*/
	ip[*n] = -ip[*n];
/*<         A(M,K) = A(K,K) >*/
	a[m + k * a_dim1] = a[k + k * a_dim1];
/*<         A(K,K) = T >*/
	a[k + k * a_dim1] = t;
/*<  20     CONTINUE >*/
L20:
/*<         IF (T .EQ. 0.D0) GO TO 80 >*/
	if (t == 0.) {
	    goto L80;
	}
/*<         T = 1.D0/T >*/
	t = 1. / t;
/*<         DO 30 I = KP1,N >*/
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<  30       A(I,K) = -A(I,K)*T >*/
/* L30: */
	    a[i__ + k * a_dim1] = -a[i__ + k * a_dim1] * t;
	}
/*<         DO 50 J = KP1,N >*/
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
/*<           T = A(M,J) >*/
	    t = a[m + j * a_dim1];
/*<           A(M,J) = A(K,J) >*/
	    a[m + j * a_dim1] = a[k + j * a_dim1];
/*<           A(K,J) = T >*/
	    a[k + j * a_dim1] = t;
/*<           IF (T .EQ. 0.D0) GO TO 45 >*/
	    if (t == 0.) {
		goto L45;
	    }
/*<           DO 40 I = KP1,N >*/
	    i__3 = *n;
	    for (i__ = kp1; i__ <= i__3; ++i__) {
/*<  40         A(I,J) = A(I,J) + A(I,K)*T >*/
/* L40: */
		a[i__ + j * a_dim1] += a[i__ + k * a_dim1] * t;
	    }
/*<  45       CONTINUE >*/
L45:
/*<  50       CONTINUE >*/
/* L50: */
	    ;
	}
/*<  60     CONTINUE >*/
/* L60: */
    }
/*<  70   K = N >*/
L70:
    k = *n;
/*<       IF (A(N,N) .EQ. 0.D0) GO TO 80 >*/
    if (a[*n + *n * a_dim1] == 0.) {
	goto L80;
    }
/*<       RETURN >*/
    return 0;
/*<  80   IER = K >*/
L80:
    *ier = k;
/*<       IP(N) = 0 >*/
    ip[*n] = 0;
/*<       RETURN >*/
    return 0;
/* ----------------------- END OF SUBROUTINE DEC ------------------------- */
/*<       END >*/
} /* dec_ */



/*<       SUBROUTINE SOL (N, NDIM, A, B, IP) >*/
/* Subroutine */ int sol_(integer *n, integer *ndim, doublereal *a, 
	doublereal *b, integer *ip)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, m;
    static doublereal t;
    static integer kb, km1, nm1, kp1;

/* VERSION REAL DOUBLE PRECISION */
/*<       INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1 >*/
/*<       DOUBLE PRECISION A,B,T >*/
/*<       DIMENSION A(NDIM,N), B(N), IP(N) >*/
/* ----------------------------------------------------------------------- */
/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
/*  INPUT.. */
/*    N = ORDER OF MATRIX. */
/*    NDIM = DECLARED DIMENSION OF ARRAY  A . */
/*    A = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
/*    B = RIGHT HAND SIDE VECTOR. */
/*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
/*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
/*  OUTPUT.. */
/*    B = SOLUTION VECTOR, X . */
/* ----------------------------------------------------------------------- */
/*<       IF (N .EQ. 1) GO TO 50 >*/
    /* Parameter adjustments */
    --ip;
    --b;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    if (*n == 1) {
	goto L50;
    }
/*<       NM1 = N - 1 >*/
    nm1 = *n - 1;
/*<       DO 20 K = 1,NM1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/*<         KP1 = K + 1 >*/
	kp1 = k + 1;
/*<         M = IP(K) >*/
	m = ip[k];
/*<         T = B(M) >*/
	t = b[m];
/*<         B(M) = B(K) >*/
	b[m] = b[k];
/*<         B(K) = T >*/
	b[k] = t;
/*<         DO 10 I = KP1,N >*/
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<  10       B(I) = B(I) + A(I,K)*T >*/
/* L10: */
	    b[i__] += a[i__ + k * a_dim1] * t;
	}
/*<  20     CONTINUE >*/
/* L20: */
    }
/*<       DO 40 KB = 1,NM1 >*/
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
/*<         KM1 = N - KB >*/
	km1 = *n - kb;
/*<         K = KM1 + 1 >*/
	k = km1 + 1;
/*<         B(K) = B(K)/A(K,K) >*/
	b[k] /= a[k + k * a_dim1];
/*<         T = -B(K) >*/
	t = -b[k];
/*<         DO 30 I = 1,KM1 >*/
	i__2 = km1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<  30       B(I) = B(I) + A(I,K)*T >*/
/* L30: */
	    b[i__] += a[i__ + k * a_dim1] * t;
	}
/*<  40     CONTINUE >*/
/* L40: */
    }
/*<  50   B(1) = B(1)/A(1,1) >*/
L50:
    b[1] /= a[a_dim1 + 1];
/*<       RETURN >*/
    return 0;
/* ----------------------- END OF SUBROUTINE SOL ------------------------- */
/*<       END >*/
} /* sol_ */



/*<       SUBROUTINE DECH (N, NDIM, A, LB, IP, IER) >*/
/* Subroutine */ int dech_(integer *n, integer *ndim, doublereal *a, integer *
	lb, integer *ip, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal t;
    static integer na, nm1, kp1;

/* VERSION REAL DOUBLE PRECISION */
/*<       INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J,LB,NA >*/
/*<       DOUBLE PRECISION A,T >*/
/*<       DIMENSION A(NDIM,N), IP(N) >*/
/* ----------------------------------------------------------------------- */
/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A HESSENBERG */
/*  MATRIX WITH LOWER BANDWIDTH LB */
/*  INPUT.. */
/*     N = ORDER OF MATRIX A. */
/*     NDIM = DECLARED DIMENSION OF ARRAY  A . */
/*     A = MATRIX TO BE TRIANGULARIZED. */
/*     LB = LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED, LB.GE.1). */
/*  OUTPUT.. */
/*     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U . */
/*     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
/*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
/*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
/*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
/*           SINGULAR AT STAGE K. */
/*  USE  SOLH  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
/*  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N). */
/*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

/*  REFERENCE.. */
/*     THIS IS A SLIGHT MODIFICATION OF */
/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
/*     C.A.C.M. 15 (1972), P. 274. */
/* ----------------------------------------------------------------------- */
/*<       IER = 0 >*/
    /* Parameter adjustments */
    --ip;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *ier = 0;
/*<       IP(N) = 1 >*/
    ip[*n] = 1;
/*<       IF (N .EQ. 1) GO TO 70 >*/
    if (*n == 1) {
	goto L70;
    }
/*<       NM1 = N - 1 >*/
    nm1 = *n - 1;
/*<       DO 60 K = 1,NM1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/*<         KP1 = K + 1 >*/
	kp1 = k + 1;
/*<         M = K >*/
	m = k;
/*<         NA = MIN0(N,LB+K) >*/
/* Computing MIN */
	i__2 = *n, i__3 = *lb + k;
	na = min(i__2,i__3);
/*<         DO 10 I = KP1,NA >*/
	i__2 = na;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<           IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I >*/
	    if ((d__1 = a[i__ + k * a_dim1], abs(d__1)) > (d__2 = a[m + k * 
		    a_dim1], abs(d__2))) {
		m = i__;
	    }
/*<  10     CONTINUE >*/
/* L10: */
	}
/*<         IP(K) = M >*/
	ip[k] = m;
/*<         T = A(M,K) >*/
	t = a[m + k * a_dim1];
/*<         IF (M .EQ. K) GO TO 20 >*/
	if (m == k) {
	    goto L20;
	}
/*<         IP(N) = -IP(N) >*/
	ip[*n] = -ip[*n];
/*<         A(M,K) = A(K,K) >*/
	a[m + k * a_dim1] = a[k + k * a_dim1];
/*<         A(K,K) = T >*/
	a[k + k * a_dim1] = t;
/*<  20     CONTINUE >*/
L20:
/*<         IF (T .EQ. 0.D0) GO TO 80 >*/
	if (t == 0.) {
	    goto L80;
	}
/*<         T = 1.D0/T >*/
	t = 1. / t;
/*<         DO 30 I = KP1,NA >*/
	i__2 = na;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<  30       A(I,K) = -A(I,K)*T >*/
/* L30: */
	    a[i__ + k * a_dim1] = -a[i__ + k * a_dim1] * t;
	}
/*<         DO 50 J = KP1,N >*/
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
/*<           T = A(M,J) >*/
	    t = a[m + j * a_dim1];
/*<           A(M,J) = A(K,J) >*/
	    a[m + j * a_dim1] = a[k + j * a_dim1];
/*<           A(K,J) = T >*/
	    a[k + j * a_dim1] = t;
/*<           IF (T .EQ. 0.D0) GO TO 45 >*/
	    if (t == 0.) {
		goto L45;
	    }
/*<           DO 40 I = KP1,NA >*/
	    i__3 = na;
	    for (i__ = kp1; i__ <= i__3; ++i__) {
/*<  40         A(I,J) = A(I,J) + A(I,K)*T >*/
/* L40: */
		a[i__ + j * a_dim1] += a[i__ + k * a_dim1] * t;
	    }
/*<  45       CONTINUE >*/
L45:
/*<  50       CONTINUE >*/
/* L50: */
	    ;
	}
/*<  60     CONTINUE >*/
/* L60: */
    }
/*<  70   K = N >*/
L70:
    k = *n;
/*<       IF (A(N,N) .EQ. 0.D0) GO TO 80 >*/
    if (a[*n + *n * a_dim1] == 0.) {
	goto L80;
    }
/*<       RETURN >*/
    return 0;
/*<  80   IER = K >*/
L80:
    *ier = k;
/*<       IP(N) = 0 >*/
    ip[*n] = 0;
/*<       RETURN >*/
    return 0;
/* ----------------------- END OF SUBROUTINE DECH ------------------------ */
/*<       END >*/
} /* dech_ */



/*<       SUBROUTINE SOLH (N, NDIM, A, LB, B, IP) >*/
/* Subroutine */ int solh_(integer *n, integer *ndim, doublereal *a, integer *
	lb, doublereal *b, integer *ip)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k, m;
    static doublereal t;
    static integer kb, na, km1, nm1, kp1;

/* VERSION REAL DOUBLE PRECISION */
/*<       INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1,LB,NA >*/
/*<       DOUBLE PRECISION A,B,T >*/
/*<       DIMENSION A(NDIM,N), B(N), IP(N) >*/
/* ----------------------------------------------------------------------- */
/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
/*  INPUT.. */
/*    N = ORDER OF MATRIX A. */
/*    NDIM = DECLARED DIMENSION OF ARRAY  A . */
/*    A = TRIANGULARIZED MATRIX OBTAINED FROM DECH. */
/*    LB = LOWER BANDWIDTH OF A. */
/*    B = RIGHT HAND SIDE VECTOR. */
/*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
/*  DO NOT USE IF DECH HAS SET IER .NE. 0. */
/*  OUTPUT.. */
/*    B = SOLUTION VECTOR, X . */
/* ----------------------------------------------------------------------- */
/*<       IF (N .EQ. 1) GO TO 50 >*/
    /* Parameter adjustments */
    --ip;
    --b;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    if (*n == 1) {
	goto L50;
    }
/*<       NM1 = N - 1 >*/
    nm1 = *n - 1;
/*<       DO 20 K = 1,NM1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/*<         KP1 = K + 1 >*/
	kp1 = k + 1;
/*<         M = IP(K) >*/
	m = ip[k];
/*<         T = B(M) >*/
	t = b[m];
/*<         B(M) = B(K) >*/
	b[m] = b[k];
/*<         B(K) = T >*/
	b[k] = t;
/*<         NA = MIN0(N,LB+K) >*/
/* Computing MIN */
	i__2 = *n, i__3 = *lb + k;
	na = min(i__2,i__3);
/*<         DO 10 I = KP1,NA >*/
	i__2 = na;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<  10       B(I) = B(I) + A(I,K)*T >*/
/* L10: */
	    b[i__] += a[i__ + k * a_dim1] * t;
	}
/*<  20     CONTINUE >*/
/* L20: */
    }
/*<       DO 40 KB = 1,NM1 >*/
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
/*<         KM1 = N - KB >*/
	km1 = *n - kb;
/*<         K = KM1 + 1 >*/
	k = km1 + 1;
/*<         B(K) = B(K)/A(K,K) >*/
	b[k] /= a[k + k * a_dim1];
/*<         T = -B(K) >*/
	t = -b[k];
/*<         DO 30 I = 1,KM1 >*/
	i__2 = km1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<  30       B(I) = B(I) + A(I,K)*T >*/
/* L30: */
	    b[i__] += a[i__ + k * a_dim1] * t;
	}
/*<  40     CONTINUE >*/
/* L40: */
    }
/*<  50   B(1) = B(1)/A(1,1) >*/
L50:
    b[1] /= a[a_dim1 + 1];
/*<       RETURN >*/
    return 0;
/* ----------------------- END OF SUBROUTINE SOLH ------------------------ */
/*<       END >*/
} /* solh_ */


/*<       SUBROUTINE DECC (N, NDIM, AR, AI, IP, IER) >*/
/* Subroutine */ int decc_(integer *n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *ip, integer *ier)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal ti, tr;
    static integer nm1, kp1;
    static doublereal den, prodi, prodr;

/* VERSION COMPLEX DOUBLE PRECISION */
/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<       INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J >*/
/*<       DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N) >*/
/* ----------------------------------------------------------------------- */
/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION */
/*  ------ MODIFICATION FOR COMPLEX MATRICES -------- */
/*  INPUT.. */
/*     N = ORDER OF MATRIX. */
/*     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI . */
/*     (AR, AI) = MATRIX TO BE TRIANGULARIZED. */
/*  OUTPUT.. */
/*     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART. */
/*     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART. */
/*     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
/*                                                    REAL PART. */
/*     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
/*                                                    IMAGINARY PART. */
/*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
/*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
/*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
/*           SINGULAR AT STAGE K. */
/*  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
/*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

/*  REFERENCE.. */
/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
/*     C.A.C.M. 15 (1972), P. 274. */
/* ----------------------------------------------------------------------- */
/*<       IER = 0 >*/
    /* Parameter adjustments */
    --ip;
    ai_dim1 = *ndim;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *ndim;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;

    /* Function Body */
    *ier = 0;
/*<       IP(N) = 1 >*/
    ip[*n] = 1;
/*<       IF (N .EQ. 1) GO TO 70 >*/
    if (*n == 1) {
	goto L70;
    }
/*<       NM1 = N - 1 >*/
    nm1 = *n - 1;
/*<       DO 60 K = 1,NM1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/*<         KP1 = K + 1 >*/
	kp1 = k + 1;
/*<         M = K >*/
	m = k;
/*<         DO 10 I = KP1,N >*/
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<        >*/
	    if ((d__1 = ar[i__ + k * ar_dim1], abs(d__1)) + (d__2 = ai[i__ + 
		    k * ai_dim1], abs(d__2)) > (d__3 = ar[m + k * ar_dim1], 
		    abs(d__3)) + (d__4 = ai[m + k * ai_dim1], abs(d__4))) {
		m = i__;
	    }
/*<  10     CONTINUE >*/
/* L10: */
	}
/*<         IP(K) = M >*/
	ip[k] = m;
/*<         TR = AR(M,K) >*/
	tr = ar[m + k * ar_dim1];
/*<         TI = AI(M,K) >*/
	ti = ai[m + k * ai_dim1];
/*<         IF (M .EQ. K) GO TO 20 >*/
	if (m == k) {
	    goto L20;
	}
/*<         IP(N) = -IP(N) >*/
	ip[*n] = -ip[*n];
/*<         AR(M,K) = AR(K,K) >*/
	ar[m + k * ar_dim1] = ar[k + k * ar_dim1];
/*<         AI(M,K) = AI(K,K) >*/
	ai[m + k * ai_dim1] = ai[k + k * ai_dim1];
/*<         AR(K,K) = TR >*/
	ar[k + k * ar_dim1] = tr;
/*<         AI(K,K) = TI >*/
	ai[k + k * ai_dim1] = ti;
/*<  20     CONTINUE >*/
L20:
/*<         IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 80 >*/
	if (abs(tr) + abs(ti) == 0.) {
	    goto L80;
	}
/*<         DEN=TR*TR+TI*TI >*/
	den = tr * tr + ti * ti;
/*<         TR=TR/DEN >*/
	tr /= den;
/*<         TI=-TI/DEN >*/
	ti = -ti / den;
/*<         DO 30 I = KP1,N >*/
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<           PRODR=AR(I,K)*TR-AI(I,K)*TI >*/
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
/*<           PRODI=AI(I,K)*TR+AR(I,K)*TI >*/
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
/*<           AR(I,K)=-PRODR >*/
	    ar[i__ + k * ar_dim1] = -prodr;
/*<           AI(I,K)=-PRODI >*/
	    ai[i__ + k * ai_dim1] = -prodi;
/*<  30       CONTINUE >*/
/* L30: */
	}
/*<         DO 50 J = KP1,N >*/
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
/*<           TR = AR(M,J) >*/
	    tr = ar[m + j * ar_dim1];
/*<           TI = AI(M,J) >*/
	    ti = ai[m + j * ai_dim1];
/*<           AR(M,J) = AR(K,J) >*/
	    ar[m + j * ar_dim1] = ar[k + j * ar_dim1];
/*<           AI(M,J) = AI(K,J) >*/
	    ai[m + j * ai_dim1] = ai[k + j * ai_dim1];
/*<           AR(K,J) = TR >*/
	    ar[k + j * ar_dim1] = tr;
/*<           AI(K,J) = TI >*/
	    ai[k + j * ai_dim1] = ti;
/*<           IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 48 >*/
	    if (abs(tr) + abs(ti) == 0.) {
		goto L48;
	    }
/*<           IF (TI .EQ. 0.D0) THEN >*/
	    if (ti == 0.) {
/*<             DO 40 I = KP1,N >*/
		i__3 = *n;
		for (i__ = kp1; i__ <= i__3; ++i__) {
/*<             PRODR=AR(I,K)*TR >*/
		    prodr = ar[i__ + k * ar_dim1] * tr;
/*<             PRODI=AI(I,K)*TR >*/
		    prodi = ai[i__ + k * ai_dim1] * tr;
/*<             AR(I,J) = AR(I,J) + PRODR >*/
		    ar[i__ + j * ar_dim1] += prodr;
/*<             AI(I,J) = AI(I,J) + PRODI >*/
		    ai[i__ + j * ai_dim1] += prodi;
/*<  40         CONTINUE >*/
/* L40: */
		}
/*<             GO TO 48 >*/
		goto L48;
/*<           END IF >*/
	    }
/*<           IF (TR .EQ. 0.D0) THEN >*/
	    if (tr == 0.) {
/*<             DO 45 I = KP1,N >*/
		i__3 = *n;
		for (i__ = kp1; i__ <= i__3; ++i__) {
/*<             PRODR=-AI(I,K)*TI >*/
		    prodr = -ai[i__ + k * ai_dim1] * ti;
/*<             PRODI=AR(I,K)*TI >*/
		    prodi = ar[i__ + k * ar_dim1] * ti;
/*<             AR(I,J) = AR(I,J) + PRODR >*/
		    ar[i__ + j * ar_dim1] += prodr;
/*<             AI(I,J) = AI(I,J) + PRODI >*/
		    ai[i__ + j * ai_dim1] += prodi;
/*<  45         CONTINUE >*/
/* L45: */
		}
/*<             GO TO 48 >*/
		goto L48;
/*<           END IF >*/
	    }
/*<           DO 47 I = KP1,N >*/
	    i__3 = *n;
	    for (i__ = kp1; i__ <= i__3; ++i__) {
/*<             PRODR=AR(I,K)*TR-AI(I,K)*TI >*/
		prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * 
			ti;
/*<             PRODI=AI(I,K)*TR+AR(I,K)*TI >*/
		prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * 
			ti;
/*<             AR(I,J) = AR(I,J) + PRODR >*/
		ar[i__ + j * ar_dim1] += prodr;
/*<             AI(I,J) = AI(I,J) + PRODI >*/
		ai[i__ + j * ai_dim1] += prodi;
/*<  47         CONTINUE >*/
/* L47: */
	    }
/*<  48       CONTINUE >*/
L48:
/*<  50       CONTINUE >*/
/* L50: */
	    ;
	}
/*<  60     CONTINUE >*/
/* L60: */
    }
/*<  70   K = N >*/
L70:
    k = *n;
/*<       IF (DABS(AR(N,N))+DABS(AI(N,N)) .EQ. 0.D0) GO TO 80 >*/
    if ((d__1 = ar[*n + *n * ar_dim1], abs(d__1)) + (d__2 = ai[*n + *n * 
	    ai_dim1], abs(d__2)) == 0.) {
	goto L80;
    }
/*<       RETURN >*/
    return 0;
/*<  80   IER = K >*/
L80:
    *ier = k;
/*<       IP(N) = 0 >*/
    ip[*n] = 0;
/*<       RETURN >*/
    return 0;
/* ----------------------- END OF SUBROUTINE DECC ------------------------ */
/*<       END >*/
} /* decc_ */



/*<       SUBROUTINE SOLC (N, NDIM, AR, AI, BR, BI, IP) >*/
/* Subroutine */ int solc_(integer *n, integer *ndim, doublereal *ar, 
	doublereal *ai, doublereal *br, doublereal *bi, integer *ip)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, m, kb;
    static doublereal ti, tr;
    static integer km1, nm1, kp1;
    static doublereal den, prodi, prodr;

/* VERSION COMPLEX DOUBLE PRECISION */
/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<       INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1 >*/
/*<       DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N) >*/
/* ----------------------------------------------------------------------- */
/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
/*  INPUT.. */
/*    N = ORDER OF MATRIX. */
/*    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI. */
/*    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
/*    (BR,BI) = RIGHT HAND SIDE VECTOR. */
/*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
/*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
/*  OUTPUT.. */
/*    (BR,BI) = SOLUTION VECTOR, X . */
/* ----------------------------------------------------------------------- */
/*<       IF (N .EQ. 1) GO TO 50 >*/
    /* Parameter adjustments */
    --ip;
    --bi;
    --br;
    ai_dim1 = *ndim;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *ndim;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;

    /* Function Body */
    if (*n == 1) {
	goto L50;
    }
/*<       NM1 = N - 1 >*/
    nm1 = *n - 1;
/*<       DO 20 K = 1,NM1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/*<         KP1 = K + 1 >*/
	kp1 = k + 1;
/*<         M = IP(K) >*/
	m = ip[k];
/*<         TR = BR(M) >*/
	tr = br[m];
/*<         TI = BI(M) >*/
	ti = bi[m];
/*<         BR(M) = BR(K) >*/
	br[m] = br[k];
/*<         BI(M) = BI(K) >*/
	bi[m] = bi[k];
/*<         BR(K) = TR >*/
	br[k] = tr;
/*<         BI(K) = TI >*/
	bi[k] = ti;
/*<         DO 10 I = KP1,N >*/
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<           PRODR=AR(I,K)*TR-AI(I,K)*TI >*/
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
/*<           PRODI=AI(I,K)*TR+AR(I,K)*TI >*/
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
/*<           BR(I) = BR(I) + PRODR >*/
	    br[i__] += prodr;
/*<           BI(I) = BI(I) + PRODI >*/
	    bi[i__] += prodi;
/*<  10       CONTINUE >*/
/* L10: */
	}
/*<  20     CONTINUE >*/
/* L20: */
    }
/*<       DO 40 KB = 1,NM1 >*/
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
/*<         KM1 = N - KB >*/
	km1 = *n - kb;
/*<         K = KM1 + 1 >*/
	k = km1 + 1;
/*<         DEN=AR(K,K)*AR(K,K)+AI(K,K)*AI(K,K) >*/
	den = ar[k + k * ar_dim1] * ar[k + k * ar_dim1] + ai[k + k * ai_dim1] 
		* ai[k + k * ai_dim1];
/*<         PRODR=BR(K)*AR(K,K)+BI(K)*AI(K,K) >*/
	prodr = br[k] * ar[k + k * ar_dim1] + bi[k] * ai[k + k * ai_dim1];
/*<         PRODI=BI(K)*AR(K,K)-BR(K)*AI(K,K) >*/
	prodi = bi[k] * ar[k + k * ar_dim1] - br[k] * ai[k + k * ai_dim1];
/*<         BR(K)=PRODR/DEN >*/
	br[k] = prodr / den;
/*<         BI(K)=PRODI/DEN >*/
	bi[k] = prodi / den;
/*<         TR = -BR(K) >*/
	tr = -br[k];
/*<         TI = -BI(K) >*/
	ti = -bi[k];
/*<         DO 30 I = 1,KM1 >*/
	i__2 = km1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           PRODR=AR(I,K)*TR-AI(I,K)*TI >*/
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
/*<           PRODI=AI(I,K)*TR+AR(I,K)*TI >*/
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
/*<           BR(I) = BR(I) + PRODR >*/
	    br[i__] += prodr;
/*<           BI(I) = BI(I) + PRODI >*/
	    bi[i__] += prodi;
/*<  30       CONTINUE >*/
/* L30: */
	}
/*<  40     CONTINUE >*/
/* L40: */
    }
/*<  50     CONTINUE >*/
L50:
/*<         DEN=AR(1,1)*AR(1,1)+AI(1,1)*AI(1,1) >*/
    den = ar[ar_dim1 + 1] * ar[ar_dim1 + 1] + ai[ai_dim1 + 1] * ai[ai_dim1 + 
	    1];
/*<         PRODR=BR(1)*AR(1,1)+BI(1)*AI(1,1) >*/
    prodr = br[1] * ar[ar_dim1 + 1] + bi[1] * ai[ai_dim1 + 1];
/*<         PRODI=BI(1)*AR(1,1)-BR(1)*AI(1,1) >*/
    prodi = bi[1] * ar[ar_dim1 + 1] - br[1] * ai[ai_dim1 + 1];
/*<         BR(1)=PRODR/DEN >*/
    br[1] = prodr / den;
/*<         BI(1)=PRODI/DEN >*/
    bi[1] = prodi / den;
/*<       RETURN >*/
    return 0;
/* ----------------------- END OF SUBROUTINE SOLC ------------------------ */
/*<       END   >*/
} /* solc_ */



/*<       SUBROUTINE DECHC (N, NDIM, AR, AI, LB, IP, IER) >*/
/* Subroutine */ int dechc_(integer *n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *lb, integer *ip, integer *ier)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, j, k, m, na;
    static doublereal ti, tr;
    static integer nm1, kp1;
    static doublereal den, prodi, prodr;

/* VERSION COMPLEX DOUBLE PRECISION */
/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<       INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J >*/
/*<       DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N) >*/
/* ----------------------------------------------------------------------- */
/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION */
/*  ------ MODIFICATION FOR COMPLEX MATRICES -------- */
/*  INPUT.. */
/*     N = ORDER OF MATRIX. */
/*     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI . */
/*     (AR, AI) = MATRIX TO BE TRIANGULARIZED. */
/*  OUTPUT.. */
/*     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART. */
/*     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART. */
/*     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
/*                                                    REAL PART. */
/*     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
/*                                                    IMAGINARY PART. */
/*     LB = LOWER BANDWIDTH OF A (DIAGONAL NOT COUNTED), LB.GE.1. */
/*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
/*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
/*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
/*           SINGULAR AT STAGE K. */
/*  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
/*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

/*  REFERENCE.. */
/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
/*     C.A.C.M. 15 (1972), P. 274. */
/* ----------------------------------------------------------------------- */
/*<       IER = 0 >*/
    /* Parameter adjustments */
    --ip;
    ai_dim1 = *ndim;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *ndim;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;

    /* Function Body */
    *ier = 0;
/*<       IP(N) = 1 >*/
    ip[*n] = 1;
/*<       IF (LB .EQ. 0) GO TO 70 >*/
    if (*lb == 0) {
	goto L70;
    }
/*<       IF (N .EQ. 1) GO TO 70 >*/
    if (*n == 1) {
	goto L70;
    }
/*<       NM1 = N - 1 >*/
    nm1 = *n - 1;
/*<       DO 60 K = 1,NM1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/*<         KP1 = K + 1 >*/
	kp1 = k + 1;
/*<         M = K  >*/
	m = k;
/*<         NA = MIN0(N,LB+K) >*/
/* Computing MIN */
	i__2 = *n, i__3 = *lb + k;
	na = min(i__2,i__3);
/*<         DO 10 I = KP1,NA >*/
	i__2 = na;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<        >*/
	    if ((d__1 = ar[i__ + k * ar_dim1], abs(d__1)) + (d__2 = ai[i__ + 
		    k * ai_dim1], abs(d__2)) > (d__3 = ar[m + k * ar_dim1], 
		    abs(d__3)) + (d__4 = ai[m + k * ai_dim1], abs(d__4))) {
		m = i__;
	    }
/*<  10     CONTINUE >*/
/* L10: */
	}
/*<         IP(K) = M >*/
	ip[k] = m;
/*<         TR = AR(M,K) >*/
	tr = ar[m + k * ar_dim1];
/*<         TI = AI(M,K) >*/
	ti = ai[m + k * ai_dim1];
/*<         IF (M .EQ. K) GO TO 20 >*/
	if (m == k) {
	    goto L20;
	}
/*<         IP(N) = -IP(N) >*/
	ip[*n] = -ip[*n];
/*<         AR(M,K) = AR(K,K) >*/
	ar[m + k * ar_dim1] = ar[k + k * ar_dim1];
/*<         AI(M,K) = AI(K,K) >*/
	ai[m + k * ai_dim1] = ai[k + k * ai_dim1];
/*<         AR(K,K) = TR >*/
	ar[k + k * ar_dim1] = tr;
/*<         AI(K,K) = TI >*/
	ai[k + k * ai_dim1] = ti;
/*<  20     CONTINUE >*/
L20:
/*<         IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 80 >*/
	if (abs(tr) + abs(ti) == 0.) {
	    goto L80;
	}
/*<         DEN=TR*TR+TI*TI >*/
	den = tr * tr + ti * ti;
/*<         TR=TR/DEN >*/
	tr /= den;
/*<         TI=-TI/DEN >*/
	ti = -ti / den;
/*<         DO 30 I = KP1,NA >*/
	i__2 = na;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<           PRODR=AR(I,K)*TR-AI(I,K)*TI >*/
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
/*<           PRODI=AI(I,K)*TR+AR(I,K)*TI >*/
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
/*<           AR(I,K)=-PRODR >*/
	    ar[i__ + k * ar_dim1] = -prodr;
/*<           AI(I,K)=-PRODI >*/
	    ai[i__ + k * ai_dim1] = -prodi;
/*<  30       CONTINUE >*/
/* L30: */
	}
/*<         DO 50 J = KP1,N >*/
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
/*<           TR = AR(M,J) >*/
	    tr = ar[m + j * ar_dim1];
/*<           TI = AI(M,J) >*/
	    ti = ai[m + j * ai_dim1];
/*<           AR(M,J) = AR(K,J) >*/
	    ar[m + j * ar_dim1] = ar[k + j * ar_dim1];
/*<           AI(M,J) = AI(K,J) >*/
	    ai[m + j * ai_dim1] = ai[k + j * ai_dim1];
/*<           AR(K,J) = TR >*/
	    ar[k + j * ar_dim1] = tr;
/*<           AI(K,J) = TI >*/
	    ai[k + j * ai_dim1] = ti;
/*<           IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 48 >*/
	    if (abs(tr) + abs(ti) == 0.) {
		goto L48;
	    }
/*<           IF (TI .EQ. 0.D0) THEN >*/
	    if (ti == 0.) {
/*<             DO 40 I = KP1,NA >*/
		i__3 = na;
		for (i__ = kp1; i__ <= i__3; ++i__) {
/*<             PRODR=AR(I,K)*TR >*/
		    prodr = ar[i__ + k * ar_dim1] * tr;
/*<             PRODI=AI(I,K)*TR >*/
		    prodi = ai[i__ + k * ai_dim1] * tr;
/*<             AR(I,J) = AR(I,J) + PRODR >*/
		    ar[i__ + j * ar_dim1] += prodr;
/*<             AI(I,J) = AI(I,J) + PRODI >*/
		    ai[i__ + j * ai_dim1] += prodi;
/*<  40         CONTINUE >*/
/* L40: */
		}
/*<             GO TO 48 >*/
		goto L48;
/*<           END IF >*/
	    }
/*<           IF (TR .EQ. 0.D0) THEN >*/
	    if (tr == 0.) {
/*<             DO 45 I = KP1,NA >*/
		i__3 = na;
		for (i__ = kp1; i__ <= i__3; ++i__) {
/*<             PRODR=-AI(I,K)*TI >*/
		    prodr = -ai[i__ + k * ai_dim1] * ti;
/*<             PRODI=AR(I,K)*TI >*/
		    prodi = ar[i__ + k * ar_dim1] * ti;
/*<             AR(I,J) = AR(I,J) + PRODR >*/
		    ar[i__ + j * ar_dim1] += prodr;
/*<             AI(I,J) = AI(I,J) + PRODI >*/
		    ai[i__ + j * ai_dim1] += prodi;
/*<  45         CONTINUE >*/
/* L45: */
		}
/*<             GO TO 48 >*/
		goto L48;
/*<           END IF >*/
	    }
/*<           DO 47 I = KP1,NA >*/
	    i__3 = na;
	    for (i__ = kp1; i__ <= i__3; ++i__) {
/*<             PRODR=AR(I,K)*TR-AI(I,K)*TI >*/
		prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * 
			ti;
/*<             PRODI=AI(I,K)*TR+AR(I,K)*TI >*/
		prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * 
			ti;
/*<             AR(I,J) = AR(I,J) + PRODR >*/
		ar[i__ + j * ar_dim1] += prodr;
/*<             AI(I,J) = AI(I,J) + PRODI >*/
		ai[i__ + j * ai_dim1] += prodi;
/*<  47         CONTINUE >*/
/* L47: */
	    }
/*<  48       CONTINUE >*/
L48:
/*<  50       CONTINUE >*/
/* L50: */
	    ;
	}
/*<  60     CONTINUE >*/
/* L60: */
    }
/*<  70   K = N >*/
L70:
    k = *n;
/*<       IF (DABS(AR(N,N))+DABS(AI(N,N)) .EQ. 0.D0) GO TO 80 >*/
    if ((d__1 = ar[*n + *n * ar_dim1], abs(d__1)) + (d__2 = ai[*n + *n * 
	    ai_dim1], abs(d__2)) == 0.) {
	goto L80;
    }
/*<       RETURN >*/
    return 0;
/*<  80   IER = K >*/
L80:
    *ier = k;
/*<       IP(N) = 0 >*/
    ip[*n] = 0;
/*<       RETURN >*/
    return 0;
/* ----------------------- END OF SUBROUTINE DECHC ----------------------- */
/*<       END >*/
} /* dechc_ */



/*<       SUBROUTINE SOLHC (N, NDIM, AR, AI, LB, BR, BI, IP) >*/
/* Subroutine */ int solhc_(integer *n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *lb, doublereal *br, doublereal *bi, integer *
	ip)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, m, kb;
    static doublereal ti, tr;
    static integer km1, nm1, kp1;
    static doublereal den, prodi, prodr;

/* VERSION COMPLEX DOUBLE PRECISION */
/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<       INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1 >*/
/*<       DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N) >*/
/* ----------------------------------------------------------------------- */
/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
/*  INPUT.. */
/*    N = ORDER OF MATRIX. */
/*    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI. */
/*    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
/*    (BR,BI) = RIGHT HAND SIDE VECTOR. */
/*    LB = LOWER BANDWIDTH OF A. */
/*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
/*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
/*  OUTPUT.. */
/*    (BR,BI) = SOLUTION VECTOR, X . */
/* ----------------------------------------------------------------------- */
/*<       IF (N .EQ. 1) GO TO 50 >*/
    /* Parameter adjustments */
    --ip;
    --bi;
    --br;
    ai_dim1 = *ndim;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *ndim;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;

    /* Function Body */
    if (*n == 1) {
	goto L50;
    }
/*<       NM1 = N - 1 >*/
    nm1 = *n - 1;
/*<       IF (LB .EQ. 0) GO TO 25 >*/
    if (*lb == 0) {
	goto L25;
    }
/*<       DO 20 K = 1,NM1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/*<         KP1 = K + 1 >*/
	kp1 = k + 1;
/*<         M = IP(K) >*/
	m = ip[k];
/*<         TR = BR(M) >*/
	tr = br[m];
/*<         TI = BI(M) >*/
	ti = bi[m];
/*<         BR(M) = BR(K) >*/
	br[m] = br[k];
/*<         BI(M) = BI(K) >*/
	bi[m] = bi[k];
/*<         BR(K) = TR >*/
	br[k] = tr;
/*<         BI(K) = TI >*/
	bi[k] = ti;
/*<         DO 10 I = KP1,MIN0(N,LB+K) >*/
/* Computing MIN */
	i__3 = *n, i__4 = *lb + k;
	i__2 = min(i__3,i__4);
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<           PRODR=AR(I,K)*TR-AI(I,K)*TI >*/
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
/*<           PRODI=AI(I,K)*TR+AR(I,K)*TI >*/
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
/*<           BR(I) = BR(I) + PRODR >*/
	    br[i__] += prodr;
/*<           BI(I) = BI(I) + PRODI >*/
	    bi[i__] += prodi;
/*<  10       CONTINUE >*/
/* L10: */
	}
/*<  20     CONTINUE >*/
/* L20: */
    }
/*<  25     CONTINUE >*/
L25:
/*<       DO 40 KB = 1,NM1 >*/
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
/*<         KM1 = N - KB >*/
	km1 = *n - kb;
/*<         K = KM1 + 1 >*/
	k = km1 + 1;
/*<         DEN=AR(K,K)*AR(K,K)+AI(K,K)*AI(K,K) >*/
	den = ar[k + k * ar_dim1] * ar[k + k * ar_dim1] + ai[k + k * ai_dim1] 
		* ai[k + k * ai_dim1];
/*<         PRODR=BR(K)*AR(K,K)+BI(K)*AI(K,K) >*/
	prodr = br[k] * ar[k + k * ar_dim1] + bi[k] * ai[k + k * ai_dim1];
/*<         PRODI=BI(K)*AR(K,K)-BR(K)*AI(K,K) >*/
	prodi = bi[k] * ar[k + k * ar_dim1] - br[k] * ai[k + k * ai_dim1];
/*<         BR(K)=PRODR/DEN >*/
	br[k] = prodr / den;
/*<         BI(K)=PRODI/DEN >*/
	bi[k] = prodi / den;
/*<         TR = -BR(K) >*/
	tr = -br[k];
/*<         TI = -BI(K) >*/
	ti = -bi[k];
/*<         DO 30 I = 1,KM1 >*/
	i__2 = km1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           PRODR=AR(I,K)*TR-AI(I,K)*TI >*/
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
/*<           PRODI=AI(I,K)*TR+AR(I,K)*TI >*/
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
/*<           BR(I) = BR(I) + PRODR >*/
	    br[i__] += prodr;
/*<           BI(I) = BI(I) + PRODI >*/
	    bi[i__] += prodi;
/*<  30       CONTINUE >*/
/* L30: */
	}
/*<  40     CONTINUE >*/
/* L40: */
    }
/*<  50     CONTINUE >*/
L50:
/*<         DEN=AR(1,1)*AR(1,1)+AI(1,1)*AI(1,1) >*/
    den = ar[ar_dim1 + 1] * ar[ar_dim1 + 1] + ai[ai_dim1 + 1] * ai[ai_dim1 + 
	    1];
/*<         PRODR=BR(1)*AR(1,1)+BI(1)*AI(1,1) >*/
    prodr = br[1] * ar[ar_dim1 + 1] + bi[1] * ai[ai_dim1 + 1];
/*<         PRODI=BI(1)*AR(1,1)-BR(1)*AI(1,1) >*/
    prodi = bi[1] * ar[ar_dim1 + 1] - br[1] * ai[ai_dim1 + 1];
/*<         BR(1)=PRODR/DEN >*/
    br[1] = prodr / den;
/*<         BI(1)=PRODI/DEN >*/
    bi[1] = prodi / den;
/*<       RETURN >*/
    return 0;
/* ----------------------- END OF SUBROUTINE SOLHC ----------------------- */
/*<       END   >*/
} /* solhc_ */


/*<       SUBROUTINE DECB (N, NDIM, A, ML, MU, IP, IER) >*/
/* Subroutine */ int decb_(integer *n, integer *ndim, doublereal *a, integer *
	ml, integer *mu, integer *ip, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal t;
    static integer md, jk, mm, ju, md1, nm1, kp1, mdl, ijk;

/*<       REAL*8 A,T >*/
/*<       DIMENSION A(NDIM,N), IP(N) >*/
/* ----------------------------------------------------------------------- */
/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED */
/*  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU */
/*  INPUT.. */
/*     N       ORDER OF THE ORIGINAL MATRIX A. */
/*     NDIM    DECLARED DIMENSION OF ARRAY  A. */
/*     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS */
/*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND */
/*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
/*                ML+1 THROUGH 2*ML+MU+1 OF  A. */
/*     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*  OUTPUT.. */
/*     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
/*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
/*     IP      INDEX VECTOR OF PIVOT INDICES. */
/*     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O . */
/*     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE */
/*                SINGULAR AT STAGE K. */
/*  USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
/*  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1. */
/*  IF IP(N)=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO. */

/*  REFERENCE.. */
/*     THIS IS A MODIFICATION OF */
/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
/*     C.A.C.M. 15 (1972), P. 274. */
/* ----------------------------------------------------------------------- */
/*<       IER = 0 >*/
    /* Parameter adjustments */
    --ip;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *ier = 0;
/*<       IP(N) = 1  >*/
    ip[*n] = 1;
/*<       MD = ML + MU + 1 >*/
    md = *ml + *mu + 1;
/*<       MD1 = MD + 1 >*/
    md1 = md + 1;
/*<       JU = 0 >*/
    ju = 0;
/*<       IF (ML .EQ. 0) GO TO 70 >*/
    if (*ml == 0) {
	goto L70;
    }
/*<       IF (N .EQ. 1) GO TO 70 >*/
    if (*n == 1) {
	goto L70;
    }
/*<       IF (N .LT. MU+2) GO TO 7 >*/
    if (*n < *mu + 2) {
	goto L7;
    }
/*<       DO 5 J = MU+2,N >*/
    i__1 = *n;
    for (j = *mu + 2; j <= i__1; ++j) {
/*<       DO 5 I = 1,ML >*/
	i__2 = *ml;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<   5   A(I,J) = 0.D0 >*/
/* L5: */
	    a[i__ + j * a_dim1] = 0.;
	}
    }
/*<   7   NM1 = N - 1 >*/
L7:
    nm1 = *n - 1;
/*<       DO 60 K = 1,NM1 >*/
    i__2 = nm1;
    for (k = 1; k <= i__2; ++k) {
/*<         KP1 = K + 1 >*/
	kp1 = k + 1;
/*<         M = MD >*/
	m = md;
/*<         MDL = MIN(ML,N-K) + MD >*/
/* Computing MIN */
	i__1 = *ml, i__3 = *n - k;
	mdl = min(i__1,i__3) + md;
/*<         DO 10 I = MD1,MDL >*/
	i__1 = mdl;
	for (i__ = md1; i__ <= i__1; ++i__) {
/*<           IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I >*/
	    if ((d__1 = a[i__ + k * a_dim1], abs(d__1)) > (d__2 = a[m + k * 
		    a_dim1], abs(d__2))) {
		m = i__;
	    }
/*<  10     CONTINUE >*/
/* L10: */
	}
/*<         IP(K) = M + K - MD >*/
	ip[k] = m + k - md;
/*<         T = A(M,K) >*/
	t = a[m + k * a_dim1];
/*<         IF (M .EQ. MD) GO TO 20 >*/
	if (m == md) {
	    goto L20;
	}
/*<         IP(N) = -IP(N) >*/
	ip[*n] = -ip[*n];
/*<         A(M,K) = A(MD,K) >*/
	a[m + k * a_dim1] = a[md + k * a_dim1];
/*<         A(MD,K) = T >*/
	a[md + k * a_dim1] = t;
/*<  20     CONTINUE >*/
L20:
/*<         IF (T .EQ. 0.D0) GO TO 80 >*/
	if (t == 0.) {
	    goto L80;
	}
/*<         T = 1.D0/T >*/
	t = 1. / t;
/*<         DO 30 I = MD1,MDL >*/
	i__1 = mdl;
	for (i__ = md1; i__ <= i__1; ++i__) {
/*<  30       A(I,K) = -A(I,K)*T  >*/
/* L30: */
	    a[i__ + k * a_dim1] = -a[i__ + k * a_dim1] * t;
	}
/*<         JU = MIN0(MAX0(JU,MU+IP(K)),N) >*/
/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ip[k];
	i__1 = max(i__3,i__4);
	ju = min(i__1,*n);
/*<         MM = MD >*/
	mm = md;
/*<         IF (JU .LT. KP1) GO TO 55 >*/
	if (ju < kp1) {
	    goto L55;
	}
/*<         DO 50 J = KP1,JU >*/
	i__1 = ju;
	for (j = kp1; j <= i__1; ++j) {
/*<           M = M - 1 >*/
	    --m;
/*<           MM = MM - 1 >*/
	    --mm;
/*<           T = A(M,J)  >*/
	    t = a[m + j * a_dim1];
/*<           IF (M .EQ. MM) GO TO 35 >*/
	    if (m == mm) {
		goto L35;
	    }
/*<           A(M,J) = A(MM,J) >*/
	    a[m + j * a_dim1] = a[mm + j * a_dim1];
/*<           A(MM,J) = T >*/
	    a[mm + j * a_dim1] = t;
/*<  35       CONTINUE >*/
L35:
/*<           IF (T .EQ. 0.D0) GO TO 45 >*/
	    if (t == 0.) {
		goto L45;
	    }
/*<           JK = J - K >*/
	    jk = j - k;
/*<           DO 40 I = MD1,MDL >*/
	    i__3 = mdl;
	    for (i__ = md1; i__ <= i__3; ++i__) {
/*<             IJK = I - JK >*/
		ijk = i__ - jk;
/*<  40         A(IJK,J) = A(IJK,J) + A(I,K)*T >*/
/* L40: */
		a[ijk + j * a_dim1] += a[i__ + k * a_dim1] * t;
	    }
/*<  45       CONTINUE >*/
L45:
/*<  50       CONTINUE >*/
/* L50: */
	    ;
	}
/*<  55     CONTINUE >*/
L55:
/*<  60     CONTINUE >*/
/* L60: */
	;
    }
/*<  70   K = N >*/
L70:
    k = *n;
/*<       IF (A(MD,N) .EQ. 0.D0) GO TO 80 >*/
    if (a[md + *n * a_dim1] == 0.) {
	goto L80;
    }
/*<       RETURN >*/
    return 0;
/*<  80   IER = K >*/
L80:
    *ier = k;
/*<       IP(N) = 0 >*/
    ip[*n] = 0;
/*<       RETURN >*/
    return 0;
/* ----------------------- END OF SUBROUTINE DECB ------------------------ */
/*<       END >*/
} /* decb_ */



/*<       SUBROUTINE SOLB (N, NDIM, A, ML, MU, B, IP) >*/
/* Subroutine */ int solb_(integer *n, integer *ndim, doublereal *a, integer *
	ml, integer *mu, doublereal *b, integer *ip)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k, m;
    static doublereal t;
    static integer kb, md, lm, md1, nm1, imd, kmd, mdl, mdm;

/*<       REAL*8 A,B,T >*/
/*<       DIMENSION A(NDIM,N), B(N), IP(N) >*/
/* ----------------------------------------------------------------------- */
/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
/*  INPUT.. */
/*    N      ORDER OF MATRIX A. */
/*    NDIM   DECLARED DIMENSION OF ARRAY  A . */
/*    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB. */
/*    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*    B      RIGHT HAND SIDE VECTOR. */
/*    IP     PIVOT VECTOR OBTAINED FROM DECB. */
/*  DO NOT USE IF DECB HAS SET IER .NE. 0. */
/*  OUTPUT.. */
/*    B      SOLUTION VECTOR, X . */
/* ----------------------------------------------------------------------- */
/*<       MD = ML + MU + 1 >*/
    /* Parameter adjustments */
    --ip;
    --b;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    md = *ml + *mu + 1;
/*<       MD1 = MD + 1 >*/
    md1 = md + 1;
/*<       MDM = MD - 1 >*/
    mdm = md - 1;
/*<       NM1 = N - 1 >*/
    nm1 = *n - 1;
/*<       IF (ML .EQ. 0) GO TO 25 >*/
    if (*ml == 0) {
	goto L25;
    }
/*<       IF (N .EQ. 1) GO TO 50 >*/
    if (*n == 1) {
	goto L50;
    }
/*<       DO 20 K = 1,NM1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/*<         M = IP(K) >*/
	m = ip[k];
/*<         T = B(M) >*/
	t = b[m];
/*<         B(M) = B(K) >*/
	b[m] = b[k];
/*<         B(K) = T >*/
	b[k] = t;
/*<         MDL = MIN(ML,N-K) + MD >*/
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	mdl = min(i__2,i__3) + md;
/*<         DO 10 I = MD1,MDL >*/
	i__2 = mdl;
	for (i__ = md1; i__ <= i__2; ++i__) {
/*<           IMD = I + K - MD >*/
	    imd = i__ + k - md;
/*<  10       B(IMD) = B(IMD) + A(I,K)*T >*/
/* L10: */
	    b[imd] += a[i__ + k * a_dim1] * t;
	}
/*<  20     CONTINUE >*/
/* L20: */
    }
/*<  25   CONTINUE >*/
L25:
/*<       DO 40 KB = 1,NM1 >*/
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
/*<         K = N + 1 - KB >*/
	k = *n + 1 - kb;
/*<         B(K) = B(K)/A(MD,K) >*/
	b[k] /= a[md + k * a_dim1];
/*<         T = -B(K)  >*/
	t = -b[k];
/*<         KMD = MD - K >*/
	kmd = md - k;
/*<         LM = MAX0(1,KMD+1) >*/
/* Computing MAX */
	i__2 = 1, i__3 = kmd + 1;
	lm = max(i__2,i__3);
/*<         DO 30 I = LM,MDM >*/
	i__2 = mdm;
	for (i__ = lm; i__ <= i__2; ++i__) {
/*<           IMD = I - KMD >*/
	    imd = i__ - kmd;
/*<  30       B(IMD) = B(IMD) + A(I,K)*T >*/
/* L30: */
	    b[imd] += a[i__ + k * a_dim1] * t;
	}
/*<  40     CONTINUE >*/
/* L40: */
    }
/*<  50   B(1) = B(1)/A(MD,1) >*/
L50:
    b[1] /= a[md + a_dim1];
/*<       RETURN >*/
    return 0;
/* ----------------------- END OF SUBROUTINE SOLB ------------------------ */
/*<       END >*/
} /* solb_ */


/*<       SUBROUTINE DECBC (N, NDIM, AR, AI, ML, MU, IP, IER) >*/
/* Subroutine */ int decbc_(integer *n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *ml, integer *mu, integer *ip, integer *ier)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, j, k, m, md, jk, mm;
    static doublereal ti;
    static integer ju;
    static doublereal tr;
    static integer md1, nm1, kp1;
    static doublereal den;
    static integer mdl, ijk;
    static doublereal prodi, prodr;

/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<       DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N) >*/
/* ----------------------------------------------------------------------- */
/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED COMPLEX */
/*  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU */
/*  INPUT.. */
/*     N       ORDER OF THE ORIGINAL MATRIX A. */
/*     NDIM    DECLARED DIMENSION OF ARRAY  A. */
/*     AR, AI     CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS */
/*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  AR (REAL */
/*                PART) AND AI (IMAGINARY PART)  AND */
/*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
/*                ML+1 THROUGH 2*ML+MU+1 OF  AR AND AI. */
/*     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*  OUTPUT.. */
/*     AR, AI  AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
/*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
/*     IP      INDEX VECTOR OF PIVOT INDICES. */
/*     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O . */
/*     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE */
/*                SINGULAR AT STAGE K. */
/*  USE  SOLBC  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
/*  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1. */
/*  IF IP(N)=O, A IS SINGULAR, SOLBC WILL DIVIDE BY ZERO. */

/*  REFERENCE.. */
/*     THIS IS A MODIFICATION OF */
/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
/*     C.A.C.M. 15 (1972), P. 274. */
/* ----------------------------------------------------------------------- */
/*<       IER = 0 >*/
    /* Parameter adjustments */
    --ip;
    ai_dim1 = *ndim;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *ndim;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;

    /* Function Body */
    *ier = 0;
/*<       IP(N) = 1  >*/
    ip[*n] = 1;
/*<       MD = ML + MU + 1 >*/
    md = *ml + *mu + 1;
/*<       MD1 = MD + 1 >*/
    md1 = md + 1;
/*<       JU = 0 >*/
    ju = 0;
/*<       IF (ML .EQ. 0) GO TO 70 >*/
    if (*ml == 0) {
	goto L70;
    }
/*<       IF (N .EQ. 1) GO TO 70 >*/
    if (*n == 1) {
	goto L70;
    }
/*<       IF (N .LT. MU+2) GO TO 7 >*/
    if (*n < *mu + 2) {
	goto L7;
    }
/*<       DO 5 J = MU+2,N >*/
    i__1 = *n;
    for (j = *mu + 2; j <= i__1; ++j) {
/*<       DO 5 I = 1,ML >*/
	i__2 = *ml;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<       AR(I,J) = 0.D0 >*/
	    ar[i__ + j * ar_dim1] = 0.;
/*<       AI(I,J) = 0.D0 >*/
	    ai[i__ + j * ai_dim1] = 0.;
/*<   5   CONTINUE >*/
/* L5: */
	}
    }
/*<   7   NM1 = N - 1 >*/
L7:
    nm1 = *n - 1;
/*<       DO 60 K = 1,NM1 >*/
    i__2 = nm1;
    for (k = 1; k <= i__2; ++k) {
/*<         KP1 = K + 1 >*/
	kp1 = k + 1;
/*<         M = MD >*/
	m = md;
/*<         MDL = MIN(ML,N-K) + MD >*/
/* Computing MIN */
	i__1 = *ml, i__3 = *n - k;
	mdl = min(i__1,i__3) + md;
/*<         DO 10 I = MD1,MDL >*/
	i__1 = mdl;
	for (i__ = md1; i__ <= i__1; ++i__) {
/*<        >*/
	    if ((d__1 = ar[i__ + k * ar_dim1], abs(d__1)) + (d__2 = ai[i__ + 
		    k * ai_dim1], abs(d__2)) > (d__3 = ar[m + k * ar_dim1], 
		    abs(d__3)) + (d__4 = ai[m + k * ai_dim1], abs(d__4))) {
		m = i__;
	    }
/*<  10     CONTINUE >*/
/* L10: */
	}
/*<         IP(K) = M + K - MD >*/
	ip[k] = m + k - md;
/*<         TR = AR(M,K) >*/
	tr = ar[m + k * ar_dim1];
/*<         TI = AI(M,K) >*/
	ti = ai[m + k * ai_dim1];
/*<         IF (M .EQ. MD) GO TO 20 >*/
	if (m == md) {
	    goto L20;
	}
/*<         IP(N) = -IP(N) >*/
	ip[*n] = -ip[*n];
/*<         AR(M,K) = AR(MD,K) >*/
	ar[m + k * ar_dim1] = ar[md + k * ar_dim1];
/*<         AI(M,K) = AI(MD,K) >*/
	ai[m + k * ai_dim1] = ai[md + k * ai_dim1];
/*<         AR(MD,K) = TR >*/
	ar[md + k * ar_dim1] = tr;
/*<         AI(MD,K) = TI >*/
	ai[md + k * ai_dim1] = ti;
/*<  20     IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 80 >*/
L20:
	if (abs(tr) + abs(ti) == 0.) {
	    goto L80;
	}
/*<         DEN=TR*TR+TI*TI >*/
	den = tr * tr + ti * ti;
/*<         TR=TR/DEN >*/
	tr /= den;
/*<         TI=-TI/DEN >*/
	ti = -ti / den;
/*<         DO 30 I = MD1,MDL >*/
	i__1 = mdl;
	for (i__ = md1; i__ <= i__1; ++i__) {
/*<           PRODR=AR(I,K)*TR-AI(I,K)*TI >*/
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
/*<           PRODI=AI(I,K)*TR+AR(I,K)*TI >*/
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
/*<           AR(I,K)=-PRODR >*/
	    ar[i__ + k * ar_dim1] = -prodr;
/*<           AI(I,K)=-PRODI >*/
	    ai[i__ + k * ai_dim1] = -prodi;
/*<  30       CONTINUE >*/
/* L30: */
	}
/*<         JU = MIN0(MAX0(JU,MU+IP(K)),N) >*/
/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ip[k];
	i__1 = max(i__3,i__4);
	ju = min(i__1,*n);
/*<         MM = MD >*/
	mm = md;
/*<         IF (JU .LT. KP1) GO TO 55 >*/
	if (ju < kp1) {
	    goto L55;
	}
/*<         DO 50 J = KP1,JU >*/
	i__1 = ju;
	for (j = kp1; j <= i__1; ++j) {
/*<           M = M - 1 >*/
	    --m;
/*<           MM = MM - 1 >*/
	    --mm;
/*<           TR = AR(M,J) >*/
	    tr = ar[m + j * ar_dim1];
/*<           TI = AI(M,J) >*/
	    ti = ai[m + j * ai_dim1];
/*<           IF (M .EQ. MM) GO TO 35 >*/
	    if (m == mm) {
		goto L35;
	    }
/*<           AR(M,J) = AR(MM,J) >*/
	    ar[m + j * ar_dim1] = ar[mm + j * ar_dim1];
/*<           AI(M,J) = AI(MM,J) >*/
	    ai[m + j * ai_dim1] = ai[mm + j * ai_dim1];
/*<           AR(MM,J) = TR >*/
	    ar[mm + j * ar_dim1] = tr;
/*<           AI(MM,J) = TI >*/
	    ai[mm + j * ai_dim1] = ti;
/*<  35       CONTINUE >*/
L35:
/*<           IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 48 >*/
	    if (abs(tr) + abs(ti) == 0.) {
		goto L48;
	    }
/*<           JK = J - K >*/
	    jk = j - k;
/*<           IF (TI .EQ. 0.D0) THEN >*/
	    if (ti == 0.) {
/*<             DO 40 I = MD1,MDL >*/
		i__3 = mdl;
		for (i__ = md1; i__ <= i__3; ++i__) {
/*<             IJK = I - JK >*/
		    ijk = i__ - jk;
/*<             PRODR=AR(I,K)*TR >*/
		    prodr = ar[i__ + k * ar_dim1] * tr;
/*<             PRODI=AI(I,K)*TR >*/
		    prodi = ai[i__ + k * ai_dim1] * tr;
/*<             AR(IJK,J) = AR(IJK,J) + PRODR >*/
		    ar[ijk + j * ar_dim1] += prodr;
/*<             AI(IJK,J) = AI(IJK,J) + PRODI >*/
		    ai[ijk + j * ai_dim1] += prodi;
/*<  40         CONTINUE >*/
/* L40: */
		}
/*<             GO TO 48 >*/
		goto L48;
/*<           END IF >*/
	    }
/*<           IF (TR .EQ. 0.D0) THEN >*/
	    if (tr == 0.) {
/*<             DO 45 I = MD1,MDL >*/
		i__3 = mdl;
		for (i__ = md1; i__ <= i__3; ++i__) {
/*<             IJK = I - JK >*/
		    ijk = i__ - jk;
/*<             PRODR=-AI(I,K)*TI >*/
		    prodr = -ai[i__ + k * ai_dim1] * ti;
/*<             PRODI=AR(I,K)*TI >*/
		    prodi = ar[i__ + k * ar_dim1] * ti;
/*<             AR(IJK,J) = AR(IJK,J) + PRODR >*/
		    ar[ijk + j * ar_dim1] += prodr;
/*<             AI(IJK,J) = AI(IJK,J) + PRODI >*/
		    ai[ijk + j * ai_dim1] += prodi;
/*<  45         CONTINUE >*/
/* L45: */
		}
/*<             GO TO 48 >*/
		goto L48;
/*<           END IF >*/
	    }
/*<           DO 47 I = MD1,MDL >*/
	    i__3 = mdl;
	    for (i__ = md1; i__ <= i__3; ++i__) {
/*<             IJK = I - JK >*/
		ijk = i__ - jk;
/*<             PRODR=AR(I,K)*TR-AI(I,K)*TI >*/
		prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * 
			ti;
/*<             PRODI=AI(I,K)*TR+AR(I,K)*TI >*/
		prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * 
			ti;
/*<             AR(IJK,J) = AR(IJK,J) + PRODR >*/
		ar[ijk + j * ar_dim1] += prodr;
/*<             AI(IJK,J) = AI(IJK,J) + PRODI >*/
		ai[ijk + j * ai_dim1] += prodi;
/*<  47         CONTINUE >*/
/* L47: */
	    }
/*<  48       CONTINUE >*/
L48:
/*<  50       CONTINUE >*/
/* L50: */
	    ;
	}
/*<  55     CONTINUE >*/
L55:
/*<  60     CONTINUE >*/
/* L60: */
	;
    }
/*<  70   K = N >*/
L70:
    k = *n;
/*<       IF (DABS(AR(MD,N))+DABS(AI(MD,N)) .EQ. 0.D0) GO TO 80 >*/
    if ((d__1 = ar[md + *n * ar_dim1], abs(d__1)) + (d__2 = ai[md + *n * 
	    ai_dim1], abs(d__2)) == 0.) {
	goto L80;
    }
/*<       RETURN >*/
    return 0;
/*<  80   IER = K >*/
L80:
    *ier = k;
/*<       IP(N) = 0 >*/
    ip[*n] = 0;
/*<       RETURN >*/
    return 0;
/* ----------------------- END OF SUBROUTINE DECBC ------------------------ */
/*<       END >*/
} /* decbc_ */



/*<       SUBROUTINE SOLBC (N, NDIM, AR, AI, ML, MU, BR, BI, IP) >*/
/* Subroutine */ int solbc_(integer *n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *ml, integer *mu, doublereal *br, doublereal *
	bi, integer *ip)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k, m, kb, md, lm;
    static doublereal ti, tr;
    static integer md1, nm1;
    static doublereal den;
    static integer imd, kmd, mdl, mdm;
    static doublereal prodi, prodr;

/*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
/*<       DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N) >*/
/* ----------------------------------------------------------------------- */
/*  SOLUTION OF LINEAR SYSTEM, A*X = B , */
/*                  VERSION BANDED AND COMPLEX-DOUBLE PRECISION. */
/*  INPUT.. */
/*    N      ORDER OF MATRIX A. */
/*    NDIM   DECLARED DIMENSION OF ARRAY  A . */
/*    AR, AI TRIANGULARIZED MATRIX OBTAINED FROM DECB (REAL AND IMAG. PART). */
/*    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*    BR, BI RIGHT HAND SIDE VECTOR (REAL AND IMAG. PART). */
/*    IP     PIVOT VECTOR OBTAINED FROM DECBC. */
/*  DO NOT USE IF DECB HAS SET IER .NE. 0. */
/*  OUTPUT.. */
/*    BR, BI SOLUTION VECTOR, X (REAL AND IMAG. PART). */
/* ----------------------------------------------------------------------- */
/*<       MD = ML + MU + 1 >*/
    /* Parameter adjustments */
    --ip;
    --bi;
    --br;
    ai_dim1 = *ndim;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *ndim;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;

    /* Function Body */
    md = *ml + *mu + 1;
/*<       MD1 = MD + 1 >*/
    md1 = md + 1;
/*<       MDM = MD - 1 >*/
    mdm = md - 1;
/*<       NM1 = N - 1 >*/
    nm1 = *n - 1;
/*<       IF (ML .EQ. 0) GO TO 25 >*/
    if (*ml == 0) {
	goto L25;
    }
/*<       IF (N .EQ. 1) GO TO 50 >*/
    if (*n == 1) {
	goto L50;
    }
/*<       DO 20 K = 1,NM1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/*<         M = IP(K) >*/
	m = ip[k];
/*<         TR = BR(M) >*/
	tr = br[m];
/*<         TI = BI(M) >*/
	ti = bi[m];
/*<         BR(M) = BR(K) >*/
	br[m] = br[k];
/*<         BI(M) = BI(K) >*/
	bi[m] = bi[k];
/*<         BR(K) = TR >*/
	br[k] = tr;
/*<         BI(K) = TI >*/
	bi[k] = ti;
/*<         MDL = MIN(ML,N-K) + MD >*/
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	mdl = min(i__2,i__3) + md;
/*<         DO 10 I = MD1,MDL >*/
	i__2 = mdl;
	for (i__ = md1; i__ <= i__2; ++i__) {
/*<           IMD = I + K - MD >*/
	    imd = i__ + k - md;
/*<           PRODR=AR(I,K)*TR-AI(I,K)*TI >*/
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
/*<           PRODI=AI(I,K)*TR+AR(I,K)*TI >*/
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
/*<           BR(IMD) = BR(IMD) + PRODR >*/
	    br[imd] += prodr;
/*<           BI(IMD) = BI(IMD) + PRODI >*/
	    bi[imd] += prodi;
/*<  10     CONTINUE >*/
/* L10: */
	}
/*<  20     CONTINUE >*/
/* L20: */
    }
/*<  25     CONTINUE >*/
L25:
/*<       DO 40 KB = 1,NM1 >*/
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
/*<         K = N + 1 - KB >*/
	k = *n + 1 - kb;
/*<         DEN=AR(MD,K)*AR(MD,K)+AI(MD,K)*AI(MD,K) >*/
	den = ar[md + k * ar_dim1] * ar[md + k * ar_dim1] + ai[md + k * 
		ai_dim1] * ai[md + k * ai_dim1];
/*<         PRODR=BR(K)*AR(MD,K)+BI(K)*AI(MD,K) >*/
	prodr = br[k] * ar[md + k * ar_dim1] + bi[k] * ai[md + k * ai_dim1];
/*<         PRODI=BI(K)*AR(MD,K)-BR(K)*AI(MD,K) >*/
	prodi = bi[k] * ar[md + k * ar_dim1] - br[k] * ai[md + k * ai_dim1];
/*<         BR(K)=PRODR/DEN >*/
	br[k] = prodr / den;
/*<         BI(K)=PRODI/DEN >*/
	bi[k] = prodi / den;
/*<         TR = -BR(K) >*/
	tr = -br[k];
/*<         TI = -BI(K) >*/
	ti = -bi[k];
/*<         KMD = MD - K >*/
	kmd = md - k;
/*<         LM = MAX0(1,KMD+1) >*/
/* Computing MAX */
	i__2 = 1, i__3 = kmd + 1;
	lm = max(i__2,i__3);
/*<         DO 30 I = LM,MDM >*/
	i__2 = mdm;
	for (i__ = lm; i__ <= i__2; ++i__) {
/*<           IMD = I - KMD >*/
	    imd = i__ - kmd;
/*<           PRODR=AR(I,K)*TR-AI(I,K)*TI >*/
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
/*<           PRODI=AI(I,K)*TR+AR(I,K)*TI >*/
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
/*<           BR(IMD) = BR(IMD) + PRODR >*/
	    br[imd] += prodr;
/*<           BI(IMD) = BI(IMD) + PRODI >*/
	    bi[imd] += prodi;
/*<  30       CONTINUE >*/
/* L30: */
	}
/*<  40     CONTINUE >*/
/* L40: */
    }
/*<         DEN=AR(MD,1)*AR(MD,1)+AI(MD,1)*AI(MD,1) >*/
    den = ar[md + ar_dim1] * ar[md + ar_dim1] + ai[md + ai_dim1] * ai[md + 
	    ai_dim1];
/*<         PRODR=BR(1)*AR(MD,1)+BI(1)*AI(MD,1) >*/
    prodr = br[1] * ar[md + ar_dim1] + bi[1] * ai[md + ai_dim1];
/*<         PRODI=BI(1)*AR(MD,1)-BR(1)*AI(MD,1) >*/
    prodi = bi[1] * ar[md + ar_dim1] - br[1] * ai[md + ai_dim1];
/*<         BR(1)=PRODR/DEN >*/
    br[1] = prodr / den;
/*<         BI(1)=PRODI/DEN >*/
    bi[1] = prodi / den;
/*<  50   CONTINUE >*/
L50:
/*<       RETURN >*/
    return 0;
/* ----------------------- END OF SUBROUTINE SOLBC ------------------------ */
/*<       END >*/
} /* solbc_ */



/*<       subroutine elmhes(nm,n,low,igh,a,int) >*/
/* Subroutine */ int elmhes_(integer *nm, integer *n, integer *low, integer *
	igh, doublereal *a, integer *int__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, m;
    static doublereal x, y;
    static integer la, mm1, kp1, mp1;


/*<       integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1 >*/
/*<       real*8 a(nm,n) >*/
/*<       real*8 x,y >*/
/*<       real*8 dabs >*/
/*<       integer int(igh) >*/

/*     this subroutine is a translation of the algol procedure elmhes, */
/*     num. math. 12, 349-368(1968) by martin and wilkinson. */
/*     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971). */

/*     given a real general matrix, this subroutine */
/*     reduces a submatrix situated in rows and columns */
/*     low through igh to upper hessenberg form by */
/*     stabilized elementary similarity transformations. */

/*     on input: */

/*      nm must be set to the row dimension of two-dimensional */
/*        array parameters as declared in the calling program */
/*        dimension statement; */

/*      n is the order of the matrix; */

/*      low and igh are integers determined by the balancing */
/*        subroutine  balanc.      if  balanc  has not been used, */
/*        set low=1, igh=n; */

/*      a contains the input matrix. */

/*     on output: */

/*      a contains the hessenberg matrix.  the multipliers */
/*        which were used in the reduction are stored in the */
/*        remaining triangle under the hessenberg matrix; */

/*      int contains information on the rows and columns */
/*        interchanged in the reduction. */
/*        only elements low through igh are used. */

/*     questions and comments should be directed to b. s. garbow, */
/*     applied mathematics division, argonne national laboratory */

/*     ------------------------------------------------------------------ */

/*<       la = igh - 1 >*/
    /* Parameter adjustments */
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --int__;

    /* Function Body */
    la = *igh - 1;
/*<       kp1 = low + 1 >*/
    kp1 = *low + 1;
/*<       if (la .lt. kp1) go to 200 >*/
    if (la < kp1) {
	goto L200;
    }

/*<       do 180 m = kp1, la >*/
    i__1 = la;
    for (m = kp1; m <= i__1; ++m) {
/*<        mm1 = m - 1 >*/
	mm1 = m - 1;
/*<        x = 0.0d0 >*/
	x = 0.;
/*<        i = m >*/
	i__ = m;

/*<        do 100 j = m, igh >*/
	i__2 = *igh;
	for (j = m; j <= i__2; ++j) {
/*<           if (dabs(a(j,mm1)) .le. dabs(x)) go to 100 >*/
	    if ((d__1 = a[j + mm1 * a_dim1], abs(d__1)) <= abs(x)) {
		goto L100;
	    }
/*<           x = a(j,mm1) >*/
	    x = a[j + mm1 * a_dim1];
/*<           i = j >*/
	    i__ = j;
/*<   100   continue >*/
L100:
	    ;
	}

/*<        int(m) = i >*/
	int__[m] = i__;
/*<        if (i .eq. m) go to 130 >*/
	if (i__ == m) {
	    goto L130;
	}
/*    :::::::::: interchange rows and columns of a :::::::::: */
/*<        do 110 j = mm1, n >*/
	i__2 = *n;
	for (j = mm1; j <= i__2; ++j) {
/*<           y = a(i,j) >*/
	    y = a[i__ + j * a_dim1];
/*<           a(i,j) = a(m,j) >*/
	    a[i__ + j * a_dim1] = a[m + j * a_dim1];
/*<           a(m,j) = y >*/
	    a[m + j * a_dim1] = y;
/*<   110   continue >*/
/* L110: */
	}

/*<        do 120 j = 1, igh >*/
	i__2 = *igh;
	for (j = 1; j <= i__2; ++j) {
/*<           y = a(j,i) >*/
	    y = a[j + i__ * a_dim1];
/*<           a(j,i) = a(j,m) >*/
	    a[j + i__ * a_dim1] = a[j + m * a_dim1];
/*<           a(j,m) = y >*/
	    a[j + m * a_dim1] = y;
/*<   120   continue >*/
/* L120: */
	}
/*    :::::::::: end interchange :::::::::: */
/*<   130   if (x .eq. 0.0d0) go to 180 >*/
L130:
	if (x == 0.) {
	    goto L180;
	}
/*<        mp1 = m + 1 >*/
	mp1 = m + 1;

/*<        do 160 i = mp1, igh >*/
	i__2 = *igh;
	for (i__ = mp1; i__ <= i__2; ++i__) {
/*<           y = a(i,mm1) >*/
	    y = a[i__ + mm1 * a_dim1];
/*<           if (y .eq. 0.0d0) go to 160 >*/
	    if (y == 0.) {
		goto L160;
	    }
/*<           y = y / x >*/
	    y /= x;
/*<           a(i,mm1) = y >*/
	    a[i__ + mm1 * a_dim1] = y;

/*<           do 140 j = m, n >*/
	    i__3 = *n;
	    for (j = m; j <= i__3; ++j) {
/*<   140      a(i,j) = a(i,j) - y * a(m,j) >*/
/* L140: */
		a[i__ + j * a_dim1] -= y * a[m + j * a_dim1];
	    }

/*<           do 150 j = 1, igh >*/
	    i__3 = *igh;
	    for (j = 1; j <= i__3; ++j) {
/*<   150      a(j,m) = a(j,m) + y * a(j,i) >*/
/* L150: */
		a[j + m * a_dim1] += y * a[j + i__ * a_dim1];
	    }

/*<   160   continue >*/
L160:
	    ;
	}

/*<   180 continue >*/
L180:
	;
    }

/*<   200 return >*/
L200:
    return 0;
/*    :::::::::: last card of elmhes :::::::::: */
/*<       end >*/
} /* elmhes_ */

