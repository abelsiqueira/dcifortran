!     ------------------------------------------------------------------
!     PROGRAM : DCI.
!     MODULE  : DCICG.
!     VERSION : 2.0.
!     DATE    : June, 2005.
!     CONTENTS: subroutine NAprojCG
!               subroutine NewStepCG
!               subroutine MsolveCG
!
!     External routines called by this module:
!
!     From dciblas : MVprod, MtVprod, vfill, vfills, ifill, icopy,
!                    vcmpr, icmpr, axpys, dots.
!     From the blas: dcopy, dscal, daxpy, ddot.
!     ------------------------------------------------------------------

      SUBROUTINE NAprojCG(m, n, prcndiag, goodx0, maxitcg, epscg, 
     $                    mu, A, Airow, Apcol, Prec, tmp1, tmp2,
     $                    tmp3, tmp4, tmp5, tmp, LimTmp, tmpmax,
     $                    r, Pr, flag)

C     This routine computes the projection of r onto N(A).
C     Pr = r - A'*inv(A*A')*A*r.
C     If A*A' is nonsingular, then mu should be set to zero. 
C     When A*A' is singular, mu may be used to regularize the problem, 
C     so (A*A'+mu*I) is used instead of A*A'. 
C     On output, tmp stores -inv(A*A')*A*r. When r is the gradient 
C     of the objective function of DCI, tmp is equal to lambda, the
C     vector of Lagrange multipliers. LimTmp is a logical parameter 
C     that indicates if tmp is to be bounded. If LimTmp = .TRUE., 
C     tmp = max(min(tmp, tmpmax), -tmpmax). flag is set to a value 
C     greater than zero when the CG method fails.

      IMPLICIT NONE

C     Routine parameters.

      INTEGER m, n, prcndiag, maxitcg, flag
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      REAL*8  epscg, mu, tmpmax
      REAL*8  A(*)
      REAL*8  tmp(m), tmp1(m), tmp2(m), tmp3(m), tmp4(m)
      REAL*8  Prec(*)
      REAL*8  r(n), Pr(n), tmp5(n)
      LOGICAL goodx0, LimTmp

C     Internal variables.

      INTEGER i, itercg
      REAL*8  relres


      IF (m.GT.0) THEN
        CALL MVprod(m, n, A, r, Pr, Airow, Apcol, 0.0D0)
        CALL MsolveCG(m, n, maxitcg, epscg, mu, A, Airow, Apcol,  
     $                prcndiag, Prec, tmp1, tmp2, tmp3, tmp4, tmp5, 
     $                Pr, goodx0, tmp, relres, itercg, flag)

C       If (flag.NE.0) WRITE(*,*)'maxitcg, flag: ',maxitcg,' ',flag

        CALL dscal(m, -1.0D0, tmp, 1)
        IF (LimTmp) THEN
          DO i = 1,m
            IF (tmp(i).GT.tmpmax) tmp(i) = tmpmax
          ENDDO
          DO i = 1,m
            IF (tmp(i).LT.(-tmpmax)) tmp(i) = -tmpmax
          ENDDO
        ENDIF
        CALL MtVprod(n, A, tmp, Pr, Airow, Apcol) 
        CALL daxpy(n, 1.0D0, r, 1, Pr, 1)
      ELSE
        CALL dcopy(n, r, 1, Pr, 1)
        flag = 0
      ENDIF

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE NewStepCG(m, n, prcndiag, maxitcg, epscg, mu,
     $                     A, Airow, Apcol, Prec, tmp1, tmp2,
     $                     tmp3, tmp4, tmp5, tmp, r, s, flag)

C     This routine computes s = -A'*inv(A*A'+mu*I)*r.
C     tmp(m),tmp1(m),...,tmp4(m) and tmp5(n) are temporary vectors.
C     On output, flag > 0 indicates that the CG method has failed.

      IMPLICIT NONE

C     Routine parameters.

      INTEGER m, n, prcndiag, maxitcg, flag
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      REAL*8  epscg, mu
      REAL*8  A(*)
      REAL*8  r(m), tmp(m), tmp1(m), tmp2(m), tmp3(m), tmp4(m)
      REAL*8  Prec(*)
      REAL*8  s(n), tmp5(n)

C     Internal variables.

      INTEGER itercg
      REAL*8  relres

c$$$C     Local variables user to compute the norm of the residual.
c$$$
c$$$      REAL*8  rnorm, bnorm
c$$$
c$$$C     Functions called when computing the norm of the residual.
c$$$
c$$$      REAL*8  ResNorm, dnrm2

      IF (m.GT.0) THEN

        CALL dcopy(m, r, 1, s, 1)
        CALL MsolveCG(m, n, maxitcg, epscg, mu, A, Airow, Apcol,  
     $                prcndiag, Prec, tmp1, tmp2, tmp3, tmp4, tmp5, 
     $                s, .FALSE., tmp, relres, itercg, flag)

c$$$        rnorm = ResNorm(m, n, A, Airow, Apcol, tmp, r)
c$$$        bnorm = dnrm2(m, r, 1)
c$$$        IF ((rnorm/(1.0D0+bnorm)).GT.1.0D-5) THEN
c$$$          WRITE(*,*)'NAprojCh warning: inacurate solution of A*Atx = b'
c$$$          WRITE(*,*)'Residual norm   : ', rnorm
c$$$          WRITE(*,*)'MSolveCG relres : ', relres
c$$$          WRITE(*,*)'MSolveCG iter   : ', itercg
c$$$          WRITE(*,*)'MSolveCG flag   : ', flag
c$$$          WRITE(*,*)'MSolveCG mu     : ', mu
c$$$          STOP
c$$$        ENDIF

        CALL dscal(m, -1.0D0, tmp, 1)
        CALL MtVprod(n, A, tmp, s, Airow, Apcol)

      ELSE

        CALL vfill(n, 0.0D0, s, 1)
        flag = 0

      ENDIF

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE MsolveCG(m, n, maxit, tol, mu, A, Airow, Apcol,  
     $                    prcndiag, Prec, p, z, r, Ap, y, b, goodx0,  
     $                    x, relres, iter, flag)

C     This routine solves the system (A.A'+mu*I).x = b by the 
C     preconditioned conjugate gradient method. Usually, mu is a 
C     small penalty parameter used to avoid the numerical problems 
C     that arise when A.A' is singular. mu may be set to zero if
C     A.A' is nonsingular.
C
C     Input parameters:
C
C     m        Number of rows in A.
C     n        Number of columns in A.
C     maxit    Maximum number of iterations
C     tol      Stopping tolerance.
C     mu       penalty parameter as explained above.
C     A, Airow, Apcol Matrix A given by columns in sparse format.
C     prcndiag Total number of (upper and lower) diagonals of the band
C              preconditioner (zero means no preconditioning at all). 
C     Prec     Vector that stores the band preconditioner. 
C     b        RHS vector.
C     goodx0   Logical variable that indicates if x stores a good 
C              starting vector (goodx0=.TRUE.) or x = 0 is to be used.
C     x        Starting vector (if goodx0=.TRUE.).
C
C     Temporary real vectors passed as parameters.
C
C     p(m), z(m), r(m), Ap(m), y(n).
C
C     Output parameters:
C
C     x        Approximate solution.
C     relres   ||b - (A*A'+mu*I)*x|| / ||b||.
C     iter     Total number of iterations spent to obtain x.
C     flag     Output flag:
C              If flag = 0, the CG method converged.
C              If flag = 1, maxit iterations were performed.
C              If flag = 2, the preconditioner is very ill-conditioned.
C              If flag = 3, alpha = 0.
C              If flag = 4, (r'*z) is out of bounds or alpha = Inf.

      IMPLICIT NONE

C     Routine parameters.

      INTEGER m, n, prcndiag, maxit, iter, flag
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      REAL*8  tol, mu, relres
      REAL*8  b(m), x(m)
      REAL*8  A(*)
      REAL*8  Prec(*)
      LOGICAL goodx0

C     Temporary vectors passed as parameters.

      REAL*8  p(m), z(m), r(m), Ap(m)
      REAL*8  y(n)

C     Internal variables.

      INTEGER prcnband
      REAL*8  bnorm, rnorm, rtol, rtz, rtzold, alpha, beta
      REAL*8  dnrm2, ddot

C     Copying prcndiag onto prcnband (prcnband may vary according
C     to the availability of a good preconditioner).

      prcnband = prcndiag

C     Checking if b = 0, so x = 0 is a solution.

      bnorm = dnrm2(m, b, 1)
      IF (bnorm.LT.1.0D-20) THEN
        CALL vfill(m, 0.0D0, x, 1)
        flag = 0
        relres = bnorm
        iter = 0
        RETURN
      ENDIF

C     Computing the initial residual.

      IF (goodx0) THEN

C       A good starting point was given.

        CALL dcopy(m, b, 1, r, 1)
        CALL MtVprod(n, A, x, y, Airow, Apcol)
        CALL MVprod(m, n, A, y, Ap, Airow, Apcol, 0.0D0)
        IF (mu.NE.0.0D0) CALL daxpy(m, mu, x, 1, Ap, 1)
        CALL daxpy(m, -1.0D0, Ap, 1, r, 1)
        rnorm = dnrm2(m, r, 1)

      ELSE

C       Using x0 = 0.

        CALL dcopy(m, b, 1, r, 1)
        CALL vfill(m, 0.0D0, x, 1)
        rnorm = bnorm

      ENDIF

      rtol = bnorm*tol
      rtz = 0.0D0
      flag = 0
      iter = 0
      
      DO WHILE ((iter.LE.maxit).AND.(rnorm.GT.rtol).AND.(flag.EQ.0))

C       Preconditioning.     

        IF (prcnband.GT.0) THEN
          CALL dciprec(m, n, mu, A, Airow, Apcol, prcnband, Prec, r, z) 
          IF (dnrm2(m, z, 1).GT.1.0D99) THEN
            flag = 2
            relres = rnorm/bnorm
            RETURN
          ENDIF
        ENDIF

C       Computing r'*z.  

        rtzold = rtz
        IF (prcnband.GT.0) THEN
          rtz = ddot(m, r, 1, z, 1)
        ELSE
          rtz = rnorm**2
        ENDIF
        IF ((DABS(rtz).LT.1.0D-50).OR.(DABS(rtz).GT.1.0D99)) THEN
          flag = 4
          relres = rnorm/bnorm
          RETURN
        ENDIF

C       Updating p.

        IF (iter.EQ.0) THEN
          IF (prcnband.GT.0) THEN
            CALL dcopy(m, z, 1, p, 1)
          ELSE
            CALL dcopy(m, r, 1, p, 1)
          ENDIF
        ELSE
          beta = rtz/rtzold
          CALL dscal(m, beta, p, 1)
          IF (prcnband.GT.0) THEN
            CALL daxpy(m, 1.0D0, z, 1, p, 1)
          ELSE
            CALL daxpy(m, 1.0D0, r, 1, p, 1)
          ENDIF
        ENDIF

C       Computing (A*A'+mu*I)*p.

        CALL MtVprod(n, A, p, y, Airow, Apcol)
        CALL MVprod(m, n, A, y, Ap, Airow, Apcol, 0.0D0)
        IF (mu.NE.0.0D0) CALL daxpy(m, mu, p, 1, Ap, 1)

C       Computing alpha.

        alpha = rtz/ddot(m, p, 1, Ap, 1)
        IF (DABS(alpha).LT.1.0D-50) THEN
          flag = 1
          relres = rnorm/bnorm
          RETURN
        ELSE IF (DABS(alpha).GT.1.0D99) THEN
          flag = 4
          relres = rnorm/bnorm
          RETURN
        ENDIF

C       Updating x and r.

        CALL daxpy(m, alpha, p, 1, x, 1)
        CALL daxpy(m, -alpha, Ap, 1, r, 1)
        iter = iter + 1

      ENDDO

C     Adjusting flag if the maximum number of iterations was exceeded.

      IF (iter.GT.maxit) flag = 1
      relres = rnorm/bnorm

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE dciprec(m, n, mu, A, Airow, Apcol, 
     $                   prcnband, Prec, b, x)

C     This routine solves the system Mx = b using the CG method.
C     M is a band preconditioner build from A*A'+mu*I. If prcavail 
C     is .FALSE., then the Cholesky factor of the decomposed 
C     preconditioner is returned in L. Otherwise, L is supposed
C     to computed in a previous call to this routine.
C     CAUTION: prcnband is changed inside this routine. 

      IMPLICIT NONE

C     Routine parameters.

      INTEGER m, n, prcnband
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      REAL*8  mu
      REAL*8  b(m), x(m)
      REAL*8  A(*)
      REAL*8  Prec(*)

C     Variable passed by COMMON.

      INTEGER prcavail

C     Internal variable.

      INTEGER flag

C     COMMON statements.

      COMMON / dciprecavail / prcavail

      flag = 0

C     Changed from .NOT.prcavail to prcavail.EQ.0
      IF (prcavail.EQ.0) THEN

C       Computing the preconditioner M.

        CALL formband(m, n, prcnband, mu, A, Airow, Apcol, Prec)
        CALL dcopy(m, Prec(MIN((prcnband-1)/2,m-1)*m+1), 1, x, 1)
        CALL bandchol(m, prcnband, Prec, flag)

C       Using the diagonal preconditioner if a larger band didn't work.

        IF ((flag.NE.0).AND.(prcnband.GT.1)) THEN
          prcnband = 1
          CALL dcopy(m, x, 1, Prec, 1) 
          CALL bandchol(m, prcnband, Prec, flag)
        ENDIF

C       Discarding the preconditioner if things go wrong.

        IF (flag.NE.0) THEN
          CALL dcopy(m, b, 1, x, 1)
          prcnband = 0
        ELSE
          prcavail = .TRUE.
        ENDIF

      ENDIF

C     Solving M*x = b.
 
      IF (prcnband.GT.0) THEN
        CALL precsolve(m, prcnband, Prec, x, b)
      ENDIF

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE formband(m, n, prcndiag, mu, A, Airow, Apcol, Prec)

C     This routine builds the band matrix Prec that is used as the 
C     preconditioner for the CG method. prcndiag is the total number
C     of diagonals in Prec. The (i,j)-element of the preconditioner
C     is stored in Prec((p+j-i)*m+j), where p is the (semi)bandwidth.
C     The same value corresponds to the (j,i)-element, since Prec is 
C     built from A*A'+mu*I and, therefore, is symmetric. A is given by 
C     columns, so "outer" products are used to compute Prec.

      IMPLICIT NONE

C     Routine parameters

      INTEGER m, n, prcndiag
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      REAL*8  mu
      REAL*8  A(*)
      REAL*8  Prec(*)

C     Internal variables.

      INTEGER i, j, k, p, row, col, imax
      REAL*8  v, maxAAt

C     Functions called by this routine.

      INTEGER idamax

C     Computing the (semi)bandwidth.

      p = MIN((prcndiag-1)/2,m-1)

C     Filling M entries with zeros.

      CALL vfill(m*(p+1), 0.0D0, Prec, 1)

C     Computing lower triangular part of M.

      DO k = 1, n
        DO i = Apcol(k), Apcol(k+1)-1
          row = Airow(i)
          v = A(i)
          DO j = Apcol(k), i
            col = Airow(j)
            IF ((row-col).LE.p) THEN
              Prec((p+col-row)*m+col) = Prec((p+col-row)*m+col) + v*A(j)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

C     Adding a small value to the diagonal of Band, trying to
C     turn this matrix into a positive definite one. The value
C     added is greater or equal to mu.

      imax = idamax(m, Prec(m*p+1), 1)
      maxAAt = DMAX1(DMIN1(1.0D-8,Prec(m*p+imax)*1.0D-10), mu)

      DO i = m*p+1,m*(p+1)
        Prec(i) = Prec(i) + maxAAt
      ENDDO

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE bandchol(m, prcndiag, Prec, flag)

C     This routine computes the Cholesky decomposition of a band matrix 
C     stored by diagonals in the vector Prec. Prec is supposed to be
C     previously generated by routine formband. The (lower triangular)
C     Cholesky factor is stored over Prec. prcndiag is the total number
C     of (upper and lower) diagonals of the preconditioner matrix.

      IMPLICIT NONE

C     Routine parameters.

      INTEGER m, prcndiag, flag
      REAL*8  Prec(*)

C     Internal variables.

      INTEGER p, j, k, lbd
      REAL*8  fac, maxfac

C     Defining maxfac.

      maxfac = 1.0D-99

C     Computing the (semi)bandwidth.

      p = MIN((prcndiag-1)/2,m-1)

C     Computing the Cholesky factor.

      DO j = 1, m
        DO k = MAX(1, j-p), j-1
          lbd = MIN(k+p, m)-j+1
          fac = -Prec((p+k-j)*m+k)
          CALL daxpy(lbd, fac, Prec((p+k-j-lbd+1)*m+k), m, 
     $               Prec((p-lbd+1)*m+j), m)
        ENDDO
        lbd = MIN(j+p, m)-j+1
        fac = Prec(p*m+j)
        IF (maxfac.LT.fac) THEN
          maxfac = fac
        ELSE IF (fac.LT.(1.0D-10*maxfac)) THEN
          IF (fac.LE.-1.0D-16) THEN
            flag = j
            RETURN
          ELSE
            fac = fac + 1.01D-10*maxfac
          ENDIF
        ENDIF
        CALL dscal(lbd, 1.0D0/SQRT(fac), Prec((p-lbd+1)*m+j), m)
      ENDDO

      flag = 0

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE precsolve(m, prcndiag, Prec, x, b)

C     This routine solves the system M.x = b, where M is a band matrix 
C     whose Cholesky factor is stored in the vector Prec. prcndiag is 
C     the total number of (upper and lower) diagonals of M.

      IMPLICIT NONE

C     Routine parameters.

      INTEGER m, prcndiag
      REAL*8  x(m), b(m)
      REAL*8  Prec(*)

C     Internal variables.

      INTEGER i, j, p
      REAL*8  t
      REAL*8  ddot

C     Computing the (semi)bandwidth.

      p = MIN((prcndiag-1)/2,m-1)

C     Copying b onto x.

      CALL dcopy(m, b, 1, x, 1)

C     Solving L.y = b.

      DO i = 1, m

C       Computing y(i).

        x(i) = x(i)/prec(m*p+i)

C       Updating the remaining variables.

C       CALL daxpy(MIN(p, m-i), -x(i), prec((p-1)*m+i), -m, x(i+1), 1)
        DO j = 1, MIN(p, m-i)
          x(i+j) = x(i+j) - x(i)*prec((p-j)*m+i)
        ENDDO

      ENDDO

C     Solving L'.x = y.

      x(m) = x(m)/prec(m*(p+1))
      DO i = m-1, 1, -1

C       Computing x(i).
 
C       t = ddot(MIN(p, m-i), prec((p-1)*m+i), -m, x(i+1), 1)
        t = 0.0D0
        DO j = 1, MIN(p, m-i)
          t = t + prec((p-j)*m+i)*x(i+j)
        ENDDO

        x(i) = (x(i)-t)/prec(m*p+i)

      ENDDO

      RETURN
      END
