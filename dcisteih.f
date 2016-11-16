!     ------------------------------------------------------------------
!     PROGRAM : DCI.
!     MODULE  : DCISteih.
!     VERSION : 2.0.
!     DATE    : June, 2000.
!     CONTENTS: subroutine dcisteih
!
!     External routines called by this module:
!
!     From dcicute : Hprod. 
!     From dcichol : NAproj.
!     From dciblas : vfill.
!     From the blas: dcopy, dscal, ddot, daxpy.
!     ------------------------------------------------------------------

      SUBROUTINE dcisteih(m, n, x, lambda, g, gp, A, Airow, Apcol, 
     $                    UseChol, Diag, Mval, Mirow, Mpcol, Mcnel, 
     $                    perm, iperm, DenseRow, delta, eps1, eps2, 
     $                    eps3, maxit, epscga, epscgr, maxitcg, 
     $                    prcndiag, s, qs, gts, iter, flag, Prec, 
     $                    tmp, tmp1, tmp2, tmp3, tmp4, tmp5, r, p, 
     $                    u, v, snew, nHprod, GotH)
C    $                    , Hprod)

C     This routine computes the horizontal step by
C     using Steihaug's algorithm.
C
C     The horizontal step is an approximate solution for:
C
C     min  0.5*s'*H*s + g'*s
C     s.t. A*s = 0
C          ||s|| < delta
C
C     Input parameters:
C
C     m,n            number of constraints and variables.
C     x              current solution (used to compute the Hessian H).
C     lambda         current Lagrangian multipliers (used to compute H).
C     g              vector that defines the objective function.
C     gp             initial residual vector.
C     A,Airow,Apcol  matrix that defines the linear constraints.
C     UseChol        logical parameter that indicates if the Cholesky
C                    decomposition of A*A' is to be computed (.TRUE.)
C                    or if the CG method should be used to solve a 
C                    system with this matrix (.FALSE.).
C     Diag, Mval, Mirow, Mpcol, Mcnel, perm, iperm, DenseRow
C                    Data structure for the Cholesky decompostion of
C                    A*A'. See the dcichol module for a description
C                    of these parameters.
C     delta          trust region radius.
C     eps1           relative stopping tolerance. The algorithm stops 
C                    when ||r||^2 <= eps1*(||r0||^2). Suggestion: 0.01.
C     eps2           absolute stopping tolerance. The algorithm stops 
C                    when ||r||^2 <= eps2. Suggestion: 1e-10.
C     eps3           relative tolerance used to identify a direction of 
C                    negative curvature.
C     maxit          maximum number of iterations.
C     epscga         Absolute stopping tolerance used when solving the 
C                    system (A'*A)x = y by the CG method.
C     epscgr         Relative stopping tolerance used when solving the 
C                    system (A'*A)x = y by the CG method.
C     maxitcg        Maximum number of iterations of the CG algorithm 
C                    (used to solve the system (A'*A)x = y).
C     epsnegc        Tolerance
C     prcndiag       Number of (upper and lower) diagonals of the band
C                    preconditioner used when solving (A'*A)x = y by the
C                    CG method (zero means no preconditioning at all).
C     Prec           Vector that stores the band preconditioner. 
C     nHprod         total number of Hessian-vector products.
C     GotH           variable that indicates if the Hessian of the
C                    Lagrangian was already computed at (x, lambda).
C
C     Output parameters:
C
C     s              horizontal step.
C     qs             objective function value at s.
C     gts            g'*s.
C     iter           total number of CG iterations.
C     flag           output flag. 
C                    If flag = 0, a solution s such that ||s||<delta 
C                       was found.
C                    If flag = 1, a direction of negative curvature 
C                       was found. ||s||=delta.
C                    If flag = 2, the minimizer is outside 
C                       the trust region. ||s||=delta.
C                    If flag = 3, the maximum number of iterations 
C                       was taken.
C     nHprod         see description above.
C     GotH           see description above.
C
C     Temporary vectors passed as parameters:
C
C     r(n), p(n), u(n), v(n), snew(n),
C     tmp(m), tmp1(m), tmp2(m), tmp3(m), tmp4(m), tmp5(n).

      IMPLICIT NONE

      EXTERNAL Hprod

C     Subroutine parameters.

      INTEGER m, n, maxit, iter, flag, nHprod, maxitcg, prcndiag
      INTEGER DenseRow
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      INTEGER Mirow(*)
      INTEGER Mpcol(m+1)
      INTEGER Mcnel(m), perm(m), iperm(m)
      REAL*8  delta, eps1, eps2, eps3, qs, epscga, epscgr, gts
      REAL*8  x(n), g(n), gp(n), s(n), r(n), p(n), u(n), v(n)
      REAL*8  snew(n), tmp5(n)
      REAL*8  tmp(m), tmp1(m), tmp2(m), tmp3(m), tmp4(m)
      REAL*8  lambda(m), Diag(m)
      REAL*8  A(*)
      REAL*8  Mval(*)
      REAL*8  Prec(*)
      LOGICAL GotH, UseChol

C     Local variables.

      INTEGER naflag
      REAL*8  alpha, beta, gamma, theta, theta0, thetanew, delta2
      REAL*8  zeta, zeta2, ptp, sts, gtp, stp, stu, stsnew, qsnew
      REAL*8  epscg, mu, Lbdmax
      LOGICAL LimLbd, goodx0

C     Functions called by the routine.

      REAL*8  dnrm2, ddot

      delta2 = delta*delta
      CALL vfill(n, 0.0D0, s, 1)
      qs = 0.0D0
      sts = 0.0D0

C     gp is the initial residual vector r.

      CALL dcopy (n, gp, 1, r, 1)
      CALL dscal(n, -1.0D0, r, 1)
      CALL dcopy(n, r, 1, p, 1)
      theta0 = ddot(n, r, 1, r, 1)
      theta  = theta0
      gts    = 0.0D0
      iter   = 0

      DO WHILE ((theta.GT.eps2).AND.(theta.GT.(eps1*theta0)).AND.
     &         (iter.LE.maxit))

        iter = iter+1
        CALL CutestHprod(n, m, x, lambda, p, u, GotH)
        nHprod = nHprod + 1
        gamma = ddot(n, p, 1, u, 1)
        stu = ddot(n, s, 1, u, 1)
        gtp = ddot(n, g, 1, p, 1)
        ptp = ddot(n, p, 1, p, 1)

        IF (gamma.LE.(eps3*ptp)) THEN

C         We have found a direction of negative curvature.

          stp = ddot(n, s, 1, p, 1)
          zeta = DSQRT(stp*stp + (delta2 - sts)*ptp)
          zeta2 = (-stp - zeta)/ptp
          zeta = (-stp + zeta)/ptp
          CALL dcopy(n, s, 1, snew, 1) 
          CALL daxpy(n, zeta2, p, 1, snew, 1)
          qsnew = qs + zeta2*stu + 0.5D0*zeta2*zeta2*gamma + zeta2*gtp
          CALL daxpy(n, zeta, p, 1, s, 1)
          qs = qs + zeta*stu + 0.5D0*zeta*zeta*gamma + zeta*gtp

          IF (qsnew.LT.qs) THEN
            CALL dcopy(n, snew, 1, s, 1)
            qs  = qsnew
            gts = gts + zeta2*gtp;
          ELSE
            gts = gts + zeta*gtp;
          ENDIF

          flag = 1
          RETURN

        ENDIF

        alpha = theta/gamma
        CALL dcopy(n, s, 1, snew, 1)
        CALL daxpy(n, alpha, p, 1, snew, 1)
        stsnew = ddot(n, snew, 1, snew, 1)

        IF (stsnew.GT.delta2) THEN

C         Keeping the step inside the trust region boundary.

          stp = ddot(n, s, 1, p, 1)
          zeta = (-stp + DSQRT(stp*stp + (delta2 - sts)*ptp))/ptp
          CALL daxpy(n, zeta, p, 1, s, 1)
          gts = gts+zeta*gtp;
          qs = qs + zeta*stu + 0.5D0*zeta*zeta*gamma + zeta*gtp
          flag = 2
          RETURN

        ENDIF

        CALL daxpy(n, -alpha, u, 1, r, 1)
        CALL dcopy(n, r, 1, v, 1)

        IF (UseChol) THEN

          LimLbd = .FALSE.
          lbdmax = 1.0D20
          CALL NAprojCh(m, n, A, Airow, Apcol, Diag, Mval, Mirow, 
     $                  Mpcol, Mcnel, perm, iperm, DenseRow, tmp, 
     $                  LimLbd, lbdmax, v, r, naflag)

        ELSE

          epscg = DMAX1(epscga, epscgr*DMIN1(1.0D0, SQRT(theta)))
          mu = 1.0D-10
          LimLbd = .FALSE.
          goodx0 = .FALSE.
          lbdmax = 1.0D20
          CALL NAprojCG(m, n, prcndiag, goodx0, maxitcg, epscg,
     $                  mu, A, Airow, Apcol, Prec, tmp1, tmp2, 
     $                  tmp3, tmp4, tmp5, tmp, LimLbd, lbdmax,
     $                  v, r, naflag)

        ENDIF

        thetanew = ddot(n, r, 1, r, 1)
        beta = thetanew/theta
        qs = qs + alpha*stu + 0.5D0*alpha*alpha*gamma + alpha*gtp
        CALL dscal(n, beta, p, 1)
        CALL daxpy(n, 1.0D0, r, 1, p, 1) 
        CALL dcopy(n, snew, 1, s, 1)
        gts = gts+alpha*gtp;
        sts = stsnew
        theta = thetanew

      ENDDO

      IF (iter.GT.maxit) THEN
        flag = 3
      ELSE
        flag = 0
      ENDIF

      RETURN
      END

