!     ------------------------------------------------------------------
!     PROGRAM : DCI.
!     MODULE  : DCIDog.
!     VERSION : 2.0.
!     DATE    : June, 2000.
!     CONTENTS: subroutine dcitrust.
!
!     External routines called by this module:
!
!     From dcicute : Constr. 
!     From dcichol : NewStep.
!     From dciblas : MtVprod, MVprod.
!     From the blas: dcopy, dscal, ddot, daxpy, dnrm2.
!     ------------------------------------------------------------------

      SUBROUTINE dcitrust(n, m, xc, hnorm, A, Airow, Apcol, FreshA, 
     &                    delta, infdelta, kappa1, kappa2, kappa3, 
     &                    kappa4, UseChol, epscg, maxitcg, prcndiag,  
     &                    Prec, Diag, Mval, Mirow, Mpcol, Mcnel, perm,
     &                    iperm, DenseRow, tmp1, tmp2, tmp3, g, d, dcp,  
     &                    dn, h, x, hc, hcnorm, nConstr, iout) 
C    &                    , Constr)

C     This routine computes a trust region vertical step. This step
C     is the solution of the problem:
C
C     min  0.5*(||h(x)+A*d||^2 - ||h(x)||^2)
C     s.t. ||d|| <= delta.
C
C     The algorithm uses a variant of the dogleg method described 
C     by M. Powell in "A hybrid method for nonlinear equations", in 
C     Numerical Methods for Nonlinear Algebraic Equations, 
C     P. Rabinowitz, ed. 1970. 
C 
C     Input parameters:
C
C     m, n          number of constraints and variables.
C     xc            current approximation for the solution.
C     hc            constraint vector, h(x).
C     hnorm         two-norm of h(x).
C     A,Airow,Apcol Jacobian of the constraints.
C     delta         trust region radius.
C     infdelta      smallest admissible value for delta.
C     kappa1        parameter used to decide if the step is to be 
C                   accepted (suggestion: 1e-4).
C     kappa2        parameter used to reduce delta (suggestion: 0.25).
C     kappa3        parameter used to decide if delta is to be 
C                   increased (suggestion: 0.7).
C     kappa4        parameter used to increase delta (suggestion: 2.5).
C     UseChol       logical parameter that indicates if the Cholesky
C                   decomposition of A*A' is to be computed (.TRUE.)
C                   or if the CG method should be used to solve a 
C                   system with this matrix (.FALSE.).
C     epscg         Stopping tolerance used when solving the system 
C                   (A'*A)x = y by the CG method.
C     maxitcg       Maximum number of iterations of the CG algorithm 
C                   (used to solve the system (A'*A)x = y).
C     prcndiag      Number of (upper and lower) diagonals of the band
C                   preconditioner used when solving (A'*A)x = y by the
C                   CG method (zero means no preconditioning at all). 
C     Prec          Vector that stores the band preconditioner. 
C     nConstr       total number of constraint evaluations.
C     Constr        A function that computes the constraint vector.
C     FreshA        parameter that defines if we are using the exact
C                   Jacobian at x or an approximation of this matrix.
C     Diag, Mval, Mirow, Mpcol, Mcnel, perm, iperm, DenseRow
C                   Data structure for the Cholesky decompostion of
C                   A*A'. See the dcichol module for a description
C                   of these parameters.
C
C     Output parameters:
C
C     xc            new point.
C     hc            new constraint vector, h(xc).
C     hcnorm        two-norm of h(xc).
C     nConstr       see description above.
C     delta         see description above.
C     Diag, Mval, Mirow, Mpcol, Mcnel, perm, iperm, DenseRow
C                   Data structure for the Cholesky decompostion of
C                   A*A'. See the dcichol module for a description
C                   of these parameters.
C     iout          Exit condition:
C                     1: d = dcp (because dn was discarded).
C                     2: d = dn.
C                     3: ||dcp|| <= ||d|| = delta <= ||dn||.
C                     4: ||d|| = delta <= ||dcp||.
C                     5: failure: delta <= infdelta.
C                     6: ||dcp|| is too small.
C
C     Temporary vectors passed as parameters:
C
C     tmp1(m), tmp2(m), tmp3(m), g(n), d(n), dcp(n), dn(n), h(m), x(n).

      IMPLICIT NONE
 
      EXTERNAL Constr

C     Subroutine parameters.

      INTEGER m, n, nConstr, maxitcg, prcndiag, DenseRow, iout
      INTEGER Airow(*)
      INTEGER Mirow(*)
      INTEGER Apcol(n+1)
      INTEGER Mpcol(m+1)
      INTEGER Mcnel(m), perm(m), iperm(m)
      REAL*8  hnorm, hcnorm, delta, epscg
      REAL*8  infdelta, kappa1, kappa2, kappa3, kappa4
      REAL*8  xc(n), g(n), d(n), dcp(n), dn(n), x(n)
      REAL*8  hc(m), h(m), Diag(m), tmp1(m), tmp2(m), tmp3(m)
      REAL*8  A(*)
      REAL*8  Mval(*)
      REAL*8  Prec(*)
      LOGICAL FreshA, UseChol

C     Internal variables.

      INTEGER naflag, iter
      REAL*8  alpha, gnorm, dcpnorm, dnnorm, dctd, dctn, deltaold
      REAL*8  Ared, Pred, dnorm
      LOGICAL dnavail

C     Functions called by the routine.

      REAL*8  ddot, dnrm2

C     Computing the Cauchy step.

      CALL dcopy(m, hc, 1, h, 1)
      CALL MtVprod(n, A, h, g, Airow, Apcol) 
      gnorm = dnrm2(n, g, 1)
      IF (gnorm.LT.1.0D-20) THEN
        hcnorm = hnorm
        iout = 6
        RETURN
      ENDIF
      CALL MVprod(m, n, A, g, d, Airow, Apcol, 0.0D0)
      alpha = DMIN1(gnorm**2/(dnrm2(m, d, 1)**2), delta/gnorm)
      CALL dcopy(n, g, 1, dcp, 1)
      CALL dscal(n, -alpha, dcp, 1)
      dcpnorm = DABS(alpha)*gnorm

C     Setting some initial values.

      CALL dcopy(n, xc, 1, x, 1)
      dnavail  = .FALSE.
      dnnorm   = 0.0D0
      Ared     = 0.0D0
      Pred     = 1.0D0
      deltaold = delta
      delta    = delta/kappa2
      dnorm    = delta
      iout     = 0
      iter     = 0

C     Main loop (run the loop only once if FreshA is false).

      DO WHILE ((Ared.LT.(kappa1*Pred)).AND.(FreshA.OR.(iter.LT.2)))

        iter = iter + 1
        delta = kappa2*dnorm

        IF (dcpnorm.LT.delta) THEN

C         Cauchy step is inside the ball. 

          IF (.NOT.dnavail) THEN

C           Computing the Newton step if it is not available.

            IF (UseChol) THEN
              CALL NewStepCh(m, n, A, Airow, Apcol, Diag, Mval, Mirow,
     $                       Mpcol, Mcnel, perm, iperm, DenseRow, tmp1,
     $                       h, dn, naflag)

            ELSE
              CALL NewStepCG(m, n, prcndiag, maxitcg, epscg, 0.0D0,
     $                       A, Airow, Apcol, Prec, tmp1, tmp2,  
     $                       tmp3, Diag, g, d, h, dn, naflag)
            ENDIF

            dnavail = .TRUE.

C           Discarding dn if NewStep has failed.

            IF (naflag.GT.1) THEN
              dnnorm = 0.0D0
            ELSE
              dnnorm = dnrm2(n, dn, 1)
            ENDIF

          ENDIF

          IF (dnnorm.LE.dcpnorm) THEN

C           Discarding dn if its norm is too small.

            CALL dcopy(n, dcp, 1, d, 1)
            dnorm = dcpnorm
            iout = 1

          ELSE IF (dnnorm.LE.delta) THEN

C           Using dn if it is inside the ball.

            CALL dcopy(n, dn, 1, d, 1)
            dnorm = dnnorm
            iout = 2

          ELSE

C           Finding a point in the dogleg path (d = dn - dcp).

            dctn = ddot(n, dcp, 1, dn, 1)
            dnorm = DSQRT(dnnorm**2 - 2.0D0*dctn + dcpnorm**2)
            dctd = dctn - dcpnorm**2
            alpha = (-dctd + DSQRT(dctd**2 - (dnorm**2)*
     $              (dcpnorm**2 - delta**2)))/(dnorm**2)
            CALL dcopy(n, dcp, 1, d, 1)
            CALL dscal(n, 1.0D0-alpha, d, 1)
            CALL daxpy(n, alpha, dn, 1, d, 1)
            dnorm = delta
            iout = 3

          ENDIF

        ELSE

C         The Cauchy step is outside the ball.

          CALL dcopy(n, dcp, 1, d, 1)
          CALL dscal(n, delta/dcpnorm, d, 1)
          dnorm = delta
          iout = 4

        ENDIF

C       Computing the new point, xc, and h(xc).

        CALL dcopy(n, x, 1, xc, 1)
        CALL daxpy(n, 1.0D0, d, 1, xc, 1)
        CALL CuterConstr(m, n, xc, hc)
        nConstr = nConstr + 1
        hcnorm = dnrm2(m, hc, 1)

C       Computing Ared and Pred.

        Ared = hnorm**2 - hcnorm**2
        IF (iout.EQ.2) THEN
          Pred = hnorm**2
        ELSE
          CALL MVprod(m, n, A, d, g, Airow, Apcol, 0.0D0)
          Pred = - dnrm2(m, g, 1)**2 - 2.0D0*ddot(m, h, 1, g, 1)
        ENDIF

      ENDDO

      IF ((Ared.LT.(kappa1*Pred)).AND.(.NOT.FreshA)) THEN

C       The algorithm has failed, A must be recomputed.
  
        CALL dcopy(n, x, 1, xc, 1)
        CALL dcopy(m, h, 1, hc, 1)
        hcnorm = hnorm;
        delta  = deltaold;
        iout   = 5;

      ELSE IF (Ared.GE.(kappa3*Pred)) THEN

C       Increasing Delta.

        delta = DMAX1(kappa4*dnorm, delta)

      ENDIF

      RETURN
      END
