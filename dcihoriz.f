!     ------------------------------------------------------------------
!     PROGRAM : DCI.
!     MODULE  : DCIHoriz.
!     VERSION : 2.0.
!     DATE    : June, 2000.
!     CONTENTS: subroutine dcihoriz
!
!     External routines called by this module:
!
!     From dcicute : Fun, Constr, Hprod. 
!     From dcichol : NewStep.
!     From dcisteih: dcisteih.
!     From the blas: dcopy, daxpy, ddot, dnrm2.
!     ------------------------------------------------------------------

      SUBROUTINE dcihoriz(m, n, xc, lambda, h, g, gp, A, Airow, Apcol, 
     *                    UseChol, Diag, Mval, Mirow, Mpcol, Mcnel, 
     *                    perm, iperm, DenseRow, Lc, hnorm, delta, csih,
     *                    eps1, eps2, eps3, maxit, epscg2a, epscg2r, 
     *                    epscg4a, epscg4r, maxitcg, prcndiag, rho, 
     *                    zeta1, zeta2, zeta3, alphar, alphai, alphas, 
     *                    eta1, eta2, eta3, xnew, fnew, hnew, Lnew, 
     *                    snorm, DLh, Prec, tmp, tmp1, tmp2, tmp3, tmp4, 
     *                    tmp5, tmp6, tmp7, tmp8, tmp9, ssoc, s, 
     *                    nSteih, nrej, nFun, nConstr, nHprod, nSoc)
C    *                    , Fun, Constr, Hprod)

C     This routine computes the horizontal step by using Steihaug's 
C     algorithm. The horizontal step is an approximate solution for:
C
C     min  0.5*s'*H*s + g'*s
C     s.t. A*s = 0
C          ||s|| <= delta.
C
C     Input parameters:
C
C     m, n          number of constraints and variables.
C     xc            current approximation for the solution.
C     lambda        current Lagrange multipliers.
C     h             h(xc).
C     g             gradient of the Lagrangian.
C     gp            projected gradient of the Lagrangian.
C     A,Airow,Apcol Jacobian of the constraints.
C     UseChol       logical parameter that indicates if the Cholesky
C                   decomposition of A*A' is to be computed (.TRUE.)
C                   or if the CG method should be used to solve a 
C                   system with this matrix (.FALSE.).
C     Diag, Mval, Mirow, Mpcol, Mcnel, perm, iperm, DenseRow
C                   Data structure for the Cholesky decompostion of
C                   A*A'. See the dcichol module for a description
C                   of these parameters.
C     Lc            Lagrangian at (xc, lambda).
C     hnorm         two-norm of h(xc).
C     delta         trust region radius.
C     csih          Stopping tolerance for the infeasibility.
C     eps1          relative stopping tolerance. The algorithm stops 
C                   when r'*v <= eps1*(r0'*v0). Suggestion: 0.01.
C     eps2          absolute stopping tolerance. The algorithm stops 
C                   when r'*v <= eps2. Suggestion: 1e-10.
C     eps3          relative tolerance used in the Steihaug's algorithm
C                   to identify a direction of negative curvature.
C                   Suggestion: 1e-8.
C     maxit         maximum no. of iterations of Steihaug's algorithm.
C     epscg2a       Absolute stopping tolerance used by the CG method 
C                   when projecting the residual in DCISteih. 
C     epscg4a       Absolute stopping tolerance used by the CG method 
C                   when computing the second order correction.
C     epscg2r       Relative stopping tolerance used by the CG method 
C                   when projecting the residual in DCISteih. 
C     epscg4r       Relative stopping tolerance used by the CG method 
C                   when computing the second order correction.
C     maxitcg       Maximum number of iterations of the CG algorithm 
C                   (used to solve the system (A'*A)x = y).
C     prcndiag      Number of (upper and lower) diagonals of the band
C                   preconditioner used when solving (A'*A)x = y by the
C                   CG method (zero means no preconditioning at all). 
C     Prec          Vector that stores the band preconditioner. 
C     rho           radius of the trust region cylinder.
C     zeta1         parameter used to decide if the step is to be
C                   accepted and also if a 2nd. order correction
C                   is to be added to the step (suggestion: 2.0).
C     zeta2         parameter used to decide if a 2nd. order correction
C                   is to be added to the step (suggestion: 1.0).
C     zeta3         parameter used to decide if a 2nd. order correction
C                   is to be added to the step (suggestion: 0.1).
C     alphar        parameter used to reduce delta (suggestion: 0.5).
C     alphai        parameter used to increase delta (suggestion: 2).
C     alphas        parameter used to reduce delta when Ared/Pred < 0 
C                   (suggestion: 0.0625).
C     eta1          parameter used to decide if delta is to be 
C                   reduced (suggestion: 1e-4).
C     eta2          parameter used to decide if delta is to be 
C                   increased (suggestion: 0.5).
C     eta3          parameter used to compute the reduction of delta 
C                   when AredO/AredF < 0 (suggestion: 0.1).
C     Fun           A function that computes the objective function.
C     Constr        A function that computes the constraint vector.
C     HProd         A function that computes the product of the Hessian
C                   of the lagrangian by a vector.
C     nFun          total number of function evaluations. 
C     nConstr       total number of constraint evaluations.
C     nHprod        total number of Hessian-vector products.
C
C     Output parameters:
C
C     xnew          new point.
C     fnew          f(xnew).
C     hnew          h(xnew).
C     Lnew          Lagrangian at (xnew, lambda).
C     hnorm         two-norm of h(xnew).
C     DLh           Lnew - Lc.
C     nSteih        total number of CG iterations.
C     snorm         two-norm of the step.
C     nRej          number of rejected steps.
C     nFun          see description above.
C     nConstr       see description above.
C     nHprod        see description above.
C     nSoc          number of second order corrections.
C
C     Temporary vectors passed as parameters:
C
C     tmp, tmp1, tmp2, tmp3, tmp4, tmp5, 
C     tmp6, tmp7, tmp8, tmp9, s, ssoc.

      IMPLICIT NONE

      EXTERNAL Fun, Constr, Hprod

C     Subroutine parameters.

      INTEGER m, n, maxit, nRej, nFun, nConstr, nHprod, nSteih, nSoc
      INTEGER maxitcg, prcndiag, DenseRow
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      INTEGER Mirow(*)
      INTEGER Mpcol(m+1)
      INTEGER Mcnel(m), perm(m), iperm(m)
      REAL*8  Lc, hnorm, delta, csih, eps1, eps2, eps3, rho 
      REAL*8  epscg2a, epscg2r, epscg4a, epscg4r, zeta1, zeta2, zeta3
      REAL*8  alphar, alphai, alphas, eta1, eta2, eta3
      REAL*8  fnew, Lnew, hnnew, DLh, snorm
      REAL*8  xc(n), g(n), gp(n), xnew(n), s(n), ssoc(n)
      REAL*8  tmp5(n), tmp6(n), tmp7(n), tmp8(n), tmp9(n)
      REAL*8  tmp(m), tmp1(m), tmp2(m), tmp3(m), tmp4(m)
      REAL*8  lambda(m), h(m), hnew(m), Diag(m)
      REAL*8  Prec(*)
      REAL*8  A(*)
      REAL*8  Mval(*)
      LOGICAL UseChol

C     Internal variables.

      INTEGER iter, flag
      REAL*8  qs, epscg, gts, asoc, alphat
      LOGICAL GotH, first

C     Functions called by the routine.

      INTEGER idamax
      REAL*8  CuterFun, ddot, dnrm2

C     Setting up a few parameters.

      GotH   = .FALSE.
      first  = .TRUE.
      nSteih = 0
      nRej   = 0
      DLh    = -1.0D0
      qs     = 2.0D0*DLh/eta1
      hnnew  = 2.0D0*rho + 1.0D0
      nSoc   = 0

C     Main loop. (Notice that DLh and qs are negative!)                 

      DO WHILE ((hnnew.GT.(zeta1*rho)).OR.(DLh.GT.(eta1*qs)))

C       Computing the horizontal step by Steihaug's algorithm
C       (ssoc is used as a temporary vector).

        CALL dcisteih(m, n, xc, lambda, g, gp, A, Airow, Apcol, UseChol, 
     $                Diag, Mval, Mirow, Mpcol, Mcnel, perm, iperm, 
     $                DenseRow, delta, eps1, eps2, eps3, maxit, epscg2a, 
     $                epscg2r, maxitcg, prcndiag, s, qs, gts, iter, 
     $                flag, Prec, tmp, tmp1, tmp2, tmp3, tmp4, tmp5, 
     $                tmp6, tmp7, tmp8, tmp9, ssoc, nHprod, GotH)
C    $                , Hprod)
        nSteih = nSteih + iter

c$$$        IF (flag.EQ.3) THEN
c$$$          WRITE(*,*)'dcisteih warning: the maximum number of',
c$$$     $              ' iterations was taken' 
c$$$        ENDIF

C       Computing xnew and h(xnew).

        CALL dcopy(n, xc, 1, xnew, 1)
        CALL daxpy(n, 1.0D0, s, 1, xnew, 1)
        CALL CuterConstr(m, n, xnew, hnew)
        nConstr = nConstr + 1
        hnnew = dnrm2(m, hnew, 1)
        IF ((flag.EQ.1).OR.(flag.EQ.2)) THEN
          snorm = delta
        ELSE
          snorm = dnrm2(n, s, 1)
        ENDIF

C       Computing a second order correction.

!       IF (((hnnew.GT.DMIN1(0.9D0*zeta1*rho,zeta1*hnorm+zeta2*rho)).OR.
        IF (((hnnew.GT.DMIN1(zeta1*rho,zeta1*hnorm+zeta2*rho)).OR.
     $     ((hnorm.LE.csih).AND.(hnnew.GT.DMAX1(csih, (zeta3*hnorm)))))
     $     .AND.first) THEN
 
C         Second order correction.

          IF (UseChol) THEN
            CALL NewStepCh(m, n, A, Airow, Apcol, Diag, Mval, Mirow,
     $                     Mpcol, Mcnel, perm, iperm, DenseRow, tmp1,
     $                     hnew, ssoc, flag)
          ELSE
            epscg = DMAX1(epscg4a, epscg4r*DMIN1(1.0D0, hnnew))
            CALL NewStepCG(m, n, prcndiag, maxitcg, epscg, 0.0D0,
     $                     A, Airow, Apcol, Prec, tmp1, tmp2,  
     $                     tmp3, tmp4, tmp5, tmp, hnew, ssoc, flag)
          ENDIF

c$$$          IF (flag.NE.0) THEN
c$$$            WRITE(*,*)'DCIHoriz warning: NewStep has failed.'
c$$$          ENDIF

          CALL daxpy(n, 1.0D0, s, 1, ssoc, 1)
          asoc = DABS(ssoc(idamax(n, ssoc, 1)))
          IF (asoc.GT.delta) THEN
            CALL dscal(n, delta/asoc, ssoc, 1)
          ENDIF
          CALL dcopy(n, xc, 1, xnew, 1)
          CALL daxpy(n, 1.0D0, ssoc, 1, xnew, 1)
          CALL CuterConstr(m, n, xnew, hnew)
          nConstr = nConstr + 1
          hnnew = dnrm2(m, hnew, 1)
          nSoc = nSoc + 1
 
        ENDIF

C       Computing the expected reduction of the Lagrangian.

        fnew = CuterFun(n, xnew)
        nFun = nFun + 1
        Lnew = fnew + ddot(m, lambda, 1, hnew, 1)
        DLh = Lnew - Lc

C       Adjusting the trust region radius.
C       (Notice that DLh and qs are negative!)

        IF ((hnnew.GT.(zeta1*rho)).OR.(DLh.GT.(eta1*qs))) THEN
          IF (DLh.GT.0) THEN 
C           Reducing delta drastically because DLh is negative.
            alphat = (1.0D0 - eta3)*gts/
     $               ((1.0D0 - eta3)*(Lc + gts) - Lnew + eta2*qs)
            delta = DMIN1(alphar*snorm, DMAX1(alphat, alphas)*delta)
          ELSE 
C           Reducing delta.
            delta = snorm*alphar;
          ENDIF
          nRej = nRej + 1
        ELSE IF (DLh.LE.(eta2*qs)) THEN
C         Increasing delta.
          delta = DMAX1(delta, snorm*alphai)
        ENDIF

        first = .FALSE.

      ENDDO

      hnorm = hnnew

      RETURN
      END
