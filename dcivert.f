!     ------------------------------------------------------------------
!     PROGRAM : DCI.
!     MODULE  : DCIVert.
!     VERSION : 2.0.
!     DATE    : June, 2000.
!     CONTENTS: subroutine dcivert
!
!     External routines called by this module:
!
!     From dcicute : Grad, Jacob. 
!     From dcichol : choldec, NAproj.
!     From dcidog  : dcitrust.
!     From the blas: dcopy, dnrm2.
!     ------------------------------------------------------------------

      SUBROUTINE dcivert(m, n, x, h, A, Airow, Apcol, Diag, Mval, Mirow,
     $                   Mpcol, Mcnel, perm, iperm, permavail, DenseRow, 
     $                   lincon, csih, csig, phi1, phi2, nfailv, rhomax, 
     $                   maxrest, hnorm, engp, delta, infdelta, kappa1, 
     $                   kappa2, kappa3, kappa4, thetar, bfgsupd, cls1, 
     $                   cls2, alphab, stepmax, nstpmax, epscg1a, 
     $                   epscg1r, epscg3a, epscg3r, maxitcg, prcndiag, 
     $                   Prec, UseChol, EpsChol, DenseWin, DenseFrc, 
     $                   CholUpd, lbdmax, iter, tmp1, tmp2, tmp3, tmp4, 
     $                   tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, 
     $                   tmp12, tmp13, dn, xc, g, lambda, gp, gpnorm, 
     $                   ngp, lnorm, rho, nRest, nbfgs, nConstr, nGrad, 
     $                   nJacob, flag)
C    $                   , Constr, Grad, Jacob)

C     This routine finds the vertical step. It also computes vectors 
C     g and gp used in the horizontal step. The purpose of the vertical 
C     step is to find a point x such that ||h(x)|| is confined to a 
C     cylinder with a predetermined radius. This is sometimes called a 
C     "restoration". The step is based on the sucessive application of 
C     Powell's dogleg method to solve the trust region problem:
C
C     min  0.5||h(x)+As||^2 - 0.5||h(x)||^2
C     s.t. ||d|| <= delta.
C
C     until ||h(x)|| is sufficiently reduced.
C
C     Input parameters:
C
C     m, n          number of constraints and variables.
C     x             current approximation for the solution.
C     h             constraint vector, h(x).
C     A,Airow,Apcol Jacobian of the constraints.
C     Diag, Mval, Mirow, Mpcol, Mcnel, perm, iperm, DenseRow
C                   Data structure for the Cholesky decompostion of
C                   A*A'. See the dcichol module for a description
C                   of these parameters.
C     permavail     parameter that indicates if the AMD algorithm has
C                   already been used to compute the row and column
C                   permutations of A*A'.
C     lincon        parameter that indicates if all of the constraints
C                   are linear (.TRUE.).
C     csih          tolerance for the infeasibility at the solution.
C     phi1          parameter used to compute rho.
C     phi2          parameter used to call dcibfgs. The BFGS algorithm
C                   is called whenever ||h(xc)|| > phi2*||h(x)|| after 
C                   nfailv iterations (suggestion: 0.95).
C     nfailv        maximum number of unsuccessful iterations. See the 
C                   definition of phi2 above.
C     rhomax        largest admissible value for rho. 
C     maxrest       maximum number of iterations allowed.
C     hnorm         two-norm of h(x).
C     engp          an approximate value for ||gp||/(||g||+1).
C     delta         dcitrust trust region radius.
C     infdelta      smallest admissible value for delta.
C     kappa1        dcitrust parameter used to decide if the step is to 
C                   be accepted (suggestion: 1e-4).
C     kappa2        dcitrust parameter used to reduce delta 
C                   (suggestion: 0.25).
C     kappa3        dcitrust parameter used to decide if delta is to be 
C                   increased (suggestion: 0.7).
C     kappa4        dcitrust parameter used to increase delta 
C                   (suggestion: 2.5).
C     thetar        parameter used to decide if A is to be recomputed 
C                   between two successive dogleg steps. A is rebuild 
C                   whenever hnorm/hnold > thetar (suggestion: 0.9).
C     bfgsupd       parameter that indicates how many BFGS hessian 
C                   updates are to be used to compute the vertical step. 
C                   The BFGS update is not used only if bfgsupd == 0, 
C                   otherwise at most bfgsupd iterations of the BFGS 
C                   algorithm are computed (suggestion: 5).
C     cls1          BFGS parameter used in the sufficient decrease 
C                   condition (suggestion: 1.0e-4). 
C     cls2          BFGS parameter used in the curvature condition 
C                   (suggestion: 0.1).
C     alphab        BFGS parameter used to increase the step in
C                   the extrapolation phase of the algorithm 
C                   (suggestion: 2).
C     stepmax       BFGS parameter. Largest step allowed in the 
C                   line search (suggestion: 5).
C     nstpmax       BFGS parameter. Maximum number of step length 
C                   reductions allowed per call of the steplength 
C                   function (suggestion: 10). 
C     epscg1a       Absolute stopping tolerance used by the CG method 
C                   when computing Lambda.
C     epscg3a       Absolute stopping tolerance used by the CG method 
C                   when computing the Newton step in DCITrust. 
C     epscg1r       Relative stopping tolerance used by the CG method 
C                   when computing Lambda.
C     epscg3r       Relative stopping tolerance used by the CG method 
C                   when computing the Newton step in DCITrust. 
C     maxitcg       Maximum number of iterations of the CG algorithm 
C                   (used to solve the system (A'*A)x = y).
C     prcndiag      Number of (upper and lower) diagonals of the band
C                   preconditioner used when solving (A'*A)x = y by the
C                   CG method (zero means no preconditioning at all). 
C     UseChol       logical parameter that indicates if the Cholesky
C                   decomposition of A*A' is to be computed (.TRUE.)
C                   or if the CG method should be used to solve a 
C                   system with this matrix (.FALSE.).
C     EpsChol       parameter that indicates if A*A' is singular.
C                   A*A' is considered singular if its diagonal
C                   contains an element lower than EpsChol
C                   (suggestion: 1e-16).
C     DenseWin      dimension of the leading submatrix of A*A' that 
C                   is treated as a dense matrix. If DenseWin is 
C                   greater or equal to n, then a dense factorization 
C                   is computed (suggestion: 200).
C     DenseFrc      DenseFrc is the ratio between the number of nonzero 
C                   elements in one row and the size of the matrix that 
C                   is used to define the dense window when computing
C                   the Cholesky decomposition of A*A'. If the number 
C                   of nonzero elements in row k is greater than 
C                   n*DenseFrc, the remaining (n-k+1)x(n-k+1) submatrix 
C                   of A*A' is treated as a dense matrix. If DenseFrc
C                   <= 0, then a dense factorization is computed. If
C                   DenseWin<=1 and DenseFrc>=1, then no dense window 
C                   is generated. See the definition of DenseWin above.
C     CholUpd       logical variable that indicates if the Cholesky
C                   decomposition of A*A' is to be updated (.TRUE.)
C                   instead of recomputed (.FALSE.).
C     lbdmax        upper limit for the Lagrange multipliers.
C     iter          current iteration
C     Constr        a function that computes the constraint vector.
C     Jacob         a function that computes the Jacobian of the 
C                   constraints.
C     Grad          a function that computes the gradient of f(x).
C     nConstr       total number of constraint evaluations.
C     nGrad         total number of gradient evaluations.
C     nJacob        total number of Jacobian evaluations.
C
C     Output parameters:
C
C     xc            new point.
C     h             constraint vector, h(xc).
C     g             gradient of the objective function.
C     lambda        vector of Lagrange multipliers.
C     gp            projected gradient.
C     A             Jacobian of the constraints.
C     Diag, Mval, Mirow, Mpcol, Mcnel, perm, iperm, DenseRow
C                   Data structure for the Cholesky decompostion of
C                   A*A'. See the dcichol module for a description
C                   of these parameters.
C     Prec          Vector that stores the band preconditioner. 
C     hnorm         two-norm of h(xc).
C     gpnorm        ||gp(xc)||.
C     ngp           ||gp(xc)||/(||g(xc)||+1).
C     lnorm         ||lambda||.
C     delta         dcitrust trust region radius.
C     rho           trust cylinder radius.
C     nRest         number of restoration steps.
C     nbfgs         number of BFGS iterations.
C     nConstr       see description above.
C     nGrad         see description above.
C     nJacob        see description above.
C     flag          Output flag. 
C                     0: the vertical step was computed successfully.
C                     1: an error has occurred.
C
C     Temporary vectors passed as parameters:
C
C     tmp1(m), tmp2(m),
C     tmp3(max(m, bfgsupd)), tmp4(max(m, bfgsupd)),
C     dn(n), tmp5(n), tmp6(n), tmp7(n), tmp8(n), tmp9(m), tmp10(n),
C     tmp11(n*(bfgsupd+1)), tmp12(n*(bfgsupd+1)), tmp13(n*(bfgsupd+1)).
                         
      IMPLICIT NONE

      EXTERNAL Constr, Grad, Jacob

C     Subroutine parameters.

      INTEGER m, n, maxrest, maxitcg, prcndiag, DenseRow
      INTEGER nfailv, bfgsupd, nstpmax, DenseWin, iter
      INTEGER nRest, nbfgs, nConstr, nGrad, nJacob, flag
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      INTEGER Mirow(*)
      INTEGER Mpcol(m+1)
      INTEGER Mcnel(m), perm(m), iperm(m)
      REAL*8  csih, csig, phi1, phi2, rhomax, hnorm, gpnorm, ngp
      REAL*8  lnorm, delta, engp, infdelta, kappa1, kappa2, kappa3
      REAL*8  kappa4, rho, epscg1a, epscg1r, epscg3a, epscg3r, thetar
      REAL*8  cls1, cls2, alphab, stepmax, EpsChol, DenseFrc, lbdmax
      REAL*8  x(n), xc(n), g(n), gp(*)
      REAL*8  dn(n), tmp5(n), tmp6(n), tmp7(n), tmp8(n), tmp10(n)
      REAL*8  tmp1(m), tmp2(m), tmp3(m), tmp4(m), tmp9(m)
      REAL*8  tmp11(n*(bfgsupd+1)), tmp12(n*(bfgsupd+1))
      REAL*8  tmp13(n*(bfgsupd+1))
      REAL*8  h(m), lambda(m), Diag(m)
      REAL*8  A(*)
      REAL*8  Mval(*)
      REAL*8  Prec(*)
      LOGICAL lincon, UseChol, CholUpd, permavail

C     Local variables.

      INTEGER naflag, trflag, oldAcnt, fail, iout, ibfgs
      REAL*8  rcyl, hnold, epscg, mu, dnnorm
      LOGICAL Aavail, gavail, LimLbd, dnavail

C     Functions called by the routine.

      REAL*8  dnrm2

      ngp = engp
      rho = DMIN1(phi1*rhomax*ngp, 0.75D0*rhomax)
      IF ((rho.LT.csih).AND.(hnorm.GT.(100.0D0*csih))) THEN
        rho = 0.1D0*hnorm;
      ENDIF
      nRest = 0
      nbfgs = 0
      hnold = hnorm
      CALL dcopy(n, x, 1, xc, 1)
      oldAcnt = 1
      flag = 0
      fail = 0
      IF (iter.EQ.1) THEN
        Aavail = .TRUE.
        gavail = .TRUE.
      ELSE
        Aavail = .FALSE.
        gavail = .FALSE.
      ENDIF

      IF ((hnorm.LE.rho).AND.(.NOT.Aavail)) THEN
C INSERIR NOVO CRITÉRIO AQUI (USANDO UM EXPOENTE)

        IF (m.GT.0) THEN

C         Recomputing g and A.

          IF (.NOT.lincon) THEN
            CALL CuterJacob(m, n, xc, A, Airow, Apcol)
            nJacob = nJacob + 1
          ENDIF
          Aavail = .TRUE.

          CALL CuterGrad(n, xc, g)
          nGrad = nGrad + 1
          gavail = .TRUE.

          IF ((UseChol).AND.(.NOT.lincon)) THEN

C           Computing the Cholesky decomposition of A*A'.

            CALL choldec(m, n, A, Airow, Apcol, Diag, Mval, Mirow, 
     $                   Mpcol, Mcnel, perm, iperm, EpsChol, DenseWin,
     $                   DenseFrc, DenseRow, CholUpd, permavail)

C           Computing gp and lambda.

            LimLbd = .TRUE.
            CALL NAprojCh(m, n, A, Airow, Apcol, Diag, Mval, Mirow, 
     $                    Mpcol, Mcnel, perm, iperm, DenseRow, lambda, 
     $                    LimLbd, lbdmax, g, gp, naflag)

          ELSE

C           Computing gp and lambda.

            epscg = DMAX1(epscg1a, epscg1r*gpnorm)
            mu = 0.0D0
            LimLbd = .TRUE.
            CALL NAprojCG(m, n, prcndiag, .TRUE., maxitcg, epscg,
     $                    mu, A, Airow, Apcol, Prec, tmp1, tmp2, 
     $                    tmp3, tmp4, tmp5, lambda, LimLbd, lbdmax,
     $                    g, gp, naflag)

          ENDIF
 
c$$$          IF (naflag.NE.0) THEN
c$$$            WRITE(*,*)'DCIVert warning: NAproj has failed.'
c$$$          ENDIF

        ELSE

          CALL CuterGrad(n, xc, g)
          nGrad = nGrad + 1

          CALL dcopy(n, g, 1, gp, 1)

        ENDIF

C       Recomputing ngp and rho.

        gpnorm = dnrm2(n, gp, 1)
        ngp = gpnorm/(dnrm2(n, g, 1) + 1.0D0)
        rho = DMIN1(phi1*rhomax*ngp, 0.75D0*rhomax);
        IF ((rho.LT.csih).AND.(hnorm.GT.(100.0D0*csih))) THEN
          rho = 0.1D0*hnorm;
        ELSE IF (ngp.LE.(5.0D0*csig)) THEN
          rho = csih
        ENDIF

        oldAcnt = 0

      ENDIF

C     Main loop.

      DO WHILE ((hnorm.GT.rho).AND.(nrest.LE.maxrest).AND.(flag.EQ.0))

        DO WHILE ((hnorm.GT.rho).AND.(nrest.LE.maxrest)
     &            .AND.(flag.EQ.0))

C         Computing the new point.

          nrest = nrest + 1

          IF (.NOT.UseChol) THEN
             epscg = DMAX1(epscg3a, epscg3r*DMIN1(1.0D0, hnorm))
          ENDIF

          CALL dcitrust(n, m, xc, hnold, A, Airow, Apcol, Aavail, delta,
     &                  infdelta, kappa1, kappa2, kappa3, kappa4,
     &                  UseChol, epscg, maxitcg, prcndiag, Prec, Diag,  
     &                  Mval, Mirow, Mpcol, Mcnel, perm, iperm,
     &                  DenseRow, tmp1, tmp2, tmp3, tmp5, tmp6, tmp7, 
     &                  tmp8, tmp9, tmp10, h, hnorm, nConstr, trflag) 
C    &                  , Constr)
          gavail = .FALSE.


C         Counting the number of successive unsuccessful steps.

          IF (hnorm.GT.(phi2*hnold)) THEN
            fail = fail + 1;
          ELSE
            fail = 0;
          ENDIF

          IF (fail.GE.nfailv) THEN

C           The trust region algorithm is not working. 
C           Calling the BFGS routine.

            IF (.NOT.lincon) THEN
              CALL CuterJacob(m, n, xc, A, Airow, Apcol)
              nJacob = nJacob + 1
            ENDIF
            Aavail = .TRUE.

            IF ((UseChol).AND.(.NOT.lincon)) THEN

              CALL choldec(m, n, A, Airow, Apcol, Diag, Mval, Mirow, 
     $                     Mpcol, Mcnel, perm, iperm, EpsChol, DenseWin,
     $                     DenseFrc, DenseRow, CholUpd, permavail)
              CALL NewStepCh(m, n, A, Airow, Apcol, Diag, Mval, Mirow, 
     $                       Mpcol, Mcnel, perm, iperm, DenseRow, tmp1,
     $                       h, dn, naflag)

            ELSE

              mu = 0.0D0
              CALL NewStepCG(m, n, prcndiag, maxitcg, epscg, mu,
     $                       A, Airow, Apcol, Prec, tmp1, tmp2,
     $                       tmp3, tmp4, tmp5, tmp9, h, dn, naflag)

            ENDIF

            dnnorm = dnrm2(n, dn, 1)
            dnavail = .TRUE.
            ibfgs = 0
            CALL dcibfgs(m, n, xc, h, A, Airow, Apcol, hnorm, dnavail, 
     $                   dn, dnnorm, rho, bfgsupd, cls1, cls2, alphab, 
     $                   stepmax, nstpmax, tmp5, tmp6, tmp7, tmp8,
     $                   tmp10, tmp3, tmp4, tmp11, tmp12, tmp13, ibfgs, 
     $                   nConstr, nJacob, iout)
C    $                   , Constr, Jacob);            
            nbfgs = nbfgs + ibfgs

            Aavail = .TRUE.
            IF ((UseChol).AND.(.NOT.lincon)) THEN
              CALL choldec(m, n, A, Airow, Apcol, Diag, Mval, Mirow, 
     $                     Mpcol, Mcnel, perm, iperm, EpsChol, DenseWin,
     $                     DenseFrc, DenseRow, CholUpd, permavail)
            ENDIF

            IF ((iout.EQ.1).OR.(iout.GT.2)) THEN

C             The restoration has failed.

              flag = 1

            ENDIF

          ELSE IF (((hnorm.GT.(thetar*hnold)).AND.(oldAcnt.GT.0)).OR.
     $             (oldAcnt.GT.5).OR.(iout.EQ.5)) THEN

C           Recomputing A, L, and perm because hnorm was 
C           not sufficiently reduced.

            IF (.NOT.lincon) THEN
              CALL CuterJacob(m, n, xc, A, Airow, Apcol)
              nJacob = nJacob + 1
            ENDIF
            Aavail = .TRUE.
            oldAcnt = 0

            IF ((UseChol).AND.(.NOT.lincon)) THEN
              CALL choldec(m, n, A, Airow, Apcol, Diag, Mval, Mirow, 
     $                     Mpcol, Mcnel, perm, iperm, EpsChol, DenseWin,
     $                     DenseFrc, DenseRow, CholUpd, permavail)
            ENDIF

          ELSE

C           Trying to use the same A once more.

            oldAcnt = oldAcnt + 1
            Aavail = .FALSE.

          ENDIF  

          hnold = hnorm
          delta = DMAX1(delta, infdelta)

        ENDDO

C       Recomputing A and its Cholesky decomposition, if necessary.

        IF (.NOT.Aavail) THEN
          IF (.NOT.lincon) THEN
            CALL CuterJacob(m, n, xc, A, Airow, Apcol)
            nJacob = nJacob + 1
          ENDIF
          Aavail = .TRUE.
          IF ((UseChol).AND.(.NOT.lincon)) THEN
            CALL choldec(m, n, A, Airow, Apcol, Diag, Mval, Mirow, 
     $                   Mpcol, Mcnel, perm, iperm, EpsChol, DenseWin,
     $                   DenseFrc, DenseRow, CholUpd, permavail)
          ENDIF
        ENDIF

        IF (.NOT.gavail) THEN

C         Computing g.

          CALL CuterGrad(n, xc, g)
          nGrad = nGrad + 1

C         Computing gp and lambda.

          IF (UseChol) THEN
            LimLbd = .TRUE.
            CALL NAprojCh(m, n, A, Airow, Apcol, Diag, Mval, Mirow, 
     $                    Mpcol, Mcnel, perm, iperm, DenseRow, lambda, 
     $                    LimLbd, lbdmax, g, gp, naflag)
          ELSE
            epscg = DMAX1(epscg1a, epscg1r*gpnorm)
            mu = 0.0D0
            LimLbd = .TRUE.
            CALL NAprojCG(m, n, prcndiag, .TRUE., maxitcg, epscg,
     $                    mu, A, Airow, Apcol, Prec, tmp1, tmp2, 
     $                    tmp3, tmp4, tmp5, lambda, LimLbd, lbdmax,
     $                    g, gp, naflag)
          ENDIF

c$$$         IF (naflag.NE.0) THEN
c$$$           WRITE(*,*)'DCIVert warning: NAproj has failed.'
c$$$         ENDIF

C         Computing ngp.

          gpnorm = dnrm2(n, gp, 1)
          ngp = gpnorm/(dnrm2(n, g, 1) + 1.0D0)

        ENDIF

C       Computing rho and the two-norms of Lambda and sn.

        IF (rho.GT.(2.0D0*rhomax*ngp)) THEN
          rho = DMIN1(phi1*rhomax*ngp, 
     $          DMAX1(1.0D-4*rhomax*ngp, 0.75D0*rhomax))
        ELSE
          rho = DMAX1(rho, DMIN1(phi1*rhomax*ngp, 
     $          DMAX1(1.0D-4*rhomax*ngp, 0.75D0*rhomax)))
        ENDIF

      ENDDO

      lnorm = dnrm2(m, lambda, 1)
      CALL dcopy(n, xc, 1, tmp5, 1)
      CALL daxpy(n, -1.0D0, x, 1, tmp5, 1)

      RETURN
      END
