!     ------------------------------------------------------------------
!     PROGRAM : DCI.
!     MODULE  : DCI.
!     VERSION : 2.1.
!     DATE    : March, 2006.
!     CONTENTS: subroutine DCI
!
!     External routines called by this module:
!
!     From dcicute : Fun, Constr, Grad, Jacob, Hprod.
!     From dciio   : ReadParameters, InitialValues, KPrint, LastPrint,
!                    dbgprnt
!     From dcivert : dcivert.
!     From dcihoriz: dcihoriz.
!     From the blas: ddot.
!     ------------------------------------------------------------------


      SUBROUTINE DCI(m, n, lincon, x, lambda, fc, h, iter, ierr)
C    $               Fun, Constr, Grad, Jacob, Hprod)

C     1) Purpose:
C
C     The purpose of this subroutine is to solve the equality
C     constrained nonconvex nonlinear programming problem
C
C
C                          minimize   f(x)
C                          subject to h(x)=0
C
C     where x \in R^n, f:R^n->R, h:R^n-R^m, m<=n, and f and h are
C     twice continuously differentiable.
C
C
C     2) Input parameters.
C
C     a) Functions and subroutines.
C
C     REAL*8 FUNCTION Fun(n,x). This function must return f(x). n
C     is the dimension of x.
C
C     SUBROUTINE Constr(m,n,x,h). This subroutine must return a 
C     m-dimensional vector h that contains h(x). The dimension of
C     x is given by n.
C
C     SUBROUTINE Grad(n,x,g). This subroutine must return g, a n-
C     dimensional vector that contains the gradient of f at x.
C
C     SUBROUTINE Jacob(m,n,x,A,Airow,Apcol). This subroutine must
C     return the triple {A, Airow, Apcol} containing the nonzero
C     elements of the Jacobian of the constraints. m is the number
C     of constraints and n is the number of variables, so the
C     Jacobian is a mxn real matrix.
C     
C     SUBROUTINE Hprod(m,n,x,lambda,u,v). This subroutine must
C     return a vector v = H.u, where H is the Hessian of the
C     Lagrangian of the problem. The n components of the current
C     approximate solution are stored in x, while the m Lagrange
C     multipliers are given by lambda. 
C
C     b) Integer parameters:
C
C     n       dimension of x.
C     m       dimension of h(x).
C
C     c) Logical parameter:
C
C     lincon  parameter that indicates if all of the constraints
C             are linear or affine (.TRUE.).
C
C     d) Real vectors:
C
C     x       initial guess for the solution of the problem.
C     lambda  initial approximation for the Lagrange multipliers.
C
C                               
C     3) Output parameters.
C
C     a) Integer parameters:
C
C     iter    number of iterations taken.       
C     ierr    exit code.
C             0: the projected gradien is small enough.
C             1: the maximum trust cylinder radius is small enough.
C             2: the maximum number of iterations was reached;             
C             3: the maximum number of restorations was reached.
C             4: the restoration has failed (delta is too small).
C             5: the norm of the step was smaller than minstep*
C                xnorm for maxssmll successive iterations.
C
C     b) Real parameters:
C
C     fc      f(x).
C
C     c) Real vectors:
C
C     x       approximate solution.        
C     h       h(x).
C     lambda  Lagrange multipliers.
C
C  
C     5) Stopping criteria.
C
C     This routine stops when one of the following criteria is
C     satisfied (please refer to the "Internal variables" section
C     below for a description of the variables used here):
C
C     a) ||gp(x)|| < csig and ||h(x)|| < csih;
C     b) rhomax < rhomin;
C     c) iter > maxit;
C     d) nrest > maxrest;
C     e) vdelta < infdelta;
C     e) snorm < minstep*xnorm for maxssmll successive iterations.
C
C     6) Parameters read from disk.
C
C     a) Integer parameters.
C
C     maxrest  maximum number of restorations allowed.
C     maxitst  maximum number of Steihaug iterations allowed.
C              See also the relitst parameter below.
C     minitst  Lower limit for maxitst. See the relitst parameter below.
C     maxitcg  Maximum number of iterations of the CG algorithm 
C              (used to solve the system (A'*A)x = y).
C              See also the relitcg parameter below.
C     minitcg  Lower limit for maxitcg (see the relitcg parameter below).
C     maxit    maximum number of iterations allowed.
C     maxssmll maximum number of iterations with a small step.
C              Suggestion: 5.
C     prcndiag Number of diagonals of A'*A included in the preconditioner
C              used to solve the system (A'*A)x = y by the CG method.
C     prnt     a variable that indicates how frequently the
C              subroutine will print data about the iteration.
C              If prnt < 0, the subroutine prints nothing;
C              if prnt = 0, only a summary is printed at the end
C                           of the subroutine;
C              if prnt > 0, a summary is printed after prnt
C                           iterations.
C     solprnt  a variable that indicates if the solution of
C              the problem is to be saved.
C              If solprnt < 0, the subroutine saves nothing;
C              if solprnt = 0, only a summary is saved at the end
C                           of the subroutine;
C              if solprnt > 0, a summary and the values of all of
C                           the variables and constraints are saved.
C     tabprnt  a variable that indicates if a one-line summary
C                           about the solution is to be saved.
C              If solprnt < 0, the subroutine saves nothing;
C              if solprnt >=0, the ole-line summary is saved.
C     nfailv   maximum number of unsuccessful vertical iterations. 
C              See the definition of phi2 below.
C     bfgsupd  dcivert parameter that indicates how many BFGS hessian 
C              updates are to be used to compute the vertical step. 
C              The BFGS update is not used only if bfgsupd == 0, 
C              otherwise at most bfgsupd iterations of the BFGS 
C              algorithm are computed (suggestion: 5).
C     nstpmaxb dcivert parameter. Maximum number of step length 
C              reductions allowed per call of the steplength 
C              function inside the BFGS algorithm (suggestion: 10). 
C     DenseWin dimension of the leading submatrix of A*A' that 
C              is treated as a dense matrix when the Cholesky 
C              decomposition of A*A' is computed. If DenseWin is 
C              greater or equal to n, then a dense factorization 
C              is used (suggestion: 200).
C
C     b) Real parameters.
C
C     relitst  maximum number of Steihaug iterations relative to the
C              size of the problem. The maximum number of iterations 
C              is set to max(minitst, min(maxitst, n*relitst).
C     relitcg  Maximum number of iterations of the CG algorithm 
C              relative to the size of the constraint vector.
C              The maximum number of CG iterations is set to 
C              max(minitcg, min(maxitcg, m*relitcg).
C     csih     tolerance for the infeasibility at the solution.
C     csig     tolerance for the projected gradient.
C     rhomin   tolerance for the size of rhomax.
C     phi1     dcivert parameter used to compute rho.
C     phi2     dcivert parameter used to call dcibfgs. The BFGS 
C              algorithm is called whenever ||h(xc)|| > phi2*||h(x)|| 
C              after nfailv iterations of the usual vertical step
C              (suggestion: 0.95).
C     kappa1   dcitrust parameter used to decide if the step is to be 
C              accepted (suggestion: 1e-4).
C     kappa2   dcitrust parameter used to reduce vdelta 
C              (suggestion: 0.25).
C     kappa3   dcitrust parameter used to decide if vdelta is to be 
C              increased (suggestion: 0.7).
C     kappa4   dcitrust parameter used to increase vdelta 
C              (suggestion: 2.5).  
C     epsst1   relative stopping tolerance for Steihaug's algorithm. 
C              Suggestion: 0.01.
C     epsst2   absolute stopping tolerance for Steihaug's algorithm. 
C              The value defined by the user is multiplied by csig 
C              at the beginning of the algorithm. Suggestion: 1e-6.
C     epsst3   relative tolerance used in the Steihaug's algorithm
C              to identify a direction of negative curvature.
C              Suggestion: 1e-8.
C     epscg1a  Absolute stopping tolerance used by the CG method when 
C              computing Lambda in DCIVert. The value defined by the 
C              user is multiplied by csig at the beginning of the
C              algorithm. Suggestion: 1e-2.
C     epscg2a  Absolute stopping tolerance used by the CG method when 
C              projecting the residual in DCISteih. The value defined 
C              by the user is multiplied by csih at the beginning of 
C              the algorithm. Suggestion: 1e-5.
C     epscg3a  Absolute stopping tolerance used by the CG method when 
C              computing the Newton step in DCITrust. The value 
C              defined by the user is multiplied by csih at the 
C              beginning of the algorithm. Suggestion: 1e-2.
C     epscg4a  Absolute stopping tolerance used by the CG method when 
C              computing the second order correction in DCIHoriz. 
C              The value defined by the user is multiplied by csih at 
C              the beginning of the algorithm. Suggestion: 1e-2.
C     epscg1r  Relative stopping tolerance used by the CG method when 
C              computing Lambda in DCIVert. Suggestion: 1e-5.
C     epscg2r  Relative stopping tolerance used by the CG method when 
C              projecting the residual in DCISteih. Suggestion: 1e-8.
C     epscg3r  Relative stopping tolerance used by the CG method when 
C              computing the Newton step in DCITrust. Suggestion: 1e-3.
C     epscg4r  Relative stopping tolerance used by the CG method when 
C              computing the second order correction in DCIHoriz. 
C              Suggestion: 1e-5.
C     zeta1    dcihoriz parameter used to decide if the step is to be
C              accepted (hnnew.LE.(zeta1*rho)) and also if a 2nd. order 
C              correction is to be added to the step (suggestion: 2.0).
C     zeta2    dcihoriz parameter used to decide if a 2nd. order 
C              correction is to be added to the step (suggestion: 1.0).
C     zeta3    dcihoriz parameter used to decide if a 2nd. order
C              correction is to be added to the step (suggestion: 0.1).
C     alphar   dcihoriz parameter used to reduce hdelta 
C              (suggestion: 0.5).
C     alphai   dcihoriz parameter used to increase hdelta 
C              (suggestion: 2).
C     alphas   dcihoriz parameter used to reduce delta when 
C              Ared/Pred < 0 (suggestion: 0.0625).
C     alphab   dcivert parameter used to increase the step in
C              the extrapolation phase of the BFGS algorithm 
C              (suggestion: 2).
C     eta1     dcihoriz parameter used to decide if hdelta is to be 
C              reduced (suggestion: 1e-4).
C     eta2     dcihoriz parameter used to decide if hdelta is to be 
C              increased (suggestion: 0.7).
C     eta3     dcihoriz parameter used to compute the reduction of 
C              delta when AredO/AredF < 0 (suggestion: 0.1).
C     delta0   lower limit for vdelta and hdelta at iteration 0.
C     mindelta smallest value allowed for hdelta and vdelta at
C              the beginning of an iteration.
C     infdelta smallest value for delta. If delta < infdelta, the
C              algorithm is stopped with ierr = 4. The value defined 
C              by the user is multiplied by csih at the beginning of 
C              the algorithm. Suggestion: 1e-10.
C     minstep  smallest step size allowed (see the ierr variable above).
C              The value defined by the user is multiplied by 
C              min(csih,csig) at the beginning of the algorithm. 
C              Suggestion: 1e-3.
C     thetar   dcivert parameter used to decide if A is to be recomputed 
C              between two successive dogleg steps. A is rebuild 
C              whenever hnorm/hnold > thetar (suggestion: 0.9).
C     cls1     dcivert parameter used in the sufficient decrease 
C              condition of the BFGS algorithm (suggestion: 1.0e-4). 
C     cls2     dcivert parameter used in the curvature condition of the
C              BFGS algorithm (suggestion: 0.1).
C     stepmaxb dcivert parameter. Largest step allowed in the line
C              search done within the BFGS algorithm (suggestion: 5).
C     EpsChol  dcivert parameter that indicates if A*A' is singular.
C              A*A' is considered singular if its diagonal contains an 
C              element lower than EpsChol (suggestion: 1e-16).
C     DenseFrc dcivert parameter. Ratio between the number of nonzero 
C              elements in one row and the size of the matrix that is 
C              used to define the dense window when computing the
C              Cholesky decomposition of A*A'. If the number of nonzero
C              elements in row k is greater than n*DenseFrc, the 
C              remaining (n-k+1)x(n-k+1) submatrix of A*A' is treated 
C              as a dense matrix. If DenseFrc <= 0, then a dense 
C              factorization is computed. If DenseWin <= 1 and DenseFrc
C              >= 1, then no dense window is generated. See the 
C              definition of DenseWin above.
C     lbdmax   upper limit for the Lagrange multipliers.
C
C     b) Logical parameters.
C
C     UseChol  parameter that indicates if the Cholesky decomposition of
C              A*A' is to be computed (TRUE) or if the CG method should
C              be used to solve a system with this matrix (FALSE).
C     UpdChol  logical variable that indicates if the Cholesky
C              decomposition of A*A' is to be updated (TRUE) instead
C              of recomputed (FALSE).
C
C     7) Local variables.
C
C     a) Integer variables:
C
C     nFun     total number of function evaluations. 
C     nConstr  total number of constraint evaluations.
C     nGrad    total number of gradient evaluations.
C     nJacob   total number of Jacobian evaluations.
C     nHprod   total number of Hessian-vector products.
C     tSoc     total number of second order corrections used.
C     tbfgs    total number of BFGS iterations used.
C     tRest    total number of restorations.
C     tSteih   total number of Steihaug iterations.
C     tRej     total number of rejected steps.
C     nSoc     number of second order corrections at each iteration.
C     nbfgs    number of BFGS iterations at each main iteration.
C     nRest    number of restorations at each iteration.
C     nsteih   number of Steihaug iterations at each iteration. 
C     nRej     number of rejected steps at each iteration.
C     vflag    dcivert output flag.
C     itssmll  number of successive iterations with a too
C              small step.
C     DenseRow first row of the dense window of the Cholesky factor
C              of A*A'. See the AAtfact routine in the dcichol module.
C
C     b) Integer vectors:
C
C     Airow    row indices of the nonzero elements of A.
C     Apcol    pointer to the first nonzero element of each
C              row stored in A.
C     Mirow, Mpcol, Mcnel, perm, iperm
C              Data structure for the Cholesky decompostion of
C              A*A'. See the dcichol module for a description
C              of these parameters.
C
C     c) Real variables:
C
C     hnormc   two-norm of h(xc) (used only for printing information
C              about the current iteration).
C     hnorm    two-norm of h(x).
C     hnormi   inf-norm of h(x).
C     Ln       Lagrangian at (x, lambda) (after the horizontal step).
C     Lc       Lagrangian at (xc, lambda) (after the vertical step).
C     Lref     vaue of Lc at the last iteration in which 
C              DLv>(-0.5*DLh).
C     DLh      Ln - Lc (horizontal change of the lagrangian).
C     DLv      Lc - Ln (vertical change of the Lagrangian). 
C     gpnorm   ||gp||.
C     gpnormi  inf-norm of gp.
C     ngp      ||gp||/(||g||+1).
C     engp     an approximate value for ngp.
C     fx       value of f(x) after the horizontal step.
C     snorm    two-norm of the horizontal step.
C     xnorm    two-norm of the current point.
C     lnorm    two-norm of the vector of Lagrange multipliers.
C     rho      trust cylinder radius.
C     rhomax   Upper bound for rho.
C     hdelta   dcihoriz trust region radius.
C     vdelta   dcitrust trust region radius.
C     gphinorm ||At*h|| at the solution.
C
C     d) Real vectors:
C
C     g        gradient of f(x).
C     gp       projected gradient.
C     xc       approximate solution after a restoration.
C     tmpi     temporary vector (i=1,...,7).
C     A        nonzero elements of the Jacobian of the constraints.
C     Diag, Mval
C              Data structure for the Cholesky decompostion of
C              A*A'. See the dcichol module for a description
C              of these parameters.
C
C     e) Logical variable:
C
C     permavail parameter that indicates if the AMD algorithm has
C              already been used to compute the row and column
C              permutations of A*A'.


      IMPLICIT NONE

      EXTERNAL Fun, Constr, Grad, Jacob, Hprod

      INTEGER nmax, mmax, amax, lmax, bmax, cmax

C     Defining vector dimensions. 
C     nmax     maximum number of primal variables.
C     mmax     maximum number of constraints.
C     amax     maximum number of nonzero elements in the Jacobian
C              of the constraints.
C     lmax     maximum number of nonzero elements in the band
C              preconditioner used when applying the CG method to 
C              solve a system.
C     bmax     maximum number of BFGS hessian updates used to 
C              compute the vertical step.
C     cmax     maximum number of nonzero elements in the Cholesky
C              factor of A*A'.
C     Caution: nmax must be greater or equal to mmax.

      PARAMETER (nmax = 60000, mmax = 40000, amax = 1000000)
      PARAMETER (lmax = 210000, bmax = 200, cmax = 9000000)

C     Parameters.

      INTEGER m, n, iter, ierr
      REAL*8  fc
      REAL*8  x(n)
      REAL*8  h(m), lambda(m)
      LOGICAL lincon

C     Parameters read from disk.

      INTEGER maxrest, maxitst, minitst, maxitcg, minitcg, maxit
      INTEGER maxssmll, prnt, solprnt, tabprnt, dbgprnt, prcndiag
      INTEGER nfailv, bfgsupd, nstpmaxb, DenseWin
      REAL*8  csih, csig, rhomin, phi1, phi2, EpsChol, lbdmax
      REAL*8  kappa1, kappa2, kappa3, kappa4, DenseFrc
      REAL*8  epsst1, epsst2, epsst3, epscg1a, epscg1r
      REAL*8  epscg2a, epscg2r, epscg3a, epscg3r, epscg4a, epscg4r
      REAL*8  zeta1, zeta2, zeta3, alphar, alphai, alphab, alphas
      REAL*8  thetar, cls1, cls2, eta1, eta2, eta3, stepmaxb
      REAL*8  delta0, mindelta, infdelta, minstep, relitst, relitcg
      LOGICAL UseChol, UpdChol
      CHARACTER*200 solpath, tabpath

C     Local variables.

      INTEGER nFun, nConstr, nGrad, nJacob, nHprod, nSoc, nbfgs
      INTEGER nRest, nSteih, nRej, tRest, tSteih, tRej, tSoc, tbfgs
      INTEGER vflag, itssmll, DenseRow
      INTEGER Airow(amax)
      INTEGER Apcol(nmax)
      INTEGER Mirow(cmax)
      INTEGER Mpcol(mmax+1)
      INTEGER Mcnel(mmax), perm(mmax), iperm(mmax)
      REAL*4  eltime
      REAL*8  hnorm, hnormc, hnormi, Ln, Lc, Lref, DLh, DLv
      REAL*8  gpnorm, gpnormi, gphinorm, fx, snorm, xnorm, lnorm 
      REAL*8  rho, rhomax, ngp, engp, vdelta, hdelta
      REAL*8  g(nmax), gp(nmax), xc(nmax)
      REAL*8  tmp5(nmax), tmp6(nmax), tmp7(nmax), tmp8(nmax)
      REAL*8  tmp9(nmax), tmp10(nmax)
      REAL*8  tmp1(mmax), tmp2(mmax), tmp3(mmax), tmp4(mmax), tmp(mmax)
      REAL*8  Diag(mmax)
      REAL*8  tmp11(nmax*(bmax+1)), tmp12(nmax*(bmax+1))
      REAL*8  tmp13(nmax*(bmax+1))
      REAL*8  A(amax)
      REAL*8  Prec(lmax)
      REAL*8  Mval(cmax)
      LOGICAL permavail

C     LOCAL CONVERGENCE DEBUGGING
      REAL*8  ngpxk, ngpxck
      INTEGER naflag, locali
      LOGICAL LimLbd, forcegrad, gavail
C     It is possible to use vector Mval to store Prec.

C     Functions called by the routine.

      INTEGER idamax
      REAL*4  secnds
      REAL*8  CutestFun, ddot, dnrm2

C     OPEN(unit=555, file='dci.lixo', status='new')

C     LOCAL CONVERGENCE
      common / dciforcegrad / forcegrad

C     Setting up constants.

      CALL ReadParameters(maxrest, maxitst, minitst, maxitcg, minitcg, 
     &                    maxit, maxssmll, prnt, solprnt, tabprnt,
     &                    dbgprnt,
     &                    nfailv, bfgsupd, nstpmaxb, DenseWin, prcndiag, 
     &                    solpath, tabpath, relitst, relitcg, csih, 
     &                    csig, rhomin, phi1, phi2, kappa1, kappa2,  
     &                    kappa3, kappa4, thetar, epsst1, epsst2,  
     &                    epsst3, epscg1a, epscg1r, epscg2a, epscg2r,  
     &                    epscg3a, epscg3r, epscg4a, epscg4r, EpsChol,  
     &                    DenseFrc, lbdmax, zeta1, zeta2, zeta3, alphar,  
     &                    alphai, alphas, alphab, eta1, eta2, eta3,   
     &                    delta0, mindelta, infdelta, minstep, cls1,  
     &                    cls2, stepmaxb, UseChol, UpdChol)

C     Starting the stopwatch.

      eltime = secnds(0.0)
      
C     Setting initial values.
      
      CALL InitialValues(m, n, iter, itssmll, trest, tsteih, trej, tSoc,
     *                   tbfgs, nFun, nConstr, nGrad, nJacob, nHprod, 
     *                   nSoc, nbfgs, ierr, vflag, delta0, relitst, 
     *                   relitcg, maxitst, minitst, maxitcg, minitcg, 
     *                   csig, csih, epsst2, epscg1a, epscg1r, epscg2a,  
     *                   epscg3a, epscg4a, infdelta, minstep, fc, hnorm, 
     *                   hnormi, snorm, Ln, Lref, DLh, DLv, gpnorm, 
     *                   gpnormi, ngp, engp, rhomax, rho, vdelta, 
     *                   hdelta, x, h, g, gp, lambda, A, Airow, Apcol, 
     *                   Diag, Mval, Mirow, Mpcol, Mcnel, perm, iperm, 
     *                   permavail, EpsChol, DenseWin, DenseFrc, 
     *                   DenseRow, UseChol, UpdChol, lbdmax, prcndiag, 
     *                   phi1, Prec, tmp1, tmp2, tmp3, tmp4, tmp5)
C    *                   Fun, Grad, Constr, Jacob)

C     Main loop.

C     LOCAL CONVERGENCE
C     DO locali = 1,n
C       x(locali) = 0.0D0
C     ENDDO
      forcegrad = .TRUE.
      ngpxk = 0
      ngpxck = 0

      DO WHILE (((hnormi.GT.csih).OR.
     *          ((gpnormi.GT.csig).AND.(ngp.GT.(csig*1.0D-2))))
     *          .AND.(iter.LE.maxit).AND.(tRest.LE.maxrest)
     *          .AND.(itssmll.LE.maxssmll).AND.(vflag.EQ.0)
     *          .AND.(rhomax.GE.rhomin))

        iter = iter + 1

C       Increasing vdelta and hdelta if they are smaller than mindelta.

        hdelta = DMAX1(hdelta, mindelta)
        vdelta = DMAX1(vdelta, mindelta)
C       IF (vdelta.GT.(50.0D0*hdelta)) vdelta = 50.0D0*hdelta
C       IF (hdelta.GT.(50.0D0*vdelta)) hdelta = 50.0D0*vdelta
C       vdelta = 10.0D0*vdelta

C       Computing the vertical step.

        CALL dcivert(m, n, x, h, A, Airow, Apcol, Diag, Mval, Mirow,
     $               Mpcol, Mcnel, perm, iperm, permavail, DenseRow, 
     $               lincon, csih, csig, phi1, phi2, nfailv, rhomax, 
     $               maxrest, hnorm, engp, vdelta, infdelta, kappa1, 
     $               kappa2, kappa3, kappa4, thetar, bfgsupd, cls1, 
     $               cls2, alphab, stepmaxb, nstpmaxb, epscg1a, epscg1r, 
     $               epscg3a, epscg3r, maxitcg, prcndiag, Prec, UseChol, 
     $               EpsChol, DenseWin, DenseFrc, UpdChol, lbdmax, iter, 
     $               tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, 
     $               tmp9, tmp10, tmp11, tmp12, tmp13, tmp, xc, g, 
     $               lambda, gp, gpnorm, ngp, lnorm, rho, nRest, nbfgs, 
     $               nConstr, nGrad, nJacob, vflag)
C    $               , Constr, Grad, Jacob)

        hnormc = hnorm
        tRest = tRest + nRest
        tbfgs = tbfgs + nbfgs

C       LOCAL CONVERGENCE
        ngpxck = gpnorm

C       Computing f at xc.

        fc = CutestFun(n, xc)

C       Checking for convergence.

        hnormi = DABS(h(idamax(m, h, 1)))
        gpnormi = DABS(gp(idamax(n, gp, 1)))
        IF (((hnormi.GT.csih).OR.
     $      ((gpnormi.GT.csig).AND.(ngp.GT.(csig*1.0D-2))))
     $      .AND.(vflag.EQ.0).AND.(rhomax.GE.rhomin)) THEN

C         Computing the new Lagrangian and its Hessian.

          Lc = fc + ddot(m, lambda, 1, h, 1)
          DLv = Lc - Ln

C         Updating rhomax and Lref.

          IF (DLv.GE.(0.5D0*(Lref - Ln))) THEN
            rhomax = 0.5D0*rhomax
          ENDIF
          IF (DLv.GT.(-0.5D0*DLh)) THEN
            Lref = Lc
          ENDIF

C         Computing the horizontal step.

          CALL dcihoriz(m, n, xc, lambda, h, g, gp, A, Airow, Apcol, 
     *                  UseChol, Diag, Mval, Mirow, Mpcol, Mcnel, perm,
     *                  iperm, DenseRow, Lc, hnorm, hdelta, csih, 
     *                  epsst1, epsst2, epsst3, maxitst, epscg2a, 
     *                  epscg2r, epscg4a, epscg4r, maxitcg, prcndiag, 
     *                  rho, zeta1, zeta2, zeta3, alphar, alphai, 
     *                  alphas, eta1, eta2, eta3, x, fx, h, Ln, snorm, 
     *                  DLh, Prec, tmp, tmp1, tmp2, tmp3, tmp4, tmp5, 
     *                  tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, nSteih, 
     *                  nrej, nFun, nConstr, nHprod, nSoc)
C    *                  , Fun, Constr, Hprod)

          tSoc = tSoc + nSoc
          tSteih = tSteih + nSteih
          tRej = tRej + nRej

C         Estimating ngp at the new point. 

          IF (nsteih.GT.0) THEN
            engp = DABS(DLh)/(DABS(fx - fc) + snorm)
C b2        engp = ngp 
C b3        engp = SQRT(DABS(DLh))/(SQRT(DABS(fx - fc)) + snorm) !b3
C b4        engp = SQRT(DABS(DLh))/(SQRT(DABS(fx - fc)) + 1.0D0) !b4
C b5        engp = SQRT(DABS(DLh))/(SQRT(DABS(fx - fc)) + 1.0D-8) !b5
C b6        engp = DMIN1(SQRT(DABS(DLh))/(SQRT(DABS(fx-fc))+1.0D0),ngp)     
C b7        engp = DMAX1(SQRT(DABS(DLh))/(SQRT(DABS(fx-fc))+1.0D0),ngp)     
C b8        engp = (SQRT(DABS(DLh))/(SQRT(DABS(fx-fc))+1.0D0)+ngp)/2.0D0
C b9        engp = SQRT(DABS(DLh))/(SQRT(DABS(fx - fc)) + 2.0D0) !b9
          ELSE
            engp = ngp
          ENDIF

          xnorm = dnrm2(n, x, 1)
          IF (snorm.LT.(xnorm*minstep)) THEN
            itssmll = itssmll + 1
          ELSE
            itssmll = 0
          ENDIF

        ELSE

          CALL dcopy(n, xc, 1, x, 1)

        ENDIF

C       LOCAL CONVERGENCE 
C       CALCULAR gp
        LimLbd = .TRUE.
        gavail = .FALSE.
        CALL CutestGrad(n, x, g)
        CALL NAprojCh(m, n, A, Airow, Apcol, Diag, Mval, Mirow,
     $                mpcol, Mcnel, perm, iperm, DenseRow, lambda,
     $                LimLbd, lbdmax, g, gp, naflag)
        ngpxk = dnrm2(n, gp, 1)
C       IMPRIMIR
        if (dbgprnt.gt.0) THEN
          PRINT *,"|gp(xk)|/|gp(xck)| = ", ngpxk/ngpxck
        endif

C       Printing some information about the current iteration.          
          
        IF (prnt.GT.0) THEN
          IF (MOD(iter,prnt).EQ.0) THEN      
            CALL KPrint(iter, nSteih, nRest, nRej, nSoc, nbfgs, fc, 
     *                  hnormc, hnorm, gpnorm, lnorm, rho, rhomax, 
     *                  snorm, vdelta, hdelta, Ln, DLh, DLv)
          ENDIF                          
        ENDIF

      ENDDO

C     Defining the stopping criterion.

      IF (vflag.GT.0) THEN
        ierr = 4
      ELSE IF (rhomax.LT.rhomin) THEN
        ierr = 1
      ELSE IF (iter.GT.maxit) THEN
        ierr = 2
      ELSE IF (tRest.GT.maxrest) THEN
        ierr = 3
      ELSE IF (itssmll.GT.maxssmll) THEN
        ierr = 5
      ELSE
        ierr = 0
      ENDIF

C     Stoppinng the stopwatch.

      eltime = secnds(eltime)

C     Printing some data before returning.
      
      IF (prnt.GE.0) THEN
        CALL MtVprod(n, A, h, tmp5, Airow, Apcol)
        gphinorm = dnrm2(n, tmp5, 1)
        CALL LastPrint(ierr, iter, tSteih, tRest, tRej, nFun, nConstr,
     *                 nGrad, nJacob, tSoc, tbfgs, fc, hnormi, gpnormi,  
     *                 gphinorm, rho, rhomax, snorm, vdelta, hdelta, 
     *                 Ln, DLh, DLv, eltime)
      
      ENDIF                                      

      IF (solprnt.GE.0) THEN
        CALL SaveSolution(solprnt, solpath, m, n, x, lambda, fc, h,
     *                    hnormi, iter, ierr, eltime, tRej, tRest, nFun,
     *                    nConstr, nGrad, nJacob, nHprod, tSoc, tbfgs, 
     *                    tSteih, maxrest, maxitst, maxitcg, maxit, 
     *                    maxssmll, prcndiag, csih, csig, phi1, phi2, 
     *                    epsst1, epsst2, epsst3, epscg1a, epscg1r, 
     *                    epscg2a, epscg2r, epscg3a, epscg3r, epscg4a, 
     *                    epscg4r, eta1, eta2, kappa1, kappa2, kappa3, 
     *                    kappa4, zeta1, zeta2, zeta3, alphar, alphai, 
     *                    mindelta, infdelta, minstep)
      ENDIF                                      
      
      IF (tabprnt.GE.0) THEN
        CALL SaveTable(tabpath, m, ierr, iter, tRej, tRest, nFun, 
     *                 nJacob, nHprod, tSoc, tbfgs, tSteih, eltime, fc, 
     *                 hnormi)
      ENDIF                                      

C     CLOSE(555)
      
      RETURN
      END
