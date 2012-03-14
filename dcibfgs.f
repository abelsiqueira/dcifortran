!     ------------------------------------------------------------------
!     PROGRAM : DCI.
!     MODULE  : DCIBFGS.
!     VERSION : 2.0.
!     DATE    : March, 2006.
!     CONTENTS: subroutine dcibfgs
!               subroutine Hiprod (internal)
!               subroutine linesearch (internal)
!               subroutine zoom (internal)
!               function interpolate (internal)
!
!     External routines called by this module:
!
!     From dcicute : Constr, Jacob. 
!     From dciblas : MtVprod.
!     From the blas: daxpy, dcopy, ddot, dnrm2, dscal.
!     ------------------------------------------------------------------

      SUBROUTINE dcibfgs(m, n, x, h, A, Airow, Apcol, hnorm, davail,
     &                   dn, dnorm, rho, maxit, c1, c2, alphab, stepmax, 
     &                   nstpmax, g, s, xold, gold, sold, alpha, csi,
     &                   d, u, y, iter, nConstr, nJacob, iout)
C    &                   , Constr, Jacob)

C     This routine reduces the infeasibility of the DCI problem by 
C     applying the BFGS algorithm (with the line search proposed by 
C     Moré e Thuente) to the nonlinear least squares problem 
C     min  0.5||h(x)||^2.
C
C     Input parameters:
C
C     m, n     number of constraints and variables.
C     x        current approximation for the solution.
C     h        constraint vector, h(x).
C     A        Jacobian of the constraints.
C     hnorm    two-norm of h(x).
C     davail   parameter that indicates if an initial direction is 
C              available. If davail = false, we use d(1) = s(1) = -g.
C     dn       Initial (Gauss-Newton) direction vector.
C     dnorm    two-norm of the initial ditection vector (if it is 
C              available).
C     rho      trust cylinder radius.
C     maxit    maximum number of iterations allowed.
C     c1       parameter used in the sufficient decrease condition
C              (suggestion: 1.0e-4). 
C     c2       parameter used in the curvature condition 
C              (suggestion: 0.1).
C     alphab   parameter used to increase the step in the extrapolation 
C              phase of the algorithm (suggestion: 2.0).
C     stepmax  largest step allowed (suggestion: 5.0).
C     nstpmax  maximum number of step length reductions allowed per 
C              call of the steplength function (suggestion: 10). 
C     constr   A function that computes the constraint vector.
C     jacob    A function that computes the Jacobian of the constraints.
C     nConstr  total number of constraint evaluations.
C     nJacob   total number of Jacobian evaluations.
C
C     Output parameters:
C
C     x        new point.
C     h        constraint vector, h(x).
C     hnorm    two-norm of h(x).
C     A        Jacobian of the constraints.
C     iter     number of BFGS steps.
C     nConstr  see description above.
C     nJacob   see description above.
C     iout     stopping flag:
C              if iout = 0, hnorm <= rho.
C              if iout = 1, the zoom function has failed.
C              if iout = 2, iter = maxit.
C              if iout = 3, norm(g) is too small.
C              if iout = 4, the step become too small.
C
C     Temporary vectors passed as parameters.
C
C     g(n), s(n), xold(n), gold(n), sold(n), alpha(maxit), csi(maxit),
C     d(n*(maxit+1)), u(n*maxit), y(n*(maxit+1)).

      IMPLICIT NONE

      EXTERNAL Constr, Jacob

C     Subroutine parameters.

      INTEGER m, n, maxit, nstpmax, iter, nConstr, nJacob, iout
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      REAL*8  hnorm, dnorm, rho, c1, c2, alphab, stepmax
      REAL*8  x(n), dn(n), g(n), s(n), xold(n), gold(n), sold(n)
      REAL*8  alpha(maxit), csi(maxit)
      REAL*8  d(n*(maxit+1)), u(n*maxit), y(n*(maxit+1))
      REAL*8  h(m)
      REAL*8  A(*)
      LOGICAL davail

C     Internal variables:

      INTEGER i
      REAL*8  angle, f, gnorm, gtd, lambda, snorm, xnorm

C     Functions called by the routine.

      REAL*8  ddot, dnrm2

C     Setting up some variables.

      f = 0.5D0*hnorm**2 
      CALL MtVprod(n, A, h, g, Airow, Apcol) 
      gnorm = dnrm2(n, g, 1)
      iter = 0
      iout = 0

C     Returning from here if the initial point is feasible.

      IF (hnorm.LT.rho) RETURN

C     Returning from here if the gradient is too small.

      if (gnorm.LT.1.0D-10) THEN
        iout = 3
        RETURN
      ENDIF

C     Creating the data structure for the BFGS algorithm.
C     s = H*g, but we use the identity matrix as the initial H.

      CALL dcopy(n, g, 1, s, 1)
      CALL dscal(n, -1.0D0, s, 1)
      snorm = gnorm
      IF (davail) THEN
        angle = ddot(n, s, 1, dn, 1)/(dnorm*snorm)
        IF (angle.LT.1.0D-4) THEN
            CALL dcopy(n, s, 1, d, 1)
        ELSE
            CALL dcopy(n, dn, 1, d, 1)
        ENDIF
      ELSE
        CALL dcopy(n, s, 1, d, 1) 
      ENDIF

C     Performing the first line search.

      gtd = ddot(n, g, 1, d, 1)
      CALL dcopy(n, x, 1, xold, 1)
      CALL dcopy(n, g, 1, gold, 1)

      CALL linesearch(m, n, xold, d, f, gtd, c1, c2, alphab, 
     &                stepmax, nstpmax, x, h, A, Airow, Apcol, 
     &                g, lambda, nConstr, nJacob, iout)

C     y(:,iter+1) = g - gold
      CALL dcopy(n, g, 1, y(iter*n+1), 1)
      CALL daxpy(n, -1.0D0, gold, 1, y(iter*n+1), 1)
C     d(:,iter+1) = x - xold
      CALL dcopy(n, x, 1, d(iter*n+1), 1)
      CALL daxpy(n, -1.0D0, xold, 1, d(iter*n+1), 1)
C     dnorm = norm(d(:,iter+1))
      dnorm = dnrm2(n, d(iter*n+1), 1)
      hnorm = dnrm2(m, h, 1)
      xnorm = dnrm2(n, x, 1)

C     Main loop.

      DO WHILE ((hnorm.GT.rho).AND.(iter.LT.maxit).AND.(iout.EQ.0).AND. 
     &          (dnorm.GE.(xnorm*1.0D-17)))

        iter = iter + 1

C       Calculando alpha(iter).
  
C       alpha(iter) = d(:,iter)'*y(:,iter)
        alpha(iter) = ddot(n, d((iter-1)*n+1), 1, y((iter-1)*n+1), 1)
        IF (alpha(iter).LT.1.0D-12) alpha(iter) = 1.0D-12

C       Calculando u(iter) e csi(iter).
  
        CALL dcopy(n, s, 1, sold, 1)
C       s = H*g, but we use the identity matrix as the initial H.
        CALL dcopy(n, g, 1, s, 1)
        DO i = 1,iter-1
          CALL Hiprod(i, n, alpha, csi, d, u, y, g, s)
        ENDDO
C       u(:,iter) = s + sold
        CALL dcopy(n, s, 1, u((iter-1)*n+1), 1)
        CALL daxpy(n, 1.0D0, sold, 1, u((iter-1)*n+1), 1)
C       csi(iter) = 1+y(:,iter)'*u(:,iter)/alpha(iter)
        csi(iter) = 1.0D0 + ddot(n,y((iter-1)*n+1),1,u((iter-1)*n+1),1)/
     &              alpha(iter)
  
C       Calculando s = -H.g e d(iter+1).
  
        CALL Hiprod(iter, n, alpha, csi, d, u, y, g, s)
        CALL dscal(n, -1.0D0, s, 1)
  
C       Calculando o novo ponto e os novos valores das restrições.
  
        gtd = ddot(n, g, 1, s, 1)
        CALL dcopy(n, x, 1, xold, 1)
        CALL dcopy(n, g, 1, gold, 1)
        CALL linesearch(m, n, xold, s, f, gtd, c1, c2, alphab, 
     &                  stepmax, nstpmax, x, h, A, Airow, Apcol, 
     &                  g, lambda, nConstr, nJacob, iout)

C       Calculando y(iter+1) e d(iter+1).
  
C       y(:,iter+1) = g - gold
        CALL dcopy(n, g, 1, y(iter*n+1), 1)
        CALL daxpy(n, -1.0D0, gold, 1, y(iter*n+1), 1)
C       d(:,iter+1) = x - xold
        CALL dcopy(n, x, 1, d(iter*n+1), 1)
        CALL daxpy(n, -1.0D0, xold, 1, d(iter*n+1), 1)
C       dnorm = norm(d(:,iter+1))
        dnorm = dnrm2(n, d(iter*n+1), 1)
        hnorm = dnrm2(m, h, 1)
        xnorm = dnrm2(n, x, 1)
  
      ENDDO

      iter = iter + 1
      IF (hnorm.LE.rho) THEN
        iout = 0
      ELSE IF (iter.GT.maxit) THEN
        iout = 2
      ELSE IF (dnorm.LT.(xnorm*1.0D-10)) THEN
        iout = 4
      ENDIF
  
      RETURN
      END ! dcibfgs.

C     ------------------------------------------------------------------
      
      SUBROUTINE Hiprod(i, n, alpha, csi, d, u, y, v, Hv)

      IMPLICIT NONE

C     Subroutine parameters.

      INTEGER i, n 
      REAL*8  alpha(*), csi(*)
      REAL*8  d(n), u(n), y(n), v(n), Hv(n)

C     Internal variables

      REAL*8  beta, rho

C     Function called by the routine.

      REAL*8  ddot

      beta = ddot(n, d((i-1)*n+1), 1, v, 1)/alpha(i)
      rho = ddot(n, y((i-1)*n+1), 1, Hv, 1)/alpha(i)
      CALL daxpy(n, (csi(i)*beta-rho), d((i-1)*n+1), 1, Hv, 1)
      CALL daxpy(n, -beta, u((i-1)*n+1), 1, Hv, 1)

      RETURN
      END ! Hiprod.

C     ------------------------------------------------------------------
      
      SUBROUTINE linesearch(m, n, x0, d, f, gtd0, c1, c2, alpha, 
     &                      lbdmax, itmax, x, h, A, Airow, Apcol, 
     &                      g, lambda, nConstr, nJacob, iout)
C    &                      , constr, jacob)
C     This function computes a new point x = x0 + lambda*d that satisfies 
C     the strong Wolfe conditions. f0 = phi(0) = f(x0), gtd0 = phi'(0) = 
C     g(x0)'*d = h(x0)'*A(x0)*d, c1 is a parameter used in the sufficient 
C     decrease condition (suggestion: 1e-4). c2 is a parameter used in
C     the curvature condition (suggestion: 0.9). alpha is a parameter 
C     used to increase lambda (suggestion: 2). lbdmax is the upper limit
C     for lambda (suggestion: 5). itmax is the maximum number of step
C     lengths to be tried in the zoom function (suggestion: 10).

      IMPLICIT NONE

      EXTERNAL constr, jacob

C     Subroutine parameters.

      INTEGER m, n, itmax, nConstr, nJacob, iout
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      REAL*8  f, gtd0, c1, c2, alpha, lbdmax, lambda
      REAL*8  x0(n), d(n), x(n), g(n)
      REAL*8  h(m)
      REAL*8  A(*)

C     Internal variables

      REAL*8  fold, lbdold, gtd, f0
      LOGICAL first, ishi

C     Function called by the routine.

      REAL*8  ddot

C     Setting up some variables.

      f0     = f
      fold   = f
      lbdold = 0.0D0
      lambda = 1.0D0
      lbdmax = DMAX1(1.0D0, lbdmax)
      first  = .TRUE.
      iout   = 0

C     Main loop.
  
      DO WHILE (.TRUE.)
  
C       Computing phi(lambda).
  
        CALL dcopy(n, x0, 1, x, 1)
        CALL daxpy(n, lambda, d, 1, x, 1)
        CALL CuterConstr(m, n, x, h)
        nConstr = nConstr + 1
        f = 0.5D0*ddot(m, h, 1, h, 1)
  
C       Zooming if the step is too large and need to be reduced.
    
        IF ((f.GT.(f0 + c1*lambda*gtd0)).OR.
     &      ((f.GE.fold).AND.(.NOT.first))) THEN

          ishi = .TRUE.
          CALL zoom(m, n, lbdold, lambda, c1, c2, x0, d, f0, gtd0,  
     &              fold, f, ishi, itmax, lambda, x, h, f, A, 
     &              Airow, Apcol, g, nConstr, nJacob, iout)
          EXIT
        ENDIF

C       Computing phi'(lambda).
    
        CALL CuterJacob(m, n, x, A, Airow, Apcol)
        nJacob = nJacob + 1
        CALL MtVprod(n, A, h, g, Airow, Apcol) 
        gtd = ddot(n, g, 1, d, 1)
    
C       Verifying if the curvature condition is satisfied.
    
        IF (DABS(gtd).LE.(-c2*gtd0)) THEN
          EXIT
        ENDIF
      
C       Zooming if the slope of phi is positive at x(lambda).
    
        IF (gtd.GE.0.0D0) THEN
          ishi = .FALSE.
          CALL zoom(m, n, lambda, lbdold, c1, c2, x0, d, f0, gtd0,  
     &              f, fold, ishi, itmax, lambda, x, h, f, A, 
     &              Airow, Apcol, g, nConstr, nJacob, iout)
          EXIT
        ENDIF
    
        IF (lambda.EQ.lbdmax) THEN
          EXIT
        ENDIF

        lbdold = lambda
        fold = f
        lambda = DMIN1(alpha*lambda, lbdmax)
        first = .FALSE.
    
      ENDDO

      RETURN
      END ! linesearch.
  
C     ------------------------------------------------------------------
      
      SUBROUTINE zoom(m, n, lbdlo0, lbdhi0, c1, c2, x0, d, f0, gtd0, 
     &                flo0, fhi0, ishi0, itmax, lambda, x, h, f, A, 
     &                Airow, Apcol, g, nConstr, nJacob, iout)
C    &                , constr, jacob)  

C     This function computes a steplength lambda, between lbdlo0 and 
C     lbdhi0, such that x0 + lambda*d satisfies the strong Wolfe 
C     conditions. f0 = phi(0) = f(x0), gtd0 = phi'(0) = g(x0)'*d = 
C     h(x0)'*A(x0)*d, flo0 = phi(lbdlo) and fhi0 = phi(lbdhi). ishi0
C     is a parameter that must be set to 'true' if lbdhi0 is the last 
C     computed value of lambda. itmax is the maximum number of step 
C     lengths to be tried. c1 ans c2 are input parameters for the 
C     linesearch routine. Notice that lbdlo0 and lbdhi0 need to to be 
C     ordered.

      IMPLICIT NONE

      EXTERNAL constr, jacob

C     Subroutine parameters.

      INTEGER m, n, itmax, nConstr, nJacob, iout
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      REAL*8  lbdlo0, lbdhi0, c1, c2, f0, gtd0, flo0, fhi0, lambda, f
      REAL*8  x0(n), d(n), x(n), g(n)
      REAL*8  h(m)
      REAL*8  A(*)
      LOGICAL ishi0

C     Internal variables

      INTEGER iter
      REAL*8  gtd, diff, lbdlo, lbdhi, flo, fhi
      LOGICAL first, enough, ishi

C     Function called by the routine.

      REAL*8  ddot, interpolate

      iter   = 1
      enough = .FALSE.
      flo    = flo0
      fhi    = fhi0
      lbdlo  = lbdlo0
      lbdhi  = lbdhi0
      ishi   = ishi0
  
      DO WHILE ((iter.LE.itmax).AND.(.NOT.enough))
  
C       Interpolating to find a new lambda.
    
        first = ((iter.EQ.1).OR.(lbdlo.EQ.0.0D0).OR.(lbdhi.EQ.0.0D0))
        IF (ishi) THEN
          lambda = interpolate(f0, gtd0, fhi, flo, lbdhi, lbdlo, first)
        ELSE
          lambda = interpolate(f0, gtd0, flo, fhi, lbdlo, lbdhi, first)
        ENDIF
    
C       Keeping lambda between lbdlo and lbdhi.
    
        diff = DABS(lbdlo - lbdhi)
        lambda = DMAX1(DMIN1(lbdlo, lbdhi) + 0.01D0*diff,
     &                 DMIN1(lambda, DMAX1(lbdlo, lbdhi) - 0.01D0*diff))
    
C       Evaluating phi(lambda).
    
        CALL dcopy(n, x0, 1, x, 1)
        CALL daxpy(n, lambda, d, 1, x, 1)
        CALL CuterConstr(m, n, x, h)
        nConstr = nConstr + 1
        f = 0.5D0*ddot(m, h, 1, h, 1)
  
        IF ((f.GT.(f0 + c1*lambda*gtd0)).OR.(f.GE.flo)) THEN

C         The sufficient decrease condition was not satisfied.
    
          lbdhi = lambda
          fhi = f
          ishi = .TRUE.
          iter = iter + 1
      
        ELSE
      
C         Computing phi'(lambda).
    
          CALL CuterJacob(m, n, x, A, Airow, Apcol)
          nJacob = nJacob + 1
          CALL MtVprod(n, A, h, g, Airow, Apcol) 
          gtd = ddot(n, g, 1, d, 1)

          IF (DABS(gtd).LE.(-c2*gtd0)) THEN

C           Exiting if the curvature condition is satisfied.
      
            enough = .TRUE.

          ELSE

C           Keeping phi'(lbdlo)*(lbdhi - lbdlo) < 0.
      
            IF (gtd*(lbdhi-lbdlo).GE.0.0D0) THEN
              lbdhi = lbdlo
              fhi = flo
            ENDIF
      
            lbdlo = lambda
            flo = f
            ishi = .FALSE.
            iter = iter + 1

          ENDIF
      
        ENDIF
    
      ENDDO
  
      IF (iter.GT.itmax) THEN

C       The sufficient decrease condition was not satisfied.
C       Recomputing phi'.
    
        IF (ishi) THEN
          CALL CuterJacob(m, n, x, A, Airow, Apcol)
          nJacob = nJacob + 1
          CALL MtVprod(n, A, h, g, Airow, Apcol) 
          gtd = ddot(n, g, 1, d, 1)
        ENDIF
    
C       The algorithm has failed.
    
        iout = 1
    
      ELSE
    
C       The Wolfe conditions were satisfied.
    
        iout = 0
    
      ENDIF

      RETURN
      END ! zoom

C     ------------------------------------------------------------------

      REAL*8 FUNCTION interpolate(f0, fl0, fm1, fm2, lm1, lm2, first)

C     This function performs a (quadratic or cubic) interpolation
C     to obtain lambda, the steplength used in the zoom procedure.
C     f0 = phi(0), fl0 = phi'(0), lm1 = lambda(k-1), lm2 = lambda(k-2),
C     fm1 = phi(lm1) and fm2 = phi(lm2). first is a parameter that 
C     must be set to 'true' only if lm2 and fm2 are not available.

      IMPLICIT NONE

C     Subroutine parameters.

      REAL*8  f0,  fl0, fm1, fm2, lm1, lm2
      LOGICAL first

C     Internal variables.

      REAL*8  c, d, delta, lambda

      IF (first) THEN

C       Using a quadratic function to obtain lambda.
    
        lambda = -0.5D0*fl0*lm1**2/(fm1-f0-fl0*lm1)
  
      ELSE

C       Using a cubic function to obtain lambda.
    
        c = (-(lm2/lm1**2)*(fm1-f0-fl0*lm1)+
     &       (lm1/lm2**2)*(fm2-f0-fl0*lm2))/(lm1-lm2)
        d = ((fm1-f0-fl0*lm1)/lm1**2-(fm2-f0-fl0*lm2)/lm2**2)/(lm1-lm2)
        delta = c**2 - 3.0D0*fl0*d
        IF (delta.GE.0.0D0) THEN
          lambda = (-c + DSQRT(delta))/(3.0D0*d)
        ELSE
          lambda = -0.5D0*fl0*lm1**2/(fm1-f0-fl0*lm1)
        ENDIF    
  
      ENDIF
  
C     Discarding lambda if it is too close or too far from the last value.
  
      IF ((DABS(lambda-lm1).LT.(1.0D-3*lm1)).OR.
     &    (lambda.LT.(1.0D-3*lm1))) THEN
        lambda = 0.5D0*lm1
      ENDIF

      interpolate = lambda

      RETURN
      END ! interpolate
