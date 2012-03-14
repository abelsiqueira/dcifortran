!     ------------------------------------------------------------------
!     PROGRAM : DCI.
!     MODULE  : DCIIO.
!     VERSION : 2.0.
!     DATE    : April/1998.
!     CONTENTS: subroutine ReadParameters.
!               subroutine InitialValues.
!               subroutine KPrint.
!               subroutine LastPrint.
!               subroutine SaveSolution.
!               subroutine SaveTable.
!
!     External routines called by this module:
!
!     From dcicute : Fun, Grad, Constr, Jacob. 
!     From dcichol : choldec.
!     From dciblas : vexpn, MtVprod.
!     From the blas: dcopy, ddot, daxpy, dnrm2, idamax
!     ------------------------------------------------------------------

      SUBROUTINE ReadName(un, name, val, stat) 

      IMPLICIT NONE

      INTEGER        un, stat, pos, l
      CHARACTER*10   name
      CHARACTER*200  val, line
      
      READ(unit=un, fmt=555, iostat=stat) line
 555  FORMAT(A200)

      IF (stat.EQ.0) THEN
        l = LEN(line)
        pos = 1
        DO WHILE (line(pos:pos).EQ.' ')
          pos = pos+1
        ENDDO
        line = line(pos:l)
        l = l-pos+1
        pos = INDEX(line,' ')
        name = line(1:pos-1)
        val = line((pos+1):l)
      ENDIF

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE ReadPath(val, path) 

      IMPLICIT NONE

      INTEGER        pos, l
      CHARACTER*200  val, path
      
      l = LEN(val)
      pos = 1
      DO WHILE (val(pos:pos).EQ.' ')
        pos = pos+1
      ENDDO
      val = val(pos:l)
      l = l-pos+1
      pos = INDEX(val,' ')
      path = val(1:pos-1)

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE ReadParameters(maxrest, maxitst, minitst, maxitcg, 
     &                          minitcg, maxit, maxssmll, prnt, solprnt, 
     &                          tabprnt, nfailv, bfgsupd, nstpmaxb, 
     &                          DenseWin, prcndiag, solpath, tabpath, 
     &                          relitst, relitcg, csih, csig, rhomin,  
     &                          phi1, phi2, kappa1, kappa2, kappa3,  
     &                          kappa4, thetar, epsst1, epsst2, epsst3,  
     &                          epscg1a, epscg1r, epscg2a, epscg2r,  
     &                          epscg3a, epscg3r, epscg4a, epscg4r,  
     &                          EpsChol, DenseFrc, lbdmax, zeta1, zeta2,  
     &                          zeta3, alphar, alphai, alphas, alphab,   
     &                          eta1, eta2, eta3, delta0, mindelta, 
     &                          infdelta, minstep, cls1, cls2, stepmaxb,
     &                          UseChol, UpdChol)

      IMPLICIT NONE

      INTEGER       maxrest, maxitst, minitst, maxitcg, minitcg, maxit
      INTEGER       maxssmll, prnt, solprnt, tabprnt, PrcnDiag
      INTEGER       stat, nfailv, bfgsupd, nstpmaxb, DenseWin
      REAL*8        csih, csig, rhomin, phi1, phi2, delta0, stepmaxb
      REAL*8        kappa1, kappa2, kappa3, kappa4
      REAL*8        epsst1, epsst2, epsst3
      REAL*8        epscg1a, epscg1r, epscg2a, epscg2r, epscg3a
      REAL*8        epscg3r, epscg4a, epscg4r, zeta1, zeta2, zeta3 
      REAL*8        alphar, alphai, alphas, alphab, eta1, eta2, eta3
      REAL*8        mindelta, infdelta, minstep, relitst, relitcg
      REAL*8        thetar, EpsChol, DenseFrc, lbdmax, cls1, cls2 
      LOGICAL       UseChol, UpdChol
      CHARACTER*200 solpath, tabpath
      CHARACTER*10  name
      CHARACTER*200 val

!     Setting initial values just in case the dci.spc file 
!     is not available.

      prnt     = 1
      solpath  = '.'
      tabpath  = '.'
      solprnt  = 0
      tabprnt  = 1
      maxrest  = 20000
      maxitst  = 100000
      minitst  = 100
      maxitcg  = 100000
      minitcg  = 100
      maxit    = 2000
      maxssmll = 5
      PrcnDiag = 3
      nfailv   = 3
      bfgsupd  = 10
      nstpmaxb = 50
      DenseWin = 200
      relitst  = 1.0D+1
      relitcg  = 1.0D+1
      csih     = 1.0D-12
      csig     = 1.0D-12
      rhomin   = 1.0D-8
      phi1     = 1.0D+0
      phi2     = 9.5D-1
      kappa1   = 1.0D-4
      kappa2   = 2.5D-1
      kappa3   = 7.0D-1
      kappa4   = 2.5D+0
      epsst1   = 1.0D-2
      epsst2   = 1.0D-14
      epsst3   = 1.0D-8
      epscg1a  = 1.0D-2
      epscg1r  = 1.0D-5
      epscg2a  = 1.0D-5
      epscg2r  = 1.0D-8
      epscg3a  = 1.0D-2
      epscg3r  = 1.0D-3
      epscg4a  = 1.0D-2
      epscg4r  = 1.0D-5
      zeta1    = 2.0D+0
      zeta2    = 1.0D+0
      zeta3    = 5.0D+0
      alphar   = 2.5D-1
      alphai   = 2.5D+0
      alphas   = 6.25D-2
      alphab   = 2.0D+0
      eta1     = 1.0D-4
      eta2     = 7.0D-1
      eta3     = 1.0D-1
      mindelta = 1.0D-4
      infdelta = 1.0D-10
      minstep  = 1.0D-3
      delta0   = 1.0D+5
      thetar   = 9.0D-1
      EpsChol  = 1.0D-16
      DenseFrc = 7.0D-1
      lbdmax   = 1.0D+6
      cls1     = 1.0D-4
      cls2     = 1.0D-1
      stepmaxb = 5.0D+0
      UseChol  = .TRUE.
      UpdChol  = .FALSE.


!     Reading the parameters.

      OPEN(unit=40, file='dci.spc', status='old', iostat=stat)

      IF (stat.EQ.0) THEN

        CALL ReadName(40, name, val, stat) 

        DO WHILE (stat.EQ.0)
          
          IF      (name.EQ.'MAXREST   ') THEN
            READ(val, fmt=*, iostat=stat) maxrest
          ELSE IF (name.EQ.'MAXITST   ') THEN
            READ(val, fmt=*, iostat=stat) maxitst
          ELSE IF (name.EQ.'MINITST   ') THEN
            READ(val, fmt=*, iostat=stat) minitst
          ELSE IF (name.EQ.'MAXITCG   ') THEN
            READ(val, fmt=*, iostat=stat) maxitcg
          ELSE IF (name.EQ.'MINITCG   ') THEN
            READ(val, fmt=*, iostat=stat) minitcg
          ELSE IF (name.EQ.'MAXIT     ') THEN
            READ(val, fmt=*, iostat=stat) maxit
          ELSE IF (name.EQ.'MAXSSMLL  ') THEN
            READ(val, fmt=*, iostat=stat) maxssmll
          ELSE IF (name.EQ.'PRCNDIAG  ') THEN
            READ(val, fmt=*, iostat=stat) prcndiag
          ELSE IF (name.EQ.'PRNT      ') THEN
            READ(val, fmt=*, iostat=stat) prnt
          ELSE IF (name.EQ.'SOLPRNT   ') THEN
            READ(val, fmt=*, iostat=stat) solprnt
          ELSE IF (name.EQ.'TABPRNT   ') THEN
            READ(val, fmt=*, iostat=stat) tabprnt
          ELSE IF (name.EQ.'NFAILV    ') THEN
            READ(val, fmt=*, iostat=stat) nfailv
          ELSE IF (name.EQ.'BFGSUPD   ') THEN
            READ(val, fmt=*, iostat=stat) bfgsupd
          ELSE IF (name.EQ.'NSTPMAXB  ') THEN
            READ(val, fmt=*, iostat=stat) nstpmaxb
          ELSE IF (name.EQ.'DENSEWIN  ') THEN
            READ(val, fmt=*, iostat=stat) densewin
          ELSE IF (name.EQ.'SOLPATH   ') THEN
            CALL ReadPath(val, solpath)
          ELSE IF (name.EQ.'TABPATH   ') THEN
            CALL ReadPath(val, tabpath)
          ELSE IF (name.EQ.'RELITST   ') THEN
            READ(val, fmt=*, iostat=stat) relitst
          ELSE IF (name.EQ.'RELITCG   ') THEN
            READ(val, fmt=*, iostat=stat) relitcg
          ELSE IF (name.EQ.'CSIG      ') THEN
            READ(val, fmt=*, iostat=stat) csig
          ELSE IF (name.EQ.'CSIH      ') THEN
            READ(val, fmt=*, iostat=stat) csih
          ELSE IF (name.EQ.'RHOMIN    ') THEN
            READ(val, fmt=*, iostat=stat) rhomin
          ELSE IF (name.EQ.'PHI1      ') THEN
            READ(val, fmt=*, iostat=stat) phi1
          ELSE IF (name.EQ.'PHI2      ') THEN
            READ(val, fmt=*, iostat=stat) phi2
          ELSE IF (name.EQ.'KAPPA1    ') THEN
            READ(val, fmt=*, iostat=stat) kappa1
          ELSE IF (name.EQ.'KAPPA2    ') THEN
            READ(val, fmt=*, iostat=stat) kappa2
          ELSE IF (name.EQ.'KAPPA3    ') THEN
            READ(val, fmt=*, iostat=stat) kappa3
          ELSE IF (name.EQ.'KAPPA4    ') THEN
            READ(val, fmt=*, iostat=stat) kappa4
          ELSE IF (name.EQ.'EPSCG1A   ') THEN
            READ(val, fmt=*, iostat=stat) epscg1a
          ELSE IF (name.EQ.'EPSCG1R   ') THEN
            READ(val, fmt=*, iostat=stat) epscg1r
          ELSE IF (name.EQ.'EPSCG2A   ') THEN
            READ(val, fmt=*, iostat=stat) epscg2a
          ELSE IF (name.EQ.'EPSCG2R   ') THEN
            READ(val, fmt=*, iostat=stat) epscg2r
          ELSE IF (name.EQ.'EPSCG3A   ') THEN
            READ(val, fmt=*, iostat=stat) epscg3a
          ELSE IF (name.EQ.'EPSCG3R   ') THEN
            READ(val, fmt=*, iostat=stat) epscg3r
          ELSE IF (name.EQ.'EPSCG4A   ') THEN
            READ(val, fmt=*, iostat=stat) epscg4a
          ELSE IF (name.EQ.'EPSCG4R   ') THEN
            READ(val, fmt=*, iostat=stat) epscg4r
          ELSE IF (name.EQ.'EPSST1    ') THEN
            READ(val, fmt=*, iostat=stat) epsst1
          ELSE IF (name.EQ.'EPSST2    ') THEN
            READ(val, fmt=*, iostat=stat) epsst2
          ELSE IF (name.EQ.'EPSST3    ') THEN
            READ(val, fmt=*, iostat=stat) epsst3
          ELSE IF (name.EQ.'ZETA1     ') THEN
            READ(val, fmt=*, iostat=stat) zeta1
          ELSE IF (name.EQ.'ZETA2     ') THEN
            READ(val, fmt=*, iostat=stat) zeta2
          ELSE IF (name.EQ.'ZETA3     ') THEN
            READ(val, fmt=*, iostat=stat) zeta3
          ELSE IF (name.EQ.'ALPHAR    ') THEN
            READ(val, fmt=*, iostat=stat) alphar
          ELSE IF (name.EQ.'ALPHAI    ') THEN
            READ(val, fmt=*, iostat=stat) alphai
          ELSE IF (name.EQ.'ALPHAB    ') THEN
            READ(val, fmt=*, iostat=stat) alphab
          ELSE IF (name.EQ.'ALPHAS    ') THEN
            READ(val, fmt=*, iostat=stat) alphas
          ELSE IF (name.EQ.'ETA1      ') THEN
            READ(val, fmt=*, iostat=stat) eta1
          ELSE IF (name.EQ.'ETA2      ') THEN
            READ(val, fmt=*, iostat=stat) eta2
          ELSE IF (name.EQ.'ETA3      ') THEN
            READ(val, fmt=*, iostat=stat) eta3
          ELSE IF (name.EQ.'MINDELTA  ') THEN
            READ(val, fmt=*, iostat=stat) mindelta
          ELSE IF (name.EQ.'INFDELTA  ') THEN
            READ(val, fmt=*, iostat=stat) infdelta
          ELSE IF (name.EQ.'DELTA0    ') THEN
            READ(val, fmt=*, iostat=stat) delta0
          ELSE IF (name.EQ.'THETAR    ') THEN
            READ(val, fmt=*, iostat=stat) thetar
          ELSE IF (name.EQ.'CLS1      ') THEN
            READ(val, fmt=*, iostat=stat) cls1
          ELSE IF (name.EQ.'CLS2      ') THEN
            READ(val, fmt=*, iostat=stat) cls2  
          ELSE IF (name.EQ.'STEPMAXB  ') THEN
            READ(val, fmt=*, iostat=stat) stepmaxb
          ELSE IF (name.EQ.'EPSCHOL   ') THEN
            READ(val, fmt=*, iostat=stat) EpsChol
          ELSE IF (name.EQ.'DENSEFRC  ') THEN
            READ(val, fmt=*, iostat=stat) DenseFrc
          ELSE IF (name.EQ.'LBDMAX    ') THEN
            READ(val, fmt=*, iostat=stat) lbdmax
          ELSE IF (name.EQ.'MINSTEP   ') THEN
            READ(val, fmt=*, iostat=stat) minstep
          ELSE IF (name.EQ.'USECHOL   ') THEN
            READ(val, fmt=*, iostat=stat) UseChol
          ELSE IF (name.EQ.'UPDCHOL   ') THEN
            READ(val, fmt=*, iostat=stat) UpdChol
          ENDIF 

          CALL ReadName(40, name, val, stat)

        ENDDO

        CLOSE(unit=40, iostat=stat)

      ELSE

        WRITE(*,*)'File dci.spc not found. Using default parameters'

      ENDIF

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE InitialValues(m, n, iter, itssmll, trest, tsteih, trej, 
     *                         tSoc, tbfgs, nFun, nConstr, nGrad, 
     *                         nJacob, nHprod, nSoc, nbfgs, ierr, vflag, 
     *                         delta0, relitst, relitcg, maxitst, 
     *                         minitst, maxitcg, minitcg, csig, csih, 
     *                         epsst2, epscg1a, epscg1r, epscg2a,  
     *                         epscg3a, epscg4a, infdelta, minstep, fx,  
     *                         hnorm, hnormi, snorm, Ln, Lref, DLh, DLv, 
     *                         gpnorm, gpnormi, ngp, engp, rhomax, rho, 
     *                         vdelta, hdelta, x, h, g, gp, lambda, A, 
     *                         Airow, Apcol, Diag, Mval, Mirow, Mpcol, 
     *                         Mcnel, perm, iperm, permavail, EpsChol, 
     *                         DenseWin, DenseFrc, DenseRow, UseChol, 
     *                         UpdChol, lbdmax, prcndiag, phi1, Prec, 
     *                         tmp1, tmp2, tmp3, tmp4, tmp5)
C    *                         Fun, Grad, Constr, Jacob)

      IMPLICIT NONE

!     Parameters.

      INTEGER m, n, iter, itssmll, trest, tsteih, trej, tSoc, tbfgs
      INTEGER nFun, nConstr, nGrad, nJacob, nHprod, nSoc, nbfgs, ierr
      INTEGER vflag, maxitst, minitst, maxitcg, minitcg, prcndiag
      INTEGER DenseWin, DenseRow
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      INTEGER Mirow(*)
      INTEGER Mpcol(m+1)
      INTEGER Mcnel(m), perm(m), iperm(m)
      REAL*8  delta0, relitst, relitcg, csig, csih, epsst2, epscg1a
      REAL*8  epscg1r, epscg2a, epscg3a, epscg4a, infdelta, minstep, fx
      REAL*8  hnorm, hnormi, snorm, Ln, Lref, DLh, DLv
      REAL*8  gpnorm, gpnormi, ngp, engp, phi1, rhomax, rho
      REAL*8  vdelta, hdelta, EpsChol, DenseFrc, lbdmax
      REAL*8  x(n), g(n), gp(n), tmp5(n) 
      REAL*8  h(m), lambda(m), Diag(m)
      REAL*8  tmp1(m), tmp2(m), tmp3(m), tmp4(m)
      REAL*8  A(*)
      REAL*8  Mval(*)
      REAL*8  Prec(*)
      LOGICAL UseChol, UpdChol, permavail

!     Local variables.

      INTEGER naflag
      LOGICAL LimLbd
      REAL*8  epscg, mu

!     Functions called by the routine.

      INTEGER idamax
      REAL*8  CuterFun, dnrm2, ddot
      
!     Setting up some obvious initial values.

      iter      = 0 
      itssmll   = 0
      trest     = 0
      tsteih    = 0
      trej      = 0
      tSoc      = 0
      tbfgs     = 0
      nFun      = 1
      nGrad     = 1
      nHprod    = 0
      nSoc      = 0
      nbfgs     = 0
      ierr      = 0
      vflag     = 0
      Lref      = 1.0D100
      DLh       = 1.0D100
      DLv       = 0.0D0
      snorm     = 0.0D0
      permavail =.FALSE.

!     Adjusting some constants according to the precision
!     required for the solution.
      
      maxitst = MAX(minitst, MIN(maxitst, INT(DFLOAT(n)*relitst+0.5D0))) 
      maxitcg = MAX(minitcg, MIN(maxitcg, INT(DFLOAT(m)*relitcg+0.5D0)))
      epsst2  = epsst2*csig
      epscg1a = epscg1a*csig
      epscg2a = epscg2a*csih
      epscg3a = epscg3a*csih
      epscg4a = epscg4a*csih
      infdelta= infdelta*csih
      minstep = minstep*DMIN1(csih, csig)

!     Computing f(x), g(x), h(x), ||h(x)||, A(x), Ln, gp, lambda  and
!     possibly the Cholesky decomposition of A*A'.

      fx = CuterFun(n, x)

      IF (m.GT.0) THEN

        CALL CuterConstr(m, n, x, h)        
        hnorm = dnrm2(m, h, 1)
        hnormi = DABS(h(idamax(m, h, 1)))
        CALL CuterJacob(m, n, x, A, Airow, Apcol)
        CALL CuterGrad(n, x, g)
        nConstr = 1
        nJacob  = 1

        IF (UseChol) THEN

C         Computing the Cholesky decomposition of A*A'.

          CALL choldec(m, n, A, Airow, Apcol, Diag, Mval, Mirow, 
     $                 Mpcol, Mcnel, perm, iperm, EpsChol, DenseWin,
     $                 DenseFrc, DenseRow, UpdChol, permavail)

C         Computing gp and lambda.

          CALL MtVprod(n, A, lambda, gp, Airow, Apcol)
          CALL daxpy(n, 1.0D0, g, 1, gp, 1)
          naflag = 0
c$$$          LimLbd = .TRUE.
c$$$          CALL NAprojCh(m, n, A, Airow, Apcol, Diag, Mval, Mirow, 
c$$$     $                  Mpcol, Mcnel, perm, iperm, DenseRow, lambda, 
c$$$     $                  LimLbd, lbdmax, g, gp, naflag)

        ELSE

C         Computing gp and lambda.

          epscg = epscg1a
          mu = 0.0D0
          LimLbd = .TRUE.
          CALL NAprojCG(m, n, prcndiag, .TRUE., maxitcg, epscg,
     $                  mu, A, Airow, Apcol, Prec, tmp1, tmp2, 
     $                  tmp3, tmp4, tmp5, lambda, LimLbd, lbdmax,
     $                  g, gp, naflag)

        ENDIF
 
        IF (naflag.NE.0) THEN
          WRITE(*,*)'InitialValues warning: NAproj has failed.'
        ENDIF

        Ln = fx + ddot(m, lambda, 1, h, 1)

      ELSE
        CALL CuterGrad(n, x, g)
        hnorm  = 0.0D0
        hnormi = 0.0D0
        Ln = fx
        CALL dcopy(n, g, 1, gp, 1)
        nConstr = 0
        nJacob  = 0
      ENDIF

!     Computing gpnorm, ngp, rhomax, rho, vdelta and hdelta.

      gpnorm  = dnrm2(n, gp, 1)
      gpnormi = DABS(gp(idamax(n, gp, 1)))
      ngp     = gpnorm/(dnrm2(n, g, 1) + 1.0D0)
      engp    = ngp
!     rhomax  = MAX(csih, MAX(2.1D0*hnorm, 10.0D0*ngp))
      rhomax  = MAX(csih, MAX(5.1D0*hnorm, 50.0D0*ngp))
      rho     = DMIN1(phi1*rhomax*ngp, 0.75D0*rhomax);
      IF ((rho.LT.csih).AND.(hnorm.GT.(100.0D0*csih))) THEN
        rho = 0.1D0*hnorm;
      ELSE IF (ngp.LE.(5.0D0*csig)) THEN
        rho = csih
      ENDIF
      vdelta  = MAX(10.0D0*dnrm2(n, x, 1), delta0)
      hdelta  = vdelta

      RETURN      
      END

C     ------------------------------------------------------------------
      
      SUBROUTINE KPrint(iter, nsteih, nrest, nrej, nsoc, nbfgs, fx, 
     *                  hnormc, hnorm, gpnorm, lnorm, rho, rhomax, 
     *                  snorm, vdelta, hdelta, Ln, DLh, DLv)

      IMPLICIT NONE

C     Parameters.

      INTEGER iter, nsteih, nrest, nrej, nsoc, nbfgs
      REAL*8  fx, hnormc, hnorm, gpnorm, rho, rhomax, snorm, lnorm
      REAL*8  vdelta, hdelta, Ln, DLh, DLv
      
      WRITE(*,*)
      WRITE(*,*)'Iteration       : ', iter
      WRITE(*,*)
      WRITE(*,*)'CG iterations   : ', nsteih
      WRITE(*,*)'Restorations    : ', nrest
      WRITE(*,*)'Rejected steps  : ', nrej
      WRITE(*,*)'2nd. order corr.: ', nsoc
      WRITE(*,*)'BFGS iterations : ', nbfgs
      WRITE(*,*)'f(x)            : ', fx
      WRITE(*,*)'||h(xc)||       : ', hnormc
      WRITE(*,*)'||h(x)||        : ', hnorm
      WRITE(*,*)'||gp(x)||       : ', gpnorm
      WRITE(*,*)               
      WRITE(*,*)'||xnew - xc||   : ', snorm
      WRITE(*,*)'||lambda||      : ', lnorm
      WRITE(*,*)'rho             : ', rho
      WRITE(*,*)'rhomax          : ', rhomax
      WRITE(*,*)'vdelta          : ', vdelta
      WRITE(*,*)'hdelta          : ', hdelta
      WRITE(*,*)'L(x)            : ', Ln
      WRITE(*,*)'L(x) - L(xc)    : ', DLh
      WRITE(*,*)'L(xc)-L(x(k-1)) : ', DLv
      WRITE(*,*)
      WRITE(*,*) '-------- xxxxxxxxxxxxxxxxxxx --------       '

      RETURN
      END
 
C     ------------------------------------------------------------------

      SUBROUTINE LastPrint(ierr, iter, tsteih, trest, trej, nFun, 
     *                     nConstr, nGrad, nJacob, nSoc, nbfgs, fx, 
     *                     hnormi, gpnorm, gphinorm, rho, rhomax, snorm, 
     *                     vdelta, hdelta, Ln, DLh, DLv, eltime)
                
      IMPLICIT NONE

C     Parameters.

      INTEGER ierr, iter, tsteih, trest, trej
      INTEGER nFun, nConstr, nGrad, nJacob, nSoc, nbfgs
      REAL*8  fx, hnormi, gpnorm, gphinorm, rho, rhomax
      REAL*8  snorm, vdelta, hdelta, Ln, DLh, DLv
      REAL*4  eltime
      
      WRITE(*,*)
      WRITE(*,*)'DCI Summary.'
      WRITE(*,*)

C     Reporting the stopping criterion.

      IF      (ierr.EQ.0)  THEN
	WRITE(*,*)'The algorithm has converged.'
      ELSE IF (ierr.EQ.1)  THEN
	WRITE(*,*)'rhomax became too short.'
      ELSE IF (ierr.EQ.2)  THEN
	WRITE(*,*)'The maximum number of iterations was reached.'
      ELSE IF (ierr.EQ.3)  THEN
	WRITE(*,*)'The maximum number of restorations was reached.'
      ELSE IF (ierr.EQ.4)  THEN
	WRITE(*,*)'The restoration has failed.'
      ELSE IF (ierr.EQ.5)  THEN
	WRITE(*,*)'The step became too short.'
      ENDIF

C     Printing all the relevant information.

      WRITE(*,*)
      WRITE(*,*)'Iterations     : ', iter
      WRITE(*,*)'Steihaug''s iter: ', tsteih
      WRITE(*,*)'Restorations   : ', trest
      WRITE(*,*)'Rejected steps : ', trej
      WRITE(*,*)'Function evals.: ', nFun
      WRITE(*,*)'Constr.  evals.: ', nConstr
      WRITE(*,*)'Gradient evals.: ', nGrad
      WRITE(*,*)'Jacobian evals.: ', nJacob
      WRITE(*,*)'2nd order Corr.: ', nSoc
      WRITE(*,*)'BFGS iterations: ', nbfgs
      WRITE(*,*)'f(x)           : ', fx
      WRITE(*,*)'||h(x)||       : ', hnormi
      WRITE(*,*)'||gp(x)||      : ', gpnorm
      WRITE(*,*)'||At*h||       : ', gphinorm
      WRITE(*,*)'Elapsed time   : ', eltime
      WRITE(*,*)               
      WRITE(*,*)'||xnew - xc||  : ', snorm
      WRITE(*,*)'rho            : ', rho
      WRITE(*,*)'rhomax         : ', rhomax
      WRITE(*,*)'vdelta         : ', vdelta
      WRITE(*,*)'hdelta         : ', hdelta
      WRITE(*,*)'L(x)           : ', Ln
      WRITE(*,*)'L(x) - L(xc)   : ', DLh
      WRITE(*,*)'L(xc)-L(x(k-1)): ', DLv
      WRITE(*,*)
      WRITE(*,*) '-------- xxxxxxxxxxxxxxxxxxxxxxx --------       '

      RETURN
      END

C     ------------------------------------------------------------------  

      SUBROUTINE SaveSolution(solprnt, solpath, m, n, x, lambda, fx, h, 
     &                        hnormi, iter, ierr, eltime, tRej, tRest,
     &                        nFun, nConstr, nGrad, nJacob, nHprod,  
     &                        nSoc, nbfgs, tSteih, maxrest, maxitst,  
     &                        maxitcg, maxit, maxssmll, prcndiag, csih,  
     &                        csig, phi1, phi2, epsst1, epsst2, epsst3, 
     &                        epscg1a, epscg1r, epscg2a, epscg2r, 
     &                        epscg3a, epscg3r, epscg4a, epscg4r, eta1, 
     &                        eta2, kappa1, kappa2, kappa3, kappa4, 
     &                        zeta1, zeta2, zeta3, alphar, alphai, 
     &                        mindelta, infdelta, minstep)
      
      IMPLICIT NONE

      INTEGER nmax, mmax
                     
      PARAMETER (nmax=60000, mmax=40000)

C     Subroutine parameters.

      INTEGER       solprnt, m, n, iter, ierr
      INTEGER       maxrest, maxitst, maxitcg, maxit, maxssmll, prcndiag
      INTEGER       tSteih, tRej, tRest, nFun, nConstr, nGrad, nJacob
      INTEGER       nHprod, nSoc, nbfgs
      REAL*8        fx, hnormi, csih, csig, phi1, phi2
      REAL*8        epscg1a, epscg1r, epscg2a, epscg2r, epscg3a, epscg3r
      REAL*8        epscg4a, epscg4r, eta1, eta2, kappa1, kappa2, kappa3
      REAL*8        kappa4, zeta1, zeta2, zeta3, alphar, alphai
      REAL*8        mindelta, infdelta, minstep, epsst1, epsst2, epsst3
      REAL*8        x(n)
      REAL*8        Lambda(m), h(m)
      REAL          eltime
      CHARACTER*200 solpath

C     Variables passed by COMMON statements.

      INTEGER       nt, nfix
      INTEGER       ind(nmax)
      REAL*8        xfix(nmax)

C     Local variables.

      INTEGER       i, ioerr
      CHARACTER*42  prob
      CHARACTER*10  varn(nmax)
      CHARACTER*10  constrn(mmax)
      CHARACTER*255 Arq

C     Functions called by the routine.

      INTEGER       idamax

C     COMMON statements.

      COMMON / dcitotvar / nt
      COMMON / dcifixvar / nfix, ind
      COMMON / dcifixsol / xfix

C     Reading problem, variable and constraint names.

      CALL cnames(nt,m,prob,varn,constrn)

C     Opening the output file.

      WRITE(Arq,*)solpath(1:INDEX(solpath,' ')-1),'/',
     *            Prob(1:INDEX(Prob,' ')-1),'.sol'
      IF (Arq(1:1).EQ.' ') Arq = Arq(2:LEN(Arq))
      OPEN(50,file=Arq(1:INDEX(Arq,' ')-1),iostat=ioerr)
      IF (ioerr.NE.0) THEN
        WRITE(*,*)'Unable to open solution file'
        RETURN
      ENDIF

C     Expanding the solution vector if there are fixed variables.

      IF (nfix.GT.0) THEN
        CALL vexpn(x, xfix, ind, 1, n)
      ELSE
        CALL dcopy(n, x, 1, xfix, 1)
      ENDIF

C     Writing problem name.

      WRITE(50,*)
      WRITE(50,*)' CUTE problem                      = ',Prob(1:10)
      WRITE(50,*)

C     Reporting the stopping criterion.

      IF      (ierr.EQ.0)  THEN
	WRITE(50,*)' The algorithm has converged.'
      ELSE IF (ierr.EQ.1)  THEN
	WRITE(50,*)' rhomax became too short.'
      ELSE IF (ierr.EQ.2)  THEN
	WRITE(50,*)' The maximum number of iterations was reached.'
      ELSE IF (ierr.EQ.3)  THEN
	WRITE(50,*)' The maximum number of restorations was reached.'
      ELSE IF (ierr.EQ.4)  THEN
	WRITE(50,*)'The restoration has failed.'
      ELSE IF (ierr.EQ.5)  THEN
	WRITE(50,*)'The step became too short.'
      ENDIF

C     saving solution only of convergence was atained.

      IF (ierr.EQ.0) THEN
      
	WRITE(50,*)
        WRITE(50,*)' Elapsed time                       = ',eltime
	WRITE(50,*)' Objective function value           = ',fx
	WRITE(50,*)' Maximum infeasibility of solution  = ',hnormi
	WRITE(50,*)' Number of iterations               = ',Iter
	WRITE(50,*)' Number of rejected steps           = ',tRej
	WRITE(50,*)' Number of restorations             = ',tRest
	WRITE(50,*)' Number of function evaluations     = ',nFun
	WRITE(50,*)' Number of constraint evaluations   = ',nConstr
	WRITE(50,*)' Number of gradient evaluations     = ',nGrad
	WRITE(50,*)' Number of Jacobian evaluations     = ',nJacob
	WRITE(50,*)' Number of Hessian-vector products  = ',nHprod
	WRITE(50,*)' Number of Steihaug iterations      = ',tSteih
	WRITE(50,*)' Number of second order corrections = ',nSoc
	WRITE(50,*)' Number of BFGS iterations          = ',nbfgs

C       Saving parameters read from disk.

        WRITE(50,*)
        WRITE(50,*)' DCI parameters: '
        WRITE(50,*)
        WRITE(50,*)' MAXIT    = ',maxit
        WRITE(50,*)' MAXREST  = ',maxrest
        WRITE(50,*)' MAXITST  = ',maxitst
        WRITE(50,*)' MAXITCG  = ',maxitcg
        WRITE(50,*)' MAXSSMLL = ',maxssmll
        WRITE(50,*)' PRCNDIAG = ',prcndiag
        WRITE(50,*)' CSIH     = ',csih
        WRITE(50,*)' CSIG     = ',csig
        WRITE(50,*)' PHI1     = ',phi1
        WRITE(50,*)' PHI2     = ',phi2
        WRITE(50,*)' KAPPA1   = ',kappa1
        WRITE(50,*)' KAPPA2   = ',kappa2
        WRITE(50,*)' KAPPA3   = ',kappa3
        WRITE(50,*)' KAPPA4   = ',kappa4
        WRITE(50,*)' EPSST1   = ',epsst1
        WRITE(50,*)' EPSST2   = ',epsst2
        WRITE(50,*)' EPSST3   = ',epsst3
        WRITE(50,*)' EPSCG1A  = ',epscg1a
        WRITE(50,*)' EPSCG1R  = ',epscg1r
        WRITE(50,*)' EPSCG2A  = ',epscg2a
        WRITE(50,*)' EPSCG2R  = ',epscg2r
        WRITE(50,*)' EPSCG3A  = ',epscg3a
        WRITE(50,*)' EPSCG3R  = ',epscg3r
        WRITE(50,*)' EPSCG4A  = ',epscg4a
        WRITE(50,*)' EPSCG4R  = ',epscg4r
        WRITE(50,*)' ZETA1    = ',zeta1
        WRITE(50,*)' ZETA2    = ',zeta2
        WRITE(50,*)' ZETA3    = ',zeta3
        WRITE(50,*)' ALPHAR   = ',alphar
        WRITE(50,*)' ALPHAI   = ',alphai
        WRITE(50,*)' ETA1     = ',eta1
        WRITE(50,*)' ETA2     = ',eta2
        WRITE(50,*)' MINDELTA = ',mindelta
        WRITE(50,*)' INFDELTA = ',infdelta
        WRITE(50,*)' MINSTEP  = ',minstep
      
	IF (solprnt.GT.0) THEN

C         Imprimindo a solucao primal.

	  WRITE(50,*)
	  WRITE(50,*)' Primal solution: '
	  WRITE(50,*)
	  WRITE(50,*)' Variable           Value'
	  DO i=1,nt
	    WRITE(50,200)varn(i),xfix(i)
  200       FORMAT(1X,A10,D20.8)
	  ENDDO

C         Imprimindo a solucao dual.

	  WRITE(50,*)
	  WRITE(50,*)' Dual solution: '
	  WRITE(50,*)
	  WRITE(50,*)' Constraint            h                lambda'
	  DO i=1,m
	    WRITE(50,300)constrn(i),h(i),lambda(i)
  300       FORMAT(1X,A10,2D20.8)
	  ENDDO
	ENDIF  
      
      ENDIF
      
      CLOSE(50)

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE SaveTable(tabpath, m, ierr, iter, tRej, tRest, nFun, 
     $                     nJacob, nHprod, nSoc, nbfgs, tSteih, eltime, 
     $                     fx, hnormi)
      
      IMPLICIT NONE

      INTEGER nmax, mmax
                     
      PARAMETER (nmax=60000, mmax=40000)

C     Subroutine parameters.

      INTEGER       m, ierr, iter, tRej, tRest, tSteih
      INTEGER       nFun, nJacob, nHprod, nSoc, nbfgs
      REAL*8        fx
      REAL          eltime
      CHARACTER*200 tabpath

C     Variables passed by COMMON statements.

      INTEGER       nt

C     Local variables.

      INTEGER       ioerr
      REAL*8        hnormi
      CHARACTER*42  prob
      CHARACTER*10  varn(nmax)
      CHARACTER*10  constrn(mmax)
      CHARACTER*255 Arq

C     Functions called by the routine.

      INTEGER       idamax

C     COMMON statement.

      COMMON / dcitotvar / nt

C     Reading problem, variable and constraint names.

      CALL cnames(nt,m,prob,varn,constrn)

C     Opening the output file.

      WRITE(Arq,*)tabpath(1:INDEX(tabpath,' ')-1),'/',
     *            Prob(1:INDEX(Prob,' ')-1),'.out'
      IF (Arq(1:1).EQ.' ') Arq = Arq(2:LEN(Arq))
      OPEN(50,file=Arq,iostat=ioerr)
      IF (ioerr.NE.0) THEN
        WRITE(*,*)'Unable to open output file'
        RETURN
      ENDIF

C     Saving problem name, f(x*), ||h(x*)||, stopping criterion, 
C     execution time, # iterations, # rejected steps, # restorations, 
C     # function evaluations, # Jacobian evaluations, # Hessian-vector
C     products, # Steihaug iterations.

      WRITE(50,100)Prob, fx, hnormi, ierr, eltime, iter, tRej, tRest, 
     $             nFun, nJacob, nHprod, tSteih, nSoc, nbfgs 

  100 FORMAT(A10,' ',D12.5,' ',D8.1,' ',I2,' ',F8.2,' ',I5,' ',
     $       I4,' ',I4,' ',I5,' ',I6,' ',I6,' ',I5,' ',I5,' ',I5)
      
      CLOSE(50)

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE MatPrint(n, A, Airow, Apcol)

      IMPLICIT NONE

      INTEGER n
      INTEGER i, j
      INTEGER Airow(*), Apcol(*)
      REAL*8  A(*)

      OPEN(unit=36, file='matrix.m')

      WRITE(36,*)'mat = ['

      DO i = 1, n
        DO j = Apcol(i), Apcol(i+1)-1
          WRITE(36,*)Airow(j),' ', i, ' ', A(j)
        ENDDO
      ENDDO

      WRITE(36,*)'];'

      CLOSE(36)

      RETURN
      END

C     ------------------------------------------------------------------  

      SUBROUTINE TriMatPrint(n, L, Lirow, Lpcol, Lcnel)

      IMPLICIT NONE

      INTEGER n
      INTEGER i, j
      INTEGER Lpcol(*), Lirow(*), Lcnel(*)
      REAL*8  L(*)

      OPEN(unit=36, file='matL.m')

      WRITE(36,*)'matL = ['

      DO i = 1, n
        DO j = Lpcol(i), Lpcol(i)+Lcnel(i)-1
          WRITE(36,*)Lirow(j),' ', i, ' ', L(j)
        ENDDO
      ENDDO

      WRITE(36,*)'];'

      CLOSE(36)

      RETURN
      END

C     ------------------------------------------------------------------  

      SUBROUTINE VecPrint(n, v)

      IMPLICIT NONE

      INTEGER n
      INTEGER i
      REAL*8  v(*)

      OPEN(unit=37, file='vector.m')

      WRITE(37,*)'vec = ['

      DO i = 1, n
        WRITE(37,*)v(i)
      ENDDO

      WRITE(37,*)'];'

      CLOSE(37)

      RETURN
      END

C     ------------------------------------------------------------------  

      SUBROUTINE IVecPrint(n, v)

      IMPLICIT NONE

      INTEGER n
      INTEGER i
      INTEGER v(*)

      OPEN(unit=37, file='ivector.m')

      WRITE(37,*)'ivec = ['

      DO i = 1, n
        WRITE(37,*)v(i)
      ENDDO

      WRITE(37,*)'];'

      CLOSE(37)

      RETURN
      END
