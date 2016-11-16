!     ------------------------------------------------------------------
!     PROGRAM : DCI.
!     MODULE  : DCIMain.
!     VERSION : 2.0.
!     DATE    : June, 2000.
!     CONTENTS: PROGRAM    DCIMain
!
!     External routines called by this module:
!
!     From dci     : dci.
!     From dcicute : SetUp.
!
!     ------------------------------------------------------------------

      PROGRAM DCIMain

C     Test program for the DCI routine. This program reads a problem
C     from the Cute collection and calls DCI to solve it.

C     Variables used in this program:
C
C     a) Integer variables.
C
C     n
C     m
C     iter
C     ierr
C
C     b) Real variables.
C
C     fx
C
C     c) Real vectors.
C
C     x
C     h
C     lambda

      IMPLICIT NONE
    
      EXTERNAL Fun, Constr, Grad, Jacob, Hprod

      INTEGER nmax, mmax

      PARAMETER (nmax=60000, mmax=40000)

      INTEGER      m, n, iter, ierr
      REAL*8       fx
      REAL*8       x(nmax)
      REAL*8       h(mmax), lambda(mmax)
      LOGICAL      lincon
      CHARACTER*10 Prob

C     COMMON variables

C     INTEGER      nt, nfix
C     INTEGER      ind(nmax), rind(nmax)
C     REAL*8       xfix(nmax), vfix(nmax), gv(nmax)
C     LOGICAL      gavail
C     LOGICAL      FreeV(nmax)

C     COMMON statements

C     COMMON / dcigradavail / gavail
C     COMMON / dcigradient  / gv
C     COMMON / dcifixvar    / nfix, ind
C     COMMON / dcirevind    / rind
C     COMMON / dcifixsol    / xfix
C     COMMON / dcifixvect   / vfix
C     COMMON / dcitotvar    / nt
C     COMMON / dcitotconstr / m
C     COMMON / dcifixtype   / FreeV

      CALL CutestSetUp(Prob, m, n, lincon, x, lambda)

      CALL DCI(m, n, lincon, x, lambda, fx, h, iter, ierr) 
C    $         Fun, Constr, Grad, Jacob, Hprod)

      STOP
      END

