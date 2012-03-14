      PROGRAM TesteChol

      IMPLICIT NONE

      INTEGER nmax, mmax, amax, lmax
                     
      PARAMETER (nmax=30000, mmax=20000, amax=1000000, lmax = 210000)

      INTEGER m, n, iter, flag, DenseWin, DenseRow
      INTEGER Apcol(nmax+1)
      INTEGER Airow(amax)
      INTEGER Mirow(amax)
      INTEGER Mpcol(nmax+1)
      INTEGER perm(nmax), iperm(nmax), Mcnel(nmax)
      REAL*8  rnorm, relres, EpsChol, DenseFrc
      REAL*8  A(amax)
      REAL*8  Mval(amax)
      REAL*8  r(mmax), b(mmax), xs(mmax), x(mmax), tmp1(nmax)
      REAL*8  y(nmax), Diag(nmax)
      REAL*8  dnrm2
      LOGICAL upd

C     Usando a rotina Msolve para resolver o sistema 
C     (A*A')*x = b usando a decomposição de Cholesky
C     fornecida por dcichol.f.
C
C     A = [0 0  0   0  -4 -17  0  0   0  1  3  0   0 0  0;
C          0 0  0   0 -11  12  0  0  12  0  0  0   0 0  0;
C          0 0  0   0   0   0  0  0   0  0  3  0   0 2  0;
C          0 0  0   0   0  -2  0  0   0  0  0  0   0 7  0;
C          0 0 -6   0  22   0 -1  1   0 11  1 -1   0 0  0;
C          0 0  0   0   0   0  0 -8   0  0  0  0   0 0  0;
C          0 0  3 -13   0   0  7  0   0  0 16  0  -7 0  0;
C          0 0  9   0   0  13  0  0 -16  0  0  0 -14 0  0;
C          0 0  6   0  -4   0  0  0   7  0  0  8   7 0 13;
C          0 0  0   0   0   0  0  7   0  0 12  0 -12 0  0]
C
C     x = (1:10)'
C
C     b = [-1140 216 815 32 1135 -216 7352 6100 754 5793]'

      m = 10
      n = 15

      A(1)  =  -6.0D0
      A(2)  =   3.0D0
      A(3)  =   9.0D0
      A(4)  =   6.0D0
      A(5)  = -13.0D0
      A(6)  =  -4.0D0
      A(7)  = -11.0D0
      A(8)  =  22.0D0
      A(9)  =  -4.0D0
      A(10) = -17.0D0
      A(11) =  12.0D0
      A(12) =  -2.0D0
      A(13) =  13.0D0
      A(14) =  -1.0D0
      A(15) =   7.0D0
      A(16) =   1.0D0
      A(17) =  -8.0D0
      A(18) =   7.0D0
      A(19) =  12.0D0
      A(20) = -16.0D0
      A(21) =   7.0D0
      A(22) =   1.0D0
      A(23) =  11.0D0
      A(24) =   3.0D0
      A(25) =   3.0D0
      A(26) =   1.0D0
      A(27) =  16.0D0
      A(28) =  12.0D0
      A(29) =  -1.0D0
      A(30) =   8.0D0
      A(31) =  -7.0D0
      A(32) = -14.0D0
      A(33) =   7.0D0
      A(34) = -12.0D0
      A(35) =   2.0D0
      A(36) =   7.0D0
      A(37) =  13.0D0

      Airow(1)  =  5
      Airow(2)  =  7
      Airow(3)  =  8
      Airow(4)  =  9
      Airow(5)  =  7
      Airow(6)  =  1
      Airow(7)  =  2
      Airow(8)  =  5
      Airow(9)  =  9
      Airow(10) =  1
      Airow(11) =  2
      Airow(12) =  4
      Airow(13) =  8
      Airow(14) =  5
      Airow(15) =  7
      Airow(16) =  5
      Airow(17) =  6
      Airow(18) =  10
      Airow(19) =  2
      Airow(20) =  8
      Airow(21) =  9
      Airow(22) =  1
      Airow(23) =  5
      Airow(24) =  1
      Airow(25) =  3
      Airow(26) =  5
      Airow(27) =  7
      Airow(28) =  10
      Airow(29) =  5
      Airow(30) =  9
      Airow(31) =  7
      Airow(32) =  8
      Airow(33) =  9
      Airow(34) =  10
      Airow(35) =  3
      Airow(36) =  4
      Airow(37) =  9

      Apcol(1)  =  1
      Apcol(2)  =  1
      Apcol(3)  =  1
      Apcol(4)  =  5
      Apcol(5)  =  6 
      Apcol(6)  = 10
      Apcol(7)  = 14
      Apcol(8)  = 16
      Apcol(9)  = 19
      Apcol(10) = 22
      Apcol(11) = 24
      Apcol(12) = 29
      Apcol(13) = 31
      Apcol(14) = 35
      Apcol(15) = 37
      Apcol(16) = 38

      b(1)  = -1140D0
      b(2)  =   216D0
      b(3)  =   815D0
      b(4)  =    32D0
      b(5)  =  1135D0
      b(6)  =  -216D0
      b(7)  =  7352D0
      b(8)  =  6100D0
      b(9)  =   754D0
      b(10) =  5793D0

      xs(1)  =  1.0D0
      xs(2)  =  2.0D0
      xs(3)  =  3.0D0
      xs(4)  =  4.0D0
      xs(5)  =  5.0D0
      xs(6)  =  6.0D0
      xs(7)  =  7.0D0
      xs(8)  =  8.0D0
      xs(9)  =  9.0D0
      xs(10) = 10.0D0

C     Conferindo se (A*A')*xs = b.

      CALL MtVprod(n, A, xs, y, Airow, Apcol)
      CALL MVprod(m, n, A, y, r, Airow, Apcol, 0.0D0)
      CALL daxpy(m, -1.0D0, b, 1, r, 1)
      rnorm = dnrm2(m, r, 1)

      WRITE(*,*)'Conferindo se (A*At)*x = b.'
      WRITE(*,*)
      WRITE(*,*)'   r = ',rnorm
      WRITE(*,*)

      DenseFrc = -1.0D0
      DenseWin = 10000
      EpsChol = 1.0D-20
      upd = .FALSE.

      CALL choldec(m, n, A, Airow, Apcol, diag, Mval, Mirow, 
     $             Mpcol, Mcnel, perm, iperm, EpsChol, 
     $             DenseWin, DenseFrc, DenseRow, upd)


      CALL dcopy(m, b, 1, x, 1)
      CALL MsolveCh(m, diag, Mval, Mpcol, Mirow, Mcnel,
     $              perm, iperm, DenseRow, tmp1, x)


      CALL MtVprod(n, A, x, y, Airow, Apcol)
      CALL MVprod(m, n, A, y, r, Airow, Apcol, 0.0D0)
      CALL daxpy(m, -1.0D0, b, 1, r, 1)
      rnorm = dnrm2(m, r, 1)

      WRITE(*,*)'Executando Msolve.'
      WRITE(*,*)
      WRITE(*,*)'   r    = ',rnorm
      WRITE(*,*)

      STOP
      END 
