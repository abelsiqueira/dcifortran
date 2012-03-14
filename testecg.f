      PROGRAM TesteCG

      IMPLICIT NONE

      INTEGER nmax, mmax, amax, lmax
                     
      PARAMETER (nmax=30000, mmax=20000, amax=1000000, lmax = 210000)

      INTEGER m, n, maxitcg, prcndiag, iter, flag
      INTEGER Apcol(nmax+1)
      INTEGER Airow(amax)
      REAL*8  epscg, rnorm, relres, mu
      REAL*8  A(amax)
      REAL*8  Prec(lmax)
      REAL*8  r(mmax), b(mmax), xs(mmax), x(mmax)
      REAL*8  tmp1(mmax), tmp2(mmax), tmp3(mmax), tmp4(mmax)
      REAL*8  y(nmax)
      REAL*8  dnrm2
      LOGICAL prcavail, goodx0

C     Usando a rotina Msolve para resolver o sistema 
C     (A*A')*x = b usando o metodo dos gradientes conjugados
C     fornecido pela rotina Msolve de dcicg.f.
C
C     A = [ 1 -2  3  0  0  0  0  1  0  0;
C           0  3 -2  2  0  0  0  0 -1  0;
C           0  0  1  4 -1  0  0  0  0  0;
C           0  0  0  4  1  3  0  0  0 -1;
C           0  0  0  0  0  1  0  0  0  1;
C           0  0  0  0  0  5 -1 -2  0  0;
C           0  0  0  0  0  0  4  1  3  0;
C           0  0  0  0  0  0  0  0 -2  2];
C
C     x = (1:8)'
C
C     b = [-5 69 129 253 64 221 93 28]'
C
C     prcndiag fornece o numero total de diagonais de A*A' 
C     usadas pelo precondicionador.
C
C     As diagonais do precondicionador pentadiagonal sao
C     dadas pela matriz P abaixo. A primeira coluna de P
C     corresponde a diagonal -2, a terceira a diagonal 0
C     e a quinta a diagonal +2 (usando a notacao do Matlab).
C
C     P = [ 3   -12    15     0     0
C           8     6    18   -12     0
C           0    15    18     6     3
C          15     2    27    15     8
C           0     5     2     2     0
C           0    -6    30     5    15
C           0    -6    26    -6     0
C           0     0     8    -6     0]
  
      COMMON / dciprecavail / prcavail

      m = 8
      n = 10

      A(1)  =  1.0D0
      A(2)  = -2.0D0
      A(3)  =  3.0D0
      A(4)  =  3.0D0
      A(5)  = -2.0D0
      A(6)  =  1.0D0
      A(7)  =  2.0D0
      A(8)  =  4.0D0
      A(9)  =  4.0D0
      A(10) = -1.0D0
      A(11) =  1.0D0
      A(12) =  3.0D0
      A(13) =  1.0D0
      A(14) =  5.0D0
      A(15) = -1.0D0
      A(16) =  4.0D0
      A(17) =  1.0D0
      A(18) = -2.0D0
      A(19) =  1.0D0
      A(20) = -1.0D0
      A(21) =  3.0D0
      A(22) = -2.0D0
      A(23) = -1.0D0
      A(24) =  1.0D0
      A(25) =  2.0D0

      Airow(1)  = 1
      Airow(2)  = 1
      Airow(3)  = 2
      Airow(4)  = 1
      Airow(5)  = 2
      Airow(6)  = 3
      Airow(7)  = 2
      Airow(8)  = 3
      Airow(9)  = 4
      Airow(10) = 3
      Airow(11) = 4
      Airow(12) = 4
      Airow(13) = 5
      Airow(14) = 6
      Airow(15) = 6
      Airow(16) = 7
      Airow(17) = 1
      Airow(18) = 6
      Airow(19) = 7
      Airow(20) = 2
      Airow(21) = 7
      Airow(22) = 8
      Airow(23) = 4
      Airow(24) = 5
      Airow(25) = 8

      Apcol(1)  =  1
      Apcol(2)  =  2
      Apcol(3)  =  4
      Apcol(4)  =  7
      Apcol(5)  = 10
      Apcol(6)  = 12
      Apcol(7)  = 15
      Apcol(8)  = 17
      Apcol(9)  = 20
      Apcol(10) = 23
      Apcol(11) = 26

      b(1) =  -5.0D0
      b(2) =  69.0D0
      b(3) = 129.0D0
      b(4) = 253.0D0
      b(5) =  64.0D0
      b(6) = 221.0D0
      b(7) =  93.0D0
      b(8) =  28.0D0

      xs(1) =  1.0D0
      xs(2) =  2.0D0
      xs(3) =  3.0D0
      xs(4) =  4.0D0
      xs(5) =  5.0D0
      xs(6) =  6.0D0
      xs(7) =  7.0D0
      xs(8) =  8.0D0

      goodx0  = .FALSE.
      mu      = 0.0D+0
      epscg   = 1.0D-8
      maxitcg = 20

C     Conferindo se (A*A')*xs = b.

      CALL MtVprod(n, A, xs, y, Airow, Apcol)
      CALL MVprod(m, n, A, y, r, Airow, Apcol, 0.0D0)
      CALL daxpy(m, -1.0D0, b, 1, r, 1)
      rnorm = dnrm2(m, r, 1)

      WRITE(*,*)'Conferindo se (A*At)*x = b.'
      WRITE(*,*)
      WRITE(*,*)'   r = ',rnorm
      WRITE(*,*)

      prcndiag = 0
      prcavail = .FALSE.      

      CALL Msolve(m, n, maxitcg, epscg, mu, A, Airow, Apcol, prcndiag, 
     $            Prec, tmp1, tmp2, tmp3, tmp4, y, b, goodx0, x, 
     $            relres, iter, flag)

      CALL MtVprod(n, A, x, y, Airow, Apcol)
      CALL MVprod(m, n, A, y, r, Airow, Apcol, 0.0D0)
      CALL daxpy(m, mu, x, 1, r, 1)
      CALL daxpy(m, -1.0D0, b, 1, r, 1)
      rnorm = dnrm2(m, r, 1)

      WRITE(*,*)'Executando Msolve.'
      WRITE(*,*)
      WRITE(*,*)'   prec = ',prcndiag
      WRITE(*,*)'   flag = ',flag
      WRITE(*,*)'   iter = ',iter
      WRITE(*,*)'   res  = ',relres 
      WRITE(*,*)'   r    = ',rnorm
      WRITE(*,*)

      prcndiag = 1
      prcavail = .FALSE.      

      CALL Msolve(m, n, maxitcg, epscg, mu, A, Airow, Apcol, prcndiag, 
     $            Prec, tmp1, tmp2, tmp3, tmp4, y, b, goodx0, x, 
     $            relres, iter, flag)

      CALL MtVprod(n, A, x, y, Airow, Apcol)
      CALL MVprod(m, n, A, y, r, Airow, Apcol, 0.0D0)
      CALL daxpy(m, mu, x, 1, r, 1)
      CALL daxpy(m, -1.0D0, b, 1, r, 1)
      rnorm = dnrm2(m, r, 1)

      WRITE(*,*)'Executando Msolve.'
      WRITE(*,*)
      WRITE(*,*)'   prec = ',prcndiag
      WRITE(*,*)'   flag = ',flag
      WRITE(*,*)'   iter = ',iter
      WRITE(*,*)'   res  = ',relres 
      WRITE(*,*)'   r    = ',rnorm
      WRITE(*,*)

      prcndiag = 3
      prcavail = .FALSE.      

      CALL Msolve(m, n, maxitcg, epscg, mu, A, Airow, Apcol, prcndiag, 
     $            Prec, tmp1, tmp2, tmp3, tmp4, y, b, goodx0, x, 
     $            relres, iter, flag)

      CALL MtVprod(n, A, x, y, Airow, Apcol)
      CALL MVprod(m, n, A, y, r, Airow, Apcol, 0.0D0)
      CALL daxpy(m, mu, x, 1, r, 1)
      CALL daxpy(m, -1.0D0, b, 1, r, 1)
      rnorm = dnrm2(m, r, 1)

      WRITE(*,*)'Executando Msolve.'
      WRITE(*,*)
      WRITE(*,*)'   prec = ',prcndiag
      WRITE(*,*)'   flag = ',flag
      WRITE(*,*)'   iter = ',iter
      WRITE(*,*)'   res  = ',relres 
      WRITE(*,*)'   r    = ',rnorm
      WRITE(*,*)

      prcndiag = 5
      prcavail = .FALSE.      

      CALL Msolve(m, n, maxitcg, epscg, mu, A, Airow, Apcol, prcndiag, 
     $            Prec, tmp1, tmp2, tmp3, tmp4, y, b, goodx0, x, 
     $            relres, iter, flag)

      CALL MtVprod(n, A, x, y, Airow, Apcol)
      CALL MVprod(m, n, A, y, r, Airow, Apcol, 0.0D0)
      CALL daxpy(m, mu, x, 1, r, 1)
      CALL daxpy(m, -1.0D0, b, 1, r, 1)
      rnorm = dnrm2(m, r, 1)

      WRITE(*,*)'Executando Msolve.'
      WRITE(*,*)
      WRITE(*,*)'   prec = ',prcndiag
      WRITE(*,*)'   flag = ',flag
      WRITE(*,*)'   iter = ',iter
      WRITE(*,*)'   res  = ',relres 
      WRITE(*,*)'   r    = ',rnorm
      WRITE(*,*)

      STOP
      END 
