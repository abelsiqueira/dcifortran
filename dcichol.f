!     ------------------------------------------------------------------
!     PROGRAM : DCI.
!     MODULE  : DCIChol.
!     VERSION : 2.0.
!     DATE    : June, 2000.
!     CONTENTS: subroutine NAprojCh
!               subroutine NewStepCh
!               subroutine choldec
!               subroutine MsolveCh
!
!     External routines called by this module:
!
!     From dciblas : MVprod, MtVprod, vfill, vfills, ifill, icopy,
!                    vcmpr, icmpr, axpys, dots.
!     From the blas: dcopy, dscal, daxpy, ddot.
!     From amdbar  : amdbar.
!     ------------------------------------------------------------------

      SUBROUTINE NAprojCh(m, n, A, Airow, Apcol, Diag, Mval, Mirow, 
     $                    Mpcol, Mcnel, perm, iperm, DenseRow, tmp, 
     $                    LimTmp, tmpmax, r, Pr, flag)

C     This routine computes the projection of r onto N(A).
C     Pr = r - A'*(L'\(L\(A*r))), where AAt = A*A', perm = symmmd(AAt)
C     and L = chol(AAt(perm,perm))'. See AAtfact for a description of
C     Diag, Mval, Mirow, Mpcol, Mcnel, perm and iperm. m and n are
C     the number of rows and columns of a, respectively. tmp(m) is a 
C     temporary vector. On output, tmp stores -L'\(L\(A*r)). When r is 
C     the gradient of the objective function of DCI, tmp is equal to 
C     lambda, the vector of Lagrange multipliers. LimTmp is a logical 
C     parameter that indicates if tmp is to be bounded. If LimTmp = 
C     .TRUE., tmp = max(min(tmp, tmpmax), -tmpmax). flag is set to a 
C     value greater than zero when the Cholesky factorization fails.

      IMPLICIT NONE

C     Routine parameters.

      INTEGER m, n, DenseRow, flag
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      INTEGER Mpcol(m+1)
      INTEGER Mcnel(m), perm(m), iperm(m)
      INTEGER Mirow(*)
      REAL*8  tmpmax
      REAL*8  A(*)
      REAL*8  Diag(m), tmp(m)
      REAL*8  Mval(*)
      REAL*8  r(n), Pr(n)
      LOGICAL LimTmp

C     Local variable.

      INTEGER i

      flag = 0

      IF (m.GT.0) THEN
        CALL MVprod(m, n, A, r, tmp, Airow, Apcol, 0.0D0)
        CALL MsolveCh(m, Diag, Mval, Mpcol, Mirow, Mcnel,
     $                perm, iperm, DenseRow, Pr, tmp)
        CALL dscal(m, -1.0D0, tmp, 1)
        IF (LimTmp) THEN
          DO i = 1,m
            IF (tmp(i).GT.tmpmax) tmp(i) = tmpmax
          ENDDO
          DO i = 1,m
            IF (tmp(i).LT.(-tmpmax)) tmp(i) = -tmpmax
          ENDDO
        ENDIF
        CALL MtVprod(n, A, tmp, Pr, Airow, Apcol) 
        CALL daxpy(n, 1.0D0, r, 1, Pr, 1)

      ELSE
        CALL dcopy(n, r, 1, Pr, 1)
      ENDIF

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE NewStepCh(m, n, A, Airow, Apcol, Diag, Mval, Mirow, 
     $                     Mpcol, Mcnel, perm, iperm, DenseRow, tmp,
     $                     r, s, flag)

C     This routine computes s = -A'*inv(A*A')*r. flag is set to a value 
C     greater than zero when the Cholesky factorization fails.

      IMPLICIT NONE

C     Routine parameters.

      INTEGER m, n
      INTEGER DenseRow, flag
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      INTEGER Mirow(*)
      INTEGER Mpcol(m+1)
      INTEGER Mcnel(m), perm(m), iperm(m)
      REAL*8  A(*)
      REAL*8  Diag(m), tmp(m), r(m)
      REAL*8  Mval(*)
      REAL*8  s(n)

c$$$C     Local variables.
c$$$
c$$$      REAL*8  rnorm, bnorm
c$$$
c$$$C     Functions called by the routine.
c$$$
c$$$      REAL*8  ResNorm, dnrm2

      flag = 0

      IF (m.GT.0) THEN
        CALL dcopy(m, r, 1, tmp, 1)
        CALL MsolveCh(m, Diag, Mval, Mpcol, Mirow, Mcnel,
     $                perm, iperm, DenseRow, s, tmp)

c$$$        rnorm = ResNorm(m, n, A, Airow, Apcol, tmp, r)
c$$$        bnorm = dnrm2(m, r, 1)
c$$$        IF ((rnorm/(1.0D0+bnorm)).GT.1.0D-7) THEN
c$$$          WRITE(*,*)'NAprojCh warning: inacurate solution of A*Atx = b'
c$$$          WRITE(*,*)'Residual norm: ', rnorm
c$$$          WRITE(*,*)'Dense row    : ', DenseRow
c$$$          STOP
c$$$        ENDIF

        CALL dscal(m, -1.0D0, tmp, 1)
        CALL MtVprod(n, A, tmp, s, Airow, Apcol)
      ELSE
        CALL vfill(n, 0.0D0, s, 1)
      ENDIF

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE choldec(m, n, A, Airow, Apcol, diag, Mval, Mirow, 
     $                   Mpcol, Mcnel, perm, iperm, EpsChol, 
     $                   DenseWin, DenseFrc, DenseRow, upd, permavail)

C     This routine computes the Cholesky decomposition of A*A'.

      IMPLICIT NONE

      INTEGER nmax, mmax, amax, lmax
                     
      PARAMETER (nmax=60000, mmax=40000, amax=1000000, lmax = 9000000)

C     Routine parameters.

      INTEGER m, n
      INTEGER DenseWin, DenseRow
      INTEGER Airow(amax)
      INTEGER Apcol(n+1)
      INTEGER Mirow(lmax)
      INTEGER Mpcol(m+1)
      INTEGER Mcnel(m), perm(m), iperm(m)
      REAL*8  EpsChol, DenseFrc
      REAL*8  A(amax)
      REAL*8  Mval(lmax)
      REAL*8  diag(m)
      LOGICAL upd, permavail

C     Local variables.

      INTEGER ncmpa, i
      INTEGER Aicol(amax)
      INTEGER Aprow(mmax), ape(mmax), alen(mmax), anv(mmax)
      INTEGER next(mmax), head(mmax), deg(mmax), aw(mmax)
      INTEGER aiw(lmax)
      REAL*8  row(mmax)
      LOGICAL temp(mmax)

      IF (.NOT.upd) THEN
        CALL ARowInfo(m, n, Airow, Apcol, Aicol, Aprow)
        CALL AAtform(m, n, A, Airow, Apcol, Aicol, Aprow, lmax, 
     $               diag, Mval, Mirow, Mpcol, aiw, next, temp)
        CALL icopy(Mpcol(m+1)-1, Mirow, 1, aiw, 1)
        CALL icopy(m+1, Mpcol, 1, ape, 1)
        DO i = 1, m
          alen(i) = Mpcol(i+1)-Mpcol(i)
        ENDDO
        IF (.NOT.permavail) THEN
          CALL AMDBAR (m, ape, aiw, alen, lmax, ape(m+1), anv, next,
     $                 perm, head, iperm, deg, ncmpa, aw)
          permavail = .TRUE.
        ENDIF
        CALL AAtfact(m, diag, Mval, Mpcol, Mirow, Mcnel, lmax, perm, 
     $               iperm, row, ape, alen, anv, next, head,
     $               deg, temp, EpsChol, DenseWin, DenseFrc, DenseRow)
      ELSE
        CALL AAtupdt(m, n, A, Airow, Apcol, Aicol, Aprow,
     $               diag, Mval, Mirow, Mpcol, aiw, next)
        CALL AAtfact(m, diag, Mval, Mpcol, Mirow, Mcnel, lmax, perm, 
     $               iperm, row, ape, alen, anv, next, head,
     $               deg, temp, EpsChol, DenseWin, DenseFrc, DenseRow)
      ENDIF

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE ARowInfo(m, n, irow, pcol, icol, prow) 
                                                                       
C     This routine adds to the triplet {A, irow, pcol} (that defines
C     the information about the columns of A), the pair {icol, prow}
C     that provides information about the rows of the matrix.

      IMPLICIT NONE                                                     
                                                                       
C     Routine parameters.

      INTEGER m, n
      INTEGER irow(*)
      INTEGER icol(*)
      INTEGER pcol(n+1)
      INTEGER prow(m+1)
 
C     Local variables.

      INTEGER i, j

C     Storing in prow the number of elements on each row of A.

      CALL ifill(m+1, 0, prow, 1)
      DO i = 1, pcol(n+1)-1
        prow(irow(i)) = prow(irow(i))+1                               
      ENDDO                                                             

C     Accumulating prow.                                                                       

      prow(1) = prow(1)+1
      DO i = 2, m+1
        prow(i) = prow(i-1) + prow(i)                                    
      ENDDO
                                                                       
C     Storing the column of each element in icol and decreasing prow.

      DO i = n, 1, -1                                                   
        DO j = pcol(i), pcol(i+1)-1                               
          prow(irow(j)) = prow(irow(j))-1                               
          icol(prow(irow(j))) = i                                        
        ENDDO                                                           
      ENDDO                                                             

      RETURN                                                            
      END

C     ------------------------------------------------------------------

      SUBROUTINE AAtform(m, n, A, Airow, Apcol, Aicol, Aprow, lmax, 
     $                   diag, Mval, Mirow, Mpcol, pos, above, temp)
                                                                       
C     This routine computes A*A', where A is stored by columns in the
C     triplet {A, Airow, Apcol}. On exit, the main diagonal of A*A' is
C     stored in the diag vector, while the remaining elements of the
C     matrix are stored in the triplet {Mval, Mirow, Mpcol} using the 
C     CSC format. pos(n), above(m) and temp(m) are temporary vectors.

      IMPLICIT NONE                                                     

C     Routine parameters.

      INTEGER m, n, lmax
      INTEGER Airow(*), Aicol(*)
      INTEGER Apcol(n+1)
      INTEGER Aprow(m+1), Mpcol(m+1)
      INTEGER Mirow(lmax)
      INTEGER pos(n)
      INTEGER above(m)
      REAL*8  A(*)
      REAL*8  Mval(lmax)
      REAL*8  diag(m)
      LOGICAL temp(m)

C     Local variables.

      INTEGER i, j, k, cont, col, jk
      REAL*8  t
                                                                       
      DO i = 1, m
        temp(i) = .FALSE.
      ENDDO
      CALL ifill(m, 0, above, 1)
      CALL icopy(n, Apcol, 1, pos, 1)
      Mpcol(1) = 1

      DO i = 1, m

        cont = Mpcol(i)
                                                     
C       Avoiding the diagonal element of A*A'.

        temp(i) = .TRUE.
        diag(i) = 0.0D0

C       Forming the rows of A*A'.

        DO j = Aprow(i), Aprow(i+1)-1

C         Analyzing the nonzero elements in row i of A.
                                                                       
          col = Aicol(j)
          t = A(pos(col))
          DO k = pos(col), Apcol(col+1)-1
                                                                       
C           Analyzing the nonzero elements in row "col" of A'.

            jk = Airow(k)
            IF (.NOT.temp(jk)) THEN
              temp(jk) = .TRUE.
              diag(jk) = t*A(k)
              Mirow(cont) = jk
              cont = cont+1
              above(jk) = above(jk)+1
            ELSE
              diag(jk) = diag(jk) + t*A(k)
            ENDIF

          ENDDO
          pos(col) = pos(col)+1
        ENDDO

        Mpcol(i+1) = cont + above(i)
        IF (Mpcol(i+1).GT.lmax) THEN
          WRITE(*,*)'DCI failed to store L. Increase lmax.'
          STOP
        ENDIF

C       Compressing the i-th column of L.

        DO j = Mpcol(i), cont-1
          Mval(j) = diag(Mirow(j))
          temp(Mirow(j)) = .FALSE.
        ENDDO

      ENDDO

C     Copying the elements above the diagonal of L.

      DO i = 1, m
        pos(i) = Mpcol(i+1)-above(i)
      ENDDO

      DO i = 1, m-1

        DO j = Mpcol(i), Mpcol(i+1)-above(i)-1
          k = Mirow(j)
          Mval(pos(k)) = Mval(j)
          Mirow(pos(k)) = i
          pos(k) = pos(k)+1
        ENDDO

      ENDDO

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE AAtupdt(m, n, A, Airow, Apcol, Aicol, Aprow, 
     $                   diag, Mval, Mirow, Mpcol, pos, above)
                                                                       
C     This routine recomputes A*A', where A is stored by columns in 
C     the triplet {A, Airow, Apcol}. On exit, the main diagonal of A*A'
C     is stored in the diag vector, while the remaining elements of the
C     matrix are stored in the triplet {Mval, Mirow, Mpcol} using the
C     CSC format. The structure of M is supposed to be generated by a
C     previous call to subroutine AAtform. pos(n) and above(m) are 
C     temporary vectors.

      IMPLICIT NONE                                                     

C     Routine parameters.

      INTEGER m, n
      INTEGER Airow(*), Aicol(*)
      INTEGER Apcol(n+1)
      INTEGER Aprow(m+1), Mpcol(m+1)
      INTEGER Mirow(*)
      INTEGER pos(n)
      INTEGER above(m)
      REAL*8  A(*)
      REAL*8  Mval(*)
      REAL*8  diag(m)
                         
C     Local variables.
                                              
      INTEGER i, j, k, col, jk
      REAL*8  t

      CALL vfill(m, 0.0D0, diag, 1)
      CALL ifill(m, -1, above, 1)
      CALL icopy(n, Apcol, 1, pos, 1)

      DO i = 1, m
                                                                       
C       Forming the rows of A*A'.
                                                                       
        DO j = Aprow(i), Aprow(i+1)-1
                                                                       
C         Analyzing the nonzero elements in row i of A.
                                                                       
          col = Aicol(j)
          t = A(pos(col))
          DO k = pos(col), Apcol(col+1)-1
                                                                       
C           Analyzing the nonzero elements in row "col" of A'.

            jk = Airow(k)
            diag(jk) = diag(jk) + t*A(k)
            above(jk) = above(jk)+1

          ENDDO
          pos(col) = pos(col)+1
        ENDDO

C       Compressing the i-th column of L.

        DO j = Mpcol(i), Mpcol(i+1)-above(i)-1
          Mval(j) = diag(Mirow(j))
          diag(Mirow(j)) = 0.0D0
        ENDDO

      ENDDO

C     Copying the elements above the diagonal of L.

      DO i = 1, m
        pos(i) = Mpcol(i+1)-above(i)
	ENDDO

      DO i = 1, m-1

        DO j = Mpcol(i), Mpcol(i+1)-above(i)-1
          k = Mirow(j)
          Mval(pos(k)) = Mval(j)
          Mirow(pos(k)) = i
          pos(k) = pos(k)+1
        ENDDO

      ENDDO

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE AAtfact(n, Diag, M, Mprow, Micol, Mrnel, Mlen, perm, 
     $                   iperm, Md, prev, next, Mrpos, llfrst, llnext,
     $                   ticol, temp, EpsChol, DenseWin, DenseFrc, 
     $                   DenseRow)

C     This routine computes U^TDU decomposition of M = A*A'. The main
C     diagonal of matrix M must be supplied in vector Diag, while the
C     remaining elements of the matrix must be supplied by rows (or
C     columns) in the triplet {M, Micol, Mprow}. perm is the order of
C     the rows (and columns) of M that reduces the fill-ins during the
C     factorization and must be generated previoulsy by a call to
C     amdbar. iperm is the inverse permutation vector and must also be
C     generated by amdbar. n is the size of the system, while Mlen is 
C     the size of vector M. On exit, Diag contains the diagonal matrix 
C     D, while the elements of the upper triangular factor U are given 
C     by rows in the quartet {M, Mprow, Micol, Mrnel}. Matrix U may
C     have a "dense window" if DenseWin and DenseFrc are chosen 
C     carefully. On input, DenseWin is the dimension of the leading 
C     submatrix that is to be treated as a dense matrix. If DenseWin
C     is greater or equal to n, then a dense factorization is computed.
C     DenseFrc is the ratio between the number of nonzero elements in
C     one row and the size of the matrix that is used to define the
C     dense window. Thus, if the number of nonzero elements in row k
C     is greater than n*DenseFrc, the remaining (n-k+1)x(n-k+1)
C     submatrix of M is treated as a dense matrix. If DenseFrc<=0
C     then a dense factorization is computed. On the other hand, if
C     DenseWin<=1 and DenseFrc>=1, then no dense window is generated.
C     On exit, DenseRow stores the first row of the dense window. 
C     Md(n), prev(n), next(n), temp(n), ticol(n), Mrpos(n), llfrst(n) 
C     and llnext(n) are temporary vectors.

      IMPLICIT NONE

C     Routine parameters.

      INTEGER n, DenseWin, DenseRow, Mlen
      INTEGER Micol(Mlen)
      INTEGER Mprow(n+1)
      INTEGER perm(n), iperm(n), Mrnel(n), prev(n), next(n)
      INTEGER Mrpos(n), llfrst(n), llnext(n), ticol(n)
      REAL*8  EpsChol, DenseFrc
      REAL*8  M(Mlen)
      REAL*8  Diag(n), Md(n)
      LOGICAL temp(n)

C     Local variables.

      INTEGER i, j, k, l, jk, llpos, lltemp, first, last, cont

C     Adjusting DenseRow.

      DenseRow = n - DenseWin + 1
      IF (DenseRow.GT.n) DenseRow = n
      IF (DenseRow.LE.0) DenseRow = 1

C     Returning if there is only one row in A*A'.

      IF (n.EQ.1) THEN
        IF (DABS(Diag(1)).LT.EpsChol) THEN
          WRITE(*,*)'AAt is almost singular.'
          STOP
        ELSE
          RETURN
        ENDIF
      ENDIF

C     Reordering the rows (but not the columns) of M.

      CALL icopy(n+1, Mprow, 1, Mrpos, 1)
      CALL icmpr(Mrpos, Mprow, perm, 1, n)
      DO i = 1, n
        l = perm(i)
        Mrnel(i) = Mrpos(l+1)-Mrpos(l)
        IF (l.GT.1) prev(i) = iperm(l-1)
        IF (l.LT.n) next(i) = iperm(l+1)
      ENDDO
      first = iperm(1)
      last = iperm(n)

C     Setting up some variables.

      DO i = 1, n
        temp(i) = .FALSE.
      ENDDO
      CALL vcmpr(Diag, Md, perm, 1, n)
      CALL vfill(n, 0.0D0, Diag, 1)
      CALL ifill(n, 0, Mrpos, 1)
      CALL ifill(n, 0, llfrst, 1)
      CALL ifill(n, 0, llnext, 1)
      ticol(1) = 0

C     Main loop.

      DO k = 1, n

        temp(k) = .TRUE.

C       Moving the elements from M to Diag.

        Diag(k) = Md(k)
        cont = 2
        DO i = Mprow(k), Mprow(k)+Mrnel(k)-1
          jk = iperm(Micol(i))
          IF (jk.GT.k) THEN
            Diag(jk) = M(i)
            temp(jk) = .TRUE.
            ticol(cont) = jk
            cont = cont+1
          ENDIF
        ENDDO

        llpos = llfrst(k)
        DO WHILE (llpos.GT.0) 

C         Performing the Cholesky decomposition (sparse part).

          jk = Mprow(llpos)+Mrpos(llpos)
          CALL axpys(-M(jk)*Diag(llpos), M, Diag, Micol,
     $                  jk, Mprow(llpos)+Mrnel(llpos)-1, 0.0D0)
          Mrpos(llpos) = Mrpos(llpos)+1
          IF (k.LT.DenseRow) THEN
            DO j = jk, Mprow(llpos)+Mrnel(llpos)-1
              IF (.NOT.temp(Micol(j))) THEN
                temp(Micol(j)) = .TRUE.
                ticol(cont) = Micol(j)
                cont = cont+1
              ENDIF
            ENDDO
          ENDIF

C         Updating the liked list that stores the columns of U.

          lltemp = llnext(llpos)
          IF (Mrpos(llpos).LT.Mrnel(llpos)) THEN
            l = Micol(Mprow(llpos)+Mrpos(llpos))
            IF (llfrst(l).EQ.0) THEN
              llnext(llpos) = 0
              llfrst(l) = llpos
            ELSE
              llnext(llpos) = llfrst(l)
              llfrst(l) = llpos
            ENDIF
          ENDIF
          llpos = lltemp

        ENDDO

C       Checking if the row is sparse or dense.

        IF (k.LT.DenseRow) THEN
          IF (DBLE(cont-2).GE.(DenseFrc*(n-k))) THEN

C           Changing to a dense factoriz. if current row is not sparse.

            DenseRow = k

          ELSE

C           Maintaining the sparse factorization.

            DO i = 2, cont-1
              temp(ticol(i)) = .FALSE.
            ENDDO

C           Sorting the row elements according to the "Insertion Sort".

            DO i = 3, cont-1
              l = ticol(i)
              j = i
              DO WHILE (ticol(j-1).GT.l)
                ticol(j) = ticol(j-1)
                j = j-1
              ENDDO
              ticol(j) = l
            ENDDO

          ENDIF
        ENDIF

        IF (k.GE.DenseRow) THEN

C         Performing the Cholesky decomposition (dense part).

          DO j = DenseRow, k-1
            jk = Mprow(j)-1+k-j
            CALL daxpy(n-k+1, -M(jk)*Diag(j), M(jk), 1, Diag(k), 1)
          ENDDO

C         Moving the row to vector M (dense part).

          IF ((n-k).GT.(Mprow(next(k))-Mprow(k))) THEN
            CALL MoveRow(k, n-k, n, Mlen, M, Mprow, Micol, Mrnel,
     $                   prev, next, first, last, DenseRow)
          ENDIF 

          IF (DABS(Diag(k)).LT.EpsChol) THEN
            WRITE(*,*)'AAt is almost singular.'
            STOP
          ENDIF
          CALL dscal(n-k, 1.0D0/Diag(k), Diag(k+1), 1)
          CALL dcopy(n-k, Diag(k+1), 1, M(Mprow(k)), 1)
          Mrnel(k) = n-k
          CALL vfill(n-k, 0.0D0, Diag(k+1), 1)

        ELSE

C         Moving the row to vector U (sparse part).

          IF ((cont-2).GT.(Mprow(next(k))-Mprow(k))) THEN
            CALL MoveRow(k, cont-1, n, Mlen, M, Mprow, Micol, Mrnel,
     $                   prev, next, first, last, DenseRow)
          ENDIF 

          CALL vcmpr(Diag, M(Mprow(k)-1), ticol, 2, cont-1)
          IF (DABS(Diag(k)).LT.EpsChol) THEN
            WRITE(*,*)'AAt is almost singular.'
            STOP
          ENDIF
          CALL dscal(cont-2, 1.0D0/Diag(k), M(Mprow(k)), 1)
          Mrnel(k) = cont-2
          CALL vfills(0.0D0, Diag, ticol, 2, cont-1)
          CALL icopy(cont-2, ticol(2), 1, Micol(Mprow(k)), 1)

C         Updating the liked list that stores the columns of U.

          IF (Mrnel(k).GT.0) THEN
            l = Micol(Mprow(k))
            IF (llfrst(l).EQ.0) THEN
              llnext(k) = 0
              llfrst(l) = k
            ELSE
              llnext(k) = llfrst(l)
              llfrst(l) = k
            ENDIF
          ENDIF

        ENDIF

      ENDDO

C     CALL TriMatPrint(n, M, Micol, Mprow, Mrnel)
C     CALL VecPrint(n, Diag)
C     CALL IVecPrint(n, perm)
C     STOP

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE MoveRow(k, rsize, n, Mlen, M, Mprow, Micol, Mrnel,
     $                   prev, next, first, last, DenseRow)

C     This auxiliary routine moves a row of M to the end of the
C     vector and compresses M and Micol if necessary.

      IMPLICIT NONE

C     Routine parameters.

      INTEGER k, rsize, n, Mlen, first, last, DenseRow
      INTEGER Micol(Mlen)
      INTEGER Mprow(n+1)
      INTEGER Mrnel(n), prev(n), next(n)
      REAL*8  M(Mlen)

C     Local variables.

      INTEGER j, l, diff

      IF (k.NE.last) THEN

C       Moving the k-th row to the end of vectors M and Micol.

        IF (k.NE.first) next(prev(k)) = next(k)
        prev(next(k)) = prev(k)
        next(last) = k
        prev(k) = last
        next(k) = n+1
        last = k
        Mprow(k) = Mprow(prev(k)) + Mrnel(prev(k))

      ENDIF

      IF ((Mprow(k)+rsize-1).GT.Mlen) THEN

C       Compressing M and Micol. Updating Mprow.

        j = next(first)
        diff = 0
        DO WHILE (j.NE.last)
          l = Mprow(j)
          diff = diff + l - Mprow(prev(j)) - Mrnel(prev(j))
          CALL dcopy(Mrnel(j), M(l), 1, M(l-diff), 1)
          IF (j.LT.DenseRow) THEN
            CALL icopy(Mrnel(j), Micol(l), 1, Micol(l-diff), 1)
          ENDIF
          Mprow(j) = Mprow(j)-diff
          j = next(j)                
        ENDDO
        Mprow(last) = Mprow(prev(last))+Mrnel(prev(last))

      ENDIF

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE MsolveCh(n, Diag, M, Mprow, Micol, Mrnel,
     $                    perm, iperm, DenseRow, tmp, v)

C     This routine solves the system (A*A')v = b, where the U^TDU
C     decomposition of (A*A') is supposed to be generated by a
C     previous call to AAtfact. Diag must contain the diagonal 
C     matrix D, while the elements of the upper triangular factor U 
C     must be given by rows in the quartet {M, Mprow, Micol, Mrnel}. 
C     n is the dimension of the system. perm is the order of the 
C     rows (and columns) of M that reduces the fill-ins during the 
C     factorization. iperm is the inverse permutation vector. 
C     On input, v must contain the RHS vector b. On exit, v stores
C     the solution of the system. tmp(n) is a temporary vector.
C     Matrix M may have a "dense window" if the parameter DenseRow
C     is set (by AAtfact) to a value less than n. In this case, the
C     leading (m-DenseRow+1)x(m-DenseRow+1) submatrix of U is 
C     treated as a dense matrix.

      IMPLICIT NONE

C     Routine parameters.

      INTEGER n, DenseRow
      INTEGER Micol(*)
      INTEGER Mprow(n+1)
      INTEGER Mrnel(n), perm(n), iperm(n)
      REAL*8  Diag(n), tmp(n), v(n)
      REAL*8  M(*)

C     Local variables.

      INTEGER k

C     Functions called by the routine.

      REAL*8  ddot, dots

C     Moving b to tmp and changing the order of the elements.

      CALL vcmpr(v, tmp, perm, 1, n)

C     Solving Ly = b(perm) (sparse part).

      DO k = 1, DenseRow-1

        CALL axpys(-tmp(k), M, tmp, Micol, Mprow(k),
     $             Mprow(k)+Mrnel(k)-1, 0.0D0)

      ENDDO 

C     Solving Ly = b(perm) (dense part).

C     beg = Mprow(DenseRow)
      DO k = DenseRow, n-1
        CALL daxpy(n-k, -tmp(k), M(Mprow(k)), 1, tmp(k+1), 1)
C       beg = beg + n - k
      ENDDO

C     Solving Dz = y.

      DO k = 1, n  
        tmp(k) = tmp(k)/Diag(k)
      ENDDO 

C     Solving L^Tx = z (dense part).

C     beg = Mprow(DenseRow) - 1 + ((n-DenseRow+1)*(n-DenseRow))/2
      DO k = n-1, DenseRow, -1
        tmp(k) = tmp(k) - ddot(n-k, M(Mprow(k)), 1, tmp(k+1), 1)
C       beg = beg - n + k - 1
      ENDDO

C     Solving L^Tx = z (sparse part).

      DO k = DenseRow-1, 1, -1

        tmp(k) = tmp(k) - dots(M, tmp, Micol, Mprow(k), 
     $                         Mprow(k)+Mrnel(k)-1)

      ENDDO

C     Reordering the elements of x.

      CALL vcmpr(tmp, v, iperm, 1, n)

      RETURN
      END

C     ------------------------------------------------------------------

      REAL*8 FUNCTION ResNorm(m, n, A, Airow, Apcol, x, b)

C     This routine computes the norm of residual r = A*A'*x - b,
C     to check if the solution obtained by MsolveCh is correct.

      IMPLICIT NONE

      INTEGER nmax
                     
      PARAMETER (nmax=60000)

C     Routine parameters.

      INTEGER m, n
      INTEGER Airow(*)
      INTEGER Apcol(n+1)
      REAL*8  x(m), b(m)
      REAL*8  A(*)

C     Local variable.

      INTEGER i
      REAL*8  tmp(nmax), r(nmax)

C     Function called by the routine.

      REAL*8  dnrm2

C     Computing r = A*A'*x - b.

      CALL MtVprod(n, A, x, tmp, Airow, Apcol)
      CALL MVprod(m, n, A, tmp, r, Airow, Apcol, 0.0D0)
      CALL daxpy(m, -1.0D0, b, 1, r, 1)

c$$$      DO i=1,m
c$$$        WRITE(*,*)'res(',i,') = ',r(i)
c$$$      ENDDO
c$$$      WRITE(*,*)
    
      ResNorm = dnrm2(m, r, 1)

      RETURN
      END
