!     ------------------------------------------------------------------
!     PROGRAM : DCI.
!     MODULE  : DCIBlas.
!     VERSION : 2.0.
!     DATE    : June, 2000.
!     CONTENTS: Basic linear algebra subroutines.
!     ------------------------------------------------------------------

      SUBROUTINE vfill(n, alpha, v, incv)

!     Input:
!       v     - Real vector.
!       alpha - Real number.
!       n     - Dimension of v.
!       incv  - Increment.
!     Output:
!       v     - Real vector.
!     Description:
!       Fills v(1:incv:n) with alpha.

      IMPLICIT NONE

      INTEGER n, incv
      INTEGER i, nincv
      REAL*8  v(*)
      REAL*8  alpha

      IF ((n.LE.0).OR.(incv.LE.0)) RETURN

      nincv = n*incv
      DO i = 1,nincv,incv
        v(i) = alpha
      ENDDO

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE vfills(alpha, v, ind, ibeg, iend)

!     Input:
!       v          - Real vector.
!       alpha      - Real number.
!       ind        - index of v elements that are to be used.
!       ibeg, iend - first and last elements of ind that are used.
!     Output:
!       v          - Real vector.
!     Description:
!       Fills v(ind) with alpha.

      INTEGER ibeg, iend
      INTEGER i
      INTEGER ind(*)
      REAL*8  alpha
      REAL*8  v(*)

      DO i = ibeg, iend
        v(ind(i)) = alpha
      ENDDO

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE ifill(n, a, v, incv)

!     Input:
!       v     - Integer vector.
!       a     - Integer number.
!       n     - Dimension of v.
!       incv  - Increment.
!     Output:
!       v     - Integer vector.
!     Description:
!       Fills v(1:incv:n) with a.

      IMPLICIT NONE

      INTEGER n, a, incv
      INTEGER i, nincv
      INTEGER v(*)

      IF ((n.LE.0).OR.(incv.LE.0)) RETURN

      nincv = n*incv
      DO i = 1,nincv,incv
        v(i) = a
      ENDDO

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE ifills(a, v, ind, ibeg, iend)

!     Input:
!       v          - Integer vector.
!       a          - Integer number.
!       ind        - index of v elements that are to be used.
!       ibeg, iend - first and last elements of ind that are used.
!     Output:
!       v          - Integer vector.
!     Description:
!       Fills v(ind) with a.

      INTEGER ibeg, iend, a
      INTEGER i
      INTEGER v(*), ind(*)

      DO i = ibeg, iend
        v(ind(i)) = a
      ENDDO

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE icopy(n, v, incv, w, incw)

!     Input:
!       v     - Integer vector.
!       n     - Dimension of v.
!       incv  - Increment in v.
!       incv  - Increment in w.
!     Output:
!       w     - Integer vector.
!     Description:
!       Copies v onto w.

      IMPLICIT NONE

      INTEGER n, incv, incw
      INTEGER i, nincv
      INTEGER v(*), w(*)

      IF ((n.LE.0).OR.(incv.LE.0)) RETURN

      nincv = n*incv
      DO i = 1,nincv,incv
        w(i) = v(i)
      ENDDO

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE axpys(alpha, x, y, ind, ibeg, iend, zero)

!     Input:
!       x, y       - Real vectors.
!       ind        - index of y elements that are to be used.
!       ibeg, iend - first and last elements of ind and x that are used.
!       alpha      - Real number.
!       zero       - Small nonnegative number. If ABS(alpha) < zero, 
!                    nothing is done.
!     Output:
!       y          - Real vector.
!     Description:
!       y(ind) = y(ind) + alpha*x.

      INTEGER ibeg, iend
      INTEGER i
      INTEGER ind(*)
      REAL*8  alpha, zero
      REAL*8  x(*), y(*)

      IF (DABS(alpha).GT.zero) THEN
        DO i = ibeg, iend
          y(ind(i)) = y(ind(i)) + alpha*x(i)
        ENDDO
      ENDIF

      RETURN
      END

C     ------------------------------------------------------------------

      REAL*8 FUNCTION dots(v, w, ind, ibeg, iend)

!     Input:
!       v, w       - Real vectors.
!       ind        - index of y elements that are to be used.
!       ibeg, iend - first and last elements of ind and w that are used.
!     Output:
!       dots       - Real number.
!     Description:
!       dots = v'*w(ind).

      INTEGER ibeg, iend
      INTEGER i
      INTEGER ind(*)
      REAL*8  t
      REAL*8  v(*), w(*)

      t = 0.0D0
      DO i = ibeg, iend
        t = t + v(i)*w(ind(i))
      ENDDO

      dots = t

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE vcmpr(v, vcomp, ind, ibeg, iend)

!     Input:
!       v          - Real vector.
!       ind        - index of v elements that are to be used.
!       ibeg, iend - first and last elements of ind that are used.
!     Output:
!       Vcomp      - Real vector.
!     Description:
!       Creates vcomp = v[ind].

      INTEGER ibeg,iend,i
      INTEGER ind(*)
      REAL*8  v(*), vcomp(*)

      DO i=ibeg,iend
        vcomp(i) = v(ind(i))
      ENDDO

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE icmpr(v, vcomp, ind, ibeg, iend)

!     Input:
!       v          - Integer vector.
!       ind        - index of v elements that are to be used.
!       ibeg, iend - first and last elements of ind that are used.
!     Output:
!       Vcomp      - Integer vector.
!     Description:
!       Creates vcomp = v[ind].

      INTEGER ibeg,iend,i
      INTEGER ind(*), v(*), vcomp(*)

      DO i=ibeg,iend
        vcomp(i) = v(ind(i))
      ENDDO

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE vexpn(v, vexp, ind, ibeg, iend)

!     Input:
!       v          - Real vector.
!       ind        - index of v elements that are to be used.
!       ibeg, iend - first and last elements of v and ind that are used.
!     Output:
!       Vexp      - Real vector.
!     Description:
!       Creates vexp[ind] = v.

      INTEGER i, ibeg, iEnd
      INTEGER ind(*)
      REAL*8  v(*), vexp(*)

      DO i=ibeg,iend
        vexp(ind(i)) = v(i)
      ENDDO

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE MVprod(m, n, Mat, v, Mv, Mirow, Mpcol, zero)

!     Input:
!       m, n  - Number of rows and columns of Mat.
!       Mat   - Real vector of the matrix nonzero elements.
!       Mirow - Integer vector. Row indices of the elements in M.
!       MPcol - Integer vector. Position in M of the first nonzero
!               element of each column of the matrix.
!       v     - Real vector.
!       zero  - Small nonnegative number. If ABS(v(i)) < zero, this
!               element is not multiplied by the i-th column of M.
!     Output:
!       Mv    - Real vector.
!     Description:
!       Computes Mv = M*v.

      INTEGER m, n
      INTEGER i
      INTEGER Mirow(*), Mpcol(*)
      REAL*8  Mat(*), v(*), Mv(*)
      REAL*8  zero

      CALL vfill(m, 0.0D0, Mv, 1)
      DO i = 1, n
 	  CALL axpys(v(i), Mat, Mv, Mirow, Mpcol(i), Mpcol(i+1)-1, zero)
      ENDDO

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE MtVprod(n, Mat, v, Mtv, Mirow, Mpcol)

!     Input:
!       n     - Number of columns of Mat.
!       Mat   - Real vector containing the matrix nonzero elements.
!       Mirow - Integer vector. Row indices of the elements in M.
!       MPcol - Integer vector. Position in M of the first nonzero
!               element of each column of the matrix.
!       v     - Real vector.
!     Output:
!       Mtv   - Real vector.
!     Description:
!       Computes Mtv = M'*v.

      INTEGER n
      INTEGER i
      INTEGER Mirow(*), Mpcol(*)
      REAL*8  Mat(*), v(*), Mtv(*)
      REAL*8  dots

      DO i = 1, n
        Mtv(i) = dots(Mat, v, Mirow, Mpcol(i), Mpcol(i+1)-1)
      ENDDO

      RETURN
      END

C     ------------------------------------------------------------------
