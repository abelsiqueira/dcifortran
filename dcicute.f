!     ------------------------------------------------------------------
!     PROGRAM : DCICute
!     VERSION : 2.0.
!     DATE    : June, 2000.
!     CONTENTS: subroutine SetUp
!               function   Fun
!               subroutine Constr
!               subroutine Grad
!               subroutine GenAg
!               subroutine Jacob
!               subroutine Hprod
!
!     External routines called by this module:
!
!     From dciblas : vcmpr, icmpr, vexpn, ifill, vfill.
!     From the blas: dcopy.
!     From the cute: csetup, ufn, cofg, ccfg, ugr, csgr, cprod, uprod. 
!
!     ------------------------------------------------------------------

      SUBROUTINE CuterSetUp(Prob, m, n, lincon, x, lambda)

      ! Function that determines the number of variables and
      ! constraints and the number of fixed variables. It also
      ! reads initial approximations for x and lambda and
      ! checks if all of the constraints are linear (in this
      ! case, lincon is set to .TRUE.)

      IMPLICIT NONE
      
      INTEGER nmax, mmax, amax
                     
      PARAMETER (nmax=60000, mmax=40000, amax=1000000)

C     Subroutine parameters.

      CHARACTER*10 Prob
      INTEGER      m, n
      REAL*8       lambda(mmax)
      REAL*8       x(nmax)
      LOGICAL      lincon

C     Variables passed by COMMON statements.

      INTEGER      mt, nfix, nt
      INTEGER      ind(nmax), rind(nmax)
      REAL*8       xfix(nmax), vfix(nmax)
      LOGICAL      FreeV(nmax)
      LOGICAL      linconstr

C     Internal variables.

      INTEGER      i, nzj, err, ifile
      LOGICAL      tf
      LOGICAL      equatn(mmax), linear(mmax)
      REAL*8       bl(nmax), bu(nmax)
      REAL*8       cl(mmax), cu(mmax)
      CHARACTER*10 varn(nmax)
      CHARACTER*10 constrn(mmax)

C     COMMON statements.

      COMMON / dcifixvar    / nfix, ind
      COMMON / dcirevind    / rind
      COMMON / dcifixsol    / xfix
      COMMON / dcitotvar    / nt
      COMMON / dcifixtype   / FreeV
      COMMON / dcifixvect   / vfix
      COMMON / dcitotconstr / mt
      COMMON / dcilinconstr / linconstr

C     Opening OUTSDIF.d.

      ifile = 30
      OPEN(ifile, FILE='OUTSDIF.d', FORM='FORMATTED', 
     $     STATUS='OLD',IOSTAT=err)
      REWIND ifile
      IF (err.NE.0) THEN
        WRITE(*,*)'Could not open the OUTSDIF.d file.'
        STOP
      ENDIF
    
C     Obtaining the number of variables and constraints and
C     verifying if the vectors are correctly dimensioned.

      CALL cdimen(ifile, nt, m)

      IF (nt.GT.nmax) THEN
        WRITE(*,*)'gmm99cute error: Increase nmax to, at least, ',nt
        STOP
      ENDIF

      IF (m.GT.mmax) THEN
        WRITE(*,*)'gmm99cute error: Increase mmax to, at least, ',m
        STOP
      ENDIF

C     Calling the setup routine from CUTE.

      tf = .FALSE.
      CALL csetup(ifile, 6, nt, m, x, bl, bu, nmax, equatn, 
     $            linear, lambda, cl, cu, mmax, tf, tf, tf)

      CALL cdimsj(nzj)

C     Checking if amax is big enough.

      IF (nzj.GT.amax) THEN
        WRITE(*,*)'gmm99cute error: Increase amax to, at least, ',nzj
        STOP
      ENDIF

C     Reading the problem name.

      CALL cnames(nt, m, prob, varn, constrn)

C     Checking if all constraints are equations (not inequalities). 

      DO i = 1,m
        IF (.NOT.equatn(i)) THEN
          WRITE(*,*)'There are inequality constraints.'
          STOP
        ENDIF
      ENDDO   

C     Checking if all of the constraints are linear. 

      i = 1
      DO WHILE ((i.LE.m).AND.(linear(i)))
        i = i+1
      ENDDO
      IF (i.LE.m) THEN
        lincon = .FALSE.
        linconstr = .FALSE.
      ELSE
C       This option is temporarily disabled.
C       lincon = .TRUE.
C       linconstr = .TRUE.
        lincon = .FALSE.
        linconstr = .FALSE.
      ENDIF   

C     Copying m to mt, that is passed to fun and grad.

      mt = m 
      
C     Checking if all variables are free (i.e., if there are no
C     bounded variables). Fixed variables are allowed.

      nfix = 0
      n = 0
      DO i = 1,nt
        IF (DABS(bu(i)-bl(i)).LT.1.0D-15) THEN
          xfix(i) = bl(i)
          vfix(i) = 0
          nfix = nfix+1
          FreeV(i) = .FALSE.
        ELSE 
          IF ((bl(i).GT.-1.0D-19).OR.(bu(i).LT.1.0D19)) THEN
            WRITE(*,*)'Ignoring bounds of the ',varn(i),' variable.'
          ENDIF
          n = n+1
          ind(n) = i
C         rind(i) = n
          FreeV(i) = .TRUE.
        ENDIF
      ENDDO

      IF (nfix.NE.0) THEN
        WRITE(*,*)'This problem has ',nfix,' fixed variables.'
        CALL vcmpr(x, x, ind, 1, n)
      ENDIF

      RETURN
      END
      
C     ------------------------------------------------------------------
      
      REAL*8 FUNCTION CuterFun(n, x)

      ! Function that computes f(x), the value of the function
      ! being minimized.

      IMPLICIT NONE
      
      INTEGER nmax

      PARAMETER (nmax=60000)

C     Subroutine parameters.

      INTEGER n
      REAL*8  x(n)

C     Variables passed by COMMON statements.

      INTEGER nt, nfix, m
      INTEGER ind(nmax)
      REAL*8  xfix(nmax)

C     Internal variables.

      REAL*8  void, f
      LOGICAL tf

C     COMMON statements.

      COMMON / dcitotvar    / nt
      COMMON / dcifixvar    / nfix, ind
      COMMON / dcifixsol    / xfix
      COMMON / dcitotconstr / m

      IF (m.EQ.0) THEN

        IF (nfix.EQ.0) THEN
          CALL ufn(nt, x, f)
        ELSE
          CALL vexpn(x, xfix, ind, 1, n)
          CALL ufn(nt, xfix, f)
        ENDIF

      ELSE

        tf = .FALSE.
        IF (nfix.EQ.0) THEN
          CALL cofg(nt, x, f, void, tf)
        ELSE
          CALL vexpn(x, xfix, ind, 1, n)
          CALL cofg(nt, xfix, f, void, tf)
        ENDIF
 
      ENDIF

      CuterFun=f

      RETURN
      END

C     ------------------------------------------------------------------
      
      SUBROUTINE CuterConstr(m, n, x, h) 

      ! Routine that computes h(x), the constraint values at x.

      IMPLICIT NONE

      INTEGER nmax

      PARAMETER (nmax=60000)

C     Subroutine parameters.

      INTEGER m, n
      REAL*8  x(n)
      REAL*8  h(m)

C     Variables passed by COMMON statements.

      INTEGER nt, nfix
      INTEGER ind(nmax)
      REAL*8  xfix(nmax)

C     Internal variables.

      REAL*8  void
      LOGICAL tf

C     COMMON statements.

      COMMON / dcitotvar / nt
      COMMON / dcifixvar / nfix, ind
      COMMON / dcifixsol / xfix

      tf = .FALSE.
      IF (nfix.EQ.0) THEN
        CALL ccfg(nt, m, x, m, h, tf, m, nt, void, tf)
      ELSE
        CALL vexpn(x, xfix, ind, 1, n)
        CALL ccfg(nt, m, xfix, m, h, tf, m, nt, void, tf)
      ENDIF 

      RETURN
      END

C     ------------------------------------------------------------------
      
      SUBROUTINE CuterGrad(n, x, g) 

      ! Routine that computes g(x), the gradient of f(x).

      IMPLICIT NONE
      
      INTEGER nmax

      PARAMETER (nmax=60000)

C     Subroutine parameters.

      INTEGER n
      REAL*8  x(n), g(n)

C     Variables passed by COMMON statements.

      INTEGER m, nt, nfix
      INTEGER ind(nmax)
      REAL*8  gv(nmax), xfix(nmax)
      LOGICAL gavail, linconstr

C     LOCAL CONVERGENCE
      LOGICAL forcegrad

C     Internal variables.

      REAL*8  f
      LOGICAL tf

C     COMMON statements.

      COMMON / dcigradavail / gavail
      COMMON / dcigradient  / gv
      COMMON / dcifixvar    / nfix, ind
      COMMON / dcifixsol    / xfix
      COMMON / dcitotvar    / nt
      COMMON / dcitotconstr / m
      COMMON / dcilinconstr / linconstr

C     LOCAL CONVERGENCE
      COMMON / dciforcegrad / forcegrad

      IF (m.EQ.0) THEN

        IF (nfix.EQ.0) THEN
          CALL ugr(nt, x, g)
        ELSE
          CALL vexpn(x, xfix, ind, 1, n)
          CALL ugr(nt, xfix, gv)
          CALL vcmpr(gv, g, ind, 1, n)
        ENDIF

      ELSE IF (gavail) THEN

        CALL dcopy(n, gv, 1, g, 1)

      ELSE IF (linconstr) THEN

        tf = .TRUE.
        IF (nfix.EQ.0) THEN
          CALL cofg(nt, x, f, g, tf)
        ELSE
          CALL vexpn(x, xfix, ind, 1, n)
          CALL cofg(nt, xfix, f, gv, tf)
          CALL vcmpr(gv, g, ind, 1, n)
        ENDIF
 
      ELSE IF (forcegrad) THEN

        IF (nfix.EQ.0) THEN
          CALL ugr(nt, x, g)
        ELSE
          CALL vexpn(x, xfix, ind, 1, n)
          CALL ugr(nt, xfix, gv)
          CALL vcmpr(gv, g, ind, 1, n)
        ENDIF

      ELSE

        WRITE(*,*)' Cannot evaluate the gradient of f(x).'
        STOP

      ENDIF

      gavail = .FALSE.

      RETURN
      END
      
C     ------------------------------------------------------------------

      SUBROUTINE GenAg(nt, n, nfix, ind, rind, FreeV, g, J, Jirow, 
     $                 Jicol, nJ, A, Airow, Apcol, nA)
      
C     Generates A, the Jacobian of the constraints, and g, the
C     gradient of f(x), removing fixed variables.
     
      IMPLICIT NONE
      
C     Subroutine parameters.

      INTEGER nt, n, nfix, nJ, nA
      INTEGER ind(n), rind(n)
      INTEGER Jicol(*), Jirow(*), Airow(*)
      INTEGER Apcol(n+1)
      LOGICAL FreeV(n)
      REAL*8  J(*), A(*)
      REAL*8  g(n)

C     Internal variables.
 
      INTEGER i, k, ti, jk, l
      REAL*8  lA

      nA = nJ
      CALL ifill(nt+1, 0, Apcol, 1)
      CALL vfill(n, 0.0D0, g, 1)
      i = 1
      
      IF (nfix.EQ.0) THEN
      
C       Storing the number of nonzero elements on each column in Apcol
C       and creating g when there is no fixed variable.

        DO WHILE (i.LE.nA)

          IF (Jirow(i).EQ.0) THEN

C           Extracting the objective function from J.

            g(Jicol(i)) = J(i)
            J(i) = J(nA)
            Jicol(i) = Jicol(nA)
            Jirow(i) = Jirow(nA)
            nA = nA-1

          ELSE
          
C           A constraint was found.

            Apcol(Jicol(i)) = Apcol(Jicol(i))+1
            i = i+1

          ENDIF  
        ENDDO

      ELSE
      
C       Storing the number of nonzero elements on each column in Apcol,
C       creating g and removing the columns of fixed variables.

        DO WHILE (i.LE.nA)
          IF (Jirow(i).EQ.0) THEN

C           Extracting the objective function from J.

C           g(rind(Jicol(i))) = J(i)
            g(Jicol(i)) = J(i)
            J(i) = J(nA)
            Jicol(i) = Jicol(nA)
            Jirow(i) = Jirow(nA)
            nA = nA-1  

          ELSE IF (FreeV(Jicol(i))) THEN

C           The variable is not fixed.

C           Apcol(rind(Jicol(i))) = Apcol(rind(Jicol(i)))+1
            Apcol(Jicol(i)) = Apcol(Jicol(i))+1
            i = i+1

          ELSE

C           The variable is fixed.

            J(i) = J(nA)
            Jicol(i) = Jicol(nA)
            Jirow(i) = Jirow(nA)
            nA = nA-1  

          ENDIF  
        ENDDO
      ENDIF

C     Accumulating the values of Apcol.

      Apcol(1) = Apcol(1)+1
      DO i = 2,nt
        Apcol(i) = Apcol(i-1) + Apcol(i)
      ENDDO

C     Copying data from J and Jirow to A and Airow, 
C     and adjusting Apcol.

      DO i = nA, 1, -1
C       k = Apcol(rind(Jicol(i)))-1
        k = Apcol(Jicol(i))-1
C       Apcol(rind(Jicol(i))) = k
        Apcol(Jicol(i)) = k
        A(k) = J(i)
        Airow(k) = Jirow(i)
      ENDDO 
      Apcol(nt+1) = nA + 1

C     Sorting the column elements of A.

      DO k = 1, nt
        ti = Airow(Apcol(k+1))
        Airow(Apcol(k+1)) = nA + 2
        DO i = Apcol(k+1)-1, Apcol(k), -1
          l = Airow(i)
          lA = A(i)
          jk = i
          DO WHILE (Airow(jk+1).LT.l)
            Airow(jk) = Airow(jk+1)
            A(jk) = A(jk+1)
            jk = jk+1
          ENDDO
          Airow(jk) = l
          A(jk) = lA
        ENDDO
        Airow(Apcol(k+1)) = ti
      ENDDO

C     Compressing g and Apcol if there are fixed variables.

      IF (nfix.NE.0) THEN
        CALL vcmpr(g, g, ind, 1, n)
        CALL icmpr(Apcol, Apcol, ind, 1, n)
      ENDIF 
      Apcol(n+1) = nA + 1
      
      RETURN
      END

C     ------------------------------------------------------------------
      
      SUBROUTINE CuterJacob(m, n, x, A, Airow, Apcol) 
      
C     Computes A, the Jacobian of the constraints. Also computes g, 
C     the gradient of f(x). Warning: A is mxn and not nxm,as usual. 
C     m is the number of constraints. n is the number of variables.
C     Vectors A, Airow and Apcol are used to store the Jacobian of 
C     the constraints. nA is the number of nonzero elements in A.
C     Vector A stores the nonzero elements by columns. Airow 
C     contains the row index of each element in vector A. 
C     Apcol(i) indicates the position on vector A of the first 
C     nonzero element in column i of the Jacobian. Apcol(n+1)=nA.
C     The gradient of f(x) is passed by COMMON to function Grad.
C     Prcavail is a variable that indicates if a preconditioner for 
C     matrix (A*A') is available. Each time A is recomputed, prcavail
C     must be set to FALSE. This is also done here.

      IMPLICIT NONE

      INTEGER mmax, nmax, amax

      PARAMETER (nmax=60000, mmax=40000, amax=1000000)
    
C     Subroutine parameters.

      INTEGER m, n
      INTEGER Airow(amax)
      INTEGER Apcol(n+1)
      REAL*8  x(n)
      REAL*8  A(amax)

C     Variables passed by COMMON statements.

      INTEGER nfix, nt
      INTEGER ind(nmax), rind(nmax)
      REAL*8  xfix(nmax), g(nmax)  
      LOGICAL gavail, prcavail
      LOGICAL FreeV(nmax)

C     Internal variables.

      INTEGER nJ, nA
      INTEGER Jirow(amax), Jicol(amax)
      REAL*8  void
      REAL*8  J(amax)
      LOGICAL tf

C     COMMON statements.

      COMMON / dcigradavail / gavail
      COMMON / dcigradient  / g
      COMMON / dcifixvar    / nfix, ind
C     COMMON / dcirevind    / rind
      COMMON / dcifixsol    / xfix
      COMMON / dcitotvar    / nt
      COMMON / dcifixtype   / FreeV
      COMMON / dciprecavail / prcavail

      tf = .FALSE.
      IF (nfix.EQ.0) THEN
        CALL csgr(nt, m, tf, m, void, x, nJ, amax, J, Jicol, Jirow)
      ELSE
        CALL vexpn(x, xfix, ind, 1, n)
        CALL csgr(nt, m, tf, m, void, xfix, nJ, amax, J, Jicol, Jirow)
      ENDIF
      
C     Generating A and g.

      CALL GenAg(nt, n, nfix, ind, rind, FreeV, g, J, Jirow,
     $           Jicol, nJ, A, Airow, Apcol, nA)

      gavail = .TRUE.
      prcavail = .FALSE.

      IF (nA.GT.amax) THEN
        WRITE(*,*)'DCI failed to store A. Increase amax.'
        STOP
      ENDIF

      RETURN
      END
      
C     ------------------------------------------------------------------

      SUBROUTINE CuterHprod(n, m, x, Lambda, v, Hv, GotH)

C     Computes the product Hv = H*v at (x, Lambda), where H is the 
C     Hessian of the Lagrangian. 
C     m is the number of constraints. n is the number of variables. 
C     GotH is a logical variable that indicates if a call to Hprod
C     has already been made at the current point (x, Lambda).
C     GotH must be set to .FALSE. whenever the point is changed.

      IMPLICIT NONE

      INTEGER nmax

      PARAMETER (nmax=60000)     

C     Subroutine parameters.

      INTEGER m, n
      REAL*8  x(n), v(n), Hv(n)
      REAL*8  Lambda(m) 
      LOGICAL GotH

C     Variables passed by COMMON statements.

      INTEGER nfix, nt
      INTEGER ind(nmax)
      REAL*8  xfix(nmax), vfix(nmax)                        

C     COMMON statements.

      COMMON / dcifixvar   / nfix, ind
      COMMON / dcifixsol   / xfix
      COMMON / dcifixvect  / vfix
      COMMON / dcitotvar   / nt

      IF (nfix.EQ.0) THEN

        IF (m.GT.0) THEN    
          CALL cprod(nt, m, GotH, x, m, Lambda, v, Hv)
        ELSE
          CALL uprod(nt, GotH, x, v, Hv)
        ENDIF

      ELSE

        CALL vexpn(x, xfix, ind, 1, n)
        CALL vexpn(v, vfix, ind, 1, n)
        IF (m.GT.0) THEN    
          CALL cprod(nt, m, GotH, xfix, m, Lambda, vfix, Hv)
        ELSE
          CALL uprod(nt, GotH, xfix, vfix, Hv)
        ENDIF
        CALL vcmpr(Hv, Hv, ind, 1, n)

      ENDIF

      GotH = .TRUE.

      RETURN
      END


