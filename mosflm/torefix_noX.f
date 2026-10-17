C $Id: torefix.for.noX11,v 1.1 2012/03/23 15:09:50 harry Exp $ dummy torefix.for routine << DON'T EDIT THIS LINE!
C     ^^^^^^^^^^^^^^^^^^^^^^^^^ this string is essential for the
C     distribution build to work correctly. Please don't edit it in any
C     way unless you also change the file "add_refix.com".
C
      SUBROUTINE TOREFIX(NSOL,BOXOPEN,RFIXCELL,RFIXDIST)
C     ==================================================
C
      IMPLICIT NONE
C&&*&& include  ../inc/ioo.f
C
C $Id: ioo.f,v 1.9 2020/12/15 20:08:54 andrew Exp $
C
C--- awk generated include file  ioo.h
C---- START of include file ioo.h
C
C
C
C     .. Scalars in common block /IOO/ ..
      INTEGER IOUT,IUNIT,ITIN,ITOUT,INOD,INMO,IDU,NWRN,IBRIEF,ISUMMR,
     +        ICOORD,SERVERFD,SERVERFDIN,NSHUTERR,dnafd
      LOGICAL ONLINE,ONEFILE,FHEADER,BRIEF,GRAPH,IOERR,
     $     NODISPLAY,LBELL,JPGOUT,SOCKLO,LBEST,dnaout,INMOOPEN,
     $     XMLOPEN,SOCKLOPHIPROF
      CHARACTER dna_image*80,xmlfile*80
C     ..
C     .. Common block /IOO/ ..
c      COMMON /IOO/IOUT,IUNIT,ONLINE,ITIN,ITOUT,INOD,INMO,IDU,NWRN,
c     +            ONEFILE,FHEADER,BRIEF,IBRIEF,GRAPH,ISUMMR,ICOORD,
c     +            IOERR,NODISPLAY,LBELL
C     ..
C
C
      COMMON /IOO/IOUT,IUNIT,ITIN,ITOUT,INOD,INMO,IDU,NWRN,IBRIEF,
     +        ISUMMR,ICOORD,SERVERFD,SERVERFDIN,NSHUTERR, dnafd, 
     +        ONLINE,ONEFILE,FHEADER,BRIEF,GRAPH,IOERR,NODISPLAY,
     $        LBELL,JPGOUT,SOCKLO,LBEST,dnaout,INMOOPEN,XMLOPEN,
     +        SOCKLOPHIPROF
      COMMON /CIOO/dna_image,xmlfile

C&&*&& end_include  ../inc/ioo.f
      LOGICAL BOXOPEN,RFIXCELL,RFIXDIST
      INTEGER IXP,IYP,NSOL,L,IBUTTON
      INTEGER  XDLSTR,LENSTR
      CHARACTER*120 LINE
      COMMON /ERRFLG/ IERRFLG
      INTEGER IERRFLG
! This is a dummy routine distributed with MOSFLM when REFIX has not been
! requested by the site. If you want REFIX, you must e-mail Andrew Leslie
! andrew@mrc-lmb.cam.ac.uk so that he can send you a copy. All other aspects 
! of MOSFLM should work without REFIX. 
      IF (ONLINE) WRITE(ITOUT,FMT=1000)
      WRITE(LINE,FMT=1000)
 1000 FORMAT(1X,'This option is only available if you e-mail ',
     +'Andrew Leslie and ask for it. Use the new-style indexing ',
     +'instead')
      RETURN
      END
C
C
C
C
      SUBROUTINE GONSYS(ACHSE,S0L,S0G,EG,IER)
      IMPLICIT NONE
      INTEGER    IER,I
      REAL       ACHSE(3),S0L(3),S0G(3),EG(3,3)
      DO 10 I=1,3
      EG(I,2)=ACHSE(I)
10    EG(I,3)=S0L(I)
      CALL UNORM(EG,IER)
      DO 20 I=1,3
20    S0G(I)=EG(1,I)*S0L(1)+EG(2,I)*S0L(2)+EG(3,I)*S0L(3)
      RETURN
      END
C
C
C
C
      SUBROUTINE UNORM(U,IER)
      IMPLICIT NONE
C
C----  normalize matrix u(3,3) to orthogonal form
C
C
      INTEGER    IER,J
      REAL       U(3,3),P,Q
      IER=-1
      P=U(1,2)**2+U(2,2)**2+U(3,2)**2
      Q=U(1,2)*U(1,3)+U(2,2)*U(2,3)+U(3,2)*U(3,3)
      DO 1 J=1,3
1     U(J,3)=U(J,3)*P-U(J,2)*Q
      Q=SQRT(U(1,3)**2+U(2,3)**2+U(3,3)**2)
      IF (Q.LT.0.000001)GO TO 3
      P=SQRT(P)
      DO 2 J=1,3
      U(J,2)=U(J,2)/P
2     U(J,3)=U(J,3)/Q
      U(1,1)=U(2,2)*U(3,3)-U(2,3)*U(3,2)
      U(2,1)=U(3,2)*U(1,3)-U(3,3)*U(1,2)
      U(3,1)=U(1,2)*U(2,3)-U(1,3)*U(2,2)
      IER=0
3     RETURN
      END
C
C
C
C     ======================
      SUBROUTINE MATCOP(A,B)
C     ======================
      IMPLICIT NONE
C
C
C
      REAL    A(3,3),B(3,3)
      INTEGER I,J
C
C
      DO 1 I=1,3
      DO 1 J=1,3
1     B(I,J)=A(I,J)
C
C
      RETURN
      END
C
C
C
C
C     ===============================
      SUBROUTINE INVCEL(CELL,RCELL,V)
C     ===============================
      IMPLICIT NONE
C
C
C
C---- Calculate reciprocal cell parameters from unit cell constants.
C     This routine may also be used to calculate the unit cell
C
C    Parameters from the reciprocal cell constants.
C
C  CELL  - unit cell constants : a,b,c,alfa,beta,gamma           (given)
C          in angstroem and degrees
C RCELL  - calculated reciprocal cell parameters                (result)
C          in reciprocal angstroem and degrees
C   V    - calculated unit cell volume in angstroem**3 .        (result)
C          v<=0.0 indicates an error situation.
C
C
      INTEGER I,J,K
      REAL    CELL(6),RCELL(6),CA(3),SA(3),ARG,V
C
C
C
      DO 10 I=1,3
      ARG=CELL(I+3)/57.29578
      CA(I)=COS(ARG)
10    SA(I)=SIN(ARG)
C
C
      V=CELL(1)*CELL(2)*CELL(3)*
     @  SQRT(1.0+2.0*CA(1)*CA(2)*CA(3)-CA(1)**2-CA(2)**2-CA(3)**2)
      IF (V.LE.0.0)GO TO 30
C
C
      DO 20 I=1,3
      J=MOD(I,3)+1
      K=MOD(I+1,3)+1
      RCELL(I)=CELL(J)*CELL(K)*SA(I)/V
      ARG=57.29578*ACOS((CA(J)*CA(K)-CA(I))/(SA(J)*SA(K)))
      J=NINT(ARG)
      IF (ABS(ARG-J).LT.0.0009)ARG=J
20    RCELL(I+3)=ARG
C
C
30    RETURN
      END
C
C
C
C
C     =============================
      SUBROUTINE METRIC(B,CELL,IER)
C     =============================
      IMPLICIT NONE
C
C
C
C
C     Calculates cell parameters from unit cell basis vectors
C                     W.KABSCH 10-1991
C
C
C   B   - basis vectors of the unit cell (stored row-wise)       (given)
C  CELL - cell parameters (angstrom & degrees)                  (result)
C  IER  - error flag: 0:no error, -1:illegal basis vectors      (result)
C
C
C
      INTEGER    IER,I,J,K
      REAL       B(3,3),CELL(6),R
C
C
C
C
      IER=-1
      DO 10 I=1,3
10    CELL(I)=SQRT(B(I,1)*B(I,1)+B(I,2)*B(I,2)+B(I,3)*B(I,3))
      DO 20 I=2,3
      DO 20 J=1,I-1
      R=CELL(I)*CELL(J)
      IF (R.LE.0.0)GO TO 30
      R=(B(I,1)*B(J,1)+B(I,2)*B(J,2)+B(I,3)*B(J,3))/R
      IF (ABS(R).GE.1.0)GO TO 30
      K=3+6/(I*J)
20    CELL(K)=57.29578*ACOS(R)
      IER=0
30    RETURN
      END
C
C
C
C     ==========================
      SUBROUTINE RFSETMAT(RCELL,A)
C     ==========================
      IMPLICIT NONE
C
C
C
C----  calculate setting matrix from reciprocal unit cell.
C
C RCELL  - RECIPROCAL UNIT CELL PARAMETERS                       (GIVEN)
C          IN RECIPROCAL ANGSTROEM AND DEGREES.
C   A    - SETTING MATRIX IN STANDARD ORIENTATION               (RESULT)
C
C
C
C
C
      INTEGER I
      REAL    RCELL(6),A(3,3),RC(3),RS(3),ARG
      DO 10 I=1,3
      ARG=RCELL(I+3)/57.29578
      RC(I)=COS(ARG)
10    RS(I)=SIN(ARG)
      A(1,1)=RCELL(1)*RS(2)
      A(1,2)=RCELL(2)*(RC(3)-RC(1)*RC(2))/RS(2)
      A(1,3)=0.0
      A(2,1)=0.0
      A(2,2)=SQRT((RCELL(2)*RS(1))**2-A(1,2)**2)
      A(2,3)=0.0
      A(3,1)=RCELL(1)*RC(2)
      A(3,2)=RCELL(2)*RC(1)
      A(3,3)=RCELL(3)
      RETURN
      END
C
C
C
C
C     =================================
      SUBROUTINE RTUMAT(ASTV,BSTV,CSTV)
C     =================================
      IMPLICIT NONE
C
C
C       Given the setting matrix found by auto-indexing, this routine
C       operates the point group matrices on it to find the one which
C       gives the smallest misetting angles which respect to a defined
C       U matrix.
C
C       Note that the array AMAT here is a LOCAL variable.
C
C        The routine DCOSFD is used to work out the unit vector defining
C       the axis of rotation and the angle. If we find such vectors for
C       each matrix, the test matrix which gives the largest dot product
C       with the axis of the target matrix  should have the most similar
C       rotation axis. But there may be an ambiguity about its sense.
C       This can be resolved by finding the matrix with the most
C       similar angle of  rotation about the axis.
C
C
C       NOTES ON ROBUSTNESS:
C     We often find two solutions with dot product of the axis vectors
C    which are identical. But because of rounding errors in REAL*4
C    representation these errors are compounded over the matrix operations
C    and conspire to reduce the identity of these products. Here we allow
C    a tolerance to the fourth decimal place. For added robustness RFMATMUL
C    and INVERS should be double precision. At present single precision
C    seems adequate but in certain unforseen circumstances the wrong
C    solution may be chosen. Either reduce 1E-4 below or use double
C    precision routines from the NAG library.PMcL
C
C
C----  External subroutines required
C
C
C       DECOMP,INVERS,RFMATMUL,MATCOP,LAUSEY,RTOMISSET,DCOSFD,UNORM
C
C
C
C
cc        IMPLICIT NONE
C     .. Array Arguments ..
      REAL ASTV(3),BSTV(3),CSTV(3)
C     ..
C     .. Scalars in Common ..
      REAL CCOM,CCX,CCY,DPHI,F,FLAMDA,ORGX,ORGY,Q,THRESH,DXY
      INTEGER FIXF,IGROUP
      LOGICAL DCOMP,FILM,LCAMC,TARGET
      INTEGER LP,LINOUT
      LOGICAL ONLINE
C     ..
C     .. Arrays in Common ..
      REAL ACHSE,CELL,ED,S0,TARMAT
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BEST,PROD
      REAL BESTA,D,TEST,TRACE,TRACEU
      INTEGER IBEST,IER,J,J1,K,K1,K2,NEQ
C     ..
C     .. Local Arrays ..
      REAL AMAT(3,3),ANG(3,24),BINV(3,3),BMAT(3,3),DUM(3,3),EUTAR(3),
     +     EUTEST(3),I(3,3),L(3,3),LBINV(3,3),LINV(3,3),LUINV(3,3),
     +     PHI(3),RMAT(3,3),UINV(3,3),UMAT(3,3),UTAR(3,3),UTEST(3,3)
      INTEGER MLAUE(432)
C     ..
C     .. External Subroutines ..
      EXTERNAL DCOSFD,DECOMP,INVERS,LAUESY,MATCOP,RFMATMUL,RTOMISSET,
     +         UNORM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /REFCOM/IGROUP,FLAMDA,CELL(6),F,THRESH,Q,DPHI,ED(3,3),
     +       ACHSE(3),S0(3),ORGX,ORGY,FIXF,TARMAT(3,3),TARGET,DCOMP,
     +       LCAMC,CCX,CCY,CCOM,FILM,DXY,FIXCELL
      LOGICAL FIXCELL
      COMMON /RINOUT/ ONLINE,LP,LINOUT
C     ..
      SAVE
C     .. Data statements ..
C
C
      DATA I/1,0,0,0,1,0,0,0,1/
      DATA MLAUE/432*0/
      DATA BEST/1.0E-12/
C     ..
C
C
C
C
C---- Make BEST negative, as in P21 both solutions can have -ve dot product
C     so no solution is chosen !
C
CAL      BEST = 1.0E-12
      BEST = -999.0
C
      DO 10 K = 1,3
        AMAT(K,1) = ASTV(K)
        AMAT(K,2) = BSTV(K)
        AMAT(K,3) = CSTV(K)
   10 CONTINUE
C
C---- DECOMP decomposes AMAT into UMAT*BMAT
C
C          **********************
      CALL DECOMP(AMAT,UMAT,BMAT)
C          **********************
C
C---- If the target matrix was given as an A MATRIX
C     then decompose it to get UTAR
C
      IF (DCOMP) THEN
C
C            ***********************
        CALL DECOMP(TARMAT,UTAR,DUM)
C            ***********************
C
CAL        IF (ONLINE) WRITE (LP,FMT=6000) ((UTAR(K,J),J=1,3),K=1,3)
CAL 6000   FORMAT('U-MATRIX DERIVED FROM INPUT A-MATRIX',/,
CAL     +   3(1X,3F10.5,/))
CAL        WRITE (LINOUT,FMT=6000) ((UTAR(K,J),J=1,3),K=1,3)
      ELSE 
C 
C            *******************
        CALL MATCOP(TARMAT,UTAR)
C            *******************
C
      END IF
C
C          *******************
      CALL MATCOP(UTAR,TARMAT)
      CALL INVERS(UTAR,UINV,D)
      CALL INVERS(BMAT,BINV,D)
C          *******************
C
C---- Get unit vector along axis of rotation (EUTAR) and
C     angle of rotation about that axis (TRACEU) from the
C     target setting matrix UTAR
C
C          *************************
      CALL DCOSFD(UTAR,EUTAR,TRACEU)
C          *************************
C
C---- Read Permutation Matrices for the Point Group
C
C          *****************************
      CALL LAUESY(CELL,IGROUP,NEQ,MLAUE)
C          ******************************
C
C---- Main loop over the point group matrices
C
      J1 = 0
CAL      WRITE (LINOUT,FMT=6030)
CAL 6030 FORMAT(//,1X,
CAL     +     'LAUE OPERATOR    PHIX    PHIY    PHIZ')
CAL      IF (ONLINE) WRITE (LP,FMT=6030)
C
C
      DO 50 K = 1,NEQ
C
C---- Gather point group matrix number K into L(3,3)
C
        DO 30 K1 = 1,3
          DO 20 K2 = 1,3
            J1 = J1 + 1
            L(K2,K1) = MLAUE(J1)
   20     CONTINUE
   30   CONTINUE
C
C---- Calculate the UMAT given this point group matrix.
C
C            ********************
        CALL INVERS(L,LINV,D)
        CALL RFMATMUL(L,UINV,LUINV)
        CALL RFMATMUL(UMAT,LUINV,RMAT)
        CALL RTOMISSET(RMAT,PHI,1)
        CALL RFMATMUL(L,BINV,LBINV)
        CALL RFMATMUL(AMAT,LBINV,UTEST)
C            ************************
C
C---- Renormalise UTEST to reduce errors built up by matrix
C      concatenation.
C
C            ****************
        CALL UNORM(UTEST,IER)
C            ****************
C
C---- Calculate how similar this UMAT  is to the target matrix.
C
C     Get unit vector along axis of rotation for the candidate
C     matrix and the angle of rotation about this axis
C
C            *************************
        CALL DCOSFD(UTEST,EUTEST,TRACE)
C            **************************
C
        TEST = ABS(TRACE-TRACEU)
C
C
        PROD = EUTEST(1)*EUTAR(1) + EUTEST(2)*EUTAR(2) +
     +         EUTEST(3)*EUTAR(3)
C
C---- If the dot product of the axis vectors are larger than the
C     best so far, make the current the best. If they are similar
C     to the fourth decimal place, then take the one describing
C     the smaller angle of rotation
C
C---- Following code changed 9/5/95, previous version would not work
C     correctly in some circumstances (eg p21)
C
        IF (DABS(PROD-BEST).GT.1.0E-4) THEN
C
C---- This is a different solution
C
          IF (PROD.GT.BEST) THEN
            BESTA = TEST
            BEST = PROD
            IBEST = K
          END IF
C
        ELSE
C
C---- This is a similar solution
C
          IF (TEST.LT.BESTA) THEN
            BESTA = TEST
            BEST = PROD
            IBEST = K
          END IF
        END IF
cal   Following is original code
cal     IF ((PROD-BEST).GT.1.0E-4) THEN
cal       BESTA = ABS(TEST)
cal       BEST = DABS(PROD)
cal       IBEST = K
cal     ELSE IF (DABS(PROD-BEST).LE.1.0E-4 .AND.
cal  +           DABS(PROD-BEST).GT.0.0) THEN
cal       IF (TEST.LT.BESTA) THEN
cal         BESTA = ABS(TEST)
cal         BEST = DABS(PROD)
cal         IBEST = K
cal       END IF
cal     END IF
C
C
        DO 40 K1 = 1,3
          ANG(K1,K) = PHI(K1)
   40   CONTINUE
C
C
CAL        WRITE (LINOUT,FMT='(10X,I2,2x,3(f8.2),/)') 
CAL   +                     K,PHI
CAL      IF (ONLINE) WRITE (LP,FMT='(10X,I2,2x,3(f8.2),/)') K,PHI
C
CAL      WRITE (LINOUT,fmt=6010) IBEST,K,((RMAT(K1,J),J=1,3),K1=1,3) 
CAL      IF (ONLINE) WRITE(LP,fmt=6010) IBEST,K,
CAL   +        ((RMAT(K1,J),J=1,3),K1=1,3)
C 6010   FORMAT(1X,'IBEST,K ', 2I5,'  RMAT IS',/,3(1X,3F10.4,/))
   50 CONTINUE
C
C---- End of Main Loop
C
C---- Now operate the IBEST Laue group matrix on the
C     input A matrix
C
      J1 = (IBEST-1)*9
C
C
      DO 70 K1 = 1,3
        DO 60 K2 = 1,3
          J1 = J1 + 1
          L(K2,K1) = MLAUE(J1)
   60   CONTINUE
   70 CONTINUE
C
C          ******************
      CALL RFMATMUL(AMAT,L,DUM)
C          ******************
C
      DO 80 K = 1,3
        ASTV(K) = DUM(K,1)
        BSTV(K) = DUM(K,2)
        CSTV(K) = DUM(K,3)
   80 CONTINUE
C
C
      END
C
C
C
C     ========================
      SUBROUTINE DECOMP(A,U,B)
C     ========================
      IMPLICIT NONE
C
C
C
C---- Calculate U and B given A.
C
C
C
C
C---- Orientation matrix given explicitly
C
C     Set UMAT to (UB)transpose
C
C     .. Array Arguments ..
      REAL A(3,3),B(3,3),U(3,3)
C     ..
C     .. Scalars in Common ..
      REAL CCOM,CCX,CCY,DPHI,F,FLAMDA,ORGX,ORGY,Q,THRESH,DXY
      INTEGER IGROUP,FIXF
      LOGICAL DCOMP,FILM,LCAMC,TARGET
      INTEGER LP,LINOUT
      LOGICAL ONLINE
C     ..
C     .. Arrays in Common ..
      REAL ACHSE,CELL,ED,S0,TARMAT
C     ..
C     .. Local Scalars ..
      REAL D
      INTEGER I,IERR
C     ..
C     .. Local Arrays ..
      REAL CELLC(6),GMAT(3,3),RCELL(6),XMAT(3,3)
C     ..
C     .. External Subroutines ..
      EXTERNAL BMATRX,CHECKU,CLCALC,INVERS,RFMATMUL,TRANSP
C     ..
C     .. Common blocks ..
      COMMON /REFCOM/IGROUP,FLAMDA,CELL(6),F,THRESH,Q,DPHI,ED(3,3),
     +       ACHSE(3),S0(3),ORGX,ORGY,FIXF,TARMAT(3,3),TARGET,DCOMP,
     +       LCAMC,CCX,CCY,CCOM,FILM,DXY,FIXCELL
      LOGICAL FIXCELL
      COMMON /RINOUT/ ONLINE,LP,LINOUT
C     ..
      SAVE
C
C---- U = transp(A), U is returned
C
C          ***********
      CALL TRANSP(U,A)
C          ***********
C
C---- Get reciprocal metric tensor G**-1 in GMAT, = (UB)T.(UB)
C     GMAT is returned
C
C          ****************
      CALL RFMATMUL(U,A,GMAT)
C          ****************
C
C---- Get reciprocal cell dimensions RCELL
C
C          ******************
      CALL CLCALC(RCELL,GMAT)
C          ******************
C
C---- Invert to get real metric tensor G in U
C     U is returned, D is determinant
C
C          ****************
      CALL INVERS(GMAT,U,D)
C          ****************
C
C      IF (ABS(D).LE.0.000001) THEN
      IF (D.EQ.0.0) THEN
        WRITE (LINOUT,FMT=6000)
        IF (ONLINE) WRITE (LP,FMT=6000)
      ELSE
C
C---- Get real cell dimensions from metric tensor
C
C            ***************
        CALL CLCALC(CELLC,U)
C            ***************
C
        DO 10 I = 1,3
          CELLC(I) = CELLC(I)*FLAMDA
   10   CONTINUE
C
C---- Rebuild cell orthogonalization matrix B in BMAT
C
C            ****************************
        CALL BMATRX(B,RCELL,CELLC,FLAMDA)
C            ****************************
C
C---- Get matrix U = (UB).(B**-1) to UMAT
C
C            ****************
        CALL INVERS(B,XMAT,D)
        CALL RFMATMUL(A,XMAT,U)
        CALL UNORM(U,IERR)
        CALL CHECKU(U)
C            *****************
C
      END IF
C
C---- Format statements
C
 6000 FORMAT (//' !!!! Zero reciprocal cell volume from orientation ma',
     +       'trix !!!!',//)
C
C
      END
C
C
C
C     ==========================
      SUBROUTINE INVERS(A,B,DET)
C     ==========================
      IMPLICIT NONE
C
C
C
C      A    - 3 x 3 matrix                                       (given)
C      B    - 3 x 3 matrix is calculated inverse matrix of a    (result)
C      DET  - calculated determinant of given matrix a          (result)
C
C
C
      INTEGER   I,J,K
      REAL      A(3,3),B(3,3),DET
C
C
C
      DO 1 I=1,3
      J=I+1
      IF (J.GT.3)J=J-3
      K=I+2
      IF (K.GT.3)K=K-3
      B(I,1)=REAL(DBLE(A(2,J))*DBLE(A(3,K))-DBLE(A(2,K))*DBLE(A(3,J)))
      B(I,2)=REAL(DBLE(A(3,J))*DBLE(A(1,K))-DBLE(A(3,K))*DBLE(A(1,J)))
1     B(I,3)=REAL(DBLE(A(1,J))*DBLE(A(2,K))-DBLE(A(1,K))*DBLE(A(2,J)))
C
C
      DET=REAL(DBLE(A(1,1))*DBLE(B(1,1))+DBLE(A(1,2))*DBLE(B(2,1))+
     +         DBLE(A(1,3))*DBLE(B(3,1)))
      IF (DET.EQ.0.0)GO TO 3
C
C
      DO 2 I=1,3
      DO 2 J=1,3
2     B(I,J)=B(I,J)/DET
3     RETURN
      END
C
C
C
C
C     =====================================
      SUBROUTINE LAUESY(CELL,IGROUP,NEQ,ML)
C     =====================================
      IMPLICIT NONE
C
C
C
C
C
C  Obtain equivalent positions of laue-group from space-group number.
C  The space-group number for your crystal can be looked up in the
C  International tables for x-ray crystallography vol.i
C
C                 W.KABSCH    7-1986
C
C
C  CELL  - unit cell parameters in angstroem and degrees         (given)
C IGROUP - space group number as obtained from the international (given)
C          tables for x-ray crystallography vol.i
C   NEQ  - number of equivalent positions in laue group         (result)
C          derived from space-group number. a value of zero
C          indicates an illegal space-group number or illegal
C          unit cell parameters. otherwise neq assumes values
C          between 1 and 24.
C   ML   - integer array of length at least neq*9 specified     (result)
C          in the main program, to describe symmetry. all
C          friedel related positions are omitted. the largest
C          number of neq is 24, as mentioned. therefore, the
C          calling program should declare
C          integer ml(216)  to handle all possible cases!!!
C          the operators returned by this subroutine are
C          proper rotations. to obtain the complete laue-group
C          you only have to add the 'neq' friedel mates to
C          the matrices given by this routine. usually this is
C          not wanted.
C          each operator is represented by 9 consecutive
C          numbers as follows:
C          if x',y',z' is an equivalent position of x,y,z then
C          x'=x*ml(1)+y*ml(2)+z*ml(3)
C          y'=x*ml(4)+y*ml(5)+z*ml(6)
C          z'=x*ml(7)+y*ml(8)+z*ml(9)
C          the next operator occupies positions ml(10)...ml(18),
C          and so on...
C
C
C
C
      REAL      CELL(6),EPS
      INTEGER   ML(*),NEQ,IGROUP,IG,I,J,K,L,M
      INTEGER*2 MAT(288),MATN(56),NM(14),NA(14)
      PARAMETER (EPS=0.001)
      DATA MAT/ 1,3*0,1,3*0,1,  -1,3*0,-1,3*0,1,  1,3*0,-1,3*0,-1,
     1         -1,3*0,1,3*0,-1,  0,1,0,-1,4*0,1,  0,-1,0,1,4*0,1,
     2          0,-1,0,-1,4*0,-1,  0,1,0,1,4*0,-1,  0,0,1,1,3*0,1,0,
     3          0,0,1,-1,3*0,-1,0,  0,1,3*0,1,1,0,0,  0,1,3*0,-1,-1,0,0,
     4        0,0,-1,1,3*0,-1,0,  0,0,-1,-1,3*0,1,0,  0,-1,3*0,1,-1,0,0,
     5      0,-1,3*0,-1,1,0,0,  -1,4*0,-1,0,-1,0,  0,0,-1,0,-1,0,-1,0,0,
     6          -1,4*0,1,0,1,0,  0,0,-1,0,1,0,1,0,0,  1,4*0,-1,0,1,0,
     7      0,0,1,0,-1,0,1,0,0,  1,4*0,1,0,-1,0,  0,0,1,0,1,0,-1,0,0,
     8        0,-1,0,1,-1,3*0,1,  -1,1,0,-1,4*0,1,  -1,0,0,-1,1,3*0,-1,
     9         1,-1,0,0,-1,3*0,-1,  0,1,0,-1,1,3*0,1,  1,-1,0,1,4*0,1,
     $         1,0,0,1,-1,3*0,-1,  -1,1,0,0,1,3*0,-1/
      DATA MATN/1,4,3,2,9,10,11,12,13,14,15,16,
     1          5,6,7,8,17,18,19,20,21,22,23,24,
     2          1,2,5,6,3,4,7,8,
     3          1,25,26,8,27,28,  1,9,11,7,17,18,
     4          1,25,26,2,29,30,  7,31,32,8,27,28/
      DATA NM/  1,2,4,12,24, 2,4,8, 3,6, 3,6, 6,12/
      DATA NA/  0,0,0,0,0,  24,24,24, 32,32,  38,38, 44,44/
C
C
C
C----  check unit cell parameters
C
      IG=0
      NEQ=0
      DO 10 J=1,6
      IF (CELL(J).LE.0.0)GO TO 30
10    CONTINUE
C
C----  find laue group from space-group number and unit cell parameters
C
C       ..... triclinic axes
C
      IF ((IGROUP.GT.0).AND.(IGROUP.LT.3))IG=1
C
C       ..... monoclinic axes
C
      IF ((IGROUP.GT.2).AND.(IGROUP.LT.16).AND.
     @    (ABS(CELL(4)-90.0).LT.EPS))THEN
	 IF (ABS(CELL(5)-90.0).LT.EPS)IG=6
	 IF (ABS(CELL(6)-90.0).LT.EPS)IG=2
      ENDIF
C
C       ..... orthorhombic axes
C
      IF ((ABS(CELL(4)-90.0).LT.EPS).AND.(ABS(CELL(5)-90.0).LT.EPS)
     @    .AND.(ABS(CELL(6)-90.0).LT.EPS))THEN
	 IF ((IGROUP.GT.15).AND.(IGROUP.LT.75))IG=3
C
C        ..... tetragonal axes
C
	 IF (ABS(CELL(1)-CELL(2)).LT.EPS)THEN
	    IF ((IGROUP.GT.74).AND.(IGROUP.LT.89))IG=7
	    IF ((IGROUP.GT.88).AND.(IGROUP.LT.143))IG=8
C
C         ..... cubic axes
C
	    IF (ABS(CELL(1)-CELL(3)).LT.EPS)THEN
	       IF ((IGROUP.GT.194).AND.(IGROUP.LT.207))IG=4
	       IF ((IGROUP.GT.206).AND.(IGROUP.LT.231))IG=5
	    ENDIF
	 ENDIF
      ENDIF
C
C         ..... hexagonal axes
C
      IF ((ABS(CELL(1)-CELL(2)).LT.EPS).AND.(ABS(CELL(4)-90.0).LT.EPS)
     @ .AND.(ABS(CELL(5)-90.0).LT.EPS).AND.(ABS(CELL(6)-120.0).LT.EPS))
     & THEN
	 IF ((IGROUP.GT.142).AND.(IGROUP.LT.149))IG=9
	 IF ((IGROUP.GT.148).AND.(IGROUP.LT.168))IG=10
	 IF ((IGROUP.GT.167).AND.(IGROUP.LT.177))IG=13
	 IF ((IGROUP.GT.176).AND.(IGROUP.LT.195))IG=14
      ENDIF
C
C           ..... rhombohedral axes
C
      IF ((ABS(CELL(1)-CELL(2)).LT.EPS).AND.
     1    (ABS(CELL(1)-CELL(3)).LT.EPS).AND.
     2    (ABS(CELL(4)-CELL(5)).LT.EPS).AND.
     3    (ABS(CELL(4)-CELL(6)).LT.EPS))THEN
	 IF ((IGROUP.EQ.146).OR.(IGROUP.EQ.148))IG=11
	 IF ((IGROUP.EQ.155).OR.(IGROUP.EQ.160).OR.(IGROUP.EQ.161).OR.
     @       (IGROUP.EQ.166).OR.(IGROUP.EQ.167))IG=12
cdms21.09.2004    Allow space groups 1146 and 1155.
         IF ((IGROUP.EQ.1146).OR.(IGROUP.EQ.1148))IG=11
         IF ((IGROUP.EQ.1155).OR.(IGROUP.EQ.1160).OR.
     @       (IGROUP.EQ.1161).OR.(IGROUP.EQ.1166).OR.
     @       (IGROUP.EQ.1167))IG=12
      ENDIF
C
C----  get matrices of equivalent positions for laue group
C
      IF (IG.LT.1)GO TO 30
      NEQ=NM(IG)
      K=NA(IG)
      I=0
      DO 20 J=1,NEQ
      L=9*(MATN(J+K)-1)
      DO 20 M=1,9
      I=I+1
20    ML(I)=MAT(M+L)
30    RETURN
      END
C
C
C
C     ========================
      SUBROUTINE RFMATMUL(A,B,C)
C     ========================
      IMPLICIT NONE
C
C
C---- Multiply 2 3x3 matrices
C
C      C = AB
C
C     .. Array Arguments ..
C
      REAL    A(3,3),B(3,3),C(3,3)
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION AD(3,3)
      INTEGER I,J,K
C
      DO 30 I = 1,3
        DO 20 J = 1,3
          AD(I,J) = 0.0D0
          DO 10 K = 1,3
            AD(I,J) = DBLE(A(I,K))*DBLE(B(K,J)) + AD(I,J)
   10     CONTINUE
          C(I,J) = REAL(AD(I,J))
   20   CONTINUE
   30 CONTINUE

      RETURN
      END
