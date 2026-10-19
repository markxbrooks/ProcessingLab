C     Phase C: fixed-form F77 call into predict_still_cpp_ (production path).
C     Analytical stills geometry is filled explicitly (no use_fixture).
C
      PROGRAM TEST_PREDICT_STILL_F77
      IMPLICIT NONE
      INTEGER NMAX, NSPOT, IERR, NF, NS
      PARAMETER (NMAX=8192)
      REAL CELL(6)
      DOUBLE PRECISION DIST, WAVE, PX, PY, BX, BY, PHI, DMIN, ETOL
      DOUBLE PRECISION XMM(NMAX), YMM(NMAX), PHIB(NMAX), PHIE(NMAX)
      INTEGER H(NMAX), K(NMAX), L(NMAX), IFULL(NMAX)
      INTEGER EXPECT
      PARAMETER (EXPECT=656)

      EXTERNAL PREDICT_STILL_CPP

C     Analytical stills fixture (testdata/oracle/predict_still.inputs.txt)
      CELL(1) = 58.0
      CELL(2) = 58.0
      CELL(3) = 150.0
      CELL(4) = 90.0
      CELL(5) = 90.0
      CELL(6) = 90.0
      DIST = 250.0D0
      WAVE = 1.0D0
      PX = 0.172D0
      PY = 0.172D0
      NF = 1475
      NS = 1679
      BX = 0.0D0
      BY = 0.0D0
      PHI = 0.0D0
      DMIN = 3.5D0
      ETOL = 0.0025D0

      CALL PREDICT_STILL_CPP(CELL, DIST, WAVE, PX, PY, NF, NS, BX, BY,
     +     PHI, DMIN, ETOL, H, K, L, XMM, YMM, PHIB, PHIE, IFULL, NMAX,
     +     NSPOT, IERR)

      IF (IERR .NE. 0) THEN
         WRITE (*,*) 'FAIL: IERR=', IERR
         STOP 1
      END IF
      IF (NSPOT .NE. EXPECT) THEN
         WRITE (*,*) 'FAIL: NSPOT=', NSPOT, ' expected', EXPECT
         STOP 1
      END IF
      IF (ABS(CELL(1) - 58.0) .GT. 1.0E-4) THEN
         WRITE (*,*) 'FAIL: CELL(1)=', CELL(1)
         STOP 1
      END IF

      WRITE (*,'(A,I0)') 'test_predict_still_f77: OK n=', NSPOT
      STOP 0
      END
