C     Fixed-form F77-style call into PREDICT_CPP (amalgamation calling convention).
C
      PROGRAM TEST_PREDICT_CPP_HOOK_F77
      IMPLICIT NONE
      INTEGER NMAX, EXPECT
      PARAMETER (NMAX=8192, EXPECT=656)
      REAL PSI(3), CELL(6), XTOFD, WAVE, RAST, RASTY, PHI, DMIN, ETOL
      REAL XG(NMAX), YG(NMAX), PHIG(NMAX), PHIWG(NMAX)
      INTEGER*2 IHG(NMAX), IKG(NMAX), ILG(NMAX), IMG(NMAX)
      INTEGER NF, NS, NSPOT, IFLAGPR
      EXTERNAL PREDICT_CPP
C
      CELL(1) = 58.0
      CELL(2) = 58.0
      CELL(3) = 150.0
      CELL(4) = 90.0
      CELL(5) = 90.0
      CELL(6) = 90.0
      XTOFD = 25000.0
      WAVE = 1.0
      RAST = 0.172
      RASTY = 0.172
      NF = 1475
      NS = 1679
      PHI = 0.0
      DMIN = 3.5
      ETOL = 0.0025
      PSI(1) = 0.0
      PSI(2) = 0.0
      PSI(3) = 0.0
C
      CALL PREDICT_CPP(PSI, IFLAGPR, CELL, XTOFD, WAVE, RAST, RASTY,
     +     NF, NS, PHI, DMIN, ETOL, NMAX, IHG, IKG, ILG, XG, YG, PHIG,
     +     PHIWG, IMG, NSPOT)
C
      IF (IFLAGPR .NE. 0) THEN
        WRITE(*,'(A,I0)') 'FAIL: iflagpr=', IFLAGPR
        STOP 1
      END IF
      IF (NSPOT .NE. EXPECT) THEN
        WRITE(*,'(A,I0)') 'FAIL: nspot=', NSPOT
        STOP 1
      END IF
      WRITE(*,'(A,I0)') 'test_predict_cpp_hook_f77: OK n=', NSPOT
      STOP 0
      END
