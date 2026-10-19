C     Fixed-form F77-style call into PREDICT_CPP (amalgamation calling convention).
      PROGRAM TEST_PREDICT_CPP_HOOK_F77
      IMPLICIT NONE
      INTEGER NMAX, EXPECT
      PARAMETER (NMAX=8192)
      REAL PSI(3), CELL(6), AMAT(3,3)
      REAL XTOFD, WAVE, RAST, RASTY, PHI0, PHI1, DMIN, ETOL
      REAL ASTAR, CSTAR
      REAL XG(NMAX), YG(NMAX), PHIG(NMAX), PHIWG(NMAX)
      INTEGER NF, NS, NSPOT, IFLAGPR
      INTEGER*2 IHG(NMAX), IKG(NMAX), ILG(NMAX), IMG(NMAX)
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
      PHI0 = 0.0
      PHI1 = 0.0
      DMIN = 3.5
      ETOL = 0.0025
      PSI(1) = 0.0
      PSI(2) = 0.0
      PSI(3) = 0.0
      ASTAR = 1.0 / 58.0
      CSTAR = 1.0 / 150.0
      AMAT = 0.0
      AMAT(1,3) = CSTAR
      AMAT(2,2) = -ASTAR
      AMAT(3,1) = ASTAR
C
      CALL PREDICT_CPP(PSI, IFLAGPR, CELL, XTOFD, WAVE, RAST, RASTY,
     +     NF, NS, PHI0, PHI1, DMIN, ETOL, AMAT, NMAX, IHG, IKG, ILG,
     +     XG, YG, PHIG, PHIWG, IMG, NSPOT)
C
      IF (IFLAGPR .NE. 0) THEN
        WRITE(*,'(A,I0)') 'FAIL: iflagpr=', IFLAGPR
        STOP 1
      END IF
      IF (NSPOT .LE. 0) THEN
        WRITE(*,'(A,I0)') 'FAIL: nspot=', NSPOT
        STOP 1
      END IF
      WRITE(*,'(A,I0)') 'test_predict_cpp_hook_f77: OK n=', NSPOT
      STOP 0
      END
