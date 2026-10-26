      PROGRAM WRAPPER
C
C---- Gets information from a CBF file via Paul Ellis's C routines and
C     puts them into a form that FORTRAN programs can deal with. Later, 
C     this will become a callable routine for MOSFLM
C
        INTEGER ID, INDEX
        INTEGER I,IERR
        INTEGER COLUMN, ROW
        DOUBLE PRECISION PIXEL_SIZE, GAIN, WAVELENGTH, DISTANCE
        INTEGER OVERLOAD, DIMENSION(2), PRECEDENCE(2)
        CHARACTER DETECTOR
        CHARACTER DETECTOR_CHARACTER
        CHARACTER DETECTOR_ID(64)
        CHARACTER ASL_INDEX(2)
        CHARACTER DIRECTION(2), ARRAY_ID, RADIATION
        CHARACTER*80 ARGV
        INTEGER MAXPIX, MINPIX, MINX, MINY, MAXX, MAXY 
C
C---- The following is used for transferring values between the C stuff
C     and Fortran for CBF images
C
        INCLUDE '../../inc/cbfinc.h'
C
C     Key to contents of arrays (but N.B. these may well change before 
C     final implementation):
C
C     CBF_INT(1)    = number of pixels in horizontal direction 
C     CBF_INT(2)    = number of pixels in vertical direction
C     CBF_INT(3)    = overload value
C     CBF_INT(4)    = unused so far 
C     CBF_INT(5)    = unused so far 
C     CBF_INT(6)    = unused so far 
C     CBF_INT(7)    = unused so far 
C     CBF_INT(8)    = unused so far 
C     CBF_INT(9)    = unused so far 
C     CBF_INT(10)   = unused so far 
C     CBF_INT(11)   = unused so far 
C     CBF_INT(12)   = unused so far 
C     CBF_INT(13)   = unused so far 
C     CBF_INT(14)   = unused so far 
C     CBF_INT(15)   = unused so far 
C     CBF_INT(16)   = unused so far 
C
C     CBF_REAL(1)   = wavelength (Angstroms)
C     CBF_REAL(2)   = crystal to detector distance (metres)
C     CBF_REAL(3)   = pixel size in horizontal direction (metres)
C     CBF_REAL(4)   = pixel size in vertical direction (metres)
C     CBF_REAL(5)   = detector gain
C     CBF_REAL(6)   = corresponds to OMEGAFD
C     CBF_REAL(7)   = unused so far
C     CBF_REAL(8)   = unused so far
C     CBF_REAL(9)   = unused so far
C     CBF_REAL(10)  = unused so far
C     CBF_REAL(11)  = unused so far
C     CBF_REAL(12)  = unused so far
C     CBF_REAL(13)  = unused so far
C     CBF_REAL(14)  = unused so far
C     CBF_REAL(15)  = unused so far
C     CBF_REAL(16)  = unused so far
C
C     CBF_CHAR(1)   = horizontal rate of change (slow or fast)
C     CBF_CHAR(2)   = horizontal direction of change (increasing or decreasing)
C     CBF_CHAR(3)   = vertical rate of change (slow or fast) <> CBF_CHAR(1)
C     CBF_CHAR(4)   = vertical direction of change (increasing or decreasing)
C     CBF_CHAR(5)   = detector name
C     CBF_CHAR(6)   = unused so far
C     CBF_CHAR(7)   = unused so far
C     CBF_CHAR(8)   = unused so far
C     CBF_CHAR(9)   = unused so far
C     CBF_CHAR(10)  = unused so far
C     CBF_CHAR(11)  = unused so far
C     CBF_CHAR(12)  = unused so far
C     CBF_CHAR(13)  = unused so far
C     CBF_CHAR(14)  = unused so far
C     CBF_CHAR(15)  = unused so far
C     CBF_CHAR(16)  = unused so far
C
C     CBF_LOG(1)    = corresponds to INVERTX
C     CBF_LOG(2)    = unused so far
C     CBF_LOG(3)    = unused so far
C     CBF_LOG(4)    = unused so far
C     CBF_LOG(5)    = unused so far
C     CBF_LOG(6)    = unused so far
C     CBF_LOG(7)    = unused so far
C     CBF_LOG(8)    = unused so far
C     CBF_LOG(9)    = unused so far
C     CBF_LOG(10)   = unused so far
C     CBF_LOG(11)   = unused so far
C     CBF_LOG(12)   = unused so far
C     CBF_LOG(13)   = unused so far
C     CBF_LOG(14)   = unused so far
C     CBF_LOG(15)   = unused so far
C     CBF_LOG(16)   = unused so far

        WRITE(*,FMT=200)
        READ(*,FMT=100)ARGV
 100    FORMAT(80A)
 200    FORMAT(' Enter CBF filename: ')
        I = 1
        DO 1010 WHILE (argv(i:i).ne.' ')
           I= I + 1
 1010   enddo
        argv(i:i) = char(0)
        IERR = 0
        CALL CBFWRAP(IERR,CBF_INT,CBF_DOUBLE,CBF_CHAR,ARGV)
        IF(IERR.GT.0)THEN
           print*,'CBF file not found!'
        else
           print*,'CBF_INT',cbf_int
           print*,'CBF_DOUBLE',cbf_double
           print*,'CBF_CHAR:'
           do 1000 i=1,16,1
              print*,cbf_char(i)
 1000      enddo
        endif
        END
