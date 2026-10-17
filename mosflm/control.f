C---- $Id: control.f,v 1.389 2021/08/19 10:54:21 andrew Exp $
C     
C==   CONTROL ==
      SUBROUTINE CONTROL(FIRSTTIME,IFIRSTPACK,NEWGENF,GENOPEN,RPTFIRST,
     +     CELLSTR,MODE)
C     =================================================================
C     
C     test cvs
      IMPLICIT NONE
C     Last format 7880 label 960
C     POWDER BLOCK EXTEND FROM LINE 9758 TO 11123
C     
C     FIRSTTIME        Input, not changed
C     
C     Set in MAIN or MXDSPL. Initially TRUE, set FALSE
C     after STRATEGY, TESTGEN, POWDER (image display), or
C     before starting actual integration in MAIN. However,
C     it is set TRUE if repeating an entire postref segment
C     run.Controls level of initialisation of parameters,
C     whether commands are read from stored input lines or
C     input stream, referring a matrix to a previous one,
C     opening generate file and returning without printing
C     all parameters.
C     
C     IFIRSTPACK  Set and returned

C     This is the pack counter for the first image to be
C     processed in this run. Note that this is a serial
C     counter, and is NOT the image number.
C     
C     NEWGENF     Input, but can be changed.

C     Set TRUE if want to call START. If a generate file is
C     currently open, then it is initially set FALSE. If a
C     GENFILE keyword is given, it is set true. Used to
C     control closure (and subsequent opening) of both
C     generate and MTZ files.
C     
C     GENOPEN     Input, but can be changed and returned.
C     
C     TRUE if a generate file is currently open.
C     
C     RPTFIRST      Input, not changed.
C     True if repeating a multi-segment post-refinement
C     because of an excessive shift in cell
C     parameters. This controls, for example, restoring
C     refined cell parameters in preference to those from
C     MATRIX or CELL keywords, controlling reading of
C     KEYWORDS from saved list, closing generate and MTZ
C     files. Note that it is only TRUE for the first pack
C     of the first segment in the repeat run, and FALSE
C     for the first pack of all subsequent segments.
C     
C     CELLSTR
C     
C     MODE         If CONTROL is called from MXDSPL, then the value of
C     MODE determines the flow through CONTROL.
C     =0     When called from MOSFLM, normal route
C     =1     When doing a prediction after autoindexing from MXDSPL.
C      -1    When doing a prediction for new GUI during integration
C     =2     When integrating images selected in MXDSPL
C     =3     When reading keyword input via display window
C     =4     When doing multiseg post-refinement
C     =10    When doing a strategy run from window
C     
C---- This reads keyworded information on generate file name,
C     packs (and films in packs) to be processed, and
C     optional flags such as filmplot and avprofile
C     uses mike levitts parser routines
C     note limit of MAXPAX packs per generate file
C     
C     ITYP    =1 Text
C     =2 A number
C     =3 A quoted token
C     
C     Last modified 10/7/89
C     Last modified 17/8/89
C     
C---- NSER is the number of SERIAL keywords read.
C---- NPACKS is no. of packs in current serial card whereas NPACK is
C     the number of packs in total
C     ISTARTP is a pointer for the first pack of the next PROCESS/SERIAL
C     keyword.
C---- IFIRSTPACK is the pack counter for the first image in this "run"
C     IFIRST is the pack counter for the first pack in
C     current "pack" card there may be more than one
C     "pack" card for each "run"
C     
C---- NTLINE is incremented each time a line is read, unless it is a
C     @FILENAME
C     Note that NTLINE is a pointer for the NEXT line, so it is one
C     greater than the actual number of lines stored. It is initialised
C     to 1 and cannot be greater than 200.
C     
C     NRLINE is set to NTLINE after the first RUN keyword has been read
C     for
C     a multisegment post-refinement run. Thus it actually points to
C     the line AFTER the RUN keyword.
C     
C     NLINE  is a pointer to the next line to be read in from lines
C     stored in
C     array NLINE.
C
C     DEBUG keyword is used to provide debug output for particular
C     subroutines, eg 'DEBUG CONTROL'. The subkeyword is the name of
C     the subroutine producing the debug output.

C     There are also special DEBUG subkeywords for debugging other
C     aspects of the program. These are currently:

C     DEBUG LATTICES (DEBUG(88)) Debug for multiple lattice work
C                            writing XML.
C     DEBUG SOCKET   (DEBUG(89)) Debug when a SOCKET_READ is performed
C     DEBUG DSPCHOOS (DEBUG(90)  
C     DEBUG IMAG     (DEBUG(91)) Statistics on image contents from DSPSCL
C     DEBUG SUMIMAGE (DEBUG(92)) Debug on image summation (from several
C                                subroutines, including SUMIMAGES
C     DEBUG ROLLING  (DEBUG(93)) Rolling shutter mode
C     DEBUG          (DEBUG(94)) Phioffset for rolling shutter
C     DEBUG          (DEBUG(95)) Normal matrix
C     DEBUG PIXEL    (DEBUG(96)) Dumps "pick" values to logfile
C     DEBUG DSPS     (DEBUG(97)) Looks for isolated pixels with value -1 
C                                (for Eiger)
C     DEBUG CLUSTER  (DEBUG(98)) Debug finding clusters of rejected backgorund pixels
C     DEBUG FLOW     (DEBUG(99)) Provides output on the flow of the program
C                            ie called and calling subroutines including
C     DEBUG SIGMA    (DEBUG(100)) Flag for new sigma calculation (not used)
C     DBEUG IMAGE <no> (DEBUG(108) Only turn on debug after image "no".
C     DEBUG ADJACENT (DEBUG(109)   Provides information of profiles where adjacent spots
C                                 appear to be intruding into the peak region of a profile.
C
C     .. Parameters ..
C&&*&& include  ../inc/parameter.f
C
C $Id: parameter.f,v 1.34 2020/12/15 20:08:54 andrew Exp $
C
C--- awk generated include file  parameter.h
C---- START of include file parameter.h
C
C       PARAMETERS
C       IYLENGTH.. maximum number of I*2 words of data in the 
C                  "fast" (ie most rapidly changing) direction in the
C                  digitised image. This will be HALF the number of pixels
C                  for film data (each pixel is stored in one byte)
C                  but will equal the number of pixels for IP data.
C       IXWDTH...  The maximum number of "stripes" of data in the scanned image
C                  ie the number of pixels in the "slow" direction
C                  (This is the Y direction in the MOSFLM convention)
C                  Note that the array "IMAGE" used to store the image is
C                  declared as size IYLENGTH*IXWDTH I*2 words for IP data
C                  and 2*IYLENGTH*IXWDTH BYTES for film data.
C                  If this exceeds the
C                  available memory, set ixwdth=1, recompile program
C                  and use keyword "NOCORE" when running program.
C                  Note that for the POSTREF and ADDPART options,
C                  two images have to be stored in memory at once so
C                  IXWDTH should be twice the number of records in an
C                  image.
C       MAXHEAD    maximum length of image header (in 4 byte words)
C       NREFLS.... maximum number of spots per film in generate file (250,000)
C       MAXBOX.... maximum number of pixels in measurement box (3000)
C       MAXDIM.... maximum box size in either direction (pixels) (41)
C       MAXPAX.... maximum number of packs per generate file (10000)
C       MXDOV2..... maxdim/2
C       MAXBUFF... maximum size of buffer (I*2) for storing ods
C                  of active spots in subroutine meas(42,000). 
C                  Must be .GE.  2 * MAXBOX*(NNLINE-1) for subroutine process,
C                  where it is an I*4 array, everywhere else I*2.
C       MREF...... maximum number of reflections to be used in post
C                  refinement (6000)
C       NEXPAND... maximum number of expansions of the input measurement
C                  box (2)
C       NMASKS.... maximum number of different profiles (25). Note the
C                  connection between this parameter and NNLINE
C       NVECT..... maximum number of vectors for storing scanned image
C                  in filmplot (10000)
C       NIMAX..... maximum number of images to be used together in 
C                  postrefinement (NADD or WIDTH options) (400)
C       NNLINE...  maximum number of boundary lines for setting up
C                  the areas for profile fitting. The maximum possible
C                  number of standard profiles will be (NNLINE-1)**2
C                  although for a circular detector the actual number
C                  may be less than this as some boxes will lie entirely
C                  outside the detector.
C       NREJMAX... Maximum number of rejected background pixels, resulting
C                  either from overlap of adjacent spots or outliers from
C                  the background plane
C       NSPOTS...  Maximum number of found spots (for autoindexing) that can
C                  be stored (for all images). Also maximum number in 
C                  a file wriitten by IMSTILLS that can be 
C                  stored/displayed/edited. THis must be an even number
C       MCOLS....  Number of columns in output MTZ file
C       MCOLSTR..  Number of columns in output MTZ file for strategy option
C     
C       MTZ Orientation block
C        MBLENG is total length of block, MBLINT, MBLREA are numbers
C          of integers & reals
C       NDETREF    Maximum number of reflections used for detector parameter
C                  refinement (previously 62)
C       NRPAR....  Maximum number of refineable parameters for detector
C                  positional refinement (subroutine RDIST)
C       NSEGMAX..  Maximum number of segments in STRATEGY
C       NSEGMAXPR  Maximum number of segments in POST-REFINEMENT
C       MULTMAX... Maximum number of observations with same hkl in COMPLETE
C       MAXDIFF... Maximum number of different packs that a given hkl occurs on
C       NRESBIN... Maximum number of resolution bins (COMPLETE)
C
C       MAXIMG...  Maximum number of images that can be read in using the IMAGE
C                  keyword or the "read Image" menu option.
C       MXSPOT...  Maximum number of spots that can be found on one image
C                  (before rejection on spot size).
C       MXCENT...  Maximum number of active spots during spot finding 
C                  (findspots)
C
C       MGRA, NGRA... maximum number of reflections and images over which a 
C                     reflection can be spread for postrefinement.
C       MAXBGBOX    Maximum number of local background background boxes 
C                   in X and Y
C
C       MAXNPEAKS   Maximum number of local maxima within one spot
C       MAXEXCL     Maximum number of excluded areas of the same type
C                   eg reactangles, quadrilaterals, arcs, resolution shells
C
C       NLAUEG      Number of Laue Groups in indexing routines
C
C       MAXLATT     Maximum number of lattices we can find in indexing
C
C       MAXIWARN    Array limit of first element of IWARN
C
C     .. Parameters ..
      INTEGER IXWDTH
c      PARAMETER (IXWDTH=8192)
      PARAMETER (IXWDTH=12288)
      INTEGER IYLENGTH
c      PARAMETER (IYLENGTH=4096)
      PARAMETER (IYLENGTH=6144)
      INTEGER MAXHEAD
      PARAMETER (MAXHEAD=5120)
      INTEGER MAXBOX
      PARAMETER (MAXBOX=3000)
      INTEGER MAXDIM
      PARAMETER (MAXDIM=41)
      INTEGER MAXPAX
      PARAMETER (MAXPAX=10000)
      INTEGER MXDOV2
      PARAMETER (MXDOV2=MAXDIM/2)
      INTEGER NEXPAND
      PARAMETER (NEXPAND=2)
      INTEGER NNLINE
      PARAMETER (NNLINE=8)
      INTEGER NMASKS
      PARAMETER (NMASKS=49)
      INTEGER MAXBUFF
      PARAMETER (MAXBUFF=42000)
      INTEGER NREFLS
      PARAMETER (NREFLS=250000)
      INTEGER MREF
      PARAMETER (MREF=6000)
      INTEGER NVECT
      PARAMETER (NVECT=10000)
      INTEGER NIMAX
      PARAMETER (NIMAX=400)
      INTEGER NREJMAX
      PARAMETER (NREJMAX=1500)
      INTEGER NSPOTS
      PARAMETER (NSPOTS=5000)
      INTEGER MCOLS
      PARAMETER (MCOLS=18)
      INTEGER MCOLSTR
      PARAMETER (MCOLSTR=6)
      INTEGER NREFSTR
C
C---- Each reflection for strategy run needs MCOLSTR I*2 words
C     plus an I*4 word for the merging
C
      PARAMETER (NREFSTR = IXWDTH*IYLENGTH/(2+MCOLSTR))
      INTEGER MBLENG,MBLINT,MBLREA
      PARAMETER (MBLENG=185,MBLINT=29,MBLREA=156)
      INTEGER NDETREF
      PARAMETER (NDETREF=400)
      INTEGER NRPAR
      PARAMETER (NRPAR=14)
      INTEGER NSEGMAX
      PARAMETER (NSEGMAX=100)
      INTEGER NSEGMAXPR
      PARAMETER (NSEGMAXPR=100)
      INTEGER MULTMAX
      PARAMETER (MULTMAX=100)
      INTEGER MAXDIFF
      PARAMETER (MAXDIFF=100)
      INTEGER NRESBIN
      PARAMETER (NRESBIN=20)
      INTEGER MXSPOT
      PARAMETER (MXSPOT=5000)      
      INTEGER MAXIMG
      PARAMETER (MAXIMG=100)      
      INTEGER NPIXBG
      PARAMETER (NPIXBG=51)
      INTEGER MXCENT
      PARAMETER (MXCENT=500)
      INTEGER    NGRA,MGRA
      PARAMETER (NGRA=100)
      PARAMETER (MGRA=50000)
      INTEGER MAXBGBOX
      PARAMETER (MAXBGBOX=500)
      INTEGER MAXNPEAKS
      PARAMETER (MAXNPEAKS=500)
      INTEGER MAXEXCL
      PARAMETER (MAXEXCL=100)
C     ..
C     ..
      INTEGER NLAUEG
      PARAMETER (NLAUEG = 44)
      INTEGER MAXLATT
      PARAMETER (MAXLATT = 10)
      INTEGER MAXIWARN
      PARAMETER (MAXIWARN = 210)
C&&*&& end_include  ../inc/parameter.f
      INTEGER NPARM
      PARAMETER (NPARM=200)
C     ..
C     .. Scalar Arguments ..
      INTEGER IFIRSTPACK,MODE
      LOGICAL FIRSTTIME,GENOPEN,NEWGENF,RPTFIRST
      CHARACTER CELLSTR*50
C     ..
C     .. Local Scalars ..
C     
C     EMBL
      LOGICAL VERS1,VERS2,VERS3
C     EMBL
C     
      REAL THRESHF,XMMF,XMMDB,XCC,TEMP,
     +     AXRMSLIM,XRMSLIM,RADEG,PHISTART,
     +     TEMPX,FXMAX,FYMAX,
     +     VXMIN,VXMAX,VYMAX,VRMAX,THETA,ROFFMAX,DSTRMX,DSTRES,
     +     DELR,T,RPLUS,THPLUS,XSEP,YSEP,XMAXRED,YMAXRED,RMAXRED,
     +     XDMID,YDMID,REDGE,RMINRED,TORSRS,
     +     COS2TH,DET,RSDMAX,PHI1,PHI2,RMINP,DSTMINP,
     +     OMEGADD,ETAD,DIVHD,DIVVD,PCCX,OVOLSCAL,RESEX1,RESEX2,
     +     VOLSCAL,CVOL,DEFVOL,XTRUE,YTRUE,
     +     XYEXX1,XYEXY1,XYEXX2,XYEXY2,XTEST,XYMAXST,
     +     DSTMAXST,THETAST,XPAD,THISRES,PHIREF,PHIWREF,PROFPHIS,
     +     PROFPHIE,PROFWIDTH,HIGHRESP,LOWRESP
      INTEGER I,II,IADD,ID2,IDFILM,IFIRST,
     +     IPCKID,IPNT,J,JJ,K,NC,NCASS,NCH,
     +     NPROF,NPRUN,NRX,NRY,NTOK2,NXS,NYS,
     +     ICOUNT,IFAIL,NCH2,NCH3,NCH4,
     +     NCHAR,N,
     +     KEYPX,IPACKF,IPACKL,IDELAMB,IMOSAIC,ICCY,
     +     NNDIR,
     +     IBULGE,MTZPRT,
     +     INODES,
     +     NP,IPART,
     +     IIPHI,IISIZE,IAUTO,IISTART,NSTRUNO,ISTAFLG,IENDFLG,IERR,
     +     NSEGRD,NRWORK,NSERRUN,NTIMES,
     +     ISTOP,NMULTI,IMGKWD,IEXTEN,
     +     IPROKWD,NPROCRUN,IUN1,IUN2,INSIZE,
     +     IIONE,NSEGOLD,INADD,INWIDTH,INSPEED,INMONO,
     +     ISTRT2,IANGLE2,IPRINT,IPKRAT,ISLOTB,ISLOTE,SNAUTO,INSEG,
     +     INFIX,DMODE,JSEG,FIRSTIMG,LASTIMG,NCIRPT,MCIRPT,IMODE,
     $     IFILESIZE,SIXSEP,SIYSEP,SISEP,EXCLUDERING,
     +     ILATT,ILATTSOL,INDENT,ISOL,LIPROCESS,IMOSFLMX,IMOSFLMY,
     +     IIH,IIK,IIL,IMOSFLMIMG,MODEPHI,IMOSFIRSTIMG,IMOSLASTIMG
      LOGICAL ASSIGN,EFILE,OTHERS,
     +     AVPR,EXPAND,PRINTL,RRSET,ARRSET,
     +     STOPRUN,FIXEDPR,
     +     FASTH,FASTV,RESET,FINDSPOT,ROTH,ROTV,ROTANTI,ROTCLOCK,
     +     ORGLR,ORGLL,ORGUR,ORGUL,SADDPART,SSUMPART,
     +     INERR,EXTRA,
     +     ANGLES,UNPACK,FRSTWARN,
     +     INPERR,TRAPERR,NULINE,
     +     READINLINE,SDPSINDEX,NOGO,
     $     DONESEG,IFSTRAT,DOSTRAT,CKNDIR,BADKEY,NEWREAD,HARRY,
     $     REMOVEALL,LDUM,BADPICK,BADSEGM,GOODSEG,NOGENOPEN,BADFIT,
     $     SIISEP,BADLINE,BADKEYWORD,NOTVERTICAL,BADFILE,FIRSTPHIP
      CHARACTER ABC*3,GRTYPE*3,KEY*4,STRL1*7,
     +     HLPMOS*400,STRL2*10,FIXSTR*80,
     +     STR*20,STRL3*22,COMFILE*200,LINE2*80,
     +     FWORK*200,SUBKEY*4,PROFFNR*100,
     +     PROFFNW*100,PRINTOP(10)*24,SCAN*4,TEMPCH*134,KEY2*4,
     +     IDXFILE*130,LATTYP*1,LINE80*80,STATION*4,
     +     DEBUGSTR*120,SIDENT*40,KEY8*8,TARFILE*130,STR2*100,
     $     LINE200*200,longoutline*4096,CALLEDFROM*80,INOUT*7,
     +     LINE400*400,THIS*30,CALLED*30,COMMENT*100,
     +     LONGLINE*1024
C     ..
C     .. Local Arrays ..
      REAL VALUE(NPARM),FIDXY(3,2),PX72(4),PX96(4),
     +     TDELPHI(3),VALUE2(NPARM),
     +     TCELL(6),SUMAT(3,3),SBMAT(3,3),
     +     WORK(3,3),UINV(3,3),WORK2(3,3),WORK3(3,3),SAMAT(3,3),
     +     SDELPHI(3),PHILEFT(40),PHISLEFT(40),
     +     TARPHI(3),TMAT(3,3),TARCELL(6),TJUNK(3,3),
     +     PKSUMS(6),CORNER(4,2)
      INTEGER IBEG(NPARM),IDEC(NPARM),IDPROF(MAXPAX),IEND(NPARM),
     +     ITYP(NPARM),NCHPR(10),LCLASS(6,8),IFIX(NRPAR),
     +     NFLEFT(40),NLLEFT(40),IBEG2(NPARM),
     +     IEND2(NPARM),ITYP2(NPARM),IDEC2(NPARM),
     +     ISERLEFT(40),MASK(MAXBOX),IXYCPT(2,200),SIRAS(5),LINEPOINT(4)
      LOGICAL UNFIX(6),EMPTY(MAXEXCL)
      CHARACTER DUMPSTR(11)*80,VERSTR(3)*6,SABC(6)*5,FIXSTRA(NRPAR)*10

C     ..
C     .. External Functions ..
      INTEGER LENSTR
      LOGICAL VAXVMS
      EXTERNAL LENSTR,VAXVMS
C     ..
C     .. External Subroutines ..
      EXTERNAL CCPUPC,MKEYNM,MOSHLP,MPARSER,QCLOSE,UGTENV,OPENODS,
     +     START,PSTART,SETMAT,CCPDPN,MINV33,MATMUL3,ROTMAT,
     +     RTOMISSET,CELLCHK,CELLFIX,MTZINI,ASUSET,
     +     MRDSYMM,GETSPOTS,WSPOT,TOREFIX,DPSIND,DPSCHOOSE,GETBLOCK,
     +     MXD_FLU,GETSEPRAS,STARTMTZ,IMGOUT,TARGMAT,
     +     WINDIO,ESTRES,ccpdat,utime,SETUPICE,SETDETECT,KWTEMP,
     +     KWIDEN,KWMATR,KWPROC,KWTEST,KWBEAM,KWAUTOIN,KWSIMU,REFLFIND
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,COS,MOD,SIN
C     ..
C     .. Common blocks ..
C&&*&& include  ../inc/amatch.f
C
C $Id: amatch.f,v 1.2 2014/05/30 15:47:41 andrew Exp $
C
C--- awk generated include file  amatch.h
C---- START of include file amatch.h
C
C
C     .. Scalars in Common /AMATCH/ ..
C     ..
      REAL RCONV,OVRLAP,RESOL1,RESOL2,AELIMIT,ARMSLIM,
     +     SECANGLE,DAMP,TRUECCOM,ETAMAX,AWRMSLIM,MOSNEW,
     $     SLOPE,INTERCEPT
      INTEGER NSTEP,NCYCA,NPASS,N2,N3,N4,N5,N6,N7,N8,N9,N10,NBEAM,
     +        MINREFAUTO
      LOGICAL MATCH,NOCENT,NOREFINE,RMOSAIC
C     ..
C     .. Common block /AMATCH/  ..
      COMMON /AMATCH/RCONV,OVRLAP,RESOL1,RESOL2,AELIMIT,ARMSLIM,
     +     SECANGLE,DAMP,TRUECCOM,ETAMAX,AWRMSLIM,MOSNEW,
     $     SLOPE,INTERCEPT,NSTEP,NCYCA,NPASS,N2,N3,N4,N5,N6,N7,
     $     N8,N9,N10,NBEAM,MATCH,NOCENT,NOREFINE,RMOSAIC,
     +     MINREFAUTO
C     ..
C
C
C&&*&& end_include  ../inc/amatch.f
C&&*&& include  ../inc/backg.f
C
C $Id: backg.f,v 1.2 2018/01/24 14:22:08 andrew Exp $
C
C--- awk generated include file  backg.h
C---- START of include file backg.h
C
C
C
C     .. Scalars in common block /BACKG/ ..
      REAL BGFRAC,BKGADD
      INTEGER NBGMIN
C     ..
C     .. Common Block /BACKG/ ..
      COMMON /BACKG/ BGFRAC,NBGMIN,BKGADD
C     ..
C
C
C&&*&& end_include  ../inc/backg.f
C&&*&& include  ../inc/cbfinc.f
C
C---- The following are used for transferring values between the C stuff
C     and Fortran for CBF images (see cbf_info.f for assignments)
C
       DOUBLE PRECISION CBF_DOUBLE(16)
       INTEGER*4 CBF_INT4(16)
       INTEGER CBF_INT(16)
       CHARACTER*24 CBF_CHAR(16)
       LOGICAL CBF_LOG(16)
       LOGICAL FULLCBF,CBFSETUP
       COMMON /CBF_PAR/ CBF_DOUBLE,CBF_INT4,CBF_INT,CBF_CHAR,CBF_LOG,
     +                  FULLCBF,CBFSETUP
C&&*&& end_include  ../inc/cbfinc.f
C&&*&& include  ../inc/ccondata.f
C
C $Id: ccondata.f,v 1.10 2014/12/11 12:13:22 harry Exp $
C
C--- awk generated include file  ccondata.h
C---- START of include file ccondata.h
C
C
C     .. Scalars in Common /CCONDATA/ ..
      CHARACTER GENFILE*200,CCFILE*200,ODEXT*8,WAXFN*134,
     +     GTITLE*80,IDENT*40,NEWMATNAM*200,MTZNAM*200,MTZDIR*200,
     $     PGNAME*10,
     +     SPGNAM*10,ODFILE*480,XODFILE*480,SEPCHAR*1,SNEWMATNAM*200,
     +     TEMPLSTART*100,TEMPLEND*100,TEMPLSAV*100,BITMAPTYPE*12,
     $     CCP4VERSION*10
C     ..
C     ..  Arrays in Common /CCONDATA/ ..
      CHARACTER FDISK(10)*400,INLINE(10000)*400
C     ..
C     .. Common block /CCONDATA/ ..
      COMMON /CCONDATA/GENFILE,CCFILE,ODEXT,WAXFN,GTITLE,IDENT,
     +     NEWMATNAM,MTZNAM,MTZDIR,PGNAME,SPGNAM,ODFILE,XODFILE,
     +     SEPCHAR,SNEWMATNAM,TEMPLSTART,TEMPLEND,
     $     TEMPLSAV,BITMAPTYPE,CCP4VERSION,FDISK,INLINE
C     ..
C
C
C&&*&& end_include  ../inc/ccondata.f
C&&*&& include  ../inc/cconst8.f
C
C $Id: cconst8.f,v 1.1 2002/05/02 10:46:42 harry Exp $
C
C--- awk generated include file  cconst8.h
C---- START of include file cconst8.h
C
C
C     .. Arrays in Common /CCONST8/ ..
      REAL CCOMA
      INTEGER CCXA,CCYA,CBARA
C     ..
C     .. Common block /CCONST8/ ..
      COMMON /CCONST8/CCOMA(8),CCXA(8),CCYA(8),CBARA(8)
C     ..
C
C
C
C&&*&& end_include  ../inc/cconst8.f
C&&*&& include  ../inc/cell.f
C
C $Id: cell.f,v 1.4 2012/11/21 14:41:22 andrew Exp $
C
C--- awk generated include file  cell.h
C---- START of include file cell.h
C
C     CELL cell dimensions (real space)
C     RCELL reciprocal cell parameters in dimensionless rlu
C     LCELL indicator for cell parameters to be refined
C     
C     LCELL(I) = -1  parameter free
C              =  0  parameter fixed
C              >  0, = J  parameter I constrained to = parameter J
C
C     .. Arrays in Common /CELLCOM/ ..
      REAL UMAT,BMAT,GMAT,CELL,RCELL,UMATCELL,SOFTCELL,VOLUME
      INTEGER LCELL,ICRYST,NUMSPG,NLAUE
C     ..
C     .. Common Block /CELLCOM/ ..
      COMMON /CELLCOM/UMAT(3,3),BMAT(3,3),GMAT(3,3),CELL(6),RCELL(6),
     $       UMATCELL(6),SOFTCELL,VOLUME,LCELL(6),ICRYST,NUMSPG,NLAUE
C     ..
C
C
C&&*&& end_include  ../inc/cell.f
C&&*&& include  ../inc/comarray1.f
C
C $Id: comarray1.f,v 1.2 2015/08/10 11:53:22 harry Exp $
C
C
C Common bits for the big array handling for post refinement
C
C---- for passing the summed partials around different routines
C
C---- COMMON BLOCK POSTPHI
C
      REAL PHIPOST,PARTITION
      COMMON /POSTPHI/ PHIPOST,PARTITION
C&&*&& end_include  ../inc/comarray1.f
C&&*&& include  ../inc/comarray2.f
C $Id: comarray2.f,v 1.2 2015/08/10 11:53:23 harry Exp $
C
C
C Common bits for the big array handling for post refinement
C
C---- for passing the summed partials around different routines
C
      INTEGER INTPOST(2),SIGPOST(2)
C
C---- COMMON BLOCK POSTVAL
C
      COMMON /POSTVAL/ INTPOST,SIGPOST
C&&*&& end_include  ../inc/comarray2.f
C&&*&& include  ../inc/condata.f
C
C $Id: condata.f,v 1.49 2021/02/11 11:30:48 andrew Exp $
C
C--- awk generated include file condata.h 
C---- START of include file condata.h 
C 
C
C     FIXSWAP     Forces reversal of normal byte-swapping choice.
C
C---- SUMPART True if two images are to be stored in memory at the
C               same time. This will be the case if post-refinement
C               (other than POSTHOC post-refinement for off-line
C               scanners) is being performed. Note that SUMPART is NOT
C               used to determine if partials should be added across
C               images (this is ADDPART)
C
C---- NEWPREF True if allowing partials over multiple images in post-
C             refinement.
C
C---- ADDPART True if partials are to be summed across images. This is
C               only appropriate for fixed origin (ie on-line)
C               scanners where the X-ray dose is identical for
C               each image.
C     Direct beam coordinates
C             These are read into XMM(3), YMM(3) where (3) is for A,B,C packs
C             Immediately on reading in, these values are transferred
C             into XCENMM, YCENMM, XCENMMIN(J,3) where J=1->MAXPAX.
C
C             If an IMAGE keyword has been given, XCENMM(1,1)=XMM(1) and
C             YCENMM(1)=YSCAL*YMM(1). If swung out, YCENMM(1,1) is updated.
C             If a BEAM keyword has not been given, XMM,YMM, XCENMM,YCENMM,
C             XCENMMIN are all set to the middle of the image.


C             Then, when reflecting input to logfile, it AGAIN sets up
C             XCENMM,YCENMM for all images using XMM,YMM, this time
C             using YSCAL to set up YCENMM, and updating YCENMM if the
C             detctor is swung out. (label 708), DO loop 710. If image
C             is "inverted", corrects XCENMM for all images.


C             XCEN,YCEN are set in "mosflm" to XCENF, YCENF plus
C             camera constants, and XCENF,YCENF are in turn set from
C             XCENMM,YCENMM. Thus XCENMM,YCENMM are what is actually
C             used to define the beam position.

C
C     LPREF   If TRUE, return to display after refinement over entire
C             image. Only active when integrating interactively. Set by
C             toggle button in Parameter window.
C

C     LOVERLAP  If TRUE, do not call the OVERLAP subroutine

C     ALLOUT  If TRUE, all reflections (whether measured or not, or
C             classified as BADSPOTS or not, will be written to the
C             output MTZ file.
C
C     NTDIG   Number of digits in image number, returned as zero if
C             there is an error in the template. Only used if TEMPLATE
C             keyword given.
C
C     TEMPLATE is TRUE if a TEMPLATE keyword has been given.
C
C     NOLP    If TRUE, do NOT apply the Lorentz Polarisation corrections
C
C     MULTIMTZ True if writing each "block" of images to a separate MTZ file
C
C     DISPMENU True if the run was started with a IMAGE keyword. Not to be
C              confused with WINOPEN which is true if the X-window display
C              is being used, but this is not necessarily starting with
C              a "IMAGE" keyword.
c     newgui - this means that we're talking through the new GUI.
c              except that we can't have it here because the prsummary 
c              subroutine, where we REALLY want it, clashes with this common
C
C     HACTIVEMASK,LACTIVEMASK,LACTIVE - only for Rigaku images which 
C     contain a compressed "active mask" (size=HACTIVEMASK bytes)
C     at the end of the image. Compression uses a modified run length encoding
C     at present (10.01.2003).
C
C     MINBATCH Minimum allowed batch number to be written to an open
C              MTZ file (MTZ library imposes that batch numbers must 
C              increase as file is written).
c
c     MERGELATT If true, prediction lists for multiple lattices will
C               be merged for output to iMosflm
C
C     NIMGSUMREF Number of images to be summed when determining spot
C                centroids for detector parameter refinement. Only
C                works for Pilatus detectors at present.
C     ..
C     .. Arrays in common /CONDATA/ ..
      REAL XCENMM,YCENMM,THFOIL,PHIBEGA,PHIENDA,RSYM,XMM,YMM,
     +     XCENMMIN,YCENMMIN
C     ..
C     .. Arrays in common /CONDATA/ ..
      INTEGER IDPACK,NFPACK,NFIRST,ICASSET
      INTEGER IPACK1A,IPACK2A,ISERAR,PICKER
      LOGICAL AVPROF,FILMPLOT,FORCEB,FORCEC,CONTROLM,CONTROLS
C
C     .. Scalars in common /CONDATA/ ..
C     ..
      REAL RMSLIM,WRMSLIM,VLIM,XMID,SEP,THICK,PTMIN,REFREJ,
     +     YSCALIN,FULLFRAC,WTIME,OVERLAPMAX
      INTEGER NPACK,NCYC,MAXNX,MAXNY,ITHRESHF,NSIG,IRFMIN,
     +        IRFINC,IXSHIFT,IYSHIFT,LIMIT,NBLOCK,INOGEN,NFGEN,MINREF,
     +        ISERADD,NSEG,NSYM,NSYMP,NRUN,NSERTOT,IWAVE,ISWUNG,
     $        IBEAMMU,ITILT,ITWIST,IPIX,IRMINXY,ICUT,IPRCUT,ILIMIT
      INTEGER IMAT,ICELL,IUMAT,LSYMM,IDIST,IBEAM,ISCAN,NLINE,NRLINE,
     +     NTLINE,NTOTLINE,ISTARTP,IRAST,INEWMAT,ISEP,NSAVELINE,IPTIME,
     +     ISTRT,IANGLE,IBACKS,IPOLAR,IDIVH,IDIVV,IGAIN,IBLOCK,
     +     IHKLOUT,IGENF,INRES,NTDIG,DELAY,IMULTI,IMISS,IMISSMAT,
     +     IIDENT,HACTIVEMASK,ITWOTHETA,MINBATCH,ILABVIEW,INVVIEW,
     $     IINVVIEW,IHKLDIR,IBGSIG,NRLEFT,WAIT,IPARTITION,ISIGSET,
     +     NSIGSAVE,ITOL,NIMGSUMREF,IBGFRAC,THISIMG,CHUNK,ICHUNK,
     $     IPTMINSET,IYSCAL
      LOGICAL FINDCC,READCC,NOFID,INTERPOL,CONVOL,NOMEAS,
     +     USEPAR,USEOVR,SUMPART,USEBOX,POSTREF,ADDPART,
     +     FIXED,POWDER,RWEIGHT,NOREF,MULTISEG,FIXSWAP,POSTHOC,
     +     LPREF,LPINTG,DISPMENU,PRMODE,PRCELL,TEMPLATE,
     +     LOVERLAP,ALLOUT,NOLP,NEWPREF,NOBACK,MOSEST,BLKEST,NUSPOT,
     $     NUBACK,NUFIND,HARVESTREADY,HEADINFO,MULTIMTZ,ANDREW, 
     +     spostref,newgui,LACTIVEMASK,USERSPOT,HEADONLY,PHIONLY,
     $     NO_XML,TESTSOCK,PICKLOG,EXPLICIT_CELL,LSIMULATE,
     $     LONGERNAME,TIMER,MERGELATT,NULLPIXTEST,KEEPNULLPIX,CCHUNK,
     +     INVXFORCEREAD
C     ..
C real arrays, then scalars, integer arrays, scalars, logicals
      COMMON /CONDATA/
     $     XCENMM(MAXPAX,3),
     +     YCENMM(MAXPAX,3), THFOIL(3),PHIBEGA(MAXPAX),
     $     PHIENDA(MAXPAX),RSYM(4,4,96),XMM(3),YMM(3),
     $     XCENMMIN(MAXPAX),YCENMMIN(MAXPAX),
     $     RMSLIM,WRMSLIM,VLIM,XMID,SEP,THICK,PTMIN,REFREJ,
     +     YSCALIN,FULLFRAC,WTIME,OVERLAPMAX,
     $     IDPACK(MAXPAX),NFPACK(MAXPAX),NFIRST(MAXPAX),
     +     ICASSET(MAXPAX),IPACK1A(50),
     +     IPACK2A(50),ISERAR(50),PICKER(4),
     $     AVPROF(MAXPAX),FILMPLOT(MAXPAX),FORCEB(MAXPAX),
     +     FORCEC(MAXPAX),
     $     NPACK,NCYC,MAXNX,MAXNY,ITHRESHF,NSIG,IRFMIN,
     +     IRFINC,IXSHIFT,IYSHIFT,LIMIT,NBLOCK,INOGEN,NFGEN,MINREF,
     +     ISERADD,NSEG,NSYM,NSYMP,NRUN,NSERTOT,IWAVE,ISWUNG,
     $     IBEAMMU,ITILT,ITWIST,IPIX,IRMINXY,ICUT,IPRCUT,ILIMIT,
     $     IMAT,ICELL,IUMAT,LSYMM,IDIST,IBEAM,ISCAN,NLINE,NRLINE,
     +     NTLINE,NTOTLINE,ISTARTP,IRAST,INEWMAT,ISEP,NSAVELINE,
     +     IPTIME,ISTRT,IANGLE,IBACKS,IPOLAR,IDIVH,IDIVV,IGAIN,IBLOCK,
     +     IHKLOUT,IGENF,INRES,NTDIG,DELAY,IMULTI,IMISS,IMISSMAT,
     +     IIDENT,HACTIVEMASK,ITWOTHETA,MINBATCH,ILABVIEW,INVVIEW,
     $     IINVVIEW,IHKLDIR,IBGSIG,NRLEFT,WAIT,IYSCAL,
     +     FINDCC,READCC,NOFID,
     +     INTERPOL,CONVOL,NOMEAS,USEPAR,USEOVR,SUMPART,USEBOX,
     +     POSTREF,spostref, ADDPART,FIXED,POWDER,RWEIGHT,NOREF,
     +     MULTISEG,FIXSWAP,POSTHOC,LPREF,LPINTG,DISPMENU,PRMODE,
     $     PRCELL,TEMPLATE,
     +     LOVERLAP,ALLOUT,NOLP,NEWPREF,NOBACK,MOSEST,BLKEST,NUSPOT,
     $     NUBACK,NUFIND,HARVESTREADY,HEADINFO,MULTIMTZ,ANDREW,
     $     newgui,LACTIVEMASK,USERSPOT,HEADONLY,PHIONLY,NO_XML,
     $     TESTSOCK,PICKLOG,EXPLICIT_CELL,LSIMULATE,
     $     LONGERNAME,TIMER,CONTROLM,CONTROLS,IPARTITION,ISIGSET,
     +     NSIGSAVE,ITOL,MERGELATT,NIMGSUMREF,NULLPIXTEST,
     +     KEEPNULLPIX,IBGFRAC,THISIMG,CHUNK,ICHUNK,CCHUNK,IPTMINSET,
     +     INVXFORCEREAD


C&&*&& end_include  ../inc/condata.f
C&&*&& include  ../inc/debug.f
C
C $Id: debug.f,v 1.8 2021/06/03 09:41:34 andrew Exp $
C
C--- awk generated include file  debug.h
C---- START of include file debug.h
C
C
C
C     .. Arrays in common /DEBUG/ ..
      REAL XWARN
      INTEGER NDEBUG,IWARN,IARRDUMP
      LOGICAL DEBUG,LPRINT,DUMP,WARN
C
C     .. Scalars in common /DEBUG/ ..
      REAL BGRLIM
      INTEGER NDUMP,IMINDUMP,IMAXDUMP,MXDUMP,IERRLEVEL
      LOGICAL SPOT
      CHARACTER MESSAGE*512
C     
C     ..
C     .. Common Block /DEBUG/..
      COMMON /DEBUG/XWARN(MAXIWARN,100),BGRLIM,NDEBUG(200),
     +       IWARN(MAXIWARN,100),
     $       NDUMP,IMINDUMP,IMAXDUMP,MXDUMP,DEBUG(200),LPRINT(20),
     +       DUMP(30),WARN(100),SPOT,IARRDUMP(30),IERRLEVEL

      COMMON /CDEBUG/ MESSAGE
C     ..
C
C&&*&& end_include  ../inc/debug.f
C&&*&& include  ../inc/detector.f
C
C $Id: detector.f,v 1.3 2007/11/13 17:21:01 harry Exp $
C
C
C---- START of include file detector.f
C
      REAL RMINIP,RMAXIP,XMAXIP,YMAXIP,XLIMIT,RSCANIP
      COMMON /DETECTOR/RMINIP,RMAXIP,XMAXIP,YMAXIP,XLIMIT,
     $     RSCANIP
C&&*&& end_include  ../inc/detector.f
C&&*&& include  ../inc/dpsindex.f
C
C $Id: dpsindex.f,v 1.27 2014/07/17 15:50:09 tnakane Exp $ 
C
C---- START of include file dpsindex.h
C
C----- Add parameters to control DPS indexing and debug output
C      AGWL
C
C     INDNOREF   Logical that controls refinement of cell parameters in
C                autoindexing. Set FALSE by default but can be toggled
C                by keywords AUTOINDEX NOREFINE. Also if there is an
C                error retrun from CELREF it is set TRUE.
C
C     ..
C     .. arrays in common DPSDATA
c     TAKANORI: added PRIORCELL
      REAL DPSXLO(10),DPSXHI(10),CL_OM(3,3,45),PRIORCELL(6),
     $     REFLERR(NSPOTS),
     $     XORIG(NSPOTS),YORIG(NSPOTS),IORIG(NSPOTS),ISDORIG(NSPOTS),
     $     X(NSPOTS),Y(NSPOTS),Z(NSPOTS),X_SPOT(NSPOTS),
     $     Y_SPOT(NSPOTS),PHIMEAN(NSPOTS),
     $     DELSQXY(NSPOTS),deltax(nspots),deltay(nspots),
     $     phiobs(nspots),phicalc(nspots)
      INTEGER USESPOT(NSPOTS)
      INTEGER*2 IH(NSPOTS),IK(NSPOTS),IL(NSPOTS),IRESORDER(44,maxlatt)
      LOGICAL  DUMPFFT(3)
C     ..
C     .. scalars in common DPSDATA
      REAL DISTANCE,XTD,WAVELENGTH,XBEAM,YBEAM,PHIAV,DMAX,DMIN,
     $     ORGX,ORGY,XCOR,YCOR,SDCUTOFF,RMSTOLER,
     +     STEPSCAL,DMAXSCAL,CONVRG,THETA1D,THETA2D,SINTHMIN,
     +     VOLMIN,HKLDEV,HKLMIN,MAXDEV
      INTEGER NUMVECT,NUMORMAT,NUMRLP,NUMSPTDPS,NPKFFTMIN,NMULTFFT,
     +     NUMVECTREF,MAXCYC,INICYC,UNINDEXED,
     $     RINGTEST,DPSEXCL,COUNT,SOLN,CENTERING
      LOGICAL DUMPVEC,DUMPREF,DUMPGEN,THETADMP,
     +     DUMPCHOOSE,DUMP3,DUMP4,DUMP5,INDNOREF,DPSINDEX,PREREF,
     $     TRIREF,LSOL,SIGMOID,DODEN,
     +     RMSONLY,MULTILATT,RECYCLE,PRUNE,CHOOSEHKL,
     +     RFIXCELL,RFIXDIST
      CHARACTER*2 BTSTR
      COMMON /DPSDATA/ DPSXLO,DPSXHI,CL_OM,PRIORCELL,CENTERING,
     +     REFLERR,XORIG,YORIG,IORIG,
     $     ISDORIG,X,Y,Z,X_SPOT,Y_SPOT,PHIMEAN,DELSQXY,
     $     deltax,deltay,phiobs,phicalc,usespot,
     $     IH,IK,IL,IRESORDER,DUMPFFT,
     $     DISTANCE,XTD,WAVELENGTH,XBEAM,YBEAM,PHIAV,DMAX,
     $     DMIN,ORGX,ORGY,XCOR,YCOR,SDCUTOFF,RMSTOLER,
     +     STEPSCAL,DMAXSCAL,CONVRG,THETA1D,THETA2D,SINTHMIN,
     +     VOLMIN,HKLDEV,HKLMIN,MAXDEV,NUMVECT,NUMORMAT,NUMRLP,
     $     NUMSPTDPS,NPKFFTMIN,
     +     NMULTFFT,NUMVECTREF,MAXCYC,INICYC,UNINDEXED,
     $     RINGTEST,DPSEXCL,COUNT,SOLN,DUMPVEC,DUMPREF,DUMPGEN,
     $     THETADMP,
     +     DUMPCHOOSE,DUMP3,DUMP4,DUMP5,INDNOREF,DPSINDEX,PREREF,
     $     TRIREF,LSOL,SIGMOID,DODEN,
     +     RMSONLY,MULTILATT,RECYCLE,PRUNE,CHOOSEHKL,RFIXCELL,RFIXDIST
      COMMON /DPSDATAC/ BTSTR
C
C
C COMMON BLOCK USED IN PERMUTING THE CELL
C
      REAL KCELL(6)
      LOGICAL KWDCELL
      COMMON /PERMUTE/ KCELL,KWDCELL
C
C
C&&*&& end_include  ../inc/dpsindex.f
C&&*&& include  ../inc/dsplyc.f
C
C $Id: dsplyc.f,v 1.7 2015/05/27 10:36:13 harry Exp $
C
C--- awk generated include file  dsplyc.h
C---- START of include file dsplyc.h
C
C*******************************************************************
C
C  COMMON  /DSPLYC/
C
C       IMGLOW, IMGHI   low & high values of 16-bit image for scaling
C                       integer*2 to byte: IMGLOW maps to  0; 
C                       IMGHI to maximum. Note that these are not
C                       necessarily the actual limits of the data
C       IMGMAX          maximum value recorded on the image. Necessary for 
C                       iMosflm.
C       IMGMIN          Minimum value recorded on image - gaps on PADs, 
C                       necessary for iMosflm
C       JDSPWD          .LT. 0  before image window has been created
C                       = +-1 for image display that can be panned
C                       = +-2 for non-interactive image display
C       MAXDEN          highest level in colour table to fill up to
C                       must be less than ~240 - number of overlay colours
C       LDSPSG          if .true., treat image as signed, ie after dark
C                          subtraction
C                       if .false., treat image as unsigned
C       NZOOM           zoom factor for image, = 0 if no zoom
C       JYZOOM, JZZOOM  1st pixel in zoomed image
C
C----   WINOPEN Flag for whether or not window is open. Do not
C               confuse with DISPMENU (/CONDATA/)which is true if the run was
C               started with a IMAGE keyword.
C
C
C       CDSPTL          banner title
C
      REAL FTOTAL,AVERAG,IMAGERMS
      INTEGER IMGLOW, IMGHI, IMGMAX,IMGMIN,JDSPWD, MAXDEN,
     $     NZOOM, JYZOOM, JZZOOM
      INTEGER PIXEL_INTENSITIES(514)
      LOGICAL LDSPSG,WINOPEN,LPAUSE,LHELP
      COMMON /DSPLYC/ FTOTAL,AVERAG,IMAGERMS,IMGLOW, IMGHI, IMGMAX,
     *     IMGMIN,JDSPWD,MAXDEN, NZOOM, JYZOOM, JZZOOM, 
     $     PIXEL_INTENSITIES,LDSPSG,WINOPEN,LPAUSE,LHELP
C
      CHARACTER  CDSPTL*200,SHORTFN*200
      COMMON /DSPLCC/  CDSPTL,SHORTFN
C                                                           
C
C*******************************************************************


C&&*&& end_include  ../inc/dsplyc.f
C&&*&& include  ../inc/extras.f
C
C $Id: extras.f,v 1.2 2003/04/04 14:22:01 harry Exp $
C
C--- awk generated include file  extras.h
C---- START of include file extras.h
C
C
C     .. Scalars in common /EXTRAS/ ..
      INTEGER JUMPAX,NPACKS
      LOGICAL ISKIP,AFILM,BFILM,CFILM,STARTA,STARTB,STARTC,NOSTOP
C     ..
C     .. Common block /EXTRAS/ ..
      COMMON /EXTRAS/JUMPAX,NPACKS,ISKIP,AFILM,BFILM,CFILM,STARTA,
     +       STARTB,STARTC,NOSTOP
C     ..
C
C
C&&*&& end_include  ../inc/extras.f
C&&*&& include  ../inc/fid.f
C
C $Id: fid.f,v 1.7 2016/04/22 14:23:39 andrew Exp $
C
C--- awk generated include file  fid.h
C---- START of include file fid.h
C
C     XCENF,YCENF  Coordinates (in 10 micron units) of the direct beam
C                  position relative to an origin at the position of the
C                  first pixel in the digitised image.(The SCANNER
C                  coordinate frame).
C                  For IP data this comes from the defined
C                  direct beam coordinates, for film data it is the
C                  midpoint of fiducials 1 and 3.
C
C     CCX,CCY      The difference (in 10 micron units) between the
C                  refined position of the direct beam (XCEN,YCEN) and
C                  the ideal direct beam coordinates (XCENF,YCENF),
C                  in the SCANNER coordinate frame.
C                  CCX,CCY are updated in RDIST. They are in "pixels"
C                  rather than mm (but expressed in 10micron units).
C
C     DTOFD        The distance (in 10 micron units) from the crystal
C                  to the detector along a normal to the detector.
C                  For flat, unswung detectors, or swung on a 2 theta arm,
C                  this is the same as CTOFD and XTOFD. It differs for
C                  Vee shaped cassettes. Assigned in START and never changed.
C                  Only actually used in RMAXR for calculating box sizes.
C
C     IOMEGAF      Whether OMEGAF has been set by user, by program, or 
C                  by reading header
C     .. Arrays in common /FID/ ..
      REAL CCOMABC
      INTEGER FSPOS,CCXABC,CCYABC
C     ..
C     .. Scalars in common /FID/ ..
      REAL OMEGAF,OMEGAFD,CCOM,DTOFD,XCENF,YCENF,OMEGAREV
      INTEGER MM,MMDB,NFID,CCX,CCY,IYOFF,LABVIEW,IOMEGAF
      LOGICAL UPDATEOMEGA
C     ..
C     .. Common block /FID/ ..
      COMMON /FID/CCOMABC(3),OMEGAF,OMEGAFD,CCOM,DTOFD,XCENF,YCENF,
     $     FSPOS(4,2),
     +     CCXABC(3),CCYABC(3),MM,MMDB,NFID,CCX,CCY,IYOFF,LABVIEW,
     $     IOMEGAF,OMEGAREV,UPDATEOMEGA
C     ..
C
C     .. Scalars in common /FIDC/ ..
C
      CHARACTER SPOTFILEORIGIN*2,SPOTFILEAXISORDER*2,SWINGAXIS*2
      COMMON /FIDC/ SPOTFILEORIGIN,SPOTFILEAXISORDER,SWINGAXIS
C&&*&& end_include  ../inc/fid.f
C&&*&& include  ../inc/graphics.f
C
C $Id: graphics.f,v 1.1 2002/05/02 10:46:48 harry Exp $
C
C--- awk generated include file  graphics.h
C---- START of include file graphics.h
C
C
C     .. Scalars in common /GRAPHICS/ ..
      REAL GRFACT,DISPLAY
      INTEGER NGR,NGX,NGY,NHX,NHY,NLI
C     ..
C     .. Common block /GRAPHICS/ ..
      COMMON /GRAPHICS/GRFACT,DISPLAY,NGR,NGX,NGY,NHX,NHY,NLI
C     ..
C
C
C&&*&& end_include  ../inc/graphics.f
C&&*&& include  ../inc/header.f
C
C $Id: header.f,v 1.16 2021/03/04 11:41:58 andrew Exp $
C
C--- awk generated include file  header.h
C---- START of include file header.h
C     .. Scalars in common block /HEADER/
      REAL HDIST,HWAVE,HRAST,HTOR,HPHIS,HPHIE,HPHIR,HOMEGAS,
     $     HOMEGAE,HTWOTHETA,HKAPPA,HCHI,
     $     HBEAMX,HBEAMY,HCCOMEGA,HTWIST,HTILT
      INTEGER NHEAD,NTAIL,HNULLPIX,HYEAR,HDATE,HCOMPR,
     $     HBIAS,SUMBYTES,HUNDER,HDATEC,HCUTOFF,HPRCUTOFF,HHEAD
      LOGICAL USEHDR,USETAIL,USEDIST,USEWAVE,USEPHI,HDRSIZE,TESTDATE
C     ..
C     .. Arrays in common block /HEADER/
      INTEGER IHEAD,HOVER
C     ..
C     .. Common Block /HEADER/
      COMMON /HEADER/ HDIST,HWAVE,HRAST,HTOR,HPHIS,HPHIE,HPHIR,
     $     HOMEGAS,HOMEGAE,HTWOTHETA,HKAPPA,HCHI,
     $     HBEAMX,HBEAMY,HCCOMEGA,HTWIST,HTILT,
     $     IHEAD(MAXHEAD),
     $     HOVER(3),NHEAD,NTAIL,HNULLPIX,HYEAR,HDATE,
     +     HCOMPR,HBIAS,SUMBYTES,HUNDER,HDATEC,HCUTOFF,HPRCUTOFF,
     $     HHEAD,USEHDR,USETAIL,
     +     USEDIST,USEWAVE,USEPHI,HDRSIZE,TESTDATE
C     ..
      CHARACTER HCOMPRESS*3,HMONTH*3,HOSCAXIS*6,HSERIAL*12
      COMMON /CHEADER/ HCOMPRESS,HMONTH,HOSCAXIS,HSERIAL
C&&*&& end_include  ../inc/header.f
C&&*&& include  ../inc/ifpwdr.f
C $Id: ifpwdr.f,v 1.10 2019/02/28 15:20:08 andrew Exp $
C---- IFPWDR block - variables only used in IFPWDR.F which were local
C     variables in CONTROL.F
C     ..
C     .. Array Arguments ..
      REAL SAVECELL(6)
      INTEGER IDAUTO(MAXIMG),IDIMG(MAXIMG)
      LOGICAL PHISET(MAXIMG)
      CHARACTER IDENTAUTO(MAXIMG)*40,TEMPLAUTO(MAXIMG)*100
C     ..
C     .. Scalar Arguments ..
      REAL OMEGA0,RSCAN,RMINXINP,RSCANRED,RAD,MAXCELL,MINCELL
      INTEGER ICCX,NSER,NAUTO,ID1,NIMAGES,ICHECK,
     $     ITINS,ICOMM,ID,IFLAG,IPACK,IPIXY,ITOR,
     $     IXOFFSET,IYOFFSET,MODEGSR,MODEOP,MODESP,
     $     MOSIMAG,NFIRSTF,NIMAGP,NSOL,
     $     NTOK
      LOGICAL AUTOINDX,BOXOPEN,CCXRESET,CELLKEEP,DISPSET,DONERUN,
     $     DPSDONE,COMREAD,
     $     FORCEREAD,PACK,FINE,LPRNT,ROTATED,
     $     SAVIND,SYMMIN
      CHARACTER SCANNER*4,KEY6*6,LINE*400,SITE*4,XMLLINE*1024
      COMMON /IFPOWDER/SAVECELL,IDAUTO,IDIMG,PHISET,
     $     OMEGA0,RSCAN,RMINXINP,RSCANRED,RAD,MAXCELL,
     $     MINCELL,ICCX,NSER,NAUTO,ID1,NIMAGES,ICHECK,ITINS,ICOMM,ID
     $     ,IFLAG,IPACK,IPIXY,ITOR,IXOFFSET,IYOFFSET,
     $     MODEGSR,MODEOP,MODESP,MOSIMAG,NFIRSTF,NIMAGP,NSOL,
     $     NTOK,AUTOINDX,BOXOPEN,CCXRESET,CELLKEEP,DISPSET,DONERUN,
     $     DPSDONE,COMREAD,FORCEREAD,PACK,FINE,LPRNT,
     $     ROTATED,SAVIND,SYMMIN
      COMMON /IFPOWDERC/ IDENTAUTO,TEMPLAUTO,SCANNER,KEY6,LINE,
     +     SITE,XMLLINE
C&&*&& end_include  ../inc/ifpwdr.f
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
C&&*&& include  ../inc/ioomtz.f
C
C $Id: ioomtz.f,v 1.5 2011/06/27 12:46:46 andrew Exp $
C
C--- awk generated include file  ioomtz.h
C---- START of include file ioomtz.h
C
C
C     .. Scalars in common block /IOOMTZ/ ..
      INTEGER MTZOUT,BESTHKL
      LOGICAL MTZOPEN,BESTOPEN,NEWMTZ
      CHARACTER*100 BSTHKL
C     ..
C     .. Common block /IOOMTZ/ ..
      COMMON /IOOMTZ/ MTZOUT,BESTHKL,MTZOPEN,BESTOPEN,BSTHKL,NEWMTZ
C&&*&& end_include  ../inc/ioomtz.f
C&&*&& include  ../inc/ioosum.f
C
C $Id: ioosum.f,v 1.3 2010/10/14 12:47:02 harry Exp $
C
C--- awk generated include file  ioosum.h
C---- START of include file ioosum.h
C
C     .. Scalars in common block /IOOSUM/
      INTEGER NLSUM1,NLSUM2,NSUMSTART1,NSUMSTART2
C
C     .. Arrays in common block /IOOSUMC/ ..
      CHARACTER*200 LINESUM1(MAXPAX),LINESUM2(MAXPAX)
      CHARACTER*101  IOLINE(100)
C     ..
C     .. Common block /IOOSUM/ ..
      COMMON /IOOSUM/ NLSUM1,NLSUM2,NSUMSTART1,NSUMSTART2
C     ..
C     .. Common block /IOOSUMC/ ..
      COMMON /IOOSUMC/ LINESUM1,LINESUM2,IOLINE
C&&*&& end_include  ../inc/ioosum.f
C&&*&& include  ../inc/lattices.f
C
C---- START of include file lattices.f
C
C----- Parameters involved in processing of multiple lattices
C      AGWL
C 
C      NLATTICE      Total number of lattices (for which an orientation 
C                    is available)
C 
C      ILATTICE      Number of the lattice being integrated
C
C      IOVWORK       Number of additional lattices that should be taken
C                    into account when integrating current lattice
C
C      IOVLATTICE    Array storing the lattice numbers of the IOVWORK
C                    lattices that are being taken into account
C
C      ACTIVELATT    Array storing flags for each lattice, TRUE if this
C                    lattice is to be considered, FALSE if it has been
C                    deleted by the user in iMosflm
C
C      CRYSTALA      AMAT matrices for the different lattices
C
C      NUMSOL        Number of refined solutions for each lattice.
C
C      CRYSTALUBREF  Refined AMATs for different solutions and different 
C                    lattices. Note that the solution number (3rd element)
C                    is the number AFTER sorting on penalty.
C
C      UBREFSET      Logical set TRUE if DPSCHOOSE has been called so that
C                    matrices are stored in CRYSTALUBREF
C
C      CRYSTALCELL   Refined cell parameters AMATs for different solutions 
C                    and different lattices.Note that the solution number 
C                    (2nd element) is the number AFTER sorting on penalty.
C
C      NPRLATT       The lattice number being refined in a cell refinement
C                    specified by LATTICE subkeyword on POSTREF keyword and
C                    sent to iMosflm in CELL_REFINE_RESPONSE XML.
C
C      TOPSOLAR      The "best" (suggested) solution for each lattice based 
C                    only on the penalty, set in DPSIND.
C
C      BESTSOLAR     The "best" (suggested) solution for each lattice. 
C                    Initially set based on the penalty and rmsd (in 
C                    DPSCHOOSE) but can be changed subsequently by the user 
C                    via the iMosflm GUI.
C
C      BRAVLATT      The 2 character Bravais lattice type (aP, P2 etc) for
C                    each solution for each lattice. Note that the solution
C                    number here is PRIOR to sorting on penalty. The array
C                    IRESORDER in /DPSINDEX/ gives the pointer into BRAVLATT
C                    for a particlular solution NSOL (after sorting)
C                    
C      BRAVAIS       Is BRAVLATT but now after sorting on penalty
C
C      SINGLETON     Treat as multi-lattice in integration, but only write
C                    singletons to the output MTZ file (DEBUG tool)  
C
C      SUGGCELL      Cell suggested by Mosflm for each lattice
C
C      SUGGLATT      Bravais lattice for SUGGCELL
C     ..
C     .. arrays in common lattices
      REAL CRYSTALA(3,3,MAXLATT),
     +     CRYSTALU(3,3,MAXLATT),CRYSTALC(6,MAXLATT),
     +     CRYSTALUB(3,3,44,MAXLATT),THISCELL(6,44,MAXLATT),
     +     CRYSTALUBREF(3,3,44,MAXLATT),CRYSTALCELL(6,44,MAXLATT),
     $     SUGGCELL(6,MAXLATT)
      INTEGER IOVLATTICE(MAXLATT),NUMSOL(MAXLATT),TOPSOLAR(MAXLATT),
     +        BESTSOLAR(MAXLATT)
      LOGICAL ACTIVELATT(MAXLATT)
      CHARACTER*2 BravLatt(44,MAXLATT),BRAVAIS(44,MAXLATT),
     $     SUGGLATT(MAXLATT)
      
C     ..
C     .. scalars in common LATTICES
      INTEGER ILATTICE,IOVWORK,NLATTICE,NPRLATT
      LOGICAL UBREFSET,SINGLETON
C
      COMMON /LATTICES/ CRYSTALA,CRYSTALU,CRYSTALC,CRYSTALUB,
     +       THISCELL,CRYSTALUBREF,NLATTICE,ILATTICE,IOVLATTICE,
     +       IOVWORK,NUMSOL,CRYSTALCELL,SUGGCELL,NPRLATT,TOPSOLAR,
     +       BESTSOLAR,UBREFSET,ACTIVELATT,SINGLETON
      COMMON /LATTICESC/ BravLatt,BRAVAIS,SUGGLATT
C&&*&& end_include  ../inc/lattices.f
C&&*&& include  ../inc/mcs.f
C $Id: mcs.f,v 1.1 2002/05/02 10:46:56 harry Exp $
C
C
C--- awk generated include file  mcs.h
C---- START of include file mcs.h
C
C
C     .. Scalars in common block /MCS/ ..
      REAL BASEOD,G1OD,CURV
      INTEGER FILM,FILMS,CTOFD,PCKIDX,TOSPT,N1OD,XSCMIN,XSCMAX
      LOGICAL VEE,VALONGX
C     ..
C     .. Common block /MCS/ ..
      COMMON /MCS/BASEOD,G1OD,CURV,FILM,FILMS,CTOFD,PCKIDX,TOSPT,
     $            N1OD,XSCMIN,XSCMAX,VEE,VALONGX
C     ..
C     CTOFD.... This was the crystal to detector distance (10 micron units)
C               but it has now been replaced by XTOFD in common block /XY/
C
C&&*&& end_include  ../inc/mcs.f
C&&*&& include  ../inc/misc.f
C
C $Id: misc.f,v 1.3 2012/03/23 12:02:50 harry Exp $
C
C--- awk generated include file  misc.h
C---- START of include file misc.h
C
C
C
C     .. Scalars in common /MISC/ ..
      DOUBLE PRECISION STARTTIME,ENDTIME
      REAL PHIBEG,PHIEND,PHIRNG,RMIN,RMAX,WAVE
      INTEGER IPACKID,MININT,IERRFLG
C     ..
C     .. Arrays in common /MISC/ ..
      REAL DELPHI,RESANI
      INTEGER IAX
C     ..
C     .. LOGICAL
      LOGICAL ANITES

C     ..
C     .. Common Block /MISC/ ..
      COMMON /MISC/DELPHI(3),RESANI(3),STARTTIME,ENDTIME,
     $     PHIBEG,PHIEND,PHIRNG,RMIN,RMAX,
     $     WAVE,IAX(3),IPACKID,MININT,IERRFLG,ANITES
C     ..
C
C
C&&*&& end_include  ../inc/misc.f
C&&*&& include  ../inc/modify.f
C
C $Id: modify.f,v 1.1 2002/05/02 10:46:58 harry Exp $
C
C--- awk generated include file  modify.h
C---- START of include file modify.h
C
C
C     .. Arrays in common block /MODIFY/ ..
      LOGICAL MODS
C     ..
C     .. Common Block /MODIFY/ ..
      COMMON /MODIFY/MODS(30)
C     ..
C
C
C&&*&& end_include  ../inc/modify.f
C&&*&& include  ../inc/mxdinc.f
C
C $Id: mxdinc.f,v 1.1 2002/05/02 10:47:00 harry Exp $
C
C--- awk generated include file  mxdinc.h
C---- START of include file mxdinc.h
C
C
c common block for mxd routines: various parameters for XDL routines
c
c     base_width    width of base frame
c     base_height   height of base frame
c     men_x,men_y   position of menu area
c     par_x,par_y   position of parameter area
C     

C     IORDER    (R)   Order of the data in the input image 
c                with respect to two local axes ax1, ax2 e.g.
c                (xf, yf) as a number from 1 to 8.
c                    1   +ax1 slow   +ax2 fast    (+xf, +yf)
c                    2   +ax1 slow   -ax2 fast    (+xf, -yf)
c                    3   -ax1 slow   +ax2 fast    (-xf, +yf)
c                    4   -ax1 slow   -ax2 fast    (-xf, -yf)
c                    5   +ax2 slow   +ax1 fast    (+yf, +xf)
c                    6   +ax2 slow   -ax1 fast    (+yf, -xf)
c                    7   -ax2 slow   +ax1 fast    (-yf, +xf)
c                    8   -ax2 slow   -ax1 fast    (-yf, -xf)
c
c
c     JORDER    (R)   Display order with respect to the two
c                local axes (1 to 8) along the X-windows
c                axes X horizontal (left to right), Y vertical 
c                (top to bottom ) with origin at top left
c                1  +ax1 X (horiz)   +ax2 Y (vert)   (+xf, +yf)
c                2  +ax1 X (horiz)   -ax2 Y (vert)   (+xf, -yf)
c                3  -ax1 X (horiz)   +ax2 Y (vert)   (-xf, +yf)
c                4  -ax1 X (horiz)   -ax2 Y (vert)   (-xf, -yf)
c                5  +ax2 X (horiz)   +ax1 Y (vert)   (+yf, +xf)
c                6  +ax2 X (horiz)   -ax1 Y (vert)   (+yf, -xf)
c                7  -ax2 X (horiz)   +ax1 Y (vert)   (-yf, +xf)
c                8  -ax2 X (horiz)   -ax1 Y (vert)   (-yf, -xf)
c                    
c                    (2 was standard for Laue programs)
c
c     nxp_cmp       number of pixels compressed to 1 in horizontal
c                   display direction
c     nyp_cmp       number of pixels compressed to 1 in vertical
c                   display direction
c     nxdpx         number of horizontal pixels in displayed image
c     nydpx         number of vertical pixels in displayed image
c     img_horiz     image direction for horizontal display (=1 fast, 2 slow)
c     img_vert      image direction for vertical display   (=1 fast, 2 slow)
c     img_map(4)    mapping of image to display
c                    imgmap(1)   1st point on fast axis (Yms)
c                    imgmap(2)   1st point on slow axis (Zms)
c                    imgmap(3)   increment on fast axis
c                    imgmap(4)   increment on slow
c     img_x,img_y   x, y position of image
c     men_x, men_y  x, y position of menu
c     par_x, par_y  x, y position of parameter table
c     io_x, io_y    x, y position of IO area
c     not_x, not_y    x, y position of notice area (parameter table
c                     with no values)
c     pmn_x, pmn_y  x, y position of popup abort menu
c     ppb_x, ppb_y  x, y position of popup progress bar
c     busy_x, busy_y x, y position of busy window
c     img_width,img_height   width, height of image window
c     men_width, men_height  width, height of menu window  
c     par_width, par_height  width, height of parameter window
c     io_width, io_height    width, height of IO window
c     not_width, not_height  width, height of notice area
c     pmn_width, pmn_height  width, height of Abort menu
c     ppb_width              width of progress bar
c     busy_width, busy_height  width, height of busy box
c     sub_process    .true. if program is sub-process of eg clips, .false.
c                    if standalone
c
C     IFTYPE      1 = unsigned byte data
C                 2 = unsigned two-byte data (i2)
C                 3 = signed integer data
C                 4 = 'squashed i2' data (if Intensity>32767 
C                    store as 65536-Intensity/8)
      INTEGER BASE_WIDTH,BASE_HEIGHT,
     $     IMAGE_ORDER,JIMAGE_ORDER,NXP_CMP,NYP_CMP,
     $     NXDPX,NYDPX,
     $     IMG_HORIZ,IMG_VERT,IMG_MAP(4),
     $     IMG_X,IMG_Y,MEN_X,MEN_Y,PAR_X,PAR_Y,IO_X,IO_Y,
     $     NOT_X,NOT_Y,PMN_X, PMN_Y,PPB_X, PPB_Y, BUSY_X, BUSY_Y,
     $     BUSY_Y2,IFTYPE
      INTEGER IMG_WIDTH, IMG_HEIGHT, MEN_WIDTH, MEN_HEIGHT,
     $     PAR_WIDTH, PAR_HEIGHT, IO_WIDTH, IO_HEIGHT,
     $     NOT_WIDTH, NOT_HEIGHT, PMN_WIDTH, PMN_HEIGHT,
     $     PPB_WIDTH, BUSY_WIDTH, BUSY_HEIGHT
      LOGICAL SUB_PROCESS, BLANK
C

c  DON'T forget to declare function xdlstr as integer
      INTEGER XDLSTR
      EXTERNAL XDLSTR

c
c     ivhbas        view object handle for base frame
c     ivhimg        view object handle for image
c     ivhpar        view object handle for parameter table
c     ivhmen        view object handle for menu
c     ivhio         view object handle for io area
c     ivhio2        view object handle for second io area (pick option)
c     ivhio3        view object handle for third io area (output)
c     ivhnot        view object handle for notice area (parameter table)
c     ivhpmn        view object handle for popup wait menu
c     ivhppb        view object handle for popup progress bar
c     ivhblank      view object handle for blank object
c     ivhbusy       view object handle for busy object
c     ivhbusy2      view object handle for second (2 line) busy object
      INTEGER IVHBAS, IVHIMG, IVHPAR, IVHMEN, IVHIO, IVHNOT, IVHPMN,
     $     IVHPPB, IVHBLANK, IVHBUSY, IVHIO2, IVHBUSY2, IVHIO3
      PARAMETER (IVHBAS=1, IVHIMG=2, IVHPAR=3, IVHMEN=4, IVHIO=5,
     $     IVHNOT=6, IVHPMN=7, IVHPPB=8, IVHBLANK=9, IVHBUSY=10,
     $     IVHIO2=11, IVHBUSY2=12, IVHIO3=13)

c  Border  boundary between windows
      INTEGER BORDER
      PARAMETER (BORDER = 4)

c     icset   colour set number
      INTEGER ICSET
      PARAMETER (ICSET=1)
c
c     iord    axis order of image
c             = 1  xf slow, yf fast
      INTEGER IORD
      PARAMETER (IORD = 1)
c
c     minw,minh   minimum width and height for image display object
      INTEGER MINW,MINH
      PARAMETER (MINW=0, MINH=0)
c
c     ixopix, iyopix      pixel origin of displayed part
      INTEGER IXOPIX,IYOPIX
      PARAMETER (IXOPIX=1, IYOPIX=1)
c
c     ibg     background menu, =0 for none
      INTEGER IBG
      PARAMETER (IBG=0)
c
c     iovly   overlay option
      INTEGER IOVLY
      PARAMETER (IOVLY=1)
c
c     max_pixel   maximum number of pixels to use for image display:
c                 this will be used to control possible compression
c                 of the image
      INTEGER MAX_PIXEL
      PARAMETER (MAX_PIXEL = 800)
c  Image stuff
c     disp_img        true if image has been displayed
c
      LOGICAL DISP_IMG
c
c     ifd     file descriptor, = -1 for no file
      INTEGER IFD
      PARAMETER (IFD=0)
c

c  Parameter table stuff
c     max_par_col     maximum number of parameter columns
c     max_par_rows    maximum number of parameter rows
c     max_par_name    maximum length of name
c     max_par_str     maximum value length
c     par_title       = -1 no title +1 title present
c     par_font        font (3=medium)
c     par_menu        = 1 if popup menus allowed
c     disp_par        true if parameter table hsa been displayed
c
      LOGICAL DISP_PAR
      INTEGER MAX_PAR_COL, MAX_PAR_ROWS, MAX_PAR_NAME, MAX_PAR_STR,
     $        PAR_TITLE, PAR_FONT, PAR_MENU
      PARAMETER (MAX_PAR_COL = 1, MAX_PAR_ROWS = 50,
     $           MAX_PAR_NAME = 17, MAX_PAR_STR = 7,
     $           PAR_TITLE = 1, PAR_FONT = 2, PAR_MENU = 0)
C
c  Menu stuff
c     max_men_itms    maximum number of menu items
c     max_men_name    maximum number of characters in menu item
c     men_font        menu font number (3=medium)
c     men_quit_flag   = 1 to allow for quit box
c     max_men_title   maximum length of menu title
c     disp_menu       true if menu displayed
      LOGICAL DISP_MENU
      INTEGER MAX_MEN_ITMS, MAX_MEN_NAME, MEN_FONT, MEN_QUIT_FLAG,
     $        MAX_MEN_TITLE 
      PARAMETER (MAX_MEN_ITMS = 19, MAX_MEN_NAME = 19, MEN_FONT = 4,
     $     MEN_QUIT_FLAG = 1, MAX_MEN_TITLE = 9)

c
c IO area stuff
c     io_font         font number
c     disp_io         true if io area displayed
c     disp_io2        true if SECOND io area displayed
c     disp_io3        true if THIRD io area displayed
C     NSCROLL         number of pages to hold for scrolling
      LOGICAL DISP_IO,DISP_IO2,DISP_IO3
      INTEGER IO_FONT,NSCROLL
      PARAMETER (IO_FONT=2)
      PARAMETER (NSCROLL=0)


c  Notice area stuff (notice area is a parameter table with no values)
c     max_not_col     maximum number of parameter columns
c     max_not_rows    maximum number of parameter rows
c     max_not_name    maximum length of name
c     max_not_str     maximum value length
c     not_title       = -1 no title +1 title present
c     not_font        font (3=medium)
c     not_menu        = 1 if popup menus allowed
c     disp_not       true if menu displayed
      LOGICAL DISP_NOT
      INTEGER MAX_NOT_COL, MAX_NOT_ROWS, MAX_NOT_NAME, MAX_NOT_STR,
     $        NOT_TITLE, NOT_FONT, NOT_MENU
      PARAMETER (MAX_NOT_COL = 1, MAX_NOT_ROWS = 18,
     $           MAX_NOT_NAME = 23,
     $           MAX_NOT_STR = 0,
     $           NOT_TITLE = +1, NOT_FONT = 2, NOT_MENU = 0)
c
c  Menu stuff for Abort menu
c     max_pmn_itms    maximum number of menu items
c     max_pmn_name    maximum number of characters in menu item
c     pmn_font        menu font number (3=medium)
c     pmn_quit_flag   = 1 to allow for quit box
c     max_pmn_title   maximum length of menu title
      INTEGER MAX_PMN_ITMS, MAX_PMN_NAME, PMN_FONT, PMN_QUIT_FLAG,
     $        MAX_PMN_TITLE 
      PARAMETER (MAX_PMN_ITMS = 0, MAX_PMN_NAME = 12, PMN_FONT = 4,
     $     PMN_QUIT_FLAG = 1, MAX_PMN_TITLE = 15)

c   Stuff for progress bar
c     ppb_font   font number (3=medium)
c     ppb_colr   colour
      INTEGER PPB_FONT, PPB_COLR
      PARAMETER (PPB_FONT = 4, PPB_COLR=3)

c   Busy box
      INTEGER BUSY_FONT
      PARAMETER (BUSY_FONT = 3)
c
c len_dialog    length of dialog box
      INTEGER LEN_DIALOG
      PARAMETER (LEN_DIALOG = 60)

c
c  Overlays
c
c  Circles
c     cir_ivec    vector set number for circles
c     cir_colr    colour for circles
c     cir_iovl    overlay number for circles
c     cir_symb   symbol for centre
      INTEGER CIR_IVEC, CIR_COLR, CIR_IOVL, CIR_SYMB
      PARAMETER (CIR_IVEC=3, CIR_COLR=6, CIR_IOVL=2, CIR_SYMB=4)

C
C---- Residual vectors
C     IRV_VEC Vector set number
C     IRV_COL Colour  (1 red, 2 yellow, 3 green)
C     IRV_OVL Overlay number
      INTEGER IRV_VEC, IRV_COL, IRV_OVL
      PARAMETER (IRV_VEC=4, IRV_COL=1, IRV_OVL=2)
c  Measure
c
c     mes_colr    colour for measure crosses
      INTEGER MES_COLR
      PARAMETER (MES_COLR=1)

c   Vertical Crosses
c     cross_iovl    overlay number for circles
c     cross_symb   symbol for centre
      INTEGER CROSS_IOVL, CROSS_SYMB
      PARAMETER (CROSS_IOVL=2, CROSS_SYMB=3)
c  Boxes
c     box_iovl    overlay number for boxes
c     box_symb   symbol for centre
      INTEGER BOX_IOVL, BOX_SYMB, BOX_VEC
      PARAMETER (BOX_VEC=2,BOX_IOVL=2, BOX_SYMB=13)
c  Circles
c     CIRC_iovl    overlay number for circles
c     CIRC_symb   symbol for centre
      INTEGER CIRC_IOVL, CIRC_VEC
      PARAMETER (CIRC_VEC=5,CIRC_IOVL=2)
c   Crosses
c     cross_iovl    overlay number for circles
c     cross_symb   symbol for centre
      INTEGER XCROSS_IOVL, XCROSS_SYMB
      PARAMETER (XCROSS_IOVL=2, XCROSS_SYMB=13)
C
C---- common block last of all
C
      COMMON /MXDCMI/ BASE_WIDTH,BASE_HEIGHT,
     $     IMAGE_ORDER,JIMAGE_ORDER,NXP_CMP,NYP_CMP,
     $     NXDPX,NYDPX,
     $     IMG_HORIZ,IMG_VERT,IMG_MAP,
     $     IMG_X,IMG_Y,MEN_X,MEN_Y,PAR_X,PAR_Y,IO_X,IO_Y,
     $     NOT_X,NOT_Y,PMN_X, PMN_Y,PPB_X, PPB_Y, BUSY_X, BUSY_Y,
     $     BUSY_Y2, IFTYPE,
     $     IMG_WIDTH, IMG_HEIGHT, MEN_WIDTH, MEN_HEIGHT,
     $     PAR_WIDTH, PAR_HEIGHT, IO_WIDTH, IO_HEIGHT,
     $     NOT_WIDTH, NOT_HEIGHT, PMN_WIDTH, PMN_HEIGHT,
     $     PPB_WIDTH, BUSY_WIDTH, BUSY_HEIGHT,
     $     SUB_PROCESS,BLANK,DISP_IMG,DISP_PAR,DISP_MENU,
     $     DISP_IO,DISP_IO2,DISP_IO3,DISP_NOT
  
C&&*&& end_include  ../inc/mxdinc.f
C&&*&& include  ../inc/myprof.f
C
C $Id: myprof.f,v 1.3 2011/07/19 11:18:20 andrew Exp $
C
C--- awk generated include file  myprof.h
C---- START of include file myprof.h
C
C
C     .. Arrays in Common Block /MYPROF/ ..
      REAL PRSCALE,XPMAX,YPMAX,XLINE,YLINE,WPROFL,WPRSUMS,PRCUT
      INTEGER ISIZE,NRFBOX,IPROFL,IBOX,IXBOX,IYBOX,NPFIRST
      LOGICAL BOX
C     ..
C     .. Scalars in Common Block /MYPROF/
      REAL TOL,BGPKRAT,FRACREJ,BADTOL,TOLMIN,RECLEVEL
      INTEGER NXLINE,NYLINE,NUMBOX,IBOUND,ITRIM,NOVERLAP,MINPKSIZ
      LOGICAL HIGHRES,LOWRES,PRSET,OFFDET,PROPT,FIXBOX,LINESET,
     +        PROPTCEN,RECOVER,NOFIXBOX,PREXICE
C     ..
C     .. Common Block /MYPROF/ ..
      COMMON /MYPROF/PRSCALE(NMASKS,2),XPMAX(NMASKS),YPMAX(NMASKS),
     $       XLINE(NNLINE),YLINE(NNLINE),
     +       WPROFL(MAXBOX,NMASKS),WPRSUMS(MAXBOX,NMASKS),TOL,
     +       BGPKRAT,FRACREJ,BADTOL,TOLMIN,RECLEVEL,
     $       ISIZE(NMASKS,2),NRFBOX(NMASKS),
     +       IPROFL(MAXBOX,NMASKS+1),IBOX(0:NNLINE,0:NNLINE),
     +       IXBOX(NMASKS,2),IYBOX(NMASKS,2),
     +       NPFIRST(NNLINE-1),NXLINE,NYLINE,NUMBOX,IBOUND,ITRIM,
     $       NOVERLAP,MINPKSIZ,BOX(NMASKS),HIGHRES,LOWRES,PRSET,
     +       OFFDET,PROPT,FIXBOX,LINESET,PROPTCEN,RECOVER,NOFIXBOX,
     +       PREXICE,PRCUT
C     ..
C
C
C&&*&& end_include  ../inc/myprof.f
C&&*&& include  ../inc/ori.f
C
C $Id: ori.f,v 1.3 2018/04/06 12:43:44 andrew Exp $
C
C--- awk generated include file  ori.h
C---- START of include file ori.h
C
C     XCEN,YCEN    Coordinates (in 10 micron units) of the direct beam
C                  position relative to an origin at the position of the
C                  first pixel in the digitised image.(The SCANNER
C                  coordinate frame). These parameters are refined for
C                  each image. 
C
C     XCEN0,YCEN0  Coordinates of direct beam position at zero swing angle.
C                  (Needed for pxtomm conversion for swung detectors)
C                  These values are assigned on the basis of input direct
C                  beam coordinates, corrected for swing angle if necessary.
C                  They are not (currently) updated during refinement.
C
C     XOFF,YOFF    Distance between centre of detector and direct beam.
C
C     ..
C     .. Arrays in common /ORI/ ..
      LOGICAL FIXPAR
C
C     .. Scalars in common block /ORI/ ..
      REAL COSOM0,SINOM0,XTOFRA,YSCAL,TILT,TWIST,BULGE,
     +     VBNEG,VBPOS,VTILT,VTWIST,VVERT,TOFF,BULGE2,PBULGE,ROFF,
     +     RDTOFF,RDROFF,ROFFPHI,TOFFPHI,JUNK3,XCEN,YCEN,XCEN0,YCEN0,
     +     XOFF,YOFF,TILTMAT,TWISTMAT,DETNOR,DPSIX,SMTWIST,SMTILT,
     +     SMOMEGA0,SMYSCAL,SMTOFF,SMROFF,SMXTOFRA
      INTEGER CBAR,NODES,NPHI,NSM1,NSM2
      LOGICAL RESETCCOM,NOSMOOTH,DAMPREF
C     ..
C     .. Common Block /ORI/ ..
      COMMON /ORI/COSOM0,SINOM0,XTOFRA,YSCAL,TILT,TWIST,BULGE,
     $       VBNEG,VBPOS,VTILT,VTWIST,VVERT,TOFF,BULGE2,PBULGE,ROFF,
     $       RDTOFF,RDROFF,ROFFPHI,TOFFPHI,JUNK3,XCEN,YCEN,XCEN0,
     +       YCEN0,XOFF,YOFF,TILTMAT,TWISTMAT,DPSIX,DETNOR(3),CBAR,
     +       NODES,NPHI,FIXPAR(NRPAR),
     $       RESETCCOM,SMTWIST,SMTILT,SMOMEGA0,SMYSCAL,SMTOFF,SMROFF,
     +       SMXTOFRA,NSM1,NSM2,NOSMOOTH,DAMPREF
C     ..
C
C
C&&*&& end_include  ../inc/ori.f
C&&*&& include  ../inc/over.f
C
C $Id: over.f,v 1.1 2002/05/02 10:47:03 harry Exp $
C
C--- awk generated include file  over.h
C---- START of include file over.h
C
C
C     .. Scalars in Common Block /OVER/ ..
C
C     ..
C     .. Scalars in common /OVER/ ..
      INTEGER MINDTX,MINDTY,NPOVL,IXSEP,IYSEP
C     ..
C     .. ARRAYS in common /OVER/ ..
      INTEGER HKLPOVL
C     ..
C     .. Common Block /OVER/ ..
      COMMON /OVER/MINDTX,MINDTY,NPOVL,IXSEP,IYSEP,HKLPOVL(NREFLS/2)
C     ..
C
C
C&&*&& end_include  ../inc/over.f
C&&*&& include  ../inc/params.f
C
C $Id: params.f,v 1.1 2002/05/02 10:47:04 harry Exp $
C
C--- awk generated include file  params.h
C---- START of include file params.h
C
C
C     .. Scalars in common block /PARAMS/ ..
      INTEGER NSDR
C     ..
C     .. Common Block /PARAMS/ ..
      COMMON /PARAMS/NSDR
C     ..
C
C
C
C&&*&& end_include  ../inc/params.f
C&&*&& include  ../inc/parm1.f
C
C $Id: parm1.f,v 1.4 2021/08/19 10:53:54 andrew Exp $
C
C--- awk generated include file  parm1.h
C---- START of include file parm1.h
C
C     SCALEINT  Can be set by keyword SCALE, used to scale up
C               intensities (and sd's) if they are numerically small
C               eg as in MicroED. Default 1.0, but 10.0 for MicroED
C
C     SCALEUP   Can be set by keyword SCALEUP and is used internally
C               to scale up intensities and sigmas to avoid truncation
C               or rounding errors when dealing with very small intensities
C               (less than 10), intorduced to help avoid BIAS in intensities
C               that should be close to zero. Default 100.0
C
C     .. Scalars in common block /PARM1/ ..
      REAL BGSIG,EFAC,EFACSQ,BGFREJ,GRADMAX,GRADMAXR,BGRAT,PKRAT,
     +     RESD,RESDLOW,SCALEINT,SCALEUP
      INTEGER CUTOFF,NOVPIX,NDSTART,NDTOT,IXDMIN,IXDMAX,IYDMIN,IYDMAX,
     +     IDMIN,IDMAX,IGRADMAX,IGRADMAXR
      LOGICAL DUMPSPOT,BADPLOT,DUMPALL,BADPLOT2
C     ..
C     .. Common Block /PARM1/ ..
      COMMON /PARM1/BGSIG,EFAC,EFACSQ,BGFREJ,GRADMAX,GRADMAXR,
     +     BGRAT,PKRAT,RESD,RESDLOW,CUTOFF,NOVPIX,
     +     NDSTART,NDTOT,IXDMIN,IXDMAX,IYDMIN,IYDMAX,
     +     IDMIN,IDMAX,IGRADMAX,IGRADMAXR,DUMPSPOT,
     $     BADPLOT,DUMPALL,BADPLOT2,SCALEINT,SCALEUP
C     ..
C
C
C&&*&& end_include  ../inc/parm1.f
C&&*&& include  ../inc/parm2.f
C
C $Id: parm2.f,v 1.3 2021/03/15 11:24:31 andrew Exp $
C
C--- awk generated include file  parm2.h
C---- START of include file parm2.h
C
C
C     .. Scalars in common block /PARM2/ ..
      REAL PRBGSIG,RMSBGPR,DISCRIM,PKWDLIM1,PKWDLIM2,PKWDLIM3,PKWDOUTL,
     +     FBGNULL,FPKNULL,FBGEDGE,FPKEDGE
      INTEGER ISDRATIO,PRCUTOFF,IPLOT,NRFMIN,NHKLD,IOUTL1,IOUTL2
      LOGICAL PROFILE,PROCES,SAVEFILE,CHANGEMASK,PRBFILM,PRCFILM,
     +        PRREAD,PRSAVE,WEIGHT,PRPART,USEOVRLD,USEDGE,VARPRO,
     +        WTPROFILE,DISCRIMINATE,PUPDATE,PKONLY,DENSE,PRFULLS,
     +        DECONV,CLEAN
C     ..
C     .. Arrays in common blocks /PEL/ and /PELC/ ..
      INTEGER IHD
C     ..
C     .. Common Block /PARM2/ ..
      COMMON /PARM2/PRBGSIG,RMSBGPR,DISCRIM,PKWDLIM1,PKWDLIM2,
     $       PKWDLIM3,PKWDOUTL,FBGNULL,FPKNULL,FBGEDGE,FPKEDGE,
     $       ISDRATIO,PRCUTOFF,IPLOT,NRFMIN,NHKLD,IOUTL1,IOUTL2,
     $       IHD(3,50),
     $       PROFILE,PROCES,SAVEFILE,CHANGEMASK,PRBFILM,PRCFILM,
     +       PRREAD,PRSAVE,WEIGHT,
     +       PRPART,USEOVRLD,USEDGE,VARPRO,
     +       WTPROFILE,DISCRIMINATE,
     +       PUPDATE,PKONLY,DENSE,
     +       PRFULLS,DECONV,CLEAN
C     ..
C
C
C&&*&& end_include  ../inc/parm2.f
C&&*&& include  ../inc/pel.f
C     **** This is the IMAGE PLATE version of this common block ****
C
C
C---- START of include file pel.h
C
C
C     .. Scalars in common blocks /PEL/ and /PELC/ ..
      INTEGER IBA,IPOINT,ISTART
      LOGICAL INCORE
C     ..
C     .. Arrays in common blocks /PEL/ and /PELC/ ..
      LOGICAL RDSTRIP
      INTEGER*2 BOXOD,IMAGE,IMAGE2
C     ..
C     .. Common Blocks /PEL/ and /PELC/ ..
      COMMON /PEL/IBA,IPOINT,INCORE,RDSTRIP(IXWDTH),BOXOD(MAXBOX),
     +            ISTART
      COMMON /PELC/ IMAGE(IYLENGTH*IXWDTH),IMAGE2(4400*4400)
C     ..
C
C
C---- END of include file pel.h
C
C&&*&& end_include  ../inc/pel.f
C&&*&& include  ../inc/phiprofile.f
C
C
C     .. Scalars in Common Block /PHIPROFILE/ ..
      INTEGER IHPROF,IKPROF,ILPROF,ISTARTLAST,IENDLAST,NBITPROF,IMGPAD,
     +        NPHIPROFIMG,NPHIPROFMISS
      REAL ETASAVE,FRACREAL,RESSAVE,RESLOWSAVE,DSTMIN2PHIP,DSTMAX2PHIP
      LOGICAL PHIPROF,RESTRESL,REPREDICT,STARTTRUNC,ENDTRUNC
C     ..
C     .. Arrays in common /PHIPROFILE/
      INTEGER INTPROFILE,IDPROFILE,IMGMISSING
      REAL FRACPROFILE
C
C     .. Common Block /PHIPROFILE/ ..
        COMMON /PHIPROFILE/ IHPROF,IKPROF,ILPROF,ISTARTLAST,IENDLAST,
     +              NBITPROF,IMGPAD,NPHIPROFIMG,NPHIPROFMISS,ETASAVE,
     +              FRACREAL,RESSAVE,RESLOWSAVE,DSTMIN2PHIP,DSTMAX2PHIP,
     +              PHIPROF,RESTRESL,REPREDICT,STARTTRUNC,ENDTRUNC,
     +              INTPROFILE(NGRA),IDPROFILE(NGRA),IMGMISSING(NGRA),
     +              FRACPROFILE(NGRA)
C&&*&& end_include  ../inc/phiprofile.f
C&&*&& include  ../inc/postchk.f
C
C $Id: postchk.f,v 1.5 2015/05/11 14:19:27 andrew Exp $
C
C--- awk generated include file  postchk.h
C---- START of include file postchk.h
C
C     DELPHIV   stores the missets for the first NADD images only,
C               so that once NADD images have been processed their
C               refined missets can be saved for writing to the
C               summary file.
C
C     NREFMOS   Minimum number of reflections for mosaicity refinement
C               to be carrid out (default 100)
C
C     .. Arrays in Common /POSTCHK/
      REAL SDBEAM,SDCELL,SDDELPHI,DELPHIV,SHIFT
      LOGICAL FCELL,VALSET
C
C     .. Scalars in Common /POSTCHK/ ..
C     ..
      REAL PRRES1,PRRES2,RESIDMAX,SDFAC,SHIFTMAX,SHIFTFAC,
     +     ANGWIDTH,CELLSHIFT,FRACMIN,FRACMAX,FRCSHIFT
      INTEGER PRNS,NREFPR,NADD,NRPT,IPRINTP,NPRMIN,NREFMOS
      LOGICAL USEBEAM,REFCELL,RESCALE,SMOOTHPHI
C     ..
C     .. Common block /POSTCHK/  ..
      COMMON /POSTCHK/ SDBEAM(6),SDCELL(6),SDDELPHI(2*NIMAX),
     +             DELPHIV(3*NIMAX),SHIFT(3),PRRES1,PRRES2,RESIDMAX,
     +             SDFAC,SHIFTMAX,SHIFTFAC,ANGWIDTH,CELLSHIFT, 
     $             FRACMIN,FRACMAX,FRCSHIFT,
     +             PRNS,NREFPR,NADD,NRPT,IPRINTP,
     +             NPRMIN,FCELL(6),USEBEAM,REFCELL,RESCALE,
     +             VALSET(MAXPAX),NREFMOS,SMOOTHPHI
C     ..
C
C
C&&*&& end_include  ../inc/postchk.f
C&&*&& include  ../inc/praccum.f
C
C $Id: praccum.f,v 1.1 2002/05/02 10:47:08 harry Exp $
C
C--- awk generated include file  praccum.h
C---- START of include file praccum.h
C
C
C     ..
C     .. Arrays in common block /PRACCUM/ ..
      REAL PRDATA
C
C     .. Scalars in common block /PRACCUM/ ..
      LOGICAL ACCUMULATE,FIRSTPASS,SECONDPASS,THIRDPASS,FIRSTFILM,
     +        NOTREAD
C     ..
C     .. Common Block /PRACCUM/ ..
      COMMON /PRACCUM/PRDATA(4,50),ACCUMULATE,FIRSTPASS,SECONDPASS,
     +       THIRDPASS,FIRSTFILM,NOTREAD
C     ..
C
C
C&&*&& end_include  ../inc/praccum.f
C&&*&& include  ../inc/precession.f
C
C $Id: precession.f,v 1.1 2002/05/02 10:47:08 harry Exp $
C
C--- awk generated include file  precession.h
C---- START of include file precession.h
C
C
C     .. Scalars in common block /PRECESSION/ ..
      REAL D1,D2,PHIPREC,PSIPREC,PCTOFD,XLAMBDA
      LOGICAL PRECESS
C     ..
C     .. Common Block /PRECESSION/ ..
      COMMON /PRECESSION/D1,D2,PHIPREC,PSIPREC,PCTOFD,XLAMBDA,PRECESS
C     ..
C
C
C&&*&& end_include  ../inc/precession.f
C&&*&& include  ../inc/ras.f
C
C $Id: ras.f,v 1.1 2002/05/02 10:47:10 harry Exp $
C
C--- awk generated include file  ras.h
C---- START of include file ras.h
C
C
C     .. Scalars in common block /RAS/ ..
      INTEGER NEWRAS,MINT
C     ..
C     .. Arrays in common block /RAS/ ..
      INTEGER IRAS
      REAL VARAS
C     ..
C     .. Common Block /RAS/ ..
      COMMON /RAS/VARAS(5),IRAS(5),NEWRAS,MINT
C     ..
C
C
C&&*&& end_include  ../inc/ras.f
C&&*&& include  ../inc/reeke.f
C
C $Id: reeke.f,v 1.17 2015/05/11 14:19:27 andrew Exp $
C
C--- awk generated include file  reeke.h
C---- START of include file reeke.h
C
C     TOR     For synchrotron sources, degree of polarisation of the beam
C     IMONO   Collimation flag for polarisation.
C             = 0   Pinhole or mirrors
C             = 1   Graphite Monochromator
C             = 2   Synchrotron, use TOR
C     NWMAX   Maximum reflection width in images
C     DSTMAX  dimensionless rlu, = WAVE/RES where RES is the maximum
C             resolution in Angstrom.
C     DSTMIN  dimensionless rlu, = WAVE/DMAX where DMAX is the maximum
C              Bragg spacing
C     S0      radiation source vector (-1,0,0), may be modified by MU values
C     SOURCE  radiation source vector - basis for s0
C     GONIO_angle:  machine angles for generalized goniostats
C     GONIO_angle_SET: fixed offset angles for 3- and 4-circle goniostats
C                      i.e. currently 54.74 deg for Bruker CHI, 50.0deg for
C                      Nonius and Oxford KAPPA 
C     VERTICAL Logical set TRUE if rotation axis is vertical, only
C              affects polarization correction for synchrotron data
C     NSMOOTH  Number of images to use in smoothing refined values of
C              mosaic spread and missetting angles (if being smoothed).
C     ..
C     .. Arrays in common block /REEKE/ ..
      REAL RMC,AMAT,SCAXIS,S0,s0unit,CAMCON,BEAMMU
C     ..
C     .. Scalars in common block /REEKE/ ..
      REAL X1,X2,Y1,Y2,Z1,Z2,XYS,DSTAR2,DSTMAX,DSTPL,DSTPL2,
     +     DIVH,DIVV,DELAMB,ETA,DELCOR,DSTMAXS,TOR,DSTMIN,WMAX,
     +     ETAFRAC,ETASM,
     $     GON_CHI,GON_PHIS,GON_KAPPA,GON_OMEGAS,GON_2THETA,
     $     GON_PHIE,GON_OMEGAE,
     +     GON_CHI_SET,GON_KAPPA_SET,
     +     BLKSIZE
      INTEGER ISYN,IMONO,NWMAX,IPAD,NSMOOTH,IGON_CHI,IGON_PHI,
     $     IGON_KAPPA,IGON_OMEGA,IGON_2THETA,
     +     IGON_CHI_SET,IGON_KAPPA_SET,IHOSCAXIS
C     ..
C     .. logicals in common block /REEKE/ ..
      LOGICAL LOGETA,NUREEK,LCAMCON,GENGON,FOURCIRCLE,THREECIRCLE,
     +        VERTICAL
C     ..
C     .. Common Block /REEKE/ ..
      COMMON /REEKE/RMC(3,3),AMAT(3,3),SCAXIS(3),S0(3),s0unit(3),
     $     CAMCON(8),BEAMMU(2),X1,X2,Y1,Y2,
     +     Z1,Z2,XYS,DSTAR2,DSTMAX,DSTPL,DSTPL2,DIVH,DIVV,DELAMB,
     +     ETA,DELCOR,DSTMAXS,TOR,DSTMIN,WMAX,ETAFRAC,ETASM,
     $     GON_CHI,GON_PHIS,GON_KAPPA,GON_OMEGAS,GON_2THETA,
     $     GON_PHIE,GON_OMEGAE,
     +     GON_CHI_SET,GON_KAPPA_SET,ISYN,IMONO,
     $     NWMAX,IPAD,
     $     NSMOOTH,IGON_CHI,IGON_PHI,
     $     IGON_KAPPA,IGON_OMEGA,IGON_2THETA,
     +     IGON_CHI_SET,IGON_KAPPA_SET,IHOSCAXIS,LOGETA,NUREEK,LCAMCON,
     $     GENGON,FOURCIRCLE,THREECIRCLE,
     +     BLKSIZE,VERTICAL
C     ..
C
C
C&&*&& end_include  ../inc/reeke.f
C&&*&& include  ../inc/reprt.f
C--- awk generated include file  reprt.h
C---- START of include file reprt.h
C
C
C
C     .. Scalars in common block /REPRT/ ..
      REAL AVBGRATIO,AVSIG,RFULL,RPART,AVPKRATIO,PKRMAXI,BGRMAXI,
     +     RFACOV,SDRATOV,SDMON,RESCUT,STHCUT
      INTEGER NREF,NOFR,NOLO,MAXBSI,MINBSI,NEDGE,NBOX,NBZERO,NNEG,
     +        NBAD,NBGRJ,NEDGE1,NPARTEND,NSPOVL,NRSYM,NHALF,NSUMPART,
     +        NOVRLAP,MSPOT,NLATTOVRLAP,NWIDE,NOUTSIDE,NBADRES,NINCUSP,
     +        NEDGE2
      LOGICAL PKACCEPT
C     ..
C     .. Arrays in common block /REPRT/ ..
      REAL RATIO,AVSD,AVSDP,AVSIG1,AVSIG2,AVINTI1,AVINTI2,AVPRI1,
     +     AVPRI2,RMSDELI1,RMSDELI2,ABSDELI1,ABSDELI2,AVPRSIG1,
     +     AVPRSIG2,AVDELSIG1,AVDELSIG2,PKRATIO,RSIGVSM,
     +     FIOVSDP,FIOVSDS,PIOVSDP,PIOVSDS,DBIN
      INTEGER IRANGE,IANAL,IANALF,NRFLS1,NRFLS2,MEANDELI1,MEANDELI2,
     +        NBGRHIST,IVSM,NIVSM,NRESPF,NRESPP,NRESSF,NRESSP,
     +        IRESPF,IRESPP,IRESSF,IRESSP,ISDRESPF,ISDRESPP,
     +        ISDRESSF,ISDRESSP
C     ..
C     .. Common Block /REPRT/ ..
      COMMON /REPRT/RATIO(10),AVSD(10),AVSDP(10),AVSIG1(10,2),
     $       AVSIG2(NMASKS,2),AVINTI1(10,2),AVINTI2(NMASKS,2),
     $       AVPRI1(10,2),AVPRI2(NMASKS,2),RMSDELI1(10,2),
     +       RMSDELI2(NMASKS,2),ABSDELI1(10,2),ABSDELI2(NMASKS,2),
     +       AVPRSIG1(10,2),AVPRSIG2(NMASKS,2),AVDELSIG1(10,2),
     +       AVDELSIG2(NMASKS,2),PKRATIO(10),RSIGVSM(13),FIOVSDP(9),
     $       FIOVSDS(9),PIOVSDP(9),PIOVSDS(9),DBIN(9),
     $       AVBGRATIO,AVSIG,RFULL,RPART,AVPKRATIO,PKRMAXI,BGRMAXI,
     $       RFACOV,SDRATOV,SDMON,RESCUT,STHCUT,
     $       IRANGE(9),IANAL(10),IANALF(10),NRFLS1(10,2),
     +       NRFLS2(NMASKS,2),MEANDELI1(10,2),MEANDELI2(NMASKS,2),
     +       NBGRHIST(32),IVSM(13),NIVSM(13),
     +       NRESPF(9),NRESPP(9),NRESSF(9),NRESSP(9),
     +       IRESPF(9),IRESPP(9),IRESSF(9),IRESSP(9),ISDRESPF(9),
     +       ISDRESPP(9),ISDRESSF(9),ISDRESSP(9),NREF,NOFR,NOLO,
     $       MAXBSI,MINBSI,NEDGE,NBOX,NBZERO,NNEG,NBAD,
     +       NBGRJ,NEDGE1,NPARTEND,NSPOVL,
     +       NRSYM,NHALF,NSUMPART,PKACCEPT,NOVRLAP,MSPOT,
     +       NLATTOVRLAP,NWIDE,NOUTSIDE,NBADRES,NINCUSP,NEDGE2
C     ..
C
C
C&&*&& end_include  ../inc/reprt.f
C&&*&& include  ../inc/restart.f
C
C $Id: restart.f,v 1.2 2005/09/20 16:41:57 andrew Exp $
C
C---- START of include file restart.h
C
C---- Saves image ID, IDENT and phi values in case of an abort or 
C     when current run has finished, determines which image will be
C     displayed with GUI at this point.
C
C     .. Scalars in common block /RESTART/ ..
      REAL RESTPHIB, RESTPHIE
      INTEGER RESTID,RNTDIG
      CHARACTER RESTIDENT*40,RTEMPLSTART*100,RTEMPLEND*100
C     ..
C     .. Arrays in common block /RESTART/ ..
C
C     .. Common Block /RESTART/ ..
      COMMON /RESTART/ RESTPHIB,RESTPHIE,RESTID,RNTDIG
C
C
      COMMON /RESTARTC/ RESTIDENT,RTEMPLSTART,RTEMPLEND
C     ..
C     ..
C
C
C&&*&& end_include  ../inc/restart.f
C&&*&& include  ../inc/rfs.f
C
C $Id: rfs.f,v 1.2 2007/07/12 13:25:57 andrew Exp $
C
C--- awk generated include file  rfs.h
C---- START of include file rfs.h
C
C
C
C     .. Scalars in common block /RFS/ ..
      REAL FDIST,ESTART,RMSRES,WRMSRES,WESTART
      INTEGER MAXREF,MAXX,MAXY,NRS,MAXR,NREJS,NCENTREF
C     ..
C     .. Arrays in common block /RFS/ ..
      REAL XRS,YRS,WXRS,WYRS
      INTEGER RRS,IHKLR
C     ..
C     .. Common Block /RFS/ ..
      COMMON /RFS/ XRS(NDETREF),YRS(NDETREF),WXRS(NDETREF),
     +       WYRS(NDETREF),FDIST,ESTART,RMSRES,WRMSRES,
     $       WESTART,RRS(NDETREF),IHKLR(3,NDETREF),MAXREF,MAXX,
     +       MAXY,NRS,MAXR,NREJS,NCENTREF
C     ..
C
C
C&&*&& end_include  ../inc/rfs.f
C&&*&& include  ../inc/savall.f
C
C $Id: savall.f,v 1.5 2018/05/09 14:48:18 andrew Exp $
C
C---- START of include file savall.h
C
C      NSAVIMG      Number of images used in autoindexing
C      ISAVIMG      Array storing image numbers used in last autoindexing
C      NSAVSEG      Number of segments used in last postref run
C      ISFIRST      Array storing image number of the first image
C                   in all segments of last postref run.
C      SAVMATSR     Indicates whether last matrix was determined by
C                   autoindexing (Autoindexing) or by postref (Post refinement)
C      SAVMATNAM    Name of the matrix file
C      SAVENAM      Name of savefile
C      SVSCN        SCANNER keyword
C      SVSITE       SITE keyword
C      RES          High resolution limit given on RESOL keyword
C      RESLOW       Low resolution limit given on RESOL keyword
C      RES2SD       estimate of high resolution limit from spot finding
C      IISCN        Set TRUE if a DETECTOR keyword has been read that specifies
C                   a particular detector type.
C      GENERICDET   Set TRUE if a generic detector specification (ROTATION,
C                   ORIGIN, FAST and RECTANGULAR/CIRCULAR) has been given.
C
C
C     .. Scalars in common block /SAVALL/ ..
      REAL RES,RESLOW
      INTEGER NSAVIMG,NDIR,NSAVSEG
      CHARACTER SAVMATSTR*80,SAVMATNAM*80,SAVENAM*80,SVSCN*80,
     +          SVSITE*80
      LOGICAL IISCN,IISITE,IIWAVE,IIDIV,IIDISP,IINULL,IIRAST,IISEP,
     +        IIOVER,IIPIX,IIBACK,IIRES,GENERICDET
C     ..
C     .. Arrays in common block /SAVALL/ ..
C
      INTEGER ISAVIMG(MAXIMG),ISFIRST(100)
C     .. Common Block /SAVALL/ ..
      COMMON /SAVALL/ RES,RESLOW,
     +       ISAVIMG,ISFIRST,NSAVIMG,NDIR,NSAVSEG,
     +       IISCN,IISITE,IIWAVE,IIDIV,IIDISP,IINULL,IIRAST,IISEP,
     +        IIOVER,IIPIX,IIBACK,IIRES,GENERICDET
C
C
      COMMON /SAVALLC/ SAVMATSTR,SAVMATNAM,SAVENAM,SVSCN,SVSITE
C     ..
C
C
C&&*&& end_include  ../inc/savall.f
C&&*&& include  ../inc/scn.f
C
C $Id: scn.f,v 1.5 2014/11/18 14:05:40 andrew Exp $
C
C--- awk generated include file  scn.h
C---- START of include file scn.h
C
C     SCNSZ   is pixel size (in microns) divided by 25
C     RAST    pixel size in slow direction in mm
C     FACT    multiplying by FACT converts from 10 micron units 
C             (the standard unit internal to the program) into pixels
C     IYLEN   The number of pixels in the Y (fastest changing) direction
C             in the digitised image.
C     NREC    The number of pixels in the X (slow) direction in the 
C             digitised image.
C     NWORD   The number of I*2 words in the Y direction.
C     NBYTE   The number of bytes in Y direction = NWORD/2
C     NHBYTE  Number of bytes in header
C
C     ICURR   When several images are stored in a single file, ICURR is the
C             pointer to the first record of the current image in the direct
C             access file (only implemented for film data)
C     NEXTRA  The number of additional (unused) bytes padding the end of
C             each record in image file
C     BYTSWAP True if the "endedness" (Big-endian/Little-endian) of the 
C             machine that MOSFLM is running on is different to that of the
C             machine on which the image was written. This is determined by
C             looking at the value of NXPIX in the header record of the
C             image file (subroutine GETHDR)
C
C     LACTIVE True if using Rigaku active mask, False if ignoring it
C     REALLYFORCEREAD True if we want to display an image where we have
C                     changed the activemask
C
C     ROWREADT Time in seconds that the detector takes to read a row of data,
C              required for rolling shutter mode readout of Tietz detector
C              used for microED.
C     ROTSPEED Spindle rotation rate (degrees per second),required for rolling
C              shutter mode readout of Tietz detector used for microED. 
C
C
C     .. Scalars in common block /SCN/ ..
      REAL FACT,SCNSZ,RAST,RASTY,ROWREADT,ROTSPEED
      INTEGER NBYTE,NREC,NWORD,IYLEN,ICURR,NEXTRA,NHBYTE
      LOGICAL BYTSWAP,LACTIVE,REALLYFORCEREAD
C     ..
C     .. Common Block /SCN/ ..
      COMMON /SCN/ FACT,SCNSZ,RAST,RASTY,ROWREADT,ROTSPEED,
     +     NBYTE,NREC,NWORD,IYLEN,
     +     ICURR,NEXTRA,NHBYTE,BYTSWAP,LACTIVE,
     $     REALLYFORCEREAD
C&&*&& end_include  ../inc/scn.f
C&&*&& include  ../inc/spots.f
C
C $Id: spots.f,v 1.4 2010/05/14 15:42:51 harry Exp $
C
C--- awk generated include file  spots.h
C---- START of include file spots.h
C
C     XSPT,YSPT,ISPT store spot coordinates (in mm) and intensity
C
C     These arrays hold the complete list of spots found on ALL images.
C     Pointers to the first and last spot in each image are kept in arrays
C
C     XSPT....  X coordinates (in mm) of spots in total spot list for
C               all images.
C
C     YSPT....  Y coordinates (in mm) of spots in total spot list for
C               all images.
C
C     ISPT....  Intensity of spots in total spot list for
C               all images.
C
C     ISDSPT..  Standard deviations of intensities
C
C     INDX...   Stores a pointer into arrays XSPT,YSPT,ISPT for the I'th
C               displayed spot (only spots above current threshold are
C               displayed).
C
C     ISTIMG(IMAG)...  Pointer to first spot for image IMAG in arrays
C                      XSPT,YSPT,ISPT
C     IENDIMG(IMAG)...  Pointer to last spot for image IMAG in arrays
C                      XSPT,YSPT,ISPT
C
C     IXSPT,IYSPT...   Spot coordinates in display pixels for current image
C     SELECT....       True for images selected for autoindexing
C     SPOTFND....      True for images on which spots have been found
C     ITHRESH          I/sig(I) threshold for spots to be used in autoindexing
C     MAXNREF          Maximum number of spots to be used per image in 
C                      autoindexing.
C
C     .. Scalars in Common Block /SPOTS/ ..
        REAL XF,YF,XFID1,YFID1,XFID2,YFID2,XFID3,YFID3,PHIIMG,
     +       LPIXEL,LYSCALE,LOMEGA,LXCEN,LYCEN
        INTEGER NSPT,ISPOT,ITHRESH,MAXTHRESH,MAXNREF,IRSCAL,LIXLEN,
     +       LIYLEN,INVFLAG,LISWUNG
        CHARACTER*400 SPTNAM,NSPTNAM

C     ..
C     .. Arrays in common /SPOTS/ ..
        REAL XSPT(NSPOTS),YSPT(NSPOTS)
        INTEGER ISPT(NSPOTS),IXSPT(NSPOTS),IYSPT(NSPOTS),
     +          ISTIMG(MAXIMG),IENDIMG(MAXIMG),ISDSPT(NSPOTS)
        INTEGER*2 INDX(NSPOTS)
        LOGICAL SELECT(MAXIMG),SPOTFND(MAXIMG)
C     ..
C     .. Common Block /SPOTS/ ..
        COMMON /SPOTS/ XF,YF,XFID1,YFID1,XFID2,YFID2,XFID3,YFID3,
     +       XSPT,YSPT,ISPT,IXSPT,IYSPT,INDX,ISPOT,NSPT,
     +       ITHRESH,MAXTHRESH,MAXNREF,PHIIMG,IRSCAL,LPIXEL,
     +       LYSCALE,LOMEGA,LIXLEN,LIYLEN,LXCEN,LYCEN,
     +       INVFLAG,LISWUNG,ISTIMG,IENDIMG,SELECT,ISDSPT,
     $       SPOTFND
C
        COMMON /CSPOTS/ SPTNAM,NSPTNAM

C&&*&& end_include  ../inc/spots.f
C&&*&& include  ../inc/spots2.f
C
C $Id: spots2.f,v 1.22 2020/04/29 10:10:15 andrew Exp $
C
C--- awk generated include file  spots.h
C---- START of include file spots.h
C
C---- These stores variables associated with spot-finding routines.
C     Note that the final coordinates are kept in XSPT,YSPT in /SPOTS/
C
C     XSPOT...  X coordinates (in mm) of all spots found on current
C               image, before rejections based on spot size
C     YSPOT...  Y coordinates (in mm) of all spots found on current
C               image, before rejections based on spot size
C     INSPOT... Intensities of all spots found on current
C               image, before rejections based on spot size
C     ISIGSP... Standard deviations of intensities
C     NOIMG     The pack or image number.
C     PHI..     Phi value of the mid-point of the image
C     NIMAG     Total number of images read in (for autoindexing)
C     NEWSPT    If TRUE, new format spot coordinate file, otherwise old
C               format.
C     MEDWXSPOT Median spot size in pixels in X determined in PICKSPOTS
C     MEDWYSPOT Median spot size in pixels in Y determined in PICKSPOTS
C     IRSTRT    Number of the image on which to restart after an abort
C
C     SPXMIN    Minimum X coord (relative to direct beam position) for
C               spots
C     SPYMIN    Minimum Y coord (relative to direct beam position) for
C               spots
C     NSPTD     Number of displayed spots
C     THRESH    Threshold for pixel to be considered part of a spot,
C               is (background + THRESH*sigma)
C     NSPK      Number of maxima in a spot
C     MAXPKS    Maximum number of local maxima in a spot
C     RMSFAC    For multiple peaks treated as a single spot, the rms
C               variation over the pixels in the spot must be more than
C               RMSFAC times the variation expected from counting statistics
C               if all pixels had the background value.
C     MAXYEXT   Maximum extent of a "spot" in Y direction.
C     ANISOMAX  Maximum allowed ratio of spot size in X and Y
C
C     The following are for local background determination
C
C     LOCBG     Logical, TRUE if using local background (default)
C     NXBOX     Size of box for local background determination in X
C     NYBOX     Size of box for local background determination in Y
C     BGLOCAL   Array dimension (MAXBGBOX,MAXBGBOX,4) square holding 
C               local background values.
C---- BGLOCAL(I,J,1) ...   Mean background for box I,J
C---- BGLOCAL(I,J,2) ...   Poissonian sigma (uses GAIN) for box I,J
C---- BGLOCAL(I,J,3) ...   Rms variation in counts for box I,J
C---- BGLOCAL(I,J,4) ...   Average +ve difference between mean background
C                          of box I,J and its neighbours, expressed as
C                          a multiple of the rms variation of box I,J
C                          Used to identify boxes in shadow.
C     IWTBG     Array dimension (MAXBGBOX,MAXBGBOX) holding weight for
C               each local background area.
C     EXCLDICE  True if spots in resolution shells of ice rings are to be 
C               excluded
C     AUTORESLN Set resolution limit to 4A if spots between 4-6 are much
C               weaker than those between 6A and infinity.
C     AUTORES1  The lower resolution bin boundary (default 6A) for AUTORESLN
C     AUTORES2  The higher resolution bin boundary (default 4A) for AUTORESLN
C     AUTOFAC1  The intensity factor (default 0.2) for AUTORESLN
C     AUTOFAC2  The reflection number factor (default 0.25) for AUTORESLN
C     AUTORING  Automatically exclude radial bins that have too many spots
C               eg due to ice rings.
C     RINGFACT  Rings with more than RINGFACT*(average number of spots) are
C               rejected as potential ice/powder rings. Default 4.0.
C     NUMICE    Number of spots rejected as potentailly being ice spots
C     NTOOFEWPIX Number of spots rejected as having too few pixels.
C     NIMGSUMSPOT Number of images to sum when searching for spots for indexing.
C     NRESOLU   Number of spots rejected by automatic resolution reduction to 4.5
C
C
C     .. Scalars in Common Block /SPOTS2/ and /SPOTS2C/ ..
      REAL CUTPIXMAX,CUTPIXMIN,CUTWXMAX,CUTWXMIN,CUTWYMAX,CUTWYMIN,
     +     XSPLIT,YSPLIT,YOFFSET,XOFFSET,THRESH,RMINSP,RMAXSP,
     +     RMINSRCH,RMAXSRCH,SCALSRCH,SPXMIN,SPYMIN,THRESHMIN,
     +     THRESHMAX,RMSFAC,ANISOMAX,AUTORES1,AUTORES2,AUTOFAC1,
     +     AUTOFAC2,RINGFACT,RES2SD
      INTEGER NBINR,NBINT,IBINLIM,NPIXMIN,NPIXMAX,IIMAG,NIMAG,
     +     MEDWXSPOT,MEDWYSPOT,NSEARCH,ITHSET,ISAFE,IRSTRT,
     +     NSPTD,MAXPKS,MAXYEXT,NXBOX,NYBOX,IPIXMIN,IRMSFAC,
     +     NUMICE,NTOOFEWPIX,NIMGSUMSPOT,NRESOLU
      LOGICAL RADX,RADY,NEWSPT,DOFIND,LOCBG,EXCLDICE,
     +     AUTORESLN,AUTORING,AUTOTHRESH,FOUND(20)
C     ..
C     .. Arrays in common /SPOTS2/  and /SPOTS2C/ ....
      INTEGER NPIX(MXSPOT),IWXSPOT(MXSPOT),IWYSPOT(MXSPOT),
     $     IORDER(MXSPOT),INSPOT(MXSPOT),NOIMG(MAXIMG),ISDSPOT(MXSPOT),
     $     NSPK(MXSPOT), IWTBG(MAXBGBOX,MAXBGBOX), AVBG(IXWDTH)
      REAL BGOD(IXWDTH),PHI(MAXIMG),PHISAVE(MAXIMG),
     $     XSPOT(MXSPOT),YSPOT(MXSPOT),RESOLU(MXSPOT),
     $     RSPOT(MXSPOT),PHISTIM(MAXIMG),PHIRNGA(MAXIMG),
     +     BGODA2(21,21,3)
C     ..
C     .. Common Block /SPOTS2/  and /SPOTS2C/ ....
      COMMON /SPOTS2/
C     .. real scalars
     $     CUTPIXMAX,CUTPIXMIN,CUTWXMAX,CUTWXMIN,CUTWYMAX,CUTWYMIN,
     +     XSPLIT,YSPLIT,YOFFSET,XOFFSET,THRESH,RMINSP,RMAXSP,
     +     RMINSRCH,RMAXSRCH,SCALSRCH,SPXMIN,SPYMIN,THRESHMIN,
     +     THRESHMAX,RMSFAC,ANISOMAX,AUTORES1,AUTORES2,AUTOFAC1,
     +     AUTOFAC2,RINGFACT,RES2SD,
C     .. integer scalars
     $     NBINR,NBINT,IBINLIM,NPIXMIN,NPIXMAX,IIMAG,NIMAG,
     +     MEDWXSPOT,MEDWYSPOT,NSEARCH,ITHSET,ISAFE,IRSTRT,
     +     NSPTD,MAXPKS,MAXYEXT,
     +     NXBOX,NYBOX,
     +     IPIXMIN,IRMSFAC,NUMICE,NTOOFEWPIX,NIMGSUMSPOT,NRESOLU,
C     -- logicals 
     $     RADX,RADY,NEWSPT,DOFIND,LOCBG,EXCLDICE,AUTORESLN,
     $     AUTORING,AUTOTHRESH,FOUND,
C     .. integer arrays
     $     NPIX,IWXSPOT,IWYSPOT,IORDER,INSPOT,NOIMG,ISDSPOT,NSPK,
     +        IWTBG,AVBG,
C     .. real arrays
     $     BGOD,PHI,XSPOT,YSPOT,RESOLU,RSPOT,PHISTIM,PHIRNGA,
     $     BGODA2,PHISAVE

      REAL BGLOCAL(MAXBGBOX,MAXBGBOX,4)
       
      COMMON /SPOTS3/BGLOCAL

      CHARACTER IMGTEMPL*100,IMGFN*100,IMGNUM*9
      COMMON /SPOTS2C/ IMGTEMPL,IMGFN,IMGNUM
C&&*&& end_include  ../inc/spots2.f
C&&*&& include  ../inc/strat.f
C--- awk generated include file  strat.h
C---- START of include file strat.h
C
C---- Stores variables for use in STRATEGY option
C     IROTAX is the axis closest to the rotation axis
C     PHIROTAX is the angle this axis makes with the rotation axis
C
C     PHIZONE is the phi value at which axis "IZONEAX" is along the
C     X-ray beam. (IZONEAX = 1 is a, =2 is b, =3 is c)
C     PHIPAD  is the rotation to be added to PHILAUE to ensure generation
C             of all unique data, and will depend on angle between the
C             unique axis and the rotation axis.
C     AUTANOM if true tries to maximise number of anomalous pairs
C     CELLSCAL is the scale factor applied to the cell edges to
C              speed up the calculation.
C     SHRUNK   is TRUE if the cell has been scaled by CELLSCAL
C     ISTRUN   is used to increment input phi angles by multiples of
C              360 degrees, so that phi values in different parts can
C              be recognised. It starts at zero.
C     FIRSTRAT Starts as TRUE (set in MOSDATA) and is set to FALSE once
C              the MTZ file has been opened for a STRATEGY run.
C              Reset TRUE after EXITing from STRATGEY prompt.
C
C     .. Scalars in common /STRAT/ ..
      INTEGER NSEGM,IPCKCUR,NSTRAT,NUNIQ,NSTRUN,ISTRUN,NSEGAUTO,
     +        IZONEAX,IROTAX,NLAST,NNPACKS,NLASTPACK
      LOGICAL STRATEGY,AUTO,SIZESET,AUTANOM,SHRUNK,FIRSTRAT,
     +         NEWSTRAT,WAITINP,OFFPHI
      REAL ROTAUTO,PHILAUE,PHIPAD,PHIZONE,PHIROTAX,CELLSCAL
C     ..
C     .. Arrays in common /STRAT/ ..
      REAL PHIST,PHIFIN,PHIINC,PHIADD,PHISEGA
      INTEGER IFIRSTONE
C     ..
C     .. Common Block /STRAT/ ..
      COMMON /STRAT/ PHIST(NSEGMAX),PHIFIN(NSEGMAX),PHIINC(NSEGMAX),
     +       PHIADD(NSEGMAX),PHISEGA(NSEGMAX),ROTAUTO,PHILAUE,PHIPAD,
     $       PHIZONE,PHIROTAX,CELLSCAL,IFIRSTONE(NSEGMAX),
     $       NSEGM,IPCKCUR,NSTRAT,NUNIQ,
     +       NSTRUN,ISTRUN,NSEGAUTO,
     +       IZONEAX,IROTAX,NLAST,NNPACKS,NLASTPACK,
     +       STRATEGY,AUTO,SIZESET,AUTANOM,SHRUNK,
     +       FIRSTRAT,NEWSTRAT,WAITINP,
     +       OFFPHI
C     ..
C
C


      
C&&*&& end_include  ../inc/strat.f
C&&*&& include  ../inc/sys.f
C
C $Id: sys.f,v 1.1 2002/05/02 10:47:21 harry Exp $
C
C--- awk generated include file  sys.h
C---- START of include file sys.h
C
C
C     .. Scalars in Common Block /SYS/ ..
      INTEGER ISYS
C     ..
C     .. Arrays in Common Block /SYS/  ..
      INTEGER KSYS
C     ..
C     ..  Common Block /SYS/ ..
      COMMON /SYS/ISYS,KSYS(3)
C     ..
C
C
C&&*&& end_include  ../inc/sys.f
C&&*&& include  ../inc/test.f
C--- awk generated include file  test.h
C---- START of include file test.h
C
C    Detector limits, used in SPTEST. All are in 10 micron units.
C    In the virtual detector frame, Y is parallel to the rotation axis
C    and X is orthogonal to the rotation axis.
C
C    In the scanner frame, X is parallel to the slow direction in the
C    image, Y is the fast direction.
C
C    RSCANSQ      Radius**2 of the scanned region of the plate (Mar scanners)
C                 This is different to the other limits because it is
C                 applied to the actual scanner coordinates rather than
C                 the virtual detector coordinates. All other limits apply
C                 to virtual detector coordinates.
C    XSCAN        Maximum X coordinate in detector frame (10 micron units)
C    YSCAN        Maximum Y coordinate in detector frame (10 micron units)
C    RSCANX,RSCANY If non-zero, these are the coordinates (10 micron units)
C                 for the centre of the circle radius RSCAN defining the
C                 useable area of the detector.
C    XSCAN,YSCAN,RSCAN are maximum coords relative to the physical centre of 
C                 the detector (ie the mid-point of the image), NOT wrt the 
C                 direct beam position.
C    RMINX,RMINY If non-zero, these are the coordinates (10 micron units)
C                 for the centre of the circle radius RMIN defining the
C                 useable area of the detector.
C    RMXSQD       Maximum radius (squared) in virtual detector frame
C    RMNSQD       Minimum radius (squared) in virtual detector frame
C    XMIN         Minimum X coordinate in  detector frame
C    XMAX         Maximum X coordinate in  detector frame
C    YMIN         Minimum Y coordinate in  detector frame
C    YMAX         Maximum Y coordinate in  detector frame
C    NXYEXC       Number of excluded areas set by user (rectangles, arcs)
C    NXYEXCDET    Number of detector specific excluded rectangles (Pilatus)
C    NOTILES      If set true, ignore all tile boundaries in Pilatus images
C    ICEWIDTH     Width of the exclusion shells for ice, this fraction
C                 of the actual resolution is used as the half-width
C    ICEWIDTHPREF Width of the exclusion shells for ice for post
C                 refinement, this fraction
C                 of the actual resolution is used as the half-width
C    PREFEXICE    Logical, true if reflections in ice rings to be
C                 excluded from post-refinement
C
C     .. Scalars in common block /TEST/ ..
      REAL RMXSQD,RMNSQD,XMIN,XMAX,YMIN,YMAX,RSCANSQ,RSCANX,RSCANY,
     +     RMINX,RMINY,XSCAN,YSCAN,CIRRAD,ICEWIDTH,ICEWIDTHPREF
      INTEGER NEXCL,NXYEXC,NXYEXCDET,IQUADNO,IARCNO,NEXCL2
      LOGICAL PREFEXICE,NOTILES
C
C     .. Arrays in common block /TEST/
      REAL RESEXL,RESEXH,XYEXC,CNREXC,EXCARC,CIRCEN,PHILIM,RES2EXL,
     +     RES2EXH,RES3EXL,RES3EXH,XYEXCSP
      INTEGER IEXCARC
      LOGICAL LEXCARC
C     ..
C     .. Common Block /TEST/ ..
      COMMON /TEST/RESEXL(MAXEXCL),RESEXH(MAXEXCL),XYEXC(4,MAXEXCL),
     +     CNREXC(8,MAXEXCL),
     $     EXCARC(5,MAXEXCL),CIRCEN(2),PHILIM(2),RMXSQD,RMNSQD,XMIN,
     +     XMAX,YMIN,YMAX,RSCANSQ,RSCANX,RSCANY,RMINX,RMINY,XSCAN,
     $     YSCAN,CIRRAD,IEXCARC(6,MAXEXCL),NEXCL,NXYEXC,IQUADNO,
     $     IARCNO,LEXCARC(MAXEXCL),NXYEXCDET,RES2EXL(MAXEXCL),
     $     RES2EXH(MAXEXCL),NEXCL2,ICEWIDTH,ICEWIDTHPREF,PREFEXICE,
     +     RES3EXL(MAXEXCL),RES3EXH(MAXEXCL),NOTILES,XYEXCSP(4,MAXEXCL)
C     ..
C
C

C&&*&& end_include  ../inc/test.f
C&&*&& include  ../inc/tgen.f
C
C $Id: tgen.f,v 1.2 2014/01/20 09:20:15 andrew Exp $
C
C--- awk generated include file  tgen.h
C---- START of include file tgen.h
C
C---- Stores variables for use in TESTGEN option
C
C     .. Arrays in common /TGEN/ ..
c      REAL XOVER
c      INTEGER ISTATS
C
C     .. Scalars in common /TGEN/ ..
      REAL XOVER,PHSTART,PHEND,PHSTEP,OSCMIN,OSCMAX,SDIVH,SDIVV,
     +     SDELCOR,SDELAMB,SETA,OSCANG,PCMAX,SBLKSIZE
      INTEGER ISTATS
      LOGICAL TESTGEN,TESTRAT
C     ..
C     ..
C     .. Common Block /TGEN/ ..
      COMMON /TGEN/XOVER(2),PHSTART,PHEND,PHSTEP,OSCMIN,OSCMAX,
     +     SDIVH,SDIVV,SDELCOR,SDELAMB,SETA,OSCANG,PCMAX,
     $     SBLKSIZE,ISTATS(MAXPAX,3),TESTGEN,TESTRAT

C&&*&& end_include  ../inc/tgen.f
C&&*&& include  ../inc/tiltlog.f
C
C $Id: tiltlog.f,v 1.1 2002/05/02 10:47:23 harry Exp $
C
C---- awk generated include file tiltlog.h
C---- START of include file tiltlog.h
C
C    for things connected with the new definitions of TILT and TWIST 
C
C     .. Scalars in common block /TILTLOG/ ..
      LOGICAL NUTWIST
      COMMON /TILTLOG/NUTWIST
C&&*&& end_include  ../inc/tiltlog.f
C&&*&& include  ../inc/trigparam.f
C 
C $Id: trigparam.f,v 1.3 2008/07/31 11:46:02 harry Exp $
C
C
C     Trigonometric parameters - set here to highly precise values, don't 
C     calculate in code; values from Clemens Vonrhein
C
C     DTOR.. convert degrees to radians
C     RTOD.. convert radians to degrees
C     PI..   value of PI
      REAL PI
      PARAMETER   ( PI = 3.1415926536 )
      REAL          TWOPI
      PARAMETER   ( TWOPI = 6.2831853072 )
      REAL RTOD
      PARAMETER   ( RTOD = 57.2957795131 )
      REAL          DTOR
      PARAMETER   ( DTOR = 0.0174532925 )
C DOUBLE PRECISION...
      DOUBLE PRECISION PI_DOUBLE
      PARAMETER ( PI_DOUBLE = 3.1415926535897932384626433832795029D0 )
      DOUBLE PRECISION TWOPI_DOUBLE
      PARAMETER ( TWOPI_DOUBLE = 
     $     6.2831853071795864769252867665590058D0 )
      DOUBLE PRECISION RTOD_DOUBLE
      PARAMETER ( RTOD_DOUBLE = 57.2957795130823208768D0 )
      DOUBLE PRECISION DTOR_DOUBLE
      PARAMETER ( DTOR_DOUBLE = 0.01745329251994329576D0 )
c

c
C&&*&& end_include  ../inc/trigparam.f
C&&*&& include  ../inc/xy.f
C
C $Id: xy.f,v 1.1 2002/05/02 10:47:25 harry Exp $
C
C--- awk generated include file  xy.h
C---- START of include file xy.h
C
C     .. Scalars in common block /XY/ ..
      REAL XTOFD,SINV,COSV,TANV,TWOTHETA
      INTEGER ICASS
C     ..
C     .. Common Block /XY/ ..
      COMMON /XY/XTOFD,SINV,COSV,TANV,TWOTHETA,ICASS
C     ..
C
C     XTOFD....  Crystal to detector distance in 10 micron units. Read from
C                keyworded input and never changed.
C
C                Spot positions are calculated in S/R XYSPOT (Called from
C                REEK) and are for an "ideal" detector at a distance of XTOFD.
C                These are converted into pixel positions in S/R MMTOPX
C                which applies the multiplicative factor XTOFRA to allow
C                for refinement of the distance. XTOFRA is the parameter
C                that is actually refined (in RDIST), rather than XTOFD.
C                The refined distance that is printed in the logfile is
C                actually XTOFRA*XTOFD
C
C     ICASS....  Indicates detector type:
C                0     Flat film
C                1     Vee shaped cassettes
C                2     FAST detector (only used in TESTGEN mode of OSCGEN)
C                3     Swung out FAST (ditto)
C                4     IP detector
C     TWOTHETA   Detector swing angle (degrees)
C&&*&& end_include  ../inc/xy.f
C&&*&& include  ../inc/virbat.f
C
C $Id: virbat.f,v 1.2 2003/11/28 12:34:49 harry Exp $
C
C
C----  include file header  for virtual batches in post-refinement
C
       INTEGER NIVB,NVIRBAT,NADDMISSET,IVIRBAT
       LOGICAL NUPR_INT
       COMMON /VIRBAT/ NIVB,NVIRBAT,NADDMISSET,
     +                 IVIRBAT,NUPR_INT
C&&*&& end_include  ../inc/virbat.f
C&&*&& include  ../inc/gui.f
c     gui.h
c     16th April 2002
c     $Id: gui.f,v 1.19 2018/09/19 13:48:26 andrew Exp $
c     
c     This defines the variables and common blocks which are used by the
c     new gui control routine. This will also be used in xdisp and control,
c     to switch the new gui on!
c     
c     

      logical gui_switch, refine_cell(6), refine_detector(11),
     +     autocomplete, autoreturn, newmatstore, sumimgforspots,
     +     sumimgforref
      integer nargs
      integer imageno, nstart, nimgsumdisp
      parameter (nargs = 200)

c     refine_detector logical block index:
c     1: beam x 2: beam y 3: ccomega 4: yscale
c     5: distance 6: tilt 7: twist 8: roff 9: toff
c     
c     refine_cell logical block index:
c     1: a 2: b 3: c
c     4: alpha 5: beta 6: gamma
c     
c     
c     

      integer keep_phi_segments, keep_phi_start(20), 
     +     keep_phi_end(20)

      real storemeanmult

      common /gui/ gui_switch, imageno, refine_cell, refine_detector,
     +     autocomplete, nstart, autoreturn, newmatstore, nimgsumdisp,
     +     sumimgforspots,sumimgforref

      common /keep_phi/ keep_phi_segments, keep_phi_start, 
     +     keep_phi_end, storemeanmult


c     refine_cell(i) is true if you wish to refine cell(i)
C&&*&& end_include  ../inc/gui.f
C&&*&& include  ../inc/spottrap.f
c     spottrap.f
c     global bits and bobs for trapping blank images
c     
c     $Id: spottrap.f,v 1.1 2002/05/20 09:40:31 graeme Exp $
c     

      integer nspt_old, nspt_new, loop_count

      common /spottrap/ nspt_old, nspt_new, loop_count
C&&*&& end_include  ../inc/spottrap.f
C     
C     Extra common blocks for IP
C&&*&& include  ../inc/lmb.f
C
C $Id: lmb.f,v 1.17 2021/03/30 09:16:05 andrew Exp $
C
C--- awk generated include file  lmb.h
C---- START of include file lmb.h
C ODSCAL is to scale image numbers into the range 0-255 for plotting
C average profiles.
C  GAIN should be equal to the overall gain (image counts per X-ray photon)
C of the system, and is used to evaluate standard deviations based on
C counting statistics assuming independent pixels (ie point spread function
C less than pixel size).
C INVERTX true if image is inverted in the slow (X) direction when read in.
C ISCAL is to scale final integrated intensities and sigmas
C IDIVIDE is the adc offset
C ICONST is a constant to be added to all pixels (normally zero) to
C        allow processing of images with zero pixel values in the scanned
C        area.
C IMGP is true if working with image plate data
C SPIRAL is true for scanners with a spiral readout (Mar, DIP2000)
C ORTHOG is true for scanners with orthogonal scan (FUJI, RAXIS, MD)
C NULLPIX is the value of pixels within the image but not in the active
C area of the detector
C TILED is true if there are inactive areas of the detector within the
C       inscribing circle or square.
C
C NTILEX Number of tiles in X direction
C NTILEY Number of tiles in Y direction
C TILEX  X Coordinates of the midpoints of the null areas between tiles
C TILEY  Y Coordinates of the midpoints of the null areas between tiles
C TILEWX Width of the null areas between tiles in X
C TILEWY Width of the null areas between tiles in Y


C MACHINE and MODEL denote the type of detector.
C DECTRIS_TYPE can be PILATUS, EIGER1 or EIGER2. It is used when setting up
C the positions of the gaps between the tiles.
C
C MACHINE is used in the following subroutines:
C GETBLK, GETHDR, INTPXL, INTPXL2, MXDSPL, OPENODS, PROCESS, PUTPXL
C Currently coded types are:
C
C         MACHINE        MODEL               COMMENTS
C         =======        =====               ========
C         MAR            180, 300, 345
C         RAXI(S)        RAXISII, RAXISIV
C         MD
C         FUJI
C         CCD1           CCD1                Princetown CCD at CHESS,Tiff format
C         CCD2           CCD ?               ESRF CCD detector
C         ADSC           QUAD1               ADSC 2x2 CCD detector
C
C TRUSTED whether or not the image header information can be trusted for 
C         automatic processing
C MICROED indicates whether or not this is electron diffraction data,
C        base on the wavelength being less than 0.05A
C
C
C     .. Scalars in Common Block /LMB/ ..
        REAL ODSCAL,GAIN,LOGA,LOGB
        INTEGER ISCAL,IDIVIDE,ICONST,NULLPIX,NTILEX,NTILEY,
     $       RAXISCOMP,IIDIVIDE,COMPRESSION,GAPVAL,IINVERTX
        LOGICAL IMGP,INVERTX,SPIRAL,ORTHOG,CIRCULAR,TILED,SETADC,
     $       NEWRAXIS,BRUKNEW,TRUSTED,MINICBF,HDF5FILE,GAPMASK,
     $       WRONGGAPVALUE,MICROED
        CHARACTER MACHINE*4,MODEL*8,LDETECTOR*20,DECTRIS_TYPE*10,
     +       CMACHINE*40

C     ..
C     .. Arrays in common /LMB/
        INTEGER TILEX,TILEY,TILEWX,TILEWY
C
C     .. Common Block /LMB/ ..
        COMMON /LMB/ ODSCAL,GAIN,LOGA,LOGB,ISCAL,IDIVIDE,ICONST,
     $       NULLPIX,NTILEX,NTILEY,RAXISCOMP,IIDIVIDE,COMPRESSION,
     $       GAPVAL,IINVERTX,TILEX(20),TILEY(20),TILEWX(20),TILEWY(20),
     $       IMGP,INVERTX,SPIRAL,ORTHOG,CIRCULAR,
     $       TILED,SETADC,NEWRAXIS,BRUKNEW,TRUSTED,MINICBF,HDF5FILE,
     $       GAPMASK,WRONGGAPVALUE,MICROED
        COMMON /LMBC/ MACHINE,MODEL,LDETECTOR,DECTRIS_TYPE,CMACHINE

C&&*&& end_include  ../inc/lmb.f
C     INCLUDE 'MOSLIB( disks )'
c     -harvest
C&&*&& include  ../inc/mharvest.f
c      INTEGER MCOLS
C     .. Scalars in common block /CHARVEST/ ..
      CHARACTER HVERSION*100,ProjectName*64,Crystalname*64,
     +          DataSetName*64,HBEAMLINE*10,Precipitant*80,
     +          XNAME_COLS(MCOLS+MAXLATT*3)*64,
     +          DNAME_COLS(MCOLS+MAXLATT*3)*64
C     .. Scalars in common block /IHARVEST/ ..
      LOGICAL USECWD,PNAMEgiven,XNAMEgiven,DNAMEgiven,DOHARVEST
      REAL    PHhar,Hartemp,
     +     HARETA(MAXPAX),HARDIV(MAXPAX),HARDIH(MAXPAX)
      CHARACTER*3 MKDIRMODE,CHMODMODE
      INTEGER KHFLMS
C     ..
C     .. Common block /CHARVEST/ ..
      COMMON /CHARVEST/HVERSION,ProjectName,Crystalname,
     +                 DataSetName,HBEAMLINE,Precipitant,
     +                 XNAME_COLS,DNAME_COLS
C     .. Common block /IHARVEST/ ..
      COMMON /IHARVEST/USECWD,PHhar,Hartemp,KHFLMS,
     +                 PNAMEgiven,XNAMEgiven,DNAMEgiven,DOHARVEST,
     +                 HARETA,HARDIV,HARDIH,MKDIRMODE,CHMODMODE
C     ..
C&&*&& end_include  ../inc/mharvest.f
c     -harvest
      real error
      common /fudge/ error(3)
C     ..
C     .. Equivalences ..
      EQUIVALENCE (IRAS(1),NXS), (IRAS(2),NYS), (IRAS(3),NC)
      EQUIVALENCE (IRAS(4),NRX), (IRAS(5),NRY)
C     ..
      SAVE
C     .. Data statements ..               
      DATA DUMPSTR/'Each reflection',
     +     'Every profile',
     +     '*** DUMP ODS not implemented, use DUMP SPOT ***',
     +     'Reflections with fraction of background pixels rejected > ',
     +     'Profile fitted overloads',
     +     'Reflections containing pixels with value NULLPIX',
     +     'Final profiles only (not summed pixel values)',
     +     'A list of spot coordinates will also be written',
     +     'Fully recorded reflections',
     +     'If number of peak pixels rejected is more than',
     +     'If number of background pixels rejected is more than'/
      DATA VERS1/.FALSE./,VERS2/.FALSE./,VERS3/.FALSE./
      DATA VERSTR/'(V1.0)','(V2.0)','(V3.0)'/
      DATA STRL1/'PARTIAL'/,STRL2/'OVERLOADED'/,
     +     STRL3/'PARTIAL AND OVERLOADED'/
      DATA SABC/'A','B','C','ALPHA','BETA','GAMMA'/
      DATA FIXSTRA/'XCEN','YCEN','OMEGA0','YSCAL','XTOFRA',
     +     'TILT','TWIST','ROFF','TOFF','RDTOFF','RDROFF',
     +     ' ',' ',' '/
      DATA PRINTOP/' Finding fiducials, ',' Refinement, ',' ',
     +     ' Intensity statistics ',' ',' ',' ',' ',' ',' '/
      DATA NCHPR/20,13,1,22,1,1,1,1,1,1/
C     ..
C     .. Data statements ..
C     .. Sizes for red scanner
      DATA XMAXRED/11000/,YMAXRED/11000/,RMAXRED/11000/
C     AL      DATA XMAXIP/9000/,YMAXIP/9000/,RMAXIP/9000/,RSCANIP/9000/
      DATA FXMAX/5800.0/,FYMAX/5800.0/,RMINRED/1000/
      DATA VXMIN/500.0/,VXMAX/8200.0/,VYMAX/5800.0/,VRMAX/8800.0/
C     
C---- Flags for crystal classes
C     
      DATA LCLASS/-1,-1,-1,-1,-1,-1,
     +     -1,-1,-1, 0,-1, 0,
     +     -1,-1,-1, 0, 0, 0,
     +     -1, 1,-1, 0, 0, 0,
     +     -1, 1,-1, 0, 0, 0,
     +     -1, 1,-1, 0, 0, 0,
     +     -1, 1, 1, 0, 0, 0,
     +     -1, 1, 1,-1, 4, 4/
C     AL this was original, it is wrong   +            -1, 0, 0,-1, 4, 4
C     /
C     
C---- Parameters DIVH,DIVV,DELAMB,DELCOR FOR BEAM LINES PX7.2 & PX9.6
C     
      DATA PX72/0.15,0.05,0.0015,0.0/,PX96/0.06,0.02,0.0015,0.0/
      DATA VOLSCAL/1.0/,OVOLSCAL/1.0/
      DATA FRSTWARN/.TRUE./
      DATA SDPSINDEX/.FALSE./,
     $     NOGO/.FALSE./,DONESEG/.FALSE./,HARRY/.FALSE./
      data xsep/0.0/,ysep/0.0/
      DATA EMPTY /MAXEXCL * .FALSE./
      DATA INSEG / 0 /,INFIX / 0 /
      DATA NOTVERTICAL/.FALSE./
      DATA FIRSTPHIP/.TRUE./      
C     ..
C
      THIS = 'CONTROL'
      indent = 1 
      CALLEDFROM = 'CONTROL'
      IF(DEBUG(52))THEN
        WRITE(IOUT,FMT=7185)FIRSTTIME,IFIRSTPACK,NEWGENF,GENOPEN,
     +       RPTFIRST,MODE,CELLSTR,GUI_SWITCH,SOCKLO,COMREAD,
     +       NRLEFT,XMLOPEN,BADKEYWORD,BADLINE
        IF(ONLINE)WRITE(ITOUT,FMT=7185)FIRSTTIME,IFIRSTPACK,NEWGENF,
     +       GENOPEN,RPTFIRST,MODE,CELLSTR,GUI_SWITCH,SOCKLO,COMREAD,
     +       NRLEFT,XMLOPEN,BADKEYWORD,BADLINE
      END IF
 7185 FORMAT(80('*'),/,'At the beginning of CONTROL, the arguments',
     $     ' have these values:',/,'FIRSTTIME: ',L1,' IFIRSTPACK: ',I5,
     $     ' NEWGENF: ',L1,' GENOPEN: ',L1,' RPTFIRST: ',L1,/,
     $     ' MODE: ',I4,/,' CELLSTR: ',A50,/,80('*'),/,1X,'GUI_SWITCH',
     +      L2,' SOCKLO',L2,' COMREAD',L2,' NRLEFT',I3,' XMLOPEN',L2,
     +     ' BADKEYWORD',L2,' BADLINE',L2)
      HEADONLY = .FALSE.
C
C---- IFILESIZE must always be set to 0 when calling control, is only set
C     by a QQINQ in either OPENODS or GETDETECT
C
      IFILESIZE = 0
      DMODE = 1
      IFSTRAT = .FALSE.
      BADLINE = .FALSE.
      BADKEYWORD = .FALSE.
      INVXFORCEREAD = .FALSE.
C     
C---- Special route if call is from MXDSPL
C     
      IF (WINOPEN) CALL MXD_FLU(I)
C     
C---- Still need to do the following initialisations
C 
      IF(.NOT.GOODSEG)THEN
        IF (COMREAD.AND..NOT.DPSDONE)THEN
          CLOSE (UNIT=ICOMM)
        ENDIF
        COMREAD = .FALSE.
      ENDIF
      GOODSEG = .FALSE.
      DELAY = 0
      STOPRUN = .FALSE.
      FORCEREAD = .FALSE.
      INPERR = .FALSE.
      TRAPERR = ((MODE.EQ.3).OR.(MODE.EQ.10))
      READINLINE = .FALSE.
C
C--- Reset MODE of -1 to 1 but set flag to not open new generate file
C
      NOGENOPEN = .FALSE.
      IF (MODE.EQ.-1) THEN
        MODE = 1
        NOGENOPEN = .TRUE.
      END IF
C     
C---- Counts the number of PROCESS keywords for each RUN keyword
C     
      NPROCRUN = 0
      NULINE = .TRUE.
C     
C---- Do jump depending on mode in two steps to keep sgi compiler happy !
C     
      IF ((MODE.GT.0).AND.(MODE.LT.10)) GOTO 32
C     
C     
      IMGKWD = 0

      IPROKWD = 0
      FINDSPOT = .FALSE.
      IF(SDPSINDEX)THEN
        NODISPLAY = .TRUE.
      ELSE
        NODISPLAY = .FALSE.
      ENDIF
C
C---- GUI_SWITCH switches input to NEOCTRL - "long form" input for some 
C     commands related to the new GUI.
C     GUI_SWITCH may be true if SOCKLO is true
C     NODISPLAY should always be true if SOCKLO is true, since we _never_ 
C     want the old gui when using the socket.
C
      IF(GUI_SWITCH.OR.SOCKLO)THEN 
        NODISPLAY = .TRUE.
      ENDIF
      INDNOREF = .FALSE.
      AUTOINDX = .FALSE.
      RESET = .FALSE.
      NOREF = .FALSE.
      RRSET = .FALSE.
      ARRSET = .FALSE.
      DISPSET = .FALSE.
      POWDER = .FALSE.                       
      OTHERS = .TRUE.
      AVPR = .TRUE.
      EXPAND = .FALSE.                               
      PRINTL = .FALSE.
      MAXNREF = NSPOTS
C     
C---- Set NEWGENF FALSE providing the generate file has actually been
C     opened. If using POWDER option to start, the generate file will
C     not yet have actually been opened
C     
      IF (GENOPEN) NEWGENF = .FALSE.
      FIXEDPR = .FALSE.
C     
C---- If repeating a multisegment post-refinement from scratch, save 
C     refined cell
C     
      IF (RPTFIRST) THEN
        DO 2 I = 1,6
          SAVECELL(I) = CELL(I)
 2      CONTINUE
      END IF
C     
C---- Input channel for reading keywords from a command file, default set ITIN
C     to unit 5 if not reading from a command file
C     
      ICOMM = 4
      IF(.NOT.COMREAD)THEN
        ITIN = 5
        ITINS = ITIN
      ENDIF
C     
C---- Initialise some parameters and save cell for stratgey option
C     Whys is this needed if TESTGEN ?
C     This is for STRATEGY runs with more than one PART, as CONTROL is
C     called again for the second part.
C     FIRSTRAT is set TRUE in MOSDATA, and is only set FALSE when the
C     MTZ file
C     has been opened and the symmetry operators obtained.
C     It is reset to TRUE when EXIT is type at the STRATEGY prompt or
C     when
C     End-of-file is reached at STRATEGY prompt.
C     
      IF ((STRATEGY.OR.TESTGEN).AND.(.NOT.FIRSTRAT)) THEN
c     -harvest
        DOHARVEST = .false.
c     -harvest
        ISTRUN = ISTRUN + 1
        IMOSAIC = 2
        IDIVH = 2
        IDIVV = 2
        DO 3 I = 1,6
          SAVECELL(I) = CELL(I)
 3      CONTINUE
      END IF
      IAUTO = 0
      IISTART = 0
      IHKLOUT = 0
chrp 16062008      IGENF = 0
C     
C---- Reset NSER to zero, UNLESS this is a multi seg run
C     
      IF (.NOT.MULTISEG) THEN
        NSER = 0
        NPACK = 0
        ISTARTP = 1
      END IF
C     
C     
C---- See if there are any runs still needing to be done. This happens
C     when
C     more than one PROCESS/SERIAL keyword has been given for a single
C     RUN 
C     keyword,but the images are not abutting in PHI and so must be 
C     processed as separate runs.
C     
      IF (NRLEFT.GT.0) THEN
        NRWORK = NRWORK + 1
        NRLEFT = NRLEFT - 1
        NSER = NSER + 1
C     
C---- restore saved values of SUMPART, ADDPART. These may have been
C     reset
C     (eg if only one image in first SERIAL run, or because SUMPART is
C     set
C     false for the last image being processed in any run.
C     
        SUMPART = SSUMPART
        ADDPART = SADDPART
        POSTREF = SPOSTREF
C     
        IPACKF = NFLEFT(NRWORK)
        IPACKL = NLLEFT(NRWORK)
        ISERADD = ISERLEFT(NRWORK)
        ISERAR(NSER) = ISERADD
        IPACK1A(NSER) = IPACKF
        IPACK2A(NSER) = IPACKL
        PHIRNG = PHILEFT(NRWORK)
        PHISTART = PHISLEFT(NRWORK)
        NPACKS = IPACKL - IPACKF + 1
C     
C---- NPACKS is no. of packs in current serial card whereas NPACK is
C     the number of packs in total
C     
        NPACK = NPACK + NPACKS
        IF (DEBUG(52)) THEN
          WRITE(IOUT,FMT=7180) NPACK,NSER,NPACKS,IFIRSTPACK,ISTARTP,
     +         ISERADD,IPACKF,IPACKL
          IF (ONLINE) WRITE(ITOUT,FMT=7180) NPACK,NSER,NPACKS,
     +         IFIRSTPACK,ISTARTP,ISERADD,IPACKF,IPACKL
 7180     FORMAT(1X,'Control: NPACK=',I3,' NSER=',I2,' NPACKS=',I3,
     +         ' IFIRSTPACK=',I4,' ISTARTP=',I3,' ISERADD=',I6,
     $         ' First image = ',I7,', Last image = ',I7)
        END IF
C     
        IF (NPACK.GT.MAXPAX) THEN
          WRITE(IOUT,FMT=6021) MAXPAX
          WRITE(ITOUT,FMT=6021) MAXPAX
          STOP
 6021     FORMAT(//,1X,'*** FATAL ERROR ***',/,1X,'Maximum number of',
     +         ' images is',I4,/,1X,'To change this change parameter',
     +         ' MAXPAX (PARAMETER (MAXPAX=10000)',/,1X,
     +         ' with a global edit and recompile.')
        ELSE
          J = 0
C     
C---- Note that IDPACK, PHIBEGA, PHIENDA are used in MAIN to set up
C     start and end oscillation angles for image IDPACK
C     
          DO 6 I = ISTARTP,NPACK
            J = J + 1
            IF (I.EQ.ISTARTP) THEN
              IDPACK(I) = IPACKF
              PHIBEGA(I) = PHISTART
            ELSE
              IDPACK(I) = IPACKF + J - 1
              PHIBEGA(I) = ((J-1)*PHIRNG) + PHISTART
            END IF
            PHIENDA(I) = PHIBEGA(I) + PHIRNG

C     AL***** Need to update this for film
            NFPACK(I) = 1
            NFIRST(I) = 1
 6        CONTINUE
C     
C     
C---- Ready for next serial card
C     
          ISTARTP = NPACK + 1
        END IF
C     
C---- Now want to process these images
C     
        KEY = 'RUN '
        GOTO 63
      ELSE
        NRWORK = 0
      END IF
C     
C---- Always set  NIMAG to zero
C     
      NIMAG = 0
      NIMAGES = 0
C     
C---- Always set IBLOCK and IPARTITION to zero
C     
      IBLOCK = 0
      IPARTITION = 0
      NAUTO = 0
      DO 7 I = 1,MAXIMG
        IDENTAUTO(I) = ' '
        TEMPLAUTO(I) = ' '
        PHISET(I) = .FALSE.
 7    CONTINUE
C     
C---- Skip initialisation if not the first time
C     
      IF (.NOT.FIRSTTIME) GO TO 40
C     
      IF (RPTFIRST) NLINE = 1
      MINBATCH = -999
      NRWORK = 0
      NRLEFT = 0
C     
C---- Do not reset NTLINE if this is a repeat of multiseg refinement.
C     
      IF (.NOT.RPTFIRST) NTLINE = 1   
C     
C---- Do not reset NSAVELINE if RPTFIRST true, in case run is aborted
C     in which case NTLINE is reset to NSAVELINE which will be 1.
C     
      IF (.NOT.RPTFIRST) NSAVELINE = NTLINE
      IF (DEBUG(52)) THEN
        WRITE(IOUT,FMT=7670) NTLINE, NSAVELINE
        IF (ONLINE) WRITE(ITOUT,FMT=7670) NTLINE, NSAVELINE
 7670   FORMAT(1X,'In CONTROL, NTLINE=',I5,' NAVELINE=',I5)
      END IF
      NRUN = 0
      NSTRUN = 0
      NSEG = 0
      NFID = 0
      ISTRUN = 0
      ROTATED = .FALSE.
      WAVE = 1.5418
C
C---- Assign a default min radius for spotfinding to be 2mm
C
      RMINSP = 2.0
      RMAXSP = 0.0
      TOR = 0.0
      INRES = 0
      RESLOW = 0.0
      DSTMAX = 0.0
      NEXCL = 0
      DPSEXCL = 0
      NSEGM = 0
      NSEGRD = 0
      DO 5 I = 1,NSEGMAX
        PHIADD(I) = 0.0
 5    CONTINUE
      CELLKEEP = .FALSE.
C     
C---- Default multiplier for postrefinement maximum residual
C     
      RSDMAX = 1.0
C     
C---- Default polarisation for synchrotron sources. This value is for
C     the SRS (UPDATED 29.06.2002 by Harry)
C     
      TORSRS = 0.95
C     
C---- Default maximum reflection width (degrees)
C     
      WMAX = 5.0
      MULTISEG = .FALSE.
      IF (RPTFIRST) MULTISEG = .TRUE.
      MTZNAM = 'HKLOUT'
      IF (.NOT.RPTFIRST) NEWMATNAM = 'NEWMAT'
C     
C---- If repeating a multiseg run, do NOT reinitialise all these
C     
      IF (.NOT.RPTFIRST) THEN
C     YSCALIN = 0.0
        ROFF = 0.0
        TOFF = 0.0
        RDTOFF = 0.0
        RDROFF = 0.0
        ROFFPHI = 0.0
        TOFFPHI = 0.0
        NODES = 0
        NPHI = 0
        TILT = 0
        ITILT = 0
        ITWIST = 0
        CCOM = 0.0
        CCX = 0
        CCY = 0
        CBAR = 0
      END IF
C     
      IGAIN = 0
      ISERADD = 0
      IPIX = 0
      IPIXY = 0
      IYSCAL = 0
      IBEAM = 0
      IRAST = 0
      ISWUNG = 0
      IIDENT = 0
      NHEAD = 0
      ITOL = 0
      IPOLAR = 0
      ITOR = 0
      IBACKS = 0
      INSIZE = 0
      ISIGSET = 0
      IPTMINSET = 0
      IEXTEN = 0
      INMONO = 0
      IMULTI = 0
      IVIRBAT = 0
      IPKRAT = 0
      ILIMIT = 0
      IBGSIG = 0
      IBGFRAC = 0
C     
C---- Set default flags for positional parameter refinement
C     Parameters 9-11 are for MAR IP data only and are set below
C     once the MACHINE type (Mar or RAXIS) is known
C     
C     Refineable parameters are:
C     1 XCEN
C     2 YCEN
C     3 OMEGA0
C     4 YSCAL
C     5 XTOFRA   (Crystal to detector distance multiplier)
C     6 TILT            
C     7 TWIST
C     8 ROFF for Image Plate, BULGE for film
C     9 TOFF
C     10 RDTOFF
C     11 RDROFF
C     
      DO 4 I = 1,NRPAR
        IFIX(I) = 0
        FIXPAR(I) = .FALSE.
        IF (I.GT.7) FIXPAR(I) = .TRUE.
C     
C---- Default is to fix RDROFF,RDTOFF
C     
        IF (I.GE.10) IFIX(I) = 1
 4    CONTINUE
      IFTYPE = 1
      IF (IMGP) THEN
        NHEAD = 1
        MACHINE = 'MAR '
        MODEL = 'M300'
        COMPRESSION = 1
        ODEXT = 'image'
        HDRSIZE = .TRUE.
        USEHDR = .TRUE.
        USETAIL = .FALSE.
        USEDIST = .TRUE.
        USEWAVE = .TRUE.
        USEPHI = .TRUE.
        IFTYPE = 4
        RAST = 0.150
        INVERTX = .TRUE.
        OMEGAFD = 90.0
        omegaf = omegafd * dtor
C     
C---- Default is to fix ROFF for IP data, but if MACHINE is MAR it will
C     be turned on below
C     
        FIXPAR(8) = .TRUE.
      END IF
      WAIT = 0.0
      XLIMIT = 0.0
      NREC = 0
      IYLEN = 0
      NEXTRA = 0
      XMM(1) = 0.0
      YMM(1) = 0.0
      XMAX = 0.0
      YMAX = 0.0
      XMIN = 0.0
      YMIN = 0.0
      RMAX = 0.0
      RMIN = 0.0
      RSCAN = 0.0
      RSCANX = 0.0
      RSCANY = 0.0
      XSCAN = 0.0
      YSCAN = 0.0
      RMINX = 0.0
      RMINXINP = 0.0
      RMINY = 0.0
      N1OD = 0
      G1OD = 0.0
      BASEOD = 0.0
      CURV = 0.0
      XMMF = 0.0
      XMMDB = 0.0
      XOFF = 0.0
      YOFF = 0.0
C     
C---- If repeating a multiseg run, do NOT reinitialise all these
C     
      IF (.NOT.RPTFIRST) THEN
        XTOFRA = -999.0
        YSCAL = -999.0
      END IF
      YSCALIN = 1.0
      GAIN = 1.0
      EFAC = 0.07
      TRUECCOM = 0.0
      IF (IMGP) ONEFILE = .TRUE.
      IF (IMGP) THEN
        ICASS = 4
      ELSE
        ICASS = 0
      END IF
      NSER = 0
      IMISS = 0
      IMISSMAT = 0
      ICELL = 0
      IMAT = 0
      IUMAT = 0
      ISEP = 0
      LSYMM = 0
      NDIR = 0
      IMOSAIC = 0      
      IDIVH = 0
      IDIVV = 0
      IDIST = 0
      IWAVE = 0
      ISYN = 0
      ICUT = 0
      IPRCUT = 0
      INODES = 0
      ISCAN = 0
      ICCX = 0
      ICCY = 0
      IXOFFSET = 0
      IYOFFSET = 0
      ITHSET = 0
      IPIXMIN = 0
      IRMSFAC = 0
C     
C---- pointer to starting paack for serial card
C     
      ISTARTP = 1               
C     
C---- Initialise profile parameters
C     
      NXLINE = 0
      NYLINE = 0
      HIGHRES = .FALSE.
      LOWRES = .FALSE.
      LINESET = .FALSE.
C     
C---- Initialise default processing mode for IP data
C     
C     This is:
C     1) Use ADDPART for partial addition. (ADDPART OFF to suppress)
C     2) Use profile fitting. Cannot override this.
C     3) Use post-refinement, mode "SINGLE" for trigonal and higher 
C     symmetry, otherwise WIDTH 10 degrees. The symmetry is obtained
C     by testing LCELL(2) which will be -ve for symmetry lower than
C     trigonal.
C     (POSTREF OFF to suppres post-refinement)
C     ********* The following not yet implemented **********
C     4) Refine and use an isotropic beam parameter, with added safety
C     margin. (POSTREF USEBEAM OFF to turn of use of refined values).
C     The value used should be obtained through a recursive filter
C     to provide additional stability,eg current value + 0.2*shift
C     5) Do an AUTOMATCH to check orientation and refine mosaic spread
C     for first image. Report error if shift is more than 0.1 degrees
C     (set by keyword ERRLIM). This must somehow be reset if the user
C     really wants to use the AUTOMATCH option to provide refinement.
C     6) 
C     
      IF (IMGP) THEN
        ADDPART = .FALSE.
C
C---- Summing of partials across images is now deprecated. Also the new style
C     post-refinement (uses partials over many images) is incompatible with
C     SUMPART, so change default SUMPART from TRUE to FALSE. 30/9/05
C
        SUMPART = .FALSE.
        POSTREF = .TRUE.
C     
C---- Set flags to indicate if post refinement mode (single or width)
C     has been set, cell parameters fixed or unfixed via keywords.
C     
        PRMODE = .FALSE.
        PRCELL = .FALSE.
        PRNS = 1
      END IF
C     AL     MATCH = .TRUE.
C     AL     RCONV = 0.5
C     AL     USEBEAM = .TRUE.
C     AL   END IF
C     
      DO 10 I = 1,200
        DEBUG(I) = .FALSE.
        NDEBUG(I) = 100
        IF (I.LE.100) THEN
          WARN(I) = .FALSE.
          DO 11 J = 1,MAXIWARN
            IWARN(J,I) = 0
            IF (J.GT.20) GOTO 11
            XWARN(J,I) = 0.0
 11       CONTINUE
          IF (I.LE.50) THEN
            ISERAR(I) = 0
          END IF
          IF (I.LE.30) THEN
            MODS(I) = .FALSE.
            DUMP(I) = .FALSE.
            IARRDUMP(I) = 10
          END IF
        END IF
 10   CONTINUE
C
      DO 12 I = 1,6
        UNFIX(I) = .FALSE.
        FCELL(I) = .FALSE.
 12   CONTINUE
C     
      DUMPSPOT = .FALSE.
      DUMPALL = .FALSE.
      NDSTART = 1
      NDTOT = NREFLS
      IXDMIN = 1
      IXDMAX = 2000
      IYDMIN = 1
      IYDMAX = 2000
C     
      DO 20 I = 1,20
        LPRINT(I) = .FALSE.
        IF (ONLINE .AND. (I.LE.2)) LPRINT(I) = .TRUE.
 20   CONTINUE
C     
C---- Print profiles by default
C     
      LPRINT(11) = .TRUE.
C     
C     
      DO 30 I = 1,MAXPAX
        AVPROF(I) = .TRUE.
        FORCEB(I) = .FALSE.
        FORCEC(I) = .FALSE.
        FILMPLOT(I) = .FALSE.
        ICASSET(I) = 0
 30   CONTINUE
      GENFILE(1:8) = '________'
      NOMEAS = .FALSE.
      ALLOUT = .FALSE.
      PROCES = .FALSE.
      SPOT = .FALSE.
      FIRSTPASS = .FALSE.
      SECONDPASS = .FALSE.
      THIRDPASS = .FALSE.
      PRECESS = .FALSE.
      NPACK = 0
C     
C---- Number of packs to be used in postrefinement
C     
      IF (.NOT.RPTFIRST) NADD = 0
C     
C---- Number of packs in a BLOCK of data (IP data only)
C     
      IF (IMGP) THEN
        NBLOCK = 10
      ELSE
        NBLOCK = 0
      END IF 
C     
C---- Put jumps depending on MODE here
C     
 32   CONTINUE
      IF (MODE.EQ.1) THEN
        GOTO 770
      END IF
      IF ((MODE.EQ.2).OR.(MODE.EQ.4)) THEN
C     
C---- Restore values of ADDPART,SUMPART,POSTREF
C     
        ADDPART = SADDPART
        SUMPART = SSUMPART
        POSTREF = SPOSTREF
        IF (MODE.EQ.2) NSER = 0
C     
C---- If previous run was a single image integration, POSTREF will be
C     FALSE,
C     so must reset this
C     
        IF (MODE.EQ.4) THEN
          POSTREF = .TRUE.
          SPOSTREF = .TRUE.
          SUMPART = .TRUE.
          SSUMPART = .TRUE.
        END IF
        DO 33 I = 1,6
          UNFIX(I) = .FALSE.
 33     CONTINUE
        GOTO 768
      END IF
      IF (MODE.EQ.3) THEN
        GOTO 50
      END IF
C     
C     
C---- IFIRSTPACK is the pack counter for the first pack in this "run"
C     IFIRST is the pack counter for the first pack in
C     current "pack" card there may be more than one
C     "pack" card for each "run"
C     
C---- Assign the first pack for this "run"
C     
 40   CONTINUE
      IFIRSTPACK = NPACK + 1
C     
 50   IF (COMREAD) GO TO 55
C
C---- Set IERR=0 to allow multiple calls to PICKIT 
C     (which sets IERR=50 on exit from IFPWDR) AGWL 20/5/15
C
      IERR = 0
C     
C---- If using menu input, read input lines from the window
C     
      IF ((MODE.EQ.3).OR.(MODE.EQ.10)) THEN
C     
C---- Check for error in input of previous value (if any)
C     
        IF (IOERR.OR.INPERR) THEN
          LINE = ' '
          WRITE(LINE,FMT=7474)
 7474     FORMAT('Error in input, please repeat.')
          CALL MXDWIO(LINE,2)
          INPERR = .FALSE.
        END IF
        LINE = ' '
        WRITE(LINE,FMT=7470)
 7470   FORMAT('MOSFLM => ')
        CALL MXDWIO(LINE,0)
        CALL MXDRIO(LINE)
        NCH = LENSTR(LINE)
        IF (NCH.GT.0) THEN
          WRITE(IOUT,FMT=7472) LINE(1:NCH)
          IF (ONLINE.AND..NOT.SOCKLO) WRITE(ITOUT,FMT=7472)LINE(1:NCH)
        END IF
 7472   FORMAT(1X,'MOSFLM => ',A)
C     
C---- Decode this line.
C     
C     ******************************************
        CALL MPARSE(LINE,IBEG,IEND,ITYP,VALUE,IDEC,NTOK)
C     ******************************************
        IF (NTOK.EQ.0) GOTO 50
        GOTO 61
      END IF
      IF (ONLINE.AND..NOT.SOCKLO) WRITE (ITOUT,FMT=6000)
      IF (BRIEF) WRITE (IBRIEF,FMT=6000)
 6000 FORMAT (1X,'MOSFLM => ',$)
C     
C---- If this is the second or subsequent run of a multisegment 
C     post-refinement, get the control cards from array INLINE instead
C     of reading them from input stream.
C---- Also if repeating a multisegment refinement from the beginning
C     (with
C     a new cell) read the keywords from those stored in INLINE.
C     
C     However, there is a problem if we have exited from IMAGE display
C     mode, because FIRSTTIME is now FALSE, and if we do STRATEGY with
C     more
C     than one PART it will incorrectly try reading from INLINE after
C     the
C     first STRATEGY keyword has been given !
C     
C     Also a problem when two successive POSTREF SEGMENT 1 commands
C     have been given. Again, it tries to read input from INLINE 
C     immediately after encountering the second POSTREF keyword
C     
C     
 55   CONTINUE
      IF ((MULTISEG.AND.(NRUN.LE.NSEG).AND.
     +     (.NOT.FIRSTTIME).AND.(.NOT.WAITINP)).OR.
     +     (MULTISEG.AND.RPTFIRST).OR.
     +     (STRATEGY.AND.(ISTRUN.LT.NSTRUN).AND.(.NOT.FIRSTTIME)
     +     .AND.(MODE.NE.10).AND.(.NOT.WAITINP).and.(ISTRUN.GT.0))) THEN
C     this test was .and..not.socklo, but it seems to be necessary sometimes 
C     to do it for socklo runs.
        READINLINE = .TRUE.
C     
        LINE = INLINE(NLINE)
        NLINE = NLINE + 1
        IF(NLINE.GT.NTOTLINE)NTOTLINE = NLINE
        IF (DEBUG(52)) THEN
          I = LENSTR(LINE)
          IF (I.NE.0) THEN
            WRITE(IOUT,FMT=6549) NLINE-1,LINE(1:LENSTR(LINE))
            IF (ONLINE) WRITE(ITOUT,FMT=6549) NLINE-1,
     +           LINE(1:LENSTR(LINE))
 6549       FORMAT(1X,'Reading line',I3,' from stored input: ',A)
          ELSE
            WRITE(IOUT,FMT=6547) NLINE-1
            IF (ONLINE) WRITE(ITOUT,FMT=6547) NLINE-1
          END IF
        END IF
 6547   FORMAT(1X,'Reading line',I3,' from stored input: ')
        IF (NLINE.GT.1000) THEN
          WRITE(IOUT,FMT=6001)NLINE
          IF (ONLINE) WRITE(ITOUT,FMT=6001)NLINE
 6001     FORMAT(//,1X,'**** ERROR ***',/,1X,'Error in reading ',
     +         'stored lines of input; ',I4,' lines to be read, but',
     $         ' maximum is 1000')
          STOP
        END IF
C     
C---- Decode this line.
C     
C     ******************************************
        CALL MPARSE(LINE,IBEG,IEND,ITYP,VALUE,IDEC,NTOK)
C     ******************************************
C---- Test for zero tokens (as on a comment card)
C     
        IF (NTOK.EQ.0) GOTO 55
C     
C     
C---- If AUTOINDEX was present in command file, substitute it with
C     a MATRIX keyword giving the refined matrix.
C     
        KEY6 = LINE(IBEG(1) :IEND(1))
        CALL CCPUPC(KEY6)
        IF (KEY6.EQ.'AUTOIN') THEN
          NOGO = .TRUE.
          WRITE(LINE,FMT=7500) NEWMATNAM(1:LENSTR(NEWMATNAM))
 7500     FORMAT('MATRIX ',A)
C     ******************************************
          CALL MPARSE(LINE,IBEG,IEND,ITYP,VALUE,IDEC,NTOK)
C     ******************************************
        END IF
C     
C---- If MOSAIC ESTIMATE was present in command file, substitute it
C     with the value determined by the process.
C     
        KEY = LINE(IBEG(1):IEND(1))
        CALL CCPUPC(KEY)
        IF(KEY.EQ.'MOSA')THEN
          KEY = LINE(IBEG(2):IEND(2))
          CALL CCPUPC(KEY)
          IF(KEY.EQ.'ESTI')THEN
            NOGO = .TRUE.
            WRITE(LINE,FMT=7501)2.0*ETA/DTOR
 7501       FORMAT('MOSAIC ',F7.3)
C     ******************************************
            CALL MPARSE(LINE,IBEG,IEND,ITYP,VALUE,IDEC,NTOK)
C     ******************************************
          ENDIF
        ENDIF
        IF(((KEY(1:2).EQ.'GO').OR.(KEY(1:3).EQ.'RUN')).AND.NOGO)THEN
          WRITE(IOUT,FMT=6009)KEY(IBEG(1):IEND(1))
          IF(ONLINE)WRITE(ITOUT,FMT=6009)KEY(IBEG(1):IEND(1))
 6009     FORMAT(1X,'This ',A,' keyword will be ignored')
          GOTO 61
        ENDIF
C     
C---- now write the command lines
C     
        WRITE(IOUT,FMT=7472) LINE(1:LENSTR(LINE))
        IF (ONLINE) WRITE(ITOUT,FMT=7472) LINE(1:LENSTR(LINE))
        IF ((KEY6(1:4).EQ.'CELL').AND.RPTFIRST) THEN
          WRITE(IOUT,FMT=7502)
          IF (ONLINE) WRITE(ITOUT,FMT=7502)
        END IF
        
 7502   FORMAT(1X,'This CELL keyword will be ignored.')
        GOTO 61
      END IF
C     
C---- Read next keyword
C     
 60   CONTINUE
c     socket
      IF(SOCKLO.AND..NOT.XMLOPEN)THEN
C
C---- if previous line was bad input...
C
        IF(BADKEYWORD)THEN
          XMLLINE = ' '
 8500     FORMAT('<?xml version="1.0"?>',
     $         '<!DOCTYPE keyword_error>',
     $         '<keyword_error>',
     $         '<status>',
     $         '<code>ok</code>',
     $         '</status>',
     $         '<keyword>',A,'</keyword>',
     $         '<message>is not recognised</message>',
     $         '</keyword_error>')
          WRITE(XMLLINE,FMT=8500)LINE(IBEG(1):IEND(1))
          CALL WRITE_SOCKET_LENGTH(SERVERFD,
     $         LENSTR(XMLLINE),XMLLINE)
        ENDIF
        IF(BADLINE)THEN
          XMLLINE = ' '
 8520     FORMAT('<?xml version="1.0"?>',
     $         '<!DOCTYPE keyword_error>',
     $         '<keyword_error>',
     $         '<status>',
     $         '<code>ok</code>',
     $         '</status>',
     $         '<keyword>',A,'</keyword>',
     $         '<command>',A,'</command>',
     $         '</keyword_error>')
          WRITE(XMLLINE,FMT=8520)KEY,LINE(1:LENSTR(LINE))
          CALL WRITE_SOCKET_LENGTH(SERVERFD,
     $         LENSTR(XMLLINE),XMLLINE)
        ENDIF
        
c     socket
        LINE200 = '<done>'
c     socket
C
C---- Print for program flow
C  
        CALLED = 'WRITE_SOCKET_LENGTH'
        COMMENT = 'To write <done> '
        IF (DEBUG(99).AND.DEBUG(52)) 
     +              CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
        CALL WRITE_SOCKET_LENGTH(SERVERFD,6,LINE200)
c     socket
        LINE200 = ' '
c     socket
c     using ntok because this will cause eof to be registered
C
C---- Print for program flow
C  
        CALLED = ' *** READ_SOCKET ***'
        COMMENT = 'using NTOK'
        IF (DEBUG(99).AND.DEBUG(52)) 
     +              CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
C
        CALL READ_SOCKET(SERVERFD,LINE400, ntok)
C Socket is broken
        IF(NTOK.EQ.-2)THEN
          NSHUTERR = 54
          CALL SHUTDOWN('CONTROL',NSHUTERR)
        ENDIF
        IF (DEBUG(89)) THEN
          WRITE(IOUT,FMT=7790) LINE400(1:120)
          IF (ONLINE) WRITE(ITOUT,FMT=7790) LINE400(1:120)
 7790     FORMAT(1X,'Socket read in CONTROL1:',A)
        END IF
c     socket
        LINE = LINE400
c     socket
        if(ntok .gt. -1) then
          write(iout,fmt=7472)line(1:LENSTR(LINE))
           CALL MPARSE(LINE,IBEG,IEND,ITYP,VALUE,IDEC,NTOK)
        end if
c     socket
      ELSE
C     ******************************************************
        CALL MPARSER(ITIN,IOUT,LINE,IBEG,IEND,ITYP,VALUE,
     $       IDEC,NTOK)
C     ******************************************************
C     socket
      ENDIF
      BADLINE = .FALSE.
      BADKEYWORD = .FALSE.
C     
C---- eof ?
C     
      IF (NTOK.EQ.-1) THEN
        IF (COMREAD) THEN
          COMREAD = .FALSE.
          ITIN = ITINS
          CLOSE (UNIT=ICOMM)
          GO TO 50
        ELSE
          GO TO 750
        END IF
      ELSE IF (NTOK.EQ.0) THEN
        GOTO 50
      END IF
C     
 61   continue
      IF (COMREAD) WRITE (ITOUT,FMT=7472) LINE(1:MIN(IEND(NTOK),120))
C     
C     
C---- first 4 chars
C     
      KEY = LINE(IBEG(1) :IEND(1))
      KEY6 = LINE(IBEG(1) :IEND(1))
C     
C---- convert to upper case
C     
C     ***********
 63   CALL CCPUPC(KEY)
      CALL CCPUPC(KEY6)
C     ***********
C---- Don't store name of indirect file, as commands will already
C     have been stored.
C     Also, do not store lines if they are being read from INLINE
C     in this call, as things can get out of step !
C     
C     hrp 12022002 NEVER store STRATEGY lines; it will confuse
C     subsequent 
C     postrefinement or integration runs; also don't store RUN or GO
C     which 
C     follows an AUTOINDEX DPS run
C     
      IF(KEY.EQ.'STRA')THEN
         IFSTRAT = .TRUE.
         DPSINDEX = .FALSE.
      ENDIF
      IF ((KEY(1:1).NE.'@').AND.(.NOT.READINLINE)) THEN
        IF(IFSTRAT.OR.(DPSINDEX.AND.((KEY.EQ.'RUN ').OR.
     $       (KEY.EQ.'GO  '))))THEN
          IF(IFSTRAT.AND.(KEY.EQ.'GO  ').OR.
     $         (KEY.EQ.'RUN ').AND..NOT.DPSINDEX)THEN 
            IFSTRAT = .FALSE.
            DOSTRAT = .TRUE.
          ENDIF
        ELSE
C
C---- Problems with new GUI, it sends commands to read Phi values
C     for all images, can be hundreds of them, so get array bounds
C     error. With new GUI, probably do not need to save these lines at all
C     Try this out
C
          IF (.NOT.SOCKLO) THEN
            INLINE(NTLINE) = LINE
            NTLINE = NTLINE + 1
            IF(NTLINE.GT.NTOTLINE)NTOTLINE = NTLINE
            IF (DEBUG(52)) THEN
              WRITE(IOUT,FMT=7650) NTLINE-1,
     +                 INLINE(NTLINE-1)(1:LENSTR(INLINE(NTLINE-1)))
              IF (ONLINE) WRITE(ITOUT,FMT=7650) NTLINE-1,
     +                 INLINE(NTLINE-1)(1:LENSTR(INLINE(NTLINE-1)))
            END IF
          END IF
C     hrp 12022002
C---- FIRSTTIME gets set .false. by strategy
C     
          IF(DOSTRAT)THEN
            FIRSTTIME = .TRUE.
            NRUN = 0
            DOSTRAT = .FALSE.
          ENDIF
        ENDIF
      END IF
 7650 FORMAT(1X,'Stored line',I4,' in CONTROL is: ',A)
C     
C     AL ****** OSCGEN INPUT STARTS HERE
C     AL
C     
C---- TITLe
C     
      IF (KEY.EQ.'TITL') THEN
C     
C---- Must trap a "null" TITLE
C     
        IF (NTOK.GE.2) THEN
          IF (ITYP(2).EQ.3) THEN
C     
C---- Quoted token
C     
            GTITLE = LINE(IBEG(2) :IEND(2))
            NCHAR = IDEC(2)
          ELSE
            GTITLE = LINE(IBEG(2) :IEND(NTOK))
            NCHAR = IEND(NTOK) - IBEG(2) + 1
          END IF
        END IF
c     -harvest
C     
C---- UCWD   use current working directory  Flag only
C     ======
C     if present file is opened in current working directory
C     default is o/p file to
C     $HOME/DepositFiles/ProjectName/DataSetName.ProgramName
C     
      ELSE IF ((KEY.EQ.'UCWD').OR.(KEY6.EQ.'USECWD')) THEN
        USECWD = .true.
C
C---- RSIZ set max width of harvest file line - this is just a place
C     -holder, default is 80 characters, but I don't know what mharvest
C     does...
      ELSE IF (KEY.EQ.'RSIZ')THEN
        CONTINUE
C     
C---- PNAME   PROJECTNAME [pname]  <string>
C     ===========
C     
C     if given with DATASET then harvest will o/p a file
C     this project_name should be always used for the
C     one structure determination
C     no default
C     
      ELSE IF (KEY.EQ.'PNAM') THEN
C     
C---- Must trap a "null" Project name
C     
        IF (NTOK.GE.2) THEN
          IF (ITYP(2).EQ.3) THEN
C     
C---- Quoted token
C     
            PROJECTNAME = LINE(IBEG(2) :IEND(2))
            NCHAR = IDEC(2)
          ELSE
            PROJECTNAME = LINE(IBEG(2) :IEND(NTOK))
            NCHAR = IEND(NTOK) - IBEG(2) + 1
          END IF
        END IF
        PNAMEgiven = .true.
C     
C---- XNAME   CRYSTALNAME [xname]  <string>
C     ===========
C     
C     CRYSTALNAME is written into the MTZ header, along
C     with PROJECTNAME and DATASETNAME. 
C     The pair CRYSTALNAME/DATASETNAME define the data structure.
C     CRYSTALNAME is not currently used in data harvesting. 
C     No default
C     
      ELSE IF (KEY.EQ.'XNAM') THEN
C     
C---- Must trap a "null" crystal name
C     
        IF (NTOK.GE.2) THEN
          IF (ITYP(2).EQ.3) THEN
C     
C---- Quoted token
C     
            CRYSTALNAME = LINE(IBEG(2) :IEND(2))
            NCHAR = IDEC(2)
          ELSE
            CRYSTALNAME = LINE(IBEG(2) :IEND(NTOK))
            NCHAR = IEND(NTOK) - IBEG(2) + 1
          END IF
        END IF
        XNAMEgiven = .true.
C     
C---- DNAME  DATASETNAME  [dname]  <string>
C     ===========
C     
C     if given with PROJECT then harvest will o/p a file
C     this dataset_name is the name of one of the diffraction
C     data sets used in a particular project
C     no default
C     
      ELSE IF (KEY.EQ.'DNAM') THEN
        DATASETNAME = LINE(IBEG(2) :IEND(2))
        DNAMEgiven = .true.
C     
C---- XDETAILS    PH value  TEMP value  PRECIPITANT   string
C     ===========
C     
      ELSE IF (KEY.EQ.'XDET') THEN
        ICOUNT = 1
 8888   ICOUNT = ICOUNT + 1
        IF (ICOUNT.GT.NTOK) GOTO 50
        SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
        CALL CCPUPC(SUBKEY)
        IF (SUBKEY(1:2).EQ.'PHDC') THEN
          ICOUNT = ICOUNT + 1
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
          PHhar = VALUE(ICOUNT)
        ELSE IF (SUBKEY(1:3).EQ.'TKDC') THEN
          ICOUNT = ICOUNT + 1
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
          Hartemp = VALUE(ICOUNT)
        ELSE IF (SUBKEY(1:4).EQ.'CRYS') THEN
          ICOUNT = ICOUNT + 1
          IF (ITYP(ICOUNT).EQ.3) THEN
            PRECIPITANT = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
          ELSE
            PRECIPITANT = LINE(IBEG(ICOUNT) :IEND(NTOK))
          END IF
        END IF
        GOTO 8888
c     -harvest
C     
C---- IDENtifier
C     
      ELSE IF (KEY.EQ.'IDEN') THEN
C     
C---- Must trap a "null" filename
C     
        CALL KWIDEN(IBEG,IEND,ITYP,VALUE,
     $       IDEC,NTOK,NPARM,
     +       INERR,LINE)
        IF (INERR) THEN
          IF (ONLINE) THEN
            IF (COMREAD) THEN
              COMREAD = .FALSE.
              ITIN = ITINS
              CLOSE (UNIT=ICOMM)
            END IF
            BADLINE = .TRUE.
            GOTO 50
          END IF
          NSHUTERR = 1
          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
        END IF

C     
C---- TEMPlate
C     
      ELSE IF (KEY.EQ.'TEMP') THEN

        CALL KWTEMP(IBEG,IEND,ITYP,VALUE,
     $         IDEC,NTOK,NPARM,MODE,
     +         INERR,LINE)
C
C---- Check for NeXus file (.nxs) as mosflm will not currently handle this.
C     It will send a "trapped_error" XML response or STOP once it has processed
C     the "go" keyword.
C
        NCH = LENSTR(TEMPLSTART)
        BADFILE = (TEMPLSTART(NCH-2:NCH).EQ.'nxs')
C
C---- CHUNk
C
      ELSE IF (KEY.EQ.'CHUN')THEN
C     .. value of CHUNK to an integer based on the oscillation width of 
C     .. the image and the mosaic spread, so obviously cannot be done until 
C     .. the images have been read in and a value for the mosaic spread 
C     .. has been obtained. However, once CHUNK has been set by 'AUTO', it 
C     .. not be allowed to change it again - or should it?
        IF(NTOK.NE.2)THEN
          IF(ONLINE)WRITE(ITOUT,FMT=7520)LINE(1:LENSTR(LINE))
          WRITE(IOUT,FMT=7520)LINE(1:LENSTR(LINE))
 7520     FORMAT(//,1X,'*** ERROR ***',/,' CHUNK keyword can only ',
     $         'be used with a single argument, either an ',/,
     $         ' integer or ''AUTO'', which lets Mosflm try to work ',
     $         'out a suitable ',/,' value for itself; ',/,
     $               6X,'Line: ''',A,''' ignored')
        ELSE
          DO ICOUNT = 2,NTOK
            IF(ITYP(ICOUNT).EQ.1)THEN
              CALL CCPUPC(LINE(IBEG(ICOUNT):IEND(ICOUNT)))
              IF(LINE(IBEG(ICOUNT):IBEG(ICOUNT)+3).EQ.'AUTO')THEN
                CCHUNK = .TRUE.
                IF(ONLINE)WRITE(ITOUT,FMT=7522)
                WRITE(IOUT,FMT=7522)
 7522           FORMAT(' Calculate number of images to chunk ',
     $               'together, based on start and end PHI ',/,
     $               ' of the first image and the mosaicity')
              ELSE
                IF(ONLINE)WRITE(ITOUT,FMT=7520)LINE(1:LENSTR(LINE))
                WRITE(IOUT,FMT=7520)LINE(1:LENSTR(LINE))
              ENDIF
            ENDIF
            IF(ITYP(ICOUNT).EQ.2)THEN
              IF(VALUE(ICOUNT).EQ.1.AND.CHUNK.EQ.1)THEN
C---- A "CHUNK 1" command is sent by iMosflm if reading a saved session file
C     where the images have not been "chunked", so this is not an error and 
C     the warning has therefore been commented out. 
C
C                IF(ONLINE)WRITE(ITOUT,FMT=7524)
C                WRITE(IOUT,FMT=7524)
 7524           FORMAT(//,1X,'*** ERROR ***',/,
     $               ' Chunk size is set to 1, but this is the ',
     $               'default anyway so I don''t know',/,
     $               ' why this keyword pair has been given:',/,
     $               6X,'Keyword ''CHUNK'' ignored')
              ELSEIF(VALUE(ICOUNT).LE.0)THEN
                IF(ONLINE)WRITE(ITOUT,FMT=7526)LINE(1:LENSTR(LINE))
                WRITE(IOUT,FMT=7526)LINE(1:LENSTR(LINE))
 7526           FORMAT(//,1X,'*** ERROR ***',/,
     $               ' Chunk size is set to a value which makes ',
     $               'no sense:',/,
     $               6X,'Line: ''',A,''' ignored')
              ELSE
                CHUNK = INT(VALUE(ICOUNT))
                IF(ONLINE)WRITE(ITOUT,FMT=7528)CHUNK
                WRITE(IOUT,FMT=7528)CHUNK
 7528           FORMAT(' Images will be grouped in blocks of ',I2,
     $               ' virtual images for processing')
              ENDIF
            ENDIF
          ENDDO
          ICHUNK = 1
        ENDIF
        
C     
C---- NEWMAT filename for output matrix from postref segment
C     
      ELSE IF (KEY.EQ.'NEWM') THEN
C     
C---- Must trap a "null" filename
C     
        IF (NTOK.EQ.1) THEN
          WRITE(IOUT,FMT=7060) LINE(IBEG(1):IEND(1))
          IF (ONLINE) WRITE(ITOUT,FMT=7060) LINE(IBEG(1):IEND(1))
          GOTO 50
        END IF
        IF (ITYP(2).EQ.1) THEN
          NEWMATNAM = LINE(IBEG(2) :IEND(2))
          INEWMAT = 2
        END IF
C
C---- WMAT matrix name; write out the current matrix to a file
C
      ELSE IF (KEY.EQ.'WMAT') THEN
        IF (NTOK.EQ.1) THEN
          WRITE(IOUT,FMT=7060) LINE(IBEG(1):IEND(1))
          IF (ONLINE) WRITE(ITOUT,FMT=7060) LINE(IBEG(1):IEND(1))
          GOTO 50
        END IF
        IF (ITYP(2).EQ.1) THEN
          NEWMATNAM = LINE(IBEG(2) :IEND(2))
          INEWMAT = 2
          CALL CCPDPN (3,NEWMATNAM,'UNKNOWN','F',80,IFAIL)
 7099     FORMAT(3(3F12.8/),(3F12.3/),3(3F12.8/),6F12.4,/,
     $         3F12.4,/,'SYMM ',A)
          WRITE(3,FMT=7099) ((AMAT(I,J),J=1,3),I=1,3),
     +         (DELPHI(J),J=1,3),
     +         ((UMAT(I,J),J=1,3),I=1,3),CELL,
     +         (DELPHI(J),J=1,3),SPGNAM
          CLOSE (UNIT=3)
        ENDIF
C     
C---- MATRix x x x x x x x x x
C     
      ELSE IF (KEY.EQ.'MATR') THEN
          CALL KWMATR(IBEG,IEND,ITYP,VALUE,
     $       IDEC,NTOK,NPARM,RPTFIRST,FIRSTTIME,SYMMIN,
     +       INERR,LINE)
        IF (INERR) THEN
          IF (ONLINE) THEN
            IF (COMREAD) THEN
              COMREAD = .FALSE.
              ITIN = ITINS
              CLOSE (UNIT=ICOMM)
            END IF
            BADLINE = .TRUE.
            GOTO 50
          END IF
          NSHUTERR = 3
          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
        END IF
        EXPLICIT_CELL = .TRUE.
C
C---- SUBMATRIX <lattice_number> <solution_number> <amatrix>
C
      ELSE IF (KEY.EQ.'SUBM') THEN
C
C---- No checks for errors (only sent by imosflm)
C
C     
C     ************************************
          CALL MKEYNM(2,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          ILATT = NINT(VALUE(2))
          ILATTSOL = NINT(VALUE(3))
C
C---- Need to assign BESTSOLAR if processing is done in
C     batch without any indexing, as BESTSOLARR is used in GENERATE
C
          IF (BESTSOLAR(ILATT).EQ.0) BESTSOLAR(ILATT)=ILATTSOL
C
C---- Need to set NLATTICE in case of batch integration with no indexing done
C
          NLATTICE = MAX(NLATTICE,ILATT)
C     ************************************
          CALL MKEYNM(9,4,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C
          N = 3     
          DO I = 1,3
            DO J = 1,3
              N = N + 1
              CRYSTALA(I,J,ILATT) = VALUE(N)
              CRYSTALUBREF(I,J,ILATTSOL,ILATT) = VALUE(N)
            ENDDO
          ENDDO        
C     
C---- TARGET matrix file
C     
      ELSE IF (KEY.EQ.'TARG') THEN
C     
C---- Must trap a "null" filename
C     
        IF (NTOK.EQ.1) THEN
          WRITE(IOUT,FMT=7060) LINE(IBEG(1):IEND(1))
          IF (ONLINE) WRITE(ITOUT,FMT=7060) LINE(IBEG(1):IEND(1))
          GOTO 50
        END IF
C     
C---- Check if second token is alphanumeric, if so it should be:
C     EITHER
C     1)the file written by IDXREF containing AMAT and missetting
C     angles
C     OR
C     2) LATTICE <NLATT> SOLUTION <NSOL> MATRIX <9 values of AMAT>
C     
C     
C     
        IF (ITYP(2).EQ.1) THEN
          TARFILE = LINE(IBEG(2) :IEND(2))
          KEY8 = TARFILE(1:8)
          CALL CCPUPC(KEY8)
C
C---- If this is LATTICE, then this is a special case for multiple
C     lattices processed with iMosflm. 
C     The passed amatrix must be used as a target for those solutions
C     with the same Bravais lattice type in other lattices, providing 
C     they are the "suggested" solution.
C
          IF (KEY8(1:7).EQ.'LATTICE') THEN


C---- No further checks for subkeywords as this command only sent by
C     imosflm.
C
C     ************************************
            CALL MKEYNM(1,3,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
            ILATT = NINT(VALUE(3))
C     
C     ************************************
            CALL MKEYNM(1,5,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
            ILATTSOL = NINT(VALUE(5))
            WRITE(IOUT,FMT=7820) ILATT, ILATTSOL
            IF (ONLINE) WRITE(ITOUT,FMT=7820) ILATT, ILATTSOL
 7820       FORMAT(/,1X,'target matrix',
     +               ' assigned to lattice',I3,' solution',I3)
C
C     ************************************
            CALL MKEYNM(9,7,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
            N = 6     
            DO I = 1,3
              DO J = 1,3
                N = N + 1
                TMAT(I,J) = VALUE(N)
              ENDDO
            ENDDO
            CALL REVAMAT(ILATT,ILATTSOL,TMAT)
            GOTO 50        
          END IF
C
C---- Test that the MATRIX keyword has already been given
C
          IF (IMAT.NE.1) THEN
            WRITE(IOUT,FMT=7510)
            IF (ONLINE) WRITE(ITOUT,FMT=7510)
 7510       FORMAT(/,/,1X,'***** FATAL ERROR *****',/,1X,
     +         'the MATRIX keyword MUST be given BEFORE the TARGET',
     +         ' keyword.')
            STOP
          END IF
C
          IFAIL = 1  
C     
C     **************************************
          IF (BRIEF) THEN
            CALL CCPDPN (-3,TARFILE,'OLD','F',80,IFAIL)
          ELSE
            CALL CCPDPN (3,TARFILE,'OLD','F',80,IFAIL)
          END IF
C     **************************************
C     
C---- Trap file open failure
C     
          IF (IFAIL.LT.0) THEN
            WRITE(IOUT,FMT=7371)
            IF (ONLINE) WRITE(ITOUT,FMT=7371)
 7371       FORMAT(1X,'**** ERROR ****',/,1X,
     +           'Cannot open TARGET matrix file.')
            IMAT = 0
            IF (ONLINE) THEN
              IF (COMREAD) THEN
                COMREAD = .FALSE.
                ITIN = ITINS
                CLOSE (UNIT=ICOMM)
              END IF
              BADLINE = .TRUE.
              GOTO 50
            END IF
            NSHUTERR = 6
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF
          READ (3,FMT=6004,END=79) ((TMAT(I,J),J=1,3),I=1,3),
     +         (TARPHI(I),I=1,3)
          READ(3,FMT=6004,END=79) ((TJUNK(I,J),J=1,3),I=1,3)
          READ(3,FMT=6003,END=79) TARCELL
          CLOSE (UNIT=3)
C     
C---- Check the missets are zero
C     
          IF ((TARPHI(1).NE.0).OR.(TARPHI(2).NE.0).OR.
     +         (TARPHI(3).NE.0)) THEN
            WRITE(IOUT,FMT=7516)
            IF (ONLINE) WRITE(ITOUT,FMT=7516)
 7516       FORMAT(/,/,1X,'***** FATAL ERROR *****',/,1X,
     +           'The missets in the target matrix must be zero')
            STOP
          END IF
          CALL TARGMAT(AMAT,TMAT,TARCELL)
          WRITE(IOUT,FMT=7518) ((TMAT(I,J),J=1,3),I=1,3),
     +         ((AMAT(I,J),J=1,3),I=1,3)
          IF (ONLINE) WRITE(ITOUT,FMT=7518) 
     +         ((TMAT(I,J),J=1,3),I=1,3),((AMAT(I,J),J=1,3),I=1,3)
 7518     FORMAT(1X,'Input target matrix',3(1X,3F10.6/),/,
     +         1X,'AMAT after transformation',3(1X,3F10.6/))
        ELSE
          WRITE(IOUT,FMT=7512)
          IF (ONLINE) WRITE(ITOUT,FMT=7512)
 7512     FORMAT(/,/,1X,'***** FATAL ERROR *****',/,1X,
     +         'The name of the file for the TARGET matrix must',
     +         ' be given.')
          STOP
        END IF
        GOTO 50
C     
C     End of file error
C     
 79     WRITE(IOUT,FMT=7514) TARFILE
        IF (ONLINE) WRITE(ITOUT,FMT=7514) TARFILE
 7514   FORMAT (//,1X,'******* End of file when reading matrix and ',
     +       'missetting angles from ',A,'  *******')

C     
C---- UMAT x x x x x x x x x
C     
      ELSE IF (KEY.EQ.'UMAT') THEN
C     
C---- Must trap a "null" filename
C     
        IF (NTOK.EQ.1) THEN
          WRITE(IOUT,FMT=7060) LINE(IBEG(1):IEND(1))
          IF (ONLINE) WRITE(ITOUT,FMT=7060) LINE(IBEG(1):IEND(1))
          BADLINE = .TRUE.
          GOTO 50
        END IF
C     
C---- First check that this is a second (or greater) segemnt of a
C     multiple segment post-refinement run, when MATRIX must be given
C     rather than UMAT
C     
        IF (NSEG.GT.1) THEN
          WRITE(IOUT,FMT=7045)
          IF (ONLINE) WRITE(ITOUT,FMT=7045)
 7045     FORMAT(1X,'**** ERROR ****',/,1X,'Must give a MATRIX',
     +         ' keyword rather than UMAT for multiple ',
     +         'serial keyword runs')
          NSHUTERR = 7
          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
        END IF
C     
C---- Check if second token is alphanumeric, if so it should
C     be the file written by REFIX or IDXREF containing AMAT and
C     missetting
C     angles followed by UMAT, CELL and missetting angles
C     
        IUMAT = 1
        IF (ITYP(2).EQ.1) THEN
          IDXFILE = LINE(IBEG(2) :IEND(NTOK))
          IFAIL = 0      
C     
C     **************************************
          IF (BRIEF) THEN
            CALL CCPDPN (-3,IDXFILE,'OLD','F',80,IFAIL)
          ELSE
            CALL CCPDPN (3,IDXFILE,'OLD','F',80,IFAIL)
          END IF
C     **************************************
C     
          READ (3,FMT=6004,END=83) ((UMAT(I,J),J=1,3),I=1,3),
     +         (TDELPHI(I),I=1,3),
     +         ((UMAT(I,J),J=1,3),I=1,3)
          READ (3,FMT=6003,END=83) (TCELL(I),I=1,6)
          READ (3,FMT=6004,END=83) (TDELPHI(I),I=1,3)
 6004     FORMAT (3F12.6)
 6003     FORMAT (6F12.6)
C     
C---- Put missets in range -180 to 180, transfer missets and cell
C     to proper locations providing these have not been specified on
C     MISSETS or CELL cards
C     
          DO 82 I=1,6
            IF (ICELL.EQ.0) CELL(I) = TCELL(I)
            IF (I.LE.3) THEN
              IF (IMISS.EQ.0) DELPHI(I) = TDELPHI(I)
              IF (DELPHI(I).GT.180.00) DELPHI(I) = 
     $             DELPHI(I) - 360.00
            ENDIF
 82       CONTINUE
C     
C     
C     Check if cell parameters have already been specified
          IF ((ICELL.EQ.1).AND.(.NOT.RPTFIRST)) THEN
            WRITE(IOUT,FMT=7011) 
            IF (ONLINE) WRITE(ITOUT,FMT=7011) 
 7011       FORMAT (' ** BEWARE ** Cell parameters read from CELL',
     +           ' card will overwrite those read',/,1X,
     $           'from file')
          END IF
          ICELL = 1
C     
C     Check if missetts have been read from MISSET card
C     
          IF (IMISS.EQ.1) THEN
            WRITE(IOUT,FMT=7013) 
            IF (ONLINE) WRITE(ITOUT,FMT=7013) 
 7013       FORMAT (' ** BEWARE ** missetting angles read from ',
     +           'MISSET  card will overwrite those read',/,1X,
     $           'from the file')
          END IF
          IMISSMAT = 1
          CLOSE (UNIT=3)
        ELSE
C     
C     
C     ************************************
          CALL MKEYNM(9,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
          N = 1
C     
C     
          DO 100 I = 1,3
            DO 90 J = 1,3
              N = N + 1
              UMAT(I,J) = VALUE(N)
 90         CONTINUE
 100      CONTINUE
C     
        END IF
        GOTO 50
C     
C---- End of file in read
C     
 83     WRITE(IOUT,FMT=7514) IDXFILE
        IUMAT = 0
        IF (ONLINE) THEN
          WRITE(ITOUT,FMT=7514) IDXFILE
          IF (COMREAD) THEN
            COMREAD = .FALSE.
            ITIN = ITINS
            CLOSE (UNIT=ICOMM)
          END IF
          BADLINE = .TRUE.
          GOTO 50
        END IF
        NSHUTERR = 8
        CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
C     
C---- ANGLes x x x (RESET)
C     or MISSet x x x (RESET)
C     
      ELSE IF(KEY.EQ.'ANGL'.OR.KEY.EQ.'MISS') THEN
C     
C     ************************************
        CALL MKEYNM(3,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
C     Check if missetts have been read from UMAT or AMAT file
        IF (IMISSMAT.EQ.1) THEN
          WRITE(IOUT,FMT=7015) 
          IF (ONLINE) WRITE(ITOUT,FMT=7015) 
 7015     FORMAT (' ** BEWARE ** These missetting angles will ',
     +         'overwrite those read from the UMAT',/,1X,
     $         'or MATRIX file')
        END IF
        IMISS = 1
C     
        DO 110 I = 1,3
          DELPHI(I) = VALUE(I+1)
C     
C---- Put missets in range -180 to 180
C     
          IF (DELPHI(I).GT.180.00) DELPHI(I) = DELPHI(I) - 360.00
 110    CONTINUE
C     
C---- Check if these are to be incorporated into the U matrix
C     
        IF (NTOK.EQ.5) THEN
          IF(LINE(IBEG(5):IEND(5)).EQ.'RESET')
     +       RESET = .TRUE.
        ENDIF
C     
C---- DISTance x
C     
      ELSE IF (KEY.EQ.'DIST'.AND.(ITYP(2).EQ.2)) THEN
C     
C     ************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
        XTOFD = VALUE(2)*100.0
        CAMCON(3) = VALUE(2)
        CALL SETDIS(ITILT,ITWIST,1)
C     RADEG = 18000.0/3.14159
C     IF (XTOFD.NE.0) FDIST = 1.0/ (XTOFD*RADEG)
C     TILT = ITILT*FDIST
C     TWIST = ITWIST*FDIST
C     
C---- Input in mms - store in 10 microns
C     
        IDIST = 2
C     
C---- CRYStal...only used to specify use of rhombohedral rather than
C     hexagonal
C     setting for rhombohedral spacegroups
C     
      ELSE IF (KEY.EQ.'CRYS') THEN
C     
C---- The SYMMETRY keyword MUST be given AFTER the CRYST keyword
C     
        IF (LSYMM.GT.0) THEN
          WRITE(IOUT,FMT=7460)
          IF (ONLINE) WRITE(ITOUT,FMT=7460)
 7460     FORMAT(/,1X,'***** FATAL ERROR *****',/,1X,'The CRYST ',
     +         'keyword MUST be given BEFORE the SYMMETRY keyword.')
          NSHUTERR = 9
          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
        END IF
        I = 2
        IF (I.LE.NTOK) THEN
          LATTYP = LINE(IBEG(I) :IEND(I))
C     
C     **********
          CALL CCPUPC(LATTYP)
C     **********
C     
C---- CRYStal RHOMbohedral
C     
          IF (LATTYP.EQ.'R') THEN
            ICRYST = 8
C     
C---- Set cell refinement flags
C     
            DO 112 I = 1,6
              LCELL(I) = LCLASS(I,ICRYST)
 112        CONTINUE
            WRITE(IOUT,7462)
            IF (ONLINE) WRITE(ITOUT,7462)
 7462       FORMAT(/,1X,'***** IMPORTANT *****',/,1X,'A new entry',
     +           ' must be given in "syminfo.lib" for a rhombohedra',
     +           'l cell',/,1X,'if not using the hexagonal settin',
     +           'g.',/,1X,'eg 1146 for R3 using the rhombohedral',
     +           ' cell.',/,1X,'Then give the symmetry as: SYMMET',
     +           'RY 1146.'/,1X,'If autoindexing without supplyin',
     +           'g the cell and spacegroup, select the triclinic',
     $           /,1X,'solution that is',
     +           ' closest to the rhombohedral cell, but give the',
     +           ' spacegroup as 1146.')
          ELSE
            WRITE(IOUT,FMT=6470)
            IF (ONLINE) WRITE(ITOUT,FMT=6470)
          END IF
 6470     FORMAT(/,1X,'**** ERROR *****',/,1X,'Only allowed CRYSTAL',
     +         ' class is RHOMBOHEDRAL which specifies use',
     +         /,1X,'of r hombohedral rather than hexagonal cell',
     +         ' for rhombohedral spacegroups',/,1X,'Use SYMMETRY',
     +         ' keyword to specify the crystal symmetry.')
        END IF
C     
C     
C---- WAVElength x 
C     
      ELSE IF (KEY.EQ.'WAVE' ) THEN
C     
C     ************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
        IIWAVE = .TRUE.
        WAVE = VALUE(2)
C     
C---- Wavelength in Angstrom
C     
        IWAVE = 2
C
C---- trap for silly values of high resolution edge
C
        IF(RES.LT.WAVE/2.0.AND.INRES.EQ.1)THEN
 8000     FORMAT(/,/,' **** WARNING **** ',/,1X,'With the current ',
     $         'wavelength of ',F8.6,'A the maximum resolution ',
     $         'given ',/,1X,'(',F8.6,'A) is too high. It will be ',
     $         'set to the maximum possible',/,1X,'resolution for ',
     $         'this wavelength of ',F8.6,/,/)
          WRITE(IOUT,FMT=8000)WAVE,RES,WAVE/2.0
          IF(ONLINE)WRITE(ITOUT,FMT=8000)WAVE,RES,WAVE/2.0
          RES = WAVE/2.0
          
        ENDIF
        
C     
C---- DSTArmax
C     
      ELSE IF (KEY.EQ.'DSTA') THEN
C     
C     ************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
        DSTMAX = VALUE(2)
C     
C---- RESOlution
C     
      ELSE IF (KEY.EQ.'RESO') THEN
c
C---- *** IMPORTANT ***
C     This keyword is also read by KWRESO so if any changes are made 
C     here, they need to be made in KWRESO as well.
C     
C---- Check for subkeyword
C   
        ICOUNT = 1
C
C---- check for 'CHEC' subkeyword, plus co-ordinates
C
        SUBKEY = LINE(IBEG(2):IEND(2))
        CALL CCPUPC(SUBKEY)
        IF((NTOK.EQ.4).AND.SUBKEY(1:4).EQ.'CHEC')THEN
          II = VALUE(3)/RAST
          IF(INVERTX)II=NREC-(VALUE(3)/RAST)
          JJ = VALUE(4)/RAST
          CALL CRESOL(II,JJ,XTOFD*0.01,TWOTHETA,WAVE,THISRES)
          write (xmlline, fmt=8600)thisres
 8600     format('<?xml version="1.0"?>',
     $         '<!DOCTYPE resolution_response>',
     $         '<resolution_response>',
     $         '<status>',
     $         '<code>ok</code>',
     $         '</status>',
     $         '<resolution>',F12.4,'</resolution>',
     $         '</resolution_response>')
          CALL WRITE_SOCKET_LENGTH(SERVERFD,
     $         LENSTR(XMLLINE),XMLLINE)
          goto 50
        endif
        IIRES = .TRUE.
C     
C---- Set DSTMAX to zero so that it does not use a previously defined
C     lower resolution limit
C     
        DSTMAX = 0.0
        IF (ITYP(2).EQ.1) GOTO 152
        ANITES = .FALSE.
C     
C---- Only high resolution given
C     
        IF (((NTOK.EQ.2).AND.(ITYP(2).EQ.2)).OR.
     +       ((NTOK.GT.2).AND.(ITYP(3).NE.2))) THEN
C     ************************************
          CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          IF(ABS(VALUE(2)).GT.1E-2)THEN
            INRES = 1
            RES = VALUE(2) 
          ELSE
            INRES = 0
            RES = 0.0
          ENDIF
          ICOUNT = ICOUNT + 1
C     
C---- Both high and low resolution given
C     
        ELSE
C     ************************************
          CALL MKEYNM(2,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          INRES = 1
          RESLOW = VALUE(2)
          RES = VALUE(3)
          IF (RES.GT.RESLOW) THEN
            TEMPX = RES
            RES = RESLOW
            RESLOW = TEMPX
          END IF
          IF(RES.LT.1E-2)THEN
            INRES = 0
            RES = 0.0
          ENDIF
          ICOUNT = ICOUNT + 2
        END IF
C
C---- trap for silly values of high resolution edge
C
        IF(WAVE.GT.1E-3.AND.RES.LT.WAVE/2.0.AND.INRES.EQ.1)THEN
          WRITE(IOUT,FMT=8000)WAVE,RES,WAVE/2.0
          IF(ONLINE)WRITE(ITOUT,FMT=8000)WAVE,RES,WAVE/2.0
          RES = MAX (RES, WAVE/2.0)
        ENDIF

        IF (ICOUNT.EQ.NTOK) GOTO 156
C     
C---- Other subkeywords given
C     
 152    ICOUNT = ICOUNT + 1
        SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     **************
        CALL CCPUPC(SUBKEY)
C     **************
        IF (SUBKEY.EQ.'CUTO') THEN
          ICOUNT = ICOUNT + 1
C     ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          RESCUT = VALUE(ICOUNT)
        ELSE IF (SUBKEY.EQ.'EXCL') THEN
          ICOUNT = ICOUNT + 1
          IF(ITYP(ICOUNT).EQ.1)THEN
C
C----     Next token is a string
C
            SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
            CALL CCPUPC(SUBKEY)
            IF(SUBKEY.EQ.'NONE')THEN
              NEXCL = 0
              DO 7145 I=1,MAXEXCL
                RESEXL(I) = 0.0
                RESEXH(I) = 0.0
 7145         ENDDO
            ELSEIF(SUBKEY(1:3).EQ.'ICE')THEN
              IF(ITYP(ICOUNT+1).EQ.2)THEN
                EXCLUDERING = VALUE(ICOUNT+1)
                ICOUNT = ICOUNT + 1
              ELSE
                EXCLUDERING = 13
              ENDIF
CcHRP corrected 27.01.2011 - listed strong to weak (roughly). Top 5 are much
C     stronger than the rest. Then four moderately strong. XPAD can be used
C     to modify the width of the resolution shells.
C
              XPAD = 0.005
              WRITE(IOUT,FMT=9876) EXCLUDERING
              IF (ONLINE) WRITE(ITOUT,FMT=9876) EXCLUDERING
 9876         FORMAT(1X,'CALLING SETUPICE EXCLUDERING IS',I3)
              CALL SETUPICE(RESEXL,RESEXH,EXCLUDERING,NEXCL,XPAD)
            ENDIF
          ELSE
C
C----     Next token is a number
C
            NEXCL = NEXCL + 1
            IF (NEXCL.GT.MAXEXCL) THEN
              WRITE(IOUT,FMT=7150) MAXEXCL
              IF (ONLINE) WRITE(ITOUT,FMT=7150) MAXEXCL
 7150         FORMAT(1X,'*** ERROR ***',/,1X,'A maximum of',I3,
     +             ' resolution ranges can be excluded.',/,1X,'The ',
     +             'remainder will be ignored')
              GOTO 156
            END IF
C
C---- Check if excluded resolution RANGE has been given with
C     an intervening "TO"
C
            IF(ITYP(ICOUNT+1).EQ.1) THEN
C
C---- There is a "TO", read the two numbers separately
C     ************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
              RESEX1 = VALUE(ICOUNT)
              ICOUNT = ICOUNT + 2
C     ************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
              RESEX2 = VALUE(ICOUNT)
            ELSE            
C
C---- There is not a "TO", read the two numbers at same time
C     ************************************
              CALL MKEYNM(2,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
              RESEX1 = VALUE(ICOUNT)
              RESEX2 = VALUE(ICOUNT+1)
              ICOUNT = ICOUNT + 1
            END IF
            IF (RESEX2.GT.RESEX1) THEN
              TEMPX = RESEX1
              RESEX1 = RESEX2
              RESEX2 = TEMPX
            END IF
            RESEXL(NEXCL) = RESEX1
            RESEXH(NEXCL) = RESEX2
            IF (ICOUNT.GE.NTOK) THEN
              WRITE(IOUT,FMT=7146) (RESEXH(J),RESEXL(J),J=1,NEXCL)
              IF (ONLINE) WRITE(ITOUT,FMT=7146) (RESEXH(J),RESEXL(J),
     +                                J=1,NEXCL)
 7146         FORMAT(/,1X,'**** Information ****',/,1X,'The following '
     $             ' resolution ranges will be excluded:',/,
     $             (1x,6(F6.3,' to',F6.3,',   ')))
            ENDIF
          ENDIF
C     
C---- Anistropic resolution limits for a*, b*, c*.
C     
        ELSE IF (SUBKEY.EQ.'ANIS') THEN
          INRES = 1
          ICOUNT = ICOUNT + 1
C     
C---- must supply three numbers
C     
C     ************************************
          CALL MKEYNM(3,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          RESANI(1) = VALUE(ICOUNT)
          RESANI(2) = VALUE(ICOUNT+1)
          RESANI(3) = VALUE(ICOUNT+2)
          RES = MIN(RESANI(1),RESANI(2))
          RES = MIN(RES,RESANI(3))
          ANITES = .TRUE.
          ICOUNT = ICOUNT + 2
C     
C---- We might want a low resolution limit if anisotropic resolution
C     limits set.
C     
        ELSE IF (SUBKEY.EQ.'LOW') THEN
C          INRES = 1
          ICOUNT = ICOUNT + 1
C     
C---- should supply three numbers, but we'll assume low resolution
C     limit is isotropic
C     
C     ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          RESLOW = VALUE(ICOUNT)
C
C---- if the user (or the new GUI) has given both low and high resolution limits
C     on separate input lines, we need to make sure they are the right way round!
C
          IF (RES.GT.RESLOW) THEN
            TEMPX = RES
            RES = RESLOW
            RESLOW = TEMPX
          END IF
C     
C---- else the subkey is wrong; tell the user
C     
        ELSE
          WRITE (IOUT,FMT=6130) SUBKEY
          IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
        END IF
        IF (ICOUNT.LT.NTOK) GOTO 152
 156    CONTINUE
C     
C---- SYNChrotron
C     
      ELSE IF (KEY.EQ.'SYNC') THEN
        INMONO = 1
        ISYN = 1
        IMONO = 2
        TOR = TORSRS
        I = 1
 150    CONTINUE
        I = I + 1
        IF (I.GT.NTOK) THEN
          GO TO 170
        ELSE
C     
C---- Skip if no more tokens on line
C     
          KEY2 = LINE(IBEG(I) :IEND(I))
C     
C     **********
          CALL CCPUPC(KEY2)
C     **********
C     
C---- SYNC PX72
C     
          IF (KEY2.EQ.'PX72') THEN
            Hbeamline='DLAB-PX7.2'
            GO TO 160
C     
C---- SYNC PX96
C     
          ELSE IF (KEY2.NE.'PX96') THEN
C     
C---- SYNC DIVH
C     
            IF (KEY2.EQ.'DIVH') THEN
              I = I + 1
C     
C     ************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
              DIVHD = VALUE(I)
              DIVH = 0.5*DTOR*DIVHD
              IDIVH = 1
              WARN(26) = .FALSE.
C     
C---- SYNC POLARISATION...defines degree of polarisation of beam
C     
            ELSE IF (KEY2.EQ.'POLA') THEN
C     
              IPOLAR = 1
              ITOR = 1
              I = I + 1
C     ************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
              TOR = VALUE(I)
C     
C---- SYNC DIVV
C     
            ELSE IF (KEY2.EQ.'DIVV') THEN
              I = I + 1
C     
C     ************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
              DIVVD = VALUE(I)
              DIVV = 0.5*DTOR*DIVVD
              IDIVV = 1
C     
C---- SYNC DELCOR
C     
            ELSE IF (KEY2.EQ.'DELC') THEN
              I = I + 1
C     
C     ************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
              DELCOR = VALUE(I)
            ELSE
              WRITE (IOUT,FMT=6130) KEY2
              IF (ONLINE) WRITE (ITOUT,FMT=6130) KEY2
            END IF
C     
            GO TO 150
          END IF
        END IF
C     
C---- Set parameters for px9.6 station
C     
        KEYPX = 1
        DIVHD = PX96(1)
        DIVVD = PX96(2)
        DIVH = 0.5*DTOR*DIVHD
        DIVV = 0.5*DTOR*DIVVD
        IDIVH = 1
        IDIVV = 1
        DELAMB = PX96(3)
        DELCOR = PX96(4)
        Hbeamline='DLAB-PX9.6'
        GO TO 170
C     
C---- Set parameters for px7.2 station
C     
 160    KEYPX = 1
        DIVHD = PX72(1)
        DIVVD = PX72(2)
        DIVH = 0.5*DTOR*DIVHD
        DIVV = 0.5*DTOR*DIVVD
        IDIVH = 1
        IDIVV = 1
        DELAMB = PX72(3)
        DELCOR = PX72(4)
        Hbeamline='DLAB-PX7.2'
C     
 170    CONTINUE
C     
C---- SERIAL....serial no. of packs and phi ranges
C     
      ELSE IF (KEY.EQ.'SERI') THEN
C     
C     
C---- Not allowed if inputting keywords via menu
C     
        IF (MODE.EQ.3) THEN
          WRITE(IOUT,FMT=7260)
          WRITE(ITOUT,FMT=7260)
          LINE = ' '
          WRITE(LINE,FMT=7476)
 7476     FORMAT('Cannot use this keyword from menu.')
          CALL MXDWIO(LINE,2)
          GOTO 50
        END IF
        IPROKWD = 1
        WRITE(IOUT,FMT=7171)
        IF (ONLINE) WRITE(ITOUT,FMT=7171)
 7171   FORMAT(/,1X,'***** WARNING *****',
     +       /,1X,'***** WARNING *****',
     +       /,1X,'***** WARNING *****',/,1X,
     +       'The SERIAL keyword is now obsolete, please use the',
     +       ' PROCESS keyword',/,1X,'eg PROCESS 1 to 20 START 5.0',
     +       ' ANGLE 0.5 BLOCK 5 ADD 1000')
        STOP
C     
C
C     
C---- PROCESS (new version of SERIAL)  no. of packs and phi ranges
C     eg PROCESS 1 TO 20 START 5 ANGLE 1.0
C     or PROCESS 1 20 START 5 OSC 1.0
C     
      ELSE IF (KEY.EQ.'PROC') THEN
        IF(MULTIMTZ.AND.MTZOPEN) THEN
          MTZOPEN = .FALSE.
          MTZPRT = 1
C     *********************
          if (debug(52))then
            write(iout,fmt=7740)
            if (online) write(itout,fmt=7740)
 7740       format(1x,'**closing MTZ file in control 1 ***')
          end if
          CALLED = 'LWCLOS'
          COMMENT = 'To close MTZ file '
          IF (DEBUG(99)) 
     +        CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)

          CALL LWCLOS(MTZOUT,MTZPRT)
          IF(LBEST.AND.BESTOPEN)THEN
            CLOSE(BESTHKL)
            BESTOPEN = .FALSE.
          ENDIF
C     *********************
        END IF
        MINBATCH = -999
        NRUN = 0

C     
C     
C---- First check that an IMAGE keyword has not been given
C     
        IF (IMGKWD.GT.0.AND..NOT.SOCKLO) THEN
          WRITE(IOUT,FMT=7310)
          IF (ONLINE) WRITE(ITOUT,FMT=7310)
 7310     FORMAT(/,1X,'***** ERROR *****',/,1X,'IMAGE and PROCESS ',
     +         'keywords must NOT be given in the same "run".',/,1X,
     +         'This keyword will be ignored')
          GOTO 50
        END IF
C     
C     
c     DPSINDEX = .not.dpsdone
        DPSINDEX = .false.
        NOGO = .FALSE.
        IPROKWD = 1
        ISTRT = 0
        IANGLE = 0
        ISERADD = 0
        INERR = .FALSE.
C     .. HRP 25.05.2006
        INSEG = 0
        INFIX = 0
C     
C---- Not allowed if inputting keywords via menu
C     
        IF (MODE.EQ.3) THEN
          WRITE(IOUT,FMT=7260)
          WRITE(ITOUT,FMT=7260)
 7260     FORMAT(1X,'**** ILLEGAL INPUT ***',/,1X,'Use the',
     +         ' "Integrate" menu option')
          LINE = ' '
          WRITE(LINE,FMT=7476)
          CALL MXDWIO(LINE,2)
          GOTO 50
        END IF
c     -harvest
chrp19082003 - turned on/off with keyword 'HARVEST'
C        DOHARVEST = .true.
c     -harvest
        CALL KWPROC(IBEG,IEND,ITYP,VALUE,
     $       IDEC,NTOK,IFIRSTPACK,NPARM,NSER,NSERRUN,IPACKF,IPACKL,
     +       INERR,BADKEY,NEWREAD,NFLEFT,NLLEFT,PHILEFT,PHISLEFT,
     +       ISERLEFT,NPROCRUN,LINE)
CAL     +       ISERLEFT,NRLEFT,NPROCRUN,MTZOPEN,LINE)
        IF (INERR.OR.BADKEY) THEN
          IF (ONLINE) THEN
            IF (COMREAD) THEN
              COMREAD = .FALSE.
              ITIN = ITINS
              CLOSE (UNIT=ICOMM)
            END IF
            BADLINE = .TRUE.
            GOTO 50
          END IF
          NSHUTERR = 10
          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
        END IF
        IF (NEWREAD) GOTO 50

C     
C---- INTEnsity
C     
      ELSE IF (KEY.EQ.'INTE') THEN
C     
C     ************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
        MINT = NINT(VALUE(2))
C     
C---- SEPAration x x
C     
      ELSE IF (KEY.EQ.'SEPA') THEN
C     
        IISEP = .TRUE.
        IF ((NTOK.GE.3).AND.(ITYP(2).EQ.2)) THEN
C     ************************************
          CALL MKEYNM(2,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
C---- Input in mms - stored as an integer in 10 microns units
C     
          IF(VALUE(2).GT.1E-2.AND.VALUE(3).GT.1E-2)THEN
            IXSEP = NINT(VALUE(2)*100.0)
            IYSEP = NINT(VALUE(3)*100.0)
            ISEP = 2
          ELSE
C
C---- if SEPARATION 0 0 we treat it as no keyword given
            IXSEP = 0
            IYSEP = 0
            ISEP = 0
            IISEP = .FALSE.
          ENDIF
          ICOUNT = 3
C     
C---- Check for subkeywords TRIM, OVERLAP, CLOSE, NOTCLOSE, UPDATE, AUTO
C     
        ELSE
          ICOUNT = 1
        END IF
C     
        IF (NTOK.GT.1) THEN
 174      ICOUNT = ICOUNT + 1
          IF (ICOUNT.GT.NTOK) GOTO 175
          SUBKEY = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     *********** 
          CALL CCPUPC(SUBKEY)
C     ***********
C     
          IF (SUBKEY.EQ.'TRIM') THEN
            ICOUNT = ICOUNT + 1
            IF(IISEP)THEN
C     
C     ************************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
C     
              ITRIM = NINT(VALUE(ICOUNT))
            ENDIF
          ELSE IF (SUBKEY.EQ.'OVER') THEN
            ICOUNT = ICOUNT + 1
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
C     
            NOVERLAP = NINT(VALUE(ICOUNT))
C
C----- Set CLOSE options
C     
          ELSE IF (SUBKEY.EQ.'CLOS') THEN
            DENSE = .TRUE.
            PKONLY = .TRUE.
C
C----- Unnet CLOSE options with NOTCLOSE
C     
C
          ELSE IF (SUBKEY.EQ.'NOTC') THEN
            DENSE = .FALSE.
            PKONLY = .FALSE.
C
C---- UPDATE, Start with this value but allow updating
C     
          ELSE IF (SUBKEY.EQ.'UPDA') THEN
            ISEP = 3
C
C---- AUTO, determine separation from scratch
C     
          ELSE IF (SUBKEY.EQ.'AUTO') THEN
            ISEP = 0
C     
C     Not recognised
C     
          ELSE
            IF (TRAPERR) INPERR = .TRUE.
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
            IISEP = .FALSE.
          END IF
          IF (ICOUNT.LT.NTOK) GOTO 174
        END IF
 175    CONTINUE
C     
C---- If reading keyword input for menu, convert to MINDTX,MINDTY
C     
        IF (MODE.EQ.3) THEN
C---- Convert spot separations (IXSEP,IYSEP) into "ideal detector"
C     coordinate
C     frame, as the spot coordinates (generate file coords) are in this
C     frame
C     
          MINDTX = NINT(ABS(IXSEP*COSOM0 + IYSEP*SINOM0))
          MINDTY = NINT(ABS(IYSEP*COSOM0 + IXSEP*SINOM0))
        END IF
C     
C---- RASTer x x x x x x
C     
      ELSE IF (KEY.EQ.'RAST') THEN
C     
        IRAST = 1
        IIRAST = .TRUE.
C     ************************************
        CALL MKEYNM(5,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
        DO 176 I = 1,5
          IRAS(I) = NINT(VALUE(I+1))
 176    CONTINUE
C     
        DO 177 I = 1,5
          IF (I.EQ.3) THEN
            IF (IRAS(I).LE.0) THEN
              IRAS(I) = 1
              WRITE(IOUT,FMT=6620) (IRAS(K),K=1,5)
              IF (ONLINE) WRITE(ITOUT,FMT=6620)  (IRAS(K),K=1,5)
              IF (BRIEF) WRITE(IBRIEF,FMT=6620)  (IRAS(K),K=1,5)
            END IF
          END IF
 6620     FORMAT(1X,'Corner parameter (NC) must be greater than zero'
     +         ,', new raster parameters',5I4)


          IF (IRAS(I).LT.0) THEN
            IRAS(I) = 1
            WRITE(IOUT,FMT=6621) (IRAS(K),K=1,5)
            IF (ONLINE) WRITE(ITOUT,FMT=6621)  (IRAS(K),K=1,5)
            IF (BRIEF) WRITE(IBRIEF,FMT=6621)  (IRAS(K),K=1,5)
          END IF
 6621     FORMAT(1X,'*** Negative raster parameters not permitted',
     +         ', new raster parameters',5I4)
 177    CONTINUE
C     
C---- Check this gives more than one peak pixel
C     
        CALL SETMASK(MASK,IRAS)
        CALL SETSUMS(MASK,IRAS,PKSUMS)
        IF (PKSUMS(5).LT.2.0) THEN
          WRITE(IOUT,FMT=7620) (IRAS(K),K=1,5),PKSUMS(5)
          IF (ONLINE) WRITE(ITOUT,FMT=7620) (IRAS(K),K=1,5),PKSUMS(5)
 7620     FORMAT(/,/,1X,'***** ERROR *****',/,1X,'The current ',
     +         'raster parameters (',5I4,') give only',F3.0,
     +         ' peak pixels.',/,1X,'There must be at least 2.')
          IF (ONLINE) THEN
            IF (COMREAD) THEN
              COMREAD = .FALSE.
              ITIN = ITINS
              CLOSE (UNIT=ICOMM)
            END IF
            BADLINE = .TRUE.
            GOTO 50
          END IF
          NSHUTERR = 15
          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
        END IF
        IF (NTOK.EQ.7) THEN
C     ************************************
          CALL MKEYNM(1,7,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          IF ((VALUE(7).LT.1.0).OR.(VALUE(7).GT.50)) THEN
            WRITE(IOUT,FMT=6460)
            IF (ONLINE) WRITE(ITOUT,FMT=6460)
 6460       FORMAT(1X,'**** WARNING ****',/,1X,'Pixel size must now',
     +           ' be given in mm using a PIXEL keyword',/,1X,
     +           'The value given here has been ignored')
            IF (ONLINE) THEN
              IF (COMREAD) THEN
                COMREAD = .FALSE.
                ITIN = ITINS
                CLOSE (UNIT=ICOMM)
              END IF
              BADLINE = .TRUE.
              GOTO 50
            END IF
            NSHUTERR = 16
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF
        END IF
C     
C---- PIXEL...pixel size (slow direction in image) in mm, optionally
C     followed by pixel size in fast direction
C     
      ELSE IF (KEY.EQ.'PIXE') THEN
        IPIX = 1
        IIPIX = .TRUE.
C     
C     ************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
        RAST = VALUE(2)
        IF (RAST.GT.1.0) THEN
          WRITE(IOUT,FMT=6461) RAST
          IF (ONLINE) THEN
            WRITE(ITOUT,FMT=6461) RAST
            IF (COMREAD) THEN
              COMREAD = .FALSE.
              ITIN = ITINS
              CLOSE (UNIT=ICOMM)
            END IF
            BADLINE = .TRUE.
            GOTO 50
          END IF
          NSHUTERR = 17
          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
        END IF
 6461   FORMAT(1X,'**** WARNING ****',/,1X,'Pixel size must now',
     +       ' be given in mm, NOT microns. Input value of',
     +       F8.2,'mm is unreasonable')       
C     
C---- Check for pixel size in fast direction, if present, use it to
C     define YSCAL
C     
        IF (NTOK.EQ.3) THEN     
          IPIXY = 1
C     ************************************
          CALL MKEYNM(1,3,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
          RASTY = VALUE(3)
          YSCAL = RAST/RASTY
          YSCALIN = YSCAL
C     
C---- Check if YSCAL already defined on DISTORTION keyword
C     
          IF (IYSCAL.NE.0) THEN
            WRITE(IOUT,FMT=6462) YSCAL
            IF (ONLINE) WRITE(ITOUT,FMT=6462) YSCAL
            IF (BRIEF) WRITE(IBRIEF,FMT=6462) YSCAL
 6462       FORMAT(1X,'***** WARNING *****',/,1X,'YSCAL calculated',
     +           ' from the ratio of the pixel sizes in the slow a',
     +           'nd ',/,1X,'fast directions (',F6.4,') will super',
     +           'cede the value given on by DISTORTION YSCAL keyw',
     $           'ords')
          END IF
        END IF
C     
C---- DISPersion x
C     
      ELSE IF (KEY.EQ.'DISP') THEN
C     ************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
        IIDISP = .TRUE.
        DELAMB = VALUE(2)
        IDELAMB = 1
C     
C---- LIMIts
C     
      ELSE IF (KEY.EQ.'LIMI') THEN
        I = 1
 178    CONTINUE
        I = I + 1
        IF (I.LE.NTOK) THEN
C     
C---- Skip if no more tokens on line
C     
          KEY2 = LINE(IBEG(I) :IEND(I))
C     
C     **********
          CALL CCPUPC(KEY2)
C     **********
C     
C---- LIMIts RSCAN... radius of scanned circle cented on the middle
C     of the image
C     
C     
          IF (KEY2.EQ.'RSCA') THEN
            I = I + 1
C     
C     ************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
            RSCAN = VALUE(I)*100.0
C     
C---- If in MODE=3 (ie keyword input from MXDSPL) set RSCANSQ
C     
            IF (MODE.EQ.3) RSCANSQ = RSCAN*RSCAN
C     
C---- LIMIts XSCAN... maximum X coordinate in digitised image
C     
C     
          ELSE IF (KEY2.EQ.'XSCA') THEN
            I = I + 1
C     
C     ************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
            XSCAN = VALUE(I)*100.0
C     
C---- LIMIts YSCAN... maximum Y coordinate in digitised image
C     
C     
          ELSE IF (KEY2.EQ.'YSCA') THEN
            I = I + 1
C     
C     ************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
            YSCAN = VALUE(I)*100.0
C     
C---- LIMIts CENTRE coordinates for the centre of the circle of radius
C     RSCAN defining the useable part of the image
C     
C     
          ELSE IF (KEY2.EQ.'CENT') THEN
            I = I + 1
C     
C     ************************************
            CALL MKEYNM(2,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
            RSCANX = VALUE(I)*100.0
            I = I + 1
            RSCANY = VALUE(I)*100.0
C     
C     
C---- LIMIts RMAX
C     
          ELSE IF (KEY2.EQ.'RMAX') THEN
            I = I + 1
C     
C     ************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
            RMAX = VALUE(I)*100.0
C     
C---- LIMIts RMIN
C     
          ELSE IF (KEY2.EQ.'RMIN') THEN
            I = I + 1
C     
C     ************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
            RMIN = VALUE(I)*100.0
C     
C---- LIMIts RCENTRE coordinates for the centre of the circle of radius
C     RMIN defining the useable part of the image
C     
C     
          ELSE IF (KEY2.EQ.'RCEN') THEN
            I = I + 1
C     
C     ************************************
            CALL MKEYNM(2,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
            RMINXINP = VALUE(I)*100.0
            RMINX = RMINXINP
            I = I + 1
            RMINY = VALUE(I)*100.0
            IRMINXY = 1
C     
C---- LIMIts XMIN
C     
          ELSE IF (KEY2.EQ.'XMIN') THEN
            I = I + 1
C     
C     ************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
            XMIN = VALUE(I)*100.0
C     
C---- LIMIts XMAX
C     
          ELSE IF (KEY2.EQ.'XMAX') THEN
            I = I + 1
C     
C     ************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
            XMAX = VALUE(I)*100.0
C     
C---- LIMIts YMIN
C     
          ELSE IF (KEY2.EQ.'YMIN') THEN
            I = I + 1
C     
C     ************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
            YMIN = VALUE(I)*100.0
C     
C---- LIMITS YMAX
C     
          ELSE IF (KEY2.EQ.'YMAX') THEN
            I = I + 1
C     
C     ************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
            YMAX = VALUE(I)*100.0
C     
C---- LIMIts EXCLude
          ELSE IF (KEY2.EQ.'EXCL'.OR.KEY2.EQ.'RECT') THEN
            I = I + 1
            NXYEXC = NXYEXC + 1
            IF (NXYEXC.GT.MAXEXCL) THEN
              WRITE(IOUT,FMT=7151) MAXEXCL,NXYEXCDET
              IF (ONLINE) WRITE(ITOUT,FMT=7151) MAXEXCL,NXYEXCDET
 7151         FORMAT(1X,'*** ERROR ***',/,1X,'A maximum of',I3,
     +             ' rectangular regions can be excluded.',/,1X,
     +             'of which',I3,' are used by detector specific',
     +             ' areas.',/,1X,
     +             'The remainder will be ignored')
              GOTO 178
            END IF
C     ************************************
            CALL MKEYNM(4,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
            XYEXX1 = VALUE(I)
            XYEXY1 = VALUE(I+1)
            XYEXX2 = VALUE(I+2)
            XYEXY2 = VALUE(I+3)
            I = I + 3
            IF (XYEXX2.LT.XYEXX1) THEN
              TEMPX = XYEXX1
              XYEXX1 = XYEXX2
              XYEXX2 = TEMPX
            END IF
            IF (XYEXY2.LT.XYEXY1) THEN
              TEMPX = XYEXY1
              XYEXY1 = XYEXY2
              XYEXY2 = TEMPX
            END IF
            XYEXC(1,NXYEXC) = 100.0*XYEXX1
            XYEXC(2,NXYEXC) = 100.0*XYEXY1
            XYEXC(3,NXYEXC) = 100.0*XYEXX2
            XYEXC(4,NXYEXC) = 100.0*XYEXY2
            WRITE (IOUT,FMT=7155) XYEXX1,XYEXY1,XYEXX2,XYEXY2
            IF (ONLINE) WRITE (ITOUT,FMT=7155) XYEXX1,XYEXY1,
     +           XYEXX2,XYEXY2
 7155       FORMAT (1X,'Exclude region with corners',F7.2,',',F7.2,
     +           ' and',F7.2,',',F7.2)
C     
C---- LIMIts QUADrilateral (EXCLUDE implied)
          ELSE IF (KEY2.EQ.'QUAD') THEN
            I = I + 1
            IQUADNO = IQUADNO + 1
            IF (IQUADNO.GT.MAXEXCL) THEN
              WRITE(IOUT,FMT=7156) MAXEXCL
              IF (ONLINE) WRITE(ITOUT,FMT=7156) MAXEXCL
 7156         FORMAT(1X,'*** ERROR ***',/,1X,'A maximum of',I3,
     +             ' quadrilateral regions can be excluded.',/,1X,
     +             'The remainder will be ignored')
              GOTO 178
            END IF
C     ************************************
            CALL MKEYNM(8,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
            CORNER(1,1) = 100 * VALUE(I)
            CORNER(1,2) = 100 * VALUE(I+1)
            CORNER(2,1) = 100 * VALUE(I+2)
            CORNER(2,2) = 100 * VALUE(I+3)
            CORNER(3,1) = 100 * VALUE(I+4)
            CORNER(3,2) = 100 * VALUE(I+5)
            CORNER(4,1) = 100 * VALUE(I+6)
            CORNER(4,2) = 100 * VALUE(I+7)
            CNREXC(1,IQUADNO) = 100 * VALUE(I)
            CNREXC(2,IQUADNO) = 100 * VALUE(I+1)
            CNREXC(3,IQUADNO) = 100 * VALUE(I+2)
            CNREXC(4,IQUADNO) = 100 * VALUE(I+3)
            CNREXC(5,IQUADNO) = 100 * VALUE(I+4)
            CNREXC(6,IQUADNO) = 100 * VALUE(I+5)
            CNREXC(7,IQUADNO) = 100 * VALUE(I+6)
            CNREXC(8,IQUADNO) = 100 * VALUE(I+7)
            CALL QUADTL(1,IQUADNO,CORNER,0.0,0.0,IERR)
            WRITE(IOUT,FMT=7157)(VALUE(J),J=I,I+7)
            IF(ONLINE)WRITE(ITOUT,FMT=7157)(VALUE(J),J=I,I+7)
 7157       FORMAT(1X,'Exclude quadrilateral with corners',
     +           3('(',F7.2,',',F7.2,')'),' and (',F7.2,',',F7.2,')')
            I = I + 7
C     
C---- LIMIts ARC  (EXCLUDE implied) - requires x,y (centre), radius(?) 
C     (all mm), phihigh, philow in degrees, plus inside or outside 
C     - 6 arguments
          ELSE IF (KEY2(1:3).EQ.'ARC') THEN
            I = I + 1
            IARCNO = IARCNO + 1
            IF (IARCNO.GT.MAXEXCL) THEN
              WRITE(IOUT,FMT=7156) MAXEXCL
              IF (ONLINE) WRITE(ITOUT,FMT=7158) MAXEXCL
 7158         FORMAT(1X,'*** ERROR ***',/,1X,'A maximum of',I3,
     +             ' arcs can be excluded.',/,1X,
     +             'The remainder will be ignored')
              GOTO 178
            END IF
C     ************************************
            CALL MKEYNM(5,I,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     ..
C     .. need to convert to PIXELS only if we know RAST & YSCAL - late in GUI 
C        code
                 EXCARC(1,IARCNO) = VALUE(I)
                 EXCARC(2,IARCNO) = VALUE(I+1)
                 EXCARC(3,IARCNO) = VALUE(I+2)
                 EXCARC(4,IARCNO) = VALUE(I+3) * DTOR
                 EXCARC(5,IARCNO) = VALUE(I+4) * DTOR

C     ..
C     .. six points in pixels
                 IEXCARC(1,IARCNO) = -999
                 I = I + 5
                 KEY2 = LINE(IBEG(I):IEND(I))
                 CALL CCPUPC(KEY2)
C     assume INSIDE default for moment 
                 LEXCARC(IARCNO) = .FALSE.
                 IF(KEY2.EQ.'INSI')LEXCARC(IARCNO) = .TRUE.
                 IF(KEY2.EQ.'OUTS')LEXCARC(IARCNO) = .FALSE.
                 I = I + 1
C
C---- LIMIts REMOve - removes selected conditions set by EXCLude (RECTangle), 
C     QUADrilateral and ARC above - this is necessary for the new GUI and
C     may also be necessary for automation. Number refers to slot.
          ELSE IF (KEY2.EQ.'REMO') THEN
            REMOVEALL = .FALSE.
            I = I + 1
            KEY2 = LINE(IBEG(I) :IEND(I))
            I = I + 1
C     
C           *****************
            CALL CCPUPC(KEY2)
C           *****************
C     .. ALL ..
            IF(KEY2(1:3).EQ.'ALL')REMOVEALL = .TRUE.
C
C---- If LIMITS REMOVE EVERYTHING given, then even remove detector
C     specific exclusions. There is no way to reverse this.
C
            IF(KEY2.EQ.'EVER') THEN
              NXYEXCDET = 0
              NOTILES = .TRUE.
            END IF
C
C     .. RECTangle .. this code is very similar to that in 'mxmask' so any changes here 
C     ..              should also be made there.
            IF(KEY2.EQ.'RECT'.OR.REMOVEALL)THEN
              DO 9015 II = NXYEXCDET+1,MAXEXCL
                EMPTY(II) = .FALSE.
 9015         ENDDO
              IF(REMOVEALL)THEN
                ISLOTB = NXYEXCDET + 1
                ISLOTE = NXYEXC
              ELSE
                ISLOTB = I
                ISLOTE = NTOK
              ENDIF
C     .. go through list
              DO 9020 II=ISLOTB,ISLOTE,1
C     .. zero chosen slots
                IF(REMOVEALL)THEN
                  J = II
                ELSE
                  J = NINT(VALUE(II))
                ENDIF
                XYEXC(1,J) = 0.0
                XYEXC(2,J) = 0.0
                XYEXC(3,J) = 0.0
                XYEXC(4,J) = 0.0
                EMPTY(J) = .TRUE.
 9020         ENDDO
              J = 1
              K = NXYEXCDET + 1
              DO 9040 II=NXYEXCDET+1,NXYEXC,1
C     .. zero remaining used slots
                IF(EMPTY(II))THEN
                  EMPTY(II) = .FALSE.
                ELSE
                  XYEXC(1,K) = XYEXC(1,II)
                  XYEXC(2,K) = XYEXC(2,II)
                  XYEXC(3,K) = XYEXC(3,II)
                  XYEXC(4,K) = XYEXC(4,II)
                  K = K + 1
                ENDIF
 9040         ENDDO
              DO 9050 II=K,NXYEXC,1
C     .. zero slots at end
                XYEXC(1,II) = 0.0
                XYEXC(2,II) = 0.0
                XYEXC(3,II) = 0.0
                XYEXC(4,II) = 0.0
                EMPTY(II) = .TRUE.
 9050         ENDDO
              NXYEXC = NXYEXC - NTOK + 3
              IF(REMOVEALL)NXYEXC = NXYEXCDET
              IF(NXYEXC.GT.0)THEN
                WRITE(IOUT,FMT=16070)NXYEXC
                IF(ONLINE)WRITE(ITOUT,FMT=16070)NXYEXC
16040           FORMAT('Slot ',I2,': corners (',F5.1,',',F5.1,'), (',
     $               F5.1,',',F5.1,')')
16070           FORMAT('Remaining ',I2,' slots contain: ')
                DO 9060 II = 1,NXYEXC,1
                  WRITE(IOUT,FMT=16040)II,(XYEXC(J,II)*0.01,J=1,4)
                  IF(ONLINE)WRITE(ITOUT,FMT=16040)II,
     $                 (XYEXC(J,II)*0.01,J=1,4)
CAL                  WRITE(IOUT,FMT=16040)II,(XYEXCSP(J,II)*0.01,J=1,4)
CAL                  IF(ONLINE)WRITE(ITOUT,FMT=16040)II,
CAL     $                 (XYEXCSP(J,II)*0.01,J=1,4)
 9060           ENDDO
              ENDIF
            ENDIF
C     .. QUADrilateral ..
            IF(KEY2.EQ.'QUAD'.OR.REMOVEALL)THEN
              DO 9105 II = 1,10
                EMPTY(II) = .FALSE.
 9105         ENDDO
              IF(REMOVEALL)THEN
                ISLOTB = 1
                ISLOTE = IQUADNO
              ELSE
                ISLOTB = I
                ISLOTE = NTOK
              ENDIF
C     .. go through list
              DO 9120 II=ISLOTB,ISLOTE,1
C     .. zero chosen slots
                IF(REMOVEALL)THEN
                  J = II
                ELSE
                  J = NINT(VALUE(II))
                ENDIF
                DO 9110 K = 1,8
                  CNREXC(K,J) = 0.0
 9110           ENDDO
                EMPTY(J) = .TRUE.
 9120         ENDDO
              K = 1
              DO 9140 II=1,IQUADNO,1
C     .. zero remaining used slots
                IF(EMPTY(II))THEN
                  EMPTY(II) = .FALSE.
                ELSE
                  DO 9130 J=1,8
                    CNREXC(J,K) = CNREXC(J,II)
 9130             ENDDO
                  K = K + 1
                ENDIF
 9140         ENDDO
              DO 9150 II=K,NXYEXC,1
C     .. zero slots at end
                DO 9145 J = 1,8
                  CNREXC(J,II) = 0.0
 9145           ENDDO
                EMPTY(II) = .TRUE.
 9150         ENDDO
              IQUADNO = IQUADNO - NTOK + 3
              IF(REMOVEALL)IQUADNO = 0
              IF(IQUADNO.GT.0)THEN
                WRITE(IOUT,FMT=16070)IQUADNO
                IF(ONLINE)WRITE(ITOUT,FMT=16070)IQUADNO
                DO 9160 II = 1,IQUADNO,1
                  WRITE(IOUT,FMT=16140)II,(CNREXC(J,II)*0.01,J=1,8)
                  IF(ONLINE)WRITE(ITOUT,FMT=16140)II,
     $                 (CNREXC(J,II)*0.01,J=1,8)

16140             FORMAT('Slot ',I2,': corners ',
     $                 3('(',F5.1,',',F5.1,'), '),
     $                 '(',F5.1,',',F5.1,')')
 9160           ENDDO
              ENDIF
            ENDIF
C     .. ARC ..
            IF(KEY2(1:3).EQ.'ARC'.OR.REMOVEALL)THEN
              DO 9205 II = 1,10
                EMPTY(II) = .FALSE.
 9205         ENDDO
              IF(REMOVEALL)THEN
                ISLOTB = 1
                ISLOTE = IARCNO
              ELSE
                ISLOTB = I
                ISLOTE = NTOK
              ENDIF
C     .. go through list
              DO 9220 II=1,NTOK,1
C     .. zero chosen slots
                IF(REMOVEALL)THEN
                  J = II
                ELSE
                  J = NINT(VALUE(II))
                ENDIF
                DO 9210 K = 1,5
                  EXCARC(K,J) = 0.0
                  IEXCARC(K,J) = 0
 9210           ENDDO
                IEXCARC(6,J) = 0
                LEXCARC(J) = .FALSE.
                EMPTY(J) = .TRUE.
 9220         ENDDO
              K = 1
              DO 9240 II=1,IARCNO,1
C     .. zero remaining used slots
                IF(EMPTY(II))THEN
                  EMPTY(II) = .FALSE.
                ELSE
                  DO 9230 J = 1,5
                    EXCARC(J,K) = EXCARC(J,II)
                    IEXCARC(J,K) = IEXCARC(J,II)
 9230             ENDDO
                  IEXCARC(6,K) = IEXCARC(6,II)
                  LEXCARC(K) = LEXCARC(II)
                  K = K + 1
                ENDIF
 9240         ENDDO
              DO 9250 II=K,IARCNO,1
C     .. zero slots at end
                DO 9245 J = 1,5
                  EXCARC(J,II) = 0.0
                  IEXCARC(J,II) = 0
 9245           ENDDO
                IEXCARC(6,II) = 0
                LEXCARC(II) = .FALSE.
                EMPTY(II) = .TRUE.
 9250         ENDDO
              IARCNO = IARCNO - NTOK + 3
              IF(REMOVEALL)IARCNO = 0
              IF(IARCNO.GT.0)THEN
                WRITE(IOUT,FMT=16070)IARCNO
                IF(ONLINE)WRITE(ITOUT,FMT=16070)IARCNO
C     
C---- loop round slots
C
                DO 9260 II = 1,IARCNO,1
                  IF(LEXCARC(I))THEN
                    INOUT = 'Outside'
                  ELSE
                    INOUT = 'Inside'
                  ENDIF
                  WRITE(IOUT,FMT=16240)II,(EXCARC(J,II),J=1,3),
     $                 (EXCARC(J,II)*RTOD,J=4,5),INOUT(1:LENSTR(INOUT))
                  IF(ONLINE)
     $                 WRITE(ITOUT,FMT=16240)II,(EXCARC(J,II),J=1,3),
     $                 (EXCARC(J,II)*RTOD,J=4,5),INOUT(1:LENSTR(INOUT))
16240             FORMAT('Slot ',I2,': centre (',F6.1,',',F6.1,'),',
     $                 ' radius ',F6.1,' angular limits ',2F5.1,1X,A)
 9260           ENDDO
              ENDIF
            ENDIF
          ELSE
C     
C---- Not recognised
C     
            IF (TRAPERR) INPERR = .TRUE.
            WRITE (IOUT,FMT=6130) KEY2
            IF (ONLINE) WRITE (ITOUT,FMT=6130) KEY2
          END IF
C     
          GO TO 178
        END IF

C     
C---- DIVErgence x [x]
C     
      ELSE IF (KEY.EQ.'DIVE') THEN
        IIDIV = .TRUE.
        N = 1
        IF (NTOK.EQ.3) N = 2
C     
C     ************************************
        CALL MKEYNM(N,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
        DIVHD = VALUE(2)
        IDIVH = 1
        IDIVV = 1
        IF (N.EQ.2) THEN
C     
C---- If only one number set vertical div
C     equal to the hor div
C     
          DIVVD = VALUE(3)
        ELSE
          DIVVD = DIVHD
        END IF
        DIVH = 0.5*DTOR*DIVHD
        DIVV = 0.5*DTOR*DIVVD
        WARN(26) = .FALSE.
c     
c---- dummy keyword for allowing a write for a jpeg file of the image 
C     REMOVE before distribution
c     
      ELSE IF (KEY.EQ.'JPEG')THEN
        JPGOUT = .TRUE.
      ELSE IF (KEY.EQ.'NOJP')THEN
        JPGOUT = .FALSE.
c     
c---- keyword for new, generalized REEK calculations
C     
c     
      ELSE IF (KEY.EQ.'NURE')THEN
        NUREEK = .TRUE.
      ELSE IF (KEY.EQ.'OLDR')THEN
        NUREEK = .FALSE.
c     
c---- keyword for new-style definition of TILT and TWIST
C     REMOVE before distribution
c     
      ELSE IF (KEY.EQ.'NEWT')THEN
        NUTWIST = .TRUE.
        CALL SETDIS(ITILT,ITWIST,1)
      ELSE IF (KEY.EQ.'OLDT')THEN
        NUTWIST = .FALSE.
        CALL SETDIS(ITILT,ITWIST,1)
C     
C---- MOSAIC
C     
      ELSE IF (KEY.EQ.'MOSA') THEN
C     
C---- Check if next keyword is a number or token
C     
        ICOUNT = 2
        DO 1000 WHILE (ICOUNT.LE.NTOK)
          IF (ITYP(ICOUNT).EQ.1) THEN
            SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
            CALL CCPUPC(SUBKEY)
            ICOUNT = ICOUNT + 1
            IF (SUBKEY.EQ.'ESTI')THEN
              NIMAGES = 1
              POWDER = .TRUE.
              ICOUNT = 3
              MOSEST = .TRUE.
              BLKEST = .FALSE.
C     
C---- stuff cribbed from AUTOMATCH (vide infra)
C     
              OTHERS = .TRUE.
              MATCH = .TRUE.
C     
C---- to set this up as an integration rather than post-refinement run
C     
              MULTISEG = .FALSE.
              NSEG = 1
C     
C---- if subkeyword was 'ESTI' - remaining numeric token (if it exists)
C     will be an image number. IMAGE, PHI, TO all assumed here.
C     
c              IF (ITYP(ICOUNT).EQ.2) THEN
c                IF(MOSIMAG.NE.0)NAUTO = 0
c                MOSIMAG = NINT(VALUE(ICOUNT))
c                ICOUNT = ICOUNT + 1
c                IF(.NOT.AUTOINDX.OR.NAUTO.LT.1)THEN
c                  NAUTO = 1
c                  ID = MOSIMAG
c                  NOIMG(1) = MOSIMAG
c                  IDPACK(1) = MOSIMAG
c                  NRUN = 1
c                  NODISPLAY = .TRUE.
C     
c                ENDIF
C     
C---- check for 'IMAGE','PHI' and 'TO' keywords and 
C     
c              ELSE
                DO 1020 WHILE (ICOUNT.LE.NTOK)
C     ..
C     .. no image keyword
                  IF (ITYP(ICOUNT).EQ.2) THEN
                    IF(MOSIMAG.NE.0)NAUTO = 0
                    MOSIMAG = NINT(VALUE(ICOUNT))
                    ICOUNT = ICOUNT + 1
                    IF(.NOT.AUTOINDX.OR.NAUTO.LT.1)THEN
                      NAUTO = 1
                      NOIMG(1) = MOSIMAG
                      IDPACK(1) = MOSIMAG
                      NRUN = 1
                      NODISPLAY = .TRUE.
                    ENDIF
                  ELSEIF(ITYP(ICOUNT).EQ.1)THEN
                    SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
                    ICOUNT = ICOUNT + 1
                    CALL CCPUPC(SUBKEY)
                    IF(SUBKEY.EQ.'IMAG')THEN
                      IF(ITYP(ICOUNT).EQ.2)THEN
                        IF(MOSIMAG.NE.0)NAUTO = 0
                        MOSIMAG = NINT(VALUE(ICOUNT))
                        ICOUNT = ICOUNT + 1
                        IF(.NOT.AUTOINDX.OR.NAUTO.LT.1)THEN
                          NAUTO = 1
                          NOIMG(1) = MOSIMAG
                          IDPACK(1) = MOSIMAG
                          NRUN = 1
                          NODISPLAY = .TRUE.
                        ENDIF
                      ENDIF
                    ELSEIF(SUBKEY(1:3).EQ.'PHI')THEN
                      IF(ITYP(ICOUNT).EQ.2)THEN
                        PHIBEG = VALUE(ICOUNT)
                        PHISTIM(1) = PHIBEG
                        ICOUNT = ICOUNT + 1
                        IANGLE = 1
                        ISTRT = 1
                      ENDIF
                    ELSEIF(SUBKEY(1:2).EQ.'TO')THEN
                      IF(ITYP(ICOUNT).EQ.2)THEN
                        PHIEND = VALUE(ICOUNT)
                        ICOUNT = ICOUNT + 1
                        IANGLE = 1
                      ENDIF
                    ENDIF
                  ENDIF
 1020           ENDDO 
                PHI(1) = 0.5 * (PHIBEG + PHIEND)
c              ENDIF
C     
C---- 
C     
            ENDIF
          ELSEIF (SUBKEY.EQ.'BLOC')THEN
C     
C---- blocksize given in microns, converted to Angstroms for use 
C     internally
C     
            BLKSIZE = 10000.0*VALUE(ICOUNT)
            ICOUNT = ICOUNT + 1
            IF (BLKSIZE.LE.0.0) BLKSIZE = 1000000.0
          ELSE
C     
C---- second token is a number (the mosaicity)
C     
            ETAD = VALUE(ICOUNT)
            ICOUNT = ICOUNT + 1
            ETA = 0.5*DTOR*ETAD
C     
C---- but if mosaicity estimation has already been done, use that value
C     instead
C     
            IF(.NOT.READINLINE)LOGETA = .FALSE.
            IF(LOGETA)ETA = SETA
            IMOSAIC = 1
            WARN(26) = .FALSE.
          END IF
 1000   ENDDO
C     
C---- BLOCKSIZE Mosaic block size in microns
C     
      ELSE IF (KEY.EQ.'BLOC') THEN
C     
C     ************************************
          CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          BLKSIZE = 10000.0*VALUE(2)
          IF (BLKSIZE.LE.0.0) BLKSIZE = 1000000.0
C
C---- CONTROL sets flags for automatically controlling mosaic spread
C             and spot separation (to avoid overlaps).
C
      ELSE IF (KEY.EQ.'CONT') THEN
        CONTROLM = .FALSE.
        CONTROLS = .FALSE.
        IF (NTOK.GT.1) THEN
          ICOUNT = 1
 180      ICOUNT = ICOUNT + 1
          SUBKEY = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C          *********** 
          CALL CCPUPC(SUBKEY)
C          ***********
C     
          IF (SUBKEY.EQ.'MOSA') THEN
            CONTROLM = .TRUE.
          ELSE IF (SUBKEY.EQ.'MAXO') THEN
            CONTROLM = .TRUE.
            ICOUNT = ICOUNT + 1
C     
C     
C---- MAXOVERLAP Maximum allowed overlap as a fraction of all
C                reflections.
C     If the predicted fraction of overlapped spots is
C     greater than this, the minimum spot separation or
C     mosaic spread will be adjusted to reach this value.
C
C     ************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
            OVERLAPMAX = VALUE(ICOUNT)
            IF (OVERLAPMAX.GT.1.0) OVERLAPMAX = 0.01*OVERLAPMAX
          ELSE IF (SUBKEY.EQ.'SEPA') THEN
            CONTROLS = .TRUE.
          END IF
          IF (ICOUNT.LT.NTOK) GOTO 180
        END IF
C     
C     
C---- CELL
C     
      ELSE IF (KEY.EQ.'CELL') THEN
        EXPLICIT_CELL = .TRUE.
        CELLKEEP = .FALSE.
        ICOUNT = 2
        ANGLES = .TRUE.
C     
C---- Check if next keyword is a number or token
C     
        IF (ITYP(2).EQ.1) THEN
          SUBKEY = LINE(IBEG(2):IEND(2))
          CALL CCPUPC(SUBKEY)
          ICOUNT = 3
          IF (SUBKEY.EQ.'KEEP') THEN
            CELLKEEP = .TRUE.
          ELSE
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
            IF (ONLINE) THEN
              IF (COMREAD) THEN
                COMREAD = .FALSE.
                ITIN = ITINS
                CLOSE (UNIT=ICOMM)
              END IF
              BADLINE = .TRUE.
              GOTO 50
            END IF
            NSHUTERR = 18
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF        
        END IF
        IF (NTOK.LT.7) THEN
          ANGLES = .FALSE.
C     ************************************
          CALL MKEYNM(3,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          CELL(4) = 90.0
          CELL(5) = 90.0
          CELL(6) = 90.0
          IF ((NUMSPG.GE.143).AND.(NUMSPG.LE.194)) CELL(6) = 120.0
          K = 3
        ELSE
C     ************************************
          CALL MKEYNM(6,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          K = 6
        END IF
C     
C---- Test for error in read
C     
        IF (IOERR) THEN
          IF (ONLINE) THEN
            IF (COMREAD) THEN
              COMREAD = .FALSE.
              ITIN = ITINS
              CLOSE (UNIT=ICOMM)
            END IF
            GOTO 50
          END IF
          STOP
        END IF
C     
C     Check if cell parameters have already been read from UMAT file
        IF (((IUMAT.EQ.1).OR.(IMAT.EQ.1)).AND.(.NOT.RPTFIRST)) THEN
          WRITE(IOUT,FMT=7017) 
          IF (ONLINE) WRITE(ITOUT,FMT=7017) 
 7017     FORMAT (' ** BEWARE ** Cell parameters given on CELL card',
     +         ' will overwrite those read from',/,1X,'file given',
     +         ' as UMAT or MATRIX')
        END IF
C     
        IF (VALUE(ICOUNT).GT.1.0) THEN
          ICELL = 1
          DO 184 I = 1,K
            CELL(I) = VALUE(I+ICOUNT-1)
 184      CONTINUE
        ELSE
          ICELL = -1
          DO 186 I = 1,K
            RCELL(I) = VALUE(I+ICOUNT-1)
 186      CONTINUE
        END IF
C     
C---- if we are using DPS indexing in background, we put the contents of
C     
C     CELL into KCELL; overrides KEEP keyword so comes before it.
C     
        IF(DPSINDEX)THEN
          DO 1861 I=1,6,1
            KCELL(I) = CELL(I)
 1861     ENDDO
        END IF
C     
C---- Test for KEEP keyword after cell
C     
        IF (ITYP(NTOK).EQ.1) THEN
          SUBKEY = LINE(IBEG(NTOK):IEND(NTOK))
          CALL CCPUPC(SUBKEY)
          IF (SUBKEY.EQ.'KEEP') THEN
            CELLKEEP = .TRUE.
          ELSE
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
            IF (ONLINE) THEN
              IF (COMREAD) THEN
                COMREAD = .FALSE.
                ITIN = ITINS
                CLOSE (UNIT=ICOMM)
              END IF
              BADLINE = .TRUE.
              GOTO 50
            END IF
            NSHUTERR = 19
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF        
        END IF        

C     
C---- STRATEGY data collection option
C     
      ELSE IF (KEY.EQ.'STRA') THEN
C     
C     
c     -harvest
        TESTRAT = .FALSE.
        DOHARVEST = .false.
c     -harvest
C
C--- Increas DEFVOL to avoid rounding errors giving higher anom completeness
C    than overall completeness. AGWL 17/7/09
C
C        DEFVOL = 400000
        DEFVOL = 4000000
        INERR = .FALSE.
        STRATEGY = .TRUE.
        NEWSTRAT = .FALSE.
        NSEGM = NSEGM + 1
        NSEGRD = NSEGRD + 1
        INSPEED = 0
        IF (NSEGM.GT.NSEGMAX) THEN
          WRITE(IOUT,FMT=6471) NSEGMAX
          IF (ONLINE) WRITE(ITOUT,FMT=6471) NSEGMAX
 6471     FORMAT(//,1X,'***** FATAL ERROR *****',/,1X,
     +         'Only',I5,' segments allowed in STRATEGY option')
          STOP
        END IF
        IFIRSTONE(NSEGM) = 0
        ICOUNT = 1
        IIPHI = 0
        IISIZE = 0
        SIZESET = .FALSE.
        IF (NTOK.EQ.1) GOTO 192
C     
 187    ICOUNT = ICOUNT + 1
        SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     **************
        CALL CCPUPC(SUBKEY)
C     **************
C---- START
C     
        IF (SUBKEY.EQ.'STAR') THEN
C     
          IISTART = 1
          IF (IAUTO.EQ.1) THEN
            WRITE(IOUT,FMT=6486)
            IF (ONLINE) WRITE(ITOUT,FMT=6486)
            INERR = .TRUE.
          END IF
 6486     FORMAT(1X,'**** ERROR *****',/,1X,'Cannot have STRATEGY',
     +         ' AUTO keyword and a STRATEGY START keyword ',
     +         ' in the same run.')
C     
          ICOUNT = ICOUNT + 1
          IIPHI = 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          IF (IOERR) INERR = .TRUE.
          PHIST(NSEGM) = VALUE(ICOUNT) + ISTRUN*360.0
C     
C---- Values must be integral
C     
          TEMPX = REAL(NINT(PHIST(NSEGM)))
          IF (ABS(PHIST(NSEGM)-TEMPX).GT.0.01) THEN
            WRITE(IOUT,FMT=6473) PHIST(NSEGM),TEMPX
            IF (ONLINE) WRITE(ITOUT,FMT=6473) PHIST(NSEGM),TEMPX
            PHIST(NSEGM) = TEMPX
          END IF
 6473     FORMAT(1X,'*** WARNING *** Phi values must be integral',
     +         ' input value',F7.1,' reset to',F7.1)
C     **************
C     
C---- END phi
C     
        ELSE IF (SUBKEY.EQ.'END') THEN
          ICOUNT = ICOUNT + 1
C     
C     ************************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
          IF (IOERR) INERR = .TRUE.
          PHIFIN(NSEGM) = VALUE(ICOUNT) + ISTRUN*360.0
C     
C---- Values must be integral
C     
          TEMPX = REAL(NINT(PHIFIN(NSEGM)))
          IF (ABS(PHIFIN(NSEGM)-TEMPX).GT.0.01) THEN
            WRITE(IOUT,FMT=6473) PHIFIN(NSEGM),TEMPX
            IF (ONLINE) WRITE(ITOUT,FMT=6473) PHIFIN(NSEGM),TEMPX
            PHIFIN(NSEGM) = TEMPX
          END IF
C     
C---- STEP phi
C     
        ELSE IF (SUBKEY.EQ.'STEP') THEN
          ICOUNT = ICOUNT + 1
C     
C     ************************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
          IF (IOERR) INERR = .TRUE.
          PHIINC(NSEGM) = ABS(VALUE(ICOUNT))
C     
C---- Values must be integral
C     
          TEMPX = REAL(NINT(PHIINC(NSEGM)))
          IF (ABS(PHIINC(NSEGM)-TEMPX).GT.0.01) THEN
            WRITE(IOUT,FMT=6473) PHIINC(NSEGM),TEMPX
            IF (ONLINE) WRITE(ITOUT,FMT=6473) PHIINC(NSEGM),TEMPX
            PHIINC(NSEGM) = TEMPX
          END IF
C     
C---- FIRSTPACK... Only appropriate for unique option. required if more
C     than one oscgen run is necessary to generate all
C     reflection data (eg using different crystal
C     orientations). in this case, for each run of oscgen
C     firstpack should be set to the final pack number
C     of the preceeding run +1.
C     
        ELSE IF (SUBKEY.EQ.'FIRS') THEN
C     
          ICOUNT = ICOUNT + 1
C     ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
          IF (IOERR) INERR = .TRUE.
          IFIRSTONE(NSEGM) = NINT(VALUE(ICOUNT))
C     
C---- PARTS or RUNS ... number of different parts (eg with different
C     crystal 
C     orientation)
C     
        ELSE IF ((SUBKEY.EQ.'RUNS').OR.(SUBKEY.EQ.'PART')) THEN
C     
          ICOUNT = ICOUNT + 1
C     ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          IF (IOERR) INERR = .TRUE.
          NSTRUNO = NSTRUN
          NSTRUN = NINT(VALUE(ICOUNT))
C     
C---- Check for consistency with previous input value (if any)
C     
          IF ((NSTRUNO.NE.0).AND.(NSTRUN.NE.NSTRUNO)) THEN
            WRITE(IOUT,FMT=6488) NSTRUN,NSTRUNO
            IF (ONLINE) WRITE(ITOUT,FMT=6488) NSTRUN,NSTRUNO
 6488       FORMAT(1X,'**** ERROR *****',/,1X,'RUNS has been',
     +           ' given as',I3,' while on a previous STRATEGY',
     +           ' keyword it was given as',I3)
            STOP
          END IF
C
C---- If running from new GUI and there is more than one part, set WAITINP
C     TRUE to avoid it trying to read further input from stored lines of input
C     after this STRATEGY keyword. 9dec05 AL
C
          IF ((NSTRUN.GT.1).AND.(SOCKLO)) WAITINP = .TRUE.
          IF (SOCKLO) NRUN = 0
C  
C---- AUTO
C     
        ELSE IF (SUBKEY.EQ.'AUTO') THEN
C     
          AUTO = .TRUE.
          IAUTO = 1
          IF (IISTART.EQ.1) THEN
            WRITE(IOUT,FMT=6486)
            IF (ONLINE) WRITE(ITOUT,FMT=6486)
            INERR = .TRUE.
          END IF
C     
C---- ANOM (Maximise anomalous pairs)
C     
        ELSE IF (SUBKEY.EQ.'ANOM') THEN
          AUTANOM = .TRUE.
C     
C---- NOTANOM (switch off maximise anomalous pairs)
C     
        ELSE IF (SUBKEY.EQ.'NOTA') THEN
          AUTANOM = .FALSE.
C     
C---- ROTATE (AUTO MODE)
C     
        ELSE IF (SUBKEY.EQ.'ROTA') THEN
C     
          ICOUNT = ICOUNT + 1
C     ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          IF (IOERR) INERR = .TRUE.
          ROTAUTO = NINT(VALUE(ICOUNT))
C     
C---- SEGMENTS (AUTO MODE)
C     
        ELSE IF (SUBKEY.EQ.'SEGM') THEN
C     
          ICOUNT = ICOUNT + 1
C     ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          IF (IOERR) INERR = .TRUE.
          NSEGAUTO = NINT(VALUE(ICOUNT))
          IF (NSEGAUTO.GT.4) THEN
            WRITE(IOUT,FMT=6377)
            IF (ONLINE) WRITE(ITOUT,FMT=6377) 
 6377       FORMAT(/,1X,'***** ERROR *****',/,1X,'Maximum number',
     +           ' of segments is 4. Segments reset to 4')
            NSEGAUTO = 4
          END IF
C     
C---- SIZES of SEGMENTS (AUTO MODE)
C     
        ELSE IF (SUBKEY.EQ.'SIZE') THEN
C     
          SIZESET = .TRUE.
 190      IISIZE = IISIZE + 1
          ICOUNT = ICOUNT + 1
C     ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          IF (IOERR) INERR = .TRUE.
          PHISEGA(IISIZE) = NINT(VALUE(ICOUNT))
          IF (ICOUNT.LT.NTOK) THEN
            IF (ITYP(ICOUNT+1).EQ.2) GOTO 190
          END IF
C     
C---- SPEEDUP (speedup factor)
C     
        ELSE IF (SUBKEY.EQ.'SPEE') THEN
C     
          INSPEED = 1
          IF (ICOUNT.LT.NTOK) THEN
            IF (ITYP(ICOUNT+1).EQ.2) THEN
              ICOUNT = ICOUNT + 1
C     ************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
              IF (IOERR) INERR = .TRUE.
              VOLSCAL = VALUE(ICOUNT)
            ELSE
              VOLSCAL = 5
            END IF
          ELSE
            VOLSCAL = 5
          END IF
C     
C---- Check for consistency
C     
          IF ((NSEGM.GT.1).AND.(NSTRUN.GT.1)) THEN
            IF (VOLSCAL.NE.OVOLSCAL) THEN
              WRITE(IOUT,FMT=7140) OVOLSCAL
              IF (ONLINE) WRITE(ITOUT,FMT=7140) OVOLSCAL
 7140         FORMAT(//,1X,'***** ERROR *****',/,1X,'You cannot ',
     +             'specify different SPEEDUP factors for ',
     $             'different runs',/,1X,'The original value of',
     $             F5.1,' has been restored')
              VOLSCAL = OVOLSCAL
            END IF
          END IF
          OVOLSCAL = VOLSCAL
C     
        ELSE IF (SUBKEY.EQ.'VOLU') THEN
          ICOUNT = ICOUNT + 1
C     ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          IF (IOERR) INERR = .TRUE.
          DEFVOL = VALUE(ICOUNT)
          
        ELSE IF (SUBKEY(1:3).EQ.'ALT') THEN
          NEWSTRAT = .TRUE.
        ELSE IF (SUBKEY.EQ.'OLD') THEN
          NEWSTRAT = .FALSE.
C
C---- run testgen automatically with range(s) found in this STRATEGY run
C     - useful for automated runs, also for Mosflm Server & expert system.
C
        ELSE IF (SUBKEY.EQ.'TEST') THEN
          TESTRAT = .TRUE.
          XOVER(1) = -10.0
C
C---- TESTGEN, PHSTART, PHEND, XOVER(1) _must_ be set in COMPLETE
C
C---- new S/R to allow DNA control over parameters
C
        ELSE IF(TESTRAT.AND.ICOUNT.LE.NTOK)THEN
          CALL KWTEST(LINE,ICOMM,ITINS,
     $         COMREAD,IERR)
          ICOUNT = NTOK
          IF(IERR.EQ.50)IERR = 0
C     Not recognised
C     
        ELSE
          IF (TRAPERR) INPERR = .TRUE.
          WRITE (IOUT,FMT=6130) SUBKEY
          IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          INERR = .TRUE.
          BADLINE = .TRUE.
        END IF
        IF (ICOUNT.LT.NTOK) GOTO 187
C     
C---- Trap an error in input
C     
        IF (INERR) THEN
          INERR = .FALSE.
          WRITE(IOUT,FMT=7160)
          IF (ONLINE) WRITE(ITOUT,FMT=7160)
 7160     FORMAT(1X,'*** Because of input error, this line has',
     +         ' been ignored ***')
          STRATEGY = .FALSE.
          NSEGM = NSEGM - 1
          NSEGRD = NSEGRD - 1
          IF (ONLINE) THEN
            IF (COMREAD) THEN
              COMREAD = .FALSE.
              ITIN = ITINS
              CLOSE (UNIT=ICOMM)
            END IF
            BADLINE = .TRUE.
            GOTO 50
          END IF
          STOP
        END IF
C     
C---- If no START keyword given, assume AUTO mode
C     
 192    IF (IIPHI.EQ.0) THEN
          AUTO = .TRUE.
          IAUTO = 1
        END IF
C     
C---- Phistart must be in range 0 to 359. Don't check if AUTO mode and
C     phistart has not been given
C     
        IF (AUTO.AND.(IIPHI.EQ.0)) GOTO 189
        IF ((PHIST(NSEGM)-ISTRUN*360).LT.0) THEN
          PHIST(NSEGM) = PHIST(NSEGM) + 360.0
          PHIFIN(NSEGM) = PHIFIN(NSEGM) + 360.0
          PHIADD(NSEGM) = 360.0
        END IF
        IF ((PHIST(NSEGM)-ISTRUN*360).GT.360) THEN
          PHIST(NSEGM) = PHIST(NSEGM) - 360.0
          PHIFIN(NSEGM) = PHIFIN(NSEGM) - 360.0
          PHIADD(NSEGM) = -360.0
        END IF
C     
C---- If no STEP given, use 5 degrees
C     
 189    IF (PHIINC(NSEGM).EQ.0) PHIINC(NSEGM) = 5.0
        IF (PHIFIN(NSEGM).LT.PHIST(NSEGM)) THEN
          TEMPX = PHIST(NSEGM)
          PHIST(NSEGM) = PHIFIN(NSEGM)
          PHIFIN(NSEGM) = TEMPX
        END IF
C     
        IF (NTOK.EQ.1) THEN
          WRITE(IOUT,FMT=7560)
          IF (ONLINE) WRITE(ITOUT,FMT=7560)
 7560     FORMAT('Strategy option will be run in automatic mode.',
     +         /,1X,'Speedup factor will be calculated ',
     +         'automatically (use keyword SPEEDUP to set ',
     +         'explicitly).',/,1X,'Type GO to continue, or',
     +         ' ABORT to stop strategy run.')
          IF (MODE.EQ.10) THEN
            LINE = ' '
            WRITE(LINE,FMT=7562)
 7562       FORMAT('Running strategy in default mode. Type GO ',
     +           'to continue or ABORT to stop.')
            CALL MXDWIO(LINE,2)
          END IF
        END IF
      ELSE IF ((KEY.EQ.'ABOR').AND.STRATEGY) THEN
        IFSTRAT = .FALSE.
        STRATEGY = .FALSE.
        NSEGM = NSEGM - 1
        NSEGRD = NSEGRD - 1
        WRITE(IOUT,FMT=7566)
        IF (ONLINE) WRITE(ITOUT,FMT=7566)
 7566   FORMAT(1X,'STRATEGY run aborted.')
        IF (MODE.EQ.10) THEN
          POWDER = .TRUE.
          CALL MXDCIO(1,0,0,0,0)
          RETURN
        END IF
C     
C     
C---- Ensure final value is an integral number of steps from start
C     *** NO LONGER NECESSARY ***
C     AL          I = NINT((PHIFIN(NSEGM)-PHIST(NSEGM))/PHIINC(NSEGM))
C     AL          PHIFIN(NSEGM) = PHIST(NSEGM) + I*PHIINC(NSEGM)
C     
C---- Symmetry, compulsory
C     
      ELSE IF (KEY.EQ.'SYMM') THEN
C     
C---- Trap symmetry given as 0
C     
        IF ((NTOK.EQ.2).AND.(ITYP(2).EQ.2).AND.
     +       (NINT(VALUE(2)).EQ.0)) THEN
          NUMSPG = 0
          LSYMM = 0
C
C---- next is necessary for multiple indexing runs in same job _if_
C     symmetry is set to 0
C
          CALL MATZERR(CELL,6)
          CALL CLEAR(AMAT)
          CALL CLEAR(BMAT)
          CALL IDNMAT(UMAT)
          IMAT = 0
          ICELL = 0
          IUMAT = 0
          GOTO 141
        END IF
cdms21.09.2004  Convert rhombohedral R# to H#, unless CRYSTAL RHOMBOHEDRAL
c               has been given
C
        IF (ITYP(2).NE.2) THEN
           IF (LINE(IBEG(2):IBEG(2)).EQ.'R') THEN
              IF (ICRYST.NE.8) LINE(IBEG(2):IBEG(2)) = 'H'
           ENDIF
        ENDIF
C     
        SYMMIN = .TRUE.
        LSYMM = 1
        CALL MRDSYMM(2,LINE,IBEG,IEND,ITYP,VALUE,NTOK,
     +       SPGNAM,NUMSPG,PGNAME,NSYM,NSYMP,RSYM,INERR)
        IF (INERR) THEN
          WRITE(IOUT,FMT=7142) NUMSPG,SPGNAM
          IF (ONLINE) WRITE(ITOUT,FMT=7142) NUMSPG,SPGNAM
 7142     FORMAT(1X,'*** ERROR *** Spacegroup number',I5,' name ',A,
     +         ' not found in SYMOP library')
          IF (ONLINE) THEN
            IF (COMREAD) THEN
              COMREAD = .FALSE.
              ITIN = ITINS
              CLOSE (UNIT=ICOMM)
            END IF
            BADLINE = .TRUE.
            GOTO 50
          END IF
            NSHUTERR = 20
          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
        END IF
        IF((NUMSPG.EQ.146).OR.(NUMSPG.EQ.148).OR.(NUMSPG.EQ.155).OR.
     +       (NUMSPG.EQ.160).OR.(NUMSPG.EQ.161).OR.(NUMSPG.EQ.166).OR.
     +       (NUMSPG.EQ.167)) THEN
          WRITE(IOUT,FMT=6472)
          IF (ONLINE) WRITE(ITOUT,FMT=6472)
 6472     FORMAT(1X,'**** Warning ****',/,1X,'For rhombohedral ',
     +         'spacegroups a hexagonal cell is used by default.',
     +         /,80('='))
        END IF
        IF((NUMSPG.EQ.1146).OR.(NUMSPG.EQ.1148).OR.(NUMSPG.EQ.1155).OR.
     +       (NUMSPG.EQ.1160).OR.(NUMSPG.EQ.1161).OR.(NUMSPG.EQ.1166)
     +       .OR.(NUMSPG.EQ.1167)) THEN
C
C---- set ICRYST here for R setting.
C
          ICRYST = 8

          WRITE(IOUT,FMT=6475)
          IF (ONLINE) WRITE(ITOUT,FMT=6475)
cdms21.09.2004   Change message, since primitive rhombohedral is now OK.
 6475     FORMAT(1X,'**** Warning ****',/,1X,'For rhombohedral ',
     +         'spacegroups a hexagonal cell is used by default.',
     +         /,1X,'You have requested a primitive rhombohedral ',
     +         'cell.',/,80('='))
        END IF
        IF (NUMSPG.LT.3) THEN
C     
C---- Triclinic (but trap case where input is CRYST R  then SYMM 0)
C     
          IF (ICRYST.NE.8) ICRYST = 1
        ELSE IF (NUMSPG.LT.16) THEN
C     
C---- Monoclinic
C     
          ICRYST = 2
        ELSE IF (NUMSPG.LT.75) THEN
C     
C---- Orthorhombic
C     
          ICRYST = 3
        ELSE IF (NUMSPG.LT.143) THEN
C     
C---- Tetragonal
C     
          ICRYST = 4
        ELSE IF (NUMSPG.LT.168) THEN
C     
C---- Trigonal, but allow for rhombohedral settings (CRYST keyword)
C     
          IF (ICRYST.NE.8) ICRYST = 5
        ELSE IF (NUMSPG.LT.195) THEN
C     
C---- Hexagonal
C     
          ICRYST = 6
        ELSE IF (NUMSPG.LE.230) THEN
C     
C---- Cubic
C     
          ICRYST = 7
        END IF
C     
C---- Set cell refinement flags
C     
        DO 140 I = 1,6
          LCELL(I) = LCLASS(I,ICRYST)
 140    CONTINUE
C     
C---- Get lattice type from spacegroup name
C     
        LATTYP = SPGNAM(1:1)    
C     
C     TAKANORI: something is very wrong with MRDSYMM!
C       C222 => C222, B2 => B2, A2 => A2, but
C       B222 turns into C222 and A222 into C222!!
C
c       print*, 'SPGNAME = ', SPGNAM, 'LATTYP = ', LATTYP
        ISYS = 2
        IF (LATTYP.EQ.'A') THEN
          KSYS(1) = 0
          KSYS(2) = 1
          KSYS(3) = 1
        ELSE IF (LATTYP.EQ.'B') THEN
          KSYS(1) = 1
          KSYS(2) = 0
          KSYS(3) = 1
        ELSE IF (LATTYP.EQ.'C') THEN
          KSYS(1) = 1
          KSYS(2) = 1
          KSYS(3) = 0
        ELSE IF (LATTYP.EQ.'I') THEN
          ISYS = 2
          KSYS(1) = 1
          KSYS(2) = 1
          KSYS(3) = 1
        ELSE IF (LATTYP.EQ.'R') THEN
C     
C---- Allow for choice of rhombohedral cell
C     
          IF (ICRYST.EQ.8) THEN
            ISYS = 0
          ELSE
C     
C---- hexagonal setting
C     
            ISYS = 3
            KSYS(1) = -1
            KSYS(2) = 1
            KSYS(3) = 1
            LATTYP = 'H'
          END IF
        ELSE IF (LATTYP.EQ.'H') THEN
          ISYS = 3
          KSYS(1) = -1
          KSYS(2) = 1
          KSYS(3) = 1
        ELSE IF (LATTYP.EQ.'F') THEN
          KSYS(1) = 1
          KSYS(2) = 1
          KSYS(3) = 1
          ISYS = 4
        ELSE IF (LATTYP.EQ.'P') THEN
          ISYS = 0
          KSYS(1) = 1
          KSYS(2) = 1
          KSYS(3) = 1
        END IF
C     
C---- If a CELL keyword has been given already, impose symmetry
C     constraints
C     on the cell.
C     
        IF ((.NOT.ANGLES).AND.((ICRYST.EQ.5).OR.(ICRYST.EQ.6)))
     +       CELL(6) = 120.0
C     
 141    CONTINUE
C     
      ELSE IF (KEY.EQ.'TEST') THEN
        IF(HARRY)THEN
          CALL KWTEST(LINE,ICOMM,ITINS,
     $         COMREAD,IERR)
          IF(IERR.EQ.50)THEN
            IERR = 0
            BADLINE = .TRUE.
            GOTO 50
          ENDIF
        ELSE
C     
C---- TESTGEN option
C    
c     -harvest
          DOHARVEST = .false.
c     -harvest
          TESTGEN = .TRUE.
          OSCANG = 0.0
          ISTAFLG = 0
          IENDFLG = 0
          ICOUNT = 1
          IF (NTOK.EQ.1) GOTO 185
 188      ICOUNT = ICOUNT + 1
          SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     **************
          CALL CCPUPC(SUBKEY)
C     **************
          IF (SUBKEY.EQ.'STAR') THEN
            ISTAFLG = 1
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            PHSTART = VALUE(ICOUNT)
            IF (ABS(PHSTART-NINT(PHSTART)).GT.0.01) THEN
              WRITE(IOUT,FMT=6474) NINT(PHSTART)
              IF (ONLINE) WRITE(ITOUT,FMT=6474) NINT(PHSTART)
 6474         FORMAT(1X,'*** WARNING ***',/,1X,'Phi values and step',
     +             ' must be integers, nearest integer',I5,' taken')
            END IF
          ELSE IF (SUBKEY.EQ.'END') THEN
          IENDFLG = 1
          ICOUNT = ICOUNT + 1
C     
C     ************************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
          PHEND = VALUE(ICOUNT) 
          IF (ABS(PHEND-NINT(PHEND)).GT.0.01) THEN
            WRITE(IOUT,FMT=6474) NINT(PHEND)
            IF (ONLINE) WRITE(ITOUT,FMT=6474) NINT(PHEND)
          END IF
        ELSE IF (SUBKEY.EQ.'STEP') THEN
          ICOUNT = ICOUNT + 1
C     
C     ************************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
          PHSTEP = VALUE(ICOUNT) 
          IF (ABS(PHSTEP-NINT(PHSTEP)).GT.0.01) THEN
            WRITE(IOUT,FMT=6474) NINT(PHSTEP)
            IF (ONLINE) WRITE(ITOUT,FMT=6474) NINT(PHSTEP)
          END IF
        ELSE IF (SUBKEY.EQ.'ANGL') THEN
          ICOUNT = ICOUNT + 1
C     
C     ************************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
          OSCANG = VALUE(ICOUNT) 
        ELSE IF (SUBKEY.EQ.'MINO') THEN
          ICOUNT = ICOUNT + 1
C     
C     ************************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
          OSCMIN = VALUE(ICOUNT) 
        ELSE IF (SUBKEY.EQ.'MAXO') THEN
          ICOUNT = ICOUNT + 1
C     
C     ************************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
          OSCMAX = VALUE(ICOUNT) 
        ELSE IF (SUBKEY.EQ.'OVER') THEN
          ICOUNT = ICOUNT + 1
C     
C     ************************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
          XOVER(1) = VALUE(ICOUNT)
C     
C     Not recognised
C     
        ELSE
          IF (TRAPERR) INPERR = .TRUE.
          WRITE (IOUT,FMT=6130) SUBKEY
          IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
        END IF
        IF (ICOUNT.LT.NTOK) GOTO 188
C     
C---- Test all values have been given
C     
 185    IF ((ISTAFLG.EQ.0).OR.(IENDFLG.EQ.0)) THEN
          WRITE(IOUT,FMT=7130)
          IF (ONLINE) WRITE(ITOUT,FMT=7130)
          IF (ONLINE) THEN
            IF (COMREAD) THEN
              COMREAD = .FALSE.
              ITIN = ITINS
              CLOSE (UNIT=ICOMM)
            END IF
            BADLINE = .TRUE.
            GOTO 50
          END IF
          NSHUTERR = 21
          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
        END IF
 7130   FORMAT(1X,'*** ERROR ***',/,1X,'START and END keywords',
     +       ' must be given',/,1X,'Full list of possible keywords',
     +       /,1X,'TESTGEN START 0 END 90 STEP 5 OVERLAP 4 MINOSC',
     +       ' 0.3 MAXOSC 4',/,1X,'This will test phi values from',
     +       ' 0 to 90 in steps of 5 degrees. At each phi value',/,1X,
     +       'The oscillation angle giving less than 4% overlapped',
     +       ' spots will be determined',/,1X,'providing this is',
     +       ' between 0.3 and 4 degrees')
        ENDIF
C     .. end of KWTEST block ..
C     
C     
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**********************************************************************
C---- Genfile ?
C     
      ELSE IF (KEY.EQ.'GENF') THEN
C     
C---- Must trap a "null" filename
C     
        IF (NTOK.EQ.1) THEN
          WRITE(IOUT,FMT=7060) LINE(IBEG(1):IEND(1))
          IF (ONLINE) WRITE(ITOUT,FMT=7060) LINE(IBEG(1):IEND(1))
          BADLINE = .TRUE.
          GOTO 50
        END IF
C     
C---- quoted token
C     
        IF (ITYP(2).EQ.3) THEN
          GENFILE = LINE(IBEG(2) :IEND(2))
        ELSE
          GENFILE = LINE(IBEG(2) :IEND(NTOK))
        END IF
C     
C---- Append .gen if not specified
C     find how many non-blank characters in genfile
C     
        NCH = LENSTR(GENFILE)
C     
C     
        DO 200 I = 1,NCH
          IF (GENFILE(I:I).EQ.'.') GO TO 202
 200    CONTINUE 
C     
C     
        GENFILE = GENFILE(1:NCH)//'.gen'
C     
 202    NEWGENF = .TRUE.
        IGENF = 1
C     
C---- HKLOUT...output mtz file
C     
      ELSE IF (KEY.EQ.'HKLO') THEN
C     
C---- Must trap a "null" filename
C     
        IF (NTOK.EQ.1) THEN
          WRITE(IOUT,FMT=7060) LINE(IBEG(1):IEND(1))
          IF (ONLINE) WRITE(ITOUT,FMT=7060) LINE(IBEG(1):IEND(1))
 7060     FORMAT(/,1X,'*** ERROR ***',/,1X,'Must supply a parameter',
     +         ' for keyword: ',A)
          BADLINE = .TRUE.
          GOTO 50
        END IF
C     
C---- quoted token
C     
        IHKLOUT = 1
        IF (ITYP(2).EQ.3) THEN
          MTZNAM = LINE(IBEG(2) :IEND(NTOK))
        ELSE
          MTZNAM = LINE(IBEG(2) :IEND(2))
        END IF
C     
C---- Trap case where only subkeyword MULTIPLE has been given and NO
C     MTZ filename, OR MULTIPLE has been given BEFORE the MTZ filename
C     
        STR2 = MTZNAM
        CALL CCPUPC(STR2)
        IF ((STR2(1:8).EQ.'MULTIPLE').OR.
     +       (STR2(1:10).EQ.'NOMULTIPLE')) THEN
          IF (STR2(1:8).EQ.'MULTIPLE') MULTIMTZ = .TRUE.
          IF (STR2(1:8).EQ.'NOMULTIPLE') MULTIMTZ = .FALSE.
          IF (NTOK.EQ.2) THEN
            IHKLOUT = 0
            MTZNAM = 'HKLOUT'
            GOTO 207
          ELSE IF (NTOK.EQ.3) THEN
            MTZNAM = LINE(IBEG(3):IEND(3))
C     
C-----Set flag so it does not check 3rd token below
C     
            J = 999
          END IF
        END IF
C     
C---- Append .mtz if not specified
C     find how many non-blank characters in MTZNAM
C     
        NCH = LENSTR(MTZNAM)
C     
C     
        DO 204 I = 1,NCH
          IF (MTZNAM(I:I).EQ.'.') GO TO 206
 204    CONTINUE 
C     
C     
        MTZNAM = MTZNAM(1:NCH)//'.mtz'
C     
 206    CONTINUE
C     
C---- Test for MULTI subkeyword as 3rd token
C     
        IF ((NTOK.EQ.3).AND.(J.NE.999)) THEN
          SUBKEY = LINE(IBEG(3):IEND(3))
          CALL CCPUPC(SUBKEY)
          IF (SUBKEY.EQ.'MULT') THEN
            MULTIMTZ = .TRUE.
          ELSE IF (SUBKEY.EQ.'NOMU') THEN
            MULTIMTZ = .FALSE.
          END IF
        END IF
 207    J = 0
C
C---- MTZ directory
C     
      ELSE IF (KEY.EQ.'MTZD') THEN
C     
C---- Must trap a "null" filename
C     
        IF (NTOK.EQ.1) THEN
          WRITE(IOUT,FMT=7060) LINE(IBEG(1):IEND(1))
          IF (ONLINE) WRITE(ITOUT,FMT=7060) LINE(IBEG(1):IEND(1))
          BADLINE = .TRUE.
          GOTO 50
        END IF
C     
C---- quoted token
C     
        IHKLDIR = 1
        IF (ITYP(2).EQ.3) THEN
          MTZDIR = LINE(IBEG(2) :IEND(NTOK))
        ELSE
          MTZDIR = LINE(IBEG(2) :IEND(2))
        END IF
C     
C---- If running under Unix, check for a "/" at end of directory
C     specification, and add one if not present.
C     
          IF (.NOT.VAXVMS()) THEN
            NCH = LENSTR(MTZDIR)
            IF (MTZDIR(NCH:NCH).NE.'/') THEN
              MTZDIR = MTZDIR(1:NCH)//'/'
            END IF
          END IF

C - also need some checks to see if it exists etc...
C     
    
C
C---- MTZ directory
C     
      ELSE IF (KEY.EQ.'CLOS') THEN
        IF (MTZOPEN) THEN
          MTZOPEN = .FALSE.
          MTZPRT = 1
          CALLED = 'LWCLOS'
          COMMENT = 'To close MTZ file '
          IF (DEBUG(99)) 
     +        CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
C     *********************
          CALL LWCLOS(MTZOUT,MTZPRT)
          if (debug(52))then
            write(iout,fmt=7742)
            if (online) write(itout,fmt=7742)
 7742       format(1x,'**closing MTZ file in control 2 ***')
          end if
          IF(LBEST.AND.BESTOPEN)THEN
            CLOSE(BESTHKL)
            BESTOPEN = .FALSE.
          ENDIF
C     *********************
        END IF
        MINBATCH = -999
        NRUN = 0
        
C---- Site
C     
      ELSE IF (KEY.EQ.'SITE') THEN
C     
C     
        IF (NTOK.LE.1) THEN
          WRITE(IOUT,FMT=6006)
 6006     FORMAT(' *** ERROR *** No value given for Key_Word SITE ')
          IF (ONLINE) WRITE(ITOUT,FMT=6006)
          IF (ONLINE) THEN
            IF (COMREAD) THEN
              COMREAD = .FALSE.
              ITIN = ITINS
              CLOSE (UNIT=ICOMM)
            END IF
            BADLINE = .TRUE.
            GOTO 50
          ELSE
            NSHUTERR = 22
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF
        END IF
C     
C     
        SITE = LINE(IBEG(2):IEND(2))
        IISITE = .TRUE.
        SVSITE = LINE(IBEG(2):IEND(NTOK))
C     
C     *************
        CALL CCPUPC(SITE)
C     *************
C---- Test against known sites
C     
        IF (SITE.EQ.'EMBL'.OR.SITE(1:3).EQ.'LMB'.OR.SITE.EQ.'DLAB'
     +       .OR.SITE.EQ.'IMPC'.OR.SITE.EQ.'CHES'.OR.SITE.EQ.'SSRL'.OR.
     +       SITE.EQ.'ALS') THEN
          CONTINUE
        ELSE         
          WRITE(IOUT,FMT=6007) SITE
          IF (ONLINE) WRITE(ITOUT,FMT=6007) SITE
          IF (ONLINE) THEN
           IF (COMREAD) THEN
              COMREAD = .FALSE.
              ITIN = ITINS
              CLOSE (UNIT=ICOMM)
            END IF
            BADLINE = .TRUE.
            GOTO 50
          ELSE
            NSHUTERR = 23
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF
        END IF
 6007   FORMAT(1X,'The site "',A,'" is not known, must be one of:',
     +       /,1X,'EMBL,LMB,DLAB,IMPC,CHESS,SSRL,ALS')
C     
C---- If EMBL, test for further keyword, defaults to SCR3
C     
        IF (SITE.EQ.'EMBL') THEN
          ISCAN = 1
          MACHINE = 'MAR '
          COMPRESSION = 1
          IF (NTOK.EQ.2) THEN
            SCANNER = 'SCR3'
          ELSE   
            SCANNER = LINE(IBEG(3):IEND(3))
C     *************
            CALL CCPUPC(SCANNER)
C     *************
          END IF
          IF (IEXTEN.EQ.0) ODEXT = 'corr'
C     
C---- Image plate version for Hamburg, SCR1,SCR2,SCR3 or MAR.
C     This is used to set up default beam centre and image size
C     (including
C     size of header block)
C     
          IF (SCANNER.EQ.'SCR1') THEN
            VERS1= .TRUE.
            IF(VERS2.OR.VERS3)
     +           STOP 'MUST USE ONLY ONE OF SCR1,SCR2,SCR3 CARDS '
            WRITE(IOUT,6005) VERSTR(1)
            IF (ONLINE)WRITE(ITOUT,6005)
 6005       FORMAT(1X,'*** HAMBURG scanner ',A,'***' /)
C     
            NREC = 1187
            IYLEN = 1187             
            RAST = 0.187
            XMAXIP = XMAXRED
            YMAXIP = YMAXRED
            RMAXIP = RMAXRED
            RSCANIP = RSCANRED
          ELSE IF (SCANNER.EQ.'SCR2') THEN
            VERS2= .TRUE.
            NHEAD = 0
            HDRSIZE = .FALSE.
            NREC = 1187
            IYLEN = 1187             
            RAST = 0.187
            XMAXIP = XMAXRED
            YMAXIP = YMAXRED
            RMAXIP = RMAXRED
            RSCANIP = RSCANRED
            IF(VERS1.OR.VERS3)
     +           STOP 'MUST USE ONLY ONE OF SCR1,SCR2,SCR3 CARDS '
            WRITE(IOUT,6005) VERSTR(2)
            IF (ONLINE)WRITE(ITOUT,6005) VERSTR(2)
C     
          ELSE IF (SCANNER.EQ.'SCR3') THEN
            VERS3= .TRUE.
            NHEAD = 0
            HDRSIZE = .FALSE.
            NREC = 1187
            IYLEN = 1187             
            RAST = 0.187
            XMAXIP = XMAXRED
            YMAXIP = YMAXRED
            RMAXIP = RMAXRED
            RMINIP = RMINRED
            RSCANIP = RSCANRED
            IF(VERS1.OR.VERS2)
     +           STOP 'MUST USE ONLY ONE OF SCR1,SCR2,SCR3 CARDS '
            WRITE(IOUT,6005) VERSTR(3)
            IF (ONLINE)WRITE(ITOUT,6005) VERSTR(3)
C     
C---- Mar Research scanner at EMBL (add any required code here)
          ELSE IF (SCANNER.EQ.'MAR') THEN
            CONTINUE
          ELSE
C     
C---- Not recognised
            IF (TRAPERR) INPERR = .TRUE.
            WRITE(IOUT,6008) SCANNER
            IF (ONLINE) WRITE(ITOUT,6008) SCANNER
            IF (ONLINE) THEN
              IF (COMREAD) THEN
                COMREAD = .FALSE.
                ITIN = ITINS
                CLOSE (UNIT=ICOMM)
              END IF
              BADLINE = .TRUE.
              GOTO 50
            ELSE       
              STOP
            END IF
 6008       FORMAT(1X,'*** ERROR ***',/,1X,'Version "',A,
     $           '" not known',/,1X,'Must be one of ',
     $           'SCR1,SCR2,SCR3')
          END IF
C     
C     LMB Prototype scanner, 1187*1187 pixels, pixel 0.150
C     
        ELSE IF (SITE.EQ.'LMB') THEN
          ISCAN = 1
          NHEAD = 0
          HDRSIZE = .FALSE.
          NREC = 1187
          IYLEN = 1187             
          RAST = 0.150
          RSCANIP = 8887
          RMAXIP = 8887
          XMAXIP = 8887
          YMAXIP = 8887
          MACHINE = 'MAR '
          COMPRESSION = 1
          IF (IEXTEN.EQ.0) ODEXT = 'pck'
C     
C---- CHESS ...eg SITE CHESS [A1 F1 F2] 
C     [FUJI [CCD [1K 2K 2KBINNED ADSC QUANTUM4]]]
C     eg SITE CHESS F1 CCD 2KBINNED
C     or SITE CHESS F1 FUJI
C     
        ELSE IF (SITE.EQ.'CHES') THEN
          IF (NTOK.LT.4) THEN
            WRITE(IOUT,FMT=7112)
            IF (ONLINE) WRITE(ITOUT,FMT=7112)
            NSHUTERR = 24
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF
 7112     FORMAT(1X,'*** ERROR ***',/,1X,'For CHESS, the station',
     +         ' and type of detector must be specified',
     +         /,1X,'eg SITE CHESS (a1 or F1 or F2) (FUJI or CCD)',
     +         'if CCD then (1K or 2K or 2KBINNED)')
C     
C---- get station
C     
          STATION = LINE(IBEG(3):IEND(3))
          CALL CCPUPC(STATION)
          Hbeamline='CHESS-' // STATION(1:lenstr(STATION))
          IF (STATION(1:2).EQ.'A1') THEN
          ELSE IF (STATION(1:2).EQ.'F1') THEN
          ELSE IF (STATION(1:2).EQ.'F2') THEN
          ELSE
            WRITE(IOUT,FMT=7111) STATION
            IF (ONLINE) WRITE(ITOUT,FMT=7111) STATION
 7111       FORMAT(1X,'*** ERROR ***',/,1X,'STATION ',A,'not ',
     +           'recognised, must be A1, F1 or F2',/,1X,'eg:',
     +           'SITE CHESS A1 CCD 2KBINNED')
            NSHUTERR = 25
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF

C     
C---- Get detector type
C     
          MODEL = LINE(IBEG(4):IEND(4))
C     
C---- Convert to upper case
C     
C     ***********
          CALL CCPUPC(MODEL)
C     ***********
          IF (MODEL(1:3).EQ.'CCD') THEN
            INVERTX = .FALSE.
            SPIRAL = .FALSE.
            CIRCULAR = .FALSE.
            ORTHOG = .TRUE.
            ISCAN = 1
            IF (STATION(1:2).EQ.'F1') THEN
              OMEGAFD = 0.0
              omegaf = omegafd * dtor
            ELSE
              OMEGAFD = 180.0
              omegaf = omegafd * dtor
            END IF
            MACHINE = 'CCD1'
            COMPRESSION = 1
            NULLPIX = 100
            RMINIP = 300
            IF (IEXTEN.EQ.0) ODEXT = 'tif'
C     
C---- Get type of CCD, can be 1K, 2K, 2KBINNED or ADSC

            IF (NTOK.LT.5) THEN
              WRITE(IOUT,FMT=7112)
              IF (ONLINE) WRITE(ITOUT,FMT=7112)
              NSHUTERR = 26
              CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
            END IF
C     
            KEY2 = LINE(IBEG(5):IEND(5))
C     ***********
            CALL CCPUPC(KEY2)
C     ***********
            IF (KEY2(1:2).EQ.'1K') THEN
              NREC = 1024
              IYLEN = 1024
              NHEAD = 0
              HDRSIZE = .FALSE.
              RAST = 0.0508
              IF (XSCAN.EQ.0) XSCAN = 2370
              IF (YSCAN.EQ.0) YSCAN = 2460
              XMAXIP = 2400
              YMAXIP = 2400
              RMAXIP = SQRT(XMAXIP**2 + YMAXIP**2)
              IF (ICUT.EQ.0) CUTOFF = 59900
              IF (IPRCUT.EQ.0) PRCUTOFF = 59900
              IF (IGAIN.EQ.0)THEN
                IGAIN = 1
                GAIN = 0.5
              ENDIF
            ELSE IF (KEY2(1:3).EQ.'2KB') THEN
              NREC = 1024
              IYLEN = 1024
              RAST = 0.0819
              NHEAD = 0
              HDRSIZE = .FALSE.
              IF (XSCAN.EQ.0) XSCAN = 3820
              IF (YSCAN.EQ.0) YSCAN = 3990
              XMAXIP = 3820
              YMAXIP = 3990
              RMAXIP = SQRT(XMAXIP**2 + YMAXIP**2)
              IF (ICUT.EQ.0) CUTOFF = 59900
              IF (IPRCUT.EQ.0) PRCUTOFF = 59900
              IF (IGAIN.EQ.0)THEN
                IGAIN = 1
                GAIN = 0.5
              ENDIF
            ELSE IF (KEY2(1:2).EQ.'2K') THEN
              NREC = 2048
              IYLEN = 2048
              NHEAD = 0
              HDRSIZE = .FALSE.
              RAST = 0.04095
              IF (XSCAN.EQ.0) XSCAN = 3820
              IF (YSCAN.EQ.0) YSCAN = 3990
              XMAXIP = 3820
              YMAXIP = 3990
              RMAXIP = SQRT(XMAXIP**2 + YMAXIP**2)
              IF (ICUT.EQ.0) CUTOFF = 14900
              IF (IPRCUT.EQ.0) PRCUTOFF = 14900
              IF (IGAIN.EQ.0)THEN
                IGAIN = 1
                GAIN = 0.5
              ENDIF
            ELSE IF (KEY2(1:2).EQ.'AD') THEN
              NREC = 1200
              IYLEN = 1200
              NHEAD = 0
              HDRSIZE = .FALSE.
              RAST = 0.0707
              IF (XSCAN.EQ.0) XSCAN = 3990
              IF (YSCAN.EQ.0) YSCAN = 3990
              XMAXIP = 3990
              YMAXIP = 3990
              RMAXIP = SQRT(XMAXIP**2 + YMAXIP**2)
              IF (ICUT.EQ.0) CUTOFF = 65500
              IF (IPRCUT.EQ.0) PRCUTOFF = 65500
              IF (IGAIN.EQ.0)THEN
                IGAIN = 1
                GAIN = 1.0
              ENDIF
              MACHINE = 'MAR '
              COMPRESSION = 1
              INVERTX = .TRUE.
              OMEGAFD = 90.0
              omegaf = omegafd * dtor
              IF(.NOT.IINULL)NULLPIX = 0
              IF (IEXTEN.EQ.0) ODEXT = 'image'
            ELSE IF (KEY2.EQ.'QUAN') THEN
              MACHINE = 'ADSC'
              COMPRESSION = 1
              MODEL = 'QUAD'
              IF (IEXTEN.EQ.0) ODEXT = 'img'
              INVERTX = .TRUE.
              SPIRAL = .FALSE.
              CIRCULAR = .FALSE.
              ORTHOG = .TRUE.
              TILED = .TRUE.
              ISCAN = 1
              OMEGAFD = 180.0
              omegaf = omegafd * dtor
              IF(.NOT.IINULL)NULLPIX = 0
              IF(IIDIVIDE.EQ.0)IDIVIDE = 0
              RMINIP = 300
              IF (XSCAN.EQ.0) XSCAN = 9400
              IF (YSCAN.EQ.0) YSCAN = 9400
              XMAXIP = 9400
              YMAXIP = 9400
              RMAXIP = SQRT(XMAXIP**2 + YMAXIP**2)
              IF (ILIMIT.EQ.0) THEN
                XLIMIT = 0.01*(XMAXIP/2)
                LIMIT = NINT(100*XLIMIT)
              END IF
              IF (ICUT.EQ.0) CUTOFF = 65500
              IF (IPRCUT.EQ.0) PRCUTOFF = 65500
              IF (IGAIN.EQ.0)THEN
                IGAIN = 1
C
C---- change default gain to 0.33, found on many detectors
C     AGWL 24/6/13
C
CAL                GAIN = 0.25
                GAIN = 0.33
              ENDIF
            ELSE
              IF (TRAPERR) INPERR = .TRUE.
              WRITE(IOUT,FMT=7108) KEY2
              IF (ONLINE) WRITE(ITOUT,FMT=7108) KEY2
 7108         FORMAT(1X,'*** ERROR ***',/,1X,'Detector ',
     +             A,' not recognised, current options are: ',
     +             '1K, 2KBINNED or 2K ')
              NSHUTERR = 27
              CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
            END IF
          ELSE IF (MODEL(1:4).EQ.'FUJI') THEN
            IF (STATION(1:2).EQ.'A1') THEN
              OMEGAFD = 0.0
              omegaf = omegafd * dtor
            ELSE
              OMEGAFD = 180.0
              omegaf = omegafd * dtor
            END IF
            ISCAN = 1
            INVERTX = .FALSE.
            NEWPREF = .TRUE.
            SPIRAL = .FALSE.
            CIRCULAR = .FALSE.
            ORTHOG = .TRUE.
            RAST = 0.10
            NHEAD = 0
            HDRSIZE = .FALSE.
            NREC = 2560
            IYLEN = 2048
            MACHINE = 'FUJI'
            COMPRESSION = 1
            LOGA = 4.0
            LOGB = 1023.0
          ELSE
            IF (TRAPERR) INPERR = .TRUE.
            WRITE(IOUT,FMT=7110) MODEL
            IF (ONLINE) WRITE(ITOUT,FMT=7110) MODEL
 7110       FORMAT(1X,'*** ERROR ***',/,1X,'Detector type ',
     +           A,' not recognised, current options are: CCD, FUJI')
            NSHUTERR = 28
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF
C     
C---- SSRL ...eg SITE SSRL ADSC
C     
        ELSE IF (SITE.EQ.'SSRL') THEN
          IF (NTOK.LT.3) THEN
            WRITE(IOUT,FMT=7480)
            IF (ONLINE) WRITE(ITOUT,FMT=7480)
            NSHUTERR = 29
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF
 7480     FORMAT(1X,'*** ERROR ***',/,1X,'For SSRL, the ',
     +         'type of detector must be specified',
     +         /,1X,'eg SITE SSRL ADSC')
C     
C---- Get detector type
C     
          MACHINE = LINE(IBEG(3):IEND(3))
C     
C---- Convert to upper case
C     
C     ***********
          CALL CCPUPC(MACHINE)
C     ***********
          IF (MACHINE.EQ.'ADSC') THEN
            COMPRESSION = 1
            MODEL = 'QUAD'
            IF (IEXTEN.EQ.0) ODEXT = 'img'
            INVERTX = .TRUE.
            SPIRAL = .FALSE.
            CIRCULAR = .FALSE.
            ORTHOG = .TRUE.
            TILED = .TRUE.
            ISCAN = 1
            OMEGAFD = 0.0
            omegaf = omegafd * dtor
            IF(.NOT.IINULL)NULLPIX = 0
            IF(IIDIVIDE.EQ.0)IDIVIDE = 0
            RMINIP = 300
            IF (XSCAN.EQ.0) XSCAN = 9400
            IF (YSCAN.EQ.0) YSCAN = 9400
            XMAXIP = 9400
            YMAXIP = 9400
            RMAXIP = SQRT(XMAXIP**2 + YMAXIP**2)
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
            IF (ICUT.EQ.0) CUTOFF = 65500
            IF (IPRCUT.EQ.0) PRCUTOFF = 65500
            IF (IGAIN.EQ.0)THEN
              IGAIN = 1
C
C---- change default gain to 0.33, found on many detectors
C     AGWL 24/6/13
C
CAL                GAIN = 0.25
              GAIN = 0.33
            ENDIF
          ELSE
            WRITE(IOUT,FMT=7481)
            IF (ONLINE) WRITE(ITOUT,FMT=7481)
 7481       FORMAT(1X,'*** ERROR ***',/,1X,'For SSRL, the ',
     +           'only type of detector allowed is ADSC',/,1X,
     +           'If using a Mar IP, do NOT use the SITE keyword.')
          END IF
C     
C---- ALS ...eg SITE ALS ADSC
C     
        ELSE IF (SITE.EQ.'ALS') THEN
          IF (NTOK.LT.3) THEN
            WRITE(IOUT,FMT=7482)
            IF (ONLINE) WRITE(ITOUT,FMT=7482)
            NSHUTERR = 30
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF
 7482     FORMAT(1X,'*** ERROR ***',/,1X,'For ALS, the ',
     +         'type of detector must be specified',
     +         /,1X,'eg SITE ALS ADSC')
C     
C---- Get detector type
C     
          MACHINE = LINE(IBEG(3):IEND(3))
C     
C---- Convert to upper case
C     
C     ***********
          CALL CCPUPC(MACHINE)
C     ***********
          IF (MACHINE.EQ.'ADSC') THEN
            COMPRESSION = 1
            MODEL = 'QUAD'
            IF (IEXTEN.EQ.0) ODEXT = 'img'
            INVERTX = .TRUE.
            SPIRAL = .FALSE.
            CIRCULAR = .FALSE.
            ORTHOG = .TRUE.
            TILED = .TRUE.
            ISCAN = 1
            OMEGAFD = 0.0
            omegaf = omegafd * dtor
            IF(.NOT.IINULL)NULLPIX = 0
            IF(IIDIVIDE.EQ.0)IDIVIDE = 0
            RMINIP = 300
            IF (XSCAN.EQ.0) XSCAN = 9400
            IF (YSCAN.EQ.0) YSCAN = 9400
            XMAXIP = 9400
            YMAXIP = 9400
            RMAXIP = SQRT(XMAXIP**2 + YMAXIP**2)
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
            IF (ICUT.EQ.0) CUTOFF = 65500
            IF (IPRCUT.EQ.0) PRCUTOFF = 65500
            IF (IGAIN.EQ.0)THEN
              IGAIN = 1
C
C---- change default gain to 0.33, found on many detectors
C     AGWL 24/6/13
C
CAL                GAIN = 0.25
              GAIN = 0.33
            ENDIF
          ELSE
            WRITE(IOUT,FMT=7483)
            IF (ONLINE) WRITE(ITOUT,FMT=7483)
 7483       FORMAT(1X,'*** ERROR ***',/,1X,'For ALS, the ',
     +           'only type of detector allowed is ADSC',/,1X,
     +           'If using a Mar IP, do NOT use the SITE keyword.')
          END IF
C     
C     
C---- Image plate version for DLAB, can be MAR (default) or RAXIS.
C     
        ELSE IF (SITE.EQ.'DLAB') THEN
          IF (NTOK.EQ.2) THEN
            MACHINE = 'MAR '
            COMPRESSION = 1
            ISCAN = 1
          ELSE   
            MACHINE = LINE(IBEG(3):IEND(3))
C     *************
            CALL CCPUPC(MACHINE)
C     *************
          END IF
C     
C---- Mar Research scanner at DLAB  (add any required code here)
          IF (MACHINE.EQ.'RAXI') THEN
            COMPRESSION = 3
            SPIRAL = .FALSE.
            CIRCULAR = .FALSE.
            ORTHOG = .TRUE.
chrp 20.08.2007             IF(.NOT.NEWRAXIS)
            INVERTX = .FALSE.
            FINE = .TRUE.
            RAST = 0.105
            ISCAN = 1
C     
C---- For the Daresbury R-axis, the rotation axis is parallel to the
C     fast scan direction, omega = 180
C     
            OMEGAFD = 180.0
            omegaf = omegafd * dtor
C     
C---- Test for keyword COARSE
C     
            IF (NTOK.GT.3) THEN
              SCAN = LINE(IBEG(4):IEND(4))
C     
C---- Convert to upper case
C     
C     ***********
              CALL CCPUPC(SCAN)
C     ***********
              IF (SCAN.EQ.'COAR') THEN
                FINE = .FALSE.
                RAST = 0.210
              ELSE IF (SCAN.EQ.'FINE') THEN
                FINE = .TRUE.
              ELSE
                IF (TRAPERR) INPERR = .TRUE.
                WRITE(IOUT,FMT=6700) SCAN
                IF (ONLINE) WRITE(ITOUT,FMT=6700) SCAN
                STOP
              END IF
            END IF
 6700       FORMAT(1X,'ERROR, Scan ',A,' not recognised',/,1X,
     +           'must be FINE or COARSE')
          ELSE IF (MACHINE.EQ.'MAR') THEN
            MACHINE = 'MAR '
            COMPRESSION = 1
            ISCAN = 1
            CONTINUE
          ELSE IF (MACHINE.EQ.'MD') THEN
            ISCAN = 1
            NHEAD = 0
            HDRSIZE = .FALSE.
            FINE = .TRUE.
            RAST = 0.08585
            SPIRAL = .FALSE.
            CIRCULAR = .FALSE.
            ORTHOG = .TRUE.
C     
C---- Test for keyword COARSE
C     
            IF (NTOK.GT.3) THEN
              SCAN = LINE(IBEG(4):IEND(4))
C     
C---- Convert to upper case
C     
C     ***********
              CALL CCPUPC(SCAN)
C     ***********
              IF (SCAN.EQ.'COAR') THEN
                FINE = .FALSE.
                RAST = 0.1720
              ELSE IF (SCAN.EQ.'FINE') THEN
                FINE = .TRUE.
              ELSE
                WRITE(IOUT,FMT=6700) SCAN
                IF (ONLINE) WRITE(ITOUT,FMT=6700) SCAN
                STOP
              END IF
            END IF

          ELSE
C     
C---- Not recognised
            IF (TRAPERR) INPERR = .TRUE.
            WRITE(IOUT,6710) MACHINE
            IF (ONLINE) WRITE(ITOUT,6710) MACHINE
            IF (ONLINE) THEN
              IF (COMREAD) THEN
                COMREAD = .FALSE.
                ITIN = ITINS
                CLOSE (UNIT=ICOMM)
              END IF
              BADLINE = .TRUE.
              GOTO 50
            ELSE       
              STOP
            END IF
 6710       FORMAT(1X,'*** ERROR ***',/,1X,'Scanner "',A,'" ',
     $           ' not known',/,1X,'Must be MAR or RAXIS')
          END IF
        END IF
C     
C---  DETECTOR or SCANNER or MACHINE
C     Subkeywords: ERROR, OFFLINE, [RAXIS, RAXIS4/RAXISIV, RAXISHTC
c                  RAXIS5/RAXISV [COARSE FINE HORIZONTAL VERTICAL]]
c                  [LIPS [HORIZONTAL VERTICAL]] ADSC, SBC1, NOIR,
C                  TAURUS, MARCCD, MAR, LMB, SMALLMAR, MARFLAT,
C                  [DIP2000/2030/2040 [ADCT VERT]] [ESRF [CCD]]
C                  [FUJI [BAS2000 BA100]] [MD [COARSE FINE]]
C                  JUPITER SATURN MERCURY A200 RIPI RIEI OXFM OXFO 
C                  [BRUKER [PROTEUM SMART]] CBF PILATUS
C                  REVERSEPHI INVERTX OMEGAF
C                  [ROTATION H/V A/C] [ORIGIN LL/LR/UL/UR]
C                  [FAST H/V] [TYPE CIRC/RECT] 
C                  ROWREAD ROTNSPEED
C     
      ELSE IF ((KEY.EQ.'MACH').OR.(KEY.EQ.'SCAN').OR.
     +       (KEY.EQ.'DETE')) THEN
c        write(*, *) 'SETTING DETECTOR TYPE'

C
C---- Only set IISCN TRUE if an actual detector type is specified.
C     AGWL 5/1/17
C
CAL        IISCN = .TRUE.
        IF (NTOK.GT.1) SVSCN = LINE(IBEG(2):IEND(NTOK))
        OMEGAREV = 0.0
        IPART = 0
        ISCAN = 1
        ICOUNT = 1
 208    ICOUNT = ICOUNT + 1
        KEY8 = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     *************
        CALL CCPUPC(KEY8)
C     *************
        SUBKEY = KEY8(1:4)
        IF (SUBKEY.EQ.'ERRO') THEN
          WRITE(IOUT,FMT=6611)
          IF (ONLINE) WRITE(ITOUT,FMT=6611)
          IF (BRIEF) WRITE(IBRIEF,FMT=6611)
 6611     FORMAT(1X,'**** ERROR ****',/,1X,'This keyword is now ',
     +         'redundant.',/,1X,'You are ',
     +         'strongly advised NOT to do this unless you',
     +         ' know what you are doing !!')
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          EFAC = VALUE(ICOUNT)
C     
        ELSE IF (SUBKEY.EQ.'OFFL') THEN
          NEWPREF = .TRUE.
C     
        ELSE  IF (SUBKEY.EQ.'RAXI') THEN
          IISCN = .TRUE.
          MACHINE = 'RAXI'
          COMPRESSION = 3
          FINE = .TRUE.
          IF (IEXTEN.EQ.0) ODEXT = 'osc'
C     
C---- Check for Raxis4 or RaxisIV
C     
          IF ((KEY8(1:6).EQ.'RAXIS4').OR.(KEY8(1:7).EQ.'RAXISIV'))
     +         THEN
            MODEL = 'RAXISIV'
            USEHDR = .TRUE.
C     
C---- Check for RaxisHTC
C     
          ELSE IF (KEY8(1:8).EQ.'RAXISHTC')THEN
            MODEL = 'RAXISHTC'
            USEHDR = .TRUE.
C     
C---- Check for Raxis5 or RaxisV
C     
          ELSE IF ((KEY8(1:6).EQ.'RAXIS5').OR.(KEY8(1:6).EQ.'RAXISV'))
     +           THEN
            MODEL = 'RAXISV'
            USEHDR = .TRUE.
          ELSE
            MODEL = 'RAXIS  '
          END IF
          CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
C     
C---- Test for keywords COARSE, FINE, HORIZONTAL (rotation axis)
C     
          IF (ICOUNT.LT.NTOK) THEN
            SCAN = LINE(IBEG(ICOUNT+1):IEND(ICOUNT+1))
C     
C---- Convert to upper case
C     
C     ***********
            CALL CCPUPC(SCAN)
C     ***********
            IF (SCAN.EQ.'COAR') THEN
              ICOUNT = ICOUNT + 1
              FINE = .FALSE.
              RAST = 0.210
              NREC = 950
              IYLEN = 950
            ELSE IF (SCAN.EQ.'FINE') THEN
              ICOUNT = ICOUNT + 1
              FINE = .TRUE.
            ELSE IF (SCAN.EQ.'HORI') THEN
              ICOUNT = ICOUNT + 1
              IF ((MODEL.EQ.'RAXISIV').OR.(MODEL.EQ.'RAXISV')
     $             .OR.(MODEL.EQ.'RAXISHTC')) THEN
                OMEGAFD = 0.0
                OMEGAF = OMEGAFD * DTOR
              ELSE
                OMEGAFD = 180.0
                OMEGAF = OMEGAFD * DTOR
              END IF
              IF (ICOUNT.LT.NTOK) THEN
                IF(LINE(IBEG(ICOUNT+1):IEND(ICOUNT+1)).EQ.'COAR')
     +               THEN
                  ICOUNT = ICOUNT + 1
                  FINE = .FALSE.
                  IF (IPIX.EQ.0) RAST = 0.210
                  NREC = 950
                  IYLEN = 950
                END IF
              END IF
            ELSE IF (SCAN.EQ.'VERT') THEN
              ICOUNT = ICOUNT + 1
C     
C---- This is the default, so do not need to set anything
C     
            END IF
          END IF
C     
C---- ESRF Large format IP
C     
        ELSE IF (SUBKEY.EQ.'LIPS') THEN
C     
C---- Test for keywords HORIZONTAL or VERTICAL (orientation of plate)
C     
          IISCN = .TRUE.
          SCAN = 'VERT'
          IF (ICOUNT.LT.NTOK) THEN
            SCAN = LINE(IBEG(ICOUNT+1):IEND(ICOUNT+1))
C     
C---- Convert to upper case
C     
C     ***********
            CALL CCPUPC(SCAN)
C     ***********
            IF (SCAN.EQ.'HORI') THEN
              ICOUNT = ICOUNT + 1
C     
            ELSE IF (SCAN.EQ.'VERT') THEN
              ICOUNT = ICOUNT + 1
C     
            ELSE
              WRITE(IOUT,FMT=7610)
 7610         FORMAT(1X,'*** ERROR ***',/,1X,' Can only give',
     +             ' HORIZONTAL or VERTICAL as subkeywords')
              IF (ONLINE) THEN
                WRITE(IOUT,FMT=7610)
                IF (COMREAD) THEN
                  COMREAD = .FALSE.
                  ITIN = ITINS
                  CLOSE (UNIT=ICOMM)
                END IF
                BADLINE = .TRUE.
                GOTO 50
              END IF
            END IF
          END IF
          IF (XSCAN.EQ.0) XSCAN = 20000
          IF (YSCAN.EQ.0) YSCAN = 40000
          XMAXIP = 20000
          YMAXIP = 40000
          MACHINE = 'LIPS'
          COMPRESSION = 1
          INVERTX = .TRUE.
          SPIRAL = .FALSE.
          CIRCULAR = .FALSE.
          ORTHOG = .TRUE.
          NEWPREF = .TRUE.
          ISCAN = 1
          IF (SCAN.EQ.'VERT') THEN
            OMEGAFD = 270.0
            omegaf = omegafd * dtor
          ELSE
            OMEGAFD = 180.0
            omegaf = omegafd * dtor
          END IF
          IF(.NOT.IINULL)NULLPIX = 0
          IF(IIDIVIDE.EQ.0)IDIVIDE = 0
          RMINIP = 300
          RMAXIP = SQRT(XMAXIP**2 + YMAXIP**2)
          IF (ILIMIT.EQ.0) THEN
            XLIMIT = 0.01*(XMAXIP/2)
            LIMIT = NINT(100*XLIMIT)
          END IF
          IF (ICUT.EQ.0) CUTOFF = 65500
          IF (IPRCUT.EQ.0) PRCUTOFF = 65500
          IF (IGAIN.EQ.0)THEN
            IGAIN = 1
            GAIN = 1.0
          ENDIF
C     
C---- ADSC CCD Detectors
C     
        ELSE IF (SUBKEY.EQ.'ADSC') THEN
          IISCN = .TRUE.
          IF(KEY8.EQ.'ADSC-PAD')THEN
            MACHINE = 'APAD'
	    MINICBF = .TRUE.
            COMPRESSION = 6
C     ..
C     .. 11.03.2015 - only model is HF-4M
            MODEL = 'HF-4M'
            IF (IEXTEN.EQ.0) ODEXT = 'cbf'
            CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
          ELSE
            MACHINE = 'ADSC'
            COMPRESSION = 1
            MODEL = 'QUAD'
            IF (IEXTEN.EQ.0) ODEXT = 'img'
            CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
          ENDIF
C     
C---- SBC1 CCD Detector (3x3 Ed Westbrook)
C     
        ELSE IF (SUBKEY.EQ.'SBC1') THEN
          IISCN = .TRUE.
          MACHINE = 'SBC1'
          COMPRESSION = 1
          IF (IEXTEN.EQ.0) ODEXT = 'img'
          CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
C     
C---- NOIR lens coupled CCD Detector (Ed Westbrook)
C     
        ELSE IF (SUBKEY.EQ.'NOIR') THEN
          IISCN = .TRUE.
          MACHINE = 'NOIR'
          COMPRESSION = 1
          IF (IEXTEN.EQ.0) ODEXT = 'img'
          CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
C     
C---- TAURUS CMOS Detector (Ed Westbrook)
C     
        ELSE IF (SUBKEY.EQ.'TAUR') THEN
          IISCN = .TRUE.
          MACHINE = 'TAUR'
          COMPRESSION = 1
          IF (IEXTEN.EQ.0) ODEXT = 'img'
          CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
C     
C---- Mar CCD Detectors
C     
        ELSE IF (SUBKEY.EQ.'MARC') THEN
          IISCN = .TRUE.
          MACHINE = 'MARC'
          COMPRESSION = 1
          MODEL = ' '
          IF (IEXTEN.EQ.0) ODEXT = 'mccd'
          CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
C     
C---- Mar (default is big scanner)
C    
        ELSE IF (SUBKEY.EQ.'MAR') THEN
c           write(*, *) 'SETTING MAR DETECTOR'
          IISCN = .TRUE.
          MACHINE = 'MAR '
          COMPRESSION = 1
          IF (IPIX.EQ.0) RAST = 0.15
          NHEAD = 1
          NTAIL = 0
          NREC = 2000
          IYLEN = 2000
          OMEGAFD = 90.0
          omegaf = omegafd * dtor
C     
C---- LMB (new big scanner)
C     
        ELSE IF (SUBKEY.EQ.'LMB') THEN
          IISCN = .TRUE.
          MACHINE = 'LMB'
          COMPRESSION = 1
          USEHDR = .FALSE.
          SPIRAL = .FALSE.
          CIRCULAR = .FALSE.
          ORTHOG = .TRUE.
          INVERTX = .TRUE.
          IF (IPIX.EQ.0) RAST = 0.16667
          NHEAD = 1
          NTAIL = 0
          NREC = 3000
          IYLEN = 3000
          OMEGAFD = 90.0
          omegaf = omegafd * dtor
          PACK = .TRUE.
          IF (IEXTEN.EQ.0) ODEXT = 'pck'
C     
C---- SmallMar 
C     
        ELSE IF (SUBKEY.EQ.'SMAL') THEN
          IISCN = .TRUE.
          MACHINE = 'MAR '
          COMPRESSION = 1
          IF (IPIX.EQ.0) RAST = 0.15
          NHEAD = 1
          NTAIL = 0
          NREC = 1200
          IYLEN = 1200
          OMEGAFD = 90.0
          omegaf = omegafd * dtor
C     
C---- Mar Flat panel
C     
        ELSE IF (SUBKEY.EQ.'MARF') THEN
          IISCN = .TRUE.
          MACHINE = 'MARF'
          COMPRESSION = 2
          MODEL = 'M555'
          OMEGAFD = 0.0
          omegaf = omegafd * dtor
          CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
C     
C---- DIP2000
C     
        ELSE IF (SUBKEY.EQ.'DIP2') THEN
          IISCN = .TRUE.
          MACHINE = 'DIP2'
          COMPRESSION = 5
c          MODEL = '16BIT'
          KEY8 = LINE(IBEG(ICOUNT):IEND(ICOUNT))
          CALL CCPUPC(KEY8)
          IF(KEY8(1:7).EQ.'DIP2000')
     $         IFILESIZE = 12501024
          IF(KEY8(1:7).EQ.'DIP2030')
     $         IFILESIZE = 18001024
          IF(KEY8(1:7).EQ.'DIP2040')
     $         IFILESIZE = 32001024
          CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
C     
C---- Check if keyword for rotation axis vertical given
C     
C     HRP from Atsushi Nakagawa             IF (ICOUNT.LT.NTOK) THEN
C     HRP from Atsushi Nakagawa                KEY2 = LINE(IBEG(ICOUNT+1
C     ):IEND(ICOUNT+1))
 7365     IF (ICOUNT.LT.NTOK) THEN
            ICOUNT=ICOUNT+1
            KEY2 = LINE(IBEG(ICOUNT):IEND(ICOUNT))
            CALL CCPUPC(KEY2)
C     HRP from Atsushi Nakagawa                IF (KEY2.EQ.'VERT') THEN
C     HRP from Atsushi Nakagawa                  OMEGAFD = 90.0
C     HRP from Atsushi Nakagawa                  ICOUNT = ICOUNT + 1
            IF (KEY2.EQ.'ADCT') THEN
              ICOUNT=ICOUNT+1
              KEY8 = LINE(IBEG(ICOUNT):IEND(ICOUNT))
              CALL CCPUPC(KEY8)
              IF(KEY8.NE.'12BIT'.AND.KEY8.NE.'16BITS'
     +             .AND.KEY8.NE.'16BITD'.AND.KEY8.NE.'16BIT') THEN
                WRITE(IOUT,FMT=7361)
                IF (ONLINE) WRITE(ITOUT,FMT=7361)
 7361           FORMAT(/,1X,'***** WARNING *****',
     +               /,1X,'***** WARNING *****',
     +               /,1X,'***** WARNING *****',
     +               /,1X,'***** WARNING *****',
     +               /,1X,'ADC has to be either 16BIT, ',
     +               /,'16BITD, 16BITS or 12BIT, or ',
     +               'taken from the tailer record')
              ELSE
                MODEL=KEY8
                WRITE(IOUT,FMT=7364) MODEL
                SETADC = .TRUE.
                IF (ONLINE) WRITE(ITOUT,FMT=7364) MODEL
 7364           FORMAT(/'*** ADCTYPE has been set to ',A)
              END IF
            ELSE IF (KEY2.EQ.'VERT') THEN
              OMEGAFD = 90.0
              OMEGAF = OMEGAFD * DTOR
            ELSE
              WRITE(IOUT,FMT=7363) KEY2
              IF (ONLINE) WRITE(ITOUT,FMT=7363) KEY2
 7363         FORMAT(1X,' UNKNOWN SUB-KEYWORD :',A)
            END IF
          ELSE
            WRITE(IOUT,FMT=7362)
            IF (ONLINE) WRITE(ITOUT,FMT=7362)
          END IF
 7362     FORMAT(/,1X,'***** WARNING *****',
     +         /,1X,'***** WARNING *****',
     +         /,1X,'***** WARNING *****',
     +         /,1X,'***** WARNING *****',
     +         /,1X,'It is assumed that the scanner has been ',
     +         'modified to have a horizontal rotation axis.',/,
     +         1X,'If this is not the case, give the extra ',
     +         'keyword VERTICAL',/,1X,'SCANNER DIP2000 VERTICAL',
     +         ', SCANNER DIP2030 VERTICAL',
     +         ' or SCANNER DIP2030 VERTICAL')
C     
C---- ESRF CCD scanner
C     
        ELSE IF (SUBKEY.EQ.'ESRF') THEN
          IISCN = .TRUE.
          Hbeamline='ESRF'
C     
C---- Get scanner type
C     
          IF (NTOK.GT.2) THEN
            ICOUNT = ICOUNT + 1
            MODEL = LINE(IBEG(3):IEND(3))
C     
C---- Convert to upper case
C     
C     ***********
            CALL CCPUPC(MODEL)
C     ***********
            IF (MODEL(1:3).EQ.'CCD') THEN
              CONTINUE
            END IF
          END IF
C     
C---- Following spatial distortion corection and non-uniformity
C     correction the
C     images are written in a pseudo Mar format, with a header record
C     that ONLY
C     contains the image size.
C     Note that the pixel size depends on the crystal to detector
C     distance,
C     so there is  no default
C     
          MACHINE = 'CCD2'
          COMPRESSION = 1
          USEHDR = .FALSE.
          SPIRAL = .FALSE.
          CIRCULAR = .FALSE.
          ORTHOG = .TRUE.
          INVERTX = .TRUE.
          NHEAD = 1
          NTAIL = 0
          NREC = 1200
          IYLEN = 1200
          OMEGAFD = 270.0
          omegaf = omegafd * dtor
          RMAXSP = 55.0
          IF (IEXTEN.EQ.0) ODEXT = 'cor'
C     
C---- Fuji offline scanner
C     
        ELSE IF (SUBKEY.EQ.'FUJI') THEN
          IISCN = .TRUE.
          HDRSIZE = .FALSE.
          USEHDR = .FALSE.
          INVERTX = .FALSE.
          ORTHOG = .TRUE.
          MACHINE = 'FUJI'
          COMPRESSION = 1
          NEWPREF = .TRUE.
          NHEAD = 0
          NREC = 2560
          IYLEN = 2048
          IF (IPIX.EQ.0) RAST = 0.10
          SPIRAL = .FALSE.
          CIRCULAR = .FALSE.
          ORTHOG = .TRUE.
          LOGA = 4.0
          LOGB = 1023.0
          IF (IEXTEN.EQ.0) ODEXT = 'fuj'
C     
C---- Test for scanner type
C     
          IF (NTOK.GT.2) THEN
            MODEL = LINE(IBEG(3):IEND(3))
            ICOUNT = ICOUNT + 1
C     
C---- Convert to upper case
C     
C     ***********
            CALL CCPUPC(MODEL)
C     ***********
            IF (MODEL(1:7).EQ.'BAS2000') THEN
              LOGA = 4.0
              LOGB = 1023.0
            ELSE IF (MODEL(1:5).EQ.'BA100') THEN
              LOGA = 1.0
              LOGB = 255.0
            ELSE
              WRITE(IOUT,FMT=6702) MODEL
              IF (ONLINE) WRITE(ITOUT,FMT=6702) MODEL
 6702         FORMAT(1X,' ***ERROR, Fuji scanner type ',A,
     +             ' not known',/,1X,'must be BAS2000 or BA100')
              STOP
            END IF
          END IF
C     
C---- Molecular Dynamics (MD)
C     
        ELSE IF (SUBKEY.EQ.'MD') THEN
          IISCN = .TRUE.
          MACHINE = 'MD'
          COMPRESSION = 4
          USEHDR = .FALSE.
          HDRSIZE = .FALSE.
          NEWPREF = .TRUE.
          SPIRAL = .FALSE.
          CIRCULAR = .FALSE.
          ORTHOG = .TRUE.
          NHEAD = 0
          FINE = .TRUE.
          IF (IPIX.EQ.0) RAST = 0.08585
C     
C---- Test for keyword COARSE
C     
          IF (NTOK.GT.2) THEN
            ICOUNT = ICOUNT + 1
            SCAN = LINE(IBEG(3):IEND(3))
C     
C---- Convert to upper case
C     
C     ***********
            CALL CCPUPC(SCAN)
C     ***********
            IF (SCAN.EQ.'COAR') THEN
              FINE = .FALSE.
              IF (IPIX.EQ.0) RAST = 0.172
            ELSE IF (SCAN.EQ.'FINE') THEN
              FINE = .TRUE.
            ELSE
              WRITE(IOUT,FMT=6700) SCAN
              IF (ONLINE) WRITE(ITOUT,FMT=6700) SCAN
              STOP
            END IF
          END IF
C     
C---- Jupiter CCD detector (also for Saturn and Mercury). Should rely on
C     active mask for edges of tiles.
C     
        ELSE IF ((SUBKEY.EQ.'JUPI').OR.(SUBKEY.EQ.'SATU').OR.
     $         (SUBKEY.EQ.'MERC').OR.(SUBKEY.EQ.'A200').OR.
     $         (SUBKEY.EQ.'RIPI').OR.(SUBKEY.EQ.'RIEI')) THEN
          IISCN = .TRUE.
          MACHINE = SUBKEY
          COMPRESSION = 3
          IF(SUBKEY.EQ.'JUPI')THEN
            MODEL = 'JUPITER'
          ELSEIF (SUBKEY.EQ.'SATU')THEN
            MODEL = 'SATURN'
          ELSEIF (SUBKEY.EQ.'MERC')THEN
            MODEL = 'MERCURY'
          ELSEIF (SUBKEY.EQ.'A200')THEN
            MODEL = 'A200'
          ELSEIF (SUBKEY.EQ.'RIPI')THEN
            MODEL = 'RIPI'
            COMPRESSION = 6
            GAIN = 0.80
          ELSEIF (SUBKEY.EQ.'RIEI')THEN
            MODEL = 'RIEI'
            COMPRESSION = 6
            GAIN = 0.80
          ENDIF
C
C---- Print for program flow
C
          CALLED = 'SETDETECT'
          COMMENT = ' to assign detector parameters for '//MODEL
          IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
          CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
C     
C---- Small Oxford CCD detector, images currently tested have been
C     converted to Mar IP format 
C     
        ELSE IF ((SUBKEY.EQ.'OXFM').OR.(SUBKEY.EQ.'OXFO')) THEN
C
          IISCN = .TRUE.
C     Oxford in Mar format probably obsolete, and certainly not autodetected!
          IF (SUBKEY.EQ.'OXFM')THEN
            MACHINE = 'OXFM'
            COMPRESSION = 1
            MODEL = 'OXM345'
            NREC = 512
            IYLEN = 512
            INVERTX = .TRUE.
            OMEGAFD = 0.0
            IF (IEXTEN.EQ.0) ODEXT = 'mar512'
            RMINIP = 300
            IF (XSCAN.EQ.0) XSCAN = 3072
            IF (YSCAN.EQ.0) YSCAN = 3072
            XMAXIP = 3072
            YMAXIP = 3072
            RMAXIP = SQRT(XMAXIP**2 + YMAXIP**2)
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
          ELSEIF (SUBKEY.EQ.'OXFO')THEN
C
C---  lots of stuff is in the image header for Oxford Sapphire images -
C     e.g. NREC, IYLEN, but we need to put something here so they're non-zero.
C
            MACHINE = 'OXFO'
            COMPRESSION = 4
            MODEL = 'SAPPHIRE'
            CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
          ENDIF
C
C---- Oxford diffraction packed files from Crysalis detector
C
        ELSE IF (SUBKEY.EQ.'OXPK') THEN
C
          IISCN = .TRUE.
          MACHINE = 'OXPK'
          CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
C     
C---- Bruker CCD detectors
C
        ELSE IF (SUBKEY.EQ.'BRUK') THEN
          IISCN = .TRUE.
          MACHINE = 'BRUK'
          COMPRESSION = 1
          IF(NTOK.GE.3)THEN
            SUBKEY = LINE(IBEG(3):IEND(3))
            CALL CCPUPC(SUBKEY)
            IF (SUBKEY.EQ.'PROT')THEN
              MODEL = 'PROTEUM'
            ELSE
              MODEL = 'SMART'
            ENDIF
          ENDIF
          CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
C
C---- CBF (Crystallographic Binary Format, binary version of imgCIF)
C     This code uses the CBFlib API called via OPENODS. Reading the 
C     header is all part of the game here.
C     
        ELSE IF (SUBKEY.EQ.'CBF ') THEN
          IISCN = .TRUE.
          MACHINE = 'CBF '
          COMPRESSION = 1
          MODEL = ''
C
C---- HDF5 file - only used for Eiger as of 02.02.2018
C     
        ELSE IF (SUBKEY.EQ.'HDF5') THEN
          MACHINE = 'HDF5'
          HDF5FILE = .TRUE.
          COMPRESSION = 6
          MODEL = ''
CAL          IF(.NOT.HEADONLY)THEN
CAL            WRITE(IOUT,FMT=7700)
CAL            IF(ONLINE)WRITE(ITOUT,FMT=7700)
CAL 7700       FORMAT('warning, warning, warning. This detector option',
CAL     $           ' has not been fully enabled yet.',/,/,
CAL     $           'You should be prepared for unexpected results.')
CAL          ENDIF
C     
C---- almost nothing should be set here for a CBF file - it should all
C     be found 
C     when reading the file itself, in OPENODS. We need to check in
C     OPENODS 
C     if either MACHINE.EQ.CBF or MACHINE.EQ.MAR which type the image is
C     , 
C     because we shouldn't have to tell the program if the image is CBF
C     type.
        ELSE IF (SUBKEY.EQ.'PILA') THEN
          IISCN = .TRUE.
	  MINICBF = .TRUE.
          MACHINE = 'PILA'
          COMPRESSION = 6
C work out in section.f          MODEL = 'PILATUS '
          IF (IEXTEN.EQ.0) ODEXT = 'cbf'
          CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
        ELSE IF (SUBKEY.EQ.'EIGE') THEN
	  MINICBF = .TRUE.
          MACHINE = 'EIGE'
          COMPRESSION = 6
C work out in section.f          MODEL = 'Eiger'
          IF (IEXTEN.EQ.0) ODEXT = 'cbf'
          CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
C
C---- PHOTON the Bruker PHOTON II or PHOTON III detector
C
        ELSE IF (SUBKEY.EQ.'PHOT') THEN
          IISCN = .TRUE.
          MACHINE = 'PHOT'
          FULLCBF = .TRUE.
c          CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
C
C---- TIME the TIMEPIX detector
C
        ELSE IF (SUBKEY.EQ.'TIME') THEN
          IISCN = .TRUE.
          MACHINE = 'TIME'
	  MINICBF = .TRUE.
          CALL SETDETECT(DMODE,INSIZE,IFILESIZE,IERR)
C
C---- ROWREAD... time (seconds) to read a row of data in rolling shutter mode
C
        ELSE IF (SUBKEY.EQ.'ROWR') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          ROWREADT  = VALUE(ICOUNT)
c          write(iout,7840) ROWREADT
c          if (online) write(itout,7840) ROWREADT
 7840     format(1x,'Row read time set to',F8.5,' seconds')
C
C---- ROTNSPEED... spindle rotation rate (degrees per second) required 
C                  for rolling shutter mode corrections 
C
        ELSE IF (SUBKEY.EQ.'ROTN') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          ROTSPEED  = VALUE(ICOUNT)
c          WRITE(IOUT,7842) ROTSPEED
c          if (online) write(itout,7842) ROTSPEED
 7842     format(1x,'Rotation speed set to',F8.4,' degrees per second')
C
C     
C     
C---- REVERSEPHI...allow normal phi direction to be reversed.
C     
        ELSE IF (SUBKEY.EQ.'REVE') THEN
          OMEGAREV = 180.0
C     
C---- INVERTX detector
C     comment out, INVERTX now has a sub-keyword to specify
C     if it is TRUE or FALSE, see below. AGWL 11/11/2016   
CAL        ELSE IF (SUBKEY.EQ.'INVE') THEN
CAL          IcounT = IcounT + 1
CAL          INVERTX = .TRUE.

C     
C---- Allow OMEGAF to be defined explicitly
C     
        ELSE IF (SUBKEY.EQ.'OMEG') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          OMEGAFD  = VALUE(ICOUNT)
          OMEGAF = OMEGAFD * DTOR
C
C---- Setting IOMEGAF=1 prevents the value given here being replaced
C     by the default value for the detector.
C
          IOMEGAF = 1
C     
C---- Allow INVERTX to be defined explicitly, must be followed by "T"
C     for TRUE, any other string will be interpreted as FALSE.
C     
        ELSE IF (SUBKEY.EQ.'INVE') THEN
          ICOUNT = ICOUNT + 1
          KEY2 = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C         ***********
          CALL CCPUPC(KEY2)
C         ***********
          IF ((KEY2(1:1).NE.'T').AND.(KEY2(1:1).NE.'F')) THEN
            WRITE(IOUT,FMT=7850) KEY2
            IF (ONLINE) WRITE(ITOUT,FMT=7850) KEY2
 7850       FORMAT(/,1X,'***** ERROR ****',/,1X,'INVERT must be',
     +      ' followed by T or F but not by:',A)
          ELSE
            INVERTX = (KEY2(1:1).EQ.'T')
            IINVERTX = 1
          END IF
C     
C---- Allow definition of new detector types using a generic description.
C     This requires the following information:
C     ROTATION axis:  Horizontal or Vertical, Clockwise or Anticlockwise
C     (Sense of rotation is when viewed from above for
C     a vertical rotation axis, or when viewed down an
C     axis running from right to left (cameraman's view)
C     for a horizontal axis
C     Cameraman's view is from behind the detector looking
C     towards the source.
C     
C     ORIGIN (first pixel): LR (lower right),LL (lower left) ,
C     UR (upper right),UL (upper left) (Cameraman's view)
C     
C     FAST direction in image: Vertical or Horizontal
C     
C     CIRCULAR or RECTANGULAR specifies the shape of the active area.
C     The physical limits of the size of the active area should be
C     given on the LIMITS keyword
C     
C     The size of the image and the number of header records should be
C     given
C     on the SIZE keyword.
C     
        ELSE IF (SUBKEY.EQ.'ROTA') THEN
          IPART = IPART + 1
          MACHINE = 'UNK'
          COMPRESSION = 1
C     
C---- Horizontal or Vertical
C     
          ICOUNT = ICOUNT + 1
          KEY2 = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     ***********
          CALL CCPUPC(KEY2)
C     ***********
          IF (KEY2(1:1).EQ.'H') THEN
            ROTH = .TRUE.
            ROTV = .FALSE.
          ELSE IF (KEY2(1:1).EQ.'V') THEN
            ROTV = .TRUE.
            ROTH = .FALSE.
          ELSE
            WRITE(IOUT,FMT=6580)
            IF (ONLINE) WRITE(ITOUT,FMT=6580)
 6580       FORMAT(1X,'*** ERROR *** ROTATION must be Horizontal',
     +           ' or Vertical followed by Clockwise or ',
     $           'Anticlockwise')
            IF (ONLINE) THEN
              IF (COMREAD) THEN
                COMREAD = .FALSE.
                ITIN = ITINS
                CLOSE (UNIT=ICOMM)
              END IF
              BADLINE = .TRUE.
              GOTO 50
            END IF
            STOP
          END IF
C     
C---- Clockwise or Anticlockwise
C     
          ICOUNT = ICOUNT + 1
          KEY2 = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     ***********
          CALL CCPUPC(KEY2)
C     ***********
          IF (KEY2(1:1).EQ.'A') THEN
            ROTANTI = .TRUE.
            ROTCLOCK = .FALSE.
          ELSE IF (KEY2(1:1).EQ.'C') THEN
            ROTCLOCK = .TRUE.
            ROTANTI = .FALSE.
          ELSE
            WRITE(IOUT,FMT=6580)
            IF (ONLINE) WRITE(ITOUT,FMT=6580)
            IF (ONLINE) THEN
              IF (COMREAD) THEN
                COMREAD = .FALSE.
                ITIN = ITINS
                CLOSE (UNIT=ICOMM)
              END IF
              BADLINE = .TRUE.
              GOTO 50
            END IF
            STOP
          END IF
C     
C---- Origin LL,LR,UL,UR
C     
        ELSE IF (SUBKEY.EQ.'ORIG') THEN
          IPART = IPART + 1
          ICOUNT = ICOUNT + 1
          KEY2 = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     ***********
          CALL CCPUPC(KEY2)
C     ***********
          IF (KEY2(1:2).EQ.'LL') THEN
            ORGLL = .TRUE.
          ELSE IF (KEY2(1:2).EQ.'LR') THEN
            ORGLR = .TRUE.
          ELSE IF (KEY2(1:2).EQ.'UL') THEN
            ORGUL = .TRUE.
          ELSE IF (KEY2(1:2).EQ.'UR') THEN
            ORGUR = .TRUE.
          ELSE
            WRITE(IOUT,FMT=6582)
            IF (ONLINE) WRITE(ITOUT,FMT=6582)
 6582       FORMAT(1X,'*** ERROR *** ORIGIN must be given as one',
     +           ' of: LL (lower left),LR (lower right),UL ',
     +           ' (upper left), or UR (upper right)')
            IF (ONLINE) THEN
              IF (COMREAD) THEN
                COMREAD = .FALSE.
                ITIN = ITINS
                CLOSE (UNIT=ICOMM)
              END IF
              BADLINE = .TRUE.
              GOTO 50
            END IF
            STOP
          END IF
C     
C---- FAST   Horizontal or Vertical
C     
        ELSE IF (SUBKEY.EQ.'FAST') THEN
          IPART = IPART + 1
          ICOUNT = ICOUNT + 1
          FASTH = .FALSE.
          FASTV = .FALSE.
          STR = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     ***********
          CALL CCPUPC(STR)
C     ***********
C     
          IF (STR(1:1).EQ.'H') THEN
            FASTH = .TRUE.
          ELSE IF (STR(1:1).EQ.'V') THEN
            FASTV = .TRUE.
          ELSE
            WRITE (IOUT,FMT=6130) STR
            IF (ONLINE) WRITE (ITOUT,FMT=6130) STR
          END IF
C     AL             WRITE(6,*),'FASTH,FASTV',FASTH,FASTV
C     
C---- TYPE detector name (used in INTPXL to convert pixel values to
C     true counts)
C     
        ELSE IF (SUBKEY.EQ.'TYPE') THEN
          ICOUNT = ICOUNT + 1
          IF (ICOUNT.LE.NTOK) THEN
            KEY2 = LINE(IBEG(ICOUNT):IEND(ICOUNT))
          ELSE
            WRITE(IOUT,FMT=7360)
            IF (ONLINE) WRITE(ITOUT,FMT=7360)
 7360       FORMAT(1X,'***** ERROR *****',/,1X,
     +           'Must specify type of detector (MAR,DIP,RAXIS',
     +           ',RAXISIV,RAXISV,RAXISHTC,FUJI,MD,CCD1,CCD2,ADSC,',
     $           'SBC1,JUPITER,UNKNOWN')
            IF (ONLINE) THEN
              IF (COMREAD) THEN
                COMREAD = .FALSE.
                ITIN = ITINS
                CLOSE (UNIT=ICOMM)
              END IF
              BADLINE = .TRUE.
              GOTO 50
            END IF
            STOP
          END IF
C     ***********
          CALL CCPUPC(KEY2)
C     ***********
          MACHINE = KEY2
C     
C---- CIRCULAR detector
C     
        ELSE IF (SUBKEY.EQ.'CIRC') THEN
          IPART = IPART + 1
          ORTHOG = .FALSE.
          SPIRAL = .TRUE.
          CIRCULAR = .TRUE.
C     
C---- RECTANGULAR detector
C     
        ELSE IF (SUBKEY.EQ.'RECT') THEN
          IPART = IPART + 1
          ORTHOG = .TRUE.
          SPIRAL = .FALSE.
          CIRCULAR = .FALSE.
        ELSE
C     
C---- Not recognised
          IF (TRAPERR) INPERR = .TRUE.
          WRITE(IOUT,6709) SUBKEY
          IF (ONLINE) WRITE(ITOUT,6709) SUBKEY
 6709     FORMAT(1X,'*** ERROR ***',/,1X,'Detector "',A,'" not known'
     $         ,/,1X,'Must be MAR, MARCCD, MARFLAT, RAXIS, MD, ',
     +         'DIP2000, ADSC, LIPS, SBC1, JUPITER,',/,1X,
     $         ' SATURN, MERCURY, A200, RIPI, RIEI, OXFORD, BRUKER or',
     $         ' NOIR',/,1X,'If defining a novel detector, keywords '
     +         'ROTATION, FAST, ORIGIN AND TYPE must all be given.')
          IF (ONLINE) THEN
            IF (COMREAD) THEN
              COMREAD = .FALSE.
              ITIN = ITINS
              CLOSE (UNIT=ICOMM)
            END IF
            BADLINE = .TRUE.
            GOTO 50
          ELSE       
            STOP
          END IF
        END IF 
C     
C---- Read additional tokens (if any) 
C     
        IF (ICOUNT.LT.NTOK) GOTO 208
C     
        IF (IPART.GT.0) THEN
C     
C---- New scanner type, check all info given
C     
          IF (IPART.LT.4) THEN
            WRITE(IOUT,FMT=6711)         
            IF (ONLINE) WRITE(ITOUT,FMT=6711)
            IF (ONLINE) THEN
              IF (COMREAD) THEN
                COMREAD = .FALSE.
                ITIN = ITINS
                CLOSE (UNIT=ICOMM)
              END IF
              BADLINE = .TRUE.
              GOTO 50
            END IF
            STOP
          END IF
 6711     FORMAT(1X,'*** ERROR ***',/,1X,'For new scanner types ',
     +         ' the ROTATION, ORIGIN, FAST AND CIRCULAR/RECTANGULAR',
     +         /,1X,' keywords MUST ALL be given')
          INVERTX = .FALSE.
          IF (ORGLR.AND.FASTV) THEN
            OMEGAFD = 90.0
            omegaf = omegafd * dtor
            INVERTX = .TRUE.
          ELSE IF (ORGLR.AND.FASTH) THEN
            OMEGAFD = 0.0
            omegaf = omegafd * dtor
          ELSE IF (ORGLL.AND.FASTV) THEN
            OMEGAFD = 90.0
            omegaf = omegafd * dtor
          ELSE IF (ORGLL.AND.FASTH) THEN
            OMEGAFD = 180.0
            omegaf = omegafd * dtor
            INVERTX = .TRUE.
          ELSE IF (ORGUR.AND.FASTV) THEN
            OMEGAFD = 270.0
            omegaf = omegafd * dtor
          ELSE IF (ORGUR.AND.FASTH) THEN
            OMEGAFD = 0.0
            omegaf = omegafd * dtor
            INVERTX = .TRUE.
          ELSE IF (ORGUL.AND.FASTV) THEN
            OMEGAFD = 270.0
            omegaf = omegafd * dtor
            INVERTX = .TRUE.
          ELSE IF (ORGUL.AND.FASTH) THEN
            OMEGAFD = 180.0
            omegaf = omegafd * dtor
          END IF
          IF (ROTH) THEN
            IF (ROTANTI) THEN
              OMEGADD = 0.0
            ELSE
              OMEGADD = 180.0
            END IF
          ELSE IF (ROTV) THEN
            IF (ROTANTI) THEN
              OMEGADD = -90.0
            ELSE
              OMEGADD = 90.0
            END IF
          END IF
          OMEGAFD = OMEGAFD + OMEGADD
          omegaf = omegafd * dtor
          GENERICDET = .TRUE. 
        END IF
        OMEGAFD = OMEGAFD + OMEGAREV
        IF (OMEGAFD.GE.360.0) OMEGAFD = OMEGAFD - 360.0
        omegaf = omegafd * dtor

C     
        IF (DEBUG(52)) THEN
          WRITE(IOUT,FMT=7330) MACHINE,MODEL,INVERTX,SPIRAL,ORTHOG,
     +         CIRCULAR,OMEGAFD
          IF (ONLINE) WRITE(ITOUT,FMT=7330) MACHINE,MODEL,INVERTX,
     +         SPIRAL,ORTHOG,CIRCULAR,OMEGAFD
 7330     FORMAT(1X,'Machine type: ',A,' Model type: ',A,' INVERTX',
     +         l3,' SPIRAL',L3,' ORTHOG',L3,'  CIRCULAR',L3,
     +         ' OMEGA',F7.1)
        END IF
C     
C     
C---- Average spot profile
C     
      ELSE IF (KEY.EQ.'AVPR') THEN
C     
C---- Totally redundant keyword.
C     
        WRITE(IOUT,FMT=6612)
        IF (ONLINE) WRITE(ITOUT,FMT=6612)
        IF (BRIEF) WRITE(IBRIEF,FMT=6612)
 6610   FORMAT(1X,'**** ERROR ****',/,1X,'This keyword is now ',
     +       'redundant. If you really want to use this',/,1X,
     +       'option, you must include the subkeyord MODIFY ',
     +       '(before numerical input, if any)',/,1X,'You are ',
     +       'strongly advised NOT to do this unless you',
     +       ' know what you are doing !!')
 6612   FORMAT(1X,'**** WARNING ****',/,1X,'This keyword is now ',
     +       'redundant, it has been ignored')
        BADLINE = .TRUE.
        GOTO 50
C     
C---- Force program to go straight into filmplot
C     
      ELSE IF (KEY.EQ.'PLOT') THEN
        IF (.NOT.ONLINE) THEN
          NWRN = NWRN + 1
          WRITE (IOUT,FMT=6026)
 6026     FORMAT (//,1X,'*** Cannot call X-windows display ',
     $         ' from batch job')
          GO TO 50
        END IF
        IF (NTOK.EQ.1) THEN
C     
C     
          DO 280 I = 1,MAXPAX
            FILMPLOT(I) = .TRUE.
 280      CONTINUE
C     
C     
          NPROF = NPACK
        ELSE 
          ICOUNT = 1
 281      ICOUNT = ICOUNT + 1
          IF (ITYP(ICOUNT).EQ.2.AND.((ICOUNT.EQ.NTOK).OR.
     +         (ITYP(ICOUNT+1).NE.2))) THEN
C     
C     *******************************************
            CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            IDFILM = NINT(VALUE(2))
C     
C     
            DO 290 I = 1,MAXPAX
              FILMPLOT(I) = (IDFILM.EQ.IDPACK(I))
 290        CONTINUE
C     
C     
            NPROF = 1
            IDPROF(1) = IDFILM
          ELSE IF (ITYP(ICOUNT).EQ.2.AND.((ICOUNT.LT.NTOK).AND.
     +           (ITYP(ICOUNT+1).EQ.2))) THEN

C     
C     *******************************************
            CALL MKEYNM(2,2,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            ID1 = NINT(VALUE(2))
            ID2 = NINT(VALUE(3))
            NPROF = 0
C     
C     
            DO 310 I = 1,MAXPAX
              IF (ID1.EQ.IDPACK(I)) THEN
                DO 300 J = I,MAXPAX
                  FILMPLOT(J) = .TRUE.
                  NPROF = NPROF + 1
                  IDPROF(NPROF) = IDPACK(J)
                  IF (IDPACK(J).EQ.ID2) GO TO 320
 300            CONTINUE
              END IF
 310        CONTINUE
C     
          ELSE
            WRITE (IOUT,FMT=6130) STR
            IF (ONLINE) WRITE (ITOUT,FMT=6130) STR
          END IF
          IF (ICOUNT.LT.NTOK) GOTO 281
        END IF
 320    IF (NPROF.EQ.NPACK) THEN
          WRITE (IOUT,FMT=6028)
 6028     FORMAT (1X,'Use interactive graphics display for all',
     $         ' images')
          IF (ONLINE) WRITE (ITOUT,FMT=6028)
        ELSE
          WRITE (IOUT,FMT=6030) (IDPROF(I),I=1,NPROF)
 6030     FORMAT(1X,'Use interactive graphics display for following ',
     +         'images', / (1X,20I7,/))
          IF (ONLINE) WRITE (ITOUT,FMT=6030) (IDPROF(I),I=1,NPROF)
        END IF
C     
C---- assignment of packs to cassettes
C     
      ELSE IF (KEY.EQ.'CASS') THEN
        IF (IMGP) THEN
          WRITE(IOUT,6017) KEY
          IF (ONLINE) WRITE(ITOUT,6017) KEY
 6017     FORMAT(/,1X,'**** ERROR ***** ',A,' not appropriate for',
     +         ' image plates, ignored')
          BADLINE = .TRUE.
          GOTO 50
        END IF
C     
C     *******************************************
        CALL MKEYNM(2,2,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
        IPCKID = NINT(VALUE(2))
        NCASS = NINT(VALUE(3))
        WRITE (IOUT,FMT=6032) IPCKID,NCASS
 6032   FORMAT (1X,'PACK',I4,' will have camera constants for ',
     +       'CASSETTE',I3,/,1X,'Succeeding packs will be assi',
     +       'gned to CASSETTES in strict order, unless redefi',
     +       'ned by another "CASSETTE" card')
        IF (ONLINE) WRITE (ITOUT,FMT=6032) IPCKID,NCASS
        ASSIGN = .FALSE.
        DO 340 I = 1,NPACK
          IF (IDPACK(I).EQ.IPCKID) THEN
            ASSIGN = .TRUE.
            DO 330 J = I,NPACK
              ICASSET(J) = NCASS
              NCASS = NCASS + 1
              IF (NCASS.GT.8) NCASS = NCASS - 8
 330        CONTINUE
          END IF
 340    CONTINUE
        IF (.NOT.ASSIGN) THEN
          IF (ONLINE) WRITE (ITOUT,FMT=6034)
 6034     FORMAT (1X,'*** The pack ID you have given is ',
     +         'not one of those declared on the PACK card,',
     $         ' PLEASE REPEAT')
          WRITE (IOUT,FMT=6034)
        END IF
C     
C---- NOMEAS... program does not integrate images
C     
      ELSE IF (KEY.EQ.'NOME') THEN
        NOMEAS = .TRUE.
C     
C---- ALLOUT... all reflection except those in cusp and outside detector
C     or resolution limits are written to output MTZ file
C     
      ELSE IF (KEY.EQ.'ALLO') THEN
        ALLOUT = .TRUE.
C     
C---- mods
C     
      ELSE IF (KEY.EQ.'MOD') THEN
C     
C     
        DO 350 I = 2,NTOK
          SUBKEY = LINE(IBEG(I) :IEND(I))
C     
C     ************
          CALL CCPUPC(SUBKEY)
C     ************
C     
          IF (SUBKEY.EQ.'OLDB') THEN
            MODS(1) = .TRUE.
          ELSE
            WRITE (IOUT,FMT=6036)
 6036       FORMAT (/,1X,'***** Modification not recognised *****')
            IF (ONLINE) WRITE (ITOUT,FMT=6036)
          END IF
 350    CONTINUE
C     
C---- Debugs
C     
      ELSE IF (KEY.EQ.'DEBU') THEN
        IF (NTOK.EQ.1) THEN
          WRITE (IOUT,FMT=6038)
 6038     FORMAT (1X,'DEBUG output will be produced in all',
     $         ' subroutines')
          IF (ONLINE) WRITE (ITOUT,FMT=6038)
C     
C     
          DO 360 I = 1,200
            DEBUG(I) = .TRUE.
 360      CONTINUE
C     
C     
        ELSE
          DEBUGSTR = ' '
C     
C     
          I = 1
 370      I = I + 1
          IF (I.GT.NTOK) GOTO 372
          SUBKEY = LINE(IBEG(I) :IEND(I))
C     
C     ************
          CALL CCPUPC(SUBKEY)
C     ************
          IF ((I.EQ.2).AND.((SUBKEY.NE.'FLOW').AND.(SUBKEY.NE.'IMAG')))
     +                                             THEN
            DEBUGSTR = LINE(IBEG(I) :IEND(I))
          ELSE
            IF ((SUBKEY.NE.'FLOW').AND.(SUBKEY.NE.'IMAG')) THEN
              NCH = LENSTR(DEBUGSTR)
              DEBUGSTR = DEBUGSTR(1:NCH)//','//LINE(IBEG(I) :IEND(I))
            END IF
          END IF
C
          IF (SUBKEY.EQ.'NONE') THEN
C
C---- Turn of all DEBUG
C
            DO 362 I = 1,200
              DEBUG(I) = .FALSE.
 362        CONTINUE
     
          ELSE IF ((SUBKEY.EQ.'MAIN').OR.(SUBKEY.EQ.'MOSF')) THEN
            DEBUG(1) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(1)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'STAR') THEN
            DEBUG(2) = .TRUE.
          ELSE IF (SUBKEY.EQ.'FIND') THEN
C     
C---- Need to differentiate FINDPACK, FINDSPOTS
C     
            KEY6 = LINE(IBEG(I):IEND(I))
C     
C---- convert to upper case
C     
C     ***********
            CALL CCPUPC(KEY6)
C     ***********
            IF (KEY6.EQ.'FINDPA') THEN
              DEBUG(3) = .TRUE.
C     
C---- test for a number
C     
              IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
                I = I + 1
C     *********************************************
                CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
                NDEBUG(3)= NINT(VALUE(I))
              END IF
            ELSE IF (KEY6.EQ.'FINDSP') THEN
              DEBUG(63) = .TRUE.
            END IF
          ELSE IF (SUBKEY.EQ.'FIDU') THEN
            DEBUG(4) = .TRUE.
          ELSE IF (SUBKEY.EQ.'SEEK') THEN
            DEBUG(5) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(5)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'CENT') THEN
            DEBUG(6) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(6)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'FLMP') THEN
            DEBUG(7) = .TRUE.
          ELSE IF (SUBKEY.EQ.'REFP') THEN
            DEBUG(8) = .TRUE.
          ELSE IF (SUBKEY.EQ.'CURR') THEN
            DEBUG(9) = .TRUE.
          ELSE IF (SUBKEY.EQ.'SPOT') THEN
            DEBUG(10) = .TRUE.
          ELSE IF (SUBKEY.EQ.'RDIS') THEN
            DEBUG(11) = .TRUE.
          ELSE IF (SUBKEY.EQ.'RVDI') THEN
            DEBUG(12) = .TRUE.
          ELSE IF (SUBKEY.EQ.'NEXT') THEN
            DEBUG(13) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(13)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'RMAX') THEN
            DEBUG(14) = .TRUE.
          ELSE IF (SUBKEY.EQ.'CHKR') THEN
            DEBUG(15) = .TRUE.
          ELSE IF (SUBKEY.EQ.'GENS') THEN
            DEBUG(16) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(16)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'MEAS') THEN
            DEBUG(17) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(17)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'PROC') THEN
            DEBUG(18) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(18)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'WRGE') THEN
            DEBUG(19) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(19)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'RASP') THEN
            DEBUG(21) = .TRUE.
          ELSE IF (SUBKEY.EQ.'MEAN') THEN
            DEBUG(22) = .TRUE.
          ELSE IF (SUBKEY.EQ.'PIXL') THEN
            DEBUG(23) = .TRUE.
          ELSE IF (SUBKEY.EQ.'BGRE') THEN
            DEBUG(24) = .TRUE.
          ELSE IF (SUBKEY.EQ.'PSTA') THEN
            DEBUG(25) = .TRUE.
          ELSE IF (SUBKEY.EQ.'PRDI') THEN
            DEBUG(26) = .TRUE.
          ELSE IF (SUBKEY.EQ.'PXYC') THEN
            DEBUG(27) = .TRUE.
          ELSE IF (SUBKEY.EQ.'PWRG') THEN
            DEBUG(28) = .TRUE.
          ELSE IF (SUBKEY.EQ.'CONV') THEN
            DEBUG(29) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(29)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'REEK') THEN
            DEBUG(29) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(29)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'AUTO') THEN
            DEBUG(30) = .TRUE.
          ELSE IF (SUBKEY.EQ.'COOR') THEN
            DEBUG(31) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(31)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'GENE') THEN
            DEBUG(32) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(32)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'EVAL') THEN
            DEBUG(33) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(33)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'INTE') THEN
C     
C---- Need to differentiate INTEG,INTEG2,INTEG3 and INTEG4
C     
            KEY6 = LINE(IBEG(I) :IEND(I))
C     
C---- convert to upper case
C     
C     ***********
            CALL CCPUPC(KEY6)
C     ***********
            IF (KEY6.EQ.'INTEG2') THEN
              DEBUG(34) = .TRUE.
            ELSE IF (KEY6.EQ.'INTEG3') THEN
              DEBUG(46) = .TRUE.
            ELSE IF (KEY6.EQ.'INTEG4') THEN
              DEBUG(89) = .TRUE.
            ELSE
              DEBUG(43) = .TRUE.
            END IF
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              IF (KEY6.EQ.'INTEG2') THEN
                NDEBUG(34)= NINT(VALUE(I))
              ELSE IF (KEY6.EQ.'INTEG3') THEN
                NDEBUG(46)= NINT(VALUE(I))
              ELSE IF (KEY6.EQ.'INTEG4') THEN
                NDEBUG(89)= NINT(VALUE(I))
              ELSE
                NDEBUG(43)= NINT(VALUE(I))
              END IF
            END IF
          ELSE IF (SUBKEY.EQ.'VARP') THEN
            DEBUG(35) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********** **********************************
C     
              NDEBUG(35)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'GETP') THEN
            DEBUG(36) = .TRUE.
          ELSE IF (SUBKEY.EQ.'POST') THEN
            DEBUG(37) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********** **********************************
C     
              NDEBUG(37)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'DSTA') THEN
            DEBUG(38) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********** **********************************
C     
              NDEBUG(38)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'IDXR') THEN
            DEBUG(39) = .TRUE.
          ELSE IF (SUBKEY.EQ.'REFR') THEN
            DEBUG(40) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********** **********************************
C     
              NDEBUG(40)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'PRSE') THEN
            DEBUG(41) = .TRUE.
          ELSE IF (SUBKEY.EQ.'BEST') THEN
            DEBUG(42) = .TRUE.
          ELSE IF (SUBKEY.EQ.'MASK') THEN
            DEBUG(44) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(44)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'OPEN') THEN
            DEBUG(45) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(45)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'PKRI') THEN
            DEBUG(47) = .TRUE.
          ELSE IF (SUBKEY.EQ.'CHEC') THEN
            DEBUG(48) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(48)= NINT(VALUE(I))
            END IF
          ELSE IF (SUBKEY.EQ.'SPRO') THEN
            DEBUG(49) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(49)= NINT(VALUE(I))
            END IF
C     
C---- WRMTZ
C     
          ELSE IF (SUBKEY.EQ.'WRMT') THEN
            DEBUG(50) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(50)= NINT(VALUE(I))
            END IF
C     
C---- RSYMM or RSYMM2
C     
          ELSE IF (SUBKEY.EQ.'RSYM') THEN
            DEBUG(51) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(51)= NINT(VALUE(I))
            END IF
C     
C---- CONTROL
C     
          ELSE IF (SUBKEY.EQ.'CONT') THEN
            DEBUG(52) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(52)= NINT(VALUE(I))
            END IF
C     
C---- NEWRMS
C     
          ELSE IF (SUBKEY.EQ.'NEWR') THEN
            DEBUG(53) = .TRUE.
C     
C---- ROTATE
C     
          ELSE IF (SUBKEY.EQ.'ROTA') THEN
            DEBUG(54) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(54)= NINT(VALUE(I))
            END IF
C     
C---- MERGHKL
C     
          ELSE IF (SUBKEY.EQ.'MERG') THEN
            DEBUG(55) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(55)= NINT(VALUE(I))
            END IF
C     
C---- TESTOVER
C     
          ELSE IF (SUBKEY.EQ.'TEST') THEN
            DEBUG(56) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(56)= NINT(VALUE(I))
            END IF
C     
C---- COMPLETE
C     
          ELSE IF (SUBKEY.EQ.'COMP') THEN
            DEBUG(57) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(57)= NINT(VALUE(I))
            END IF
C     
C---- GETREJ
C     
          ELSE IF (SUBKEY.EQ.'GETR') THEN
            DEBUG(58) = .TRUE.
C     
C---- GETMOREBG
C     
          ELSE IF (SUBKEY.EQ.'GETM') THEN
            DEBUG(59) = .TRUE.
C     
C---- ALIGN
C     
          ELSE IF (SUBKEY.EQ.'ALIG') THEN
            DEBUG(60) = .TRUE.
C     
C---- RADBG
C     
          ELSE IF (SUBKEY.EQ.'RADB') THEN
            DEBUG(61) = .TRUE.
C     
C---- PICKSPOTS
C     
          ELSE IF (SUBKEY.EQ.'PICK') THEN
            DEBUG(62) = .TRUE.
            NDEBUG(62) = 0
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(62)= NINT(VALUE(I))
            END IF

C     *** WARNING DEBUG(63) USED FOR FINDSPOTS ABOVE ***
C     
C     
C---- REFIX
C     
          ELSE IF (SUBKEY.EQ.'REFI') THEN
            DEBUG(64) = .TRUE.
C     
C---- MXDSPL
C     
          ELSE IF (SUBKEY.EQ.'MXDS') THEN
            DEBUG(65) = .TRUE.
C     
C---- TIFF
C     
          ELSE IF (SUBKEY.EQ.'TIFF') THEN
            DEBUG(66) = .TRUE.
C     
C---- ADDSPOTS
C     
          ELSE IF (SUBKEY.EQ.'ADDS') THEN
            DEBUG(67) = .TRUE.
C     
C---- MODARRAY
C     
          ELSE IF (SUBKEY.EQ.'MODA') THEN
            DEBUG(68) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(68)= NINT(VALUE(I))
            END IF
C     
C---- CBF_INFO and associated routines
C     
          ELSE IF (SUBKEY.EQ.'CBFI') THEN
            DEBUG(69) = .TRUE.
C     
C---- DPS Index
C     
          ELSE IF (SUBKEY.EQ.'DPSI') THEN
            DEBUG(70) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(70)= NINT(VALUE(I))
            END IF
C     
C---- DDDN, PRDDET and associated routines
C     
          ELSE IF (SUBKEY.EQ.'DDDN') THEN
            DEBUG(71) = .TRUE.
C     
C---- INDREF and associated routines
C     
          ELSE IF (SUBKEY.EQ.'INDR') THEN
            DEBUG(72) = .TRUE.
C     
C---- Permutate
C     
          ELSE IF (SUBKEY.EQ.'PERM') THEN
            DEBUG(73) = .TRUE.
C
C---- CELLREFSET
C
          ELSE IF (SUBKEY.EQ.'CELL') THEN
            DEBUG(74) = .TRUE.
C
C---- NEW SPOT PROFILE READING/WRITING ROUTINES
C
          ELSE IF (SUBKEY.EQ.'SPWR') THEN
            DEBUG(75) = .TRUE.
C
C---- cell refinement routines c/w indexing
C
          ELSE IF (SUBKEY.EQ.'CELR') THEN
            DEBUG(76) = .TRUE.
C
C---- extra statistics from various routines - this should be used with an
C     '.AND.DEBUG(n)' for the appropriate routines.
C
          ELSE IF (SUBKEY.EQ.'STAT') THEN
            DEBUG(77) = .TRUE.
C
C---- LOCBKG local background
C
          ELSE IF (SUBKEY.EQ.'LOCL') THEN
            DEBUG(78) = .TRUE.
C
C---- IFPWDR for command prompt findspots,indexing
C
          ELSE IF (SUBKEY.EQ.'IFPW') THEN
            DEBUG(79) = .TRUE.
C
C---- BACKSTOP for automated backstop routine
C     ***** NOTE *****1
c     NDEBUG(80) is used with DUMP keywords
C
          ELSE IF (SUBKEY.EQ.'BACK') THEN
            DEBUG(80) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(80)= NINT(VALUE(I))
            END IF

C
C---- GETSPOTS and GETSEPRAS
C
          ELSE IF (SUBKEY.EQ.'GETS') THEN
C     
C---- Need to differentiate GETSPOTS and GETSEPRAS
C     
            KEY6 = LINE(IBEG(I):IEND(I))
C     
C---- convert to upper case
C     
C     ***********
            CALL CCPUPC(KEY6)
C     ***********
            IF (KEY6(1:6).EQ.'GETSPO') THEN
              DEBUG(81) = .TRUE.
            ELSE IF (KEY6.EQ.'GETSEP') THEN
              DEBUG(85) = .TRUE.
            END IF
C
C---- CGFIT and CGFIND
C
          ELSE IF (SUBKEY.EQ.'CGFI') THEN
C     
C---- Need to differentiate CGFIT and CGFIND
C     
            KEY6 = LINE(IBEG(I):IEND(I))
C     
C---- convert to upper case
C     
C     ***********
            CALL CCPUPC(KEY6)
C     ***********
            IF (KEY6(1:5).EQ.'CGFIT') THEN
              DEBUG(82) = .TRUE.
            ELSE IF (KEY6.EQ.'CGFIND') THEN
              DEBUG(83) = .TRUE.
            END IF
C
C---- SETDETECT  
C
          ELSE IF (SUBKEY.EQ.'SETD') THEN
            DEBUG(84) = .TRUE.
C
C---- PREDICT_SPOTS
C
          ELSE IF (SUBKEY.EQ.'PRED') THEN
            DEBUG(86) = .TRUE.
C
C---- OVERLAP
C
          ELSE IF (SUBKEY.EQ.'OVER') THEN
            DEBUG(87) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(87)= NINT(VALUE(I))
            END IF
C
C---- LATTICES Not a S/R, but debugging associated with multi-lattice
C              integration
C
          ELSE IF (SUBKEY.EQ.'LATT') THEN
            DEBUG(88) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(88)= NINT(VALUE(I))
            END IF
C
C---- SOCKET Not a subroutine, but debug output when a SOCKET_READ 
C            is performed
C
          ELSE IF (SUBKEY.EQ.'SOCK') THEN
            DEBUG(89) = .TRUE.
C     
C---- DPSCHOOSE
C     
          ELSE IF (SUBKEY.EQ.'DPSC') THEN
            DEBUG(90) = .TRUE.
C     
C---- test for a number
C     
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(90)= NINT(VALUE(I))
            END IF
C
C---- STAT statistics for this image (from dspscl.f)
C
         ELSE IF (SUBKEY.EQ.'STAT') THEN
           DEBUG(91) = .TRUE.
C
C---- XYZTOHKL
C
         ELSE IF (SUBKEY.EQ.'XYZT') THEN
           DEBUG(92) = .TRUE.
C
C---- SUMIMAGE .. image summation
C
          ELSE IF (SUBKEY.EQ.'SUMI') THEN
            DEBUG(92) = .TRUE.
C
C---- ROLLING .. Rolling shutter mode
C
          ELSE IF (SUBKEY.EQ.'ROLL') THEN
            DEBUG(93) = .TRUE.
C
C---- PHIOFFSET .. Rolling shutter modifying phi for post-refinement
C
          ELSE IF (SUBKEY.EQ.'PHIO') THEN
            DEBUG(94) = .TRUE.
C
C---- INVNOR .. Normal matrix 
C
          ELSE IF (SUBKEY.EQ.'INVN') THEN
            DEBUG(95) = .TRUE.
C
C
C---- PIXEL .. Dump "pick" values to logfile so that statistics on background
C              pixel can be calculated.
C
          ELSE IF (SUBKEY.EQ.'PIXE') THEN
            DEBUG(96) = .TRUE.
C
C
C---- DSPSCL  Look for isolated pixels with value -1, used when debugging Eiger images
C
          ELSE IF (SUBKEY.EQ.'DSPS') THEN
            DEBUG(97) = .TRUE.
C     
C---- CLUSTER
C     
          ELSE IF (SUBKEY.EQ.'CLUS') THEN
            DEBUG(98) = .TRUE.
C     
C---- FLOW not a subroutine, but prints info on program flow
C     
          ELSE IF (SUBKEY.EQ.'FLOW') THEN
            DEBUG(99) = .TRUE.
C
C---- SIGMA  (not a subroutine, but flag for new sigma calculation
C
          ELSE IF (SUBKEY.EQ.'SIGM') THEN
            DEBUG(100) = .TRUE.
C
C---- EIGER  
C
          ELSE IF (SUBKEY.EQ.'EIGE') THEN
            DEBUG(101) = .TRUE.
C
C---- GETHDR
C
          ELSE IF (SUBKEY.EQ.'GETH') THEN
            DEBUG(102) = .TRUE.
C
C---- SECTION
C
          ELSE IF (SUBKEY.EQ.'SECT') THEN
            DEBUG(103) = .TRUE.
C
C---- GETGAPS
C
          ELSE IF (SUBKEY.EQ.'GETG') THEN
            DEBUG(104) = .TRUE.
C
C---- READOXIMAGE
C
          ELSE IF (SUBKEY.EQ.'READ') THEN
            DEBUG(105) = .TRUE.
C
C---- SEGMENT
C
          ELSE IF (SUBKEY.EQ.'SEGM') THEN
            DEBUG(106) = .TRUE.
C
C---- REFLFIND
C
          ELSE IF (SUBKEY.EQ.'REFL') THEN
            DEBUG(107) = .TRUE.
C
C---- IMAGE <number> Do not turn on debug until image specified has been reached.
C
          ELSE IF (SUBKEY.EQ.'IMAG') THEN
            DEBUG(108) = .TRUE.
            NDEBUG(108) = 1
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              NDEBUG(108)= NINT(VALUE(I))
            END IF
C
C---- ADJACENT debug about adjacent spots intruding in standard profile
C              code is in BESTMASK
C
          ELSE IF (SUBKEY.EQ.'ADJA') THEN
            DEBUG(109) = .TRUE.
C
C---- GETDETECT
C
          ELSE IF (SUBKEY.EQ.'GETD') THEN
            DEBUG(110) = .TRUE.
C
C---- Subroutine name not recognised
C     
          ELSE               
            IF (ONLINE) WRITE (ITOUT,FMT=6042) SUBKEY
 6042       FORMAT (1X,'*** Subroutine ',A,' not recognised')
            WRITE (IOUT,FMT=6042) SUBKEY
            BADLINE = .TRUE.
          END IF
          GOTO 370
 372      CONTINUE
          WRITE (IOUT,FMT=6040) DEBUGSTR(1:LENSTR(DEBUGSTR))
 6040     FORMAT (1X,'DEBUG output will be printed in following ',
     +         'subroutines:', /,1X,A)
          IF (ONLINE) WRITE(ITOUT,FMT=6040)
     $         DEBUGSTR(1:LENSTR(DEBUGSTR))
          IF (DEBUG(99)) THEN
            WRITE(IOUT,FMT=6041)
            IF (ONLINE) WRITE(ITOUT,FMT=6041)
 6041       FORMAT(1X,'*** Information on the program flow will be',
     +                ' provided (not complete) ***.')
          END IF
          IF (DEBUG(108)) THEN
            WRITE(IOUT,FMT=6044) NDEBUG(108)
            IF (ONLINE) WRITE(ITOUT,FMT=6044) NDEBUG(108)
 6044       FORMAT(1X,'*** Debug output will only be turned on after',
     +        ' image number',I6,' ***')
          END IF
        END IF
C
C---- NODEBUG turns off all debugging
C
      ELSE IF (KEY6.EQ.'NODEBU') THEN
        DO 373 I = 1,200
          DEBUG(I) = .FALSE.
 373    CONTINUE
        WRITE(IOUT,FMT=6043)
        IF (ONLINE) WRITE(ITOUT,FMT=6043)
 6043   FORMAT(/,1X,'*** DEBUGGING TURNED OFF FOR ALL SUBROUTINES ***')
C     
C---- Bad
C     write specified reflections to badspots file
C     
      ELSE IF (KEY.EQ.'BAD ') THEN
        DUMPSPOT = .TRUE.
        IF (NTOK.GT.1) THEN
          ICOUNT = 1
 374      ICOUNT = ICOUNT + 1
          SUBKEY = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     *********** 
          CALL CCPUPC(SUBKEY)
C     ***********       
          IF (SUBKEY(1:3).EQ.'ALL') THEN
            DUMPALL = .TRUE.
C     
C---- START Only dump reflections after NDSTART in the list sorted on
C     X scanner (pixel) coordinate
C     
          ELSE IF (SUBKEY.EQ.'STAR') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NDSTART = NINT(VALUE(ICOUNT))
            DUMPALL = .TRUE.
C     
C---- TOTAL maximum number of reflections to be dumped
C     
          ELSE IF (SUBKEY.EQ.'TOTA') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NDTOT = NINT(VALUE(ICOUNT))
            DUMPALL = .TRUE.
C     
C---- IXMIN  minimum X pixel coordinate
C     
          ELSE IF (SUBKEY.EQ.'IXMI') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            IXDMIN = NINT(VALUE(ICOUNT))
            DUMPALL = .TRUE.
C     
C---- IXMAX  maximum X pixel coordinate
C     
          ELSE IF (SUBKEY.EQ.'IXMA') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            IXDMAX = NINT(VALUE(ICOUNT))
            DUMPALL = .TRUE.
C     
C---- IYMIN  minimum Y pixel coordinate
C     
          ELSE IF (SUBKEY.EQ.'IYMI') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            IYDMIN = NINT(VALUE(ICOUNT))
            DUMPALL = .TRUE.
C     
C---- IYMAX  maximum Y pixel coordinate
C     
          ELSE IF (SUBKEY.EQ.'IYMA') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            IYDMAX = NINT(VALUE(ICOUNT))
            DUMPALL = .TRUE.
C     
C---- RESOL  Resolution limts
C     
          ELSE IF (SUBKEY.EQ.'RESO') THEN
C     *******************************************
            CALL MKEYNM(2,ICOUNT+1,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            RESDLOW = VALUE(ICOUNT+1)
            RESD = VALUE(ICOUNT+2)
            IF (RESD.GT.RESDLOW) THEN
              TEMPX = RESD
              RESD = RESDLOW
              RESDLOW = TEMPX
            END IF
            ICOUNT = ICOUNT + 2
C     
C---- IMIN  Minimum intensity
C     
          ELSE IF (SUBKEY.EQ.'IMIN') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            IDMIN = NINT(VALUE(ICOUNT))
C     
C---- IMAX  Minimum intensity
C     
          ELSE IF (SUBKEY.EQ.'IMAX') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            IDMAX = NINT(VALUE(ICOUNT))
C     
C---- HKL to dump specific reflections
C     
          ELSE IF (SUBKEY.EQ.'HKL') THEN
C     *******************************************
            CALL MKEYNM(3,ICOUNT+1,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            IF (NHKLD.LE.49) THEN
              NHKLD = NHKLD + 1
              DO 375  I = 1,3
                IHD(I,NHKLD) = NINT(VALUE(ICOUNT+I))
 375          CONTINUE
              ICOUNT = ICOUNT + 3
            END IF
            DUMPALL = .FALSE.
            NDTOT = NHKLD
          ELSE
C     
C---- Not recognised
C     
            IF (TRAPERR) INPERR = .TRUE.
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          END IF
          IF (ICOUNT.LT.NTOK) GOTO 374
          IF (DUMPALL) THEN
            WRITE (IOUT,FMT=6047) NDTOT,NDSTART,IXDMIN,IXDMAX,IYDMIN,
     +           IYDMAX,RESD,RESDLOW,IDMIN,IDMAX
            IF (ONLINE) WRITE (ITOUT,FMT=6047) NDTOT,NDSTART,IXDMIN,
     +           IXDMAX,IYDMIN,IYDMAX,RESD,RESDLOW,IDMIN,IDMAX
 6047       FORMAT(1X,'Dumping up to',I5,' reflections starting',
     +           ' at number',I5,' in list sorted on scanner',
     +           /,1X,'X coordinate and having pixel coordinates',
     +           ' between',I5,' and',I5,' in X',/,1X,'and',I5,
     +           ' and',I5,' in Y and between resolution limits',
     +           F6.2,' and',F6.2,/,1X,'and with intensities ',
     +           'between',I6,' and',I10)
          ELSE
            WRITE (IOUT,FMT=6046) ((IHD(I,J),I=1,3),J=1,NHKLD)
 6046       FORMAT (1X,'Extensive debug data will be printed for ',
     +           'the following reflections:',/,(1X,3I5))
            IF(ONLINE)WRITE(ITOUT,FMT=6046)
     $           ((IHD(I,J),I=1,3),J=1,NHKLD)
          END IF
        END IF
C     
C---- Graphics...specifies graphics device (sigma or tectronix)
C     
      ELSE IF (KEY.EQ.'GRAP') THEN
        GRTYPE = LINE(IBEG(2) :IEND(2))
C     
C     **************
        CALL CCPUPC(GRTYPE)
C     **************
C     
        IF (GRTYPE.EQ.'SIG') THEN
          NGR = 1
          IF (ONLINE) WRITE (ITOUT,FMT=6048)
 6048     FORMAT (1X,'Graphics device is SIGMA ')
          WRITE (IOUT,FMT=6048)
        ELSE IF (GRTYPE(1:2).EQ.'SG') THEN
          NGR = 6
        ELSE IF ((GRTYPE(1:3).EQ.'TEK') .OR.
     +         (GRTYPE(1:2).EQ.'TX')) THEN
          NGR = 2
C     
C     
          IF (NTOK.EQ.3) THEN
C     
C     *******************************************
            CALL MKEYNM(1,3,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NLI = NINT(VALUE(3))
          END IF
C     
C     
          IF (ONLINE) WRITE (ITOUT,FMT=6050)
 6050     FORMAT (1X,'Graphics device is TEKTRONIX ')
          WRITE (IOUT,FMT=6050)
          IF (NTOK.EQ.3) WRITE (ITOUT,FMT=6052) NLI
 6052     FORMAT (1X,'Number of lines on screen is',I4)
        ELSE IF (GRTYPE(1:1).EQ.'X') THEN
          NGR = 7
        ELSE
          IF (ONLINE) WRITE (ITOUT,FMT=6054)
 6054     FORMAT (1X,'*** UNKNOWN GRAPHICS, Must be SIG or TX (TEK)')
          WRITE (IOUT,FMT=6054)
        END IF
C     
C     
C---- PROFILE...use profile fitting for intensity measurement
C     
      ELSE IF (KEY.EQ.'NOPR') THEN
        PROFILE = .FALSE.
        ACCUMULATE = .FALSE.
        INTERPOL = .FALSE.
C     
C     
      ELSE IF (KEY.EQ.'PROF') THEN
C     
cdms21.09.2004   Be sure PROFILE is set.
         PROFILE = .TRUE.     
         IF (NTOK.GT.1) THEN
          I = 1
 380      I = I + 1
          IF (I.GT.NTOK) GO TO 390
          STR = LINE(IBEG(I) :IEND(I))
C     
C     ***********
          CALL CCPUPC(STR)
C     ***********
C     
          IF (STR(1:3).EQ.'OFF') THEN
            PROFILE = .FALSE.
            ACCUMULATE = .FALSE.
            INTERPOL = .FALSE.
          ELSE IF (STR(1:4).EQ.'PRIN') THEN
            LPRINT(11) = .TRUE.
          ELSE IF (STR(1:5).EQ.'NOPRI') THEN
            LPRINT(11) = .FALSE.
          ELSE IF (STR(1:4).EQ.'PART') THEN
            PRPART = .TRUE.
            IF (PRFULLS) THEN
              WRITE(IOUT,FMT=7050)
              IF (ONLINE) WRITE(ITOUT,FMT=7050)
 7050         FORMAT(/,1X,'*** ERROR ***',/,1X,'Cannot include',
     +             ' keywords "PARTIALS" and "FULLS", they are',
     +             ' mutually exclusive')
            NSHUTERR = 31
              CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
            END IF
          ELSE IF (STR(1:4).EQ.'FULL') THEN
            PRFULLS = .TRUE.
            IF (PRPART) THEN
              WRITE(IOUT,FMT=7050)
              IF (ONLINE) WRITE(ITOUT,FMT=7050)
            NSHUTERR = 32
              CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
            END IF
          ELSE IF (STR(1:4).EQ.'CHAN') THEN
            CHANGEMASK = .TRUE.
          ELSE IF (STR(1:4).EQ.'NOBF') THEN
            PRBFILM = .FALSE.
          ELSE IF (STR(1:4).EQ.'NOCF') THEN
            PRCFILM = .FALSE.
          ELSE IF (STR(1:4).EQ.'BFIL') THEN
            PRBFILM = .TRUE.
          ELSE IF (STR(1:4).EQ.'CFIL') THEN
            PRCFILM = .TRUE.
C     
          ELSE IF (STR(1:4).EQ.'SAVE') THEN
            PRSAVE = .TRUE.
C     
C     
C     
C---- Get the template for the profile filenames. If a "." is present,
C     strip off everything after (and including) the "."
            FWORK = ' '
            CALL UGTENV('PROFILE',FWORK)
c     -vms            inquire (file='profile',exist=efile)
            NCH = LENSTR(FWORK)
            DO 381 J=1,NCH
              IF (FWORK(J:J).EQ.'.') THEN
                NCH4 = J - 1
                GOTO 382
              END IF
 381        CONTINUE
            NCH4 = NCH
 382        PROFFNW = FWORK(1:NCH4)
C     
          ELSE IF (STR(1:4).EQ.'READ') THEN
            PRREAD = .TRUE.
C     
C---- Get template name for profiles
C     have an extension and second must add the standard extension.
C     
            FWORK = ' '
            CALL UGTENV('PROFILE',FWORK)
c     -vms            inquire (file='profile',exist=efile)
            NCH = LENSTR(FWORK)
            DO 383 J=1,NCH
              IF (FWORK(J:J).EQ.'.') THEN
                NCH3 = J - 1
                GOTO 384
              END IF
 383        CONTINUE
            NCH3 = NCH
 384        PROFFNR = FWORK(1:NCH3)//'_001.PRF'
C     
C---- Make lower case for unix systems
C     
C     **************
            CALL CCPLWC(PROFFNR)
C     **************
            INQUIRE (FILE=PROFFNR,EXIST=EFILE)
            NCH3 = LENSTR(PROFFNR)
            IF (.NOT.EFILE) THEN
              NWRN = NWRN + 1
              WRITE (IOUT,FMT=6060) PROFFNR(1:NCH3)
 6060         FORMAT (//,1X,'**** ERROR ****',/,1X,
     $             'The file containing the prof',
     +             'iles to be read does not exist',
     $             /,1X,'Filename is ',A)
              IF (ONLINE) THEN
                WRITE (ITOUT,FMT=6060) PROFFNR(1:NCH3)
                GO TO 50
              END IF
              STOP
            END IF
          ELSE IF (STR(1:4).EQ.'ACCU') THEN
            ACCUMULATE = .TRUE.
C     
          ELSE IF (STR(1:4).EQ.'INTE') THEN
            INTERPOL = .TRUE.
          ELSE IF (STR(1:4).EQ.'NOIN') THEN
            INTERPOL = .FALSE.
C     
          ELSE IF (STR(1:4).EQ.'NOAC') THEN
            ACCUMULATE = .FALSE.
          ELSE IF (STR(1:4).EQ.'NOCH') THEN
            CHANGEMASK = .FALSE.
          ELSE IF (STR(1:4).EQ.'ISDR') THEN
C     
C     *********************************************
            CALL MKEYNM(1,I+1,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
            ISDRATIO = NINT(VALUE(I+1))
            I = I + 1
          ELSE IF (STR(1:4).EQ.'BGSI') THEN
C     
C---- This should not be changed from default...force use of MODIFY 
C     subkeyword and give warning message
C     
            IF (ITYP(I+1).NE.1) THEN
C     
C---- No MODIFY keyword present
C     
              WRITE(IOUT,FMT=6614)
              IF (ONLINE) WRITE(ITOUT,FMT=6614)
              IF (BRIEF) WRITE(IBRIEF,FMT=6614)
 6614         FORMAT(1X,'*** WARNING ***. With the new program ',
     +             ' the parameter BGSIG for profiles',/,1X,
     +             ' should NOT be changed from its default.',
     $             ' If you really do want',
     +             ' do want to change it',/,1X,
     $             'you must include the subkeyword MODIFY ',
     $             'before the numeric value.',/,1X,'THIS IS ',
     +             ' *** NOT *** RECOMMENDED')
              I = I + 1
              GOTO 380
            ELSE IF (ITYP(I+1).EQ.1) THEN
              I = I + 1
              SUBKEY = LINE(IBEG(I) :IEND(I))
C     
C     ************
              CALL CCPUPC(SUBKEY)
C     ************
              IF (SUBKEY.NE.'MODI') THEN
                WRITE(IOUT,FMT=6616)
                IF (ONLINE) WRITE(ITOUT,FMT=6616)
                IF (BRIEF) WRITE(IBRIEF,FMT=6616)
 6616           FORMAT(1X,'*** INVALID SUBKEYWORD AFTER BGSIG ***',
     +               /,1X,'The only valid subkeyword is MODIFY')
                GOTO 380
              ELSE
C     *********************************************
                CALL MKEYNM(1,I+1,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
                PRBGSIG = VALUE(I+1)
                I = I + 1
              END IF
            END IF
C     
C---- CUTOFF... sets level for "overloads" when forming profiles
C     
          ELSE IF (STR(1:4).EQ.'CUTO') THEN
C     
C     *********************************************
            CALL MKEYNM(1,I+1,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
            PRCUTOFF = NINT(VALUE(I+1)) - 1
            IPRCUT = 2
            I = I + 1
C     
C---- NREF.... minimum number of reflections in a profile
C     
          ELSE IF (STR(1:4).EQ.'NREF') THEN
C     
C     *********************************************
            CALL MKEYNM(1,I+1,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
            NRFMIN = NINT(VALUE(I+1))
            I = I + 1
          ELSE IF (STR(1:4).EQ.'RMSB') THEN
C     
C     *********************************************
            CALL MKEYNM(1,I+1,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
            RMSBGPR = VALUE(I+1)
            I = I + 1
C     
C---- PLOT.... not used
C     
          ELSE IF (STR(1:4).EQ.'PLOT') THEN
C     
C     *********************************************
            CALL MKEYNM(1,I+1,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
            IPLOT = NINT(VALUE(I+1))
            I = I + 1
C     
C---- WEIGHT...use weighting in fitting profile to individual spots
C     
          ELSE IF (STR(1:4).EQ.'WEIG') THEN
            WEIGHT = .TRUE.
C     
C---- WSUM....use weighting in formation of standard profiles from
C     individual spots
C     
          ELSE IF (STR(1:4).EQ.'WSUM') THEN
            WTPROFILE = .TRUE.
          ELSE IF (STR(1:4).EQ.'NOWS') THEN
            WTPROFILE = .FALSE.
          ELSE IF (STR(1:4).EQ.'OVER') THEN
            USEOVRLD = .TRUE.
          ELSE IF (STR(1:4).EQ.'NOWE') THEN
            WEIGHT = .FALSE.
          ELSE IF (STR(1:4).EQ.'NOOV') THEN
            USEOVRLD = .FALSE.
          ELSE IF (STR(1:4).EQ.'EDGE') THEN
            USEDGE = .TRUE.
C     
C---- VARPROFILES...calculate a profile for each spot as a weighted sum
C     of neighbouring "standard" profiles
C     
          ELSE IF (STR(1:4).EQ.'VARI') THEN
            VARPRO = .TRUE.
          ELSE IF (STR(1:4).EQ.'NOVA') THEN
            VARPRO = .FALSE.
C     
C---- OPTIMISE...optimise the raster parameters (rim and corner cutoff)
C     for each standard profile. Set PROPTCEN TRUE as well, in case it
C     has previously been set FALSE by NOOPT ATALL
C     
          ELSE IF (STR(1:4).EQ.'OPTI') THEN
            PROPT = .TRUE.
            PROPTCEN = .TRUE.
          ELSE IF (STR(1:4).EQ.'NOOP') THEN
            PROPT = .FALSE.
C     
C---- Check if optimisation is also to be turned of for the central area
C     
            IF (I.LT.NTOK) THEN
              IF (ITYP(I+1).EQ.1) THEN
C     
                SUBKEY = LINE(IBEG(I+1) :IEND(I+1))
C     **************
                CALL CCPUPC(SUBKEY)
C     ***************
                IF (SUBKEY.EQ.'ATAL') THEN
                  PROPTCEN = .FALSE.
                  I = I +1
                ELSE
                  PROPTCEN = .TRUE.
                END IF
              END IF
            END IF
            
C     
C---- TOLERANCE... Used in optimisation of measurement box. Maximum
C     acceptable loss in total intensity
C     
          ELSE IF (STR(1:4).EQ.'TOLE') THEN
            ITOL = 1
C     
C---- Need to check for presence of one or two numbers. If only one,
C     then this is maximum tolerance, it two then the first is the
C     minimum tolerance and the second is the maximum
C     
            IF (((I+2).LE.NTOK).AND.(ITYP(I+2).EQ.2)) THEN
              I = I + 1
C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
              TOLMIN = VALUE(I)
            END IF
            I = I + 1
C     *********************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
            TOL = VALUE(I)
C
C---- Now check if the next subkeyword is UPDATE, in which case
C     ipmosflm is allowed to increase the TOLERANCE to avoid overlaps.
C
            IF (ITYP(I+1).EQ.1) THEN
              SUBKEY = LINE(IBEG(I+1):IEND(I+1))
              IF (SUBKEY.EQ.'UPDA') THEN
C
C---- Allow updating of TOLERANCE by CHKRAS
C     Correction. CHKRAS no longer uses ITOL. If it is set to zero here
C     then the second value (max tolerance) is ignored becuase TOL is set to 0.03
C     below (its default) if ITOL=0, so comment out next line AGWL 23/7/13
CAL                ITOL = 0
                I = I + 1
              ELSE IF (SUBKEY.EQ.'NOUP') THEN
C
C---- Do not allow updating of TOLERANCE by CHKRAS
C
                ITOL = 2
                I = I + 1
              END IF
            END IF
C     
C---- BADTOLERANCE... Used in optimisation of measurement box. Maximum
C     fraction of peak pixels that may be part of neighbouring peak
C     
          ELSE IF (STR(1:4).EQ.'BADT') THEN
            I = I + 1
C     *********************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
            BADTOL = VALUE(I)
C     
C---- BOUNDARY.... Used in optimisation of measurement box. Amount
C     (in pixels) by which boundary is to be extended after finding the
C     point which gives maximum I/sigma(I)
C     
          ELSE IF (STR(1:4).EQ.'BOUN') THEN
            I = I + 1
C     *********************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
            IBOUND = NINT(VALUE(I))
C     
C---- FIXBOX...don't optimise the overall size of box
C     
          ELSE IF (STR(1:4).EQ.'FIXB') THEN
            FIXBOX = .TRUE.
            NOFIXBOX = .FALSE.
C     
C---- NOFIXBOX...don't fix overall box size
C     
          ELSE IF (STR(1:4).EQ.'NOFI') THEN
            NOFIXBOX = .TRUE.
            FIXBOX = .FALSE.
C     
C---- RATIO... minimum ratio of number of background pixels to
C     peak pixels used in overall box size optimisation
C     
          ELSE IF (STR(1:4).EQ.'RATI') THEN
            I = I + 1
C     *********************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
            BGPKRAT = VALUE(I)
C     
C---- STOP.... When optimising overall dimensions of box, stop
C     expanding in a particular direction if a fraction greater
C     than STOP of the new background pixels are rejected
C     
          ELSE IF (STR(1:4).EQ.'STOP') THEN
            I = I + 1
C     *********************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
            FRACREJ = VALUE(I)
C     
C---- DISCRIMINATE...used to eliminate spots with strong neighbours
C     from formation of profiles
C     
          ELSE IF (STR(1:4).EQ.'DISC') THEN
            DISCRIMINATE = .TRUE.
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1

C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
              DISCRIM = VALUE(I)
            END IF
C     
C---- DECONVOLUTE to deconvolute the standard profiles and down-weight
C     peak edge pixels whihc may be overlapped by the peak from a 
C     neighbouring spot.
C     ****** NOT YET PROPERLY IMPLEMENTED ***
C     
          ELSE IF (STR(1:4).EQ.'DECO') THEN
            DECONV = .TRUE.
C     
C     
C---- WDLIM1 Controls rejection of peak pixels from profile fitting
C     Pixels are rejected if W*DEL**2 is greater than PKWDLIM1,
C     but for pixels adjacent to overlapped pixels the test 
C     is against PKWDLIM2.
C     
          ELSE IF (STR(1:6).EQ.'WDLIM1') THEN
            I = I + 1
C     *********************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
            PKWDLIM1 = VALUE(I)
C     
          ELSE IF (STR(1:6).EQ.'WDLIM2') THEN
            I = I + 1
C     *********************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
            PKWDLIM2 = VALUE(I)
C     
          ELSE IF (STR(1:6).EQ.'WDLIM3') THEN
            I = I + 1
C     *********************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
            PKWDLIM3 = VALUE(I)
C     
C---- OUTLIERS Controls rejection of peak pixels that have pixel values
C     in a 
C     specified range (IOUTL1 to IOUTL2) where the profile fit is worse
C     than PKWDOUTL. This is to deal with "speckles" due to switch-over
C     of dynamic range on Mar scanners, and could also deal with zero
C     spiral "strips" due to SCSI errors.
C     
          ELSE IF (STR(1:4).EQ.'OUTL') THEN
            I = I + 1
C     *********************************************
            CALL MKEYNM(3,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
            IOUTL1 = NINT(VALUE(I))
            IOUTL2 = NINT(VALUE(I+1))
            PKWDOUTL = VALUE(I+2)
            I = I + 2
C     
C---- PRUPDATE... if true, for first block of images go back
C     and redetermine profiles after optimising raster box parameters
C     
          ELSE IF (STR(1:4).EQ.'PRUP') THEN
            PUPDATE = .TRUE.
          ELSE IF (STR(1:5).EQ.'NOPRU') THEN
            PUPDATE = .FALSE.
C
C---- ICE, exclude reflections that may lie in ice rings from profiles
C
          ELSE IF (STR(1:3).EQ.'ICE') THEN
            PREXICE = .TRUE.
C
C---- Test for a number that gives the half-widths of the exclusion shells
C     as a fraction of the actual resolution.
C
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              I = I + 1

C     *********************************************
              CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
              ICEWIDTH = VALUE(I)
            END IF
C
C---- No ICE present, don't exclude resolution shells (unless specified
C     by RESOLUTION EXCLUDE
C
          ELSE IF (STR(1:4).EQ.'NOIC') THEN
            PREXICE = .FALSE.
C
C---- PRCUT. The fraction PRCUT reflections with the highest pixel
C     values in each standard profile area will be excluded from profile formation
C     default 0.05
C
          ELSE IF (STR(1:5).EQ.'PRCUT') THEN
            I = I + 1
C     *********************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
            PRCUT = VALUE(I)
            IF (PRCUT.GT.1.0) PRCUT = 0.01*PRCUT
C
C---- FBGNULL. Reflections with a fraction FBGNULL of background pixels having
C     a value <= NULLPIX will be rejections outright from integration
C
          ELSE IF (STR(1:4).EQ.'FBGN') THEN
            I = I + 1
C     *********************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
            IF ((VALUE(I).GT.0.0).AND.(VALUE(I).LT.1.0)) 
     +                          FBGNULL = VALUE(I)
C
C---- FPKNULL. Reflections with a fraction FPKNULL of peak pixels having
C     a value <= NULLPIX will be rejections outright from integration
C
          ELSE IF (STR(1:4).EQ.'FPKN') THEN
            I = I + 1
C     *********************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
            IF ((VALUE(I).GT.0.0).AND.(VALUE(I).LT.1.0)) 
     +                         FPKNULL = VALUE(I)
C
C---- FBGEDGE. Reflections with a fraction FBGEDGE of background pixels having
C     a value <= NULLPIX will be rejections outright from integration
C
          ELSE IF (STR(1:4).EQ.'FBGE') THEN
            I = I + 1
C     *********************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
            IF ((VALUE(I).GT.0.0).AND.(VALUE(I).LT.1.0)) 
     +                          FBGEDGE = VALUE(I)
C
C---- FPKEDGE. Reflections with a fraction FPKEDGE of peak pixels having
C     a value <= NULLPIX will be rejections outright from integration
C
          ELSE IF (STR(1:4).EQ.'FPKE') THEN
            I = I + 1
C     *********************************************
            CALL MKEYNM(1,I,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
            IF ((VALUE(I).GT.0.0).AND.(VALUE(I).LT.1.0))
     +                           FPKEDGE = VALUE(I)
C
C---- PKONLY
C
          ELSE IF (STR(1:5).EQ.'PEAKO') THEN
            PKONLY = .TRUE.
C     
          ELSE IF (STR(1:4).EQ.'HIGH') THEN
            PRSET = .TRUE.
            HIGHRES = .TRUE.
            LOWRES = .FALSE.
          ELSE IF (STR(1:4).EQ.'LOWR') THEN
            PRSET = .TRUE.
            LOWRES = .TRUE.
            HIGHRES = .FALSE.
          ELSE IF (STR(1:4).EQ.'XLIN') THEN
C     
C---- Profile boundaries (in X direction in mm wrt lower left
C     hand corner of image looking towards source)
C     
            NXLINE = 0
            PRSET = .TRUE.
            LINESET = .TRUE.
C     
C---- test for a number             
C     
 388        CONTINUE
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              NXLINE = NXLINE + 1
              IF (NXLINE.GT.NNLINE) THEN
                WRITE(IOUT,FMT=6062) NNLINE
                IF (ONLINE) WRITE(ITOUT,FMT=6062) NNLINE
 6062           FORMAT(1X,'*** fatal error ***',/,1X,'You have ',
     +               'exceeded the maximum number of boundaries(',
     +               I2,')',/,1X,'Either give fewer boundaries o',
     +               'r change parameter NNLINE in include file ',
     +               '"parameter"',/,1X,
     $               'and recompile the program')
              END IF
C     *********************************************
              CALL MKEYNM(1,I+1,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C---- Convert to 10 micron units
C     
              XLINE(NXLINE) = VALUE(I+1)*100.0
C     
C---- Check they are given in increasing order
C     
              IF (NXLINE.GT.1) THEN   
                IF (XLINE(NXLINE).LT.XLINE(NXLINE-1)) THEN
                  WRITE(IOUT,FMT=6073)
                  IF (ONLINE) WRITE(ITOUT,FMT=6073)
 6073             FORMAT(//,1X,'**** FATAL ERROR ***',/,1X,
     $                 'Coordinates of profile',
     +                 ' boundaries MUST be given in ',
     $                 'increasing magnitude')
                  IF (ONLINE) THEN
                    BADLINE = .TRUE.
                    GOTO 50
                  ENDIF
                  STOP
                END IF
              END IF
              I = I + 1
              GOTO 388
            ELSE
              GOTO 380
            END IF            
          ELSE IF (STR(1:4).EQ.'YLIN') THEN
C     
C---- Profile boundaries (in X direction in mm wrt lower left
C     hand corner of image looking towards source)
C     
            NYLINE = 0
            PRSET = .TRUE.
            LINESET = .TRUE.
C     
C---- test for a number
C     
 389        CONTINUE  
            IF (I.LT.NTOK .AND. ITYP(I+1).EQ.2) THEN
              NYLINE = NYLINE + 1
              IF (NYLINE.GT.NNLINE) THEN
                WRITE(IOUT,FMT=6062) NNLINE
                IF (ONLINE) WRITE(ITOUT,FMT=6062) NNLINE
              END IF
C     *********************************************
              CALL MKEYNM(1,I+1,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C---- Convert to 10 micron units
C     
              YLINE(NYLINE) = VALUE(I+1)*100.0
C     
C---- Check they are given in increasing order
C     
              IF (NYLINE.GT.1) THEN   
                IF (YLINE(NYLINE).LT.YLINE(NYLINE-1)) THEN
                  WRITE(IOUT,FMT=6073)
                  IF (ONLINE) WRITE(ITOUT,FMT=6073)
                  IF (ONLINE) THEN
                    BADLINE = .TRUE.
                    GOTO 50
                  ENDIF
                  STOP
                END IF
              END IF
              I = I + 1
              GOTO 389
            ELSE
              GOTO 380
            END IF            
          ELSE
            WRITE (IOUT,FMT=6130) STR
 6130       FORMAT (//,1X,'********** Sub-keyword NOT Recognised:',A)
            IF (ONLINE) WRITE (ITOUT,FMT=6130) STR
          END IF
          GO TO 380
        END IF
 390    PROFILE = .TRUE.
C     
C---- PROCESS...process a preexisting mosflm.out file ** REMOVED **
C     
C     AL      ELSE IF (KEY.EQ.'PROC') THEN
C     AL        PROCES = .TRUE.
C     AL        WRITE (IOUT,FMT=6070)
C     AL 6070 FORMAT (1X,'Going straight into PROCESS, assumes mosflm
C     .out alre',
C     AL     +       'ady exists')
C     AL        IF (ONLINE) WRITE (ITOUT,FMT=6070)
C     
C---- DUMP... additional debug output allows spot pixel values, profiles
C     and data for each reflection to be output, subject to a minimum
C     intensity if appropriate.
C     Possible subkeywords: REFL [N], SPOT, PROF, ODS, IMIN X, TOTAL N,
C     BGR, OVER, EDGE, PKREJ N, BGREJ N, SPOTLIST, HKL (3 indices),
C     RESO (low and high, either order), IMAGE (start dumping at image)
C     
      ELSE IF (KEY.EQ.'DUMP') THEN
        DUMPSPOT = .TRUE.
        IF (MXDUMP.EQ.0) MXDUMP = 100
        IF (NTOK.GT.1) THEN
          ICOUNT = 1
 331      ICOUNT = ICOUNT + 1
          STR = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C         ***********
          CALL CCPUPC(STR)
C         ***********
C     
          IF (STR(1:4).EQ.'REFL') THEN
            DUMP(1) = .TRUE.
C     
C---- test for a number N
C     
            IF (ICOUNT.LT.NTOK .AND. ITYP(ICOUNT+1).EQ.2) THEN
              ICOUNT = ICOUNT + 1
C     
C     ************************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
              NDUMP = NINT(VALUE(ICOUNT))
            END IF
C
C---- Summed pixel values for forming profiles and profiles
C     
          ELSE IF (STR(1:4).EQ.'PROF') THEN
            DUMP(2) = .TRUE.
C     
C---- Not yet implemented...use SPOT
C     
          ELSE IF (STR(1:3).EQ.'ODS') THEN
            DUMP(3) = .TRUE.
C
C---- Fraction of background pixels rejected greater than BGRLIM
C
          ELSE IF (STR(1:4).EQ.'BGRF') THEN
            DUMP(4) = .TRUE.

C     
C---- test for a number for BGRLIM
C     
            IF (ICOUNT.LT.NTOK .AND. ITYP(ICOUNT+1).EQ.2) THEN
              ICOUNT = ICOUNT + 1
C     
C     ************************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
              BGRLIM = VALUE(ICOUNT)
            END IF
C
C---- Overloads
C
          ELSE IF (STR(1:4).EQ.'OVER') THEN
            DUMP(5) = .TRUE.
C
C---- Spots where measurement box contains pixels with value NULLPIX
C
          ELSE IF ((STR(1:4).EQ.'EDGE').AND.(STR(1:7).NE.'EDGEBAD')
     +                      .AND.(STR(1:7).NE.'EDGEREJ')) THEN
            DUMP(6) = .TRUE.
C     
C-----SPOTLIST: Write a list of spot coordinates as output by IMSTILLS
C     
          ELSE IF (STR(1:5).EQ.'SPOTL') THEN
            DUMP(8) = .TRUE.
            WRITE(IOUT,FMT=7120)
            IF (ONLINE) WRITE(ITOUT,FMT=7120)
 7120       FORMAT(1X,'A spotlist will be written')
C
C---- Fully recorded reflections
C
          ELSE IF (STR(1:4).EQ.'FULL') THEN
            DUMP(9) = .TRUE.
C     
C---- PKREJ N
C     Dump spots with more than N peak pixels rejected
C     
          ELSE IF (STR(1:4).EQ.'PKRE') THEN
            DUMP(10) = .TRUE.
C     
C---- test for a number N
C     
            IF (ICOUNT.LT.NTOK .AND. ITYP(ICOUNT+1).EQ.2) THEN
              ICOUNT = ICOUNT + 1
C     
C     ************************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
              IARRDUMP(10) = NINT(VALUE(ICOUNT))
            END IF
C     
C---- DUMP BGREJ N
C     Dump spots with more than N background pixels rejected
C     
          ELSE IF (STR(1:4).EQ.'BGRE') THEN
            DUMP(11) = .TRUE.
C     
C---- test for a number N
C     
            IF (ICOUNT.LT.NTOK .AND. ITYP(ICOUNT+1).EQ.2) THEN
              ICOUNT = ICOUNT + 1
C     
C     ************************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
              IARRDUMP(11) = NINT(VALUE(ICOUNT))
            END IF
C     
C---- TOTAL total number of reflections to be dumped
C     
          ELSE IF (STR(1:4).EQ.'TOTA') THEN
C     
C---- test for a number N
C     
            IF (ICOUNT.LT.NTOK .AND. ITYP(ICOUNT+1).EQ.2) THEN
              ICOUNT = ICOUNT + 1
C     
C     ************************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
              MXDUMP = NINT(VALUE(ICOUNT))
            END IF
C     
C---- IMIN  Minimum intensity for a reflection to be dumped
C     
          ELSE IF (STR(1:4).EQ.'IMIN') THEN
C     
C---- test for a number N
C     
            IF (ICOUNT.LT.NTOK .AND. ITYP(ICOUNT+1).EQ.2) THEN
              ICOUNT = ICOUNT + 1
C     
C     ************************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
              IMINDUMP = NINT(VALUE(ICOUNT))
            END IF
C     
C---- IMAX  Maximum intensity for a reflection to be dumped
C     
          ELSE IF (STR(1:4).EQ.'IMAX') THEN
C     
C---- test for a number N
C     
            IF (ICOUNT.LT.NTOK .AND. ITYP(ICOUNT+1).EQ.2) THEN
              ICOUNT = ICOUNT + 1
C     
C     ************************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
              IMAXDUMP = NINT(VALUE(ICOUNT))
            END IF
C     
C     
C---- SPOT.... if debug is on for a particular subroutine and spot
C     is also set, then the optical densities of all spots will be
C     displayed when working online
C     
          ELSE IF (STR(1:4).EQ.'SPOT') THEN
            IF (ONLINE) WRITE (ITOUT,FMT=6064)
 6064       FORMAT (1X,'Pixel values of all spots will be dumped',
     +           ' for those subroutines which have DEBUG ',
     $           'turned on.')
            WRITE (IOUT,FMT=6064)
            SPOT = .TRUE.
C     
C---- HKL to dump specific reflections
C     
          ELSE IF (STR(1:3).EQ.'HKL') THEN
C     *******************************************
            CALL MKEYNM(3,ICOUNT+1,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            IF (NHKLD.LE.49) THEN
              NHKLD = NHKLD + 1
              DO 376  I = 1,3
                IHD(I,NHKLD) = NINT(VALUE(ICOUNT+I))
 376          CONTINUE
              ICOUNT = ICOUNT + 3
            END IF
            DUMPALL = .FALSE.
            NDTOT = NHKLD
C     
C---- RESOL  Resolution limts
C     
          ELSE IF (STR(1:4).EQ.'RESO') THEN
C     *******************************************
            CALL MKEYNM(2,ICOUNT+1,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            RESDLOW = VALUE(ICOUNT+1)
            RESD = VALUE(ICOUNT+2)
            IF (RESD.GT.RESDLOW) THEN
              TEMPX = RESD
              RESD = RESDLOW
              RESDLOW = TEMPX
            END IF
            ICOUNT = ICOUNT + 2
C
C---- IMAGE <number> Do not turn on DUMP until image specified has been reached
C
          ELSE IF (STR(1:4).EQ.'IMAG') THEN
            DUMP(12) = .TRUE.
            IARRDUMP(12) = 1
            IF (ICOUNT.LT.NTOK .AND. ITYP(ICOUNT+1).EQ.2) THEN
              ICOUNT = ICOUNT + 1
C     *********************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
C     
              IARRDUMP(12) = NINT(VALUE(ICOUNT))
            END IF
C
C---- Spots flagged as "EDGE" in output MTZ file
C
          ELSE IF (STR(1:7).EQ.'EDGEBAD') THEN
            DUMP(13) = .TRUE.
C
C---- Spots totally rejected because too many null pixels
C
          ELSE IF (STR(1:7).EQ.'EDGEREJ') THEN
            DUMP(14) = .TRUE.
C     
C---- NONE to turn off all dumping
C     
          ELSE IF (STR(1:4).EQ.'NONE') THEN
            WRITE(IOUT,FMT=7890)
            IF (ONLINE) WRITE(ITOUT,FMT=7890)
 7890       FORMAT(1X,'*** ALL DUMPING TURNED OFF ***')
            DUMPSPOT = .FALSE.
            SPOT = .FALSE.
            NHKLD = 0
            MXDUMP = 0
            DO 377 I = 1,30
              DUMP(I) = .FALSE.
 377        CONTINUE
C
          ELSE
C     
C---- Not recognised
C     
            IF (TRAPERR) INPERR = .TRUE.
            WRITE (IOUT,FMT=6130) STR(1:LENSTR(STR))
            IF (ONLINE) WRITE (ITOUT,FMT=6130) STR(1:LENSTR(STR))
          END IF
C     
C     
          IF (ICOUNT.LT.NTOK) GO TO 331
        ELSE          
C     
C     
          DO 410 I = 1,30
            DUMP(I) = .TRUE.
 410      CONTINUE
C     
C     
        END IF
C     
C     
        IF (DUMPSPOT) THEN
          WRITE (IOUT,FMT=6074)
 6074     FORMAT (/,1X,'Additional DEBUG output (voluminous) will be ',
     +       'given for:')
          IF (ONLINE) WRITE (ITOUT,FMT=6074)
        END IF
C     
C     
        DO 420 I = 1,11
          IF (DUMP(I)) THEN
            IF (I.EQ.4) THEN
              NCH = LENSTR(DUMPSTR(I))
              WRITE (IOUT,FMT=6075) DUMPSTR(I)(1:NCH),BGRLIM
 6075         FORMAT (1X,A,F6.2)
              IF (ONLINE) WRITE (ITOUT,FMT=6075) DUMPSTR(I)(1:NCH),
     +                    BGRLIM
            ELSE IF (I.GE.10) THEN
              NCH = LENSTR(DUMPSTR(I))
              WRITE (IOUT,FMT=6072) DUMPSTR(I)(1:NCH),IARRDUMP(I)
 6072         FORMAT (1X,A,I4)
              IF (ONLINE) WRITE (ITOUT,FMT=6072) DUMPSTR(I)(1:NCH),
     +                    IARRDUMP(I)
            ELSE
              WRITE (IOUT,FMT=6076) DUMPSTR(I)
 6076         FORMAT (1X,A)
              IF (ONLINE) WRITE (ITOUT,FMT=6076) DUMPSTR(I)
            END IF
C     
C     
            IF (I.EQ.1) THEN
              IF (ONLINE) WRITE(ITOUT,FMT=6077) NDUMP,MXDUMP,IMINDUMP,
     +                    IMAXDUMP
 6077         FORMAT(1X,'Every',I4,'th reflection will be dumped ,',
     +             ' up to a maximum of ',I5,' reflections',/,
     +             ' Reflection profiles will be dumped only ',
     +             'if the intensity is between',I7,' and',I7)
              WRITE(IOUT,FMT=6077) NDUMP,MXDUMP,IMINDUMP,IMAXDUMP
            END IF
C     
C     
          END IF
 420    CONTINUE
        IF (NHKLD.GT.0) THEN
          WRITE (IOUT,FMT=6046) ((IHD(I,J),I=1,3),J=1,NHKLD)
          IF(ONLINE)WRITE(ITOUT,FMT=6046)
     $           ((IHD(I,J),I=1,3),J=1,NHKLD)
        END IF
C     
C---- sdfac
C     
      ELSE IF (KEY.EQ.'SDFA') THEN
        WRITE (IOUT,FMT=6078)
 6078   FORMAT (/,1X,'***** This option is OBSOLETE ****',/,1X,
     +       'These parameters have been replaced by the SCAN',
     +       'NER INSTRUMENT FACTOR which is calculated by th',
     $       'e program')
        IF (ONLINE) WRITE (ITOUT,FMT=6078)
C     
C     
C---- Print options
C     
      ELSE IF (KEY.EQ.'PRIN') THEN
        OTHERS = .TRUE.
        PRINTL = .TRUE.
        IF (NTOK.EQ.1) THEN
          DO 440 I = 1,10
            LPRINT(I) = .TRUE.
 440      CONTINUE
        ELSE
C     
C     
          DO 450 I = 2,NTOK
            SUBKEY = LINE(IBEG(I) :IEND(I))
C     
C     ************
            CALL CCPUPC(SUBKEY)
C     ************
C     
            IF (SUBKEY.EQ.'FIDU') THEN
              LPRINT(1) = .TRUE.
            ELSE IF (SUBKEY.EQ.'REFI') THEN
              LPRINT(2) = .TRUE.
            ELSE IF (SUBKEY.EQ.'PROF') THEN
              LPRINT(3) = .TRUE.
            ELSE IF (SUBKEY.EQ.'STAT') THEN
              LPRINT(4) = .TRUE.
            ELSE
              IF (TRAPERR) INPERR = .TRUE.
              WRITE (IOUT,FMT=6090) SUBKEY
 6090         FORMAT (1X,'*** PRINT OPTION ',A,' Not recognised')
              IF (ONLINE) WRITE (ITOUT,FMT=6090) SUBKEY
            END IF   
 450      CONTINUE
        END IF
C     
C---- FOIL thicknesses
C     
      ELSE IF (KEY.EQ.'FOIL') THEN
        IF (IMGP) THEN
          WRITE(IOUT,6017) KEY
          IF (ONLINE) WRITE(ITOUT,6017) KEY
          BADLINE = .TRUE.
          GOTO 50
        END IF
        IF (NTOK.EQ.1) THEN
          THFOIL(2) = 0.15
          THFOIL(3) = 0.15
        ELSE IF (NTOK.EQ.2) THEN
C     
C     *******************************************
          CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          THFOIL(2) = VALUE(2)
          THFOIL(3) = THFOIL(2)
        ELSE
C     
C     *******************************************
          CALL MKEYNM(2,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ********************************************
C     
          THFOIL(2) = VALUE(2)
          THFOIL(3) = VALUE(3)
        END IF
C     
C     
        WRITE (IOUT,FMT=6094) (THFOIL(I),I=2,3)
 6094   FORMAT (/,1X,'The thickness of the FOILS between the "A" and ',
     +       '"B" films is',F5.2,'mm and',/,1X,
     +       'between "B" and "C" films is',F5.2,'mm')
        IF (ONLINE) WRITE (ITOUT,FMT=6094) (THFOIL(I),I=2,3)
C     
C---- convert to 10micron units
C     
        THFOIL(2) = THFOIL(2)*100.0
        THFOIL(3) = THFOIL(3)*100.0
C     
C---- SHIFT...shift measurment box so that the centre of the box is
C     no longer placed at the calculated centre of gravity
C     
      ELSE IF (KEY.EQ.'SHIF') THEN
C     
C     *******************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
        IF (NTOK.EQ.2) THEN
          IXSHIFT = NINT(VALUE(2))
        ELSE IF (NTOK.EQ.3) THEN
C     
C     *******************************************
          CALL MKEYNM(1,3,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          IXSHIFT = NINT(VALUE(2))
          IYSHIFT = NINT(VALUE(3))
        END IF
C     
C     
        IF (ONLINE) WRITE (ITOUT,FMT=6100) IXSHIFT,IYSHIFT
 6100   FORMAT (1X,'The measurement box will be displaced by',I2,' ',
     +       'PIXELS along scanner X and',I2,' pixels along scanner Y',
     +       /,1X,'from the calculated (CENTRE OF GRAVITY) position')
        WRITE (IOUT,FMT=6100) IXSHIFT,IYSHIFT
C     
C---- Findcc
C     
C     AL      ELSE IF (KEY.EQ.'FIND') THEN
C     AL        WRITE(IOUT,FMT=6450)
C     AL        IF (ONLINE) WRITE(ITOUT,FMT=6450)
 6450   FORMAT(1X,'*** This option not currently implemented ***')
C     AL        I = 0
C     AL        IF (I.EQ.0) GOTO 50
C     AL        IF (ONLINE) THEN
C     AL          FINDCC = .TRUE.
C     AL          WRITE (ITOUT,FMT=6102)
C     AL 6102 FORMAT (/,1X,'FILMPLOT will be called for all packs. after
C     refine',
C     AL     +       'ment using spots from the',/2X,'centre of the film
C     the R',
C     AL     +       'efined camera constants and distortion parameters'
C     ,/,1X,
C     AL     +       'will be written back to the GENERATE FILE',/,1X,'
C     The prog',
C     AL     +       'ram will then go on to the next PACK (only A films
C     will ',
C     AL     +       'be considered)',/,1X,'The camera constants and
C     distortion',
C     AL     +       ' parameters can be read back',/,1X,'From the
C     GENERATE FIL',
C     AL     +       'E for an "A" film if the keyword "READCC" is given'
C     
C     )
C     AL          WRITE (IOUT,FMT=6102)
C     AL        ELSE
C     AL          WRITE (IOUT,FMT=6104)
C     AL 6104 FORMAT (/,1X,'***** Cannot call FILMPLOT from a batch job
C     , card ',
C     AL     +       'ignored *****')
C     AL        END IF
C     
C---- Readcc
C     
      ELSE IF (KEY.EQ.'READ') THEN
        WRITE(IOUT,FMT=6450)
        IF (ONLINE) WRITE(ITOUT,FMT=6450)
        I = 0
        IF (I.EQ.0) THEN
          BADLINE = .TRUE.
          GOTO 50
        ENDIF
        READCC = .TRUE.
        IF (ONLINE) WRITE (ITOUT,FMT=6106)
 6106   FORMAT (/,1X,'Camera constants and distortion parameters ',
     +       'will be read from the generate file for "A" films')
        WRITE (IOUT,FMT=6106)
C     
C---- Precess
C     
      ELSE IF (KEY.EQ.'PREC') THEN
        IF (ONLINE) WRITE (ITOUT,FMT=6112)
 6112   FORMAT (/,1X,'Processing a PRECESSION PHOTOGRAPH',/,1X,
     +       'This will be treated as a film with NO FIDUCIALS')
        WRITE (IOUT,FMT=6112)
        PRECESS = .TRUE.
        NOFID = .TRUE.
        OMEGAFD = 0.0
        COSOM0 = COS(OMEGAFD)
        SINOM0 = SIN(OMEGAFD)
C     
C---- Nofid
C     
      ELSE IF (KEY.EQ.'NOFI') THEN
        OMEGAFD = 9999.0
C     
C---- alphanumeric, ie abc
C     
        IF (ITYP(2).EQ.1) THEN
          ABC = LINE(IBEG(2) :IEND(2))
C     
C     ***********
          CALL CCPUPC(ABC)
C     ***********
C     
          IF (ABC(1:1).EQ.'H') THEN
            OMEGAFD = 90.0
          ELSE IF (ABC(1:1).EQ.'V') THEN
            OMEGAFD = 0.0
          END IF
        END IF
        IF ((NTOK.EQ.1) .OR. (OMEGAFD.EQ.9999.0)) THEN
          IF (ONLINE) THEN
            WRITE (ITOUT,FMT=6118)
 6118       FORMAT (/,1X,'*** MUST Specify orientation of rotation ',
     +           'axis on scanner as "H" (HORIZONTAL) or ',
     $           '"V" (VERTICAL)')
            GO TO 50
          ELSE
            WRITE (IOUT,FMT=6118)
            STOP
          END IF
        END IF
        OMEGAF = OMEGAFD*DTOR

c        write(*, *) 'OMEGAF = ', omegaf

        OMEGA0 = OMEGAF
        COSOM0 = COS(OMEGAF)
        SINOM0 = SIN(OMEGAF)
        WRITE (IOUT,FMT=6120)
 6120   FORMAT (/,1X,'FIDUCIALS will not be used when processing ',
     $       'the film')
        IF (ONLINE) WRITE (ITOUT,FMT=6120)
        NOFID = .TRUE.
C     
C---- Interpolate
C     
      ELSE IF (KEY.EQ.'INTE') THEN
        INTERPOL = .TRUE.
C     
C---- @COMFILE.... read commands from file 'comfile.dat'
C     
      ELSE IF (KEY(1:1).EQ.'@') THEN
C     
C---- quoted token
C     
        IF (ITYP(2).EQ.3) THEN
          COMFILE = LINE(IBEG(1)+1:IEND(1))
        ELSE
          COMFILE = LINE(IBEG(1)+1:IEND(NTOK))
        END IF
C     
C---- Append .dat if not specified ***NO, NOT ANY MORE ****
C     find how many non-blank characters in comfile
C     
        NCH = LENSTR(COMFILE)
        WRITE (IOUT,FMT=6124) COMFILE(1:NCH)
 6124   FORMAT (1X,'Command filename: ',A)
        IF (ONLINE) WRITE (ITOUT,FMT=6124) COMFILE(1:NCH)
C     
C---- Test if command file exists
C     
        INQUIRE (FILE=COMFILE,EXIST=EFILE)
        IF (.NOT.EFILE) THEN
          NWRN = NWRN + 1
          WRITE (IOUT,FMT=6126) COMFILE(1:NCH)
 6126     FORMAT (//,1X,'**** ERROR ****',/,1X,'Command file ',A,
     +         ' DOES NOT EXIST')
          IF (ONLINE) THEN
            WRITE (ITOUT,FMT=6126) COMFILE(1:NCH)
            GO TO 50
          END IF
            NSHUTERR = 33
          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
        END IF
C     
C---- open command file
C     
        IFAIL = 1
C     
C     **********************************
        CALL CCPOPN(-ICOMM,COMFILE,3,1,80,IFAIL)
c     -vms        open (unit=icomm,file=comfile,status='old',readonly)
C     **********************************
C     
        COMREAD = .TRUE.
        ITINS = ITIN
        ITIN = ICOMM
        GO TO 60
C     
C---- Automatch
C     
      ELSE IF (KEY6.EQ.'AUTOMA') THEN
C     
C---- Need to differentiate AUTOMATCH and AUTOINDEX
C     
        
        OTHERS = .TRUE.
        MATCH = .TRUE.
C     
C---- Test for further keywords
C     
        IF (NTOK.GT.1) THEN
          IPNT = 1
 490      IPNT = IPNT + 1
          IF (IPNT.GT.NTOK) GO TO 500
          SUBKEY = LINE(IBEG(IPNT) :IEND(IPNT))
C     
C     ************
          CALL CCPUPC(SUBKEY)
C     ************
C     
          IF (SUBKEY.EQ.'RESO') THEN
            IPNT = IPNT + 1
C     
C     **********************************************
            CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     
            RESOL1 = VALUE(IPNT)
C     
C---- Test for second resolution limit
C     
            IF (IPNT.LT.NTOK .AND. ITYP(IPNT+1).EQ.2) THEN
              IPNT = IPNT + 1
C     
C     **********************************************
              CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     
              RESOL2 = VALUE(IPNT)
            END IF
C     
C     
            GO TO 490
          ELSE IF (SUBKEY.EQ.'RCON') THEN
            IPNT = IPNT + 1
C     
C     **********************************************
            CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     
            RCONV = VALUE(IPNT)
          ELSE IF (SUBKEY.EQ.'OVER') THEN
            IPNT = IPNT + 1
C     
C     **********************************************
            CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     
            OVRLAP = VALUE(IPNT)
          ELSE IF (SUBKEY.EQ.'NSTE') THEN
            IPNT = IPNT + 1
C     
C     **********************************************
            CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     
            NSTEP = NINT(VALUE(IPNT))
c     hrp 04012000            IF (NSTEP.GT.10) NSTEP = 10
          ELSE IF (SUBKEY.EQ.'NCYC') THEN
            IPNT = IPNT + 1
C     
C     **********************************************
            CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     
            NCYCA = NINT(VALUE(IPNT))
          ELSE IF (SUBKEY.EQ.'ANGL') THEN
            IPNT = IPNT + 1
C     
C     **********************************************
            CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     
            SECANGLE = VALUE(IPNT)
          ELSE IF (SUBKEY.EQ.'NOCE') THEN
            NOCENT = .TRUE.
          ELSE IF (SUBKEY.EQ.'NPAS') THEN
            IPNT = IPNT + 1
C     
C     **********************************************
            CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     
            NPASS = NINT(VALUE(IPNT))
          ELSE IF (SUBKEY.EQ.'DAMP') THEN
            IPNT = IPNT + 1
C     
C     **********************************************
            CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     
            DAMP = VALUE(IPNT)
          ELSE IF (SUBKEY.EQ.'RESI') THEN
            IPNT = IPNT + 1
C     
C     **********************************************
            CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     

            AXRMSLIM = VALUE(IPNT)
            ARRSET = .TRUE.
          ELSE IF (SUBKEY.EQ.'ELIM') THEN
            IPNT = IPNT + 1
C     
C     **********************************************
            CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     
            AELIMIT = 100.0*VALUE(IPNT)
          ELSE IF (SUBKEY.EQ.'CCOM') THEN
            IPNT = IPNT + 1
C     
C     **********************************************
            CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     
            TRUECCOM = VALUE(IPNT)
          ELSE IF (SUBKEY.EQ.'NOME') THEN
            NOMEAS = .TRUE.
          ELSE IF (SUBKEY.EQ.'MEAS') THEN
            NOMEAS = .FALSE.
          ELSE IF (SUBKEY.EQ.'MOSA') THEN
            RMOSAIC = .TRUE.
C     
C---- Test for maximum mosaic spread to try
C     
            IF (IPNT.LT.NTOK .AND. ITYP(IPNT+1).EQ.2) THEN
              IPNT = IPNT + 1
C     
C     **********************************************
              CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     
              ETAMAX = VALUE(IPNT)
            END IF
C     
C---- Check for turning OFF mosaic spread refinement
C     
            IF (IPNT.LT.NTOK) THEN
              KEY2 = LINE(IBEG(IPNT+1):IEND(IPNT+1))
              IF (KEY2(1:3).EQ.'OFF') RMOSAIC = .FALSE.
            END IF
          ELSE IF (SUBKEY.EQ.'BEAM') THEN
            IPNT = IPNT + 1
C     
C     **********************************************
            CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     
            NBEAM = NINT(VALUE(IPNT))
C     
C---- Test for suppression of orientation refinement, and also if we
C     want to turn it on again
C     
          ELSE IF (SUBKEY.EQ.'NORE') THEN
            NOREFINE = .TRUE.
            IF(IPNT.LT.NTOK)THEN
              IPNT = IPNT + 1
              SUBKEY = LINE(IBEG(IPNT) :IEND(IPNT))
              IF (SUBKEY(1:3).EQ.'OFF')THEN
                NOREFINE = .FALSE.
                IPNT = IPNT + 1
              ENDIF
            ENDIF
C     
C---- NREF  minimum acceptable number of reflections for positional
C     refinement in AUTOMATCH
C     
          ELSE IF (SUBKEY.EQ.'NREF') THEN
C     
            IPNT = IPNT + 1
C     *******************************************
            CALL MKEYNM(1,IPNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            MINREFAUTO = NINT(VALUE(IPNT))
C     
C---- Test for suppression of AUTOMATCH option
C     
          ELSE IF (SUBKEY(1:3).EQ.'OFF') THEN
            MATCH = .FALSE.
          ELSE        
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
            GO TO 500
          END IF
          GO TO 490
        END IF
 500    CONTINUE
C     
C     
C---- BEAM
C     Specify direct beam position for image plates or for use with 
C     NOFID option with film.
C     
      ELSE IF (KEY.EQ.'BEAM') THEN
        IBEAM = 2
         IF(HARRY)THEN
            IERR = 0
            CALL KWBEAM(LINE,MODE,COMREAD,ICOMM,IDFILM,IPROKWD,ITINS,
     $     OMEGAF,IERR)
            IF(IERR.EQ.50)THEN
              BADLINE = .TRUE.
              GOTO 50
            ENDIF
         ELSE
C IBEAM = 0, beam not set anywhere
C IBEAM = 1, set by the X-GUI
C IBEAM = 2, beam entered with BEAM keyword
C IBEAM = 3, read from image header
        IBEAM = 2
        IADD = 0
        IDFILM = 99999
        ISWUNG = 0
C     
C---- Test for OLD subkeyword for prelease 5.00 users !
C     
        IF ((NTOK.EQ.4).AND.(ITYP(4).EQ.1)) THEN
          SUBKEY = LINE(IBEG(4):IEND(4))
C     ***********
          CALL CCPUPC(SUBKEY)
C     ***********
        END IF
C     
C---- See if a pack or SWUNG_OUT subkeyword is present (will give an
C     even number of tokens)
C     
        IF (MOD(NTOK,2).EQ.0) THEN
C     
C---- Check if SWUNG_OUT given
C     
          IF (ITYP(2).EQ.1) THEN
            KEY2 = LINE(IBEG(2) :IEND(NTOK))
            CALL CCPUPC(KEY2)
            IF (KEY2.EQ.'SWUN') THEN
              ISWUNG = 1
            ELSE
              WRITE (IOUT,FMT=6130) KEY2
              IF (ONLINE) WRITE (ITOUT,FMT=6130) KEY2
            END IF
          ELSE
C     
C---- Get pack number
C     
C     *******************************************
            CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            IDFILM = NINT(VALUE(2))
          END IF
          IADD = 1
        END IF
C     
C---- Now get film centres..may be specified for
C     A, A and B, or A,B and C films.
C     
        
        IF (NTOK.EQ.3+IADD) THEN
C     
C     ************************************************
          CALL MKEYNM(2,2+IADD,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
          XMM(1) = VALUE(2+IADD)
          YMM(1) = VALUE(3+IADD)
          XMM(2) = XMM(1)
          YMM(2) = YMM(1)
          XMM(3) = XMM(1)
          YMM(3) = YMM(1)
        ELSE IF (NTOK.EQ.5+IADD) THEN
C     
C     ************************************************
          CALL MKEYNM(4,2+IADD,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
          XMM(1) = VALUE(2+IADD)
          YMM(1) = VALUE(3+IADD)
          XMM(2) = VALUE(4+IADD)
          YMM(2) = VALUE(5+IADD)
          XMM(3) = XMM(2)
          YMM(3) = YMM(2)
        ELSE IF (NTOK.EQ.7+IADD) THEN
C     
C     ************************************************
          CALL MKEYNM(6,2+IADD,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
          XMM(1) = VALUE(2+IADD)
          YMM(1) = VALUE(3+IADD)
          XMM(2) = VALUE(4+IADD)
          YMM(2) = VALUE(5+IADD)
          XMM(3) = VALUE(6+IADD)
          YMM(3) = VALUE(7+IADD)
        ELSE
          IF (ONLINE) THEN
            WRITE (ITOUT,FMT=6142)
 6142       FORMAT (//,1X,'*** ERROR IN Parameters on BEAM card ***')
            GO TO 50
          ELSE
            WRITE (IOUT,FMT=6142)
            NSHUTERR = 34
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF
        END IF
C     
C---- If no pack specified, set all film centres
C     
        IF (IDFILM.EQ.99999) THEN
          DO 520 I = 1,3
            DO 510 J = 1,MAXPAX
              XCENMM(J,I) = XMM(I)
              YCENMM(J,I) = YMM(I)
              XCENMMIN(J) = XMM(1)
              YCENMMIN(J) = YMM(1)
 510        CONTINUE
 520      CONTINUE
        ELSE
C     
C---- Find specified pack
C     
C---- First check that the PROCESS keyword has already been given
C     
          IF (IPROKWD.EQ.0) THEN
            WRITE(IOUT,FMT=7590)
            IF (ONLINE) WRITE(ITOUT,FMT=7590)
 7590       FORMAT(/,1X,'*** FATAL ERROR ***',/,1X,
     +           'If specifying different direct beam coordinates',
     +           ' for different images, the PROCESS',/,1X,'keyword',
     +           ' MUST be given before the BEAM keywords.')
            IF (ONLINE) THEN
              IF (COMREAD) THEN
                COMREAD = .FALSE.
                ITIN = ITINS
                CLOSE (UNIT=ICOMM)
              END IF
              BADLINE = .TRUE.
              GOTO 50
            END IF
            NSHUTERR = 35
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF
C     
          DO 540 J = 1,MAXPAX
            IF (IDFILM.EQ.IDPACK(J)) THEN
              DO 530 I = 1,3
                XCENMM(J,I) = XMM(I)
                YCENMM(J,I) = YMM(I)
 530          CONTINUE
              XCENMMIN(J) = XMM(1)
              YCENMMIN(J) = YMM(1)
              GO TO 550
            END IF
 540      CONTINUE
          IF (ONLINE) WRITE (ITOUT,FMT=6144) IDFILM
 6144     FORMAT (/,1X,'*** ERROR, Image',I7,' Has not been ',
     +         'specified on the PROCESS keyword ***')
          WRITE (IOUT,FMT=6144) IDFILM
          IF (ONLINE) THEN
            GO TO 50
          ELSE
            STOP
          END IF
        END IF
C     
C---- If input via window, need to define XCEN,YCEN and XCEN0, YCEN0
C     
 550    IF (MODE.EQ.3) THEN
          DO 551 I = 1,MAXPAX
            YCENMM(I,1) = YSCAL*YMM(1)
            IF ((TWOTHETA.NE.0.0).AND.(ISWUNG.EQ.0)) THEN
              XCENMM(I,1) = XCENMM(I,1) + 
     +             COS(OMEGAF)*0.01*XTOFD*TAN(TWOTHETA*DTOR)
              YCENMM(I,1) = YCENMM(I,1) + 
     +             SIN(OMEGAF)*YSCAL*0.01*XTOFD*TAN(TWOTHETA*DTOR)
            END IF
            IF (INVERTX) XCENMM(I,1) = NREC*RAST - XCENMM(I,1)
 551      CONTINUE
          XCEN = 100.0*XCENMM(1,1)
          YCEN = 100.0*YCENMM(1,1)
C          write(iout,fmt=9876) xcen
C          if (online) write(itout,fmt=9876) xcen
C 9876     format(1x,'xcen changed in control from xcenmm, now;',F8.2)
          XCEN0 = NINT(100.0*XMM(1) - 
     +         COS(OMEGAF)*XTOFD*TAN(TWOTHETA*DTOR))
          IF (INVERTX) XCEN0 = NINT(100.0*NREC*RAST-XCEN0)
          YCEN0 = NINT(100.0*YMM(1) - 
     +         SIN(OMEGAF)*XTOFD*TAN(TWOTHETA*DTOR))
C     
        END IF
        IF (.NOT.IMGP) THEN
          IF (IDFILM.EQ.99999) THEN
            WRITE (IOUT,FMT=6146) (XMM(I),YMM(I),I=1,3)
 6146       FORMAT (1X,'Direct beam coordinates for A,B, and C',
     +           ' films for all packs set to',3 (2F6.2,4X))
            IF (ONLINE) WRITE (ITOUT,FMT=6146) (XMM(I),YMM(I),I=1,3)
          ELSE
            WRITE (IOUT,FMT=6148) IDFILM, (XMM(I),YMM(I),I=1,3)
 6148       FORMAT (1X,'Direct beam coordinates for A,B, and C ',
     +           'films for pack ',I3,' set to',3 (2F6.2,4X))
            IF (ONLINE) WRITE (ITOUT,FMT=6148) IDFILM,
     +           (XMM(I),YMM(I),I=1,3)
          END IF
        END IF     
C----
C     end of "IF HARRY" block
      ENDIF
C----
C     
C     GAIN
C     
      ELSE IF (KEY.EQ.'GAIN') THEN
C     *******************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
        IGAIN = 1
        GAIN = VALUE(2)
C     
C     ADC OFFSET (IDIVIDE)
C     
      ELSE IF ((KEY.EQ.'ADCO').OR.(KEY.EQ.'OFFS')) THEN
C     *******************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
        IDIVIDE = NINT(VALUE(2))
        IIDIVIDE = 1
C     
C---- HELP library
C     
      ELSE IF (KEY.EQ.'HELP') THEN
        HLPMOS = LINE
C     
C     ***************
        CALL MOSHLP(HLPMOS)
C     ***************
C     
C     ORIEntation
C     
      ELSE IF (KEY.EQ.'ORIE') THEN
        SUBKEY = LINE(IBEG(2):IEND(2))
C     ***********
        CALL CCPUPC(SUBKEY)
C     ***********
C     
C     ORIEntation ROTAted
C     
        IF (SUBKEY.EQ.'ROTA') THEN
          ROTATED = .TRUE.
C     
C     ORIEntation STANdard- the default
C     
        ELSE IF (SUBKEY.EQ.'STAN') THEN
          ROTATED = .FALSE.
C     
C     ORIEntation ??????
C     
        ELSE
          IF (ONLINE) WRITE (ITOUT,FMT=6042) SUBKEY
          WRITE (IOUT,FMT=6042) SUBKEY
          IF (.NOT.ONLINE) STOP     
        END IF
C     
C     FILM characteristics (note that FLMPLOT is now invoked by PLOT)
C     sub keywords: ONEOD, SELWYN, BASEOD, NONLIN
C     
      ELSE IF (KEY.EQ.'FILM') THEN
        IF (NTOK.GT.1) THEN
          ICOUNT = 1
 580      ICOUNT = ICOUNT + 1
          SUBKEY = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     *********** 
          CALL CCPUPC(SUBKEY)
C     ***********
C     
          IF (SUBKEY.EQ.'BASE') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            BASEOD = VALUE(ICOUNT)
          ELSE IF (SUBKEY.EQ.'SELW') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            G1OD = VALUE(ICOUNT)
          ELSE IF (SUBKEY.EQ.'ONEO') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            N1OD = NINT(VALUE(ICOUNT))
          ELSE IF (SUBKEY.EQ.'NONL') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            CURV = VALUE(ICOUNT)
C     
C     Not recognised
C     
          ELSE
            IF (TRAPERR) INPERR = .TRUE.
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          END IF       
C     
C     
          IF (ICOUNT.LT.NTOK) GO TO 580
        END IF
C     
C     
C     SIZE image size, number of stripes of data and number of pixels
C     per stripe
      ELSE IF (KEY.EQ.'SIZE') THEN
        ISTOP = 0
        INSIZE = 1
C     *******************************************
        CALL MKEYNM(2,2,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
        NREC = NINT(VALUE(2))
        IYLEN = NINT(VALUE(3))
        IF (NREC.GT.IXWDTH/2) THEN
          WRITE(IOUT,FMT=7220) 2*NREC,IXWDTH,2*NREC
          IF (ONLINE) WRITE(ITOUT,FMT=7220) 2*NREC,IXWDTH,2*NREC
 7220     FORMAT(/,1X,'***** FATAL ERROR *****',/,1X,'The image ',
     +         'is too large. Change the parameter IXWDTH in the ',
     +         'code',/,1X,'to ',I6,' with a global edit and ',
     $         'recompile.',/,1X,
     +         'ie change all occurences of (IXWDTH=',I4,') to ',
     +         ' (IXWDTH=',I5,')')
          ISTOP = 1
        END IF
        IF (IYLEN.GT.IYLENGTH) THEN
          WRITE(IOUT,FMT=7222) 2*IYLEN,IYLENGTH,2*IYLEN
          IF (ONLINE) WRITE(ITOUT,FMT=7222) 2*IYLEN,IYLENGTH,2*IYLEN
 7222     FORMAT(/,1X,'***** FATAL ERROR *****',/,1X,'The image ',
     +         'is too large. Change the parameter IYLENGTH in the ',
     +         ' code',/,1X,'to ',I6,' with a global edit and ',
     $         'recompile.',/,1X,
     +         'ie change all occurences of (IYLENGTH=',I4,') to ',
     +         ' (IYLENGTH=',I5,')')
          ISTOP = 1
        END IF
        IF (ISTOP.NE.0) STOP
C     
C---- Check if number of header records given
C     
        IF (NTOK.GT.3) THEN
          SUBKEY = LINE(IBEG(4) :IEND(4))
          CALL CCPUPC(SUBKEY)
          IF (SUBKEY.EQ.'HEAD') THEN
C     ************************************************
            CALL MKEYNM(1,5,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            NHEAD = NINT(VALUE(5))
C     
C---- See if size of header is specified.
C     
            IF (NTOK.EQ.7) THEN
              SUBKEY = LINE(IBEG(6) :IEND(6))
              CALL CCPUPC(SUBKEY)
              IF (SUBKEY.EQ.'BYTE') THEN
C     ************************************************
                CALL MKEYNM(1,7,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
                NHBYTE = NINT(VALUE(7))
              ELSE
                WRITE (IOUT,FMT=6130) SUBKEY
                IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
              END IF
            END IF
C     
C     Not recognised
C     
          ELSE
            IF (TRAPERR) INPERR = .TRUE.
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          END IF
        END IF
C     
C     FIDUcials, subkeywords THRESHOLD (Od units), POSITION coords in mm
C     ,
C     SEARCH size of search box (mm), BEAM direct beam search size
C     
      ELSE IF (KEY.EQ.'FIDU') THEN
C     
        IF (NTOK.GT.1) THEN
          ICOUNT = 1
 585      ICOUNT = ICOUNT + 1
          SUBKEY = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     *********** 
          CALL CCPUPC(SUBKEY)
C     ***********
C     
          IF (SUBKEY.EQ.'SEAR') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
            XMMF = VALUE(ICOUNT)
          ELSE IF (SUBKEY.EQ.'THRE') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            THRESHF = VALUE(ICOUNT)
          ELSE IF (SUBKEY.EQ.'BEAM') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            XMMDB = VALUE(ICOUNT)
          ELSE IF (SUBKEY.EQ.'POSI') THEN
C     
C     Check next token is a number
C     
            NC = 0
 590        IF (ICOUNT.LT.NTOK .AND. ITYP(ICOUNT+1).EQ.2) THEN
              ICOUNT = ICOUNT + 1
C     
C     ************************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
              NC = NC + 1
              J = MOD(NC,2)
              IF (J.EQ.0) J = 2
              NFID = (NC-1)/2 + 1
              FIDXY(NFID,J) = VALUE(ICOUNT)
              GO TO 590
            ELSE   
              IF (ICOUNT.LT.NTOK) GOTO 585
            END IF          
          END IF
          IF (ICOUNT.LT.NTOK) GOTO 585
        END IF
C     
C     DISTortion parameters, sub keywords XTOFRA, YSCALE, TILT, TWIST,
C     ROFF,
C     TOFF RDROFF, RDTOFF, NODE (image plate)
      ELSE IF (KEY.EQ.'DIST') THEN
C     
C---- Ignore if this is a repeat of a multiseg post refinement (use
C     refined values instead)
C     
        IF (MULTISEG.AND.RPTFIRST) GOTO 50
C     
        IF (NTOK.GT.1) THEN
          ICOUNT = 1
 600      ICOUNT = ICOUNT + 1
          SUBKEY = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     *********** 
          CALL CCPUPC(SUBKEY)
C     ***********
C     
          IF (SUBKEY.EQ.'XTOF') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
C     
            XTOFRA = VALUE(ICOUNT)           
          ELSE IF (SUBKEY.EQ.'YSCA') THEN
            ICOUNT = ICOUNT + 1
            IYSCAL = 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
C     
C---- Check if YSCAL has already been set from ratio of pixel sizes
C     in slow and fast directions
C     
            IF (IPIXY.NE.0) THEN
              WRITE(IOUT,FMT=6464) YSCAL
              IF (ONLINE) WRITE(ITOUT,FMT=6464) YSCAL
              IF (BRIEF) WRITE(IBRIEF,FMT=6464) YSCAL
 6464         FORMAT(1X,'***** WARNING *****',/,1X,
     +             'YSCAL calculated',
     +             ' from the ratio of the pixel sizes in the slow',
     +             ' and ',/,1X,'fast directions (',F6.4,') will b',
     +             'e overwritten by the value given by DISTORTION',
     $             ' YSCAL keywords')
            END IF

            YSCAL = VALUE(ICOUNT)
            YSCALIN = YSCAL
          ELSE IF (SUBKEY.EQ.'TILT') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
C     
            ITILT = NINT(VALUE(ICOUNT))
            CALL SETDIS(ITILT,ITWIST,1)
c     RADEG = 18000.0/3.14159
c     IF (XTOFD.NE.0) FDIST = 1.0/ (XTOFD*RADEG)
c     TILT = ITILT*FDIST
          ELSE IF (SUBKEY.EQ.'TWIS') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            ITWIST = NINT(VALUE(ICOUNT))
            CALL SETDIS(ITILT,ITWIST,1)
c     RADEG = 18000.0/3.14159
c     IF (XTOFD.NE.0) FDIST = 1.0/ (XTOFD*RADEG)
c     TWIST = ITWIST*FDIST
C     
          ELSE IF (SUBKEY.EQ.'ROFF') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            ROFF = 100.0*VALUE(ICOUNT)
            IF (ABS(ROFF).GT.100.0) THEN
              WRITE(IOUT,FMT=6511) 0.01*ROFF
              IF (ONLINE) WRITE(ITOUT,FMT=6511) 0.01*ROFF
 6511         FORMAT(///,1X,'********************** WARNING ',
     $             '**************',
     +             /,1X,'Offsets must be given in mm, ',
     $             'NOT 10 micron units',
     +             /,1X,' Is the offset really',F7.2,'mm ?'//)
            END IF
C     
C---- Radially dependent radial offset
C     
          ELSE IF (SUBKEY.EQ.'RDRO') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            RDROFF = 100.0*VALUE(ICOUNT)
            IF (ABS(RDROFF).GT.100) THEN
              WRITE(IOUT,FMT=6511) 0.01*RDROFF
              IF (ONLINE) WRITE(ITOUT,FMT=6511) 0.01*RDROFF
            END IF
C     
          ELSE IF (SUBKEY.EQ.'BULG') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            IBULGE = NINT(VALUE(ICOUNT))
C     
C---- Offset (fixed tangential offset should be given in mm)
          ELSE IF (SUBKEY.EQ.'TOFF') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            TOFF = 100.0*VALUE(ICOUNT)
            IF (ABS(TOFF).GT.100.0) THEN
              WRITE(IOUT,FMT=6511) 0.01*TOFF
              IF (ONLINE) WRITE(ITOUT,FMT=6511) 0.01*TOFF
            END IF
C     
C---- Radially dependent tangential offset
C     
          ELSE IF (SUBKEY.EQ.'RDTO') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            RDTOFF = 100.0*VALUE(ICOUNT)
            IF (ABS(RDTOFF).GT.100.0) THEN
              WRITE(IOUT,FMT=6511) 0.01*RDTOFF
              IF (ONLINE) WRITE(ITOUT,FMT=6511) 0.01*RDTOFF
            END IF
C     
C---- NUmber of nodes in radially dependent ROFF, TOFF distortion.
C     
          ELSE IF (SUBKEY.EQ.'NODE') THEN
            INODES = 1
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            NODES = NINT(VALUE(ICOUNT)) - 1
            IF (NODES.EQ.0) THEN
              WRITE(IOUT,FMT=6515)
              IF (ONLINE) WRITE(ITOUT,FMT=6515)
 6515         FORMAT(1X,'*** WARNING *** A value of 1 for NODES',
     +             ' will have NO EFFECT, only values of 2',/,1X,
     +             'or more are sensible')
            END IF
C     
C---- Specify a given phase,as a multiple of pi/4. Thus sensible
C     values are 1,2,3 (0 is the default)
C     which are in fixed increments of pi/4.
C     
          ELSE IF (SUBKEY.EQ.'PHI') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            NPHI = VALUE(ICOUNT)
            IF ((NPHI.LT.0).OR.(NPHI.GT.3)) THEN
              WRITE(IOUT,FMT=6513) NPHI
              IF (ONLINE) WRITE(ITOUT,FMT=6513) NPHI
 6513         FORMAT(1X,'*** ERROR *** NPHI can only have values ',
     +             '0 to 3',/,1X,'Corresponding to starting ',
     $             'phases of NPHI*pi/4',/,1X,
     +             'NPHI will be reset to zero')
              NPHI = 0
            END IF
C     
C     Not recognised
C     
          ELSE
            IF (TRAPERR) INPERR = .TRUE.
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          END IF
          IF (ICOUNT.LT.NTOK) GOTO 600
        END IF
C     
C     CAMCON camera constants, subkeywords CCX,CCY,CCOMEGA,CBAR (mm)
      ELSE IF (KEY.EQ.'CAMC') THEN
C     
C---- Ignore if this is a repeat of a multiseg post refinement (use
C     refined values instead)
C     
        IF (MULTISEG.AND.RPTFIRST) GOTO 50
C     
        IF (NTOK.GT.1) THEN  
          ICOUNT = 1
 610      ICOUNT = ICOUNT + 1
          SUBKEY = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     *********** 
          CALL CCPUPC(SUBKEY)
C     ***********
C     
          IF (SUBKEY(1:3).EQ.'CCX') THEN
C     
C---- Discourage use of camera constants, unless this is an iMosflm job
C     in which case camera constants have to be reset when detector 
C     parameters are reset
C     
            IF (.NOT.SOCKLO) THEN
              WRITE(IOUT,FMT=7240)
              IF (ONLINE) WRITE(ITOUT,FMT=7240)
            END IF
 7240       FORMAT(/,1X,'***** WARNING *****',/,1X,
     +           /,1X,'***** WARNING *****',
     +           /,1X,'***** WARNING *****',
     +           /,1X,'***** WARNING *****',
     +           'The use',
     +           ' of camera constants CCX and CCY is to be avoi',
     +           'ded.',/,1X,'Please input the CORRECT direct b',
     +           'eam coordinates and do not rely',/,1X,'on usin',
     +           'g CCX and CCY.',/,/,1X,'In addition, do NOT us',
     +           'e CCOMEGA to allow for non-standard scanner',/,
     +           1X,'orientations (eg vertical rotation axis for',
     +           ' Mar). ',/,1X,'Use the SCANNER OMEGA keyword f',
     +           'or this. (Additional information in help libra',
     $           'ry).')
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
C     
            XCC = VALUE(ICOUNT)
            IF (ABS(XCC).GT.3.0) THEN
              WRITE(IOUT,FMT=6510) XCC
              IF (ONLINE) WRITE(ITOUT,FMT=6510) XCC
 6510         FORMAT(///,1X,
     $             '********************** WARNING **************',
     +             /,1X,'Camera constants must be in mm, ',
     $             'NOT 10 micron units',/,1X,
     +             ' Is the camera constant really',F6.1,'mm ?'//)
            END IF
            CCX = NINT(100.0*XCC)
            ICCX = 1
          ELSE IF (SUBKEY(1:3).EQ.'CCY') THEN
            IF (.NOT.SOCKLO) THEN
              WRITE(IOUT,FMT=7240)
              IF (ONLINE) WRITE(ITOUT,FMT=7240)
            END IF
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            XCC = VALUE(ICOUNT)
            IF (ABS(XCC).GT.3.0) THEN
              WRITE(IOUT,FMT=6510) XCC
              IF (ONLINE) WRITE(ITOUT,FMT=6510) XCC
            END IF
            CCY = NINT(100.0*XCC)
            ICCY = 1
          ELSE IF (SUBKEY.EQ.'CCOM') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************        
C     
            CCOM = VALUE(ICOUNT)
          ELSE IF (SUBKEY.EQ.'CBAR') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            CBAR = NINT(100.0*VALUE(ICOUNT))
C     
C     
C     Not recognised
C     
          ELSE
            IF (TRAPERR) INPERR = .TRUE.
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          END IF
C     
          IF (ICOUNT.LT.NTOK) GOTO 610
        END IF
C     
C---- OVERLOADS, subkeywords NOVER, CUTOFF. Defines "overloaded"
C     spots as those with more than NOVPIX pixels with values above
C     CUTOFF                       
C     
      ELSE IF (KEY.EQ.'OVER') THEN
        IF (NTOK.GT.1) THEN
          ICOUNT = 1
 620      ICOUNT = ICOUNT + 1
          SUBKEY = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     *********** 
          CALL CCPUPC(SUBKEY)
C     ***********
C---- CUTOFF....od gt cutoff will be treated as overload
C     
          IF (SUBKEY.EQ.'CUTO') THEN
C     
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            CUTOFF = NINT(VALUE(ICOUNT))
C     
C---- for safety !!
C     
            CUTOFF = CUTOFF - 1  
            ICUT = 2
C     
C     NOVER...sets NOVPIX, maximum number of pixels with od greater than
C     cutoff for a spot to be flagged as an overload
C     
          ELSE IF (SUBKEY.EQ.'NOVE') THEN          
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            NOVPIX=NINT(VALUE(ICOUNT))
C     
C---- Not recognised
C     
          ELSE
            IF (TRAPERR) INPERR = .TRUE.
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          END IF
C     
          IF (ICOUNT.LT.NTOK) GOTO 620
        END IF
C     
C---- ADDPartials over adjacent images
C     
      ELSE IF (KEY.EQ.'ADDP') THEN
C
C     This option is deprecated and the keyword will generate an error
C     message (and not be applied) unless a subkeyword REALLY is given.
C
        IF ((NTOK.EQ.1).OR.((NTOK.EQ.2).AND.
     +    (LINE(IBEG(2):IEND(2)).NE.'REAL'))) THEN
          WRITE(IOUT,FMT=7720)
          IF (ONLINE) WRITE(ITOUT,FMT=7720)
 7720     FORMAT(/,/,1X,'*** This keyword should not be used and has',
     +           ' been ignored ***',/,1X,
     +      'If you really need to use it, send an Email to ',
     +      'andrew@mrc-lmb.cam.ac.uk')
          GOTO 50
        END IF
C
        IF (.NOT.IMGP) THEN
          WRITE(IOUT,6019) KEY
          IF (ONLINE) WRITE(ITOUT,6019) KEY
 6019     FORMAT(/,1X,'**** ERROR ***** ',A,' not appropriate for',
     +         ' film data, ignored')
          BADLINE = .TRUE.
          GOTO 50
        END IF
        ADDPART = .TRUE.
        SUMPART = .TRUE.
        SADDPART = .TRUE.
        SSUMPART = .TRUE.
C     
C---- Check for turning ADDPART OFF
C     
        IF (NTOK.GT.1) THEN
          SUBKEY = LINE(IBEG(2):IEND(2))
C     ***********
          CALL CCPUPC(SUBKEY)
C     ***********
          IF (SUBKEY(1:3).EQ.'OFF') THEN
            ADDPART = .FALSE.
            SADDPART = .FALSE.
          ELSE
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          END IF
        END IF
C     
C---- THICKNESS... Effective thickness of active area of detector.
C     This affects the raster expansion.
C     
      ELSE IF (KEY.EQ.'THIC') THEN
C     
C---- Redundant keyword....only use if MODIFY is also specified
C     
        SUBKEY = LINE(IBEG(2):IEND(2))
C     ***********
        CALL CCPUPC(SUBKEY)
C     ***********
        IF (SUBKEY.NE.'MODI') THEN
          WRITE(IOUT,FMT=6610)
          IF (ONLINE) WRITE(ITOUT,FMT=6610)
          IF (BRIEF) WRITE(IBRIEF,FMT=6610)
          BADLINE = .TRUE.
          GOTO 50
        END IF
C     *******************************************
        CALL MKEYNM(1,3,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
        THICK = VALUE(3)
C     
C     Convert to 10 micron units
C     
        THICK = THICK * 100.0
C     
C     DIRECTORY... Logical name for directory containing images, can be
C     up to
C     ten of these
C     
      ELSE IF (KEY.EQ.'DIRE') THEN
        NNDIR = NTOK - 1
        IF(SOCKLO)NNDIR = 1
        DO 624 I = 1,NNDIR
          CKNDIR = .FALSE.
          NDIR = NDIR + 1
          IF (NDIR.GT.10) GOTO 624
C
C---- Check against any previously stored directories
C
          DO 623 J=1,NDIR-1
            IF(FDISK(J)(1:1+(IEND(I+1)-IBEG(I+1)))
     $           .EQ.LINE(IBEG(I+1):IEND(I+1)))CKNDIR = .TRUE.
 623      ENDDO
          IF(CKNDIR)THEN
            NDIR = NDIR - 1
          ELSE
            IF(SOCKLO)THEN
              FDISK(NDIR) = LINE(IBEG(I+1) :IEND(NTOK))
            ELSE
              FDISK(NDIR) = LINE(IBEG(I+1) :IEND(I+1))
            ENDIF
          ENDIF
C     
C---- If running under Unix, check for a "/" at end of directory
C     specification, and add one if not present.
C     
          IF (.NOT.VAXVMS()) THEN
            NCH = LENSTR(FDISK(NDIR))
            IF (FDISK(NDIR)(NCH:NCH).NE.'/') THEN
              FDISK(NDIR) = FDISK(NDIR)(1:NCH)//'/'
            END IF
          END IF
 624    CONTINUE
C     
C     EXTENSION.... Extension for image filenames
C     
      ELSE IF (KEY.EQ.'EXTE') THEN
        ODEXT = ' '
        IF (NTOK.GT.1) THEN
          IF (ITYP(2).EQ.3) THEN
            ODEXT = LINE(IBEG(2) :IEND(2))
          ELSE
            ODEXT = LINE(IBEG(2) :IEND(NTOK))
          END IF        
          IEXTEN = 1
        ELSE
          WRITE(IOUT,FMT=7450)
          IF (ONLINE) WRITE(ITOUT,FMT=7450)
 7450     FORMAT(1X,'***** ERROR *****',/,1X,'No extension given.')
        END IF
        NCHAR = LENSTR(ODEXT)
C     
C---- only set PACK true if extension is "pck" and not "pck1200" etc
C     
        PACK = (((ODEXT(1:3).EQ.'PCK').OR.(ODEXT(1:3).EQ.'pck')).AND.
     +       (NCHAR.EQ.3))
        IF (PACK) USEHDR = .FALSE.
C     
C     BACKGROUND Subkeywords: BGSIG BGFRAC RECOVER
C---- BGSIG
C     points in background area which deviate by more than
C     bgsig*rmsbg from the background plane will be eliminated,
C     and the background plane recalculated.
C     BGFRAC                                 
C     Only a fraction BGFRAC of the background pixels will be used
C     in the initial determination of the background plane (def 0.8)
C     RECOVER RECLEVEL
C     Get additional background pixels if there are less than 
C     RECLEVEL*NBGMIN
C     
      ELSE IF (KEY6.EQ.'BACKGR') THEN
        IF (NTOK.GT.1) THEN
          ICOUNT = 1
 630      ICOUNT = ICOUNT + 1
          SUBKEY = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     *********** 
          CALL CCPUPC(SUBKEY)
C     ***********       
C     
          IF (SUBKEY.EQ.'BGSI') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            IBGSIG = 1
            BGSIG = VALUE(ICOUNT)
          ELSE IF (SUBKEY.EQ.'BGFR') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C
            IBGFRAC = 1     
            BGFRAC = VALUE(ICOUNT)
          ELSE IF (SUBKEY.EQ.'BGAD') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            BKGADD = VALUE(ICOUNT)
          ELSE IF (SUBKEY.EQ.'RECO') THEN
            RECOVER = .TRUE.
C     
C---- Check for a number setting recovery level as a multiple
C     of the minimum number of background pixels
C     
            IF (((ICOUNT+1).LE.NTOK).AND.(ITYP(ICOUNT+1).EQ.2)) THEN
              ICOUNT = ICOUNT + 1
C     *********************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
              RECLEVEL = VALUE(ICOUNT)
            END IF
          ELSE
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          END IF
C     
          IF (ICOUNT.LT.NTOK) GOTO 630
        END IF
C     
C     REFINEMENT Subkeywords: GRADMAXR BGREJ USEBOX FIX FREE etc
C---- Spots with gradient/background greater than GRADMAXR will be
C     rejected from
C     refinement.
C     BGREJR... spots where more than a fraction BGREJR of the total
C     number of background pixels have been rejected will be rejected
C     from refinement.
C     USEBOX... Use measurement box in refinement of central region 
C     of image
C     MAXREF    Maximum number of reflections to use
C     NCENT     Number of reflections from central area to use when
C               refining with spots over entire image
C     NSM1      Number of images to be used in smoothing variable
C               parameters CCOM,YSCALE,DIST
C     NSM2      Number of images to be used in smoothing "constant"
C               parameters TILT,TWIST,ROFF,TOFF
C     SMOOTH    Turn on parameter smoothing
C     NOSMOOTH  Turn off parameter smoothing
C     INCLUDE   Include either or both partials or/and overloads
C     EXCLUDE   Turn off INCLUDE
C
      ELSE IF (KEY.EQ.'REFI') THEN
        IF (NTOK.GT.1) THEN
          ICOUNT = 1
 640      ICOUNT = ICOUNT + 1
 641      SUBKEY = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     *********** 
          CALL CCPUPC(SUBKEY)
C     ***********
C---- Weighted refinement
C     
          IF (SUBKEY.EQ.'WEIG') THEN
            RWEIGHT = .TRUE.
          ELSE IF (SUBKEY.EQ.'NOWE') THEN
            RWEIGHT = .FALSE.
C
C---- SUMIMAGES .. Number of images to add for centroid determination
C
          ELSE IF (SUBKEY.EQ.'SUMI') THEN   
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            NIMGSUMREF = VALUE(ICOUNT)
            NIMGSUMREF = MAX(NIMGSUMREF,1)
C     
C---- REJECTION    Rejection limit in refinement...
C                  if del.gt.REFREJ*SIGMA the reflection is rejected
C
          ELSE IF (SUBKEY.EQ.'REJE') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            REFREJ = VALUE(ICOUNT)
C     
C---- USEBOX/NOUSEBOX
C     
          ELSE IF (SUBKEY.EQ.'USEB') THEN
            USEBOX = .TRUE.
          ELSE IF (SUBKEY.EQ.'NOUS') THEN
            USEBOX = .FALSE.
C     
C---- BGREJECT (Max allowed fraction of background points rejected)
C     
          ELSE IF (SUBKEY.EQ.'BGRE') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            BGFREJ = VALUE(ICOUNT)
C     
C---- GRADIENT
C     
          ELSE IF (SUBKEY.EQ.'GRAD') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            GRADMAXR = VALUE(ICOUNT)
            IGRADMAXR = 1
C     
C---- ISDR...sets int/sd ratio for selection of refinement
C     spots in next
C     
          ELSE IF (SUBKEY.EQ.'ISDR') THEN
C     
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NSDR = NINT(VALUE(ICOUNT))
            IF (IMGP) THEN
              WRITE(IOUT,6017) SUBKEY
              IF (ONLINE) WRITE(ITOUT,6017) SUBKEY
              IF (ICOUNT.LT.NTOK) GOTO 640
            END IF
C     
C---- NREF  minimum acceptable number of reflections for positional
C     refinement
C     
          ELSE IF (SUBKEY.EQ.'NREF') THEN
C     
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            MINREF = NINT(VALUE(ICOUNT))
C     
C---- MAXREF  maximum number of reflections to be used for positional
C     refinement
C     
          ELSE IF (SUBKEY.EQ.'MAXR') THEN
C     
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            MAXREF = NINT(VALUE(ICOUNT))
C     
C---- NCENT  Number of reflections from central area to carry over
C            for refinement over entire image.
C     
          ELSE IF (SUBKEY.EQ.'NCEN') THEN
C     
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NCENTREF = NINT(VALUE(ICOUNT))
C     
C---- NSM1      Number of degrees of data to be used in smoothing variable
C               parameters CCOM,YSCALE,DIST
C     
          ELSE IF (SUBKEY.EQ.'NSM1') THEN
C     
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NSM1 = NINT(VALUE(ICOUNT))
C     
C---- NSM2      Number of degrees to be used in smoothing "constant"
C               parameters  TILT,TWIST,ROFF,TOFF
C     
          ELSE IF (SUBKEY.EQ.'NSM2') THEN
C     
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NSM2 = NINT(VALUE(ICOUNT))
C     
C---- SMOOTH  Turn on parameter smoothing 
C     
          ELSE IF (SUBKEY.EQ.'SMOO') THEN
C     
          NOSMOOTH = .FALSE.
C     
C---- NOSMOOTH  Turn off parameter smoothing 
C     
          ELSE IF (SUBKEY.EQ.'NOSM') THEN
C     
          NOSMOOTH = .TRUE.
C     
C---- DAMP  Turn on parameter shift damping (TILT and TWIST only)
C     
          ELSE IF (SUBKEY.EQ.'DAMP') THEN
C     
          DAMPREF = .TRUE.
C     
C---- NODAMP  Turn off parameter damping (TILT and TWIST only)
C     
          ELSE IF (SUBKEY.EQ.'NODA') THEN
C     
          DAMPREF = .FALSE.
C     
C---- RESID....sets maximum allowed residual following refinement
C     of central region (centrs) and over the whole film
C     
          ELSE IF (SUBKEY.EQ.'RESI') THEN
C     
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            XRMSLIM = VALUE(ICOUNT)
            RRSET = .TRUE.
C     
C---- CYCLES...number of refinement cycles in centrs
C     
          ELSE IF (SUBKEY.EQ.'CYCL') THEN
C     
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NCYC = NINT(VALUE(ICOUNT))
C     
C---- NSIG....rejection criterion for weak reflections in centrs
C     
          ELSE IF (SUBKEY.EQ.'NSIG') THEN
C     
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NSIG = NINT(VALUE(ICOUNT))
            NSIGSAVE = NSIG
            ISIGSET = 1
C     
C---- IMIN...rejection criterion for weak reflections in seekrs
C     
          ELSE IF (SUBKEY.EQ.'IMIN') THEN
C     
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            IRFMIN = NINT(VALUE(ICOUNT))
C     
C---- Test for second number (increment)
C     
            IF (ICOUNT.LT.NTOK .AND. ITYP(ICOUNT+1).EQ.2) THEN
              ICOUNT = ICOUNT + 1
C     
C     **********************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     **********************************************
C     
              IRFINC = NINT(VALUE(ICOUNT))
            ELSE 
              IRFINC = IRFMIN/2
            END IF
C     
C     
C     
C---- LIMIT   (used in centrs)
C     
          ELSE IF (SUBKEY.EQ.'LIMI') THEN
            ICOUNT = ICOUNT + 1
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            XLIMIT = VALUE(ICOUNT)
            LIMIT = 100*NINT(XLIMIT)
            ILIMIT = 1
C     
C---- VLIMIT   (used in centrs for vee films)
C     
          ELSE IF (SUBKEY.EQ.'VLIM') THEN
            ICOUNT = ICOUNT + 1
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
            VLIM = VALUE(ICOUNT)
            WRITE (IOUT,FMT=6116) VLIM
 6116       FORMAT (/,1X,'Maximum X coordinate for spots to be used',
     +           ' in refinement of centre of film set to',
     $           F5.1,' mm')
            IF (ONLINE) WRITE (ITOUT,FMT=6116) VLIM
            VLIM = 100*VLIM
C     
C---- FULLFRAC include partials if fraction of fully recorded reflection
C     is less
C     than FULLFRAC
C     
          ELSE IF (SUBKEY.EQ.'FULL') THEN
            ICOUNT = ICOUNT + 1
C     ************************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
            FULLFRAC = VALUE(ICOUNT)
C     
C---- INCLUDE partials or overloads in refinement
C     
          ELSE IF (SUBKEY.EQ.'INCL') THEN
 560        ICOUNT = ICOUNT + 1
            IF (ICOUNT.GT.NTOK) GO TO 570
            STR = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     ***********
            CALL CCPUPC(STR)
C     ***********
C     
            IF (STR(1:4).EQ.'PART') THEN
              USEPAR = .TRUE.
C     
C---- Test for a minimum partiality
C     
              IF (ICOUNT.LT.NTOK .AND. ITYP(ICOUNT+1).EQ.2) THEN
                ICOUNT = ICOUNT + 1
C     *********************************************
                CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
                PTMIN = VALUE(ICOUNT)
                IPTMINSET = 1
              END IF
            ELSE IF (STR(1:4).EQ.'OVER') THEN
              USEOVR = .TRUE.
            ELSE
C     
C---- If USEPAR or USEOVR set, this may be another keyword, so go
C     back to top of loop to check
C     
              IF (USEPAR.OR.USEOVR) GOTO 570
              WRITE (IOUT,FMT=6130) STR
              IF (ONLINE) WRITE (ITOUT,FMT=6130) STR
              GOTO 640
            END IF
            GO TO 560
 570        IF (ICOUNT.LE.NTOK) GOTO 641
C     
C---- EXCLUDE partials or overloads in refinement
C     
          ELSE IF (SUBKEY.EQ.'EXCL') THEN
 5601       ICOUNT = ICOUNT + 1
            IF (ICOUNT.GT.NTOK) GO TO 5701
            STR = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     ***********
            CALL CCPUPC(STR)
C     ***********
C     
            IF (STR(1:4).EQ.'PART') THEN
              USEPAR = .FALSE.
            ELSE IF (STR(1:4).EQ.'OVER') THEN
              USEOVR = .FALSE.
            ELSE
C     
C---- If USEPAR or USEOVR set, this may be another keyword, so go
C     back to top of loop to check
C     
              IF (.NOT.USEPAR.OR..NOT.USEOVR) GOTO 5701
              WRITE (IOUT,FMT=6130) STR
              IF (ONLINE) WRITE (ITOUT,FMT=6130) STR
              GOTO 640
            END IF
            GO TO 5601
 5701       IF (ICOUNT.LE.NTOK) GOTO 641
C     
C---- FREE parameters. (NOTE that all parameters except RDROFF, RDTOFF
C     are free by default).
C     
          ELSE IF (SUBKEY.EQ.'FREE') THEN
 562        ICOUNT = ICOUNT + 1
            IF (ICOUNT.GT.NTOK) GO TO 578
            STR = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     ***********
            CALL CCPUPC(STR)
C     ***********
C     
            IF (STR(1:4).EQ.'XCEN') THEN
              IFIX(1) = 0
              FIXPAR(1) = .FALSE.
            ELSE IF (STR(1:4).EQ.'YCEN') THEN
              IFIX(2) = 0
              FIXPAR(2) = .FALSE.
            ELSE IF((STR(1:4).EQ.'OMEG').OR.(STR(1:4).EQ.'CCOM'))THEN
              IFIX(3) = 0
              FIXPAR(3) = .FALSE.
            ELSE IF (STR(1:4).EQ.'YSCA') THEN
              IFIX(4) = 0
              FIXPAR(4) = .FALSE.
            ELSE IF ((STR(1:4).EQ.'XTOF')
     +             .OR.(STR(1:4).EQ.'DIST')) THEN
              IFIX(5) = 0
              FIXPAR(5) = .FALSE.
            ELSE IF (STR(1:4).EQ.'TILT') THEN
              IFIX(6) = 0
              FIXPAR(6) = .FALSE.
            ELSE IF (STR(1:4).EQ.'TWIS') THEN
              IFIX(7) = 0
              FIXPAR(7) = .FALSE.
            ELSE IF (STR(1:4).EQ.'ROFF'.OR.STR(1:4).EQ.'BULG') THEN
              IFIX(8) = 0
              FIXPAR(8) = .FALSE.
            ELSE IF (STR(1:4).EQ.'TOFF') THEN
              IFIX(9) = 0
              FIXPAR(9) = .FALSE.
            ELSE IF (STR(1:4).EQ.'RDTO') THEN
              FIXPAR(10) = .FALSE.
              IFIX(10) = 0
            ELSE IF (STR(1:4).EQ.'RDRO') THEN
              FIXPAR(11) = .FALSE.
              IFIX(11) = 0
            ELSE
C     
C---- This may be another subkeyword, go back to top of list to check
C     
              GOTO 641
            END IF
C     
C---- Get next subkeyword
C     
            GOTO 562
C     
C     
C---- FIX parameters. Can fix:
C     1 XCEN
C     2 YCEN
C     3 OMEGA0
C     4 YSCAL
C     5 XTOFRA or DIST  (Crystal to detector distance multiplier)
C     6 TILT            
C     7 TWIST
C     8 ROFF for Image Plate, BULGE for film
C     9 TOFF
C     10 RDTOFF
C     11 RDROFF
C     
          ELSE IF (SUBKEY.EQ.'FIX') THEN
 572        ICOUNT = ICOUNT + 1
            IF (ICOUNT.GT.NTOK) GO TO 578
            STR = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     ***********
            CALL CCPUPC(STR)
C     ***********
C     
            IF (STR(1:4).EQ.'XCEN') THEN
              IFIX(1) = 1
              FIXPAR(1) = .TRUE.
            ELSE IF (STR(1:4).EQ.'YCEN') THEN
              IFIX(2) = 1
              FIXPAR(2) = .TRUE.
            ELSE IF((STR(1:4).EQ.'OMEG').OR.(STR(1:4).EQ.'CCOM'))THEN
              IFIX(3) = 1
              FIXPAR(3) = .TRUE.
            ELSE IF (STR(1:4).EQ.'YSCA') THEN
              IFIX(4) = 1
              FIXPAR(4) = .TRUE.
            ELSE IF ((STR(1:4).EQ.'XTOF')
     +             .OR.(STR(1:4).EQ.'DIST')) THEN
              IFIX(5) = 1
              FIXPAR(5) = .TRUE.
            ELSE IF (STR(1:4).EQ.'TILT') THEN
              IFIX(6) = 1
              FIXPAR(6) = .TRUE.
            ELSE IF (STR(1:4).EQ.'TWIS') THEN
              IFIX(7) = 1
              FIXPAR(7) = .TRUE.
            ELSE IF (STR(1:4).EQ.'ROFF'.OR.STR(1:4).EQ.'BULG') THEN
              IFIX(8) = 1
              FIXPAR(8) = .TRUE.
            ELSE IF (STR(1:4).EQ.'TOFF') THEN
              IFIX(9) = 1
              FIXPAR(9) = .TRUE.
            ELSE IF (STR(1:4).EQ.'RDTO') THEN
              IFIX(10) = 1
              FIXPAR(10) = .TRUE.
            ELSE IF (STR(1:4).EQ.'RDRO') THEN
              IFIX(11) = 1
              FIXPAR(11) = .TRUE.
            ELSE
C     
C---- This may be another subkeyword, go back to top of list to check
C     
              GOTO 641
            END IF
C     
C---- Get next subkeyword
C     
            GOTO 572
C     
C---- Not recognised
C     
          ELSE
            IF (TRAPERR) INPERR = .TRUE.
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          END IF
C     
 578      IF (ICOUNT.LT.NTOK) GOTO 640
        END IF
C     
C     REJECTION Subkeywords: GRADMAX BGRATIO PKRATIO MINB PLOT DUMP
C---- Spots with gradient/background greater than GRADMAX, or with a
C     background
C     ratio greater than BGRATIO or a profile fit ratio greater than
C     PKRATIO or with fewer than NBGMIN background pixels after 
C     background point rejection will be flagged as BADSPOTS
C
      ELSE IF (KEY.EQ.'REJE') THEN
        IF (NTOK.GT.1) THEN
          ICOUNT = 1
 650      ICOUNT = ICOUNT + 1
          SUBKEY = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     *********** 
 651      CALL CCPUPC(SUBKEY)
C     ***********
C     
          IF (SUBKEY.EQ.'BGRA') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            BGRAT = VALUE(ICOUNT)
          ELSE IF (SUBKEY.EQ.'PKRA') THEN
C     
C---- PKRATIO... allowed combinations are:
C     PKRATIO 3.5
C     PKRATI0 3.5 ACCEPT
C     PKRATIO ACCEPT
C     
C     where ACCEPT signifies that reflections failing the
C     PKRATIO test should have the profile fitted intensity
C     set to the integrated intensity
C     
C---- Check if next token is a number
C     
            IF (ITYP(ICOUNT+1).EQ.2) THEN
              ICOUNT = ICOUNT + 1
C     *******************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
              PKRAT = VALUE(ICOUNT)
              IPKRAT = 1
C     
C---- Test if next token qualifies PKRATIO test, only allowed qualifying
C     keyword is ACCEPT
C     
              IF (ICOUNT.LT.NTOK) THEN
                ICOUNT = ICOUNT + 1
                SUBKEY = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     *********** 
                CALL CCPUPC(SUBKEY)
C     ***********
                IF (SUBKEY.EQ.'ACCE') THEN
Chrp25.11.2004 one of Marian's changes
                  PKACCEPT = .TRUE.
                ELSE
C     
C---- NOT a PKRATIO qualifier, treat as regular subkeyword
C     
                  GOTO 651
                END IF
              END IF
            ELSE
C     
C---- NOT a number following PKRATIO, test for qualifier ACCEPT
C     
              ICOUNT = ICOUNT + 1
              SUBKEY = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     *********** 
              CALL CCPUPC(SUBKEY)
C     ***********
              IF (SUBKEY.EQ.'ACCE') THEN
                PKACCEPT = .TRUE.
              ELSE
C     
C---- Not recognised
                IF (TRAPERR) INPERR = .TRUE.
                WRITE (IOUT,FMT=6130) SUBKEY
                IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
              END IF
            END IF
C     
          ELSE IF (SUBKEY.EQ.'GRAD') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            GRADMAX = VALUE(ICOUNT)
            IGRADMAX = 1
          ELSE IF (SUBKEY.EQ.'MINB') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NBGMIN = NINT(VALUE(ICOUNT))
          ELSE IF (SUBKEY.EQ.'PLOT') THEN
            BADPLOT = .TRUE.
          ELSE IF (SUBKEY.EQ.'DUMP') THEN
            BADPLOT2 = .TRUE.
          ELSE
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          END IF
C     
          IF (ICOUNT.LT.NTOK) GOTO 650
        END IF
C     
C---- POSTREF... Do Post refinement to refine missetting angles (and
C     optionally cell parameters) at the end of each pack. This requires
C     the use of profile fitting and adding in partials from the next
C     image, and therefore can only be used with IP data.   
C     
      ELSE IF (KEY.EQ.'POST') THEN
        IF (.NOT.IMGP) THEN
          WRITE(IOUT,6019) KEY
          IF (ONLINE) WRITE(ITOUT,6019) KEY
          BADLINE = .TRUE.
          GOTO 50
        END IF
C     
C---- NUPR_INT .eq. TRUE for integration run, FALSE for refinement only;
C     only used if NEWPREF.
C     
        nupr_int = .true.
        POSTREF = .TRUE.
        SUMPART = .TRUE.
        SPOSTREF = .TRUE.
        SSUMPART = .TRUE.
        INWIDTH = 0
        INADD = 0
chrp 25.05.2006        INSEG = 0
chrp 25.05.2006        INFIX = 0
        IF (NTOK.GT.1) THEN
          ICOUNT = 1
 660      ICOUNT = ICOUNT + 1
          SUBKEY = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     *********** 
          CALL CCPUPC(SUBKEY)
C     ***********       
C     
 661      CONTINUE
C     
C---- OFF to suppress post refinement
C     
          IF (SUBKEY(1:3).EQ.'OFF') THEN
            POSTREF = .FALSE.
            SPOSTREF = .FALSE.
C     
C---- POSThoc...does not assume images are in register
C     
          ELSE IF (SUBKEY.EQ.'POST') THEN
            NEWPREF = .TRUE.

C     
C---- MULTiple...allows use of partials reflections over multiple images
C     default from Version 6.2.0
C     
          ELSE IF (SUBKEY.EQ.'MULT') THEN
            IMULTI = 1
            NEWPREF = .TRUE.
C     
C---- NOMULTiple... normal postrefinement using partials over only two
C     images;
C     
          ELSE IF (SUBKEY.EQ.'NOMU') THEN
            IMULTI = 1
            NEWPREF = .FALSE.
            NUSPOT = .FALSE.
C     
C---- SINGle...only use data from one image at a time in postrefinement.
C     In this mode, cell parameters are not refined unless
C     explicitly UNFIXED
C     
          ELSE IF (SUBKEY.EQ.'SING') THEN
            NADD = 1
            PRMODE = .TRUE.
C     
C---- LATTICE  specifies the lattice number that is being refined
C     in multi-lattice processing. This is only used to specify
C     the lattice in the CELL_REFINE_RESPONSE XML to iMosflm
C     
          ELSE IF (SUBKEY.EQ.'LATT') THEN
            ICOUNT = ICOUNT + 1                                  
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NPRLATT = NINT(VALUE(ICOUNT))
C     
C---- NREF  specifies minimum number of reflections for post-refinement
C     to be carried out
C     
          ELSE IF (SUBKEY.EQ.'NREF') THEN
            ICOUNT = ICOUNT + 1                                  
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NPRMIN = NINT(VALUE(ICOUNT))
C     
C---- NRMOSAIC  specifies minimum number of reflections for post-refinement
C     of mosaicity to be carried out (NREFMOS)
C     
          ELSE IF (SUBKEY.EQ.'NRMO') THEN
            ICOUNT = ICOUNT + 1                                  
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NREFMOS = NINT(VALUE(ICOUNT))
C     
C---- PRINt  additional print output for postrefinement
          ELSE IF (SUBKEY.EQ.'PRIN') THEN
            IPRINTP = 1
C     
C---- SEGMENTS...use NSEG different non-contiguous segments for the 
C     post-refinement
C     
          ELSE IF (SUBKEY(1:3).EQ.'SEG') THEN
            NUPR_INT = .FALSE.
            ICOUNT = ICOUNT + 1                                  
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            IF (MULTISEG) THEN
              NSEGOLD = NSEG
            ELSE
c     
c     hrp 15022002
C---- for jobs from the command line, if we have just integrated an
C     image then we need to make sure FIRSTTIME is true for a multiseg
C     run, and also that NRUN = 0
c     
              FIRSTTIME = .TRUE.
              NRUN = 0
            END IF
            NSEG = NINT(VALUE(ICOUNT))
            IF (MULTISEG.AND.(NSEGOLD.NE.NSEG).AND.(NSEGOLD.NE.0))
     +           THEN
              WRITE(IOUT,FMT=7540) NSEG, NSEGOLD
              IF (ONLINE) WRITE(ITOUT,FMT=7540) NSEG, NSEGOLD
 7540         FORMAT(1X,'**** ERROR ****',/,1X,'The number of s',
     +             'egments has been given as ',I3,' but on a p',
     +             'revious POSTREF keyword',/,1X,'it was given',
     +             ' as',I3,/,1X,'It has been reset to original',
     $             ' value.')
              NSEG = NSEGOLD
            END IF
            MULTISEG = .TRUE.
            DONESEG = .TRUE.
            IF (NSEG.EQ.1) WAITINP = .TRUE.
            INSEG = 1
C     
C---- hrp 19022002
C     SPECIAL CASE; if we have two (or more) POSTREF SEGMENT 1 runs, we
C     need 
C     to make sure GENFILE is open! This addresses the symptom rather than the 
C     cause - GENFILE is closed after P/r following integration, and 
C     reflection file is also o/p at this point - later than one might
C     expect.
C     
            IF(.NOT.GENOPEN)THEN
              MTZOUT = 1
C
C---- if GENFILE hasn't been inititalized, set it to GENFILE
C
              IF(GENFILE(1:8).EQ.'________')GENFILE = 'GENFILE'
              CALLED = 'QOPEN'
              COMMENT = 'To open generate file'
              IF (DEBUG(99)) 
     +           CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
              CALL QOPEN(IUNIT,GENFILE,'UNKNOWN')
              GENOPEN = .TRUE.
            ENDIF

C     
C---- NOSEGMENT  turns off POSTREF SEGMENT
C     
          ELSE IF (SUBKEY.EQ.'NOSE') THEN
C     
C---- if we have done postrefinement in this run, then we need to set
C     NEWGENF = .true. so that STARTMTZ is called at the appropriate
C     point...
C     
            IF(DONESEG)NEWGENF = .TRUE.
            DONESEG = .FALSE.
            MULTISEG = .FALSE.
            RPTFIRST = .FALSE.
            FIRSTTIME = .TRUE.
            NSEG = 0
            NSEGOLD = NSEG
            IF(GENOPEN)THEN
C
C---- Print for program flow
C
              CALLED = 'QCLOSE'
              COMMENT = 'Close generate file when processing NOSEGMENT'
              IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
              CALL QCLOSE(IUNIT)
              GENOPEN = .FALSE.
c
chrp09052002
c
              NEWGENF = .TRUE.
            END IF
C     ******************
C     
C---- if we've just done postrefinement, MTZOPEN will be true
C     
            IF (MTZOPEN) THEN
              MTZPRT = 1
C     *********************
              if (debug(52))then
                write(iout,fmt=7744)
                if (online) write(itout,fmt=7744)
 7744           format(1x,'**closing MTZ file in control 3 ***')
              end if
              CALLED = 'LWCLOS'
              COMMENT = 'To close MTZ file '
              IF (DEBUG(99)) 
     +           CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
              CALL LWCLOS(MTZOUT,MTZPRT)
              MTZOPEN = .FALSE.
              IF(LBEST.AND.BESTOPEN)THEN
                CLOSE(BESTHKL)
                BESTOPEN = .FALSE.
              ENDIF
C     *********************
            END IF
            NRUN = 0
C     
C---- FRMIN... set FRACMIN, minimum allowed fraction for post-refinement
C     
          ELSE IF (SUBKEY.EQ.'FRMI') THEN
            ICOUNT = ICOUNT + 1                                  
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
C     
C---- FRMAX... set FRACMAX, minimum allowed fraction for post-refinement
C     
            FRACMIN = VALUE(ICOUNT)
          ELSE IF (SUBKEY.EQ.'FRMA') THEN
            ICOUNT = ICOUNT + 1                                  
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            FRACMAX = VALUE(ICOUNT)
C     
C---- ADD...use NADD images to do post-refinement. By default the number
C     is chosen to give a wedge of at least 5 degrees of data, but
C     can be set explicitly here
C     
          ELSE IF (SUBKEY.EQ.'ADD') THEN
            ICOUNT = ICOUNT + 1                                  
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NADD = NINT(VALUE(ICOUNT))
            PRMODE = .TRUE.
            INADD = 1
C     
            IF (NADD.GT.NIMAX) THEN
              WRITE(IOUT,FMT=6096) NADD,NIMAX
              IF (ONLINE) WRITE(ITOUT,FMT=6096) NADD,NIMAX
 6096         FORMAT(1X,'**** FATAL ERROR ****',/,1X,'You have ',
     +             'asked for the post-refinement to be done over',
     +             I3,/,1X,'images but this exceeds the maximum',
     $             ' allowed (',I3,/,1X,'Either reduce WIDTH or ', 
     +             'change parameter NIMAX and recompile')
              NSHUTERR = 36
              CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
            END IF
C     
C---- WIDTh... angular width for postrefinement
          ELSE IF (SUBKEY.EQ.'WIDT') THEN
            ICOUNT = ICOUNT + 1                                  
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            ANGWIDTH = VALUE(ICOUNT)
            PRMODE = .TRUE.
            INWIDTH = 1
C     
C---- FIXED... refine using a wedge of data but keep all missetting
C     angles within the wedge constant.
C     
          ELSE IF (SUBKEY.EQ.'FIXE') THEN
            FIXED = .TRUE.
C     
C---- SHIFtfac... ONLY used when doing postrefinement using data from
C     more than a single image. If on the first run of
C     postrefinement (ie on the first NADD images) the
C     shift in cell parameters is greater than SHIFTFAC
C     times the sd of that parameter, go back to the
C     first pack and repeat measurement of the first
C     NADD images. Do this up to a total of NRPT times
C     
          ELSE IF (SUBKEY.EQ.'SHIF') THEN
            ICOUNT = ICOUNT + 1                                  
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            SHIFTFAC = VALUE(ICOUNT)
C     
C---- REPEat NRPT.. maximum number of times to repeat reprocessing of
C     first
C     NADD images (see above)
C     
          ELSE IF (SUBKEY.EQ.'REPE') THEN
            ICOUNT = ICOUNT + 1                                  
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            NRPT = NINT(VALUE(ICOUNT))
C     
C---- SDFAC, only reflections with I > SDFAC*SIG(I) are used
C     
          ELSE IF (SUBKEY.EQ.'SDFA') THEN
            ICOUNT = ICOUNT + 1                                  
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            SDFAC = VALUE(ICOUNT)
C     
C---- RESOlution limits (inner and outer)
          ELSE IF (SUBKEY.EQ.'RESO') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            PRRES1 = VALUE(ICOUNT)
C     
C---- Test for second resolution limit
C     
            IF (ICOUNT.LT.NTOK .AND. ITYP(ICOUNT+1).EQ.2) THEN
              ICOUNT = ICOUNT + 1
C     *******************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
              PRRES2 = VALUE(ICOUNT)
            END IF
C     
C---- MAXResidual... terminate job if postref residual exceeds this
C     factor times (mosaic spread + beam divergence)
C     
          ELSE IF (SUBKEY.EQ.'MAXR') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            RSDMAX = VALUE(ICOUNT)
C     
C---- MAXShift...  Remeasure current image if shift EXCEEDS this limit.
C     In case of several images being used together, the
C     maximum allowed shift is set to the GREATER of
C     SHIFTMAX and SHIFTFAC*(sigma of missets)
C     
          ELSE IF (SUBKEY.EQ.'MAXS') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            SHIFTMAX = VALUE(ICOUNT)
C     
C---- GROUP...group images together for post-refinement,
C     useful if there are not many reflections on each image.
C     
          ELSE IF (SUBKEY.EQ.'GROU') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            NVIRBAT = VALUE(ICOUNT)
            IVIRBAT = 1
C     
C---- BEAM sets refined beam parameters (horizontal, vertical
C     divergences
C     and mosaic spread)
C     
          ELSE IF (SUBKEY.EQ.'BEAM') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            PRNS = NINT(VALUE(ICOUNT))
C     
C---- USEBeam...use the refined beam parameters in reflection list
C     generation
C     
          ELSE IF (SUBKEY.EQ.'USEB') THEN
            USEBEAM = .TRUE.
C     
C---- Check for turning OFF using refined beam parameters
C     
            IF (ICOUNT.LT.NTOK) THEN
              KEY2 = LINE(IBEG(ICOUNT+1):IEND(ICOUNT+1))
              CALL CCPUPC(KEY2)
              IF (KEY2(1:3).EQ.'OFF') THEN
                USEBEAM = .FALSE.
                ICOUNT = ICOUNT + 1
              END IF
            END IF
C     
C---- SMOOTH... smooth the missetting angles (the number of images will be the
C               same as for smoothing the mosaic spread.
C     
          ELSE IF (SUBKEY.EQ.'SMOO') THEN
            SMOOTHPHI = .TRUE.
C     
C---- NOSMOOTH... turn off smoothing the missetting angles
C     
          ELSE IF (SUBKEY.EQ.'NOSM') THEN
            SMOOTHPHI = .FALSE.
C     
C---- MOSSMOOTH... number of images to use in smoothing the refined 
C     mosaic spread
C     
          ELSE IF (SUBKEY.EQ.'MOSS') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            NSMOOTH = NINT(VALUE(ICOUNT))
C     
C---- MOSADD... percentage of mosaic spread to be added to refined value as a
C       safety factor.
C     
          ELSE IF (SUBKEY.EQ.'MOSA') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            ETAFRAC  = VALUE(ICOUNT)



C     
C---- UNFIx cell parameters. If only using a single image for refinement
C     (NADD=1) default is all cell parameters are not refined, but
C     refinement can be turned on using UNFIX.
C     If (NADD>1) default is to refine all variable cell parameters
C     unless they are FIXED.
C     
          ELSE IF (SUBKEY.EQ.'UNFI') THEN
            PRCELL = .TRUE.
C     
C---- Default is ALL cell parameters fixed for single image refinement,
C     ALL parameters refined for multiple image refinement.
C     
 662        ICOUNT = ICOUNT + 1
            IF (ICOUNT.LE.NTOK) THEN
C     
C---- Read next token and test against mosaic,a,b,c,alpha,beta,gamma,ALL
              SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     *********** 
              CALL CCPUPC(SUBKEY)
C     ***********       
              IF (SUBKEY.EQ.'MOSA') THEN                
                USEBEAM = .TRUE.
                GOTO 662
              END IF
C
              IF (SUBKEY.EQ.'ALL') THEN                
                DO 663 I = 1,6
                  UNFIX(I) = .TRUE.
                  FCELL(I) = .FALSE.
 663            CONTINUE
                GOTO 662
              END IF
C     
              DO 664 I = 1,6
                IF (SUBKEY.EQ.SABC(I)(1:4)) THEN
                  UNFIX(I) = .TRUE.
                  FCELL(I) = .FALSE.
                  GOTO 662
                END IF
 664          CONTINUE
C     
C---- SUBKEY not recognised, may be another subkeyword so jump back
C     First check if mosaic or any cell parameters have been unfixed.
C
              IF (USEBEAM) GOTO 661     
              DO 666 I = 1,6
                IF (UNFIX(I)) GOTO 661
 666          CONTINUE
C     
C---- None unfixed... error
              WRITE (IOUT,FMT=6230)
              IF (ONLINE) WRITE (ITOUT,FMT=6230)
 6230         FORMAT(/,1X,'**** ERROR ****, parameters to be',
     +             ' unfixed must be specified (a,b,c,alpha,beta,',
     +             'gamma, mosaic)')
            ELSE
C     
C---- ICOUNT now GT NTOK, check that some cell parameters have been
C     unfixed
C     
              DO 668 I = 1,6
                IF (UNFIX(I).OR.USEBEAM) GOTO 670
 668          CONTINUE
C     
C---- None unfixed, print error message
              WRITE (IOUT,FMT=6230)
              IF (ONLINE) WRITE (ITOUT,FMT=6230)
 670          CONTINUE
            END IF
C     
C---- FIX cell parameters or mosaic spread for multiple image 
C     post-refinement
C     
          ELSE IF (SUBKEY.EQ.'FIX') THEN
            PRCELL = .TRUE.
            INFIX = 1
C     
C     
C---- Default is ALL cell parameters fixed for single image refinement,
C     ALL parameters refined for multiple image refinement.
C     
 672        ICOUNT = ICOUNT + 1
            IF (ICOUNT.LE.NTOK) THEN
C     
C---- Read next token and test against a,b,c,alpha,beta,gamma,mosaic
              SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     *********** 
              CALL CCPUPC(SUBKEY)
C     ***********       
              IF (SUBKEY.EQ.'ALL') THEN
                DO 673 I = 1,6
                  FCELL(I) = .TRUE.
 673            CONTINUE
                GOTO 672
              END IF
C     
              IF (SUBKEY.EQ.'MOSA') THEN
                USEBEAM = .FALSE.
                GOTO 672
              END IF
C
              DO 674 I = 1,6
                IF (SUBKEY.EQ.SABC(I)(1:4)) THEN
                  FCELL(I) = .TRUE.
                  GOTO 672
                END IF
 674          CONTINUE
C     
C---- SUBKEY not recognised, may be another subkeyword so jump back
C     First check if mosaic or any cell parameters have been fixed.
C
              IF (.NOT.USEBEAM) GOTO 661     
              DO 676 I = 1,6
                IF (FCELL(I)) GOTO 661
 676          CONTINUE  
C     
C---- None fixed... error
              WRITE (IOUT,FMT=6232)
              IF (ONLINE) WRITE (ITOUT,FMT=6232)
 6232         FORMAT(/,1X,'**** ERROR ****, parameters to be',
     +             ' fixed must be specified (a,b,c,alpha,beta,',
     +             'gamma, mosaic)')
            ELSE
C     
C---- ICOUNT now GT NTOK, check that mosaic spread or some cell 
C     parameters have been fixed
C     
              IF (.NOT.USEBEAM) GOTO 680
C
              DO 678 I = 1,6
                IF (FCELL(I)) GOTO 680
 678          CONTINUE
C     
C---- None fixed, print error message
              WRITE (IOUT,FMT=6232)
              IF (ONLINE) WRITE (ITOUT,FMT=6232)
 680          CONTINUE
            END IF
C
C---- re-scaling the Normal matrix on or off
C
          ELSE IF (SUBKEY.EQ.'RESC') THEN
            ICOUNT = ICOUNT + 1
            IF( ICOUNT .LE. NTOK) THEN
              SUBKEY =  LINE(IBEG(ICOUNT):IEND(ICOUNT))
              CALL CCPUPC(SUBKEY)
              IF(SUBKEY(1:2).EQ.'ON')THEN
                RESCALE = .TRUE.
              ELSEIF(SUBKEY(1:3).EQ.'OFF')THEN
                RESCALE = .FALSE.
              ELSE
                BADLINE = .TRUE.
                GOTO 50
              ENDIF
            ELSE
              RESCALE = .TRUE.
            ENDIF
          ELSE IF (SUBKEY.EQ.'PART') THEN
            ICOUNT = ICOUNT + 1
            IF(ITYP(ICOUNT).EQ.2)THEN
              IF(VALUE(ICOUNT).GE.0.2.AND.VALUE(ICOUNT).LE.0.5)THEN
                PARTITION  = VALUE(ICOUNT)
                IPARTITION = 1
              ELSE
                WRITE(IOUT,FMT=7543)VALUE(ICOUNT)
                IF(ONLINE)WRITE(ITOUT,FMT=7543)VALUE(ICOUNT)
 7543           FORMAT(1X,'***** ERROR *****',/,1X,
     +               '***** ERROR *****',/,1X,'You supplied a value ',
     $               'type ',F4.2,' for PARTITION, but it can only ',
     +               'be in the range 0.2 - 0.5')
              ENDIF
            ENDIF
              
C
C---- ICE, exclude reflections that may lie in ice rings from post refinement
C
          ELSE IF (SUBKEY(1:3).EQ.'ICE') THEN
            PREFEXICE = .TRUE.
C
C---- Test for a number that gives the half-widths of the exclusion shells
C     as a fraction of the actual resolution.
C
            IF (ICOUNT.LT.NTOK .AND. ITYP(ICOUNT+1).EQ.2) THEN
              ICOUNT = ICOUNT + 1

C     *********************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *********************************************
              ICEWIDTHPREF = VALUE(ICOUNT)
            END IF
C
C----  Don't exclude reflection in resolution shells (unless specified
C     by RESOLUTION EXCLUDE)
C
          ELSE IF (SUBKEY(1:4).EQ.'NOIC') THEN
            PREFEXICE = .FALSE.
C
          ELSE
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          END IF
C     
          IF (ICOUNT.LT.NTOK) GOTO 660
          IF ((MULTISEG).AND.((INADD.GT.0).OR.(INWIDTH.GT.0))) THEN
            WRITE(IOUT,FMT=7542)
            IF (ONLINE) WRITE(ITOUT,FMT=7542)
 7542       FORMAT(1X,'***** ERROR *****',/,1X,
     +           '***** ERROR *****',/,1X,'Do not give ',
     +           'WIDTH or ADD subkeywords with the SEGMENT',
     $           ' subkeyword.')
          END IF
C
C---- If a SEGMENT keyword has been given, which flags cell refinement,
C     then set all cell refinement flags to allow refinement UNLESS a
C     FIX keyword has been given. This is because if a previous integration
C     run has been carried out, all cell parameters will be fixed
C     (by default). This should really be handled by an END TASK.
C
          IF ((INSEG.EQ.1).AND.(INFIX.EQ.0)) THEN
            DO 675 I = 1,6
              FCELL(I) = .FALSE.
 675        CONTINUE
          END IF
C
C---- end of IF NTOK.GT.1
        END IF
C     
C---- BRIEf mini output for online
C     
      ELSE IF (KEY.EQ.'BRIE') THEN
        BRIEF = .TRUE.
        IBRIEF = 20
C     
C---- Open output file
C     
        CALL CCPOPN(-IBRIEF,'MINIOUT',1,1,80,IFAIL)
        WRITE(IBRIEF,FMT=6202)
 6202   FORMAT(1X,'PROGRAM MOSFLM')
C     
C---- IMAGe  Allow examination of an image file. When run interactively
C     the image can be autoindexed and integrated etc via the
C     menu options. In batch, this keyword specifies images to be
C     used for autoindexing.
C     
      ELSE IF ((KEY.EQ.'IMAG').OR.(KEY.EQ.'POWD')) THEN
C
C---- Give a warning that this keyword is no longer recommended if
C     not being used with new gui - NOT here hrp 04.10.2012 - see
C     after (KEY.eq.'RUN'
c        IF (.NOT.SOCKLO) THEN
c          WRITE(IOUT,FMT=7730)
c          IF (ONLINE) WRITE(ITOUT,FMT=7730)
c        END IF
C     
C---- First check that a PROCESS keyword has not been given
C     
        IF (IPROKWD.GT.0) THEN
          WRITE(IOUT,FMT=7300)
          IF (ONLINE) WRITE(ITOUT,FMT=7300)
 7300     FORMAT(/,1X,'***** ERROR *****',/,1X,'IMAGE and PROCESS ',
     +         'keywords must NOT be given in the same "run".',/,1X,
     +         'This keyword will be ignored')
          BADLINE = .TRUE.
          GOTO 50
        END IF
C     
        IMGKWD = 1
        IANGLE = 0
C     
  
        IF(.NOT.HEADONLY)POWDER = .TRUE.
        IF (ONLINE.AND..NOT.SOCKLO) DISPMENU = .TRUE.

        IF (NTOK.LT.2) THEN
          WRITE(ITOUT,FMT=6204)
          IF (BRIEF) WRITE(IBRIEF,FMT=6204)
 6204     FORMAT(1X,'Give filename of  image as second parameter',
     +         /,1X,'eg IMAGE /hx0/andrew/catx1_001.image')
          IF (COMREAD) THEN
            COMREAD = .FALSE.
            ITIN = ITINS
            CLOSE (UNIT=ICOMM)
          END IF
          POWDER = .FALSE.
          DISPMENU = .FALSE.
          BADLINE = .TRUE.
          GOTO 50
        ELSE
C     
C---- If a template has been given, only the image number is given here
C     
          IF (TEMPLATE.OR.HDF5FILE) GOTO 764
C     
          WAXFN = ' '
          WAXFN = LINE(IBEG(2) :IEND(2))
C     
C---- Check that the filename conforms to the standard. It must be of
C     the 
C     form:  ABCDE_###.EXT or ABCDE-###.EXT where # is a number and the
C     initial part of the name "ABCDE" can be up to 40 chars long
C     
C     Search for the separator between text and number...can be "_" or "
C     -"
C     
          NCH = LENSTR(WAXFN)
          DO 758 I = NCH,1,-1
            IF ((WAXFN(I:I).EQ.'_').OR.(WAXFN(I:I).EQ.'-')) THEN
              SEPCHAR = WAXFN(I:I)
              GOTO 759
            END IF
 758      CONTINUE
C     
C---- No separator found. fatal error
C     
          WRITE(IOUT,FMT=6207)                       
          IF (ONLINE) WRITE(ITOUT,FMT=6207)                       
 6207     FORMAT(/,1X,'***** FATAL ERROR *****',/,1X,
     +         'Image filenames must be of the form ABCDE_###.ext',
     +         ' or ABCDE-###.ext where the initial string',/,1X,
     +         'can be up to 40 characters long and must be ',
     +         'separated from a 3 digit number by',/,1X,'a ',
     +         '"_" or "-", and the extension (ext) can be up',
     +         ' to 8 characters long.',/,1X,'Use the TEMPLATE',
     +         ' keyword to read other types of filename.')
          NSHUTERR = 37
          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
C     
C---- Extract directory (if any) and reset WAXFN to string stripped of
C     the directory
C     
 759      NCH2 = NCH
          DO 760 I = NCH,1,-1
            IF ((WAXFN(I:I).EQ.'/').OR.(WAXFN(I:I).EQ.']')) THEN
              NDIR = 1
              FDISK(1) = WAXFN(1:I)
              TEMPCH = WAXFN
              WAXFN = ' '
              WAXFN = TEMPCH(I+1:NCH)
              NCH2 = LENSTR(WAXFN)
              GOTO 761
            END IF
 760      CONTINUE
C     
 761      NCH = NCH2
C     
C---- Check if extension has been given, if so transfer it to ODEXT
C     and reset WAXFN to filename excluding extension
C     
          DO 684 I = NCH,1,-1
            IF (WAXFN(I:I).EQ.'.') THEN
              IF ((NCH-I).GT.8) THEN
                WRITE(IOUT,FMT=6207)
                IF (ONLINE) WRITE(ITOUT,FMT=6207)
                NSHUTERR = 38
                CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
              END IF
              ODEXT = ' '
              ODEXT = WAXFN(I+1:NCH)
              TEMPCH = WAXFN
              WAXFN = ' '
              WAXFN = TEMPCH(1:I-1)
              IMGFN = WAXFN
              J = I
              K = J-1
              GOTO 762
            END IF
 684      CONTINUE
          J = NCH
          K = J
C     
C---- Extract image template, assuming filename of form ABCD_00N.ext and
C     searching for _, then set WAXFN to filename stripped of number 00N
C     
 762      DO 685 I = J,1,-1
            IF (WAXFN(I:I).EQ.SEPCHAR) THEN
C     
C---- Get the image number as a string
C     
              IMGNUM = WAXFN(I+1:K)
              TEMPCH = WAXFN
              WAXFN = TEMPCH(1:I-1)
              IMGTEMPL = WAXFN(1:I-1)
              IDENT = WAXFN
              RESTIDENT = IDENT
              IIDENT = 1
              GOTO 692
            END IF
 685      CONTINUE
C     
C     
C---- Extract image number as a value from the string
C     
C     ******************************************
 692      NCH2 = LENSTR(IMGNUM)
          IF (NCH2.NE.3) THEN
            WRITE(IOUT,FMT=6207)                       
            IF (ONLINE) WRITE(ITOUT,FMT=6207)  
            NSHUTERR = 39                     
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF
          CALL MPARSE(IMGNUM,IBEG2,IEND2,ITYP2,VALUE2,IDEC2,NTOK2)
          CALL MKEYNM(1,1,IMGNUM,IBEG2,IEND2,ITYP2,NTOK2)
C     *******************************************
C---- Trap error in number
C     
          IF (IOERR) THEN
            BADLINE = .TRUE.
            GOTO 50
          END IF
          ID = NINT(VALUE2(1))
C---- Should never want to input more than a single IMAGE keyword, as
C     only one can be examined !
C     
C     
          RESTID = ID     
          IDIMG(1) = ID
          NCHAR = LENSTR(ODEXT)
          IF(MOSEST.AND..NOT.AUTOINDX)THEN
            MOSIMAG = ID
            NAUTO = 1
CHRP 15.05.2013            IDAUTO(1) = MOSIMAG
            NOIMG(1) = MOSIMAG
            IDPACK(1) = MOSIMAG
            NRUN = 1
            POWDER = .TRUE.
            NODISPLAY = .TRUE.
          ENDIF

C     
C---- only set PACK true if extension is "pck" and not "pck1200" etc
C     
          PACK = (((ODEXT(1:3).EQ.'PCK').OR.(ODEXT(1:3).EQ.'pck'))
     +         .AND.(NCHAR.EQ.3))
          IF (PACK) USEHDR = .FALSE.
          GOTO 767
C     
C---- Deal with TEMPLATE case, get image number
C     
C     ************************************
 764      CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          ID = NINT(VALUE(2))
C     
C---- Set up WAXFN and image number as a string as these are used to
C     build
C     filenames for spots and matrix.
C     
          IF (LENSTR(TEMPLSTART).GT.0) THEN
            WAXFN = TEMPLSTART(1:LENSTR(TEMPLSTART))
          ELSE IF (LENSTR(TEMPLEND).GT.0) THEN
            WAXFN = TEMPLEND(1:LENSTR(TEMPLEND))
          ELSE
            WAXFN = 'X'
          END IF
C     
C---- If image template ends in a ".", eg from a filename like lysox1
C     .001
C     then remove the "." from WAXFN
C     
          I = LENSTR(WAXFN)
          IF (I.GT.1) THEN
            IF (WAXFN(I:I).EQ.'.') WAXFN(I:I) = ' '
          END IF
C     
          CALL IMGMAKE(NTDIG,ID,IMGNUM)
C     
C---- Should never want to input more than a single IMAGE keyword, as
C     only one can be examined !
C     
C     
          IDIMG(1) = ID
          RESTID = ID
C     ..
C     .. wanted for HDF5 files
          THISIMG = ID
C     
C---- Get phi value (if given)
C     
 767      IF (NTOK.GT.2) THEN
            ICOUNT = 2
 687        ICOUNT = ICOUNT + 1
            SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     **************
            CALL CCPUPC(SUBKEY)
C     **************
            IF (SUBKEY.EQ.'PHI') THEN
              ICOUNT = ICOUNT + 1
C     *******************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
              IF (IOERR) THEN
                BADLINE = .TRUE.
                GOTO 50
              END IF
C     
              PHIBEG = VALUE(ICOUNT)
c     hrp08112001                  DEFPHI = .TRUE.
              ICOUNT = ICOUNT + 1
C     
C---- Check for subkeyword "TO"
C     
              IF (ITYP(ICOUNT).EQ.1) THEN
                KEY2 = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     **************
                CALL CCPUPC(KEY2)
C     **************
                IF (KEY2.NE.'TO') THEN
                  WRITE(IOUT,FMT=6592)
                  IF (ONLINE) THEN
                    WRITE(ITOUT,FMT=6592)
                    IF (COMREAD) THEN
                      COMREAD = .FALSE.
                      ITIN = ITINS
                      CLOSE (UNIT=ICOMM)
                    END IF
                    BADLINE = .TRUE.
                    GOTO 50
                  END IF
                  STOP
                END IF
 6592           FORMAT(1X,'*** ERROR ***',/,1X,'Give PHI value ',
     +               ' in the form: PHI 0.0 TO 1.0')
                ICOUNT = ICOUNT + 1
              END IF
C     
C     *******************************************
              CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
              IF (IOERR) THEN
                BADLINE = .TRUE.
                GOTO 50
              END IF
              PHIEND = VALUE(ICOUNT)
              IF (PHIEND.LT.PHIBEG) THEN
                TEMPX = PHIEND
                PHIEND = PHIBEG
                PHIBEG = TEMPX
              END IF
              PHI(1) = 0.5*(PHIBEG+PHIEND)
              PHISTIM(1) = PHIBEG
              RESTPHIB = PHIBEG
              RESTPHIE = PHIEND
              IANGLE = 1
              ISTRT = 1
C     
C---- NODISPLAY to simulate batch job when running online
C     
            ELSE IF (SUBKEY.EQ.'NODI') THEN
              NODISPLAY = .TRUE.
            ELSE
C     
C     Not recognised
C     
              IF (TRAPERR) INPERR = .TRUE.
              WRITE (IOUT,FMT=6130) SUBKEY
              IF (ONLINE) THEN
                WRITE(ITOUT,FMT=6130) SUBKEY
                IF (COMREAD) THEN
                  COMREAD = .FALSE.
                  ITIN = ITINS
                  CLOSE (UNIT=ICOMM)
                END IF
                BADLINE = .TRUE.
                GOTO 50
              END IF
              STOP
            END IF
C     
            IF (ICOUNT.LT.NTOK) GOTO 687
C     
          END IF
          IF (TEMPLATE.or.hdf5file) THEN
            I = ID
          ELSE
C     
C---- Extract image number as a value from the string
C     
C     ******************************************
            CALL MPARSE(IMGNUM,IBEG,IEND,ITYP,VALUE,IDEC,NTOK)
            CALL MKEYNM(1,1,IMGNUM,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
C---- Trap for any error in reading numbers with parser
C     
            IF (IOERR) THEN
              IF (ONLINE) THEN
                IF (COMREAD) THEN
                  COMREAD = .FALSE.
                  ITIN = ITINS
                  CLOSE (UNIT=ICOMM)
                END IF
                BADLINE = .TRUE.
                GOTO 50
              END IF
              STOP
            END IF
C     
            I = NINT(VALUE(1))
          END IF
C     
          NOIMG(1) = I
          IDPACK(1) = I
          NIMAGES = 1
          IRSTRT = 1
C     
C---- If processing interactively want to increment NIMAG here also
C     because it is use in MXDSPL
C     
          IF (ONLINE.AND.(.NOT.NODISPLAY)) NIMAG = NIMAG + 1
        END IF

c     graeme's program mods for the xkey words... DEPRECATED!
CHRP 20.02.2006
C I think this is the only one we still use...
      else if(key .eq. 'XGUI') then
        ONLINE = .TRUE.
        IOUT = 8
        ITOUT = 6
c     let's do this properly, and actually inspect the second word!
         if(ntok .EQ. 2) then
            SUBKEY = LINE(IBEG(2):IEND(2))
            CALL CCPUPC(SUBKEY)
            IF(SUBKEY .EQ. 'ON') THEN
               GUI_SWITCH = .TRUE.
               NEWGUI = .TRUE.
               NODISPLAY = .FALSE.
            ELSE
               NEWGUI = .FALSE.
               GUI_SWITCH = .FALSE.
            END IF
          ELSE
            WRITE(IOUT,FMT=9200)
            IF(ONLINE)WRITE(IOUT,FMT=9200)
 9200       FORMAT(/,' Usage: XGUI <ON/OFF>')
            BADLINE = .TRUE.
            GOTO 50
          END IF
C     
C---- FUDGE FOR ERRORS IN PHI
C     
      ELSE IF (KEY.EQ.'ERRO') THEN
        DO 686 ICOUNT = 2,4
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          ERROR(ICOUNT-1) = VALUE(ICOUNT)
 686    CONTINUE
CAL        WRITE(6,*),'ERROR',ERROR
      ELSE IF (KEY.EQ.'NORE') THEN
        NOREF = .TRUE.
        IF(NTOK.GT.1)THEN
          SUBKEY = LINE(IBEG(2) :IEND(2))
          CALL CCPUPC(SUBKEY)
          IF (SUBKEY(1:3).EQ.'OFF')THEN
            NOREF = .FALSE.
            WRITE(IOUT,FMT=7123)
            IF (ONLINE) WRITE(ITOUT,FMT=7123)
          ENDIF
        ENDIF
        IF(NOREF)THEN
          MINREF = 0
          NCYC = 1
          WRITE(IOUT,FMT=7122)
          IF (ONLINE) WRITE(ITOUT,FMT=7122)
        ENDIF
 7122   FORMAT(1X,'No positional refinement will be done')
 7123   FORMAT(1X,'Positional refinement will be performed')
C     
C---- WAIT for x minutes, y seconds for image to exist; or cycle every 5
C     
C     seconds for X minutes if X <= 0; a third parameter is added to
C     allow
C     an extra wait if there's a slow network or reading off tape.
C     
      ELSE IF (KEY.EQ.'WAIT') THEN
C     *******************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
        WAIT = VALUE(2)*60
        IF(NTOK.GT.2)WAIT=WAIT+(VALUE(3))
        IF(NTOK.EQ.4)THEN
          DELAY = MAX(INT(VALUE(4)),MIN(10,INT(WAIT/10)))
        ELSE
          DELAY = 0
        ENDIF
        WRITE(IOUT,FMT=6440) INT(WAIT/60),
     $       INT(WAIT-(60*INT(WAIT/60))),DELAY
        IF (ONLINE) WRITE(ITOUT,FMT=6440) 
     $       INT(WAIT/60),INT(WAIT-(60*INT(WAIT/60))),DELAY
 6440   FORMAT(1X,'Image file wait time set to ',I3,' minutes ',
     $       I2,' seconds with a delay of ',I2,' seconds to ',
     $       'allow for slow writing of images')
C     
C---- BIAS...add a constant value to all pixel values
C     
      ELSE IF (KEY.EQ.'BIAS') THEN
C     *******************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
        ICONST = NINT(VALUE(2))
        IF (MACHINE.EQ.'RAXI') THEN
          WRITE(IOUT,FMT=6449)
          IF (ONLINE) WRITE(ITOUT,FMT=6449)
 6449     FORMAT(1X,'***** BIAS keyword not implemented for RAXIS ',
     +         'images')
          BADLINE = .TRUE.
          GOTO 50
        END IF
        IF (ICONST.LT.0) THEN
          WRITE(IOUT,FMT=6446)
          IF (ONLINE) WRITE(ITOUT,FMT=6446)
          ICONST = 0
 6446     FORMAT(//,1X,'Negative BIAS values are NOT permitted,',
     +         ' reset to zero')
        ELSE
          WRITE(IOUT,FMT=6448) ICONST,ICONST
          IF (ONLINE) WRITE(ITOUT,FMT=6448) ICONST,ICONST
 6448     FORMAT(1X,'A constant value of',I3,' will be added to all',
     +         ' pixel values in the image',/,1X,'The ADC offset ',
     +         'will also be reset to',I3,' so that sigma estimates',
     +         ' are still correct')
          IDIVIDE = ICONST
        END IF
C     
C---- Polarisation, can be PINHOLE, MONOCHROMATOR, MIRRORS, SYNCHROTRON
C     
      ELSE IF (KEY.EQ.'POLA') THEN
        IPOLAR = 1
        IF (NTOK.LT.2) THEN
          WRITE(IOUT,FMT=6452)
          IF (ONLINE) THEN
            WRITE(ITOUT,FMT=6452)
            IF (COMREAD) THEN
              COMREAD = .FALSE.
              ITIN = ITINS
              CLOSE (UNIT=ICOMM)
            END IF
            BADLINE = .TRUE.
            GOTO 50
          END IF
          NSHUTERR = 40
          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
 6452     FORMAT(1X,'*** ERROR ***',/,1X,'The type of polarisation',
     +         ' correction must',
     +         ' be given',/,1X,'Possibilities are PINHOLE,'
     +         ,' MONOCHROMATOR MIRRORS (treated as pinhole)',
     $         ' or SYNCHROTRON',/,1X,'followed by the degree',
     $         ' of polarisation (default 0.95)')
        ELSE
C
C---- Allow a number to be specified in place of keywords
C
          IF (ITYP(2).EQ.2) THEN
C     ************************************
            CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
            TOR = VALUE(2)
            IMONO = 3
          END IF
          INMONO = 1
          SUBKEY = LINE(IBEG(2):IEND(2))
C     *********** 
          CALL CCPUPC(SUBKEY)
C     ***********
          IF ((SUBKEY.EQ.'PINH').OR.(SUBKEY.EQ.'MIRR')) THEN
            IMONO = 0
          ELSE IF (SUBKEY.EQ.'MONO') THEN
            IMONO = 1
          ELSE IF (SUBKEY.EQ.'SYNC') THEN
            IMONO = 2
C     
C---- Test for specified degree of polarisation
C     
            IF (NTOK.GT.2) THEN
C     ************************************
              CALL MKEYNM(1,3,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
              TOR = VALUE(3)
              ITOR = 1
            ELSE
              TOR = TORSRS
            END IF
          ELSE
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          END IF
        END IF
C     
C---- MAXWIDTH...Maximum reflection width (in degrees)
C     
C     
      ELSE IF (KEY.EQ.'MAXW') THEN
C     ************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
        WMAX = VALUE(2)
        IF (NTOK.EQ.4) THEN
          SUBKEY = LINE(IBEG(3):IEND(3))
          CALL CCPUPC(SUBKEY)
          IF (SUBKEY.EQ.'PAD') THEN
C     ************************************
            CALL MKEYNM(1,4,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
            IPAD = VALUE(4)
          ELSE
C     
C---- Not recognised
C     
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          END IF
        END IF
C     
C---- SDMON.... Monitor reflections if any observation differs from
C     the weighted mean by more than SDMON sd's when there are
C     multiple fully recorded observations on a single image
C     (or added partials)
C     
      ELSE IF (KEY.EQ.'SDMO') THEN
C     ************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
        SDMON = VALUE(2)
C     
C---- BSWAP...will force opposite of what the program would
C     otherwise do wrt byte-swapping
C     
      ELSE IF (KEY.EQ.'BSWA') THEN
        FIXSWAP = .TRUE.
C     
C     
C---- TWOTHETA swing angle for detector.
C     
      ELSE IF (KEY.EQ.'TWOT') THEN
C     ************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
        TWOTHETA = VALUE(2)
        CAMCON(5) = -TWOTHETA
        ITWOTHETA = 1
C     
C---- *** Keywords associated with finding spots for autoindexing follow
C     **
C     
C-----FINDSPOT
C     
      ELSE IF (KEY.EQ.'FIND') THEN
C        DOFIND = .FALSE.
        FINDSPOT = .TRUE.
        ICOUNT = 1
 800    ICOUNT = ICOUNT + 1
        IF (ICOUNT.GT.NTOK) GOTO 820
        KEY8 = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     *************
        CALL CCPUPC(KEY8)
C     *************
        SUBKEY = KEY8
C     
C---- THRESHOLD... sets the threshold for finding spots
C     
        IF (SUBKEY.EQ.'THRE') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          IF(VALUE(ICOUNT).EQ.0)THEN
             ITHSET = 0
          ELSE
             ITHSET = 1
          ENDIF
          THRESH = VALUE(ICOUNT)
C     
C---- Do not allow threshold lt 1.0
C     
          THRESH = MAX(THRESH,1.0)
C     
C---- RMIN... minimum radius for spot search
C     
        ELSE IF (SUBKEY.EQ.'RMIN') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          RMINSP = VALUE(ICOUNT)
          rmaxsp = sign(rmaxsp,rminsp)
C     
C---- RMAX... maximum radius for spot search
C     
        ELSE IF (SUBKEY.EQ.'RMAX') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          RMAXSP = VALUE(ICOUNT)
          rminsp = sign(rminsp,rmaxsp)
C     
C---- SPLIT... sets spot splitting parameters in X and Y
C     
        ELSE IF (SUBKEY.EQ.'SPLI') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(2,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     

          XSPLIT = VALUE(ICOUNT)
          ICOUNT = ICOUNT + 1
          YSPLIT = VALUE(ICOUNT)
C     
C---- NEW... write new format spots file
C     
        ELSE IF (SUBKEY.EQ.'NEW') THEN
          NEWSPT = .TRUE.
C     
C---- MINX...minimum spot size in X as a function of determined size
C     
        ELSE IF (SUBKEY.EQ.'MINX') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          CUTWXMIN = VALUE(ICOUNT)
C     
C---- MAXX...maximum spot size in X as a function of determined size
C     
        ELSE IF (SUBKEY.EQ.'MAXX') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          CUTWXMAX = VALUE(ICOUNT)
C     
C---- MINY...minimum spot size in Y as a function of determined size
C     
        ELSE IF (SUBKEY.EQ.'MINY') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          CUTWYMIN = VALUE(ICOUNT)
C     
C---- MAXY...maximum spot size in Y as a function of determined size
C     
        ELSE IF (SUBKEY.EQ.'MAXY') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          CUTWYMAX = VALUE(ICOUNT)
C     
C     
C---- MINPIX...minimum number of pixels in a spot
C     
        ELSE IF (SUBKEY.EQ.'MINP') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          NPIXMIN = NINT(VALUE(ICOUNT))
          IPIXMIN = 1
C     
C     
C---- XOFFSET... Displacement of background strip in X. If given,
C     forces background strip to be parallel to Y
C     
        ELSE IF (SUBKEY.EQ.'XOFF') THEN
          IXOFFSET = 1
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          XOFFSET = VALUE(ICOUNT)
          RADY = .TRUE.
          RADX = .FALSE.
C     
C---- YOFFSET... Displacement of background strip in Y. If given,
C     forces background strip to be parallel to X
C     
        ELSE IF (SUBKEY.EQ.'YOFF') THEN
          IYOFFSET = 1
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          YOFFSET = VALUE(ICOUNT)
          RADX = .TRUE.
          RADY = .FALSE.



C     
C---- XMIN...minimum spot X  coordinate (relative to direct beam) in mm
C     
        ELSE IF (SUBKEY.EQ.'XMIN') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          SPXMIN = VALUE(ICOUNT)
C     
C---- YMIN...minimum spot Y  coordinate (relative to direct beam) in mm
C     
        ELSE IF (SUBKEY.EQ.'YMIN') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          SPYMIN = VALUE(ICOUNT)
C     
C---- FIND - actually find spots and write the list to a file
C     
        ELSE IF (SUBKEY.EQ.'FIND') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C
C---- reset some numbers if we are doing a command line or new GUI findspots
C
          IF(SNAUTO.NE.0)THEN
            NAUTO = 0
            NIMAG = 0
            SNAUTO = 0
            DO 805 I = 1,20
              FOUND(I) = .FALSE.
 805        ENDDO
          ENDIF
          IF(NAUTO.EQ.0.OR.DPSDONE)THEN
             NAUTO = 1
          ELSE
             NAUTO = NAUTO + 1
          ENDIF
          ID = NINT(VALUE(ICOUNT))
CHRP 15.05.2013          IDAUTO(NAUTO) = ID
          IDPACK(NAUTO) = ID
CHRP 15.05.2013          NOIMG(NAUTO) = IDAUTO(NAUTO)
          NOIMG(NAUTO) = ID
          DOFIND = .TRUE.
          NODISPLAY = .TRUE.
          POWDER = .TRUE.
c          dpsindex = .false.
c          sdpsindex = .false.
C     
C---- FILE - allow the user to specify the spot filename
C     
        ELSE IF (SUBKEY.EQ.'FILE') THEN
          ICOUNT = ICOUNT + 1
C     .. allow numeric names...
c          IF(ITYP(ICOUNT).EQ.1)THEN
            SPTNAM = LINE(IBEG(ICOUNT):IEND(ICOUNT))
            USERSPOT = .TRUE.
c          ENDIF
C     
C---- PHI...gives phi range for this image
C     
        ELSE IF (SUBKEY.EQ.'PHI') THEN
          IF(DOFIND)THEN
            IF (NMULTI.GT.0) THEN
              WRITE(IOUT,FMT=7273)
              IF (ONLINE) WRITE(ITOUT,FMT=7273)
 7273         FORMAT(1X,'***** ERROR *****',/,1X,'You must only give',
     $             ' one image number on the FIND subkeyword if a ',/,1X
     $             ,'PHI subkeyword is given...see Help library.')
              NSHUTERR = 41
              CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
            END IF

            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(2,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            PHI1 = VALUE(ICOUNT)
            ICOUNT = ICOUNT + 1
            PHI2 = VALUE(ICOUNT)
            PHIBEG = PHI1
            PHIEND = PHI2
            PHISET(NAUTO) = .TRUE.
            PHI(NAUTO) = 0.5*(PHI1+PHI2)
            PHISTIM(NAUTO) = PHIBEG
            PHIRNGA(NAUTO) = ABS(PHI2-PHI1)
            PHIRNG = PHIRNGA(NAUTO)
            IANGLE = 1
          ELSE
            WRITE(IOUT,7274)
            IF (ONLINE) WRITE(ITOUT,FMT=7274)
 7274       FORMAT(/,'***** ERROR *****',/,' The PHI subkeyword is',
     $           ' only valid with the FINDSPOTS keyword if a FIND',
     $           /,' subkeyword has already been given on the SAME ',
     $           'command line',/)
            ICOUNT = ICOUNT + 2
          ENDIF
C
C---- toggle between local background method
C
        ELSE IF (SUBKEY.EQ.'RADI') THEN
          LOCBG = .FALSE.
        ELSE IF (SUBKEY.EQ.'LOCA') THEN
          LOCBG = .TRUE.
C     
C     
C---- MAXPEAKS...maximum number of local maxima within a spot
C     
        ELSE IF (SUBKEY.EQ.'MAXP') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          MAXPKS = MIN(NINT(VALUE(ICOUNT)),MAXNPEAKS)
C
C---- Specify box size for local background
C
        ELSE IF (SUBKEY.EQ.'BOX') THEN
          ICOUNT = ICOUNT + 1
C              ************************************
          CALL MKEYNM(2,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C              ************************************
          NXBOX = VALUE(ICOUNT)
          ICOUNT = ICOUNT + 1
          NYBOX = VALUE(ICOUNT)
C
C---- MAXYEXT Maximum extent of a spot in Y direction
C
        ELSE IF (KEY8.EQ.'MAXYEXT') THEN
          ICOUNT = ICOUNT + 1
C              ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C              ************************************
          MAXYEXT = VALUE(ICOUNT)
C
C---- RMSFAC Multiple peak spots must have rms variation over spot
C     greater than RMSFAC times expected rms from counting stats on a
C     uniform background.
C
        ELSE IF (SUBKEY.EQ.'RMSF') THEN
          ICOUNT = ICOUNT + 1
C              ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C              ************************************
          RMSFAC = VALUE(ICOUNT)
          IRMSFAC = 1
C
C---- ANISOMAX MaximuM ratio of spot sizes in X and Y (either way)
C
        ELSE IF (SUBKEY.EQ.'ANIS') THEN
          ICOUNT = ICOUNT + 1
C              ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C              ************************************
          ANISOMAX = VALUE(ICOUNT)
C     
C---- NOICE ... do not exclude spots in shells at resolution of ice rings
C     
        ELSE IF (SUBKEY.EQ.'NOIC') THEN
          EXCLDICE  = .FALSE.
C     
C---- ICE ...  exclude spots in shells at resolution of ice rings
C     
        ELSE IF (SUBKEY.EQ.'ICE') THEN
          EXCLDICE  = .TRUE.
C
C---- AUTORESL Automatically cut back resolution to 4A if mean intensity
C     between 6-4A is a factor "X" less than mean intensity in 6A-infinity bin
C     
        ELSE IF (KEY8.EQ.'AUTORESL') THEN
          AUTORESLN  = .TRUE.
C
C---- AUTORES1 reset the lower resolution bin boundary (default 6A) for
C     automatically cutting back resolution
C     
        ELSE IF (KEY8.EQ.'AUTORES1') THEN
          ICOUNT = ICOUNT + 1
C              ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C              ************************************
          AUTORES1 = VALUE(ICOUNT)
C
C---- AUTORES2 reset the higher resolution bin boundary (deafult 4A) for
C     automatically cutting back resolution
C     
        ELSE IF (KEY8.EQ.'AUTORES2') THEN
          ICOUNT = ICOUNT + 1
C              ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C              ************************************
          AUTORES2 = VALUE(ICOUNT)
C
C
C---- AUTOFAC1 reset the intensity factor (default 0.2) for
C     automatically cutting back resolution
C     
        ELSE IF (KEY8.EQ.'AUTOFAC1') THEN
          ICOUNT = ICOUNT + 1
C              ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C              ************************************
          AUTOFAC1 = VALUE(ICOUNT)
C
C
C---- AUTOFAC2 reset the reflection number factor (default 0.25) for
C     automatically cutting back resolution
C     
        ELSE IF (KEY8.EQ.'AUTOFAC2') THEN
          ICOUNT = ICOUNT + 1
C              ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C              ************************************
          AUTOFAC2 = VALUE(ICOUNT)
C
C---- NOAUTORESLN turn off AUTORESLN
C
        ELSE IF (KEY8.EQ.'NOAUTORE') THEN
          AUTORESLN  = .FALSE.
C
C---- AUTORING. Automatic exclusion of spots lying on powder rings
C     After diving image into 100 radial bins, if the number of spots
C     in any bin is more than "X" times the average number that bin
C     is excluded. "X" defaults to 4 but can be reset with subkeyword
C     RINGFACT
C
        ELSE IF (KEY8.EQ.'AUTORING') THEN
          AUTORING  = .TRUE.
C
C---- RINGFACT reset the factor (default 4.0) for rejecting "overpopulated"
C     resolution rings (which could be due to ice or other poser diffraction).
C     Rings with more than RINGFACT times the average number reflection are
C     rejected.
C     
        ELSE IF (KEY8.EQ.'RINGFACT') THEN
          ICOUNT = ICOUNT + 1
C              ************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C              ************************************
          RINGFACT = VALUE(ICOUNT)
C
C---- NOAUTORING
C
C
        ELSE IF (KEY8.EQ.'NOAUTORI') THEN
          AUTORING  = .FALSE.
C
C---- AUTOTHRESH. Automatically set a I/sig(I) threshold for spots
C     to be used in autoindexing based on histograms of I/sig(I).
C
        ELSE IF (KEY8.EQ.'AUTOTHRE') THEN
          AUTOTHRESH  = .TRUE.
C
C---- NOAUTOTHRESH
C
C
        ELSE IF (KEY8.EQ.'NOAUTOTH') THEN
          AUTOTHRESH  = .FALSE.
C
C
C---- SUMIMAGE Specify number of images to sum before finding spots
C     for indexing
C
        ELSE IF (KEY8.EQ.'SUMIMAGE') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          NIMGSUMSPOT = VALUE(ICOUNT)
          NIMGSUMSPOT = MAX(NIMGSUMSPOT,1)
          NIMGSUMSPOT = MIN(NIMGSUMSPOT,50)
        ELSE
C     
C     Not recognised
C     
          IF (TRAPERR) INPERR = .TRUE.
          WRITE (IOUT,FMT=6130) SUBKEY
          IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
        END IF
        IF (ICOUNT.LT.NTOK) GOTO 800
 820    CONTINUE
C     
C---- NULLPIX... The value assigned to pixels outside the active area
c                If KEEP is given after the value, this value will not
C                be changed following the gap pixel test for Pilatus
C                detectors.
C     
      ELSE IF (KEY.EQ.'NULL') THEN
C     
C     ************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
        IINULL = .TRUE.
        NULLPIX = NINT(VALUE(2))
        IF (NTOK.GT.2) THEN
          SUBKEY = LINE(IBEG(3):IEND(3))
C         **************
          CALL CCPUPC(SUBKEY)
C         **************
          KEEPNULLPIX = (SUBKEY.EQ.'KEEP')
        ELSE
          KEEPNULLPIX = .FALSE.
        END IF
C     
C---- AUTOINDEX... Autoindexing     HRP --- replaced 06.02.2004
C     
      ELSE IF (KEY6.EQ.'AUTOIN') THEN
CHRP 15.05.2013        CALL KWAUTOIN(IERR,MODE,MAXCELL,MINCELL,IDAUTO,NAUTO,
        CALL KWAUTOIN(IERR,MODE,MAXCELL,MINCELL,NAUTO,
     $       NMULTI,NIMAGES,
     $       NSOL,AUTOINDX,
     $       SDPSINDEX,DPSDONE,SAVIND,TRAPERR,INPERR,PHISET,LINE,
     $       IBEG,IEND,ITYP,VALUE,IDEC,NTOK)
        IF(IERR.NE.0)THEN
          BADLINE = .TRUE.
          IERR = 0
          GOTO 50
        ENDIF
C     
C---- SPOTSEARCH... initial search to determine best threshold and
C     optionally
C     separation and raster parameters
C     
      ELSE IF (KEY.EQ.'SPOT') THEN
        ICOUNT = 1
 860    ICOUNT = ICOUNT + 1
        IF (ICOUNT.GT.NTOK) GOTO 870
        SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     **************
        CALL CCPUPC(SUBKEY)
C     **************
        IF (SUBKEY.EQ.'RMIN') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          RMINSRCH = VALUE(ICOUNT)
          rmaxsrch = sign(rmaxsrch,rminsrch)
        ELSE IF (SUBKEY.EQ.'RMAX') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          RMAXSRCH = VALUE(ICOUNT)
          rminsrch = sign(rminsrch,rmaxsrch)
        ELSE IF (SUBKEY.EQ.'NSEA') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          NSEARCH = NINT(VALUE(ICOUNT))
        ELSE IF (SUBKEY.EQ.'SCAL') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          SCALSRCH = VALUE(ICOUNT)
        ELSE IF (SUBKEY.EQ.'SAFE') THEN
          ICOUNT = ICOUNT + 1
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          ISAFE = VALUE(ICOUNT)
        ELSE
C     
C     Not recognised
C     
          IF (TRAPERR) INPERR = .TRUE.
          WRITE (IOUT,FMT=6130) SUBKEY
          IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
        END IF
        IF (ICOUNT.LT.NTOK) GOTO 860
 870    CONTINUE
C     
C---- NONLINEARITY...Defines CURV, I = I(1.0+CURV*I)
C     
      ELSE IF (KEY.EQ.'NONL') THEN
C     ************************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************************
C     
        CURV = VALUE(2)
C     
C---- BACKSTOP (CENTRE RADIUS) Defines backstop shadow region.
C     Equivalent to LIMITS RCENTRE RMIN
C     
      ELSE IF (KEY6.EQ.'BACKST') THEN
        IIBACK = .TRUE.
        IBACKS = 1
        ICOUNT = 1
 790    ICOUNT = ICOUNT + 1
        IF (ICOUNT.LE.NTOK) THEN
C     
C---- Skip if no more tokens on line
C     
          KEY2 = LINE(IBEG(ICOUNT) :IEND(ICOUNT))
C     
C     **********
          CALL CCPUPC(KEY2)
C     **********
C     
C---- BACKSTOP RADIUS radius of backstop shadow mm
C     
C     
          IF (KEY2.EQ.'RADI') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
            RMIN = VALUE(ICOUNT)*100.0
C     
C---- BACKSTOP CENTRE... centre of backstop shadow circle (mm)
C     
C     
          ELSE IF (KEY2.EQ.'CENT') THEN
            ICOUNT = ICOUNT + 1
C     
C     ************************************
            CALL MKEYNM(2,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
            RMINXINP = VALUE(ICOUNT)*100.0
C     
C---- Check if being input from X-window interface
C     
            RMINX = RMINXINP
            IF ((MODE.EQ.3).AND.(INVERTX))
     +           RMINX =  100.0*NREC*RAST - RMINXINP
            ICOUNT = ICOUNT + 1
            RMINY = VALUE(ICOUNT)*100.0
            IRMINXY = 1
C     
C---- code for automatically finding a circular backstop - placeholder 
C     at the moment
C
          ELSE IF (KEY2.EQ.'AUTO')THEN
            ICOUNT = ICOUNT + 1
            IF(ICOUNT.LE.NTOK.and.ITYP(ICOUNT).EQ.2)THEN
              IERR = VALUE(ICOUNT)
              ICOUNT = ICOUNT + 1
            ELSE
              IERR = 1
            ENDIF
            CALL BACKSTOP(PACK,ID,IERR)
          ELSE
C     
C---- Not recognised
C     
            IF (TRAPERR) INPERR = .TRUE.
            WRITE (IOUT,FMT=6130) KEY2
            IF (ONLINE) WRITE (ITOUT,FMT=6130) KEY2
          END IF
C     
          GO TO 790
        END IF
C     
C---- Timeout delay (for X-windows) (seconds)
C     
      ELSE IF (KEY.EQ.'TIME') THEN
C     ************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
        WTIME = VALUE(2)
C     
C---- UNPACK img1 img2
C     
      ELSE IF (KEY.EQ.'UNPA') THEN
        IF (NTOK.EQ.3) THEN
C     ************************************
          CALL MKEYNM(2,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          IUN1 = VALUE(2)
          IUN2 = VALUE(3)
        ELSE IF (NTOK.EQ.2) THEN
C     ************************************
          CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          IUN1 = VALUE(2)
          IUN2 = IUN1
        ELSE
          WRITE(IOUT,FMT=7420)
          IF (ONLINE) WRITE(ITOUT,FMT=7420)
 7420     FORMAT(1X,'Must supply image numbers to be unpacked')
          STOP
        END IF
        UNPACK = .TRUE.
C     
C---- SAVE
C     
      ELSE IF (KEY.EQ.'SAVE') THEN
        IF (NTOK.EQ.2) THEN
          SAVENAM = LINE(IBEG(2):IEND(2))
        ELSE
          IF(LENSTR(WAXFN).GT.0)THEN
            SAVENAM =  WAXFN(1:LENSTR(WAXFN))//'_'//
     +           IMGNUM(1:LENSTR(IMGNUM))//'.sav'
          ELSE
            IF(LENSTR(IDENT).GT.2)THEN
              SAVENAM = IDENT(1:LENSTR(IDENT))//'.sav'
            ELSE
              SAVENAM = 'mosflm.inp'
            ENDIF
          ENDIF
        END IF
        IF (MODE.EQ.3) THEN
 960      WRITE(LINE,FMT=7600) SAVENAM(1:LENSTR(SAVENAM))
 7600     FORMAT('Name of save file (',A,') :')
          CALL MXDWIO(LINE, 0)
          CALL MXDRIO(LINE2)
C     
C---- Get filename using PARSER
C     
C     ******************************************
          CALL MPARSE(LINE2,IBEG,IEND,ITYP,VALUE,IDEC,NTOK)
C     ******************************************
          IF (NTOK.NE.0) THEN
            NCH = IEND(1) - IBEG(1) + 1
            IF (NCH.GT.70) THEN
              WRITE(LINE,FMT=7630)
 7630         FORMAT(1X,'Too many characters (max 70)')
              CALL MXDWIO(LINE, 0)
              GOTO 960
            END IF
            SAVENAM = LINE2(IBEG(1):IEND(1))
          END IF
          LINE = 'This file can be used as input to a processing job.'
          CALL MXDWIO(LINE,3)
        ELSE
          WRITE(IOUT,FMT=7602) SAVENAM(1:LENSTR(SAVENAM))
          IF (ONLINE) WRITE(ITOUT,FMT=7602) SAVENAM(1:LENSTR(SAVENAM))
 7602     FORMAT(1X,'Current processing commands will be written',
     +         ' to file:',A)
        END IF
C     
        CALL SAVEINP
C     
C---- NOLP...Do NOT apply the Lorentz Polarisation corrections, so
C     output
C     intensities are raw intensities.
      ELSE IF (KEY.EQ.'NOLP') THEN
        NOLP = .TRUE.
C     
C---- NOBACK...Do NOT subtract background from summation integration
C     intensities.
      ELSE IF (KEY.EQ.'NOBA') THEN
        NOBACK = .TRUE.
C
C---- three or four-circle goniostat
C
      ELSE IF (KEY.EQ.'GONI')THEN
         GENGON = .TRUE.
        ICOUNT = 1
 890    ICOUNT = ICOUNT + 1 
        IF(ICOUNT.GT.NTOK)GOTO 895
        SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     **************
        CALL CCPUPC(SUBKEY)
ccccc
C     mechanical setting angles for CHI or KAPPA axis -usually
C     something like 54.74 or 50.0.
        IF(SUBKEY.EQ.'SET ') THEN
          ICOUNT = ICOUNT + 1
          SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
          CALL CCPUPC(SUBKEY)
          IF(SUBKEY.EQ.'CHI ') THEN 
            ICOUNT = ICOUNT + 1 
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            GON_CHI_SET = VALUE(ICOUNT)
            IGON_CHI = 1
            IGON_CHI_SET = 1
C     .. this is probably a 3-circle like a Bruker if GON_CHI_SET is set.
            THREECIRCLE = .TRUE.
            IGON_KAPPA_SET = 0
            IF(GON_KAPPA_SET.GT.1E-3)THEN
              WRITE(IOUT,FMT=7385)GON_KAPPA_SET,GON_CHI_SET
              IF(ONLINE)WRITE(ITOUT,FMT=7385)GON_KAPPA_SET,
     $             GON_CHI_SET
 7385         FORMAT(/,3(' **** WARNING ****',/),' Previously set ',
     $             'Kappa setting angle of ',
     $             F7.3,' has been set to 0.00',/,' because the ',
     $             'Chi setting angle has been set to ',F7.3,
     $             ' with the GONIOSTAT SET ',
     $             'CHI',/,' keywords',/,'----',/,
     $             ' GONIOSTAT SET CHI and GONIOSTAT SET KAPPA ',
     $             'CANNOT be used in the same run')
              GON_KAPPA_SET = 0.0
            ENDIF
          ELSE IF(SUBKEY.EQ.'KAPP') THEN 
            ICOUNT = ICOUNT + 1 
            GON_KAPPA_SET = VALUE(ICOUNT)
            IGON_KAPPA_SET = 1
C     .. almost certainly a four circle if GON_KAPPA_SET 
            FOURCIRCLE = .TRUE.
            IGON_CHI_SET = 0
            IF(GON_CHI_SET.GT.1E-3)THEN
              WRITE(IOUT,FMT=7386)GON_CHI_SET,GON_KAPPA_SET
              IF(ONLINE)WRITE(ITOUT,FMT=7386)GON_CHI_SET,
     $             GON_KAPPA_SET
 7386         FORMAT(/,3(' **** WARNING ****',/),' Previously set ',
     $             'Chi setting angle of ',
     $             F7.3,' has been set to 0.00',/,' because the ',
     $             'Kappa setting angle has been set to ',F7.3,
     $             ' with the GONIOSTAT SET ',
     $             'KAPPA',/,' keywords',/,'----',/,
     $             ' GONIOSTAT SET KAPPA and GONIOSTAT SET CHI ',
     $             'CANNOT be used in the same run')
              GON_CHI_SET = 0.0
            ENDIF
          ENDIF

ccccc
        ELSE IF(SUBKEY.EQ.'CHI ') THEN 
          ICOUNT = ICOUNT + 1 
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          GON_CHI = VALUE(ICOUNT)
          IGON_CHI = 1
          IF(ABS(GON_KAPPA).GT.1E-3)THEN
            WRITE(IOUT,FMT=7382)GON_KAPPA,GON_CHI
            IF(ONLINE)WRITE(ITOUT,FMT=7382)GON_KAPPA,GON_CHI
 7382       FORMAT(/,3(' **** WARNING ****',/),' Previously set ',
     $           'Kappa value of ',
     $           F7.3,' has been set to 0.00 because Chi ',/,
     $           ' has been set to ',F7.3,' with the GONIOSTAT',
     $           ' CHI keywords',/,'----',/,' Kappa and Chi ',
     $           'CANNOT be used in the same run')
            GON_KAPPA = 0.0
            IGON_KAPPA = 0
          ENDIF
        ELSE IF(SUBKEY.EQ.'KAPP') THEN 
          ICOUNT = ICOUNT + 1 
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          GON_KAPPA = VALUE(ICOUNT)
          IGON_KAPPA = 1
          IF(ABS(GON_CHI).GT.1E-3)THEN
            WRITE(IOUT,FMT=7383)GON_CHI,GON_KAPPA
            IF(ONLINE)WRITE(ITOUT,FMT=7383)GON_CHI,GON_KAPPA
 7383       FORMAT(/,3(' **** WARNING ****',/),' Previously set ',
     $           'Chi value of ',
     $           F7.3,' has been set to 0.00 because Kappa ',/,
     $           ' has been set to ',F7.3,' with the GONIOSTAT ',
     $           'KAPPA keywords',/,'----',/,' Kappa and Chi ',
     $           'CANNOT be used in the same run')
            GON_CHI = 0.0
            IGON_CHI = 0
          ENDIF
          ELSE IF(SUBKEY.EQ.'OMEG') THEN 
          ICOUNT = ICOUNT + 1 
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          GON_OMEGAS = VALUE(ICOUNT)
          IGON_OMEGA = 1
C
C---- of course, gonio_2theta is the same as 'TWOTHETA'
C
        ELSE IF(SUBKEY.EQ.'TWOT') THEN 
          ICOUNT = ICOUNT + 1 
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          GON_2THETA = VALUE(ICOUNT)
          IGON_2THETA = 1
C
C---- specify which is the oscillation axis (should be omega or phi)
C
        ELSE IF(SUBKEY.EQ.'OSCA') THEN 
          ICOUNT = ICOUNT + 1 
C     *******************************************
          CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
          SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
          CALL CCPUPC(SUBKEY)
          IF ((SUBKEY.EQ.'OMEG').OR.(SUBKEY.EQ.'PHI '))THEN
            HOSCAXIS = SUBKEY
            IHOSCAXIS = 1
          ELSE
            WRITE(IOUT,FMT=7384)SUBKEY
            IF(ONLINE)WRITE(ITOUT,FMT=7384)SUBKEY
            BADLINE = .TRUE.
            GOTO 50
 7384       FORMAT(/,3(' **** WARNING ****',/),'Mosflm does not',
     $           'recognize oscillations about axes other than',
     $           ' PHI or OMEGA')
          ENDIF
        ELSE
C     
C     Not recognised
C     
          IF (TRAPERR) INPERR = .TRUE.
          WRITE (IOUT,FMT=6130) SUBKEY
          IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
        END IF
        IF (ICOUNT.LT.NTOK) GOTO 890
 895    CONTINUE
C
C---- use Rigaku active mask or not - on by default
C
      ELSE IF (KEY.EQ.'NOAC')THEN
        LACTIVE = .FALSE.
        REALLYFORCEREAD = .TRUE.
      ELSE IF (KEY.EQ.'ACTI')THEN
        LACTIVE = .TRUE.
        REALLYFORCEREAD = .TRUE.
C
C---- turn BEST output on or off
C     
      ELSE IF (KEY.EQ.'BEST')THEN
        ICOUNT = 1
 896    ICOUNT = ICOUNT + 1
        SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
        CALL CCPUPC(SUBKEY)
        IF(SUBKEY(1:2).EQ.'ON')THEN
          LBEST = .TRUE.
        ELSEIF(SUBKEY(1:3).EQ.'OFF')THEN
          IF(LBEST.AND.BESTOPEN)THEN
            CLOSE(BESTHKL)
            BESTOPEN = .FALSE.
          ENDIF
          LBEST = .FALSE.
        ENDIF
        IF(SUBKEY.EQ.'FILE')THEN
C     ..
C     .. by implication, if we name bestfile.hkl, LBEST should be true.
          IF(LBEST.AND.BESTOPEN)THEN
            CLOSE(BESTHKL)
            BESTOPEN = .FALSE.
          ENDIF
          LBEST = .TRUE.
          ICOUNT = ICOUNT + 1
          IF(ITYP(ICOUNT).EQ.1)THEN
            BSTHKL = LINE(IBEG(ICOUNT):IEND(ICOUNT))
          ELSE
            WRITE(IOUT,FMT=7570)
            IF(ONLINE)WRITE(ITOUT,FMT=7570)
 7570         FORMAT(1X,30('*'),/,' Error in name of BESTFILE.HKL')
          ENDIF
        ENDIF
        IF(ICOUNT.LT.NTOK) GOTO 896
C
C---- HEAD  read only head or tail of image (for new GUI)
C
        ELSE IF (KEY.EQ.'HEAD')THEN
C
C---- Need to set NULLPIXTEST to FALSE here for iMosflm
C
          NULLPIXTEST = .FALSE.
          HEADONLY = .TRUE.
          PHIONLY = .FALSE.
          NO_XML = .FALSE.
          NODISPLAY = .TRUE.
          DISPMENU = .FALSE.
          POWDER = .FALSE.
C
C---- option to output limited information only - PHI angles at present
C
        IF(NTOK.GT.1)THEN
          SUBKEY = LINE(IBEG(2):IEND(2))
          CALL CCPUPC(SUBKEY)
          IF(SUBKEY.EQ.'BRIE')PHIONLY = .TRUE.
        ENDIF
        USEHDR = .TRUE.
CAL31/5/05          IMGP = .TRUE.
        HDRSIZE = .TRUE.
C
        SUMPART = .FALSE.
C
C---- LOAD read only head or tail of image (for new GUI) - output no XML
C
        ELSE IF (KEY.EQ.'LOAD')THEN
          HEADONLY = .TRUE.
          PHIONLY = .FALSE.
          NO_XML = .TRUE.
          NODISPLAY = .TRUE.
          DISPMENU = .FALSE.
          POWDER = .FALSE.
C
C---- PICK keyword (for new GUI)
C
        ELSE IF (KEY.EQ.'PICK')THEN
          WRITE(*,*) ' IERR is', IERR
          BADPICK = .FALSE.
          IF(NTOK .NE. 6)THEN
            WRITE(IOUT,FMT=9100)
            IF(ONLINE)WRITE(ITOUT,FMT=9100)
            BADLINE = .TRUE.
            GOTO 50
          ELSE
            DO ICOUNT = 2,6
              IF(ITYP(ICOUNT).NE.2)BADPICK = .TRUE.
            ENDDO
            IF(BADPICK)THEN
              WRITE(IOUT,FMT=9100)
              IF(ONLINE)WRITE(ITOUT,FMT=9100)
              BADLINE = .TRUE.
              GOTO 50
            ELSE
              NAUTO = 1
              NRUN = 1 
              ID = VALUE(2)
              NSER = ID
CHRP 15.05.2013              IDAUTO(1) = ID
              IDPACK(1) = ID
              NOIMG(1) = ID
              IDPACK(1) = ID
c              IFIRSTPACK = 1
              NODISPLAY = .TRUE.
              DOFIND = .FALSE.
              DO ICOUNT = 3,6
                PICKER(ICOUNT-2) = INT(VALUE(ICOUNT))
              ENDDO
              PICKLOG = .TRUE.
              POWDER =.TRUE.
              IF(IERR.NE.0)THEN
                WRITE(IOUT,FMT=9100)
                IF(ONLINE)WRITE(ITOUT,FMT=9100)
              ENDIF
            ENDIF
          ENDIF
          GOTO 50
 9100     FORMAT(/,'Usage: PICK <img_no> <x> <y> <x_size> <y_size>')
C
C---- FITCIRCLE keyword (for new GUI) - can have up to 20 points!
C
        ELSE IF (KEY.EQ.'FITC')THEN
          DO ICOUNT=1,200
            IXYCPT(1,ICOUNT) = 0
            IXYCPT(2,ICOUNT) = 0
          ENDDO
          BADFIT = .FALSE.
          IF(ITYP(2).EQ.2)THEN
            NCIRPT = VALUE(2)
            MCIRPT = NCIRPT
            IF(NCIRPT.LT.3)BADFIT = .TRUE.
            IF((NTOK-2)/2.NE.NCIRPT.AND.MOD(NTOK-2,2).NE.NCIRPT)
     $           BADFIT = .TRUE.
          ELSE
            BADFIT = .TRUE.
          ENDIF
          IF(BADFIT)THEN
            WRITE(IOUT,FMT=9300)
            IF(ONLINE)WRITE(ITOUT,FMT=9300)
            BADLINE = .TRUE.
            GOTO 50
          ENDIF
          DO ICOUNT=4,NTOK,2
            IXYCPT(1,(ICOUNT/2)-1) = VALUE(ICOUNT - 1)
            IXYCPT(2,(ICOUNT/2)-1) = VALUE(ICOUNT)
          ENDDO
          IMODE = 5
          CALL FITCIR(NCIRPT,MCIRPT,IXYCPT,CIRCEN,CIRRAD,RMSRES,IMODE)
          IF(IMODE.EQ.-999)THEN
            WRITE(IOUT,FMT=9340)
            IF(ONLINE)WRITE(ITOUT,FMT=9340)
          ELSE
            WRITE(IOUT,FMT=9320)CIRCEN,CIRRAD,RMSRES
            IF(ONLINE)WRITE(ITOUT,FMT=9320)CIRCEN,CIRRAD,RMSRES
          ENDIF
c          IF(SOCKLO)THEN            
            XMLLINE = ' '
C            WRITE(iout,FMT=9360)NCIRPT
            WRITE(XMLLINE,FMT=9360)NCIRPT
            CALL WRITE_SOCKET_SECTION(SERVERFD,LENSTR(XMLLINE),XMLLINE)
            DO ICOUNT = 1,NCIRPT,1
              XMLLINE = ' '
C              WRITE(iout,FMT=9380)
C     $             ICOUNT,IXYCPT(1,ICOUNT),IXYCPT(2,ICOUNT)
              WRITE(XMLLINE,FMT=9380)
     $             ICOUNT,IXYCPT(1,ICOUNT),IXYCPT(2,ICOUNT)
              CALL WRITE_SOCKET_SECTION(SERVERFD,LENSTR(XMLLINE),
     $             XMLLINE)
            ENDDO
            XMLLINE = ' '
C            WRITE(iout,FMT=9400)CIRCEN,CIRRAD,RMSRES
            WRITE(XMLLINE,FMT=9400)CIRCEN,CIRRAD,RMSRES
            CALL WRITE_SOCKET_LENGTH(SERVERFD,LENSTR(XMLLINE),XMLLINE)
c          ENDIF
          GOTO 50
 9300     FORMAT(/,' Usage: FITCIRCLE <nPts> <x1> <y1> <x2> <y2> ... ',
     $         '<xn> <yn>',/,' You must specify at least 3 points!')
 9320     FORMAT(/,' centre is ',F7.2,'mm',F7.2,'mm, Radius is ',F7.2,
     $         'mm and the RMS residual is ',F5.2,'mm')
 9340     FORMAT(/,' These points do not give a circle ')
 9360     FORMAT('<?xml version="1.0"?>',
     $         '<!DOCTYPE circle_fitting_response>',
     $         '<circle_fitting_response><number_of_points>',I3,
     $         '</number_of_points><list_of_points>')
 9380     FORMAT('<point>',I3,'</point><x>',I5,'</x><y>',I5,'</y>')
 9400     FORMAT('</list_of_points>',
     $         '<circle_centre_x>',F8.2,'</circle_centre_x>',
     $         '<circle_centre_y>',F8.2,'</circle_centre_y>',
     $         '<radius>',F8.2,'</radius>',
     $         '<rms_deviation>',F5.2,'</rms_deviation>',
     $         '</circle_fitting_response>')
C
C---- SEGMENT keyword (for new GUI) - to work out the appropriate minimum
C     sensible segments for postrefinement
C
        ELSE IF (KEY.EQ.'SEGM')THEN
          BADSEGM = .FALSE.
          IF(ITYP(2).NE.2)BADSEGM = .TRUE.
          IF(ITYP(3).NE.2)BADSEGM = .TRUE.
          ICOUNT = 3
          IF(BADSEGM)THEN
            WRITE(IOUT,FMT=9420)
            IF(ONLINE)WRITE(ITOUT,FMT=9420)
            BADLINE = .TRUE.
            GOTO 50
          ENDIF
C     ..
C     .. check for angles - if we have START then we must have ANGLE
          IF(NTOK.GT.3)THEN
            ICOUNT = 4
 5000       IF(ITYP(ICOUNT).EQ.1)THEN
              SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
              CALL CCPUPC(SUBKEY)
              IF(SUBKEY.EQ.'STAR')THEN
                ICOUNT = ICOUNT + 1
                IF(ITYP(ICOUNT).EQ.2)THEN
                  PHIBEG = VALUE(ICOUNT)
                ELSE
                  WRITE(IOUT,FMT=9440)
                  IF(ONLINE)WRITE(ITOUT,FMT=9440)
 9440             FORMAT(/,' **** ERROR ****',/,' If given, these',
     $                 ' keywords must follow the following pattern:',
     $                 /,/,'     SEGMENT <first> <last>',/,' OR',/,
     $                 '     SEGMENT <first> <last> (START <x> ',
     $                 'ANGLE <y>)')
                  BADLINE = .TRUE.
                  GOTO 50
                ENDIF
              ELSEIF(SUBKEY.EQ.'ANGL')THEN
                ICOUNT = ICOUNT + 1
                IF(ITYP(ICOUNT).EQ.2)THEN
                  PHIRNG = VALUE(ICOUNT)
                ELSE
                  WRITE(IOUT,FMT=9440)
                  IF(ONLINE)WRITE(ITOUT,FMT=9440)
                  BADLINE = .TRUE.
                  GOTO 50
                ENDIF
              ENDIF
              ICOUNT = ICOUNT + 1
              IF(ICOUNT.LT.NTOK)GOTO 5000
            ENDIF
            PHIEND = PHIBEG + PHIRNG
          ENDIF
          FIRSTIMG = INT(VALUE(2))
          LASTIMG = INT(VALUE(3))
          GOODSEG = .TRUE.
          CALL SEGMENT(JSEG,NUMSPG,FIRSTIMG,LASTIMG,IERR)
c          GOODSEG = .TRUE.
          IF(IERR.NE.0)THEN
            WRITE(IOUT,FMT=9420)
            IF(ONLINE)WRITE(ITOUT,FMT=9420)
c            BADLINE = .TRUE.
            GOTO 50
          ENDIF
 9420     FORMAT(/,' Usage: SEGMENT <first> <last> (START <first> ',
     $           'ANGLE <first>)')
C     
C---- BELL turn bell on or off
C
      ELSE IF (KEY.EQ.'BELL')THEN
        LBELL = .NOT.LBELL
        ICOUNT = 2
        SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
        CALL CCPUPC(SUBKEY)
        IF(SUBKEY(1:2).EQ.'ON')LBELL = .TRUE.
        IF(SUBKEY(1:3).EQ.'OFF')LBELL = .FALSE.
C     
C---- SIMUlated spot list - this spot list is the one used for autoindexing, and
C     this keyword writes a file from a single image using the predicted spot
C     positions. Therefore, it needs a matrix first!
C
      ELSE IF (KEY.EQ.'SIMU')THEN
        CALL KWSIMU(LINE,IBEG,IEND,NPARM,IERR)
C
C---- NUSPOT switch for new style SPOTOD writing/reading
C
      ELSE IF (KEY.EQ.'NUSP')THEN
        ICOUNT = 2
        SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
        CALL CCPUPC(SUBKEY)
        IF(SUBKEY(1:2).EQ.'ON')THEN
          NUSPOT = .TRUE.
          IF(.NOT.NEWPREF)THEN
            NEWPREF = .TRUE.
            WRITE(IOUT,FMT=7387)
            IF(ONLINE)WRITE(ITOUT,FMT=7387)
 7387       FORMAT (/,' ***** Information *****',/,
     &           ' NUSPOT ON requires the NEW-STYLE  ',
     &           'the NEW-STYLE POSTREFINEMENT (POSTREF MULTI),',
     &           ' so it has been turned on.')
          ENDIF
        ENDIF
        IF(SUBKEY(1:3).EQ.'OFF')NUSPOT = .FALSE.
C
C---- NUBACK switch for new style background calculation in spotfinding
C
      ELSE IF (KEY.EQ.'NUBA')THEN
        ICOUNT = 2
        SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
        CALL CCPUPC(SUBKEY)
        IF(SUBKEY(1:2).EQ.'ON')NUBACK = .TRUE.
        IF(SUBKEY(1:3).EQ.'OFF')NUBACK = .FALSE.
C
C---- NUFIND switch for new style spotfinding - not implemented yet!
C
      ELSE IF (KEY.EQ.'NUFI')THEN
        ICOUNT = 2
        SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
        CALL CCPUPC(SUBKEY)
        IF(SUBKEY(1:2).EQ.'ON')NUFIND = .TRUE.
        IF(SUBKEY(1:3).EQ.'OFF')NUFIND = .FALSE.
C
C---- HARRY keyword for prototype code
C
      ELSE IF (KEY.EQ.'HARR')THEN
        ICOUNT = 2
        SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
        CALL CCPUPC(SUBKEY)
        IF(SUBKEY(1:2).EQ.'ON')HARRY = .TRUE.
        IF(SUBKEY(1:3).EQ.'OFF')HARRY = .FALSE.
C
C---- MINPKSIZ minimum spot size in pixels
C
      ELSE IF (KEY.EQ.'MINP') THEN     
C     ************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
        MINPKSIZ = VALUE(2)
        IF (MINPKSIZ.LT.1) MINPKSIZ = 1
C
C---- LINE2D keyword for a line across an image, to give a 2D plot of
C     intensity
C
      ELSE IF (KEY.EQ.'LINE')THEN
        IF(NTOK.NE.5)THEN
          WRITE(IOUT,FMT=50000)LINE(IEND(1)+1:LENSTR(LINE))
          IF(ONLINE)WRITE(ITOUT,FMT=50000)LINE(IEND(1)+1:LENSTR(LINE))
50000     FORMAT(/,' Usage: LINE2D <x1> <y1> <x2> <y2>, not ',/,
     $         '        LINE2D ',A)
        ELSE
          LINEPOINT(1) = 0
          LINEPOINT(2) = 0
          LINEPOINT(3) = 0
          LINEPOINT(4) = 0
          DO II = 2,5
            IF(ITYP(II).EQ.2)THEN
              LINEPOINT(II -1) = VALUE(II)
            ELSE
              WRITE(IOUT,FMT=50010)II-1,LINE(IBEG(II):IEND(II))
              IF(ONLINE)WRITE(ITOUT,FMT=50010)II-1,
     $             LINE(IBEG(II):IEND(II))
50010         FORMAT(/,' For the command LINE2D, ',
     $             ' co-ordinate ',I1,' should be a number, but ',
     $             'it is ''',A,'''')
              GOTO 50
            ENDIF
          ENDDO
          CALL LINE2D(LINEPOINT)
        ENDIF
C
C---- MTZNEW adds columns giving measurement box parameters to MTZ file
C
      ELSE IF (KEY6.EQ.'MTZNEW')THEN
        NEWMTZ = .TRUE.
C
C---- VALIdate cell with symmetry
C
      ELSE IF (KEY.EQ.'VALI')THEN
        IFLAG = 0
        CALL VALIDATE(IFLAG)
        GOTO 50
C
C---- MERGE ON/OFF  will turn on/off merging of prediction lists when
C     there are multiple lattices
C
      ELSE IF (KEY.EQ.'MERG')THEN
        IF (NTOK.EQ.1) THEN
          MERGELATT = .TRUE.
        ELSE
          ICOUNT = 2
          SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
          CALL CCPUPC(SUBKEY)
          IF(SUBKEY(1:2).EQ.'ON') MERGELATT = .TRUE.
          IF(SUBKEY(1:3).EQ.'OFF') MERGELATT = .FALSE.
        END IF
C
C---- LATTICE <N> defines which lattice is to be processed when multiple
C     lattices are present.
C     Subkeyword OVERLAP defines which lattices should be included when
C     allowing for multiple lattices, there are IOVWORK of these and their
C     lattice numbers are stored in IOVLATTICE
C
C     LATTICE DELETE <n> means delete (inactivate) lattice <n>
C     LATTICE SINGLETON means only write singletons to the output MTZ file
C     
C
      ELSE IF (KEY.EQ.'LATT')THEN
C
C---- Test for string DELETE or SINGLETON in second position
C
        IF ((ITYP(2).EQ.1)) THEN
          SUBKEY = LINE(IBEG(2):IEND(2))
          CALL CCPUPC(SUBKEY)
          IF (SUBKEY.EQ.'DELE') THEN
C     ************************************
            CALL MKEYNM(1,3,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
            I = VALUE(3)
            ACTIVELATT(I) = .FALSE.
            WRITE(IOUT,FMT=7782) I
            IF (ONLINE) WRITE(ITOUT,FMT=7782) I
 7782       FORMAT(/,1X,'Deleting lattice',I3)
          ELSE IF (SUBKEY.EQ.'SING') THEN
            SINGLETON = .TRUE.
            WRITE(IOUT,FMT=7784)
            IF (ONLINE) WRITE(ITOUT,FMT=7784)
 7784       FORMAT(/,1X,'Only singletons will be written to the MTZ',
     +               ' file.')
          END IF
        ELSE

C     ************************************
          CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     
          ILATTICE = VALUE(2)
          IOVWORK = 0
C
C---- Following 2 lines needed for batch integration when no indexing done.
C
          NLATTICE = MAX(NLATTICE,ILATTICE)
          ACTIVELATT(NLATTICE) = .TRUE.
C
C---- If MATRIX keyword has been read, set CRYSTALA to AMAT. I'm not sure that
C     this is actually used, seem to use CRYSTALUBREF instead now.
C
          IF (IMAT.GT.0) THEN
            DO I = 1,3
              DO J = 1,3
                CRYSTALA(I,J,ILATTICE) = AMAT(I,J)
              ENDDO
            ENDDO
          ENDIF
          IF (NTOK.GE.4) THEN
            IOVWORK = 0
            ICOUNT = 3
            SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
            CALL CCPUPC(SUBKEY)
            IF (SUBKEY.EQ.'OVER') THEN
 900           ICOUNT = ICOUNT + 1
C            ************************************
               CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C            ************************************
               IOVWORK = IOVWORK + 1
               IOVLATTICE(IOVWORK) = NINT(VALUE(ICOUNT))
C
C---- Following 2 lines needed for batch integration when no indexing done.
C
               NLATTICE = MAX(NLATTICE,IOVLATTICE(IOVWORK))
               ACTIVELATT(IOVLATTICE(IOVWORK)) = .TRUE.
               IF (ICOUNT.LT.NTOK) GOTO 900
               WRITE(IOUT,FMT=7780) ILATTICE,(IOVLATTICE(I),I=1,IOVWORK)
               IF (ONLINE) WRITE(ITOUT,FMT=7780) ILATTICE,
     +                     (IOVLATTICE(I),I=1,IOVWORK)
            END IF
          ELSE
C
C---- If only LATTICE <N1> is given, set all other lattices as
C     overlapping
C
            DO 902 I = 1,NLATTICE
              IF (I.NE.ILATTICE) THEN
                IOVWORK = IOVWORK + 1
                IOVLATTICE(I) = I
              END IF
 902        CONTINUE
            IF (IOVWORK.GT.1) THEN
              WRITE(IOUT,FMT=7780) ILATTICE,(IOVLATTICE(I),I=1,IOVWORK)
              IF (ONLINE) WRITE(ITOUT,FMT=7780) ILATTICE,
     +                     (IOVLATTICE(I),I=1,IOVWORK)
            ELSE
              WRITE(IOUT,FMT=7781) ILATTICE
              IF (ONLINE) WRITE(ITOUT,FMT=7781) ILATTICE
            END IF 
 7780       FORMAT(/,1X,'Dealing with lattice',I2,' allowing for ',
     +              'presence of lattice(s):',10I2)
 7781       FORMAT(/,1X,'Dealing with lattice',I2)
          END IF
        END IF
C     
C---- CLEAN reduce influence of non-predicted neighbouring spots
C
      ELSE IF (KEY.EQ.'CLEA')THEN
        CLEAN = .TRUE.
        IF (NTOK.GT.1) THEN
          ICOUNT = 2
          SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
          CALL CCPUPC(SUBKEY)
          IF (SUBKEY(1:2).EQ.'ON') THEN
            CLEAN = .TRUE.
          ELSE IF(SUBKEY(1:3).EQ.'OFF') THEN
            CLEAN = .FALSE.
          END IF
        END IF
C     
C---- HARVEST turn harvesting on
C
      ELSE IF (KEY.EQ.'HARV')THEN
        ICOUNT = 2
        SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
        CALL CCPUPC(SUBKEY)
        IF(SUBKEY(1:2).EQ.'ON')THEN
          DOHARVEST = .TRUE.
          CALL MHARVEST(0)
        ELSEIF(SUBKEY(1:3).EQ.'OFF')THEN
          DOHARVEST = .FALSE.
        ELSE
          WRITE(IOUT,FMT=9000)
          IF(ONLINE)WRITE(ITOUT,FMT=9000)
 9000     FORMAT(/,' **** ERROR ****',/,
     $         ' Usage: HARVEST ON or HARVEST OFF',/)
        ENDIF
C     
C---- NOTVERTICAL Giving this keyword will override the change in
C     polarisation correction for PILATUS detectors with OMEGA=90 or 270,
C     that applies a polarisation correction for a vertical rotation axis.
C
      ELSE IF (KEY.EQ.'NOTV')THEN
        NOTVERTICAL = .TRUE.
C     
C     SCALE factor for internal scaling of I's and SDs.
C     
      ELSE IF (KEY6.EQ.'SCALEU') THEN
C     *******************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
        SCALEUP = VALUE(2)
C     
C     SCALE factor for output intensities and sigmas.
C     
      ELSE IF (KEY.EQ.'SCAL') THEN
C     *******************************************
        CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
        SCALEINT = VALUE(2)
C
C---- ONLINE overrides result of CCPONL as this seems broken for OSX,
C           returning FALSE when it should be TRUE
C
      ELSE IF (KEY.EQ.'ONLI') THEN
         ONLINE = .TRUE.
C
C---- OFFLINE overrides result of CCPONL (or a previous ONLINE keyword)
C
      ELSE IF (KEY.EQ.'OFFL') THEN
         ONLINE = .FALSE.
C
C---- PHIPROFILE calculate a phi profile for the reflection with iMosflm
C     pixel coordinates IMOSFLMX, IMOSFLMY on image IMOSFLMIMG, which
C     has start and end phi values of PROFPHIS and PROFPHIE. The first
C     and last images in the dataset are IMOSFIRSTIMG,IMOSLASTIMG.
C     IMGPAD sets the number of images to be added to the integration range
C     at the start and end of the integration range determined by the width
C     of the reflection (to avoid having to repeat integration for other 
C     reflections).
C     PHIPRIFILE REPEAT uses the same parameters as a previous call.
C     As this command is only sent by iMosflm, subkeywords are not used
C     except for REPEAT, RESTRICT and UNRESTRICT. If RESTRICT is present
C     then the resolution limits are reset to values very slightly above 
C     and below the resolution of the reflection. UNRESTRICT will reset
C     resolution limts to the original values.
C
      ELSE IF (KEY.EQ.'PHIP') THEN
        IF (FIRSTPHIP) THEN
C
C---- Save the original resolution limits if this is the first call to PHIPROFILE
C
          RESSAVE = RES
          RESLOWSAVE = RESLOW
          FIRSTPHIP = .FALSE.
         END IF
C
C---- Need to test if the previous PHIPROFILE command had the keyword RESTRICT
C     and if so, reset the resolution limits to the original ones, otherwise
C     it will not be able to find the new reflection.
c     *** This did not work as RESTRESL set false in REFLFIND after profile has been 
C     extracted from MODARRAY, try setting REPREDICT in RESTRESL instead.
C
         IF (RESTRESL) THEN
           REPREDICT = .TRUE.
CAL           DSTMIN = WAVE/RESLOWSAVE
CAL           DSTMAX = WAVE/RESSAVE
         END IF
C
C---- check if next subkeyword is REPEAT, if true, use values
C     from the previous PHIPROFILE command. If text is not REPEAT
C     ignore completely.
C
        IF (ITYP(2).EQ.1) THEN
          KEY2 = LINE(IBEG(2):IEND(2))
          CALL CCPUPC(KEY2)
          IF (KEY2(1:4).NE.'REPE') GOTO 50
        ELSE
C     ************************************
          CALL MKEYNM(7,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
          IMOSFLMIMG = VALUE(2)
          PROFPHIS = value(3)
          PROFPHIE = value(4)
          IMOSFLMX = VALUE(5)
          IMOSFLMY = VALUE(6)
          IMOSFIRSTIMG = VALUE(7)
          IMOSLASTIMG = VALUE(8)
          IMGPAD = VALUE(9)
        END IF
C
C---- Interpret RESTrict or UNREstrict keyword which determines if
C     resolution limits are to be restricted to just above and below
C     the resolution of this reflection, or go back to original resolution
C     limits. There are 3 tokens if 'PHIPROFILE REPEAT RESTRICT/UNRESTRICT' 
C     is sent, otherwise 9.
C
        IF ((NTOK.EQ.3).OR.(NTOK.EQ.9)) THEN
          KEY2 = LINE(IBEG(NTOK):IEND(NTOK))
          CALL CCPUPC(KEY2)
          IF (KEY2(1:4).EQ.'REST') THEN
C
C---- Restrict resolution limits
C
            RESTRESL = .TRUE.
          ELSE IF (KEY2(1:4).EQ.'UNRE') THEN
C
C---- Restore original resolution limits
C
            RES = RESSAVE
            RESLOW = RESLOWSAVE
            RESTRESL = .FALSE.
          END IF
        END IF
        MODEPHI = 1
        PROFWIDTH = 0.0
C
C---- Print for program flow
C
        CALLED = 'REFLFIND'
        COMMENT = 'to find reflection indices from pixel coords'
        IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
C
        CALL REFLFIND(IMOSFLMX,IMOSFLMY,IMOSFLMIMG,PROFPHIS,PROFPHIE,
     +                 IMOSFIRSTIMG,IMOSLASTIMG,PROFWIDTH,
     +                PHIREF,PHIWREF,IPACKF,IPACKL,NSER,
     +                NSERRUN,IFIRSTPACK,MODEPHI,IERR)
C
C---- If MODEPHI has been set to 3 in REFLFIND, the images required for the
C     phi profile already been integrated and we do not want to repeat the 
C     integration run, so skip the next keyword (which is "GO"). 
C     Also skip the GO if the reflection has not been found at all, (IERR = -1)
C     or if the reflection is too wide (too many images) when IERR = -2
        IF ((MODEPHI.EQ.3).OR.(IERR.EQ.-1).OR.(IERR.EQ.-2)) THEN
C
C---- first write <done> to socket
C
          LINE200 = '<done>'
C
C---- Print for program flow
C  
          CALLED = 'WRITE_SOCKET_LENGTH'
          COMMENT = 'To write <done> before skipping blank line'
          IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
          CALL WRITE_SOCKET_LENGTH(SERVERFD,6,LINE200)
c     socket
          LINE200 = ' '
C
C---- Print for program flow
C  
          CALLED = ' *** READ_SOCKET ***'
          COMMENT = 'to skip the blank line'
          IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
C
          CALL READ_SOCKET(SERVERFD,LINE400, ntok)

C
C---- Now repeat
C
C
C---- first write <done> to socket
C
          LINE200 = '<done>'
C
C---- Print for program flow
C  
          CALLED = 'WRITE_SOCKET_LENGTH'
          COMMENT = 'To write <done> before skipping a GO keyword'
          IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
          CALL WRITE_SOCKET_LENGTH(SERVERFD,6,LINE200)
c     socket
          LINE200 = ' '
C
C---- Print for program flow
C  
          CALLED = ' *** READ_SOCKET ***'
          COMMENT = 'to skip the GO keyword'
          IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
C
          CALL READ_SOCKET(SERVERFD,LINE400, ntok)


C Socket is broken
          IF(NTOK.EQ.-2)THEN
            NSHUTERR = 54
            CALL SHUTDOWN('CONTROL',NSHUTERR)
          ENDIF
          IF (DEBUG(89)) THEN
            WRITE(IOUT,FMT=7790) LINE400(1:120)
            IF (ONLINE) WRITE(ITOUT,FMT=7790) LINE400(1:120)
          END IF
        END IF
        IF (IERR.NE.0) GOTO 50
C
C---- FORCEREAD Used to force reading of image when the INVERTX value is changed in iMosflm GUI 
C     (Experiment settings). 
C
      ELSE IF (KEY.EQ.'FORC') THEN
        INVXFORCEREAD = .TRUE.
C
C----- *** end of keywords except RUN/GO *** 
C     
C---- RUN or GO
C     
      ELSE IF (KEY.EQ.'RUN ' .OR. KEY.EQ.'GO  ') THEN
C
C---- If getting data for a phi profile, set SOCKLO false so no XML is generated.
C
        IF (PHIPROF) THEN
          SOCKLO = .FALSE.
          SOCKLOPHIPROF = .TRUE.
        END IF
C
C---- If this is not a phi profile run (ie a normal integration or other operation)
C     set ISTARTLAST,IENDLAST (the first and last images used in a phi profile
C     integration run) to -999
C
        IF (.NOT.PHIPROF) THEN
          ISTARTLAST = -999
          IEND LAST = -999
        END IF
C
C---- Give a warning that this keyword is no longer recommended if
C     not being used with new gui - NOT here hrp 04.10.2012
C
        IF (BADFILE) THEN
          IF (SOCKLO) THEN
            LINE200 = ' '
            WRITE(LINE200,FMT=7870)
 7870       FORMAT('iMosflm does not currently handle',
     +      ' NeXus (.nxs) files.__nl__Please use the',
     +      ' master.h5 file instead')
            WRITE(LONGLINE,FMT=7872) LINE200(1:LENSTR(LINE200))
            CALL WRITE_SOCKET_LENGTH(SERVERFD,
     $        LENSTR(LONGLINE),LONGLINE)
 7872       FORMAT('<?xml version="1.0"?>',
     $             '<!DOCTYPE trapped_error_response>',
     $             '<trapped_error_response>',
     $             '<type>Wrong image file type</type>',
     $             '<message>',A,'</message>',
     $             '</trapped_error_response>')
             GOTO 50
          ELSE
            WRITE(IOUT,FMT=7874)
            IF (ONLINE) WRITE(ITOUT,FMT=7874)
 7874       FORMAT(/1X,'*** FATAL ERROR ***',/,1X,'MOSFLM does',
     +      ' not currently handle NeXus file format.',/,1X,
     +      'Please use the "master.h5" file instead')
            STOP
          END IF
        END IF
        IF (DEBUG(89)) THEN
          WRITE(IOUT,FMT=7731) OMEGAFD
          IF (.NOT.ONLINE) WRITE(IOUT,FMT=7731) OMEGAFD
 7731     FORMAT(1X,'OMEGAFD=',F9.2)
        END IF
        IF (.not.newgui.and..not.headonly.and..not.findspot
     $       .and..not.autoindx.and..not.winopen.and..not.socklo
     $       .and.dispmenu) THEN
          WRITE(IOUT,FMT=7730)
          IF (ONLINE) WRITE(ITOUT,FMT=7730)
        END IF
 7730   FORMAT(/,/,/,8(1X,'*********************************',
     +       '********************************************',/),
     +       1X,'Please do not use the old mosflm GUI, use imosflm',
     +       ' instead.',/,8(1X,'*********************************',
     +       '********************************************',/))

        IF(RESIDMAX.LT.0.001)
     $       RESIDMAX = RSDMAX*2.0*(ETA + 0.5*(DIVH+DIVV))/DTOR 
        DPSDONE = SDPSINDEX
        SDPSINDEX = DPSINDEX
C
C---- CHRP 070605 need to zero spot list for indexing here - 
C
        IF(DOFIND)THEN
chrp 15.05.2013          SNAUTO = NAUTO
cagwl 12.06.13 Put this line back to get multixtal indexing working
          SNAUTO = NAUTO
          DO 10000 I=1,NSPOTS
            XSPT(I) = 0.0
            YSPT(I) = 0.0
10000     ENDDO
          DO 10010 I = 1,MAXIMG
            ISTIMG(I) = 0.0
            IENDIMG(I) = 0.0
            SELECT(I) = .FALSE.
10010     ENDDO
        ENDIF
        GOTO 950
C     
C---- End
C     
      ELSE IF ((KEY.EQ.'END').OR.(KEY.EQ.'EXIT').OR.(KEY.EQ.'STOP')
     +       .OR.(KEY.EQ.'QUIT')) THEN
        IF (MODE.EQ.3) RETURN
        IF (MODE.EQ.10) THEN
          IF (DISPMENU) THEN
            POWDER = .TRUE.
            IF (WINOPEN) CALL MXDCIO(1,0,0,0,0)
            RETURN
          END IF
        END IF
       CALL CTIMER(ENDTIME)
c       WRITE(IOUT,FMT=7800)ENDTIME-STARTTIME
c       IF(ONLINE)WRITE(ITOUT,FMT=7800)ENDTIME-STARTTIME
 7800  FORMAT(/,'This run of Mosflm took ',F9.3,' seconds',/)
        WRITE (IOUT,FMT=6196)
 6196   FORMAT (///,1X,
     $       '*********** END OF PROCESSING *****************')
        IF (ONLINE) WRITE (ITOUT,FMT=6196)
        call dna_end_output
C     
C     ******************
c     another case of failing to reset the genopen flag...
        IF (GENOPEN) then
C
C---- Print for program flow
C
          CALLED = 'QCLOSE'
          COMMENT = 'To close generate file (if open) after END keyword'
          IF (DEBUG(99)) 
     +       CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
          CALL QCLOSE(IUNIT)
          GENOPEN = .FALSE.
        END IF
C     ******************
        IF (MTZOPEN) THEN
          MTZPRT = 1
C     *********************
          if (debug(52))then
            write(iout,fmt=7746)
            if (online) write(itout,fmt=7746)
 7746       format(1x,'**closing MTZ file in control 4 ***')
          end if
          CALLED = 'LWCLOS'
          COMMENT = 'To close MTZ file '
          IF (DEBUG(99)) 
     +       CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
          CALL LWCLOS(MTZOUT,MTZPRT)
          MTZOPEN = .FALSE.
          IF(LBEST.AND.BESTOPEN)THEN
            CLOSE(BESTHKL)
            BESTOPEN = .FALSE.
          ENDIF
C     *********************
        END IF
C     
c     -harvest
C     
C        IF(HARVESTREADY)
       IF(DOHARVEST)CALL MHARVEST(1)
c     -harvest

c     close dna stuff here

c     this could be bad

       call dna_end_output
       if(dnaout) close(dnafd)
        STOP
C     
C---- ***** Not recognised *****
C     
      ELSE
        IF (TRAPERR) INPERR = .TRUE.
        IF (ONLINE) WRITE (ITOUT,FMT=6198)
        BADKEYWORD = .TRUE.
 6198   FORMAT (2X,'****** Keyword NOT RECOGNISED ***** ')
        WRITE (IOUT,FMT=6198)
C     
C---- End of control cards
C     
      END IF
C     
C---- Read next control card
C     
      GO TO 50        
C     
C---- Trap inadvertent entering of GO when doing keyword input from
C     MXDSPL
C     
 950  CONTINUE
      IF (MODE.EQ.3) RETURN
C     
C---- If simply unpacking images, do it now.
C     
      IF (UNPACK) THEN
        IF (IIDENT.EQ.0) THEN
          WRITE(IOUT,FMT=7430)
          IF (ONLINE) WRITE(IOUT,FMT=7430)
 7430     FORMAT(1X,'Must supply an IDENT keyword for unpacking')
          STOP
        END IF
        SUMPART = .FALSE.
        LIPROCESS = 4
C     
        DO 930 ID = IUN1,IUN2
C     
C-----OPENODS to UNPACK  images
C     
C     ********************************************************
          CALL OPENODS(IDENT,ID,NFIRSTF,ODEXT,FDISK,MODEOP,
     +         PACK,ODFILE,SEPCHAR,FORCEREAD,IPACK,
     +         TEMPLSTART,TEMPLEND,LIPROCESS,IERR)
C     *******************************************************
C---- Write out unpacked image
C     
          CALL IMGOUT(ID,ODFILE)
 930    CONTINUE
        STOP
      END IF
C
C---- only opening an image to read the header (or tailer)
C     
      IF(HEADONLY)THEN
        MODEOP = 2
        LIPROCESS = 0
C
C---- Print for program flow
C
        CALLED = 'OPENODS'
        COMMENT = 'Reading header only, MODEOP=2, LIPROCESS=0'
        IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
        CALL OPENODS(IDENT,ID,NFIRSTF,ODEXT,FDISK,MODEOP,
     +       PACK,ODFILE,SEPCHAR,FORCEREAD,IPACK,
     +       TEMPLSTART,TEMPLEND,LIPROCESS,IERR)
C HRP 08.12.2005
C     SECTION should _NOT_ be called here if there's been some kind of error 
C     reading the image, but Geoff has written 
C     his code so that he needs to follow this route.That's wrong, but it's 
C     much easier not to argue, and in the short term anyway, it works.
c HRP uncommented ierr test 10.12.2009
        IF(IERR.EQ.0)THEN
C
C---- Print for program flow
C
          CALLED = 'SECTION'
          WRITE(COMMENT,FMT=7860) HEADONLY, NO_XML
 7860     FORMAT(1X,'HEADONLY=',L1,'  NO_XML=',L1)
          IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)

          CALL SECTION(RSCAN,RSCANRED,RMINXINP,NIMAGES,ID,IPIXY,
     $         IXOFFSET,IYOFFSET,DISPSET,ICCX,CCXRESET,SITE,SCANNER)
        ENDIF
        HEADONLY = .FALSE.
        IMGKWD = 0
        GOTO 50
      ENDIF
C     
C---- If doing a multi segment post-refinement, only allowed one PROCESS
C     keyword per "RUN" keyword
C     
      IF (MULTISEG.AND.(NPROCRUN.GT.1)) THEN
        WRITE(IOUT,FMT=7350) NPROCRUN
        IF (ONLINE) WRITE(ITOUT,FMT=7350) NPROCRUN
 7350   FORMAT(1X,'***** FATAL ERROR *****',/,1X,
     +       'In a multi-segment post-refinement run, you must only',
     +       ' give a single "PROCESS"',/,1X,'keyword for each ',
     +       'RUN keyword.',/,1X,'For this RUN there are',I2,' PROC',
     +       'ESS keywords.')
        NSHUTERR = 44
        CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
      END IF
      NPROCRUN = 0
C     
C---- If this is a TESTGEN run, check that a SCANNER keyword has been
C     given
C     
      IF ((TESTGEN.OR.STRATEGY).AND.(ISCAN.EQ.0)) THEN
        WRITE(IOUT,FMT=6476)
        IF (ONLINE) WRITE(ITOUT,FMT=6476)
 6476   FORMAT(1X,'*** ERROR ***',/,1X,'A SCANNER keyword must be ',
     +       'given for TESTGEN and STATEGY options',/,1X,
     +       'Possibilities are MAR, SMALLMAR, RAXIS, DIP2000,',
     +       'DIP3000, FUJI, MD, ADSC, JUPITER')
        IF (COMREAD) THEN
          COMREAD = .FALSE.
          ITIN = ITINS
          CLOSE (UNIT=ICOMM)
        END IF
        BADLINE = .TRUE.
        IF (STRATEGY.AND.(NSTRUN.GT.1)) STOP
        GOTO 50
      END IF
C     
C---- If phistart and/or osc angle have NOT been set on PROCESS keyword,
C     check that they can be read from image header.
C     
C     AL       HEADINFO = ((((MACHINE.EQ.'MAR ').OR.(MACHINE.EQ.'RAXI'))
C     .AND.
C     AL     +               (NHEAD.EQ.1).AND.(.NOT.PACK)).OR.
C     AL     +              ((MACHINE.EQ.'DIP2').AND.(NTAIL.EQ.1)))
      HEADINFO = ((USEHDR.OR.USETAIL).AND.USEPHI)
      IF ((IPROKWD.GT.0).AND.(.NOT.HEADINFO)) THEN
        IF ((ISTRT.EQ.0).OR.(IANGLE.EQ.0)) THEN
          IF (ISTRT.EQ.0) THEN
            WRITE(IOUT,FMT=7197)
            IF (ONLINE) WRITE(ITOUT,FMT=7197)
          END IF
          IF ((IANGLE.EQ.0).AND.(ISTRT.NE.0)) THEN
            WRITE(IOUT,FMT=7197)
            IF (ONLINE) WRITE(ITOUT,FMT=7197)
          END IF
 7197     FORMAT(1X,' *** ERROR ***',/,1X,'With this scanner it ',
     +         'is not possible to read phi values from the image ',
     +         'header.',/,1X,'You must therefore give both',
     +         ' START and ANGLE subkeywords on the PROCESS keyword.',
     +         /,1X,'eg PROCESS 1 TO 20 ',
     +         'START 4 ANGLE 1.0')
C     
          STOPRUN = .TRUE.
        END IF
      END IF
C     
C---- If repeating a multisegment post-refinement from scratch, reset
C     cell
C     parameters to saved values, and set missets to the refined values
C     for the first image (stored in DELPHIV).
C     
      IF (RPTFIRST) THEN
        DO 688 I = 1,6
          CELL(I) = SAVECELL(I)
 688    CONTINUE
        ICELL = 1
        DO 691 I = 1,3
          DELPHI(I) = DELPHIV(I)
 691    CONTINUE
      END IF
C 
      if(.not.dpsdone)then
        NRUN = NRUN + 1
      else
        nrun = 1
      endif
c     hrp07122001      NRUN = NRUN + 1
      IF ((NRUN.GT.1).AND.(.NOT.MULTISEG)) THEN
        SUMPART = SSUMPART
        POSTREF = SPOSTREF
        ADDPART = SADDPART
      END IF
      IF (DEBUG(52)) THEN
        WRITE(IOUT,FMT=6570) ITIN,NRUN,MULTISEG,RPTFIRST,COMREAD,
     +       NSER,NSEG,ISTRUN,NSTRUN,SUMPART,
     +       ADDPART,POSTREF,NPACKS,NADD
        IF (ONLINE) WRITE(ITOUT,FMT=6570) ITIN,NRUN,MULTISEG,
     +       RPTFIRST,COMREAD,NSER,NSEG,ISTRUN,NSTRUN,
     +       SUMPART,ADDPART,POSTREF,NPACKS,NADD
      END IF
 6570 FORMAT(1X,'Read a RUN card, ITIN=',I3,' NRUN is now',I3,
     +     ' MULTISEG',L2,' RPTFIRST',L2,' COMREAD',L2,/,1X,
     +     'NSER=',I3,' NSEG=',I3,' ISTRUN=',I3,' NSTRUN=',I3,
     +     ' SUMPART',L2,' ADDPART',L2,' POSTREF',L2,/,1X,'NPACKS',
     +     I5,'  NADD',I5)
C     
C---- If input keywords were being read from a file, change input
C     channel back to normal so that raster box parameters can be
C     supplied by terminal                                
C  
      IF (COMREAD.and..not.dpsdone.and..not.dofind) ITIN = ITINS
C     
C---- If this is a multi-segment postrefinement run, read the data for
C     all the segments at this point,store each line and determine how
C     many packs are to be included in the refinement.
C     
      IF (MULTISEG.AND.(NSEG.EQ.1)) THEN
        NADD = IPACK2A(1) - IPACK1A(1)
        NRUN = 1
        NRLINE = NTLINE
        WAITINP = .FALSE.
      END IF
C     
C---- *** Read in additional lines for multi run STRATEGY or POSTREF
C     SEGMENT
C     
      IF (MULTISEG.AND.(NSER.LT.NSEG).OR.
     +     (STRATEGY.AND.(NRUN.LT.NSTRUN))) THEN

C     
C---- If this is a repeat from scratch run, all additional lines have
C     already been read and stored so don't need to do this again.
C     Initialise variables
C
        IF (RPTFIRST) THEN
          NSER = NSEG
          NRUN = 1
          NLINE = NRLINE
          GOTO 696
        END IF
C     
C---- Also, if this is after the first (or later) segment has been
C     processed
C     then again, all lines have been stored, so skip rest of this
C     
        IF (NRUN.GT.1) GOTO 696
C     
C     
        NRLINE = NTLINE
C
C---- Because we are now  using new program flow (CONTROL not called 
C     for repeating entire cycles of cell refinement) only need to 
C     store these additional keywords
C
        NTLINE = 1
        NRLINE = 1

        NRUN = 1
        NADD = IPACK2A(1) - IPACK1A(1)
 689    IF (ONLINE.AND.(.NOT.COMREAD)) THEN
          IF (MULTISEG) WRITE(ITOUT,FMT=6540) NSEG-1
 6540     FORMAT(1X,'Please supply PROCESS, RUN (IDENT  if desired)',
     +         ' keywords for the remaining',I3,' segments',/,1X,
     +         'to be used in the post-refinement')
          IF (STRATEGY) THEN
            WRITE(ITOUT,FMT=6478) NSTRUN - 1
 6478       FORMAT(1X,'Please supply remaining keywords',
     +           ' for the remaining',I3,' parts')
            IF (MODE.EQ.10) THEN
              WRITE(LINE,FMT=7490) NSTRUN - 1
 7490         FORMAT(1X,'Please supply remaining keywords',
     +             ' for the remaining',I3,' parts')
              CALL MXDWIO(LINE,2)
            END IF
          END IF
          
        END IF
C     
C---- If previously reading from a file, continue reading from that file
C     
        IF (COMREAD) ITIN = ICOMM
C     
C---- Read next line of input
C     
 690    CONTINUE 
        IF (ONLINE.AND..NOT.SOCKLO) WRITE (ITOUT,FMT=6000)
        WRITE (IOUT,FMT=6000)
        IF (BRIEF) WRITE (IBRIEF,FMT=6000)
C
C---- If using new GUI, read from socket (12/09/05)
C
        IF(SOCKLO.AND..NOT.XMLOPEN)THEN
          CALL READ_SOCKET(SERVERFD,LINE400,NTOK)
C Socket is broken
          IF(NTOK.EQ.-2)THEN
            NSHUTERR = 55
            CALL SHUTDOWN('CONTROL',NSHUTERR)
          ENDIF
          IF (DEBUG(89)) THEN
            WRITE(IOUT,FMT=7810) LINE400(1:120)
            IF (ONLINE) WRITE(ITOUT,FMT=7810) LINE400(1:120)
 7810       FORMAT(1X,'Socket read in CONTROL2:',A)
          END IF
          LINE = LINE400
        ELSE
          IF (MODE.EQ.10) THEN           
            WRITE(LINE,FMT=7470)
            CALL MXDWIO(LINE,0)
            CALL MXDRIO(LINE)
            NCH = LENSTR(LINE)
          ELSE
            READ(ITIN,FMT=6542,END=694) LINE
 6542       FORMAT(A)
          END IF
        ENDIF
C
        IF (ONLINE) WRITE(ITOUT,FMT=6542) LINE(1:LENSTR(LINE))
        WRITE(IOUT,FMT=6542) LINE(1:LENSTR(LINE))
C     
        INLINE(NTLINE) = LINE
        NTLINE = NTLINE + 1
        IF(NTLINE.GT.NTOTLINE)NTOTLINE = NTLINE
        IF (DEBUG(52)) THEN
          WRITE(IOUT,FMT=7650) NTLINE-1, 
     +                 INLINE(NTLINE-1)(1:LENSTR(INLINE(NTLINE-1)))
          IF (ONLINE) WRITE(ITOUT,FMT=7650) NTLINE-1,
     +                 INLINE(NTLINE-1)(1:LENSTR(INLINE(NTLINE-1)))
        END IF
        IF (NTLINE.GT.1000) THEN
          WRITE(IOUT,FMT=6544)
          IF (ONLINE) WRITE(ITOUT,FMT=6544)
 6544     FORMAT(//,1X,'**** ERROR ***',/,1X,'More than 1000 lines',
     +         ' of input to MOSFLM before eof')
          IF (MODE.EQ.10) THEN
            LINE = 'Fatal error, see terminal window'
            CALL MXDWIO(LINE,1)
          END IF
          STOP
        END IF
C     
C---- Decode this line to see if contains a PROCESS or RUN card
C     
C     
C     ******************************************
        CALL MPARSE(LINE,IBEG,IEND,ITYP,VALUE,IDEC,NTOK)
C     ******************************************
C     
C---- Trap blank lines
C     
        IF (NTOK.EQ.0) THEN
          NTLINE = NTLINE - 1
          GOTO 690
        END IF
C     
C---- first 4 chars
C     
        KEY = LINE(IBEG(1) :IEND(1))
C     
C---- convert to upper case
C     
C     ***********
        CALL CCPUPC(KEY)
C     ***********
        IF ((KEY.EQ.'RUN').OR.(KEY.EQ.'GO')) THEN
          IF(RESIDMAX.LT.0.001)
     $         RESIDMAX = RSDMAX*2.0*(ETA + 0.5*(DIVH+DIVV))/DTOR 
          NRUN = NRUN + 1
          IF (NPROCRUN.GT.1) THEN
            WRITE(IOUT,FMT=7350) NPROCRUN
            IF (ONLINE) WRITE(ITOUT,FMT=7350) NPROCRUN
            NSHUTERR = 45
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          END IF
C     
C---- call AUTOMATCH if mosaicity estimation is required
C     
c     hrp 07022002
          IF(IFSTRAT)then
            NRUN = NRUN - 1
            ifstrat = .false.
          endif
          IF (MULTISEG.AND.(NSER.LT.NRUN)) THEN
            WRITE(IOUT,FMT=6546) 
            IF (ONLINE) WRITE(ITOUT,FMT=6546)
 6546       FORMAT(/,1X,'In a multi-segment postrefinement run a ',
     +           '"RUN" keyword has been given without giving',/,1X,
     +           'the PROCESS keyword')
            STOP
          END IF
c     -harvest
c     cxEBI can now choose to stop if either
c     cxEBI PNAMEgiven or DNAMEgiven are false ??
c     cxEBI
C     AL        IF (.not. PNAMEgiven ) THEN
C     AL          WRITE (IOUT,FMT=8001)
C     AL 8001     FORMAT(' Error No PROTEIN NAME GIVEN by KeyWord PNAME'
C     )
C     AL          IF (ONLINE) THEN
C     AL            WRITE (IOUT,FMT=8002)
C     AL 8002       FORMAT(' Enter PNAME now !')
C     AL            GO TO 50
C     AL          END IF
C     AL          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
C     AL        END IF
C     AL        IF (.not. DNAMEgiven ) THEN
C     AL          WRITE (IOUT,FMT=8003)
C     AL 8003     FORMAT(' Error No DATA SET NAME GIVEN by KeyWord DNAME'
C     
C     )
C     AL          IF (ONLINE) THEN
C     AL            WRITE (IOUT,FMT=8004)
C     AL 8004       FORMAT(' Enter DNAME now !')
C     AL            GO TO 50
C     AL          END IF
C     AL          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
C     AL        END IF
c     -harvest
          NPROCRUN = 0
C     
        ELSE IF (KEY.EQ.'PROC') THEN
          NOGO = .FALSE.
          DPSINDEX = .FALSE.
          ISTRT2 = 0
          IANGLE2 = 0
          NSER = NSER + 1
          NPROCRUN = NPROCRUN + 1
C     
C---- First check if the "TO" has been specified
C     
          IF (ITYP(3).EQ.2) THEN
C     ************************************
            CALL MKEYNM(2,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
            ICOUNT = 3
            IPACKF = NINT(VALUE(2))
            IPACKL = NINT(VALUE(3))
          ELSE
C     ************************************
            CALL MKEYNM(1,2,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
C     ************************************
            CALL MKEYNM(1,4,LINE,IBEG,IEND,ITYP,NTOK)
C     ************************************
            ICOUNT = 4
            IPACKF = NINT(VALUE(2))
            IPACKL = NINT(VALUE(4))
          END IF
C     
C---- Check for START, OSC/ANGLE, BLOCK or FILM keywords
C     
          IF (ICOUNT.EQ.NTOK) GOTO 832
C     
 830      ICOUNT = ICOUNT + 1
          SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     **************
          CALL CCPUPC(SUBKEY)
C     **************
          IF (SUBKEY.EQ.'STAR') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            PHISTART = VALUE(ICOUNT)
            ISTRT2 = 1
          ELSE IF ((SUBKEY(1:3).EQ.'OSC').OR.
     +           (SUBKEY(1:3).EQ.'ANG')) THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
C     
            PHIRNG = VALUE(ICOUNT)
            IANGLE2 = 1
          ELSE IF (SUBKEY(1:3).EQ.'ADD') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          ELSE IF (SUBKEY.EQ.'BLOC') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
          ELSE
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
          END IF
C     
          IF (ICOUNT.LT.NTOK) GOTO 830
C     
C---- Check both START and OSC/ANGLE have been given
C     
 832      IF (((ISTRT2.EQ.0).OR.(IANGLE2.EQ.0)).AND.(.NOT.HEADINFO))
     +         THEN
            IF (ISTRT2.EQ.0) THEN
              WRITE(IOUT,FMT=7197)
              IF (ONLINE) WRITE(ITOUT,FMT=7197)
            END IF
            IF ((IANGLE2.EQ.0).AND.(ISTRT2.NE.0)) THEN
              WRITE(IOUT,FMT=7197)
              IF (ONLINE) WRITE(ITOUT,FMT=7197)
            END IF
            STOP
          END IF
C     
          IPACK1A(NSER) = IPACKF
          IPACK2A(NSER) = IPACKL
          NADD = NADD + (IPACKL - IPACKF)
          IF (DEBUG(52)) THEN
            WRITE(IOUT,FMT=6548) NTLINE-1,NSER,NRUN,NADD,NRLINE
            IF (ONLINE) WRITE(ITOUT,FMT=6548) NLINE,NSER,NRUN,NADD,
     +           NRLINE
          END IF
 6548     FORMAT(1X,'On line',I3,', found SERIAL keyword number',
     +         I3,', NRUN currently',I3,' NADD is',I3,' NRLINE is',I3)
C     
C---- Trap zero increment
C     
          IF (PHIRNG.LT.0.0) THEN
            WRITE(IOUT,FMT=7101)
            IF (ONLINE) WRITE(ITOUT,FMT=7101)
            IF (BRIEF) WRITE(IBRIEF,FMT=7101)
 7101       FORMAT(1X,'*** ERROR ***',/,1X,'The oscillation angle',
     +           ' per image must be positive',/,1X,'If the phi',
     +           ' axis is rotating in the opposite sense',/,1X,
     +           'this should be corrected by redefining the OMEGA',
     +           ' angle',/,1X,'(Currently ',F6.1,') to 180 + ',
     +           ' current value using SCANNER OMEGA',/,1X,
     +           'keywords, then use positive phi increment.')
            STOP
          END IF
C     
C---- 
C     
        ELSE IF (KEY.EQ.'STRA') THEN
           dpsdone = .false.
          NSEGRD = NSEGRD + 1
C     
C---- Need to check for a SPEEDUP subkeyword and if present ensure it
C     is the same as original value.
C     
          ICOUNT = 1
          INERR = .FALSE.
 880      CONTINUE 
          ICOUNT = ICOUNT + 1
          SUBKEY = LINE(IBEG(ICOUNT):IEND(ICOUNT))
C     **************
          CALL CCPUPC(SUBKEY)
C     **************
          IF (SUBKEY.EQ.'STAR') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            IF (IOERR) INERR = .TRUE.
          ELSE IF (SUBKEY.EQ.'END') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            IF (IOERR) INERR = .TRUE.
          ELSE IF (SUBKEY.EQ.'STEP') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            IF (IOERR) INERR = .TRUE.
          ELSE IF ((SUBKEY.EQ.'RUNS').OR.(SUBKEY.EQ.'PART')) THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            IF (IOERR) INERR = .TRUE.
          ELSE IF (SUBKEY.EQ.'ROTA') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            IF (IOERR) INERR = .TRUE.
          ELSE IF (SUBKEY.EQ.'SEGM') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            IF (IOERR) INERR = .TRUE.
          ELSE IF (SUBKEY.EQ.'SIZE') THEN
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            IF (IOERR) INERR = .TRUE.
          ELSE IF (SUBKEY.EQ.'ANOM') THEN
            AUTANOM = .TRUE.
            CONTINUE
          ELSE IF (SUBKEY.EQ.'AUTO') THEN
            CONTINUE
          ELSE IF ((SUBKEY.EQ.'SPEE').AND.(ICOUNT.LT.NTOK)) THEN
            INSPEED = 1
            ICOUNT = ICOUNT + 1
C     *******************************************
            CALL MKEYNM(1,ICOUNT,LINE,IBEG,IEND,ITYP,NTOK)
C     *******************************************
            IF (IOERR) INERR = .TRUE.
            VOLSCAL = VALUE(ICOUNT)
            IF (VOLSCAL.NE.OVOLSCAL) THEN
              IF (OVOLSCAL.NE.1) THEN
                WRITE(IOUT,FMT=7140) OVOLSCAL
                IF (ONLINE) WRITE(ITOUT,FMT=7140) OVOLSCAL
                IF (MODE.EQ.10) THEN
                  LINE = 
     $                 '** ERROR ** You cannot specify different'
                  CALL MXDWIO(LINE,1)
                  LINE = 'SPEEDUP factors for different runs.'
                  CALL MXDWIO(LINE,1)
                  LINE = 'The original value has been restored'
                  CALL MXDWIO(LINE,3)
                END IF
              ELSE
                WRITE(IOUT,FMT=7340) VOLSCAL
                IF (ONLINE) WRITE(ITOUT,FMT=7340) VOLSCAL
                OVOLSCAL = VOLSCAL
 7340           FORMAT(//,1X,'*** Speedup factor of ',F7.1,
     +               ' taken from this line.')
                IF (MODE.EQ.10) THEN
                  WRITE(LINE,FMT=7492) VOLSCAL
 7492             FORMAT(//,1X,'Speedup factor of ',F7.1,
     +                 'taken from this line.')
                  CALL MXDWIO(LINE,1)
                END IF             
              END IF             
            END IF
C     
C     
C     Not recognised
C     
          ELSE
            IF (TRAPERR) INPERR = .TRUE.
            WRITE (IOUT,FMT=6130) SUBKEY
            IF (ONLINE) WRITE (ITOUT,FMT=6130) SUBKEY
            INERR = .TRUE.
          END IF  
          IF (ICOUNT.LT.NTOK) GOTO 880
C     
C---- Trap an error in input
C     
          IF (INERR) THEN
            INERR = .FALSE.
            WRITE(IOUT,FMT=7380)
            IF (ONLINE) WRITE(ITOUT,FMT=7380)
 7380       FORMAT(1X,'*** Because of input error, this line has',
     +           ' been ignored, please give it again ***')
            IF (MODE.EQ.10) THEN
              LINE = ' '
              WRITE(LINE,FMT=7474)
              CALL MXDWIO(LINE,2)
            END IF
            NSEGRD = NSEGRD - 1
            NTLINE = NTLINE - 1
            GOTO 690
          END IF
        END IF
C     
C---- Read more control cards if required
C     
        IF (MULTISEG.AND.(NRUN.LT.NSEG)) GOTO 690
        IF (STRATEGY.AND.(NRUN.LT.NSTRUN)) GOTO 690
C     
C---- All required information has been read, reset NSER,NRUN,NLINE
C     
        NSER = 1
        
        NRUN = 1
        NLINE = NRLINE
        WAITINP = .FALSE.
        GOTO 696
C     
C---- EOF before reading all data
C     Need to differentiate between POSTREF SEGMENT and STRATEGY
C     
        
 694    IF (MULTISEG) THEN 
          WRITE(IOUT,FMT=6550) NSEG,NSER,NRUN
          IF (ONLINE) WRITE(ITOUT,FMT=6550) NSEG,NSER,NRUN
 6550     FORMAT(//,1X,'Post-refinement using',I3,' segments was ',
     +         'requested, but end-of-file was',/,1X,
     +         'encountered after reading',
     +         I3,' PROCESS/SERIAL keywords and',I3,' RUN keywords.',
     +         /,1X,'See User Guide for example input.')
        ELSE IF (STRATEGY) THEN
          WRITE(IOUT,FMT=6551) NSTRUN,NSEGRD,NRUN
          IF (ONLINE) WRITE(ITOUT,FMT=6551) NSTRUN,NSEGRD,NRUN
 6551     FORMAT(//,1X,'The STRATEGY option with',I2,' different',
     +         ' runs was requested but end-of-file',/,1X,
     +         'was encountered after reading',
     +         I3,' STRATEGY keywords and',I3,' RUN keywords',
     +         /,1X,'See Help file for example input')
          IF (MODE.EQ.10) THEN
            LINE = 'ERROR...see terminal window'
            CALL MXDWIO(LINE,1)
          END IF
        END IF
        IF (ONLINE.AND.COMREAD) THEN
          COMREAD = .FALSE.
          ITIN = ITINS
          GOTO 689
        END IF
        STOP
      END IF
C     
C---- Set up defaults if POWDER keyword given
C     
 696  continue
      IF (COMREAD.and..not.dpsdone.and..not.dofind) ITIN = ITINS
C     
C---- Trap zero length TITLE
C     
      IF (LENSTR(GTITLE).EQ.0) GTITLE = '.'
C     
C---- Trap case where AUTOINDEX keyword has been given but no PROCESS
C     keyword
C     and no images have been specified. However, if AUTOINDEX NOREFINE,
C     allow
C     to proceed.
C     
      IF ((AUTOINDX.AND.(NSER.EQ.0)).AND.(NAUTO.EQ.0).AND.
     $     (NSOL.NE.-999).AND..NOT.DOFIND) THEN
        IF(.NOT.INDNOREF)THEN
          WRITE (IOUT,FMT=7290)
          IF (ONLINE) WRITE (ITOUT,FMT=7290)
 7290     FORMAT(1X,'***** FATAL ERROR *****',/,1X,'An AUTOINDEX ',
     +         'keyword has been given without a PROCESS keyword,',
     +         /,1X,'and no IMAGES have been specified.')
          NSHUTERR = 46
          CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
        END IF
      END IF
      IF ((MOSEST.AND.(NSER.EQ.0)).AND.(NAUTO.EQ.0)) THEN
c     IF(.NOT.INDNOREF)THEN
        WRITE (IOUT,FMT=7291)
        IF (ONLINE) WRITE (ITOUT,FMT=7291)
 7291   FORMAT(1X,'***** FATAL ERROR *****',/,1X,'A MOSAIC ESTI',
     +       'MATE keyword combination has been given without a',
     $       'ny IMAGES',/,' being specified.')
          NSHUTERR = 47
        CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
c     END IF
      END IF
      IF (PICKLOG.AND.(NSER.EQ.0)) THEN
        WRITE (IOUT,FMT=72910)
        IF (ONLINE) WRITE (ITOUT,FMT=72910)
72910   FORMAT(1X,'***** FATAL ERROR *****',/,1X,'A PICK command ',
     $       'has been given without a',
     $       'ny IMAGES',/,' being specified.')
          NSHUTERR = 471
        IF(.NOT.SOCKLO)CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
      END IF

C     
      IF (POWDER) THEN
C
C
C---- Print for program flow
C
          CALLED = 'IFPWDR'
          COMMENT = 'in IF POWDER block'
          IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
          CALL IFPWDR(IFIRSTPACK,GENOPEN,IERR,IFAIL)
C
C---- Trap open error for spot file (IFAIL.LT.0)
C
          IF (SOCKLO.AND.IFAIL.LT.0) THEN
            XMLLINE = ' '
            WRITE(XMLLINE,FMT=7710)templsav(1:lenstr(templsav)),
     $           ifirstpack
            CALL WRITE_SOCKET_LENGTH(SERVERFD,LENSTR(XMLLINE),XMLLINE)
 7710       FORMAT('<?xml version="1.0"?>',
     $           '<!DOCTYPE spot_search_response>',
     $           '<spot_search_response>',
     $           '<status>',
     $           '<code>error</code>',
     $           '<message>Unable to open spot file</message>',
     $           '</status>',
     $           '<image_template>',A,'</image_template>',
     $           '<image_number>',I8,'</image_number>',
     $           '</spot_search_response>')
            GOTO 50
          END IF
CAL          write(*,*),'Return from IFPWDR, IERR ',ierr
c
c---- hmmm, next line has to be if(nser.eq.0)goto 50 for background indexing...
c
          IF(IERR.EQ.50)GOTO 50
          IF(IERR.EQ.920)GOTO 920
          IF (STRATEGY.OR.TESTGEN) GOTO 697
          RETURN
C     
C---- END of IF POWDER block
C     
      END IF
C     
C**** CHECK OSCGEN INPUT
C     
C---- Skip several of these checks if STRATEGY is being used.
C     
 920  IF (STRATEGY.OR.TESTGEN) GOTO 697
C     
C---- Check for IDENT or TEMPLATE keywords
C     
      IF ((IIDENT.EQ.0).AND.(.NOT.TEMPLATE).and.(.not.HDF5FILE)) THEN
        WRITE (IOUT,FMT=6716)
        IF (ONLINE) THEN
          WRITE (ITOUT,FMT=6716)
        END IF
 6716   FORMAT(/,1X,'****** ERROR ******',/,1X,'No IDENT or TEMPLATE',
     +       ' keywords have been given.',/,1X,
     +       '(used as a template to form image file names.)',/,/) 
        STOPRUN = .TRUE.
      END IF
C     
C---- Check for SERIAL keyword
C     
      IF (NSER.EQ.0) THEN
        WRITE (IOUT,FMT=6714)
        IF (ONLINE) THEN
          WRITE (ITOUT,FMT=6714)
        END IF
 6714   FORMAT(1X,'****** No PROCESS (SERIAL) or IMAGE keyword has',
     +       ' been given *****')
        STOPRUN = .TRUE.
      END IF
C     
C     
C---- If this multiseg run has been set up in MXDSPL all lines have been
C     stored
C     
 768  IF (MODE.EQ.4) THEN
        
        NSER = 1
        NRUN = 1
c hrp 28.11.2002 - should nrline be the same as ntline here??
c        nrline = ntline
        NLINE = NRLINE
      END IF
C     
C---- if an integration run from MXDSPL, set up NSERRUN
C     
      IF ((MODE.EQ.2))NSERRUN = 1
C     
C---- Check for RASTER keyword
C     
      IF (IRAST.EQ.0) THEN
        WRITE (IOUT,FMT=6718)
        IF (ONLINE) THEN
          WRITE (ITOUT,FMT=6718)
        END IF
 6718   FORMAT(/,1X,'**** Information ****',/,1X,'No RASTER ',
     +       'keyword has been given.',/,
     +       1X,'(Gives the starting parameters for the ',
     +       'measurement box).',/,1X,'Suitable ',
     +       'parameters will be determined automatically.',/)
C     AL          STOPRUN = .TRUE.
        IF ((IXOFFSET.EQ.0).AND.(IYOFFSET.EQ.0)) THEN
          RADX = ROTATED
          RADY = (.NOT.RADX)
        END IF
      END IF
C     
C     AL *******  START HERE IF CALLING AFTER AUTOINDEXING **
C     
C---- Check NEWPREF
C     
 770  IF (NEWPREF) SUMPART = .FALSE.
C     
C---- Check that synchrotron parameters set correctly
C     
      IF (ISYN.NE.0) THEN
C     
C---- If TOLERANCE not set, make default 0.03
C     
        IF (ITOL.EQ.0) TOL = 0.03
C     
        IF (IDELAMB.EQ.0 .AND. KEYPX.EQ.0) THEN
C     
C---- add default value here if DISPERSION hasn't been set but warn user
C     
C     anyway
C     
          IDELAMB = 1
          DELAMB = 0.0015
          WRITE (IOUT,FMT=7038)
          IF (ONLINE) THEN
            WRITE (ITOUT,FMT=7038)
            WRITE (ITOUT,FMT=7039)
          END IF
C     STOPRUN = .TRUE.
        END IF
      END IF
 7036 FORMAT (' WAVELENGTH card must be present with "SYNCH" card')
 7038 FORMAT (' DISPERSION card must be present with "SYNCH" card',/,
     $     '***** it has been set to a default value of 0.0015A but',
     $     ' this is probably wrong! *****')
 7039 FORMAT ('      You should add an appropriate value before',
     $     ' continuing')
 7040 FORMAT (' MOSAIC card must be present with "SYNCH" card')
C     
C---- If MONOCHROMATOR collimation specified, calculate the polarisation
C     factor assuming a graphite monochromator
C     
 697  IF (IMONO.EQ.1) THEN
        COS2TH = COS(2.0*ASIN(WAVE/(2*3.427)))
        TOR = (COS2TH - 1.0)/(COS2TH + 1.0)
      END IF
C     
C---- Assign crystal class based on spacegroup number, after checking
C     that spacegroup has been assigned
C     
      IF (LSYMM.EQ.0) THEN
        WRITE(IOUT,FMT=6480)
        IF (ONLINE) THEN
          WRITE(ITOUT,FMT=6480)
          IF (STRATEGY.OR.TESTGEN) THEN
            BADLINE = .TRUE.
            GOTO 50
          ENDIF
        END IF
        IF (STRATEGY.OR.TESTGEN) STOP
        STOPRUN = .TRUE.
      END IF
 6480 FORMAT(/,1X,'****** ERROR *****',/,1X,'No SYMMETRY keyword ',
     +     ' has been given.',/,1X,
     +     '(Gives space group name or number.)'/)
C     
C---- Assignment of crystal class moved from here to follow reading of 
c     SYMM keyword
C     
C     
C---- Check for presence of SEPARATION card
C     
      IF (ISEP.EQ.0) THEN
        IF (MODE.EQ.1) THEN
          WRITE (IOUT,FMT=7041)
          IF (ONLINE) WRITE (ITOUT,FMT=7041)
 7041     FORMAT (/,1X,'**** Information ****',/,1X,'No ',
     +         'SEPARATION keyword has been given, no spots will be ',
     +         'flagged as overlapping.',/,1X,'Suitable separation ',
     +         'parameters will be determined automatically',
     +         ' prior to integration.')
        ELSE
          WRITE (IOUT,FMT=7043)
          IF (ONLINE) WRITE (ITOUT,FMT=7043)
        END IF
 7043   FORMAT (/,1X,'**** Information ****',/,1X,'No SEPARATION ',
     +       'keyword has been given.',
     +       /,1X,'(Gives minimum spot separation before spots ',
     +       'are flagged as overlapping.',/,1X,'Suitable ',
     +       'parameters will be determined automatically.',/)
C     AL        STOPRUN = .TRUE.
      END IF
C     
C---- If autoindexing and not fixing cell, give message about status
C     of input cell
C     
      IF (AUTOINDX.AND.(.NOT.RFIXCELL)) THEN
        IF (CELLKEEP) THEN
          WRITE(IOUT,FMT=7400)
          IF (ONLINE) WRITE(ITOUT,FMT=7400)
        ELSE
          WRITE(IOUT,FMT=7402)
          IF (ONLINE) WRITE(ITOUT,FMT=7402)
        END IF
      END IF
 7400 FORMAT(/,1X,'**** Information ****',/,1X,
     +     '***** Information ****',/,1X,
     +     '***** Information ****',/,1X,
     +     'Because the KEEP',
     +     ' subkeyword has been given on the CELL keyword,',/,1X,
     +     'the input cell will override that determined from ',
     +     'autoindexing.',/,/)
 7402 FORMAT(/,1X,'**** Information ****',/,1X,
     +     '**** Information ****',/,1X,
     +     '**** Information ****',/,1X,'The cell derived',
     +     ' from autoindexing will override that given on the ',
     +     'CELL keyword.',/,1X,'To force the program to use ',
     +     'the input cell, add the keyword KEEP. eg:'/,1X,
     +     'CELL KEEP 74.2 74.2 35.1 90 90 90',/,/)
C     
C---- Check for BEAM keyword
C     


C     
C---- If requested (RESET on MISSET keyword) incorporate missetting
C     angles into the U- or A-matrix (which ever has been given)
C     
      IF (RESET) THEN
        CALL ROTMAT(DELPHI,WORK2,1)
        IF ((IMAT.EQ.0).AND.(IUMAT.EQ.1)) THEN
          DO 144 I = 1,3
            DO 142 J= 1,3
              WORK(I,J) = UMAT(I,J)
 142        CONTINUE
 144      CONTINUE
          CALL MATMUL3(UMAT,WORK2,WORK)
          call checku(umat)
        ELSE IF (IMAT.EQ.1) THEN
          DO 148 I = 1,3
            DO 146 J= 1,3
              WORK(I,J) = AMAT(I,J)
 146        CONTINUE
 148      CONTINUE
          CALL MATMUL3(AMAT,WORK2,WORK)
        END IF
        DO 149 I = 1,3
          DELPHI(I) = 0.0
 149    CONTINUE
      END IF
C     
C     First check that if a CELL keyword has been supplied, it obeys
C     the spacegroup symmetry
C     
      IF (ICELL.EQ.1) THEN
        IFLAG = 0
        CALL CELLCHK(ICRYST,CELL,IFLAG)
        IF ((IFLAG.NE.0).AND.(NUMSPG.GT.0)) THEN
          WRITE(IOUT,FMT=7048) NUMSPG,SPGNAM(1:LENSTR(SPGNAM)),CELL
          IF (ONLINE) WRITE(ITOUT,FMT=7048) NUMSPG,
     $         SPGNAM(1:LENSTR(SPGNAM)),CELL
 7048     FORMAT(//,1X,'*** ERROR ***',/,1X,'Space group number ',
     +         I4,' (',A,') is inconsistent with supplied cell',
     +         ' parameters',/,1X,'Cell:',6F10.3)
          STOPRUN = .TRUE.
          IF(SOCKLO)THEN
            LINE400 = ' '
            WRITE(LINE400,FMT=8020)SPGNAM(1:LENSTR(SPGNAM)),CELL
            CALL WRITE_SOCKET_LENGTH(SERVERFD,
     $           LENSTR(LINE400),LINE400)
 8020       FORMAT('<?xml version="1.0"?>',
     $           '<!DOCTYPE interface_input_response>',
     $           '<interface_input_response>',
     $           '<status>',
     $           '<code>error</code>',
     $           '</status>',
     $           '<message>Space group supplied is inconsistent ',
     $           'with the unit cell dimensions</message>',    
     $           '<symmetry>',A,'</symmetry>',
     $           '<cell>',                              
     $           '<a>',F9.4,'</a>',                              
     $           '<b>',F9.4,'</b>',                              
     $           '<c>',F9.4,'</c>',                              
     $           '<alpha>',F8.4,'</alpha>',                              
     $           '<beta>',F8.4,'</beta>',                              
     $           '<gamma>',F8.4,'</gamma>',
     $           '</cell>',
     $           '</interface_input_response>')
            
          ENDIF
        END IF
      END IF

C     
C---- If doing strategy run, calculate orientation matrix from cell
C     param
C     and U matrix if orientation matrix not given.  Otherwise calculate
C     cell
C     param from the orientation matrix. Need to do this now so cell is
C     available when opening  MTZ file.
C     
      IF (STRATEGY) THEN
        ICHECK = 1
        IF (ICELL.EQ.1) ICHECK = 0
C     
C     ************************
        CALL SETMAT(IMAT,IUMAT,ICELL,ICHECK)
C     ************************
C     
C---- Calculate default SPEEDUP from cell volume and symmetry
C     
C     
C---- Get cell volume
C     
C     ********************
        CALL CELLVOL(CELL,CVOL)
C     ********************
C     
C---- Do not calculate default speedup if this is the second or
C     subsequent part of a multipart run.
C     
        IF ((INSPEED.EQ.0).AND.(ISTRUN.EQ.0)) THEN
          IF (CVOL.GT.0) THEN
            IF (NSYMP.GT.0) VOLSCAL = CVOL*REAL(NSYMP)/DEFVOL
          END IF
          VOLSCAL = NINT(VOLSCAL)
          VOLSCAL = MAX(1.0,VOLSCAL)
          WRITE(IOUT,FMT=7550) VOLSCAL
          IF (ONLINE) WRITE(ITOUT,FMT=7550) VOLSCAL
 7550     FORMAT(1X,'Default speedup set to',F9.1)
          IF (MODE.EQ.10) THEN
            LINE = ' '
            WRITE(LINE,FMT=7568) VOLSCAL
 7568       FORMAT('SPEEDUP factor set to',F8.1)
            CALL MXDWIO(LINE,1)
          END IF
          CELLSCAL = VOLSCAL**0.333333
        END IF
C     
      END IF
C     
C---- For strategy, if this is a second or later run, the cell from the
C     orientation matrix will probably be slightly different. Impose the
C     cell from the first run on all subsequent runs.
C     
      IF (STRATEGY.AND.(ISTRUN.GE.1)) THEN
        IF (DEBUG(52)) THEN
          WRITE(IOUT,FMT=6482) SAVECELL,CELL
          IF (ONLINE) WRITE(ITOUT,FMT=6482) SAVECELL,CELL
 6482     FORMAT(1X,'Savecell',6F8.2,/,1X,'Cell',6F8.2,/,1X,
     +         'Redetermine AMAT using original cell')
        END IF
        DO 701 I = 1,6
          CELL(I) = SAVECELL(I)
 701    CONTINUE
        IMAT = 0
        ICELL = 1
        IUMAT = 1
C     ************************
        CALL SETMAT(IMAT,IUMAT,ICELL,ICHECK)
C     ************************
        IF (DEBUG(52)) THEN
          WRITE(IOUT,FMT=6484)((AMAT(I,J),J=1,3),I=1,3) 
          IF (ONLINE) WRITE(ITOUT,FMT=6484) ((AMAT(I,J),J=1,3),I=1,3) 
 6484     FORMAT(1X,3F12.6)
        END IF
C     
C---- Now reset IMAT, or a subsequent call for another STRATEGY run will
C     fail
C     
        IMAT = 1
        IUMAT = 0
      END IF
C     
C---- Set up symmetry operations for STRATEGY mode (NOTE no longer need an
C     MTZ file).
C     
      IF (STRATEGY.AND.FIRSTRAT) THEN
        FIRSTRAT = .FALSE.
        IF (NSTRUN.EQ.0) NSTRUN = 1
C     
C---- Need EPSILON for unique option
C
        IPRINT = 0
        CALL EPSLN(NSYM,NSYMP,RSYM,IPRINT)
        CALL CENTRIC(NSYM,RSYM,IPRINT)
C     
C---- Set up matrices to reduce reflections to asymmetric unit
C     
C     ******
        CALL ASUSET(SPGNAM,NUMSPG,PGNAME,NSYM,RSYM,NSYMP,NLAUE,
     +       DEBUG(52))
C     ******
      END IF
C********************************END OF OSCGEN CHECKS **************
C*******************************************************************
C     
C---- Set maximum reflection width in number of images, hard limit of
C     IPAD. However, if doing integration to get data for a phi profile
C     do not assign NWMAX here as it has been set in REFLFIND
C     

      IF ((PHIRNG.NE.0.0).AND.(.NOT.PHIPROF)) NWMAX = NINT(WMAX/PHIRNG)
      IPAD = 100
      NWMAX = MIN(NWMAX,IPAD)
C     
C---- Limit NWMAX to minimum of 2, otherwise will not integrate partials
C     under any circumstances
C     
      NWMAX = MAX(NWMAX,2)
C     
C---- If both ADDPART and POSTREF turned off, set SUMPART false
C     
      IF ((.NOT.ADDPART).AND.(.NOT.POSTREF)) SUMPART = .FALSE.
C     
C---- Save SUMPART,ADDPART,POSTREF for case of processing "left-over"
C     SERIAL
C     runs. Don't do this is this is a prediction call from MXDSPL
C     because SUMPART,ADDPART,POSTREF have already been save and set
C     false.
C     Similarly, don't do it for a STRATEGY call from MXDSPL
C     
      IF ((MODE.NE.1).AND.(MODE.NE.10)) THEN
        SSUMPART = SUMPART
        SADDPART = ADDPART
        SPOSTREF = POSTREF
      END IF
C     
C---- Skip this if getting PREDICTION from MXDSPL
C     
      IF (MODE.EQ.1) GOTO 780
C     
C---- Check genfile has been specified
C     
      IF (.NOT.(STRATEGY.OR.TESTGEN)) THEN
C     
C---- First MTZ file
C     
        IF ((IHKLOUT.EQ.0).AND.(MTZNAM.EQ.'HKLOUT')) THEN
C     
C---- First check if it has been assigned on the command line
C     
          FWORK = ' '
          CALL UGTENV('HKLOUT',FWORK)
C     
C---- Only reset MTZNAM if no HKLOUT keyword given and environment
C     variable
C     HKLOUT has been set
C     
          IF (FWORK(1:1).NE.' ') THEN
            MTZNAM = FWORK
          ELSE
C     
C---- Not set on command line, set up filename from identifier and first
C     image.
            WRITE(ABC,7320) IPACK1A(1)
 7320       FORMAT(I3.3)
            IF (LENSTR(IDENT).NE.0) 
     +           MTZNAM = IDENT(1:LENSTR(IDENT))//'_'//ABC//'.mtz'
            WRITE(IOUT,FMT=7322) MTZNAM
            IF (ONLINE) WRITE(ITOUT,FMT=7322) MTZNAM
 7322       FORMAT(1X,'**** Information ****',/,1X,
     +           'No MTZ filename given, so filename has',
     +           ' been set to: ',A)
          END IF
        END IF
C     
        IF (GENFILE(1:8).EQ.'________') THEN
          WRITE (IOUT,FMT=6160)
          IF (ONLINE) WRITE (ITOUT,FMT=6160)
 6160     FORMAT (//,1X,'**** Information ****',/,1X,'No GENFILE ',
     +         'keyword has been given.',/,
     +         1X,'It will be set to the MTZ filename with an ',
     +         'extension ".gen"')
C     
C---- Set GENFILE to same name as MTZ file but with extension ".gen"
C     
          I = LENSTR(MTZNAM)
          DO 778 J = I,1,-1
            IF (MTZNAM(J:J).EQ.'.') THEN
              K = J - 1
              GOTO 779
            END IF
 778      CONTINUE
          K = I
 779      GENFILE = MTZNAM(1:K)//'.gen'
          NEWGENF = .TRUE.
        END IF
      END IF
C     
C hrp 21012010 I think this needs to be NPACKS not NPACK here?
      NPRUN = NPACK + 1 - IFIRSTPACK

C     
C---- If only one image specified, can't do post-refinement or addpart
C     
      IF (NPRUN.EQ.1) THEN
        ADDPART = .FALSE.
        SUMPART = .FALSE.
        POSTREF = .FALSE.
        WRITE(IOUT,FMT=6490)
        IF (ONLINE) WRITE(ITOUT,FMT=6490)
        IF (BRIEF) WRITE(IBRIEF,FMT=6490)
      END IF
 6490 FORMAT(/,1X,'***** Because only one pack is to be processed',
     +     ', no post refinement or',/,1X,'addition of ',
     +     'partials will be carried out *****')
C     
C---- If using CLOSE option, ADDPART cannot be used.
C     
      IF (DENSE.AND.ADDPART) THEN
        ADDPART = .FALSE.
        WRITE(IOUT,FMT=6492)
        IF (ONLINE) WRITE(ITOUT,FMT=6492)
      END IF
 6492 FORMAT(//,1X,'***** WARNING *****',/,1X,'Because the ',
     +     'SEPARATION CLOSE option is',
     +     ' being used ADDPART has been turned off.')
      IF (NEWPREF.AND.ADDPART) THEN
        ADDPART = .FALSE.
        WRITE(IOUT,FMT=6493)
        IF (ONLINE) WRITE(ITOUT,FMT=6493)
      END IF
 6493 FORMAT(//,1X,'***** WARNING *****',/,1X,'Because the ',
     +     'new style post-refinement is',
     +     ' being used ADDPART has been turned off.',/,1X,
     +     'To use ADDPART you must give keywords POSTREF NOMULTI',
     +     ' which will restrict post-refinement',/,1X,'to ',
     +     'reflections spread over no more than 2 images.')
C     
C---- Open first image file to get image size (Mar, Raxis and Mac
C     Science only)
C     Do not do this if calling from MXDSPL
C     
      IF ((MODE.GT.0).AND.(MODE.LT.10)) GOTO 780
C     
C---- If IDENT has not been given, cannot open image file
C     
      IF ((IIDENT.EQ.0).AND.(.NOT.TEMPLATE)) GOTO 784
C     
      IF (IMGP.AND.HDRSIZE
     +     .AND.((.NOT.STRATEGY).AND.(.NOT.TESTGEN))) THEN
        MODEOP = 2
        ID1 = IDPACK(IFIRSTPACK)
        NFIRSTF = 1
        LIPROCESS = 0
C     
C---- Call OPENODS to get image size and header information when NOT in
C     POWDER mode. Note that for DIP images, the whole image is read in
C     order to get the tailer information
C     
C     Thus for DIP scanners we need to assign NWORD, NBYTE now
C     
        IF (MACHINE.EQ.'DIP2') THEN
          NWORD = IYLEN   
          NBYTE = IYLEN*2
        END IF
C
C---- Print for program flow
C
        CALLED = 'OPENODS'
        COMMENT = 'when not in POWDER mode LIPROCESS=0, MODEOP=2'
        IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)

C     ********************************************************
        CALL OPENODS(IDENT,ID1,NFIRSTF,ODEXT,FDISK,MODEOP,
     +       PACK,ODFILE,SEPCHAR,FORCEREAD,IPACK,
     +       TEMPLSTART,TEMPLEND,LIPROCESS,IERR)
C     *******************************************************
C     
cdms21.09.2004   If file does not exist, jump out. Also need to reset NIMAGES
C     to zero since it always uses the first image here
C
        IF (ID1.EQ.-999) THEN
           IF (ONLINE) THEN
              IF (COMREAD) THEN
                 COMREAD = .FALSE.
                 ITIN = ITINS
                 CLOSE (UNIT=ICOMM)
              END IF
              NIMAGES = 0
              NIMAG = 0
              BADLINE = .TRUE.
              GOTO 50
           END IF
           CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
        END IF
C
C---- If using header or tailer information, do checks now
C
        IF (USEHDR.OR.USETAIL) THEN
          IF ((IDIST.GT.0).AND.(ABS(0.01*XTOFD - HDIST).GT.0.1))
     +         THEN
            WRITE(IOUT,FMT=6713) 0.01*XTOFD, HDIST
            IF (ONLINE) WRITE(ITOUT,FMT=6713) 0.01*XTOFD, HDIST
          END IF
C     
C---- If distance not specified, set to value from header
C     
          IF ((IDIST.EQ.0).AND.(HDIST.NE.0)) THEN
            XTOFD = 100.0*HDIST
            CAMCON(3) = HDIST
            CALL SETDIS(ITILT,ITWIST,1)
c     RADEG = 18000.0/3.14159
c     IF (XTOFD.NE.0) FDIST = 1.0/ (XTOFD*RADEG)
c     TILT = ITILT*FDIST
c     TWIST = ITWIST*FDIST
            IDIST = 1
            WRITE(IOUT,FMT=6717) HDIST
            IF (ONLINE) WRITE(ITOUT,FMT=6717) HDIST
          END IF
C     
C---- Check wavelength
C     
          IF ((IWAVE.EQ.2).AND.(ABS(WAVE - HWAVE).GT.0.001)) THEN
            WRITE(IOUT,FMT=6715) WAVE,HWAVE
            IF (ONLINE) WRITE(ITOUT,FMT=6715) WAVE,HWAVE
          END IF
C     
C---- If wavelength not specified, set to value from header if non-zero
C     
          IF ((IWAVE.EQ.0).AND.(HWAVE.NE.0)) THEN
            WRITE(IOUT,FMT=6719) HWAVE
            IF (ONLINE) WRITE(ITOUT,FMT=6719) HWAVE
            WAVE = HWAVE
            IWAVE = 1
          END IF
C
C---- If CUTOFF not given by a keyword and there is a non-zero value
C     from the image header, set to the image header value. Important
C     for Eiger and Pilatus detectors where it depends on frame rate.
C
          IF ((ICUT.EQ.0).AND.(HCUTOFF.GT.0)) THEN
            WRITE(IOUT,FMT=6727) HCUTOFF
            IF (ONLINE) WRITE(ITOUT,FMT=6727) HCUTOFF
            CUTOFF = HCUTOFF
            IF (IPRCUT.EQ.0) THEN
              PRCUTOFF = CUTOFF
              IPRCUT = 1
              WRITE(IOUT,FMT=6730) HCUTOFF
              IF (ONLINE) WRITE(ITOUT,FMT=6730) HCUTOFF
            END IF
            ICUT = 1
          END IF

C     
C---- Check header pixel size
C     
          IF ((IPIX.EQ.0).AND.(HRAST.NE.0)) THEN
            WRITE(IOUT,FMT=6726) HRAST
            IF (ONLINE) WRITE(ITOUT,FMT=6726) HRAST
            RAST = HRAST
            IPIX = 1
          END IF
C     
C---- Check pixel size against input value (if given)
C     
          IF ((IPIX.EQ.1).AND.(ABS(RAST - HRAST).GT.0.001)) THEN
            WRITE(IOUT,FMT=6729) RAST,HRAST
            IF (ONLINE) WRITE(ITOUT,FMT=6729) RAST,HRAST
          END IF
C     
C---- If NULLPIX set in header (Mar CCD) and not set by keyword, use
C     the  value from header
C     
          IF ((.NOT.IINULL).AND.(HNULLPIX.GT.0)) THEN
            NULLPIX = HNULLPIX
            WRITE(IOUT,FMT=6728) NULLPIX
            IF (ONLINE) WRITE(ITOUT,FMT=6728) NULLPIX
          END IF        
C     
C---- Check that oscillation angle agrees
C     
          IF ((HPHIE-HPHIS).GT.0.0) THEN
            IF ((IANGLE.GT.0).AND.
     +           ((ABS(PHIRNG-(HPHIE-HPHIS)).GT.0.01))) THEN
              WRITE(IOUT,FMT=7230) PHIRNG,(HPHIE-HPHIS)
              IF (ONLINE) WRITE(ITOUT,FMT=7230)PHIRNG,(HPHIE-HPHIS)
 7230         FORMAT(4(1X,'***** WARNING *****',/),
     +             1X,'Input oscillation angle (',F6.2,
     +             ' degrees) does NOT agree with value in the ',
     +             'image header (',F6.2,' degrees)',/,1X,'The ',
     +             'input oscillation angle will be used.')

            END IF
          END IF
C     
C---- If oscillation angle not given, set to value from header, and set
C     up oscillation angles for all packs:
C     
          IF ((IANGLE.EQ.0).OR.(ISTRT.EQ.0)) THEN
            IF (IANGLE.EQ.0) PHIRNG = HPHIE - HPHIS
C     
C---- Check that oscillation angle from header is non-zero, if not give
C     a warning
C     
            IF (PHIRNG.EQ.0.0) THEN
              WRITE(IOUT,FMT=6725)
              IF (ONLINE) WRITE(ITOUT,FMT=6725)
            END IF
            IF (ISTRT.EQ.0) PHISTART = HPHIS
            J = 0
C     
C---- Note that PHIBEGA, PHIENDA are used in MAIN to set up
C     start and end oscillation angles.
C     
            DO 782 I = IFIRSTPACK,NPACK
              J = J + 1
              IF (I.EQ.IFIRSTPACK) THEN
                PHIBEGA(I) = PHISTART
              ELSE
                PHIBEGA(I) = ((J-1)*PHIRNG) + PHISTART
              END IF
              PHIENDA(I) = PHIBEGA(I) + PHIRNG
 782        CONTINUE
          END IF
        END IF
C
C---- End of "if using header or tailer information" block
C

        IF ((MODEL.EQ.'RAXISIV').OR.(MODEL.EQ.'RAXISHTC')) THEN
          IF (NREC.EQ.1500) THEN
            RAST = 0.2
          ELSE IF (NREC.EQ.3000) THEN
            RAST = 0.1
          ELSE IF (NREC.EQ.6000) THEN
            RAST = 0.05
          END IF
        END IF
        IF (MODEL.EQ.'RAXISV') THEN
          IF (NREC.EQ.2000) THEN
            RAST = 0.2
          ELSE IF (NREC.EQ.4000) THEN
            RAST = 0.1
          ELSE IF (NREC.EQ.8000) THEN
            RAST = 0.05
          END IF
        END IF
C
C---- END OF "IF (IMGP.AND.HDRSIZE.AND.((.NOT.STRATEGY).AND.(.NOT.TESTGEN)))" block
C
      END IF
C     
C---- Set an appropriate BLOCK size if not set explicitly
C     
      IF (IBLOCK.EQ.0) CALL GETBLOCK(NPRUN,PHIRNG,NBLOCK)
      IF (DEBUG(52)) THEN
        WRITE(IOUT,FMT=7770) NBLOCK, PHIRNG
        IF (ONLINE) WRITE(ItOUT,FMT=7770) NBLOCK, PHIRNG
 7770   FORMAT(/,1x,'After calling GETBLOCK, NBLOCK=',I4,' and',
     +         ' PHIRNG is',F5.2)
      END IF
      IF(SOCKLO)THEN
 7000    FORMAT('<?xml version="1.0"?>',
     $        '<!DOCTYPE block_size_notification SYSTEM "">',
     $        '<block_size_notification>',
     $        '<status>',
     $        '<code>ok</code>',
     $        '</status>',
     $        '<block_size>',I3,'</block_size>',
     $        '</block_size_notification>')
         WRITE(LONGOUTLINE,FMT=7000)NBLOCK
         CALL WRITE_SOCKET_LENGTH(SERVERFD,
     $        LENSTR(LONGOUTLINE),LONGOUTLINE)
      ENDIF
C     
C---- Final tests on distance, wavelength (may be read from header)
C     
 780  IF ((ISYN.NE.0).AND.(IWAVE.EQ.0)) THEN
        IF ((USETAIL).AND.
     +       ((.NOT.STRATEGY).AND.(.NOT.TESTGEN))) THEN
          WRITE(IOUT,FMT=6723)
          IF (ONLINE) WRITE(ITOUT,FMT=6723)
 6723     FORMAT(/,1X,'Wavelength will be taken from image tailer')
        ELSE
          WRITE (IOUT,FMT=7036)
          IF (ONLINE) WRITE (ITOUT,FMT=7036)
          STOPRUN = .TRUE.
        END IF
      END IF
C     
C     
C---- Check for DISTANCE keyword
C     
      IF (IDIST.EQ.0) THEN
        IF ((USETAIL).AND.
     +       ((.NOT.STRATEGY).AND.(.NOT.TESTGEN))) THEN
          WRITE(IOUT,FMT=6721)
          IF (ONLINE) WRITE(ITOUT,FMT=6721)
 6721     FORMAT(/,1X,'Distance will be taken from image tailer')
        ELSE
          WRITE (IOUT,FMT=6720)
          IF (ONLINE) WRITE (ITOUT,FMT=6720)
 6720     FORMAT(1X,'****** No DISTANCE keyword has been given *****')
          STOPRUN = .TRUE.
        END IF
      END IF
C     
C     ***** machine specific code follows *****
C     
C
C---- Print for program flow
C
        CALLED = ' '
        COMMENT = '*** Now using machine (detector) specific code ***'
        IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
C     
C---- Quantum 4 binned
C     
      IF (MACHINE .EQ. 'ADSC')THEN
        IF (abs(RAST-0.1632).le.1e-5)THEN
          nrec = 1152
          iylen = 1152
          TILEX(1) = 577
          TILEY(1) = 577
          TILEWX(1) = 2
          TILEWY(1) = 2
C     
C---- Quantum 210 unbinned
C     
        ELSEIF ((NREC.EQ.4096).AND.(IYLEN.EQ.4096))THEN
          MODEL = 'Q210'
          IF (XSCAN.EQ.9400) XSCAN = 10500
          IF (YSCAN.EQ.9400) YSCAN = 10500
          XMAXIP = 10500
          YMAXIP = 10500
          IF (ILIMIT.EQ.0) THEN
            XLIMIT = 0.01*(XMAXIP/2)
            LIMIT = NINT(100*XLIMIT)
          END IF
          TILEX(1) = 2049
          TILEY(1) = 2049
          TILEWX(1) = 8
          TILEWY(1) = 8
C     
C---- Quantum 210 binned
C     
        ELSEIF ((NREC.EQ.2048).AND.(IYLEN.EQ.2048))THEN
          MODEL = 'Q210'
          IF (XSCAN.EQ.9400) XSCAN = 10500
          IF (YSCAN.EQ.9400) YSCAN = 10500
          XMAXIP = 10500
          YMAXIP = 10500
          IF (ILIMIT.EQ.0) THEN
            XLIMIT = 0.01*(XMAXIP/2)
            LIMIT = NINT(100*XLIMIT)
          END IF
          TILEWX(1) = 4
          TILEWY(1) = 4
          TILEX(1) = 1025
          TILEY(1) = 1025
C     
C---- Quantum 315 unbinned
C     
        ELSEIF ((NREC.EQ.6144).AND.(IYLEN.EQ.6144))THEN
C
C---- image too big for arrays...
C
          IF((IYLENGTH.LT.6144).OR.(IXWDTH.LT.12288))THEN
            WRITE(IOUT,FMT=6650)
            IF(ONLINE)WRITE(ITOUT,FMT=6650)
            NSHUTERR = 52
            CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
          ENDIF
          MODEL = 'Q315'
          IF (XSCAN.EQ.9400) XSCAN = 15750
          IF (YSCAN.EQ.9400) YSCAN = 15750
          XMAXIP = 15750
          YMAXIP = 15750
          IF (ILIMIT.EQ.0) THEN
            XLIMIT = 0.01*(XMAXIP/2)
            LIMIT = NINT(100*XLIMIT)
          END IF
          TILEX(1) = 2049
          TILEY(1) = 2049
          TILEWX(1) = 8
          TILEWY(1) = 8
C     
C---- Quantum 315 binned
C     
        ELSEIF ((NREC.EQ.3072).AND.(IYLEN.EQ.3072))THEN
          MODEL = 'Q315'
          IF (XSCAN.EQ.9400) XSCAN = 15750
          IF (YSCAN.EQ.9400) YSCAN = 15750
          XMAXIP = 15750
          YMAXIP = 15750
          IF (ILIMIT.EQ.0) THEN
            XLIMIT = 0.01*(XMAXIP/2)
            LIMIT = NINT(100*XLIMIT)
          END IF
          TILEX(1) = 1025
          TILEY(1) = 1025
          TILEWX(1) = 4
          TILEWY(1) = 4
C     
C---- Quantum 270 unbinned
C     
        ELSEIF ((NREC.EQ.4168).AND.(IYLEN.EQ.4168))THEN
          MODEL = 'Q270'
          XSCAN = 13500
          YSCAN = 13500
          XMAXIP = 13500
          YMAXIP = 13500
          RMAXIP = SQRT(XMAXIP**2 + YMAXIP**2)
          IF (ILIMIT.EQ.0) THEN
            XLIMIT = 0.01*(XMAXIP/2)
            LIMIT = NINT(100*XLIMIT)
          END IF
          NTILEX = 2
          NTILEY = 2
          TILEX(1) = 2080
          TILEY(1) = 2080
          TILEWX(1) = 8
          TILEWY(1) = 8
C     
C---- Quantum 270 binned
C     
        ELSEIF ((NREC.EQ.2084).AND.(IYLEN.EQ.2084))THEN
          MODEL = 'Q270'
          XSCAN = 13500
          YSCAN = 13500
          XMAXIP = 13500
          YMAXIP = 13500
          RMAXIP = SQRT(XMAXIP**2 + YMAXIP**2)
          IF (ILIMIT.EQ.0) THEN
            XLIMIT = 0.01*(XMAXIP/2)
            LIMIT = NINT(100*XLIMIT)
          END IF
          NTILEX = 2
          NTILEY = 2
          TILEX(1) = 1040
          TILEY(1) = 1040
          TILEWX(1) = 4
          TILEWY(1) = 4

        END IF
      END IF
      
C     
C---- Calculate orientation matrix from cell param and U matrix
C     if orientation matrix not given.  Otherwise calculate cell
C     param from the orientation matrix.
C     NB Cannot do this earlier because need wavelength to convert AMAT
C     to
C     real cell parameters, and wavelength may only be read from
C     image header (or tailer)
      ICHECK = 1
      IF (ICELL.EQ.1) ICHECK = 0
C     
C     ************************
      CALL SETMAT(IMAT,IUMAT,ICELL,ICHECK)
C     ************************
C
C---- Need to reset CRYSTALUBREF because AMAT may have changed if a CELL
C     keyword was used to impose new cell parameters. Only important for
C     multilattice cases. ?? Does this need to be changed if CELL is
C     updated via iMosflm ?? Perhaps don't need test on UBREFSET ? 
C     Must not do this is MERGELATT is true, because ILATTICE will not
C     be correct if LATTICE keyword specifies anything other than
C     the first lattice. 
C
      IF ((NLATTICE.GT.1).AND.(.NOT.MERGELATT)) THEN
        ISOL = BESTSOLAR(ILATTICE)
        DO J =1,3
          DO K =1,3
            CRYSTALUBREF(J,K,ISOL,ILATTICE) = AMAT(J,K)
          ENDDO
        ENDDO
        IF (DEBUG(52)) THEN
          WRITE(IOUT,FMT=7830) ISOL,ILATTICE
          IF (ONLINE) WRITE(ITOUT,FMT=7830) ISOL,ILATTICE
 7830     FORMAT(1X,'In CONTROL, CRYSTALUBREF for solution',I3,
     +            ' lattice',I3,' set to AMAT')
        END IF
      END IF
C     
C---- If a CELL keyword had been given, then these cell parameters
C     overwrite
C     those in any MATRIX file. However, after the first "run" keyword
C     has been given, we want to be able to read in new MATRIX files and
C     use their cell parameters, so UNLESS "KEEP" has been given on the
C     cell
C     keyword, reset ICELL to 0.
C     However, do NOT want to do this if doing prediction from mxdspl
C     and
C     a cell and UMAT have been given.
C     
      IF (.NOT.CELLKEEP) ICELL = 0
      IMISS = 0
      IMISSMAT = 0
C     
C---- Check cell parameters are consistent with symmetry
C     
      IFLAG = 0
      CALL CELLCHK(ICRYST,CELL,IFLAG)
      IF ((IFLAG.NE.0).AND.(NUMSPG.GT.0)) THEN
        WRITE(IOUT,FMT=7049) NUMSPG,SPGNAM(1:LENSTR(SPGNAM)),CELL,
     $       AMAT
        IF (ONLINE) WRITE(ITOUT,FMT=7049) NUMSPG,
     $       SPGNAM(1:LENSTR(SPGNAM)),CELL,AMAT
        IF (BRIEF) WRITE(IBRIEF,FMT=7049) NUMSPG,
     $       SPGNAM(1:LENSTR(SPGNAM)),CELL,AMAT
 7049   FORMAT(//,1X,'*** ERROR ***',/,1X,'Space group number ',
     +       I4,' (',A,') is inconsistent with cell',
     +       ' parameters derived from orientation matrix.',/,1X,
     +       '(Cell:',6F10.3,')',/,'AMAT: ',3F9.6,2(/,6X,3F9.6),
     +       /,1X,'Provide consistent cell',
     +       ' parameters with a CELL keyword')
        STOPRUN = .TRUE.
          IF(SOCKLO)THEN
            LINE400 = ' '
            WRITE(LINE400,FMT=8020)SPGNAM(1:LENSTR(SPGNAM)),CELL
            CALL WRITE_SOCKET_LENGTH(SERVERFD,
     $           LENSTR(LINE400),LINE400)
          ENDIF
      END IF
C     
C---- Now strictly impose symmetry constraints on cell
C     
      SOFTCELL = 1.0
      CALL CELLFIX(CELL)
      CALL CELLFIX(RCELL)

C     
C---- By this point, must know the wavelength. If not CuKalpha or Mo,
C     assume a synchrotron source and set appropriate defaults.
C     
      IF ((ABS(WAVE-1.5418).GT.0.0019).AND.(WAVE.GT.0.05)) THEN
        IF (ABS(WAVE-0.7107).GT.0.0002) THEN
          IF (FRSTWARN) THEN
            WRITE(IOUT,FMT=7440) WAVE
            IF (ONLINE) WRITE(ITOUT,FMT=7440) WAVE
            FRSTWARN = .FALSE.
          END IF
 7440     FORMAT(/,1X,'**** WARNING ****',/,1X,'Because the wavelen',
     +         'gth (',F7.4,') is not CuKa (1.5418) or Mo (0.7107),',
     +         /,1X,'source is assumed to',
     +         ' be a synchrotron and synchrotron defaults for ',
     +         'polarisation and beam divergence',/,1X,'will be use',
     +         'd if these have not been defined explicitly ',
     +         '(SYNCH POLAR and DIVH/DIVV keywords.',/)
          ISYN = 1
          IMONO = 2
          IF (IPOLAR.EQ.0) TOR = TORSRS
C---- If TOLERANCE not set, make default 0.03
C     
          IF (ITOL.EQ.0) THEN
            TOLMIN = 0.02
            TOL = 0.03
          END IF
          IF (IDELAMB.EQ.0) THEN
            IDELAMB = 1
            DELAMB = 0.0002
          END IF
          IF (IDIVH.EQ.0) THEN
            IDIVH = 1
            DIVHD = 0.02
          END IF
          IF (IDIVV.EQ.0) THEN
            IDIVV = 1
            DIVVD = 0.02
          END IF
          DIVH = 0.5*DTOR*DIVHD
          DIVV = 0.5*DTOR*DIVVD
        END IF
      END IF
C     
C---- By this point, must know the wavelength. If less than 0.05
C     assume electron diffraction data and set beam parameters. Note
C     that detector parameters have been changed in SECTION.
C     
      IF (WAVE.LT.0.05) THEN
          IF (FRSTWARN) THEN
            WRITE(IOUT,FMT=7442) WAVE
            IF (ONLINE) WRITE(ITOUT,FMT=7442) WAVE
            FRSTWARN = .FALSE.
          END IF
 7442     FORMAT(/,1X,'**** WARNING ****',/,1X,'Because the wavelength',
     +        '(',F7.4,') is less than 0.05A the images are assumed to',
     +         /,1X,'be electron diffraction',
     +         ' and polarisation, dispersion and beam divergence ',
     +         '',/,1X,'will be set to zero ',
     +         'unless these have been defined explicitly by keywords',
     +         /,1X,'POLAR, DISP and DIVERGENCE.',/)
          ISYN = 0
          IMONO = 0
C
C---- If not set by a keyword, set polarisation, dispersion and divergence all to zero
C
          IF (IPOLAR.EQ.0) TOR = 0.0
          IF (IDELAMB.EQ.0) THEN
            IDELAMB = 1
            DELAMB = 0.0000
          END IF
          IF (IDIVH.EQ.0) THEN
            IDIVH = 1
            DIVHD = 0.00
          END IF
          IF (IDIVV.EQ.0) THEN
            IDIVV = 1
            DIVVD = 0.00
          END IF
          DIVH = 0.5*DTOR*DIVHD
          DIVV = 0.5*DTOR*DIVVD
      END IF

C     
 784  CONTINUE
      IF (IBEAM.EQ.0) THEN
         WRITE (IOUT,FMT=7047)
         IF (ONLINE) THEN
            WRITE (ITOUT,FMT=7047)
         END IF
 7047    FORMAT (/,1X,'****** ERROR *****',/,1X,'No BEAM keyword',
     +        ' (specifying coordinates ',
     +        'of direct beam position) has been given,'/,' nor ',
     $        'has a beam centre been found in the image header')
         STOPRUN = .TRUE.
      END IF
      IF (STOPRUN) THEN
C     
C---- Trap being called from MXDSPL
C     
        IF ((MODE.GT.0).AND.(MODE.LT.10)) THEN
          MODE = 99
          RETURN
        END IF
C     
        IF (ONLINE) THEN 
          COMREAD = .FALSE.
          STOPRUN = .FALSE.
          GOTO 50
        ELSE
          STOP
        END IF
      END IF
C     
C---- Set geometric limits on detector
C     
      
C     
      IF (ICASS.EQ.0) THEN
C     
C---- FLAT film detector
C     
        IF (XMAX.EQ.0.0) XMAX = FXMAX
        IF (YMAX.EQ.0.0) YMAX = FYMAX
        IF (RMAX.EQ.0.0) RMAX = MIN(SQRT(XMAX*XMAX+YMAX*YMAX),FXMAX)
        IF (RMIN.EQ.0.0) RMIN = 350.0
      ELSE IF (ICASS.EQ.1) THEN
C     
C---- VEE cassette
C     
        VEE = .TRUE.
        IF (XMIN.EQ.0.0) XMIN = VXMIN
        IF (XMAX.EQ.0.0) XMAX = VXMAX
        IF (YMAX.EQ.0.0) YMAX = VYMAX
        IF (RMIN.EQ.0.0) RMIN = 350.0
        IF (RMAX.EQ.0.0) RMAX = VRMAX
      ELSE IF (ICASS.EQ.4) THEN


C     
C---- IP detector (Mar or Rigaku)
C     
C     
C     ***** machine specific code follows *****
C---- For scanners with more than one possible image size, cannot assign
C     these limits before an image has been read (to determine the size)
C     .
C     
C---- Set defaults for Mar scanners (extended to deal with Mar345 images
C     written with "IMAGE" or "pck" formats, also Mar555 flat panel.
C     
C     Possibilities are now:
C     NREC     SIZE     PIXEL
C     3450     345.0     0.1
C     3000     300.0     0.1
C     2400     240.0     0.1
C     2300     345.0     0.15
C     2000     300.0     0.15
C     1800     180.0     0.1
C     1600     240.0     0.15
C     1200     180.0     0.15
C     3072   flat panel  0.139 (IYLEN = 2560)
C
        IF (MACHINE.EQ.'MAR ') THEN
          IF (NREC.EQ.3450) THEN
            XMAXIP = 17250
            YMAXIP = 17250
            RMAXIP = 17250
            RSCANIP = 17250
            IF (IPIX.EQ.0) RAST = 0.10
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
          ELSE IF (NREC.EQ.3000) THEN
            XMAXIP = 15000
            YMAXIP = 15000
            RMAXIP = 15000
            RSCANIP = 15000
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
            IF (IPIX.EQ.0) RAST = 0.10
          ELSE IF (NREC.EQ.2400) THEN
            XMAXIP = 12000
            YMAXIP = 12000
            RMAXIP = 12000
            RSCANIP = 12000
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
            IF (IPIX.EQ.0) RAST = 0.10
          ELSE IF (NREC.EQ.2300) THEN
            XMAXIP = 17250
            YMAXIP = 17250
            RMAXIP = 17250
            RSCANIP = 17250
            IF (IPIX.EQ.0) RAST = 0.15
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
          ELSE IF (NREC.EQ.2000) THEN
            XMAXIP = 15000
            YMAXIP = 15000
            RMAXIP = 15000
            RSCANIP = 15000
            IF (IPIX.EQ.0) RAST = 0.15
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
          ELSE IF (NREC.EQ.1800) THEN
            XMAXIP = 9000
            YMAXIP = 9000
            RMAXIP = 9000
            RSCANIP = 9000
            IF (IPIX.EQ.0) RAST = 0.10
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
          ELSE IF (NREC.EQ.1600) THEN
            XMAXIP = 12000
            YMAXIP = 12000
            RMAXIP = 12000
            RSCANIP = 12000
            IF (IPIX.EQ.0) RAST = 0.15
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
          ELSE IF (NREC.EQ.1200) THEN
            XMAXIP = 9000
            YMAXIP = 9000
            RMAXIP = 9000
            RSCANIP = 9000
            IF (IPIX.EQ.0) RAST = 0.15
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
          ELSE
C
C---- If this is a pck format file, it will have been recognised as Mar IP, but this might
C     not be the case, eg Bruker files written in pck format. Try to allow for that here.
C
            XMAXIP =  NINT(NREC*50*RAST)
            YMAXIP =  NINT(IYLEN*50*RAST)
            RSCANIP = XMAXIP
            RMAXIP = SQRT(XMAXIP**2 + YMAXIP**2)
            XSCAN = XMAXIP
            YSCAN = YMAXIP

            IF (RSCANIP.EQ.0) THEN
              WRITE(IOUT,FMT=7494) NREC
              IF (ONLINE) WRITE(ITOUT,FMT=7494) NREC
 7494         FORMAT(//,1X,'*** FATAL ERROR ***',/,1X,
     +             'Image size of',I6,' not recognised for Mar ',
     +             'Scanners.')
              IF (MODE.EQ.10) THEN
                LINE = 'ERROR...see terminal window'
                CALL MXDWIO(LINE,1)
              END IF
              NSHUTERR = 53
cal              CALL SHUTDOWN(CALLEDFROM,LENSTR(CALLEDFROM))
            END IF
          END IF

        ELSE IF (MACHINE.EQ.'MARF') THEN
C
C---- Flat panel, 3072 x 2560
C
          IF (NREC.EQ.3072) THEN
            XMAXIP = 21350
            YMAXIP = 17792
C
C---- what should RMAXIP be for a rectangular detector? diagonal to 
C     corner
C
            RMAXIP = 27792
            RSCANIP = RMAXIP
            IF (IPIX.EQ.0) RAST = 0.139
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
            CIRCULAR = .FALSE.
            ORTHOG = .TRUE.
            omegafd = 0.0
            omegaf = 0.0
          ENDIF
        ELSE IF (MACHINE.EQ.'MARC') THEN
C     
C     Set up detector limits based on pixel size,
C     nominally 80 um (for 165mm CCD) vs 64 microns for 133mm
C     
chrp 17.06.2014          IF (RAST.GT.0.070) THEN
chrp 17.06.2014            RMINIP = 300
chrp 17.06.2014            XMAXIP = 8100
chrp 17.06.2014            YMAXIP = 8100
chrp 17.06.2014            RMAXIP = 8100
chrp 17.06.2014            RSCANIP = RMAXIP
chrp 17.06.2014            IF (ILIMIT.EQ.0) THEN
chrp 17.06.2014              XLIMIT = 0.01*(XMAXIP/2)
chrp 17.06.2014              LIMIT = NINT(100*XLIMIT)
chrp 17.06.2014            END IF
chrp 17.06.2014          END IF
C
C---- Mods for new tiled 3x3 ccd ("Mar Mosaic 225?") or new fast MX225HS
C
          IF (NREC.EQ.3072.OR.NREC.EQ.2880) THEN
chrp 17.06.2014            RMINIP = 300
chrp 17.06.2014            XMAXIP = 11200
chrp 17.06.2014            YMAXIP = 11200
chrp 17.06.2014            RMAXIP = 15800
chrp 17.06.2014            IF (XSCAN.EQ.0) XSCAN = 11200
chrp 17.06.2014            IF (YSCAN.EQ.0) YSCAN = 11200
chrp 17.06.2014            RSCANIP = RMAXIP
            CIRCULAR = .FALSE.
            ORTHOG = .TRUE.
chrp 17.06.2014            IF (ILIMIT.EQ.0) THEN
chrp 17.06.2014              XLIMIT = 0.01*(XMAXIP/2)
chrp 17.06.2014              LIMIT = NINT(100*XLIMIT)
chrp 17.06.2014            END IF
C
C---- Mods for new tiled 4x4 ccd ("Mar Mosaic 325")
C
          ELSEIF (NREC.EQ.4096) THEN
chrp 17.06.2014            RMINIP = 300
chrp 17.06.2014            XMAXIP = 16250
chrp 17.06.2014            YMAXIP = 16250
chrp 17.06.2014            RMAXIP = 22980
chrp 17.06.2014            IF (XSCAN.EQ.0) XSCAN = 16250
chrp 17.06.2014            IF (YSCAN.EQ.0) YSCAN = 16250
chrp 17.06.2014            RSCANIP = RMAXIP
            CIRCULAR = .FALSE.
            ORTHOG = .TRUE.
chrp 17.06.2014            IF (ILIMIT.EQ.0) THEN
chrp 17.06.2014              XLIMIT = 0.01*(XMAXIP/2)
chrp 17.06.2014              LIMIT = NINT(100*XLIMIT)
chrp 17.06.2014            END IF
          END IF
C
C---- HRP 16th June 2014 calculate XMAXIP etc from pixel size and NREC, IYLEN;
C     over-ride hard-coded values above; use RAST here, HRAST in gethdr.f
C
          XMAXIP = INT(NREC*RAST)*50
          YMAXIP = INT(IYLEN*RAST)*50
          IF(CIRCULAR)THEN
            RMAXIP = XMAXIP
          ELSE
            RMAXIP = SQRT(XMAXIP**2+YMAXIP**2)
          ENDIF
          XSCAN = XMAXIP
          YSCAN = YMAXIP
          RSCANIP = RMAXIP
          IF (ILIMIT.EQ.0) THEN
            XLIMIT = 0.01*(XMAXIP/2)
            LIMIT = NINT(100*XLIMIT)
          END IF
        ELSE IF (MACHINE.EQ.'LMB') THEN
          XMAXIP = 25000
          YMAXIP = 25000
          RMAXIP = 25000*SQRT(2.0)
          RSCANIP = RMAXIP
          IF (ILIMIT.EQ.0) THEN
            XLIMIT = 0.01*(XMAXIP/2)
            LIMIT = NINT(100*XLIMIT)
          END IF
        ELSE IF (MACHINE.EQ.'CCD2') THEN
          XMAXIP = 7000
          YMAXIP = 9000
          RMAXIP = 8770
          RSCANIP = RMAXIP
          IF (ILIMIT.EQ.0) THEN
            XLIMIT = 0.01*(XMAXIP/2)
            LIMIT = NINT(100*XLIMIT)
          END IF
        ELSE IF (MACHINE.EQ.'RAXI') THEN
          IF ((MODEL.EQ.'RAXISIV').OR.(MODEL.EQ.'RAXISHTC')) THEN
            XMAXIP = 15000
            YMAXIP = 15000
            RMAXIP = 15000*SQRT(2.0)
            RSCANIP = RMAXIP
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
          ELSE IF (MODEL.EQ.'RAXISV') THEN
            XMAXIP = 20000
            YMAXIP = 20000
            RMAXIP = 20000*SQRT(2.0)
            RSCANIP = RMAXIP
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
          ELSE
            XMAXIP = 10000
            YMAXIP = 10000
            RMAXIP = 10000*SQRT(2.0)
            RSCANIP = RMAXIP
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
          END IF
        ELSE IF (MACHINE.EQ.'DIP2') THEN
          XMAXIP = 10000
          YMAXIP = 10000
          RMAXIP = 10000
          RSCANIP = RMAXIP
          IF (ILIMIT.EQ.0) THEN
            XLIMIT = 0.01*(XMAXIP/2)
            LIMIT = NINT(100*XLIMIT)
          END IF
          IF (NREC.EQ.3000) THEN
            XMAXIP = 15000
            YMAXIP = 15000
            RMAXIP = 15000
            RSCANIP = 15000
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
          END IF
          IF (NREC.EQ.4000) THEN
            XMAXIP = 20000
            YMAXIP = 20000
            RMAXIP = 20000
            RSCANIP = 20000
            IF (ILIMIT.EQ.0) THEN
              XLIMIT = 0.01*(XMAXIP/2)
              LIMIT = NINT(100*XLIMIT)
            END IF
          END IF
        ELSE IF (MACHINE.EQ.'FUJI') THEN
          XMAXIP = 12500
          YMAXIP = 10000
          RMAXIP = SQRT(XMAXIP**2+YMAXIP**2)
          RSCANIP = RMAXIP
          IF (ILIMIT.EQ.0) THEN
            XLIMIT = 0.01*(XMAXIP/2)
            LIMIT = NINT(100*XLIMIT)
          END IF
C     
C---- set up parameters for CBF files based on values in image header
C     
        ELSE IF (MACHINE.EQ.'CBF ') THEN
          XMAXIP = NINT(NREC*50*RAST)
          YMAXIP = NINT(IYLEN*YSCAL*50*RAST)
          RMAXIP = SQRT(XMAXIP**2+YMAXIP**2)
          RSCANIP = RMAXIP
          IF(XLIMIT.EQ.0)THEN
            XLIMIT = 0.01*(XMAXIP/2)
            LIMIT = NINT(100*XLIMIT)
          ENDIF
C     
C---- add new machines below here
C     
C     
C---- add new machines above here
C     
        END IF
        XDMID = 0.5*NREC*RAST*100.0
        YDMID = 0.5*IYLEN*RAST*100.0
C     
C---- Need to allow for swung out detector and direct beam coords given
C     for twotheta=0 rather than true twotheta
C     
        IF ((TWOTHETA.NE.0.0).AND.(ISWUNG.EQ.0)) THEN 
          XTRUE = XCENMMIN(1) + 
     +         COS(OMEGAF)*0.01*XTOFD*TAN(TWOTHETA*DTOR)
          YTRUE = YCENMMIN(1) + 
     +         SIN(OMEGAF)*YSCAL*0.01*XTOFD*TAN(TWOTHETA*DTOR)
          IF (XOFF.EQ.0.0) XOFF = XDMID - 100.0*XTRUE
          IF (YOFF.EQ.0.0) YOFF = YDMID - 100.0*YTRUE
          IF (DEBUG(52)) THEN
            WRITE(IOUT,FMT=7124) XTRUE,YTRUE
            IF (ONLINE) WRITE(ITOUT,FMT=7124) XTRUE,YTRUE
          END IF
        ELSE
 7124     FORMAT(1X,'XCENMM, YCENMM  after correction for ',
     +         'swing angle ',2F8.2)
chrp05.12.2005          IF (XOFF.EQ.0.0) XOFF = XDMID - 100.0*XCENMM(1,1)
chrp05.12.2005          IF (YOFF.EQ.0.0) YOFF = YDMID - 100.0*YCENMM(1,1)
          IF (XOFF.EQ.0.0) XOFF = XDMID - 100.0*XCENMMIN(1)
          IF (YOFF.EQ.0.0) YOFF = YDMID - 100.0*YCENMMIN(1)
        END IF
C     
        IF (XMAX.EQ.0.0) THEN
C     
C---- Note for the following to work, NREC and IYLEN must have
C     been assigned. For Mar and RAXIS they are determined from
C     the header record of the image. FOR DIP2 they are assigned
C     when the SCANNER keyword is read. This is because XOFF is
C     determined
C     from XMID which is determined by NREC. In fact this should ALWAYS
C     be
C     the case, so no need to test.
C     
C     AL          IF ((MACHINE.EQ.'MAR ').OR.(MACHINE.EQ.'RAXI').OR.
C     AL     +        (MACHINE.EQ.'DIP2').OR.(MACHINE.EQ.'CCD1').OR.
C     AL     +        (MACHINE.EQ.'CCD2')) THEN

          XMAX = XMAXIP + ABS(XOFF)

C     AL          ELSE
C     AL            XMAX = XMAXIP
C     AL          END IF
        END IF
        IF (YMAX.EQ.0.0) THEN
C     AL          IF ((MACHINE.EQ.'MAR ').OR.(MACHINE.EQ.'RAXI').OR.
C     AL     +        (MACHINE.EQ.'DIP2').OR.(MACHINE.EQ.'CCD1').OR.
C     AL     +        (MACHINE.EQ.'CCD2')) THEN
          YMAX = YMAXIP + ABS(YOFF)
C     AL          ELSE
C     AL            YMAX = YMAXIP
C     AL          END IF
        END IF
        IF (RMIN.EQ.0.0) THEN
          IF (RMINIP.GT.0) THEN
            RMIN = RMINIP
          ELSE
            RMIN = MIN(0.05*XMAX, 0.05*YMAX)
          END IF
        END IF
        IF (RMAX.EQ.0.0) THEN
C     AL          IF ((MACHINE.EQ.'MAR ').OR.(MACHINE.EQ.'RAXI').OR.
C     AL     +        (MACHINE.EQ.'DIP2')) THEN
          IF (CIRCULAR) THEN
            RMAX = RMAXIP + SQRT(XOFF**2+YOFF**2)
          ELSE
            RMAX = SQRT(XMAX**2 + YMAX**2)
          END IF 
C     AL          ELSE
C     AL            RMAX = RMAXIP
C     AL          END IF
        END IF
        IF (RSCAN.EQ.0.0) RSCAN = RSCANIP
        IF (NFGEN.EQ.3) NFGEN = 1
C     
C---- XSCAN, YSCAN are used to test if spots lie off the physical edges
C     of the detector (this is for rectangular detectors). For circular
C     detectors, RSCAN is used. XSCAN,YSCAN,RSCAN are maximum coords
C     relative to the physical centre of the detector (ie the mid-point
C     of
C     the image), NOT wrt the direct beam position.
C     
        IF (MACHINE.EQ.'RAXI') THEN
          IF ((MODEL.EQ.'RAXISIV').OR.(MODEL.EQ.'RAXISHTC')) THEN
            IF (XSCAN.EQ.0.0) XSCAN = 14975
            IF (YSCAN.EQ.0.0) YSCAN = 14975
          ELSE IF (MODEL.EQ.'RAXISV') THEN
            IF (XSCAN.EQ.0.0) XSCAN = 19975
            IF (YSCAN.EQ.0.0) YSCAN = 19975
          ELSE
            IF (XSCAN.EQ.0.0) XSCAN = 9975
            IF (YSCAN.EQ.0.0) YSCAN = 9975
          END IF
        ELSE IF (MACHINE.EQ.'FUJI') THEN
          IF (XSCAN.EQ.0.0) XSCAN = 12475
          IF (YSCAN.EQ.0.0) YSCAN = 9975
        ELSE IF (MACHINE.EQ.'CCD2') THEN
          IF (XSCAN.EQ.0.0) XSCAN = 7000
          IF (YSCAN.EQ.0.0) YSCAN = 9000
        END IF
      END IF
C     
C     
      RMNSQD = RMIN*RMIN
      RMXSQD = RMAX*RMAX
C     
C---- If BIAS keyword has been given, to add a constant to all pixel
C     values because there are zero pixel values within scanned area, then must
C     ensure that NO spots are predicted whose measurement box might
C     contain any pixels outside the scanned area, because can no longer use
C     test of a zero pixel value to check for pixels outside scanned area.
C     Assume measurement box has maximum size of 31 pixels ie 4.5mm across
C     
C---- but don't do this if we are running 'predict' from the GUI
C
      REDGE = MAX(2.25,REAL(NXS)*RAST*0.5)
      IF ((ICONST.NE.0).AND.(MODE.NE.1)) THEN
        RSCAN = RSCAN - 100*REDGE
        XSCAN = XSCAN - 100*REDGE
        YSCAN = YSCAN - 100*REDGE
      END IF
      RSCANSQ = RSCAN*RSCAN
C     
C---- Determine the radius of the reciprocal sphere DSTMAX.
C     DSTMAX is set to the minimum value of:
C     a) that obtained from the maximum radius of the film (DSTRMX)
C     b) that obtained from the resolution in Angstroms (if given)
C     (DSTRES)
C     c) DSTMAX in the input data (if given)
C     
C     
      IF (ICASS.EQ.0) THEN
C     
C---- Film case
C     
        THETA = ATAN(RMAX/XTOFD)*0.5
C     
      ELSE IF (ICASS.EQ.4) THEN
C     
C---  IP...Allow for translationally offset detector. 
C     NOTE dimensions are in 10 micron units here.
C     
        ROFFMAX = MAX(XMAX,YMAX)
C     AL        ROFFMAX = MIN(RMAX,ROFFMAX)
        ROFFMAX = MAX(RMAX,ROFFMAX)
        THETA = ATAN(ROFFMAX/XTOFD)*0.5
        IF (ABS(TWOTHETA).GT.0.0) THEN
          IF (CIRCULAR) THEN
            THETA = 0.5*(ABS(TWOTHETA)*DTOR + ATAN(RSCAN/XTOFD))
          ELSE IF (ORTHOG) THEN
            THETA = 0.5*(ABS(TWOTHETA)*DTOR + 
     +           ATAN(SQRT(XSCAN**2+YSCAN**2)/XTOFD))
          END IF
        END IF
      ELSE IF (ICASS.EQ.1) THEN
C     
C---- Vee cassette
C     
        THETA = ATAN(SQRT(3.0)/ (2.0*XTOFD/RMAX-1.0))*0.5
      END IF
      DSTRMX = SIN(THETA)*2.0
C
C---- if RESOLUTION has been given in a keyword, then set DSTMAX now.
C     
      IF (INRES.EQ.1) THEN
        IF (RES.NE.0) DSTMAX = WAVE/RES
C     
      ELSE
        IF (DSTMAX.EQ.0.0) THEN
          DSTMAX = DSTRMX
        ELSE
          DSTMAX = MIN(DSTMAX,DSTRMX)
        END IF
      END IF
C     
chrp 14.03.2003      IF (RES.NE.0) DSTRES = WAVE/RES
chrp 14.03.2003      DSTMAX = MIN(DSTMAX,DSTRES)
      DSTMAXS = DSTMAX
C     
C---- Need to allow for RMIN being defined relative to a centre other
C     than
C     the point of intersection of direct beam and detector
C     (set up by LIMITS RCENTRE.
C     
      IF ((RMINX.NE.0.0).OR.(RMINY.NE.0.0)) THEN
        RMINP = RMIN - 100.0*SQRT((XCENMM(1,1)-0.01*RMINX)**2 + 
     +       (YCENMM(1,1)-0.01*RMINY)**2)
C     
C---- RMINP negative actually means backstop shadow does not cover
C     direct beam position ! Set RMINP to a very small value so that
C     DMAX is not infinite. DMAX IS HARDWIRED TO 1000A further on in
C     this
C     S/R
C     
        IF (RMINP.LT.0.0) RMINP = 0.0
      ELSE
        RMINP = RMIN
      END IF
      THETA = ATAN(RMINP/XTOFD)*0.5
C     
C---- DSTMIN is dimensionless rlu (=WAVE/DMAX where DMAX is 
C     max real D spacing)
C     
      DSTMIN = SIN(THETA)*2.0
C     
C     
      IF (RESLOW.NE.0.0) THEN
        DSTRES = WAVE/RESLOW
        DSTMIN = MAX(DSTRES,DSTMIN)
      END IF
C     
C---- Calculate the reciprocal sphere radius DSTPL - corresponding to
C     a slightly higher resolution to be used in checking overlaps on
C     the outside of the picture
C     
      DELR = MAX(IXSEP,IYSEP)*2.0
C     
C---- If separation not given, assume spot size of 2mm
C     
      IF (DELR.EQ.0) DELR = 400.0
C     
      THETA = ASIN(DSTMAX/2.0)
      T = TAN(2.0*THETA)
C     
      IF ((ICASS.EQ.0).OR.(ICASS.EQ.4)) THEN
C     
C---- Flat cassette or IP detector
C     
        RPLUS = XTOFD*T + DELR
        THPLUS = ATAN(RPLUS/XTOFD)*0.5
C     
      ELSE IF (ICASS.EQ.1) THEN
C     
C---- Vee cassette
C     
        RPLUS = 2.0*XTOFD*T/ (SQRT(3.0)+T) + DELR
        THPLUS = ATAN(SQRT(3.0)/ (2.0*XTOFD/RPLUS-1.0))*0.5
C     
      END IF
C     
      DSTPL = SIN(THPLUS)*2.0
      DSTPL2 = DSTPL*DSTPL
C     
C---- Open the generate file and write header
C     
C---- if a new generate file, close old generate file (if any),
C     
C---- If doing a stratgey run, and this is a rectangular detector, and
C     the resolution limit is beyond the inscribed circle, give a
C     warning
C     and allow the resolution to be reset.
C     
      IF (STRATEGY.AND.(TWOTHETA.EQ.0.0)) THEN
        IF ((ORTHOG).AND.(MODE.EQ.10)) THEN
          XYMAXST = MAX(XMAX,YMAX)
          IF (XTOFD.GT.0) THETAST = ATAN(XYMAXST/XTOFD)*0.5
          DSTMAXST = SIN(THETAST)*2
          IF (DSTMAXST.LT.DSTMAX) THEN
            NULINE = .TRUE.
            WRITE(IOUT,FMT=7660) WAVE/DSTMAX,WAVE/DSTMAXST
            IF (ONLINE)WRITE(ITOUT,FMT=7660) WAVE/DSTMAX,
     $           WAVE/DSTMAXST
            WRITE(IOLINE,FMT=7660) WAVE/DSTMAX,WAVE/DSTMAXST
 7660       FORMAT('**** WARNING ****',/,
     +           'The current resolution limit (',F5.2,'A) is ',
     +           'beyond the inscribed circle',/,1X,'limit (',
     $           F5.2,'A) ie extends into the corners of the',
     +           'detector.',/,1X,'You may not be able to ge',
     +           't 100% completeness at high resolution,'/,1X,
     $           'and so the overall completeness may be low.')
            CALL WINDIO(NULINE)
          END IF
        END IF
      END IF
C     
C     
      IF (STRATEGY.OR.TESTGEN) THEN
C     
C---- If this is not the first call to CONTROL, only reset ETA etc
C     if they have been input in this call. For 2nd and later calls
C     to CONTROL IMOSAIC etc are initialied at 2
C*****No longer needed, done when input
C     AL        IF (FIRSTTIME.OR.(IMOSAIC.EQ.1)) ETA = 0.5*DTOR*ETA
C     AL        IF (FIRSTTIME.OR.(IDIVH.EQ.1)) DIVH = 0.5*DTOR*DIVH
C     AL        IF (FIRSTTIME.OR.(IDIVV.EQ.1)) DIVV = 0.5*DTOR*DIVV
        GOTO 706
      END IF
C     
C---- If there is more than one set of SERIAL/RUN keywords then we
C     want to open a new generate file for each run, but keep the same
C     MTZ file unless a new HKLOUT keyword has been given.
C     However, do not do this test if NOGENOPEN has been set TRUE,
C     which coresponds to calling for a prediction for new GUI while
C     integrating a run, as this will cause an existing genfile to be
C     closed in cases where there are 2 PROCESS keywords (NSERTOT>1)
C      AGWL 23/3/06
C 
      IF ((NEWGENF.OR.((NSERTOT.GT.1).AND.(.NOT.MULTISEG))).AND.
     +     (.NOT.NOGENOPEN)) THEN
C     
C     ******************
        IF ((.NOT.FIRSTTIME) .AND. GENOPEN) THEN
C
C---- Print for program flow
C
          CALLED = 'QCLOSE'
          COMMENT = 'Close generate file when more than 1 RUN keyword'
          IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
          CALL QCLOSE(IUNIT)
          GENOPEN = .FALSE.
        ENDIF
C     
C---- If repeating an entire multiseg run, close the genfile
C     
        IF (RPTFIRST) THEN
c     Another instance of forgetting to set the genopen flag
           IF (GENOPEN) then
C
C---- Print for program flow
C
            CALLED = 'QCLOSE'
            COMMENT = 'Close generate file, repeating a MULTISEG run'
            IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
C
              CALL QCLOSE(IUNIT)
              GENOPEN = .FALSE.
           END IF
          IF (MTZOPEN) THEN
            MTZPRT = 1
C     *********************
            if (debug(52))then
              write(iout,fmt=7748)
              if (online) write(itout,fmt=7748)
 7748         format(1x,'**closing MTZ file in control 5 ***')
            end if
            CALLED = 'LWCLOS'
            COMMENT = 'To close MTZ file '
            IF (DEBUG(99))             CALLED = 'LWCLOS'
            COMMENT = 'To close MTZ file '
            IF (DEBUG(99)) 
     +        CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
            CALL LWCLOS(MTZOUT,MTZPRT)
            MTZOPEN = .FALSE.
            IF(LBEST.AND.BESTOPEN)THEN
              CLOSE(BESTHKL)
              BESTOPEN = .FALSE.
            ENDIF
C     *********************
          END IF
        END IF

C     
        IF (PRECESS) THEN
C     
C     ***********************************
          CALL PSTART(GENFILE,GTITLE,GENOPEN,IDENT)
C     ***********************************
C     
        ELSE
C     
          NTIMES = NRUN
          IF (RPTFIRST) NTIMES = 2
          CALLED = 'START'
          WRITE(COMMENT, FMT=7880) GENOPEN
 7880     FORMAT(' To open generate file (if closed) GENOPEN is',L2)
          IF (DEBUG(99)) 
     +        CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
C     *****************************************
          CALL START(GENFILE,GTITLE,GENOPEN,IDENT,INOGEN,NTIMES,
     +               NLATTICE,NPACK)
C     *****************************************
C     
        END IF
C     
C---- If this is first run, always open MTZ file (unless multiseg post
C     refinement)but on subsequent runs, only open new MTZ file if 
C     HKLOUT keyword given.
C     
        IF (NRUN.EQ.1) THEN
C
chrp 17042002
C
        IF (.not. PNAMEgiven ) THEN
          PROJECTNAME = 'Unspecified'
          PNAMEGIVEN = .TRUE.
          WRITE (IOUT,FMT=8001)PROJECTNAME(1:LENSTR(PROJECTNAME))
          IF(ONLINE)WRITE (ITOUT,FMT=8001)
     $         PROJECTNAME(1:LENSTR(PROJECTNAME))
 8001     FORMAT(/,
     $         ' Warning: No PROTEIN NAME given by KeyWord PNAME',/,
     $         ' It has been set to "',A,'"',/)
        END IF
        IF (.not. XNAMEgiven ) THEN
          CRYSTALNAME = IDENT(1:LENSTR(IDENT))
          XNAMEGIVEN = .TRUE.
          WRITE (IOUT,FMT=8004)CRYSTALNAME(1:LENSTR(CRYSTALNAME))
          IF(ONLINE)WRITE (ITOUT,FMT=8004)
     $         CRYSTALNAME(1:LENSTR(CRYSTALNAME))
 8004     FORMAT(' Warning: No CRYSTAL NAME given by KeyWord XNAME',/,
     $         ' It has been set to "',A,'"',/)
        END IF
        IF (.NOT. DNAMEGIVEN ) THEN
          DATASETNAME = 'Unspecified'
c          CALL CCPDAT(MOSDATE)
c          CALL UTIME(MOSTIME)
c          DO 971 I=1,8
c            IF(MOSDATE(I:I).EQ.'/')MOSDATE(I:I) = '_'
c            IF(MOSDATE(I:I).EQ.' ')MOSDATE(I:I) = '0'
c 971      ENDDO
c          WRITE(DATASETNAME,FMT=8002)MOSDATE,MOSTIME
c 8002     FORMAT(A8,':',A8)
          DNAMEGIVEN = .TRUE.
          WRITE (IOUT,FMT=8003)DATASETNAME(1:LENSTR(DATASETNAME))
          IF(ONLINE)WRITE (ITOUT,FMT=8003)
     $         DATASETNAME(1:LENSTR(DATASETNAME))
 8003     FORMAT(/,
     $         ' Warning: No DATASET NAME given by KeyWord DNAME',
     $          /,' It has been set to "',A,'"',/)
        END IF
        DO 971 I=1,MCOLS+3*MAXLATT
          XNAME_COLS(I) = CRYSTALNAME
          DNAME_COLS(I) = DATASETNAME
 971    ENDDO
C
C---- open besthkl file if necessary here
C
        IF(LBEST.AND..NOT.BESTOPEN)THEN
          CALL CCPOPN(-BESTHKL,BSTHKL,1,1,80,IFAIL)
          BESTOPEN = .TRUE.
        ENDIF
          

C
chrp 17042002
C
          IF (.NOT.MULTISEG) THEN
             IF (MTZOPEN) THEN
C
C---- I see no reason why the MTZ file is closed here. If there is more than
C     one set of PROCESS/GO keywords this results in the MTZ file being
C     closed after the first set of PROCESS/GO images have been integrated
C     and the MTZ file is then overwritten by the second set !
C     Comment out the file close here. 
C     Leave the BEST closure and I'm not sure what is required.
C     AGWL 22/12/09
C
CAL                MTZPRT = 1
C     *********************
          if (debug(52))then
            write(iout,fmt=7750)nrun
            if (online) write(itout,fmt=7750)nrun
 7750       format(1x,'**NOT closing MTZ file in control 6 ***',
     +       ' NRUN=',I3)
          end if
CAL                CALL LWCLOS(MTZOUT,MTZPRT)
CAL                MTZOPEN = .FALSE.
                IF(LBEST.AND.BESTOPEN)THEN
                  CLOSE(BESTHKL)
                  BESTOPEN = .FALSE.
                ENDIF
C     *********************
             END IF
C     
C---- If writing to multiple MTZ files, set up the filename here
C     
c             IF (MULTIMTZ.AND..NOT.LONGERNAME) THEN
     
             IF (MULTIMTZ) THEN
               LONGERNAME = .TRUE.
               NCH = LENSTR(MTZNAM)
               DO 970 I = NCH,1,-1
                 IF (MTZNAM(I:I).EQ.'.') THEN
                   MTZNAM = MTZNAM(1:I-1)//'_001'//MTZNAM(I:NCH)
                   GOTO 972
                 END IF
 970           CONTINUE
            END IF
 972        CONTINUE

C
C---- Unfortunately NRUN is not being properly incremented, so with
C     two successive PROCESS/GO keywords, NRUN is still 1 for the second
C     one and so STARTMTZ is called again, deleting previous batch
C     headers. See if we can prevent this by checking for whether an 
C     MTZ file is already open.  AGWL 22 Dec 09
C
            IF (.NOT.MTZOPEN) THEN
              IF (DEBUG(52)) THEN
                write(iout,fmt=7760) MTZOPEN,NRUN
                if (online) write(itout,fmt=7760) MTZOPEN,NRUN
 7760           format(1x,'** In CONTROL call STARTMTZ for mtz file 1 ',
     +           '**  MTZOPEN',L2,' nrun=',I3)
              END IF
              CALLED = 'STARTMTZ'
              COMMENT = 'To open new MTZ file '
              IF (DEBUG(99)) 
     +          CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
              CALL STARTMTZ
            END IF
            NLSUM1 = 0
            NLSUM2 = 0
         END IF
      ELSE
          IF ((.NOT.MULTISEG).AND.(IHKLOUT.EQ.1)) THEN
            IF (MTZOPEN) THEN
              MTZPRT = 1
C     *********************
              if (debug(52))then
                write(iout,fmt=7752) NRUN
                if (online) write(itout,fmt=7752)NRUN
 7752           format(1x,'**closing MTZ file in control 7 ***',
     +             ' NRUN+',I3)
              end if
              CALLED = 'LWCLOS'
              COMMENT = 'To close MTZ file '
              IF (DEBUG(99)) 
     +          CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
              CALL LWCLOS(MTZOUT,MTZPRT)
              MTZOPEN = .FALSE.
              IF(LBEST.AND.BESTOPEN)THEN
                CLOSE(BESTHKL)
                BESTOPEN = .FALSE.
              ENDIF
C     *********************
            END IF
C     
C---- If writing to multiple MTZ files, set up the filename here
C     
            IF (MULTIMTZ.AND..NOT.LONGERNAME) THEN
              NCH = LENSTR(MTZNAM)
              LONGERNAME = .TRUE.
              DO 974 I = NCH,1,-1
                IF (MTZNAM(I:I).EQ.'.') THEN
                  MTZNAM = MTZNAM(1:I-1)//'_001'//MTZNAM(I:NCH)
                  GOTO 976
                END IF
 974          CONTINUE
            END IF
 976        CONTINUE
            if (debug(52)) then
              write(iout,fmt=7762) NRUN
              if (online) write(itout,fmt=7762) NRUN
 7762         format(1x,'** In CONTROL call STARTMTZ to open mtz file',
     +              ' 2 **  NRUN=',I3)
            end if
            CALLED = 'STARTMTZ'
            COMMENT = 'To open new MTZ file '
            IF (DEBUG(99)) 
     +        CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
            CALL STARTMTZ
CAL            NLSUM1 = 0
CAL            NLSUM2 = 0
          END IF
          NLSUM1 = 0
          NLSUM2 = 0
        END IF
C     
chrp 17.07.2006 - moved to just after calls to pstart & start        GENOPEN = .TRUE.
C     
C---- Calculate effective GAIN for film for SD calculations
C     
        IF (.NOT.IMGP) GAIN = G1OD*G1OD*N1OD/(2.0*(25.0*SCNSZ)**2)
C     
C     
C---- Set flag for orientation of vee films on scanner
C     use NREC to do this. If length of film along scanner x
C     is less than 160mm, assume length of film is around drum,
C     ie the v is parallel to scanner x.
C     
        VALONGX = (VEE .AND. (NREC.LT. (40/SCNSZ)*160))
C     
C     Convert tilt,twist,bulge etc  to correct internal
C     units.
C---- Factor to convert from 1/100th deg to distortion units
        RADEG = 100.0 * RTOD
C     
        IF (VEE) THEN
          IF (CBAR.EQ.0) CBAR = 50
          IF (XTOFD.NE.0) FDIST = 1.0/XTOFD
C     
C---- VBNEG, VBPOS, VVERT not passed, so set zero
C     

          VBNEG = 0.0
          VBPOS = 0.0
          VTILT = TILT*FDIST/RADEG
          VTWIST = TWIST*FDIST/RADEG
          VVERT = 0.0
        ELSE
          CALL SETDIS(ITILT,ITWIST,1)
c     IF (XTOFD.NE.0) FDIST = 1.0/ (XTOFD*RADEG)
c     TWIST = ITWIST*FDIST
c     TILT = ITILT*FDIST
          IF (.NOT.IMGP) THEN
            BULGE = IBULGE*FDIST
          END IF
        END IF
C     
C     
C---- write generate file header and crystal identifier to summary
C     
        WRITE (ISUMMR,FMT=6065) GTITLE(1:LENSTR(GTITLE)),
     +       IDENT(1:LENSTR(IDENT))
 6065   FORMAT (/,1X,'Generate File TITLE:   ',A,/,1X,
     +       'Crystal Identifier:    ',A)
C     
      END IF

      IF ((NPRUN.EQ.0).AND.(MODE.NE.1)) THEN
        WRITE (IOUT,FMT=6158)
        IF (ONLINE) THEN
          WRITE (ITOUT,FMT=6158)
          GO TO 50
        END IF
 6158   FORMAT (//,1X,'******* NO images have been specified on ',
     +       'PROCESS keyword *****')
        STOP
      END IF
C     
C---- Checks on supplied profile boundaries, and conflicting definitions
C     
      IF (LINESET.AND.(HIGHRES.OR.LOWRES)) THEN
        WRITE(IOUT,FMT=6067)
        IF (ONLINE) WRITE(ITOUT,FMT=6067)
        IF (BRIEF) WRITE(IBRIEF,FMT=6067)
 6067   FORMAT(/,1X,'When specifying the number of profiles (using',
     +       ' the PROFILE keyword)',/,1X,'the keywords LOWRES,',
     +       ' HIGHRES, and XLINE/YLINE are mutually exclusive')
        IF (ONLINE) GOTO 50
        STOP
      END IF
      IF (PRSET.AND.(.NOT.HIGHRES).AND.(.NOT.LOWRES)) THEN
        IF ((NXLINE.EQ.0).OR.(NYLINE.EQ.0)) THEN
          WRITE(IOUT,FMT=6066)
          IF (ONLINE) WRITE(ITOUT,FMT=6066)
 6066     FORMAT(/,1X,'***** ERROR *****',/,1X,
     +         'If profile boundaries are being supplied, ',
     +         'they must be given for both X and Y directions',/,1X,
     +         'Keywords XLINE,YLINE')
          IF (ONLINE) GOTO 50
          STOP
        END IF
      END IF
C     
C---- If boundaries have been supplied and there is only one box, turn
C     off variable profile option
C     
      IF (LINESET.AND.(NXLINE.EQ.2).AND.(NYLINE.EQ.2)) 
     +     VARPRO = .FALSE.
C     
C---- In the code below, it is assumed that for image plates the
C     defaults are independant of the site, while for film data
C     there are site dependant defaults (eg size of image, ONEOD, 
C     granularity etc
C     
C     *************************************************************
C---- Set up defaults for image plate processing 
C     *************************************************************
 706  IF (IMGP) THEN
C     
C     Default distortion parameters
        IF (YSCAL.LT.-100) YSCAL = 1.0
        IF (XTOFRA.LT.-100) XTOFRA = 1.0
C     
C---- Spiral scan distortion parameters. Turn on by default
C     for Mar and DIP2000 scanners,  unless FIX keyword given.
C     
        IF (SPIRAL) THEN
          DO 695 I = 8,11
            IF (IFIX(I).EQ.0) FIXPAR(I) = .FALSE.
 695     CONTINUE
          IF (INODES.EQ.0) NODES = 1
        ELSE
C
C---- Not a spiral image plate scanner, make sure ROFF,TOFF,RDROFF,RDTOFF
C     parameters are all fixed (they could have been set FREE by a keyword !).
C
          DO 699 I = 8,11
            FIXPAR(I) = .TRUE.
            IFIX(I) = 0
 699      CONTINUE
        END IF
C     
C---- If OVERLOAD CUTOFF has been set, but CUTOFF has NOT been changed
C     on
C     the PROFILE keyword, then set the PROFILE CUTOFF to the smaller
C     of its default value and that given on OVERLOAD CUTOFF
C     
        IF ((ICUT.NE.0).AND.(IPRCUT.EQ.0).AND.(PRCUTOFF.GT.CUTOFF))
     +       THEN
          PRCUTOFF = CUTOFF
          IPRCUT = 1
        END IF
C     
C---- Set default gain, image size,extension if not already set 
C     by keywords GAIN, SIZE, EXTENSION
C     
        IF (SITE.EQ.'LMB') THEN
          IF (IGAIN.EQ.0) THEN
            IGAIN = 1
            GAIN = 1.0
          ENDIF
          IF (NREC.EQ.0) NREC = 1187
          IF (IYLEN.EQ.0) IYLEN = 1187
          IF (ICUT.EQ.0) THEN
            CUTOFF = 31999
          ELSE
            IF (CUTOFF.GT.31999) THEN
              WRITE(IOUT,FMT=6900) CUTOFF
              IF (ONLINE) WRITE(ITOUT,FMT=6900) CUTOFF
            END IF
          END IF
          IF (IPRCUT.EQ.0) THEN
            PRCUTOFF = 31999
          ELSE
            IF (PRCUTOFF.GT.31999) THEN
              WRITE(IOUT,FMT=6902) PRCUTOFF
              IF (ONLINE) WRITE(ITOUT,FMT=6902) PRCUTOFF
            END IF
          END IF
C     
C---- Change adc offset to 16 for LMB prototype scanner
          IF (IDIVIDE.EQ.8.AND.IIDIVIDE.EQ.0) IDIVIDE = 16
          IF (ODEXT(1:1).EQ.' ') ODEXT = 'corr'
C     
C---- CHESS Fuji scanner
C     
        ELSE IF (SITE.EQ.'CHES') THEN
          IF (MACHINE.EQ.'CCD1') THEN
            IF (IGAIN.EQ.0) THEN
              IGAIN = 1
              GAIN = 0.3
            ENDIF
          ELSE IF (MACHINE.EQ.'FUJI') THEN
            IF (IGAIN.EQ.0) THEN
              IGAIN = 1
              GAIN = 1.0
            ENDIF
            IF (ICUT.EQ.0) THEN
              CUTOFF = 9999
            ELSE
              IF (CUTOFF.GT.128000) THEN
                WRITE(IOUT,FMT=6904) CUTOFF
                IF (ONLINE) WRITE(ITOUT,FMT=6904) CUTOFF
 6904           FORMAT(1X,'Cutoff of',I8,' is probably too high',
     +               ' for this type of scanner')
              END IF
            END IF
            IF (IPRCUT.EQ.0) THEN
              PRCUTOFF = 9999
            ELSE
              IF (PRCUTOFF.GT.100000) THEN
                WRITE(IOUT,FMT=6906) PRCUTOFF
                IF (ONLINE) WRITE(ITOUT,FMT=6906) PRCUTOFF
 6906           FORMAT(1X,'Profile cutoff of',I8
     $               ,' is probably too high for this'
     $               ,' type of scanner')
              END IF
            END IF
C     
C---- Change adc offset to 1 for Fuji scanners
C     
            IF(IIDIVIDE.EQ.0)IDIVIDE = 1
          END IF
C     
C---- EMBL Outstation Hamburg
C     
        ELSE IF (SITE.EQ.'EMBL') THEN
          Hbeamline='EMBL'
          IF (IGAIN.EQ.0) THEN
            IGAIN = 1 
            GAIN = 1.0
          ENDIF
C     
C---- SCR1 scanner
C     
          IF (SCANNER.EQ.'SCR1') THEN
            IF (NREC.EQ.0) NREC = 1187
            IF (IYLEN.EQ.0) IYLEN = 1187
            IF (ICUT.EQ.0) THEN
              CUTOFF = 31999
            ELSE
              IF (CUTOFF.GT.31999) THEN
                WRITE(IOUT,FMT=6900) CUTOFF
                IF (ONLINE) WRITE(ITOUT,FMT=6900) CUTOFF
 6900           FORMAT(1X,'******* WARNING *******',/,1X,
     +               'Cutoff value of',I10,' is not ',
     +               'appropriate for this type of scanner. ',
     +               'Should be 31999 *******')
              END IF
            END IF
            IF (IPRCUT.EQ.0) THEN
              PRCUTOFF = 31999
            ELSE
              IF (PRCUTOFF.GT.31999) THEN
                WRITE(IOUT,FMT=6902) PRCUTOFF
                IF (ONLINE) WRITE(ITOUT,FMT=6902) PRCUTOFF
 6902           FORMAT(1X,'******* WARNING *******',/,1X,
     +               'Profile cutoff value of',I10,' is not ',
     +               'appropriate for this type of scanner. ',
     +               'Should be 31999 *******')
              END IF
            END IF
C     
C---- SCR2 scanner
C     
          ELSE IF (SCANNER.EQ.'SCR2') THEN
            IF (NREC.EQ.0) NREC = 1187
            IF (IYLEN.EQ.0) IYLEN = 1187
            IF (ICUT.EQ.0) THEN
              CUTOFF = 31999
            ELSE
              IF (CUTOFF.GT.31999) THEN
                WRITE(IOUT,FMT=6900) CUTOFF
                IF (ONLINE) WRITE(ITOUT,FMT=6900) CUTOFF
              END IF
            END IF
            IF (IPRCUT.EQ.0) THEN
              PRCUTOFF = 31999
            ELSE
              IF (PRCUTOFF.GT.31999) THEN
                WRITE(IOUT,FMT=6902) PRCUTOFF
                IF (ONLINE) WRITE(ITOUT,FMT=6902) PRCUTOFF
              END IF
            END IF
C     
C---- SCR3 scanner
C     
          ELSE IF (SCANNER.EQ.'SCR3') THEN
            IF (NREC.EQ.0) NREC = 1187
            IF (IYLEN.EQ.0) IYLEN = 1187
            IF (ICUT.EQ.0) THEN
              CUTOFF = 31999
            ELSE
              IF (CUTOFF.GT.31999) THEN
                WRITE(IOUT,FMT=6900) CUTOFF
                IF (ONLINE) WRITE(ITOUT,FMT=6900) CUTOFF
              END IF
            END IF
            IF (IPRCUT.EQ.0) THEN
              PRCUTOFF = 31999
            ELSE
              IF (PRCUTOFF.GT.31999) THEN
                WRITE(IOUT,FMT=6902) PRCUTOFF
                IF (ONLINE) WRITE(ITOUT,FMT=6902) PRCUTOFF
              END IF
            END IF
C     
          ELSE IF (MACHINE.EQ.'MAR ') THEN
            IF (NREC.EQ.0) NREC = 1200
            IF (IYLEN.EQ.0) IYLEN = 1200
          END IF
        ELSE IF (SITE.EQ.'DLAB') THEN
          IF (MACHINE.EQ.'MAR ') THEN
            IF (IGAIN.EQ.0) THEN
              IGAIN = 1
              GAIN = 1.0
            ENDIF
            IF (NREC.EQ.0) NREC = 1200
            IF (IYLEN.EQ.0) IYLEN = 1200
          ELSE IF (MACHINE.EQ.'RAXI') THEN
            IF (IGAIN.EQ.0) THEN
              IGAIN = 1
              GAIN = 5.0
            ENDIF
C     AL            IF (NHEAD.EQ.-999) NHEAD = 1
            IF (IPRCUT.EQ.0) PRCUTOFF = 128000
            IF (ICUT.EQ.0) CUTOFF = 250000
C     
C---- YSCAL based on pixel sizes of 101.7mu in fast direction, 105mu in 
C     slow.
            YSCAL = 105.0/101.7
            YSCALIN = YSCAL
          END IF
        ELSE
C     
C---- No SITE given, set defaults for commercial Mar scanner if
C     no value supplied by keyword.
C     
C     
C     ***** machine specific code follows *****
C     
          IF (MACHINE(1:3).EQ.'MAR ') THEN
            IF (IGAIN.EQ.0) THEN
              IGAIN = 1 
              GAIN = 1.0
            ENDIF
C     AL            IF (NHEAD.EQ.-999) NHEAD = 1
            IF (NREC.EQ.0) NREC = 1200
            IF (IYLEN.EQ.0) IYLEN = 1200
C     
C---- Mar345 has dynamic range of 65K for 100mu pixel, 130K for 150mu
C     pixel. After correction, these values can be bigger.
C     
            IF (MODEL(1:4).EQ.'M345') THEN
              IF (RAST.GT.0.149) THEN
                IF (ICUT.EQ.0) CUTOFF = 150000
                IF (IPRCUT.EQ.0) PRCUTOFF = 100000
              ELSE
                IF (ICUT.EQ.0) CUTOFF = 65000
                IF (IPRCUT.EQ.0) PRCUTOFF = 50000
              END IF
            END IF

C     
C---- If MACHINE LMB has been specified set up defaults
C     
          ELSE IF (MACHINE.EQ.'LMB') THEN
            IF (IGAIN.EQ.0) THEN
              IGAIN = 1
              GAIN = 1.0
            ENDIF
            IF (IPRCUT.EQ.0) PRCUTOFF = 262000
            IF (ICUT.EQ.0) CUTOFF = 262000
            IF (NREC.EQ.0) NREC = 3000
            IF (IYLEN.EQ.0) IYLEN = 3000
            YSCAL = 1.0

C     
C---- If ESRF CCD has been specified set up defaults
C     
          ELSE IF (MACHINE.EQ.'CCD2') THEN
            IF (IGAIN.EQ.0) THEN
              IGAIN = 1
              GAIN = 1.0
            ENDIF
            IF (IPRCUT.EQ.0) PRCUTOFF = 50000
            IF (ICUT.EQ.0) CUTOFF = 65000
            YSCAL = 1.0
C     
C---- If MACHINE DIP2 has been specified set up defaults
C     
          ELSE IF (MACHINE.EQ.'DIP2') THEN
            IF (IGAIN.EQ.0) THEN
              IGAIN = 1
              GAIN = 1.0
            ENDIF
            IF (IPRCUT.EQ.0) PRCUTOFF = 62000
            IF (ICUT.EQ.0) CUTOFF = 500000
            IF (NREC.EQ.0) NREC = 2500
            IF (IYLEN.EQ.0) IYLEN = 2500
            YSCAL = 1.0
C     
C---- If MACHINE RAXIS has been specified set up defaults
C     
          ELSE IF (MACHINE.EQ.'RAXI') THEN
            IF ((MODEL.EQ.'RAXISIV').OR.(MODEL.EQ.'RAXISV')
     $           .OR.(MODEL.EQ.'RAXISHTC')) THEN
C     
C---- RAXIS IV defaults (also okay for RAXIS V)
C     
              IF (IGAIN.EQ.0) THEN
                IGAIN = 1
                GAIN = 1.0
              ENDIF
              IF (IPRCUT.EQ.0) PRCUTOFF = 250000
              IF (ICUT.EQ.0) CUTOFF = 1000000
              IF ((IYSCAL.EQ.0).AND.(IPIXY.EQ.0)) YSCAL = 1.0
              IF ((MODEL.EQ.'RAXISHTC').AND.(IYSCAL.EQ.0).AND.
     +                  (IPIXY.EQ.0)) THEN
                     YSCAL = 0.995
                     YSCALIN = YSCAL
              END IF
            ELSE
C     
C---- RAXIS II defaults
C     
              IF (IGAIN.EQ.0) THEN
                IGAIN = 1
                GAIN = 5.0
              ENDIF
              IF (IPRCUT.EQ.0) PRCUTOFF = 128000
              IF (ICUT.EQ.0) CUTOFF = 250000
C     
C---- YSCAL based on pixel sizes of 101.7mu in fast direction, 105mu in 
C     slow.
              IF ((IYSCAL.EQ.0).AND.(IPIXY.EQ.0)) YSCAL = 105.0/101.7
              YSCALIN = YSCAL
            END IF
          END IF
        END IF
C     
C---- NWORD is number of 2-byte words in a strip of ods for film data,
C     or number of pixels in a row for image plate data
C     NBYTE is the number of bytes in a strip
C     
        NWORD = IYLEN
        NBYTE = IYLEN*2
        ITHRESHF = THRESHF
C     
C***********************************END OF IMGP DEFAULTS **********
C***********************************END OF IMGP DEFAULTS **********
      ELSE
C     
C-----Site dependant defaults for film processing
C     
        IF (SITE.EQ.'LMB') THEN
          ONEFILE= .TRUE.
          FHEADER = .FALSE.
C     
C     Image size (number of records and number of pixels per record)
          IF (NREC.EQ.0) NREC = 2260
          IF (IYLEN.EQ.0) IYLEN = 2300
C     
C     
C     Set default film characteristics (if not already set by keyword)
          IF (N1OD.EQ.0) N1OD = 85
          IF (G1OD.EQ.0.0) G1OD = 3.7
          IF (BASEOD.EQ.0.0) BASEOD = 0.02
          IF (CURV.EQ.0.0) CURV = 0.07
        ELSE IF (SITE.EQ.'DLAB') THEN
          ONEFILE= .TRUE.
          FHEADER = .FALSE.
          IF (NREC.EQ.0) NREC = 2400
          IF (IYLEN.EQ.0) IYLEN = 2400
C     
C     
C     Set default film characteristics (if not already set by keyword)
          IF (N1OD.EQ.0) N1OD = 128
          IF (G1OD.EQ.0.0) G1OD = 3.7
          IF (BASEOD.EQ.0.0) BASEOD = 0.02
          IF (CURV.EQ.0.0) CURV = 0.07
        ELSE IF (SITE.EQ.'IMPC') THEN
          ONEFILE= .FALSE.
          FHEADER = .FALSE.
          IF (NREC.EQ.0) NREC = 2496
          IF (IYLEN.EQ.0) IYLEN = 2496
C     
C     
C     Set default film characteristics (if not already set by keyword)
          IF (N1OD.EQ.0) N1OD = 128
          IF (G1OD.EQ.0.0) G1OD = 2.5
          IF (BASEOD.EQ.0.0) BASEOD = 0.02
          IF (CURV.EQ.0.0) CURV = 0.07
        END IF
C     ***
C     *** End of site specific set up
C     ***
C     *************************************************************
C---- Set up defaults for film processing that are site independant
C     *************************************************************
C     
C     
C     
C     Default distortion parameters
        IF (YSCAL.LT.-100) YSCAL = 0.9985
        IF (XTOFRA.LT.-100) XTOFRA = 1.0
C     
C---- NWORD is number of 2-byte words in a strip of ods for film data,
C     or number of pixels in a row for image plate data
C     NBYTE is the number of bytes in a strip
C     
        NBYTE = IYLEN
        NWORD = IYLEN/2        
C     
C     Fiducial coordinates..these are permuted according to the
C     orientation
C     of the film. The values below are for the Enraf Nonius Arndt
C     Wonacott
C     camera             
        IF (NFID.EQ.0) THEN
          NFID = 3  
          FIDXY(1,1) = 40.0
          FIDXY(1,2) = -50.0
          FIDXY(2,1) = -40.0
          FIDXY(2,2) = -50.0
          FIDXY(3,1) = -40.0
          FIDXY(3,2) =  50.0
        END IF
C     
C     Correct for orientation  
        IF (ROTATED) THEN
          DO 700 I=1,3
            TEMP = FIDXY(I,1)
            FIDXY(I,1) = FIDXY(I,2)
            FIDXY(I,2) = -TEMP
 700      CONTINUE
        END IF
C     
C     Fiducial threshold
C     
        IF (THRESHF.EQ.0) THRESHF = 1.0
        ITHRESHF = THRESHF*N1OD

      END IF
C     
C*********END of film specific defaults *********************
C*********END of film specific defaults *********************
C     
C     Set up default fiducial search box size, convert to half width in
C     10 mu
      IF (XMMF.EQ.0.0) XMMF = 10.0
      MM = NINT(100.0*XMMF/2.0)
      MMDB = NINT(100.0*XMMDB/2)
C     
C---- Convert FIDXY to 10 micron units
C     
      DO 702 I = 1,NFID
        DO 704 J = 1,2
          FSPOS(I,J) = FIDXY(I,J)*100
C     AL          IF (FSPOS(I,J).GT.9000) THEN
C     AL            IF (ONLINE) WRITE (ITOUT,FMT=6300) I,(FIDXY(I,K),K=1
C     ,2)
C     AL            WRITE (IOUT,FMT=6300) I, (FIDXY(I,K),K=1,2)
C     AL            IF (.NOT.ONLINE) STOP
C     AL          END IF
 704    CONTINUE
 702  CONTINUE
 6300 FORMAT (/,1X,'**** FIDUCIAL ',I2,' Has coordinates off the edge',
     +     ' of the film !! (',I8,',',I8,')',/,1X,'are coordinates su',
     +     'pplied in mm ?')

C     
C---- Set limits of scanning (used in gensort)
C     
      XSCMIN = 1
      XSCMAX = NREC
C     
C---- Deal with refinement residual limits. If these have been supplied
C     on
C     a keyword, must assign them to the correct variable depending on 
C     whether weighted or weighted refinement is being done. If weighted
C     refinement, the RESID applies to the weighted residual limit etc
C     convert unweighted residual limit from mm to 10 micron units
C     
      IF (RRSET) THEN
        IF (RWEIGHT) THEN
          WRMSLIM = XRMSLIM
        ELSE
          RMSLIM = 100.0*RMSLIM
        END IF
      END IF
      IF (ARRSET) THEN
        IF (RWEIGHT) THEN
          AWRMSLIM = AXRMSLIM
        ELSE
          ARMSLIM = 100.0*AXRMSLIM
        END IF
      END IF
C     
C     Check image size
C     
      IF (NWORD.GT.IYLENGTH) THEN
        IF (ONLINE) WRITE (ITOUT,FMT=7222) NWORD,IYLENGTH,NWORD
        WRITE (IOUT,FMT=7222) NWORD,IYLENGTH,NWORD
        STOP
      ELSE IF (NREC.GT.IXWDTH/2 .AND. INCORE) THEN
        IF (ONLINE) WRITE (ITOUT,FMT=7220) 2*NREC,IXWDTH,2*NREC
        WRITE (IOUT,FMT=7220) 2*NREC,IXWDTH,2*NREC
        STOP
      END IF
C     
C     If SUMPARTIALS given, check that this is image plate and IXWDTH
C     and
C     INCORE are true.
C     
      IF (SUMPART) THEN
        IF (.NOT.IMGP) THEN
          WRITE(IOUT,6306)
          IF (ONLINE) WRITE(ITOUT,6306)
          STOP
        END IF
        IF (.NOT.INCORE) THEN
          WRITE(IOUT,6308)
          IF (ONLINE) WRITE(ITOUT,6308)
          STOP
        END IF
        IF (IXWDTH.LT.2*NREC) THEN
          WRITE(IOUT,6309) 2*NREC
          IF (ONLINE) WRITE(ITOUT,6309) 2*NREC
          STOP
        END IF
      END IF
 6306 FORMAT(/,1X,'** Can only sum partials with image plate data **')
 6308 FORMAT(//,1X,'If summing partials, the images MUST be stored',
     +     ' in memory (Keyword INCORE)')
 6309 FORMAT(//,1X,'If summing partials, IXLENGTH must be at',
     +     ' least twice the number of strips, ie',I7,/,1X,
     +     'Change parameter IXLENGTH in common block PARAMETER',
     +     ' and recompile and link program')
C
C---- Change default PKRATIO rejection criterion for Mar 3x3 mosaic - assume 
C     4x4 is the same?
C
      IF ((MACHINE.EQ.'MARC').AND.
     $     ((MODEL.EQ.'MOSAIC3').OR.(MODEL.EQ.'MOSAIC4')).AND.
     +     (IPKRAT.EQ.0)) PKRAT = 5.0
C     
C     
C---- Start reflecting input and defaults 
C     
      IF (MULTISEG.AND.(.NOT.FIRSTTIME).AND.(MODE.NE.4)) RETURN
      IF (MODE.EQ.1) RETURN
C     
C     
      WRITE(IOUT,FMT=6800) WAVE,2.0*DIVH/DTOR,2.0*DIVV/DTOR,
     +     DELAMB,DELCOR
      IF (ONLINE)WRITE(ITOUT,FMT=6800)WAVE,2.0*DIVH/DTOR,2.0*DIVV/DTOR,
     +     DELAMB,DELCOR
 6800 FORMAT(//,1X,'Beam Parameters',/,1X,'===============',/,1X,
     +     'Wavelength',F7.4,/,1X,'Beam divergence:  Horizontal',
     +     ' (DIVH)',F6.3,'   Vertical (DIVV)',F6.3,/,1X,
     +     'Wavelength dispersion (DISP)',
     +     F8.5,' Correlated Del-Lambda (DELCOR)',F8.5 )
C
C---- Set logical VERTICAL to TRUE for Pilatus and Eiger detectors that have
C     OMEGA values 90 or 270 as this corresponds to a vertical
C     rotation axis and this affects the polarisation correction.
C
      VERTICAL = (((MACHINE.EQ.'PILA').OR.(MACHINE.EQ.'EIGE').OR.
     +             (MACHINE.EQ.'HDF5'))
     +      .AND.((ABS(OMEGAFD-90.0).LT.1.0)
     +     .OR.(ABS(OMEGAFD-270.0).LT.1.0)))
      IF (NOTVERTICAL) VERTICAL = .FALSE.
      IF (IMONO.EQ.0) THEN
        WRITE(IOUT,FMT=6802)
        IF (ONLINE) WRITE(ITOUT,FMT=6802)
 6802   FORMAT(1X,'The X-ray beam is assumed to be unpolarised (this ',
     +       'is appropriate for pinhole',/,1X,'or double mirror ',
     +       'collimation)')
      ELSE IF (IMONO.EQ.1) THEN
        WRITE(IOUT,FMT=6804)
        IF (ONLINE) WRITE(ITOUT,FMT=6804)
 6804   FORMAT(1X,'The polarisation correction for a ',
     +       'graphite monochromator will be applied.')
      ELSE IF ((IMONO.EQ.2).OR.(IMONO.EQ.3))   THEN
        WRITE(IOUT,FMT=6806) TOR
        IF (ONLINE) WRITE(ITOUT,FMT=6806) TOR
 6806   FORMAT(1X,'The X-ray beam is assumed to have a degree of',
     +       ' polarisation of',F6.3)
      END IF
      IF (VERTICAL) THEN
        WRITE(IOUT,FMT=6805)
        IF (ONLINE) WRITE(ITOUT,FMT=6805)
 6805   FORMAT(1X,'A polarization correction for a vertical rotation',
     +       ' axis will be applied.')
      END IF
C     
C     
      NCH = LENSTR(IDENT)
      IF (NCH.EQ.0) NCH = 1
C     
C---- Don't reflect this on second or later segments of STRATEGY run
C     because if SPEEDUP is used, CELL and AMAT will be wrong !!
C     
      IF (ISTRUN.EQ.0) THEN
        WRITE(IOUT,FMT=6810) IDENT(1:NCH),NUMSPG,
     +       SPGNAM(1:LENSTR(SPGNAM)),
     +       CELL,RCELL,2.0*ETA/DTOR,0.0001*BLKSIZE,
     +       ((UMAT(I,J),J=1,3),I=1,3),((AMAT(I,J),J=1,3),I=1,3),DELPHI
        IF (ONLINE) WRITE(ITOUT,FMT=6810) IDENT(1:NCH),NUMSPG,
     +       SPGNAM(1:LENSTR(SPGNAM)),CELL,RCELL,
     +       2.0*ETA/DTOR, 0.0001*BLKSIZE,((UMAT(I,J),J=1,3),I=1,3),
     +       ((AMAT(I,J),J=1,3),I=1,3),DELPHI
      END IF
 6810 FORMAT (/,1X,'Crystal Parameters',/,1X,'===================',/,
     +     1X,'Crystal identifier (used as template for image file',
     +     ' names (IDENT): ',A,/,1X,'Space group (SYMMETRY)',I5,
     +     1X,'(',A,')',
     +     /,1X,'Real cell parameters (CELL):',/,1X,
     +     6F10.3,/,1X,'Reciprocal cell parameters :',/,1X,3F10.6,
     +     3F10.2,//,1X,'Mosaic spread (MOSAIC)',F6.3,' degrees',
     +     /,1X,'Mosaic block size (BLOCKSIZE)',F7.2,' microns',
     +     /,1X,'Rotation matrix U defining standard setting :',
     +     /3 (12X,3F9.5,/),/,1X,'Orientation Matrix [A], ',
     +     'Components of A*,B*,C*',/'   ',
     +     ' Along X ',3F10.6,/'    Along Y ',3F10.6,/'    Along Z ',
     +     3F10.6,/,1X,'          Misorientation Angles',/5X,3F9.3)
      
C     Set up defaults for image plate data Get rid of "IF (IMGP)" block
C     AL      IF (IMGP) THEN
      FHEADER = .FALSE.
C     
C     Do the NOFID stuff
      OMEGAF = OMEGAFD*DTOR

c      write(*, *) 'OMEGAF = ', omegaf

      OMEGA0 = OMEGAF + CCOM*DTOR
      COSOM0 = COS(OMEGA0)
      SINOM0 = SIN(OMEGA0)
      NOFID = .TRUE.
      ROTATED = .FALSE.
C     
      WRITE(IOUT,FMT=6830) RAST
      IF (ONLINE) WRITE(ITOUT,FMT=6830) RAST
 6830 FORMAT(/,1X,'Detector parameters',/,1X,'===================',
     +     /,1X,'Pixel size in the "slow" (X) direction in the ',
     +     'image (PIXEL)  ',F7.5,' mm')
C     
      IF (IPIXY.NE.0) THEN
        WRITE(IOUT,FMT=6832) RAST/YSCAL
        IF (ONLINE) WRITE(ITOUT,FMT=6832) RAST/YSCAL
      END IF
 6832 FORMAT(1X,'Pixel size in the "fast" direction (PIXEL)  ',
     +     F5.3,' mm')
C     
      WRITE(IOUT,FMT=6831) 0.01*XMIN,0.01*XMAX,0.01*YMIN,0.01*YMAX,
     +     0.01*RMIN,0.01*RMAX
      IF (ONLINE) WRITE(ITOUT,FMT=6831) 0.01*XMIN,0.01*XMAX,
     +     0.01*YMIN,0.01*YMAX,0.01*RMIN,0.01*RMAX
 6831 FORMAT(1X,'Physical limits of detector (LIMITS)',
     +     ' relative to direct beam position',/,1X,
     +     'Minimum and maximum X coordinate (XMIN,XMAX)',
     +     2F7.2,' mm',/,1X,
     +     'Minimum and maximum Y coordinate (YMIN,YMAX)',
     +     2F7.2,' mm',/,1X,'Minimum',
     +     ' and maximum radial coordinates (RMIN,RMAX)',2F7.2,' mm')
C     
C     
C     ***** machine specific code follows *****
C     
      IF (MACHINE.EQ.'MAR ') THEN
        WRITE(IOUT,FMT=6564)
        IF (ONLINE) WRITE(ITOUT,FMT=6564)
      ELSE IF (MACHINE.EQ.'RAXI') THEN
        WRITE(IOUT,FMT=6566)
        IF (ONLINE) WRITE(ITOUT,FMT=6566)
      END IF
 6564 FORMAT(1X,'For Mar scanners, X is horizontal, Y is vertical')
 6566 FORMAT(1X,'For R-axis scanners, X is vertical, Y is ',
     +     ' horizontal')
      IF (CIRCULAR) THEN
        WRITE(IOUT,FMT=6829) 0.01*RSCAN
        IF (ONLINE) WRITE(ITOUT,FMT=6829) 0.01*RSCAN
 6829   FORMAT(/,1X,'Radius of scanned area (RSCAN)',F6.1,'mm')
        IF (RSCANX.NE.0) THEN
          WRITE(IOUT,FMT=6833) 0.01*RSCANX,0.01*RSCANY
          IF (ONLINE) WRITE(ITOUT,FMT=6833) 0.01*RSCANX,0.01*RSCANY
        END IF
 6833   FORMAT(1X,'with a centre (CENTRE) at',2F8.2,'mm relative to',
     +       ' the first pixel in the image')
      ELSE IF (ORTHOG) THEN
        WRITE(IOUT,FMT=6825) 0.01*XSCAN,0.01*YSCAN
        IF (ONLINE) WRITE(ITOUT,FMT=6825) 0.01*XSCAN,0.01*YSCAN
 6825   FORMAT(1X,'Detector limits relative to the physical ',
     +       ' centre of the detector:',/,1X,'in X (XSCAN) ',F7.2,
     +       ' mm, in Y (YSCAN) ',F7.2,'mm')
      END IF
C     
      IF (RMINX.NE.0) THEN
        WRITE(IOUT,FMT=6835) 0.01*RMINX,0.01*RMINY
        IF (ONLINE) WRITE(ITOUT,FMT=6835) 0.01*RMINX,0.01*RMINY
      END IF
 6835 FORMAT(1X,'The minimum radius (RMIN) is applied with an ',
     +     'origin (RCENTRE) at',/,1X,2F8.2,'mm relative to',
     +     ' the first pixel in the image.')
C     
C     
      IF (DSTMIN.NE.0.0) THEN
        DSTMINP = WAVE/DSTMIN
      ELSE
        DSTMINP = 1000.0
C     
C---- reset DSTMIN (if zero) to correspond to 1000A
C     
        DSTMIN = WAVE/DSTMINP
      END IF
C
C---- If this is a phi profile with restricted resolution, need to 
C     modify the resolution limits. Note that the full resolution
C     limits are used in the generation algorithm (REEK) but only 
C     those within the restricted resolution range are accepted inC
C     SPTEST. This is because restricting the resolution range in
C     REEK itself lead to some rflections not being correctly generated
C     due to a problem in the algorithm.
C
      IF (RESTRESL) THEN
        HIGHRESP = SQRT(WAVE/DSTMAX2PHIP)
        LOWRESP = SQRT(WAVE/DSTMIN2PHIP)
      ELSE
        HIGHRESP = WAVE/DSTMAX
        LOWRESP = DSTMINP
      END IF
C
      WRITE(IOUT,FMT=6837) 0.01*XOFF,0.01*YOFF,0.01*XTOFD,
     +     TWOTHETA,HIGHRESP,LOWRESP
      IF (ONLINE) WRITE(ITOUT,FMT=6837) 0.01*XOFF,
     +     0.01*YOFF,0.01*XTOFD,TWOTHETA,HIGHRESP,LOWRESP
 6837 FORMAT(//,1X,'Detector translations: XOFF',F7.2,' mm   YOFF',
     $     F7.2,' mm',/,1X
     $     ,'These are the distances between the direct beam'
     $     ,' position and the centre of',/,1X,'the image.'//,1X
     $     ,'Crystal to detector distance (DIST)',F7.2,' mm',/,1X
     $     ,'Detector swing angle (TWOTHETA)',F7.2,' degrees',/,1X
     $     ,'Maximum resolution (RESOLUTION or RMAX)',F6.3,'A',/,1X
     $     ,'Low resolution cutoff (RESOL or RMIN and RCENTRE)',F9.2,
     $     'A')
      IF (NEXCL.GT.0) THEN
        WRITE(IOUT,FMT=7152) (RESEXL(I),RESEXH(I),I=1,NEXCL)
        IF (ONLINE) WRITE(ITOUT,FMT=7152) 
     +       (RESEXL(I),RESEXH(I),I=1,NEXCL)
      END IF
 7152 FORMAT(1X,'Spots lying within the following resolution',
     +     ' ranges will not be integrated:',/,
     +     5(1X,F5.2,' to ',F5.2,3X))
C     
      IF (NXYEXC.GT.0) THEN
        WRITE(IOUT,FMT=7153)
        IF (ONLINE) WRITE(ITOUT,FMT=7153)
        DO 705 I=1,NXYEXC
          X1=XYEXC(1,I)*0.01
          Y1=XYEXC(2,I)*0.01
          X2=XYEXC(3,I)*0.01
          Y2=XYEXC(4,I)*0.01
          WRITE (IOUT,FMT=7154) X1,Y1,X2,Y2
          IF (ONLINE) WRITE(ITOUT,FMT=7154) X1,Y1,X2,Y2
 705    CONTINUE
      END IF
 7153 FORMAT(1X,'Spots lying within the rectangular areas defined',
     +     ' by the following corners will not be integrated:')
 7154 FORMAT(1X,F7.2,',',F7.2,' to ',F7.2,',',F7.2,' (mm.)')
C     
      IF (RESCUT.NE.0.0) THEN
        WRITE(IOUT,FMT=6822) RESCUT
        IF (ONLINE) WRITE(ITOUT,FMT=6822) RESCUT
      END IF
 6822 FORMAT(1X,'If the mean I/sigma(I) drops below',F6.2,
     +     ' at any resolution then data beyond',/,1X,'this ',
     +     'resolution will NOT be written to MTZ file')
      WRITE(IOUT,FMT=6824) 0.01*IXSEP,0.01*IYSEP
      IF (ONLINE) WRITE(ITOUT,FMT=6824) 0.01*IXSEP,0.01*IYSEP
 6824 FORMAT(/,1X,'Minimum',
     $     ' spot separation before spots will be classified',
     $     ' as overlapping:',/,1X,'(SEPARATION) In scanner X,Y'
     $     ,' directions',2F6.2,'mm')
C     
C     
      IF (STRATEGY.OR.TESTGEN) GOTO 708
C     
      IF ((ITRIM.NE.0).OR.(NOVERLAP.NE.0)) THEN
        WRITE(IOUT,FMT=6827) 2*ITRIM,NOVERLAP
        IF (ONLINE) WRITE(ITOUT,FMT=6827) 2*ITRIM,NOVERLAP
      END IF
 6827 FORMAT(1X,'When rejecting background pixels due to overlap',
     +     ' of adjacent spots',/,1X,'the peak size of the ',
     +     'neighbouring spots will be reduced by',I2,' pixels',
     +     /,1X,'in both directions (TRIM).',/,1X,'In addition, ',
     +     'there must be at least',I2,' adjacent spots ',
     +     'overlapping',/,1X,'any given pixel before it is ',
     +     'rejected (NOVERLAP).')
C     
      IF (DENSE) THEN
        WRITE(IOUT,FMT=6840)
        IF (ONLINE) WRITE(ITOUT,FMT=6840)
      END IF
 6840 FORMAT(1X,'Pixels overlapped by neighbouring spots will be',
     +     ' determined for each spot',/,1X,'individually ',
     +     '(CLOSE subkeyword).')
C     
C     
      WRITE(IOUT,FMT=6161) GAIN,IDIVIDE,NULLPIX,CURV,
     +     0.01*THICK,XCENMMIN(1),YCENMMIN(1)
C     AL     +                       0.01*THICK,XMM(1),YMM(1)
C     
      IF (ONLINE) WRITE(ITOUT,FMT=6161) GAIN,IDIVIDE,NULLPIX,
     +     CURV,0.01*THICK,XCENMMIN(1),YCENMMIN(1)
C     AL     +                       CURV,0.01*THICK,XMM(1),YMM(1)
 6161 FORMAT(/,1X,'Detector gain (GAIN)',
     +     F5.2,/,1X,'Scanner adc offset (ADCOFFSET)',I3,/,
     +     1X,'Pixels outside active area have values of',I4,
     +     ' (NULLPIX)',/,1X,'Non-linearity correction (quadratic',
     +     ' correction)',E9.3,' (NONLINEARITY)',/,
     +     1X,'Nominal detector thickness (controls expansion',
     +     ' of measurement box) ',F4.2,' mm',/,
     +     1X,'Direct beam coordinates (BEAM) set to',
     +     2F10.3,' mm',/,1X,'Note that these coordinates are wrt ',
     +     'an origin at the first pixel',/,1X,'in the image, with',
     +     ' X the slowly changing direction and Y the',/,1X,'fast ',
     +     'direction in the image.')
C     
C     ***** machine specific code follows *****
C     
      IF (MACHINE.EQ.'MAR ') THEN
        WRITE(IOUT,FMT=6560)
        IF (ONLINE) WRITE(ITOUT,FMT=6560)
 6560   FORMAT(1X,'For the Mar scanner this  ',
     +       'is the lower right corner of the image as',/,1X,'viewe',
     +       'd from behind the detector looking towards the source.',
     +       /,1X,'The scanner X axis is horizontal and the Y axis i',
     +       's vertical')
      ELSE IF (MACHINE.EQ.'RAXI') THEN
        WRITE(IOUT,FMT=6562)
        IF (ONLINE) WRITE(ITOUT,FMT=6562)
      END IF
 6562 FORMAT(1X,'For the R-axis scanner this  ',
     +     'is the lower left corner of the image as',/,1X,'viewed ',
     +     'from behind the detector looking towards the source.',/,
     +     1X,'The scanner X axis is vertical and the Y axis is ',
     +     'horizontal')
C     
C     Set all film centres, modifying Y coordinate by YSCAL. If 
C     detector is swung out, and direct beam coordinates for twotheta=0
C     have been given, correct them for the swing angle
C     
 708  IF (ISWUNG.EQ.1) THEN
        XCEN0 = NINT(100.0*XCENMMIN(1) - 
     +       COS(OMEGAF)*XTOFD*TAN(TWOTHETA*DTOR))
        YCEN0 = NINT(100.0*YCENMMIN(1)*YSCAL - 
     +       YSCAL*SIN(OMEGAF)*XTOFD*TAN(TWOTHETA*DTOR))
      ELSE
        XCEN0 = NINT(100.0*XCENMMIN(1))
        YCEN0 = NINT(100.0*YCENMMIN(1)*YSCAL)
      END IF
C     
C---- Have to do this now rather than earlier in code because it depends
C     on
C     knowing YSCAL. Note that even if different BEAM coords were given
C     for different images, these will now be preserved.
C     
      DO 710 J=1,MAXPAX
        XCENMM(J,1)=XCENMMIN(J)
        YCENMM(J,1)=YSCAL*YCENMMIN(J)
        IF ((TWOTHETA.NE.0.0).AND.(ISWUNG.EQ.0)) THEN
          XCENMM(J,1) = XCENMM(I,1) + 
     +         COS(OMEGAF)*0.01*XTOFD*TAN(TWOTHETA*DTOR)
          YCENMM(J,1) = YCENMM(J,1) + 
     +         SIN(OMEGAF)*YSCAL*0.01*XTOFD*TAN(TWOTHETA*DTOR)
        END IF
 710  CONTINUE

C     Get rid of "IF (IMGP)" block
C     AL      END IF

C     
C---- For Mar image plate data, correct direct beam X coordinate and CCX
C     for
C     inversion of image.
C     Raster size is RAST mm
C     
      IF (INVERTX) THEN
        DO 712 I = 1,MAXPAX
          XCENMM(I,1) = NREC*RAST - XCENMMIN(I)
C     
C---- Because we have to use XCENMMIN here, must apply correction
C     for swung detectors
C     
          IF ((TWOTHETA.NE.0.0).AND.(ISWUNG.EQ.0)) THEN
            XCENMM(I,1) = XCENMM(I,1) + 
     +           COS(OMEGAF)*0.01*XTOFD*TAN(TWOTHETA*DTOR)
          END IF
 712    CONTINUE
        XCEN0 = 100.0*NREC*RAST - XCEN0
        IF (RSCANX.NE.0.0) RSCANX = 100.0*NREC*RAST - RSCANX
        IF (RMINXINP.NE.0.0) RMINX = 100.0*NREC*RAST - RMINXINP
C     
C---- Only change CCX if it was read in from input. It may have been
C     passed from the previous round of a multi-segment post refinement
C     in which case it MUST NOT be reset.
C     
        IF ((ICCX.EQ.1).AND.(.NOT.CCXRESET)) THEN
          CCX = -CCX
          CCXRESET = .TRUE.
        END IF
      END IF
C     
      IF (STRATEGY) THEN
        IF (AUTO) THEN
          WRITE(IOUT,FMT=6384)
          IF (ONLINE) WRITE(ITOUT,FMT=6384)
 6384     FORMAT(//,1X,'Automatic data collection strategy')
          IF (VOLSCAL.NE.1.0) THEN
            WRITE(IOUT,FMT=6385) VOLSCAL
            IF (ONLINE) WRITE(ITOUT,FMT=6385) VOLSCAL
 6385       FORMAT(1X,'To speed up the calculation, the cell ',
     +           'volume will be reduced by a factor of',F6.1)
            CELLSCAL = VOLSCAL**0.333333
          END IF
C     
C---- If using AUTO mode, the total rotation must cover an integral
C     number
C     of "steps" of data. Check this and modify values if necessary
C     
          IF (ROTAUTO.NE.0.0) THEN
            I = NINT(ROTAUTO/PHIINC(NSEGM))
            IF (ABS(I*PHIINC(NSEGM)-ROTAUTO).GT.0.1) THEN
              ROTAUTO = I*PHIINC(NSEGM)
              WRITE(IOUT,FMT=6376) ROTAUTO, PHIINC(NSEGM)
              IF (ONLINE) WRITE(ITOUT,FMT=6376) ROTAUTO,
     +             PHIINC(NSEGM)
 6376         FORMAT(1X,'*** WARNING ***',/,1X,'Total rotation',
     +             ' changed to',F5.0,' degrees so that it',
     +             /,1X,'is a multiple of the step size',
     +             ' (',F4.0,' degrees)')
            END IF
C     
C---- Now check there is at least one step in each segment
C     
            IF (NSEGAUTO*PHIINC(NSEGM).GT.ROTAUTO) THEN
              ROTAUTO = NSEGAUTO*PHIINC(NSEGM)
              WRITE(IOUT,FMT=6378) ROTAUTO,NSEGAUTO,PHIINC(NSEGM)
              IF (ONLINE) WRITE(ITOUT,FMT=6378) ROTAUTO,NSEGAUTO,
     +             PHIINC(NSEGM)
            END IF
 6378       FORMAT(1X,'*** WARNING ***',/,1X,'Total rotation',
     +           ' increased to',F5.0,' degrees so that there',
     +           /,1X,'is at least one step in each of the',I3,
     +           ' segments requested',/,1X,'Current step size is',
     +           F4.0,' degrees (change with STEP keyword)')
          END IF
          IF (ROTAUTO.NE.0) THEN
            WRITE(IOUT,FMT=6386) ROTAUTO
            IF (ONLINE) WRITE(ITOUT,FMT=6386) ROTAUTO
 6386       FORMAT(1X,'The angular rotation will be limited',
     +           ' to ',F5.1,' degrees.')
          END IF
          IF (NSEGAUTO.GT.1) THEN
            IF (SIZESET) THEN
              WRITE(IOUT,FMT=6390) NSEGAUTO,
     +             (PHISEGA(I),I=1,NSEGAUTO)
              IF (ONLINE) WRITE(ITOUT,FMT=6390) NSEGAUTO,
     +             (PHISEGA(I),I=1,NSEGAUTO)
            ELSE
              WRITE(IOUT,FMT=6388) NSEGAUTO
              IF (ONLINE) WRITE(ITOUT,FMT=6388) NSEGAUTO
            END IF
 6388       FORMAT(1X,'The rotation will be split up into',I3,
     +           ' segments of approximately equal size')
 6390       FORMAT(1X,'The rotation will be split up into',I3,
     +           ' segments with sizes (degrees):',8F6.1)
          END IF
C     
        ELSE
C     
C---- NOT in AUTO mode
C     
          WRITE(IOUT,FMT=6370) NSEGM
          IF (ONLINE) WRITE(ITOUT,FMT=6370) NSEGM
 6370     FORMAT(//,1X,'Data collection strategy ',/,1X,
     +         '========================',/,1X,'Reflections',
     +         ' will be generated in',I4,' segments as listed below',
     +         /,1X,'   Phi start     Phi end   step size')
          K = 0
          DO 714 I = 1,NSEGM
            J = NINT(PHIST(I))/360 + 1
            IF (J.NE.K) THEN
              WRITE(IOUT,FMT=6372) J
              IF (ONLINE) WRITE(ITOUT,FMT=6372) J
 6372         FORMAT(1X,'Run number',I3)
              K = J
            END IF
            PHI1 = PHIST(I) - (J-1)*360 - PHIADD(I)
            PHI2 = PHIFIN(I) - (J-1)*360 - PHIADD(I)
            WRITE(IOUT,FMT=6374) PHI1,PHI2,PHIINC(I)
            IF (ONLINE) WRITE(ITOUT,FMT=6374) PHI1,PHI2,PHIINC(I)
 6374       FORMAT(1X,3F12.1)
 714      CONTINUE
          IF (VOLSCAL.NE.1.0) THEN
            WRITE(IOUT,FMT=6385) VOLSCAL
            IF (ONLINE) WRITE(ITOUT,FMT=6385) VOLSCAL
            CELLSCAL = VOLSCAL**0.333333
          END IF
        END IF
C     
C---- If called from dispay window, which had MODE=10, reset to 0
C     
        IF (MODE.EQ.10) MODE = 0
C     
C---- Convert spot separations (IXSEP,IYSEP) into "ideal detector"
C     coordinate
C     frame, as the spot coordinates (generate file coords) are in this
C     frame
C     
        MINDTX = NINT(ABS(IXSEP*COSOM0 + IYSEP*SINOM0))
        MINDTY = NINT(ABS(IYSEP*COSOM0 + IXSEP*SINOM0))
        RETURN
      END IF
C     
      IF (TESTGEN) THEN
C     
C---- Convert spot separations (IXSEP,IYSEP) into "ideal detector"
C     coordinate
C     frame, as the spot coordinates (generate file coords) are in this
C     frame
C     
        MINDTX = NINT(ABS(IXSEP*COSOM0 + IYSEP*SINOM0))
        MINDTY = NINT(ABS(IYSEP*COSOM0 + IXSEP*SINOM0))
C     
C---- Get cell parameters
C     
        ICHECK = 1
        IF (ICELL.EQ.1) ICHECK = 0
C     
C     ************************
        CALL SETMAT(IMAT,IUMAT,ICELL,ICHECK)
C     ************************
        RETURN
      END IF
C     
C---- PRINT INFO
C     
C     Removed processing pack.... from here
C*********In here
      IF (NDIR.EQ.0) THEN
        WRITE(IOUT,FMT=6381) ODEXT
        IF (ONLINE) WRITE(ITOUT,FMT=6381) ODEXT
      ELSE
        WRITE(IOUT,FMT=6380) NDIR,
     +       (FDISK(NNDIR)(1:LENSTR(FDISK(NNDIR))),NNDIR = 1,NDIR)
        WRITE(IOUT,FMT=6383) ODEXT
        IF (ONLINE) THEN
          WRITE(ITOUT,FMT=6380) NDIR,
     +         (FDISK(NNDIR)(1:LENSTR(FDISK(NNDIR))),NNDIR = 1,NDIR)
          WRITE(ITOUT,FMT=6383) ODEXT
        END IF
      END IF
 6380 FORMAT(/,1X,'Images will be read from the following',
     +     I3,' directories',/,1X,'(DIRECTORY...up to 10 can be',
     +     ' given):',/,(1X,A))
 6383 FORMAT(1X,'The filename extension (EXTENSION) is:',A)
 6381 FORMAT(/,1X,'Images will be read from the local directory',
     +     /,1X,'The filename extension is:',1X,A)
C     
C     
      IF (EFAC.GT.-900.0) THEN
        WRITE (IOUT,FMT=6152) EFAC
 6152   FORMAT (/,1X,'Scanner instrument error factor set to ',F6.3)
        IF (ONLINE) WRITE (ITOUT,FMT=6152) EFAC
        IF (PROFILE) THEN
          WRITE(IOUT,FMT=6155)
          IF (ONLINE) WRITE(ITOUT,FMT=6155)
        END IF
 6155   FORMAT(1X,'This value will override that calculated by the',
     +       ' program')
      ELSE
        IF (PROFILE) THEN
          WRITE (IOUT,FMT=6153)
 6153     FORMAT (/,1X,'Scanner instrument error factor will be ',
     +         'calculated by the program')
          IF (ONLINE) WRITE (ITOUT,FMT=6153)
        ELSE
          WRITE(IOUT,FMT=6157)
          IF (ONLINE) WRITE(ITOUT,FMT=6157)
 6157     FORMAT(/,1X,'Scanner error cannot be calculated without',
     +         ' profile fitting',/,1X,'so scanner error instrument',
     +         ' factor set to 0.0')
        END IF
      END IF
C     
C---- If NOFID is set, check that BEAM has been set for all packs
C     

c     gw hack to prevent beam complaints...
c     graeme used 'if(socklo)' Harry suggests 'if(nofid)', because this
C     also has to work with background runs
C 
c     this will be removed unless it proves that it is necessary...

c     if(nofid) then
c     do 1701 i = 1, 3
c     do 1702 j = 1, maxpax
c     xcenmm(j,i) = xmm(i)
c     ycenmm(j,i) = ymm(i)
c     xcenmmin(j) = xmm(1)
c     ycenmmin(j) = ymm(1)
c     1702       end do
c     1701    end do
c     end if
      IF (NOFID) THEN                       
        IFLAG = 0
        XTEST = 0.0
        IF (INVERTX) XTEST = NREC*RAST
        DO 720 I = IFIRSTPACK,NPACK
          IF ((XCENMM(I,1).EQ.XTEST).AND.(YCENMM(I,1).EQ.0.0)) THEN
            IF (ONLINE) WRITE (ITOUT,FMT=6162) IDPACK(I)
 6162       FORMAT (1X,'*** NO BEAM Parameters ',
     $           'supplied for pack',I4,' ***')
            WRITE (IOUT,FMT=6162) IDPACK(I)
            IFLAG = 1
          END IF
 720    CONTINUE
        IF (IFLAG.EQ.1) GO TO 50
      END IF
C     
C     Print film characteristics, image size, distortion parameters
C     fiducial coordinates
      IF (.NOT.IMGP) THEN
        WRITE(IOUT,FMT=6310) G1OD,BASEOD,CURV,N1OD,0.01*THICK,NFID,
     +       XMMF,THRESH,(I,(FIDXY(I,J),J=1,2),I=1,NFID)
        IF (ONLINE) WRITE(ITOUT,FMT=6310) G1OD,BASEOD,CURV,N1OD,
     +       0.01*THICK,NFID,XMMF,THRESH,
     $       (I,(FIDXY(I,J),J=1,2),I=1,NFID)
 6310   FORMAT(/,1X,'Film characteristics:',/,1X,'Selwyn granularity',
     +       F5.1,'   Base od',F5.2,'   Non linearity ',F6.4,/,1X,
     +       'Number of grey levels corresponding to an Od of 1.0 is',
     +       I4,/,'Film thickness (this affects expansion of ',
     +       ',measurement box)',F5.2,' mm'/,1X,
     +       'There are',I2,' fiducials. Search box size is ',F4.1,
     +       'mm and the threshold is',F4.1,' od',/,1X,'Fid. ',
     +       'coords (mm) :',3('  Fid',I2,2F6.1))
C     
      END IF
C     
      WRITE(IOUT,FMT=6315) XMMDB,ITHRESHF
      IF (ONLINE) WRITE(ITOUT,FMT=6315) XMMDB,ITHRESHF
 6315 FORMAT(1X,'Direct beam search box size (FIDUCIAL BEAM)(mm)',
     +     F4.1,/,1X,'Threshold for locating direct beam',
     +     ' (FIDUCIAL THRESHOLD',I7)
C     
      WRITE(IOUT,FMT=6320) NREC,IYLEN
      IF (ONLINE) WRITE(ITOUT,FMT=6320) NREC,IYLEN
 6320 FORMAT(1X,'Image consists of',I6,' stripes each of',I6,' pixels',
     +     ' (SIZE)')
      IF (NHEAD.EQ.0) THEN
        WRITE(IOUT,FMT=6322)
        IF (ONLINE) WRITE(ITOUT,FMT=6322)
 6322   FORMAT(1X,'There is NO header information in image file')
      ELSE                                            
        WRITE(IOUT,FMT=6324) NHEAD
        IF (ONLINE) WRITE(ITOUT,FMT=6324) NHEAD
 6324   FORMAT(1X,'Header information is contained in the first',I2,
     +       ' records of the file')
      END IF
C     
      IF (IMGP) THEN
        IF (YSCALIN.EQ.0.0) YSCALIN = YSCAL
        IF (SPIRAL) THEN
          WRITE(IOUT,FMT=6325) XTOFRA,YSCAL,ITILT,ITWIST,0.01*ROFF,
     +       0.01*TOFF,0.01*RDROFF,0.01*RDTOFF
          IF (ONLINE) WRITE(ITOUT,FMT=6325) XTOFRA,YSCAL,ITILT,ITWIST,
     +       0.01*ROFF,0.01*TOFF,0.01*RDROFF,0.01*RDTOFF
          IF ((NODES.NE.0).AND.(.NOT.FIXPAR(10)).AND.(.NOT.FIXPAR(11)))
     +       THEN
            WRITE(IOUT,FMT=6327) NODES+1
            IF (ONLINE) WRITE(ITOUT,FMT=6327) NODES+1
 6327     FORMAT(1X,'A radially dependent ROFF and TOFF will also', 
     +         ' be refined.',/,1X,'The dependence will be ',
     +         'sinusoidal, with ',I2,' nodes (zero values)',
     +         'between',/,1X,'the centre and outside of scan')     
            IF (NPHI.NE.0) THEN
              WRITE(IOUT,FMT=6329) NPHI
              IF (ONLINE) WRITE(ITOUT,FMT=6329) NPHI
 6329       FORMAT(1X,'The phase of the radially dependent distor',
     +           'tion at zero radius will',/,1X,'be',I2,'xPI/4')
            END IF
          END IF
        ELSE
C
C---- Non spiral image plates, CCDs etc
C
          WRITE(IOUT,FMT=6328) XTOFRA,YSCAL,ITILT,ITWIST
          IF (ONLINE) WRITE(ITOUT,FMT=6328) XTOFRA,YSCAL,ITILT,ITWIST
        END IF
 6325     FORMAT(/,1X,'Distortion parameters (DISTORTION):',/,1X,
     +       'Crystal to detector distance scalar (XTOFRA)',F7.4,/,
     +       1X,'Dividing factor relating pixel size in fast scan ',
     +       'direction',/,1X,'to that in slow direction ( YSCAL)',
     +       F7.4,/,1X,'Detector Tilt (TILT)',I5,'   Twist (TWIST)',
     +       I5,/,1X,'Radial offset (mm) (ROFF)',
     +       F6.2,'  Tangential offset (mm) (TOFF)',F7.2,/,1X,
     +       'Radially dependent ROFF (RDROFF)',F6.2,'mm,',
     +       'Radially dependent TOFF (RDTOFF)',F6.2,'mm')
 6328     FORMAT(/,1X,'Distortion parameters (DISTORTION):',/,1X,
     +         'Crystal to detector distance scalar (XTOFRA)',F7.4,/,
     +         1X,'Dividing factor relating pixel size in fast scan ',
     +         'direction',/,1X,'to that in slow direction ( YSCAL)',
     +         F7.4,/,1X,'Detector Tilt (TILT)',I5,'   Twist (TWIST)',
     $         I5)

      ELSE
C
C---- Film !!
C
        WRITE(IOUT,FMT=6330) XTOFRA,YSCAL,ITILT,ITWIST,IBULGE
        IF (ONLINE) WRITE(ITOUT,FMT=6330) XTOFRA,YSCAL,ITILT,ITWIST,
     +       IBULGE
 6330   FORMAT(1X,'Distortion parameters (DISTORTION):',/,1X,
     +       'Crystal to detector distance scalar (XTOFRA)',F7.4,
     +       'Dividing factor relating pixel size in fast scan ',
     +       'direction to that in slow direction ( YSCAL)',F7.4,
     +       /,1X,'Detector Tilt',I5,'   Twist',I5,'   Bulge',I5)
      END IF
C     
      PCCX = CCX
      IF (INVERTX) PCCX = -CCX
      WRITE(IOUT,FMT=6335) 0.01*PCCX,0.01*CCY,CCOM
      IF (ONLINE) WRITE(ITOUT,FMT=6335) 0.01*PCCX,0.01*CCY,CCOM
 6335 FORMAT(1X,'Camera constants (mm and degrees): CCX',F7.2,'  CCY',
     +     F7.2,'  CCOM',F8.3)
      ROTATED = ((ABS(OMEGAF/DTOR).LT.1.0).OR.((ABS(OMEGAF/DTOR-180.0))
     +     .LT.1.0))
C     
C---  Set direction of radial background  strip
C     
      IF ((IXOFFSET.EQ.0).AND.(IYOFFSET.EQ.0)) THEN
        RADX = ROTATED
        RADY = (.NOT.RADX)
      END IF
      IF (ROTATED) THEN
        WRITE(IOUT,FMT=6340)
        IF (ONLINE) WRITE(ITOUT,FMT=6340)
      ELSE
        WRITE(IOUT,FMT=6350)
        IF (ONLINE) WRITE(ITOUT,FMT=6350)
      END IF
 6340 FORMAT(/,1X,'The oscillation axis is parallel to the fast ',
     +     'axis in the image (scanner Y)')
 6350 FORMAT(/,1X,'The oscillation axis is parallel to the slow ',
     +     'axis in the image (scanner X)')
C     
C---- If online and findcc or precess is set, then
C     set filmplot true for all packs
C     
      IF (ONLINE .AND. (FINDCC.OR.PRECESS)) THEN
        DO 730 I = 1,MAXPAX
          FILMPLOT(I) = .TRUE.
 730    CONTINUE
      END IF
C     
C---- If precession photo, increase value of limit and display
C     (if set to defaults)
C     
      IF (PRECESS) THEN
        IF (LIMIT.EQ.2500) LIMIT = 4000                          
        IF (DISPLAY.EQ.25.0) DISPLAY = 40.0
      END IF
C     
C---- Print parameters associated with refinement
C     
C     
C---- Set default value of LIMIT to half max X coordinate
C     
      IF (ILIMIT.EQ.0.AND.XLIMIT.EQ.0.0) THEN
        XLIMIT = 0.5*XMAX*0.01
        LIMIT = 100*XLIMIT
      END IF
      IF (USEBOX.AND.(ISIGSET.EQ.0)) NSIG = 20
C
C---- Skip printing of refinement parameters if doing a phi profile integration
C
      IF (PHIPROF) GOTO 980
C
      WRITE(IOUT,FMT=6352) MAXREF,XLIMIT,MINREF,NSIG
      IF (ONLINE) WRITE(ITOUT,FMT=6352) MAXREF,XLIMIT,MINREF,NSIG
 6352 FORMAT(/,1X,'Parameters affecting refinement (REFINE keyword)',
     +     /,1X,'================================================',
     +     /,1X,'Initial refinement uses up to',I4,' (MAXREF) ',
     +     'reflections within',F6.1,
     +     'mm (LIMIT)',/,1X,'of the image centre. A minimum of',I4,
     +     ' reflections (NREF) are required.',/,1X,'The spots must',
     +     ' have an rms',
     +     ' variation in pixel values more than ',I2,' times ',/,1X,
     +     'that expected for a uniform background (NSIG).')
      IF (USEPAR .AND. .NOT.USEOVR) THEN
        WRITE (IOUT,FMT=6154) STRL1
 6154   FORMAT (1X,A,' Reflections will be used for positional ',
     +       'refinement (INCLUDE PARTIALS)',/,
     +       1X,'(You may wish to increase the allo',
     +       'wed maximum residual (keyword RESID)')
        IF (ONLINE) WRITE (ITOUT,FMT=6154) STRL1
        IF (ADDPART) THEN
          WRITE(IOUT,FMT=6163)
          IF (ONLINE) WRITE(ITOUT,FMT=6163)
 6163     FORMAT(1X,'Partials will be summed before using them for',
     +         ' refinement')
        ELSE
          WRITE(IOUT,FMT=6159) PTMIN
          IF (ONLINE) WRITE(ITOUT,FMT=6159) PTMIN
 6159     FORMAT(1X,'(Only partials greater than',F5.2,' recorded',
     +         ' will be used (INCLUDE PARTIALS)',F5.1,')')
        END IF
      ELSE IF (.NOT.USEPAR .AND. USEOVR) THEN
        WRITE (IOUT,FMT=6154) STRL2
        IF (ONLINE) WRITE (ITOUT,FMT=6154) STRL2
      ELSE IF (USEPAR .AND. USEOVR) THEN
        WRITE (IOUT,FMT=6154) STRL3
        IF (ONLINE) WRITE (ITOUT,FMT=6154) STRL3
        IF (SUMPART) THEN
          WRITE(IOUT,FMT=6163)
          IF (ONLINE) WRITE(ITOUT,FMT=6163)
        ELSE
          WRITE(IOUT,FMT=6159) PTMIN
          IF (ONLINE) WRITE(ITOUT,FMT=6159) PTMIN
        END IF
      ELSE
        WRITE (IOUT,FMT=6156)
 6156   FORMAT (1X,'PARTIALS and OVERLOADS will be rejected for posi',
     +       'tional refinement unless',/,1X,'there are not enough ',
     +       'fully recorded non-overloaded reflections.',/,1X,
     +       '(Use INCLUDE PARTIALS or',
     +       ' INCLUDE OVERLOADS to use these spots by default.)')
        IF (ONLINE) WRITE (ITOUT,FMT=6156)
      END IF
      IF (USEBOX) THEN
        WRITE(IOUT,FMT=6430)
        IF (ONLINE) WRITE(ITOUT,FMT=6430)
      END IF
 6430 FORMAT(1X,'The measurement box will be used in the ',
     +     'initial determination of the centre',/,1X,'of gravity',
     +     ' of spots in the central region (USEBOX)')
      IF (RWEIGHT) THEN
        WRITE(IOUT,FMT=6520) NCYC,WRMSLIM
        IF (ONLINE) WRITE(ITOUT,FMT=6520) NCYC,WRMSLIM
 6520   FORMAT(/,1X,'Following ',I3,' cycles of refinement (CYCLES) ',
     +       'if the weighted residual',/,1X,'exceeds',F5.1,
     +       ' (RESID) then processing will be abandoned.')
      ELSE
        WRITE(IOUT,FMT=6354) NCYC,0.01*RMSLIM
        IF (ONLINE) WRITE(ITOUT,FMT=6354) NCYC,0.01*RMSLIM
 6354   FORMAT(/,1X,'Following ',I3,' cycles of refinement (CYCLES) ',
     +       'if the rms residual exceeds '/,1X,F6.1,' mm ',
     +       '(RESID) then processing will be abandoned.')
      END IF
      IF (ONLINE) WRITE (ITOUT,FMT=6098) IRFMIN,IRFINC
 6098 FORMAT (1X,'For the outer regions of the image, spots with',
     +     ' I/(sigma(I)) > CUTOFF will be',/,1X,'selected ',
     +     'where CUTOFF is',I4,' for the outermost bins',
     +     ' and is incremented by',/,1X,I3,' for each bin',
     +     ' working towards centre of image (IMIN)')
      WRITE (IOUT,FMT=6098) IRFMIN,IRFINC
      IF (.NOT.IMGP) THEN     
        IF (ONLINE) WRITE (ITOUT,FMT=6058) NSDR
 6058   FORMAT (/,1X,'Intensity/SD Ratio (ISDR) for selection of ',
     +       'refinement spots for B and C',/,1X,'films  set to',I4)
        WRITE (IOUT,FMT=6058) NSDR
      END IF
      WRITE(IOUT,FMT=6356) GRADMAXR,BGFREJ
      IF (ONLINE) WRITE(ITOUT,FMT=6356) GRADMAXR,BGFREJ
 6356 FORMAT(1X,'For spots selected for refinement, the maximum ',
     +     ' allowed value for',/,1X,'(background gradient)/',
     +     '(average background is',F6.3,' (GRADIENT) and the',/,1X,
     +     'maximum allowed fraction of rejected background ',
     +     'pixels is',F5.2,' (BGREJECT)')
      IF (RWEIGHT) THEN
        WRITE(IOUT,FMT=6357)
        IF (ONLINE) WRITE(ITOUT,FMT=6357)
      END IF
 6357 FORMAT(1X,'Reflections will be weighted by the estimated',
     +     ' error in the position',/,1X,'of their centre ',
     +     'of gravity (WEIGHT)')
C     
C---- Make up list of parameters that have been explicitly FIXED for
C     the positional refinement (XCEN,YCEN,OMEGA0,YSCAL etc)
      NP = 7
C     
C     ***** machine specific code follows *****
C     
      IF (IMGP.AND.SPIRAL) NP = 11
      FIXSTR = ' '
      DO 734 I = 1,NP
        IF (FIXPAR(I)) THEN
          IF (I.EQ.8.AND..NOT.IMGP) FIXSTRA(I) = 'BULGE'
          NCH = LENSTR(FIXSTR)
          IF (NCH.EQ.0) NCH = 1
          NCH2 = LENSTR(FIXSTRA(I))
          IF (NCH.EQ.1) THEN
            FIXSTR = FIXSTRA(I)(1:NCH2)
          ELSE
            FIXSTR = FIXSTR(1:NCH)//','//FIXSTRA(I)(1:NCH2)
          END IF
          FIXEDPR = .TRUE.
        END IF
 734  CONTINUE
      IF (FIXEDPR) THEN
        WRITE(IOUT,FMT=6359) FIXSTR
        IF (ONLINE) WRITE(ITOUT,FMT=6359) FIXSTR
 6359   FORMAT(/,1X,'The following parameters will be fixed during',
     +       ' the positional refinement',/,1X,'(FIX; ',
     +       'use FREE to allow refinement of parameters',
     +       ' that are fixed by default):',/,1X,A)
      END IF
      IF (NOSMOOTH) THEN
        WRITE(IOUT,FMT=6361)
        IF (ONLINE) WRITE(ITOUT,FMT=6361)
 6361   FORMAT(/,1X,'No smoothing of detector parameters will be',
     +          ' applied (use SMOOTH to change).')
      ELSE
        WRITE(IOUT,FMT=6363) NSM1,NSM2
        IF (ONLINE) WRITE(ITOUT,FMT=6363) NSM1,NSM2
 6363   FORMAT(/,1X,'The refined detector parameters will be ',
     +     'smoothed (use NOSMOOTH to stop',/,1X,'smoothing). The',
     +     ' smoothing will be done over',I4,' degrees (NSM1)',
     +        ' for CCOMEGA,',/,1X,'YSCALE and DISTANCE and over ',
     +        I4,' (NSM2) degrees for all other parameters.')
      END IF
      IF (DAMPREF) THEN
        WRITE(IOUT,FMT=6365)
        IF (ONLINE) WRITE(ITOUT,FMT=6365)
      END IF
 6365 FORMAT(1X,'Shifts in TILT and TWIST parameters will be damped',
     +       ' according to their estimated errors.')
C
      IF (NIMGSUMREF.GT.1) THEN
        WRITE(IOUT,FMT=6364) NIMGSUMREF
        IF (ONLINE) WRITE(ITOUT,FMT=6364) NIMGSUMREF
 6364   FORMAT(/,1X,I3,' images will be summed when determining',
     +          ' reflection centroids (SUMIMAGES).')
      END IF
C     
C---- POSTREFINEMENT
C     

 980  IF ((POSTREF).AND.(.NOT.PHIPROF)) THEN
        IF (DEBUG(52)) THEN
          WRITE(IOUT,FMT=6401) NPACKS,NADD,PRCELL,GUI_SWITCH,MULTISEG,
     +          (LCELL(I),I=1,6),(FCELL(I),I=1,6),
     +          (REFINE_CELL(I),I=1,6)
          IF (ONLINE) WRITE(ITOUT,FMT=6401) NPACKS,NADD,PRCELL,
     +          GUI_SWITCH,MULTISEG,(LCELL(I),I=1,6),(FCELL(I),I=1,6),
     +          (REFINE_CELL(I),I=1,6)
 6401     FORMAT(1X,'NPACKS',I4,'NADD=',I3,' PRCELL',L2,' GUI_SWITCH',
     +      L2,' MULTISEG',L2,'  LCELL',6I3,/,1X,'FCELL:',6L2,
     +      '  REFINE_CELL:',6L2)
        END IF
        IF (.NOT.PROFILE) THEN
          PROFILE = .TRUE.
          WRITE(IOUT,FMT=6400)
          IF (ONLINE) WRITE(ITOUT,FMT=6400)
 6400     FORMAT(/,1X,'Postrefinement option requires that profile',
     +         ' fitting is used, so it has been turned on')
        END IF
C     
C---- Set up appropriate mode of post-refinement based on crystal
C     symmetry.
C     For trigonal or higher symmetry used SINGLE and refine cell and
C     missets
C     For Orthorhombic or lower use WIDTH 10 and refine cell and missets
C     .
C     
C     Do NOT set default if ADD or SINGLE or WIDTH has been given on
C     POSTREF card, or if multisegment post-refinement is being used
C     
        IF ((.NOT.PRMODE).AND.(.NOT.MULTISEG)) THEN
          IF (LCELL(2).EQ.-1) THEN
C     
C---- Orthorhombic or lower
C     
            ANGWIDTH = 10
          ELSE
C     
C---- Trigonal or higher
C     
            NADD = 1
          END IF
          IF(NEWPREF)THEN 
c     
c     should this really be NIVB = ....?
C     
            NADD = MAX(NIVB,(INT(ETA*2.0/(phirng*DTOR))
     $           +1))
          END IF
        END IF
C     
C---- Set default refined cell parameters for SINGLE case
C     If user has not explicitly fixed/unfixed parameters then
C     if trigonal or higher unfix all refineable params
C     
        IF ((NADD.EQ.1).AND.(.NOT.PRCELL).AND.
     +       (LCELL(2).NE.-1)) THEN
          DO 683 I = 1,6
            IF (LCELL(I).EQ.-1) UNFIX(I) = .TRUE.
 683      CONTINUE  
        END IF
C     
C---- Reset cell refinement flags if in single image mode, default
C     is to fix all cell parameters unless explicitly UNFIXED
C     23sep05 Change this, default for single image mode now
C     the same as if many images used.
CAL        IF (NADD.EQ.1) THEN
CAL          DO 682 I = 1,6
CAL            FCELL(I) = .TRUE.
CAL            IF (UNFIX(I)) FCELL(I) = .FALSE.
CAL 682      CONTINUE  
CAL        END IF 

c     Things to control the setting of fcell in the instance when
c     we're using integrate_setup - aha, but what if we're using 
c     refine_setup? this can be a pain! -> set refine_cell and lcell?
C     "refine_cell" below is never assigned a value. Therefore this code
C     is not useful and has been commented out. agwl 10/10/12
C
CAL        if(gui_switch) then
CAL           do i = 1, 6
CAL              fcell(i) = .true.
CAL              if(lcell(i) .eq. -1) then
c     we're free to refine cell(i)
CAL                 if(refine_cell(i)) then
CAL                    fcell(i) = .false.
CAL                 end if
CAL              end if
CAL           end do
CAL        end if

c     Things to control the refinement of the detector parameters when
c     we're using the new GUI interface - oops - maybe this will also
c     affect the postrefinement -> will have to set refine_detector
c     to something useful too.
CAL   This code does not do anything useful, comment it all out now. 
CAL   AGWL 10/10/12
CAL        if(gui_switch) then
CAL           do i = 1, 11
CAL Comment this out until we find from Geoff what is supposed to happen here.
CAL              fixpar(i) = .not. refine_detector(i)
CAL           end do
CAL        end if
           
C     
C---- Now make list of cell parameters that are to be refined.
C     
C     LCELL(I) = -1  parameter free
C     =  0  parameter fixed
C     .gt. 0, = J  parameter I constrained to = parameter J
C     
        REFCELL = .FALSE.
        CELLSTR = ' '
        DO 732 I = 1,6
          IF ((.NOT.FCELL(I)).AND.(LCELL(I).LT.0)) THEN
            IF (REFCELL) THEN
              NCH = LENSTR(CELLSTR)
              NCH2 = LENSTR(SABC(I))
              CELLSTR = CELLSTR(1:NCH)//','//SABC(I)(1:NCH2)
            ELSE
              CELLSTR = SABC(I)
            END IF
            REFCELL = .TRUE.
          END IF
 732    CONTINUE
        IF (DEBUG(52)) THEN
          WRITE(IOUT,FMT=6403) NADD,PRCELL,GUI_SWITCH,MULTISEG,
     +          (LCELL(I),I=1,6),(FCELL(I),I=1,6),
     +          (REFINE_CELL(I),I=1,6),REFCELL
          IF (ONLINE) WRITE(ITOUT,FMT=6403) NADD,PRCELL,GUI_SWITCH,
     +          MULTISEG,(LCELL(I),I=1,6),(FCELL(I),I=1,6),
     +          (REFINE_CELL(I),I=1,6),REFCELL
 6403     FORMAT(1X,'NADD=',I3,' PRCELL',L2,' GUI_SWITCH',L2,
     +      ' MULTISEG',L2,'  LCELL',6I3,/,1X,'FCELL:',6L2,
     +      '  REFINE_CELL:',6L2,'  REFCELL',L2)
        END IF
C     
C---- Set up maximum allowed residual from input mosaic spread, beam
C     divergence (in generate file). Need to convert eta etc from
C     half-widths in radians.
C     
        RESIDMAX = RSDMAX*2.0*(ETA + 0.5*(DIVH+DIVV))/DTOR 
C     
C---- For post refinement using several images, if angular WIDTH was
C     specified rather than an explicit number of images (ADD) then
C     convert the angular width to number of images here.
C     
        IF (NADD.EQ.0) THEN
          NADD = NINT(ANGWIDTH/PHIRNG)
C     
C---- If all cell parameters have been fixed, then make the default NADD
C     =1
C     so that the mosaic spread is refined after every image rather than
C     after 10 degrees
C     
          IF (.NOT.REFCELL) NADD = 1
          IF (NADD.GT.NIMAX) THEN
            WRITE(IOUT,FMT=6432) NADD,NIMAX
            IF (ONLINE) WRITE(ITOUT,FMT=6432) NADD,NIMAX
            IF (BRIEF) WRITE(IBRIEF,FMT=6432) NADD,NIMAX
 6432       FORMAT(1X,'**** FATAL ERROR ****',/,1X,'You have ',
     +           'asked for the post-refinement to be done over',
     +           I3,/,1X,'images but this exceeds the maximum a',
     +           'llowed (',I3,/,1X,'Either reduce WIDTH or cha',
     +           'nge parameter NIMAX and recompile')
            IF(SOCKLO)THEN
              LINE400 = ' '
              WRITE(LINE400,FMT=7100)NADD,NIMAX
              CALL WRITE_SOCKET_LENGTH(SERVERFD,
     $             LENSTR(LINE400),LINE400)
 7100         FORMAT('<?xml version="1.0"?>',
     $             '<!DOCTYPE cell_refine_response>',
     $             '<cell_refine_response>',
     $             '<status><code>error</code><message>',
     $             'Requested number of images for postrefinement is ',
     $             I4,', but maximum allowed is ',I4,
     $             '</message></status></cell_refine_response>')
              GOTO 50
              
            ELSE
              STOP
            ENDIF
          END IF
        END IF
C     
        IF (NADD.EQ.1) THEN
          IF (NPACKS.GT.1) THEN
            WRITE(IOUT,FMT=6410)
            IF (ONLINE) WRITE(ITOUT,FMT=6410)
 6410       FORMAT(/,1X,'POST REFINEMENT'/,1X,'===============',
     +           /,1X,'Post refinement will be used to refine ',
     +           'the missetting angles after each',/,1X,'imag',
     +           'e (POSTREF, use POSTREF OFF to prevent post ',
     $           'refinement')
C     

            IF (MULTISEG) THEN
              WRITE(IOUT,FMT=6415) NADD,NSEG,
     +             NEWMATNAM(1:LENSTR(NEWMATNAM))
              IF (ONLINE) WRITE(ITOUT,FMT=6415) NADD,NSEG,
     +             NEWMATNAM(1:LENSTR(NEWMATNAM))
              IF (NADD.GT.NIMAX) THEN
                WRITE(IOUT,FMT=6097) NADD,NIMAX
                IF (ONLINE) WRITE(ITOUT,FMT=6097) NADD,NIMAX
                STOP   
              END IF
            END IF
C     
            WRITE(IOUT,FMT=6409) NPRMIN
            IF (ONLINE) WRITE(ITOUT,FMT=6409) NPRMIN
            IF (REFCELL) THEN
              WRITE(IOUT,FMT=6412) CELLSTR
              IF (ONLINE) WRITE(ITOUT,FMT=6412) CELLSTR
              IF (LCELL(2).EQ.-1) THEN
                WRITE(IOUT,FMT=6240)
                IF (ONLINE) WRITE(ITOUT,FMT=6240)
 6240           FORMAT(1X,' *** BEWARE ***',
     +               /,1X,' CELL REFINEMENT MAY BE UNSTABLE ',
     +               'FOR CRYSTAL SYMMETRY LOWER THAN TRIGONAL')
              END IF
 6412         FORMAT(1X,'In addition, the following cell paramet',
     +             'ers will be refined:',/,1X,A)
            ELSE
              WRITE(IOUT,FMT=6419)
              IF (ONLINE) WRITE(ITOUT,FMT=6419)
 6419         FORMAT(1X,'All cell parameters will be fixed')
            END IF
          ELSE
C     
C---- Not enough images for refinement
C     
            WRITE(IOUT,FMT=6421)
            IF (ONLINE) WRITE(ITOUT,FMT=6421)
 6421       FORMAT(1X,'POST REFINEMENT (POSTREF)'/,1X,
     +           '===============',/,1X,
     +           'No post-refinement will be done as only one',
     +           ' image is being processed')
          END IF
        ELSE
C     
C---- NADD > 1
C     
          IF ((NPACKS.LT.(NADD+1)).AND.(.NOT.MULTISEG).AND.REFCELL)
     +         THEN
            WRITE(IOUT,FMT=6429) NPACKS,NADD
            IF (ONLINE) WRITE(ITOUT,FMT=6429) NPACKS,NADD
 6429       FORMAT(1X,'POST REFINEMENT (POSTREF)'/,1X,
     +           '==============='
     +           ,/,1X,'**** WARNING ****'
     +           ,/,1X,'**** WARNING ****'
     +           ,/,1X,'**** WARNING ****'
     +           ,/,1X,'**** WARNING ****'
     +           ,/,1X,'**** WARNING ****'
     +           ,/,1X,'Cell parameters will ',
     +           ' NOT be refined because only',I3,' images are to',
     +           ' be',/,1X,'processed (PROCESS keyword) but',I3,
     +           ' images are required for cell parameter',/,1X,
     +           ' refinement (NADD or WIDTH subkeywords on POSTRE',
     +           'F keyword).',/,1X,'Because post-refinement uses ',
     +           'partially recorded reflections at the end of',/,
     $           1X,'one image and the start of the next, it is ne',
     +           'cessary to process ONE MORE',/,1X,
     $           'image than the number to be used in the',
     +           ' refinement',//,1X,'The crystal orientation will',
     +           ' be refined for every image',/)
          ELSE
            IF (REFCELL) THEN
              WRITE(IOUT,FMT=6411)
              IF (ONLINE) WRITE(ITOUT,FMT=6411)
            ELSE
              WRITE(IOUT,FMT=7390)
              IF (ONLINE) WRITE(ITOUT,FMT=7390)
            END IF
 6411       FORMAT(/,1X,'POST REFINEMENT (POSTREF)'/,1X,
     +           '===============',/,1X,
     +           'Post refinement will be used to refine cell ',
     +           'parameters and missetting angles.',/,
     +           1X,'Use POSTREF OFF to prevent post refinement.')
 7390       FORMAT(/,1X,'POST REFINEMENT (POSTREF)'/,1X,
     +           '===============',/,1X,
     +           'Post refinement will be used to refine ',
     +           'missetting angles.',/,
     +           1X,'(Use POSTREF OFF to prevent post refinement.)')
            WRITE(IOUT,FMT=6409) NPRMIN
            IF (ONLINE) WRITE(ITOUT,FMT=6409) NPRMIN
 6409       FORMAT(1X,'Refinement will only be carried out if th',
     +           'ere are more than',I5,' reflections',/,1X,'sele',
     +           'cted for the refinement (NREF).')
            IF (REFCELL) THEN
              WRITE(IOUT,FMT=6417) CELLSTR
              IF (ONLINE) WRITE(ITOUT,FMT=6417) CELLSTR
 6417         FORMAT(1X,'The following cell parameters',
     +             ' will be refined ',/,1X,'(Use FIX keyword to',
     +             ' fix individual parameters or FIX ALL to fix',
     +             /,1X,'all of them):',/,1X,A)
            ELSE
              WRITE(IOUT,FMT=6419)
              IF (ONLINE) WRITE(ITOUT,FMT=6419)
            END IF
C     
C---- Set REFCELL true as this is used in writing summary information
C     This is NOT correct. REFCELL is set correctly above and should not
C     be reset here. AGWL 25/07//07  
CAL            REFCELL = .TRUE.
C
            IF (MULTISEG) THEN
              WRITE(IOUT,FMT=6415) NADD,NSEG,
     +             NEWMATNAM(1:LENSTR(NEWMATNAM))
              IF (ONLINE) WRITE(ITOUT,FMT=6415) NADD,NSEG,
     +             NEWMATNAM(1:LENSTR(NEWMATNAM))
              IF (NADD.GT.NIMAX) THEN
                WRITE(IOUT,FMT=6097) NADD,NIMAX
                IF (ONLINE) WRITE(ITOUT,FMT=6097) NADD,NIMAX
            IF(SOCKLO)THEN
              LINE400 = ' '
              WRITE(LINE400,FMT=7100)NADD,NIMAX
              CALL WRITE_SOCKET_LENGTH(SERVERFD,
     $             LENSTR(LINE400),LINE400)
              GOTO 50
            ELSE
              STOP
            ENDIF
 6097           FORMAT(1X,'**** FATAL ERROR ****',/,1X,'You have ',
     +               'asked for the post-refinement using a total',
     +               ' of',I3,/,1X,'images but this exceeds the m',
     +               'aximum allowed (',I3,/,1X,'Either reduce WI',
     +               'DTH or change parameter NIMAX and recompile')
              END IF
 6415         FORMAT(1X,'Data from ',I3,' (ADD) images in',I3,
     +             ' (SEGMENT) different segments will be combined',
     +             /,1X,'for use in post-refinement. Each segment ',
     +             'must be specified on a separate PROCESS',/,1X,
     +             'keyword followed by a RUN keyword.',/,1X,'The ',
     +             'images will NOT be integrated',/,1X,
     +             'The final orientation matrix and cell will',
     +             ' be written to file:'/,1X,A,/)
            ELSE
              IF (REFCELL) THEN
                WRITE(IOUT,FMT=6413) NADD
                IF (ONLINE) WRITE(ITOUT,FMT=6413) NADD
 6413         FORMAT(1X,'Data from the',I3,' (ADD or WIDTH) previ',
     +             'ous images will be ',
     +             'added together',/,1X,'and used to refine the ',
     +             'current parameters.'/,1X,'However, the ',
     +             'missetting angles will be refined after every',
     +             ' image'/,1X,'until enough images have been ',
     +             'processed to allow cell parameter refinement')
              ELSE
                WRITE(IOUT,FMT=6431)
                IF (ONLINE) WRITE(ITOUT,FMT=6431)
 6431           FORMAT(1X,'Missetting angles will be refined.')
                IF (SMOOTHPHI) THEN
                  WRITE(IOUT,FMT=6433) NSMOOTH
                  IF (ONLINE) WRITE(ITOUT,FMT=6433) NSMOOTH
 6433             FORMAT(1X,'Refined missetting angles will be ',
     +                      'smoothed over',I3,' images.')
                END IF
              END IF
            END IF
          END IF
        END IF
        IF (NPACKS.GT.NADD) THEN
          IF (PRNS.EQ.0) THEN
            WRITE(IOUT,FMT=6414)
            IF (ONLINE) WRITE(ITOUT,FMT=6414)
          ELSE IF (PRNS.EQ.1) THEN
            IF (USEBEAM) THEN
              WRITE(IOUT,FMT=6416) NSMOOTH,ETAFRAC
              IF (ONLINE) WRITE(ITOUT,FMT=6416) NSMOOTH,ETAFRAC
            ELSE
              WRITE(IOUT,FMT=6422) NSMOOTH,ETAFRAC
              IF (ONLINE) WRITE(ITOUT,FMT=6422) NSMOOTH,ETAFRAC
            END IF
          ELSE IF (PRNS.EQ.2) THEN
C     
C---- If anisotropic divergence is to be refined the initial values must
C     be non-zero. If there are zero, transfer mosaic spread to
C     diveregences
C     and if still zero, stop.
C     
            IF ((DIVH+DIVV).EQ.0.0) THEN
              IF (ETA.GT.0) THEN
                DIVH = ETA
                DIVV = DIVH
                ETA = 0
                ETAD = 0
                DIVHD = 2.0*DIVH/DTOR
                DIVVD = 2.0*DIVV/DTOR
                WRITE(IOUT,FMT=7250) DIVHD
                IF (ONLINE) WRITE(ITOUT,FMT=7250) DIVHD
              ELSE
                DIVHD = 0.1
                DIVH = 0.5*DIVHD*DTOR
                DIVVD = DIVHD
                DIVH = DIVV
                WRITE(IOUT,FMT=7252) DIVHD
                IF (ONLINE) WRITE(ITOUT,FMT=7252) DIVHD
              END IF
            END IF
            WRITE(IOUT,FMT=6418)
            IF (ONLINE) WRITE(ITOUT,FMT=6418)
          END IF
          IF (USEBEAM.AND.(PRNS.GT.0)) THEN
            WRITE(IOUT,FMT=6420)
            IF (ONLINE) WRITE(ITOUT,FMT=6420)
          END IF
          IF (PRNS.NE.0) THEN
            WRITE(IOUT,FMT=6423)
            IF (ONLINE) WRITE(ITOUT,FMT=6423)
          END IF
 6423     FORMAT(1X,'Use POSTREF BEAM 0 to turn off beam parameter',
     +         ' refinement')
 6414     FORMAT(1X,'No beam parameters will be refined (BEAM).')
 6416     FORMAT(1X,'The mosaic spread will be refined.',/,1X,
     +           'The refined value will be smoothed over',I3,
     +           ' images (MOSSMOOTH).',/,1X,'A "safety" factor ',
     +           'equal to',F5.2,' times the current value will be',
     +            ' added (MOSADD).',/,1X,'To stop smoothing, set',
     +               ' MOSSMOOTH to one.')
 6422     FORMAT(1X,'The mosaic spread will be refined BUT the input',
     +            ' value will actually be used.',/,1X,
     +           'The refined value will be smoothed over',I3,
     +           ' images (MOSSMOOTH).',/,1X,'A "safety" factor ',
     +           'equal to',F5.2,' times the current value will be',
     +            ' added (MOSADD).',/,1X,'To stop smoothing, set',
     +               ' MOSSMOOTH to one.')
 6418     FORMAT(1X,'Horizontal and vertical beam divergences ',
     +         'will be refined (BEAM).')
 6420     FORMAT(1X,'The refined mosaic spread will be used ',
     +         ' during integration.',/,1X,
     +        '(Use POSTREF FIX MOSAIC to use the input mosaic',
     +           ' spread instead of the refined one.)')
        END IF
        WRITE(IOUT,FMT=6424) SDFAC,RESIDMAX
        IF (ONLINE) WRITE(ITOUT,FMT=6424) SDFAC,RESIDMAX
 6424   FORMAT(1X,'Reflections with I .GT.',F4.1,' sigma will be',
     +       ' used in refinement (SDFAC)',/,
     +       1X,'If the refinement residual is greater than',F5.2,
     +       ' processing will be abandoned.',/,1X,'(Controlled by',
     +       ' MAXRESID, limit is MAXRESID*(EPS+MEAN(DIVH,DIVV))).')
        IF (NADD.EQ.1) THEN
          WRITE(IOUT,FMT=6425) SHIFTMAX,SHIFTFAC
          IF (ONLINE) WRITE(ITOUT,FMT=6425) SHIFTMAX,SHIFTFAC
 6425     FORMAT(1X,'If the rms change in missetting angles is ',
     +         'greater than',F5.2,' degrees',/,1X,'(MAXSHIFT), or',
     +         ' if the change in cell parameters is more than',
     +         /,1X,F5.1,' times the estimated standard deviation',
     +         ' (SHIFTFAC)',/,1X,'the ',
     +         'image will be reprocessed with an updated ',
     +         'reflection list.')
        ELSE IF (.NOT.MULTISEG) THEN
          WRITE(IOUT,FMT=6427) SHIFTMAX,SHIFTFAC,NRPT
          IF (ONLINE) WRITE(ITOUT,FMT=6427) SHIFTMAX,SHIFTFAC,NRPT
 6427     FORMAT(1X,'If the rms change in ',
     +         'missetting angles is greater than',F5.2,
     +         ' degrees',/,1X,'(MAXSHIFT) then that image will be ',
     +         'reprocessed',/,1X,'If, on the first refinement of ',
     +         'cell parameters, any cell parameter changes',
     +         ' by',/,1X,'more than',F5.1,' times its estimated sd',
     +         ' (SHIFTFAC) then ',
     +         'the processing will be',/,1X,'restarted using ',
     +         'the updated parameters',/,1X,'This will be done a',
     +         ' maximum of ',I2,' (REPEAT) times')
        END IF
        WRITE(IOUT,FMT=6249) PARTITION, 1.0-PARTITION
        IF (ONLINE) WRITE(ITOUT,FMT=6249) PARTITION,1.0-PARTITION
 6249   FORMAT(1X,'Partials will be partitioned into two parts ',
     +     'with fractional',/,1X,'intensities',F5.2,' and',F5.2,
     +     ' (PARTITION)')
        
        IF (PREFEXICE) THEN
          WRITE(IOUT,FMT=6248) ICEWIDTHPREF
          IF (ONLINE) WRITE(ITOUT,FMT=6248) ICEWIDTHPREF
 6248     FORMAT(/,1X,'Reflections that lie within shells centred on',
     +       ' ice rings will be excluded.',/,1X,
     +       'The half-width of the shells is set to ',
     +         F6.4,' times the actual resolution (ICE).')
        END IF
C     
C---- End of IF (POSTREF) block
C     
      END IF
C     
C     
      IF ((PROFILE).AND.(.NOT.PHIPROF)) THEN
C     
C---- Check non-compatible keywords
C     
        IF (PRBFILM .OR. PRCFILM .AND. (.NOT.PRREAD)) THEN
          IF (.NOT.ACCUMULATE) THEN
            IF (ONLINE) WRITE (ITOUT,FMT=6172)
 6172       FORMAT (//,1X,'**** WARNING ****',/,1X,'if using the ',
     +           'B or C films to form the standard profiles, then',
     $           /,1X,'ACCUMULATE must  be turned on. this has be',
     $           'en done.')
            WRITE (IOUT,FMT=6172)
            ACCUMULATE = .TRUE.
          END IF
C     
          IF (PRCFILM) PRBFILM = .TRUE.
          DO 740 I = 1,MAXPAX
            FORCEB(I) = .TRUE.
            IF (PRCFILM) FORCEC(I) = .TRUE.
 740      CONTINUE
        END IF
C     
        IF (PRREAD) THEN 
          ACCUMULATE = .FALSE.
          FIRSTPASS = .FALSE.
          PRBFILM = .FALSE.
          PRCFILM = .FALSE.
        END IF                   
C     
        IF (MULTISEG) GOTO 745
C     
        WRITE (IOUT,FMT=6174) 
 6174   FORMAT (//,1X,'PROFILE FITTING'/,1X,'===============',/,1X,
     +       '(Use PROFILE OFF to suppress profile fitting)')
C     
C---- If partials have been included in refinement and ADDPART is
C     not being used, then include them in profile formation UNLESS
C     PROFILE FULLS has been specified.
C     
        IF ((USEPAR).AND.(.NOT.ADDPART).AND.(.NOT.PRFULLS).AND.
     +       (.NOT.PRPART)) THEN
          PRPART = .TRUE.
          WRITE(IOUT,FMT=6175)
 6175     FORMAT(/,/,1X,'***** WARNING *****',/,1X,
     +         '***** WARNING *****',/,1X,
     +         '***** WARNING *****',/,1X,
     +         '***** WARNING *****',/,1X,
     +         '***** WARNING *****',/,1X,
     +         'Because you have requested inclusion ',
     +         'of partials in refinement (REFINEMENT ',/,1X,
     +         'INCLUDE PARTIALS) partials will ALSO be used ',
     +         'in forming the standard profiles.',/,1X,'This',
     +         ' is equivalent to including keywords: ',
     +         'PROFILE PARTIALS',/,1X,
     +         'If you do NOT want to include partials',
     +         ' in forming profiles, use keywords: ',/,1X,
     +         'PROFILE FULLS')
        END IF
C     
        IF (PRREAD) THEN
          WRITE (IOUT,FMT=6190) PROFFNR(1:NCH3)
 6190     FORMAT (1X,'PROFILES for this run will be read from logi',
     +         'cal file PROFILE (READ)'/,1X,'Filename: ',A)
        ELSE
          IF (PRSET) THEN
            IF (LINESET) THEN
              WRITE(IOUT,FMT=6185) (0.01*XLINE(I),I=1,NXLINE)
              WRITE(IOUT,FMT=6187) (0.01*YLINE(I),I=1,NYLINE)
 6185         FORMAT(1X,'The coordinates of the lines ',
     +             '(in mm) defining the standard areas are:',/,1X,
     +             'In X direction',9(1X,F6.1))
 6187         FORMAT(1X,'In Y direction',9(1X,F6.1))
            ELSE IF (HIGHRES) THEN
              WRITE(IOUT,FMT=6177)
            ELSE IF (LOWRES) THEN
              WRITE(IOUT,FMT=6179)
            END IF
          ELSE
            IF (WAVE/DSTMAX.LT.2.5) THEN
              WRITE(IOUT,FMT=6181)
            ELSE
              WRITE(IOUT,FMT=6191)
            END IF
          END IF
 6177     FORMAT(1X,'High resolution profile binning (21',
     +         ' standard profiles) has been selected')
 6179     FORMAT(1X,'Low resolution profile binning (9',
     +         ' standard profiles) has been selected')
 6181     FORMAT(1X,'The default profile binning for this resolution',
     +         ' giving 25 standard profiles',/,1X,'will be ',
     +         'used (Default is 25 bins above 2.5A, 9 below.)',
     +         /,1X,'To change this use XLINES/YLINES to define',
     +         ' the areas you want')
 6191     FORMAT(1X,'The default profile binning for this resolution',
     +         ' giving 9 standard profiles',/,1X,'will be ',
     +         'used (Default is 25 bins above 2.5A, 9 below)',
     +         /,1X,'To change this use XLINES/YLINES to define',
     +         ' the areas you want')
C     
          IF (PROPT) THEN
            WRITE(IOUT,FMT=6171) TOLMIN,TOL,IBOUND
          ELSE IF (PROPTCEN) THEN
            WRITE(IOUT,FMT=6601) TOLMIN,TOL,IBOUND
          END IF
          
 6171     FORMAT(1X,'The measurement box parameters will be optim',
     +         'ised for every standard',/,1X,'profile independen',
     +         'tly (OPTIMISE). The',
     +         ' parameters used in the optimisation are:',/,1X,
     +         'Toleance (TOLERANCE) Minimum ',F5.3,' Maximum ',
     +         F5.3,'  Boundary (BOUNDARY)',I3,
     +         /,1X,'See the help library for a description of th',
     +         'ese parameters',/,1X,'(To switch off optimisation',
     +         ' for all standard profiles use NOOPT, to turn it ',
     +         'off',/,1X,'for the average central spot profile a',
     +         'lso use NOOPT ATALL)')
 6601     FORMAT(1X,'The measurement box parameters will be optim',
     +         'ised for the central region',/,1X,'only (NOOPTIMI',
     +         'SE).  The',
     +         ' parameters used in the optimisation are:',/,1X,
     +         'Tolerance (TOLERANCE) Minimum ',F5.3,' Maximum ',
     +         F5.3,'  Boundary (BOUNDARY)',I3,
     +         /,1X,'See the help library for a description of th',
     +         'ese parameters',/,1X,'(To switch off optimisation',
     +         ' completely use NOOPTIMISE ATALL)')
          IF (FIXBOX) THEN
            WRITE(IOUT,FMT=6600)
          ELSE                  
            WRITE(IOUT,FMT=6602) BGPKRAT,FRACREJ
          END IF
 6600     FORMAT(1X,'The overall dimensions of the box will not',
     +         ' be altered in the optimisation (FIXBOX)')
 6602     FORMAT(1X,'The overall dimensions will be set to achi',
     +         'eve a minimum ratio of background',/,1X,'to pea',
     +         'k pixels of', F5.2,' (RATIO)',/,1X,'The expansi',
     +         'on in either direction to achieve this will be',/,
     +         1X,'halted when more than',F5.2,' of the pixels ',
     +         'are rejected (STOP)',/,1X,'(Use FIXBOX to switc',
     +         'h off optimisation of OVERALL dimensions of box)')
          WRITE (IOUT,FMT=6176) ISDRATIO,PRBGSIG,NOVPIX,PRCUTOFF,
     +         NRFMIN,RMSBGPR                
 6176     FORMAT (1X,'Reflections with peak pixel intensity less th',
     +         'an',I3,' times (ISDR) the rms',/,1X,'variation in b',
     +         'ackground density will be rejected from the STANDAR',
     +         'D PROFILES',/,1X,'For the standard profiles, backgr',
     +         'ound points which lie more than',/,1X,F5.1,' sigma ',
     +         '(BGSIG) from the best least squares plane will be r',
     +         'ejected',/,1X,'Spots with more than',I3,' (NOVPIX u',
     +         'nder OVERLOADS) pixel values greater',/,1X,'than',I12,
     $         ' (CUTOFF',
     +         ') will be treated as OVERLOADS and excluded.',
     +         /,1X,'The rejection criteria applied to each profile',
     +         ' are:',/,1X,'Minimum number of reflections con',
     +         'tributing to profile',I4,' (NREF)',/,1X,
     +         'Maximum rms variation in background (after scaling',
     +         ' peak to 255)',F5.1,' (RMSBG)')
          IF (PUPDATE) WRITE(IOUT,FMT=6173)
 6173     FORMAT(/,1X,'For the first block of images, the profiles ',
     +         'will be redetermined',/,1X,'after the raster box pa',
     +         'rameters have been optimised.',/,1X,'For succeeding',
     +         ' blocks, the optimised parameters of the previous b',
     +         'lock',/,1X,'will be used (use PRUPDATE/NOPRUPDATE).')
          IF (INTERPOL) WRITE(IOUT,FMT=6122)
 6122     FORMAT (/,1X,'Pixel values will be interpolated to place ',
     +         'CALCULATED spot position',/,1X,'in centre of measur',
     +         'ement box (INTERPOLATE).')
          IF (CHANGEMASK) WRITE (IOUT,FMT=6178)
 6178     FORMAT (/,1X,'The background mask for each measurement bo',
     +         'x will be updated on the basis of',/,1X,'background',
     +         ' points rejected from the standard profile for that',
     +         ' box.',/,1X,'The updated mask will be used in the ',
     +         ' integration pass (CHANGEMASK).')
          IF (ACCUMULATE) THEN
            IF (PRCFILM) THEN
              WRITE (IOUT,FMT=6184)
 6184         FORMAT (/,1X,'Profiles will be accumulated from the',
     $             ' A,B and C films of all the packs specified',/,
     $             1X,'prior to measuring all the films in each',
     $             ' pack',/,1X,'in order to achieve this',
     +             ', B and C films will be treated as A films ',
     $             'for positional refinement',/,2X,'that is, st',
     $             'arting with reflections fr',
     +             'om the central region',/,1X,
     +             'and then using the whole film')
            ELSE IF (PRBFILM) THEN
              WRITE (IOUT,FMT=6182)
 6182         FORMAT (/,1X,'Profiles will be accumulated from the',
     $             ' A and B films of all the packs specified',/,
     $             1X,'prior to measuring first all the A+B films',
     $             ' and then the subsequent films in',/,1X,
     +             'each pack',/,1X,'in order to achieve this,',
     $             ' B films will be treated as A films for ',
     $             'positional refinement',/,2X,
     +             'that is, starting with reflections from ',
     $             'the central region'
     +             ,/,1X,'and then using the whole film')
            ELSE
              IF (IMGP) THEN
                WRITE (IOUT,FMT=6183) NBLOCK
              ELSE
                WRITE (IOUT,FMT=6180)
              END IF
 6180         FORMAT (/,1X,'Profiles will be accumulated from the ',
     $             'A films of all the packs specified',/,1X,
     $             'Prior to measuring first all the A films ',
     $             'and then the subsequent films in',/,1X,'each',
     +             ' pack (ACCUMULATE).')
 6183         FORMAT (/,1X,'Profiles will be accumulated over',I3,
     +             ' images (BLOCK subkeyword on PROCESS line)',/,
     +             1X,'prior to integration (ACCUMULATE).')
            END IF              
          END IF
          IF (LPRINT(11)) WRITE (IOUT,FMT=6186)
 6186     FORMAT (/,1X,'All profiles will be printed (PRINT).')
          IF (PRSAVE) WRITE (IOUT,FMT=6188) PROFFNW(1:NCH4)
 6188     FORMAT (/,1X,'The PROFILES will be written to ',A,
     +         '_00n.prf',/,1X,' where n is 1,2,3 etc for ',
     +         'each succeeding block of data (SAVE).')
          IF (WEIGHT) WRITE (IOUT,FMT=6192)
 6192     FORMAT (/,1X,'The least-squares fit of the standard profil',
     +         'es to individual spot pixel',/,1X,'values will be',
     +         ' weighted (WEIGHT).')
          IF (PKONLY) WRITE(IOUT,FMT=6502)
 6502     FORMAT(1X,'The profile will be fitted to the peak',
     +         ' pixels only, rather than fitting',/,1X,
     +         'a plane plus the scaled profile to the entire box.')
          IF (WTPROFILE.AND.(.NOT.PRPART)) WRITE(IOUT,FMT=6199)
 6199     FORMAT(1X,'In the formation of the standard profiles, each',
     +         ' reflection will',/,1X,'be weighted using counting',
     +         ' statistics (WSUM).')
          IF (VARPRO) WRITE (IOUT,FMT=6193)
 6193     FORMAT (/,1X,'A separate profile will be calculated for ',
     +         'every reflection as a',/,1X,'weighted sum of ',
     +         'neighbouring profiles (VARIABLE).')
          IF (PRPART) THEN
            WRITE (IOUT,FMT=6194)
            WTPROFILE = .FALSE.
          END IF
 6194     FORMAT (/,1X,'** PARTIALS will be included in forming the ',
     +         'standard PROFILES **',/,1X,
     +         'In this case it is not valid',
     +         ' to use weighting in the formation of',/,1X,'the ',
     +         'standard profiles so this has been suppressed')
          IF (DISCRIMINATE) WRITE(IOUT,FMT=6201) DISCRIM
 6201     FORMAT(/,1X,'Spots for which the highest pixel value (af',
     +         'ter background subtraction)',/,1X,'is less than',
     +         F6.1,' times greater than the largest pixel value',
     +         ' in',/,1X,'the background region will be eliminate',
     +         'd from the standard profiles.')
          WRITE(IOUT,FMT=6203) PKWDLIM1,PKWDLIM2
 6203     FORMAT(/,1X,'For fully recorded reflections, individual ',
     +         'peak pixels will be rejected',/,1X,'from profile f',
     +         'itting if their fit to the scaled profile',
     +         ' deviates by more',/,1X,'than',F7.1,' times the ex',
     +         'pected error (WDLIM1).',/,1X,'For ALL reflections,',
     +         ' for pixels adjacent to overlapped pixels, the',/,
     +         1X,'rejection factor is',F7.1,' (WDLIM2).')
          IF (PKWDOUTL.NE.0) WRITE(IOUT,FMT=6205) IOUTL1,IOUTL2,
     +         PKWDOUTL
 6205     FORMAT(1X,'In addition, peak pixels with values between ',
     +         I6,' and',I6,' which deviate',/,1X,'by more than',
     +         F7.1,' times the expected error will be rejected.')
          IF (PREXICE) WRITE(IOUT,FMT=6206) ICEWIDTH
 6206     FORMAT(/,1X,'When forming the profiles, spots that lie ',
     +         'within shells centred on ice',/,1X,'rings will be ',
     +         'excluded. The half-width of the shells is set to ',
     +         F6.4,/,1X,'times the actual resolution (ICE).')
          WRITE(IOUT,FMT=6208) PRCUT
 6208     FORMAT(/,1X,'When forming the profiles, a fraction',F6.3,
     +         ' of the reflections that have',/,1X,'the highest',
     +         ' pixel values in each profile will be excluded, ',
     +         'to eliminate',
     +         /,1X,'the effect of zingers and hot pixels (PRCUT).')
        END IF
C     
        WRITE(IOUT,FMT=6210) NOVPIX,CUTOFF
 6210   FORMAT(//,1X,'REFLECTION INTEGRATION'/,1X,'================',
     +       '======'/,1X,'OVERLOADS (OVERLOAD):',/,1X,
     +       'Any reflection with',
     +       ' more than',I3,' (NOVER) pixels ',
     +       'with a',/,1X,'value greater than ',I12,' (CUTOFF)',
     +       ' will be flagged as overloads.') 
        IF (USEOVRLD) WRITE(IOUT,FMT=6195)
 6195   FORMAT (/,1X,'The intensity of overloaded reflections will ',
     +       'be estimated by profile fitting.')
C     
        WRITE(IOUT,FMT=6220) BGFRAC,BGSIG
 6220   FORMAT(/,1X,'Background evaluation (BACKGROUND):',/,1X,
     +       'A fraction',F5.2,' (BGFRAC)',
     +       ' of the background pixels',/,1X,'will be used in the ',
     +       'initial determination of the plane constants',/,1X,
     +       'Pixels deviating by more than',F5.1,' (BGSIG) sigma ',
     +       'from this initial',/,1X,'plane will be rejected.')
C     
        IF (RECOVER) WRITE(IOUT,FMT=6630) NINT(RECLEVEL*NBGMIN),
     +       NINT(RECLEVEL*NBGMIN)
 6630   FORMAT(/,1X,'If there are fewer than',I3,' background ',
     +       'pixels remaining after rejecting',/,1X,'those overlapp',
     +       'ed by neighbouring spots, background pixels with',/,1X,
     +       'the lowest values in the standard profile will be ',
     +       'included until',I3,' background',/,1X,
     +       'pixels are obtained. This number is RECLEVEL*NBGMIN ',
     +       'where RECLEVEL is set by',/,1X,'keywords BACKGROUND ',
     +       'RECOVER RECLEVEL and NBGMIN by REJECTION',
     +       ' MINBG NBGMIN')
C     
        WRITE(IOUT,FMT=6631) GRADMAX,NBGMIN,BGRAT,PKRAT
 6631   FORMAT(/,1X,'REJECTION CRITERIA (REJECTION):',/,1X,
     +       'Reflections for which the (background gradient)/',
     +       '(average background) is',/,1X,'greater than',F6.3,
     +       ' will be rejected',
     +       ' (GRADMAX)',/,1X,'Minimum number of background ',
     +       'pixels remaining after rejection of'/,1X,'outliers',I4,
     +       ' (MINB)',/,1X,
     +       'Maximum BGRATIO ',F5.1,' (BGRATIO)',/,1X,'Maximum ',
     +       'PKRATIO ',F5.1,' (PKRATIO)')
        IF (PKACCEPT) THEN
          WRITE(IOUT,FMT=6223)
        ELSE
          WRITE(IOUT,FMT=6221)
        END IF
        WRITE(IOUT,FMT=6225)
 6221   FORMAT(1X,'Reflections failing the PKRATIO test will ',
     +       'be rejected (use REJECT PKRATIO ACCEPT',/,1X,
     +       ' to keep the summation integration value.)')
 6223   FORMAT(1X,'For reflections which fail the PKRATIO test,',
     +       ' the profile fitted intensity',/,1X,'and sd will be ',
     +       'replaced by the summation integration intensity and',
     +       ' sd.')
 6225   FORMAT(/,1X,'Rejected reflections will be written to the',
     +       ' MTZ file, but will be',/,1X,'rejected by AIMLESS/SCALA',
     +       ' by default. Use the ACCEPT keywords in',/,1X,
     +       ' AIMLESS/SCALA to include these measurements.')
C     
C---- Test if PLOT option requested without specifying SCANNER (leaves
C     EFAC set to -999). Note BADPLOT is also set true by DUMPSPOT so
C     must exclude this possibility.
C     
        IF (BADPLOT.AND.(.NOT.DUMPSPOT)) THEN
          IF (EFAC.LT.-900) THEN
            WRITE(IOUT,FMT=6222)
          ELSE
            WRITE(IOUT,FMT=6224)
          END IF
        END IF
 6222   FORMAT(/,1X,'**** WARNING ****',/,1X,'Pixel values for ',
     +       '"badspots" have been requested (PLOT)',/,1X,'but the ',
     +       'scanner error has not been assigned (SCANNER keyword)',
     +       /,1X,'To get pixel values the scanner error MUST be ',
     +       'assigned a value')
 6224   FORMAT(1X,'Pixel values of "badspots" will be printed')
C     
C     
        IF (USEDGE) WRITE(IOUT,FMT=6197)
 6197   FORMAT (/,1X,'The intensity of reflections with up to half ',
     +       'the pixels outside the',/,1X,'scanned area will be ',
     +       'estimated by profile fitting')
        IF (ADDPART) WRITE(IOUT,FMT=6360)          
 6360   FORMAT(/,1X,'PARTIALS:',/,1X,'Partials spanning two image',
     +       's will be summed to give the fully recorded',/,1X,
     +       'reflection (ADDPART... Use ADDPART OFF to suppress)')
        IF (IANGLE.EQ.1) WRITE(IOUT,FMT=6362) NWMAX
 6362   FORMAT(/,1X,'Partial reflections spanning more than',I3,
     +       ' images will NOT be integrated',/,1X,
     $       '(use MAXWIDTH to change)')
C     
C---- Now the same ouput to channel ITOUT
C     
        IF (ONLINE) THEN
          WRITE (ITOUT,FMT=6174)
          IF ((USEPAR).AND.(.NOT.ADDPART).AND.(.NOT.PRFULLS).AND.
     +         (.NOT.PRPART)) WRITE(ITOUT,FMT=6175)
          IF (PRREAD) THEN
            WRITE (ITOUT,FMT=6190)  PROFFNR(1:NCH3)       
          ELSE
            IF (PRSET) THEN
              IF (LINESET) THEN
                WRITE(ITOUT,FMT=6185) (0.01*XLINE(I),I=1,NXLINE)
                WRITE(ITOUT,FMT=6187) (0.01*YLINE(I),I=1,NYLINE)
              ELSE IF (HIGHRES) THEN
                WRITE(ITOUT,FMT=6177)
              ELSE IF (LOWRES) THEN
                WRITE(ITOUT,FMT=6179)
              END IF
            ELSE
              IF (WAVE/DSTMAX.LT.2.5) THEN
                WRITE(ITOUT,FMT=6181)
              ELSE
                WRITE(ITOUT,FMT=6191)
              END IF
            END IF
            IF (PROPT) THEN
              WRITE(ITOUT,FMT=6171) TOLMIN,TOL,IBOUND
            ELSE IF (PROPTCEN) THEN
              WRITE(ITOUT,FMT=6601) TOLMIN,TOL,IBOUND
            END IF
            IF (FIXBOX) THEN
              WRITE(ITOUT,FMT=6600)
            ELSE                  
              WRITE(ITOUT,FMT=6602) BGPKRAT,FRACREJ
            END IF
            WRITE (ITOUT,FMT=6176) ISDRATIO,PRBGSIG,NOVPIX,PRCUTOFF,
     +           NRFMIN,RMSBGPR
            IF (PUPDATE) WRITE(ITOUT,FMT=6173)
            IF (INTERPOL) WRITE(ITOUT,FMT=6122)  
            IF (CHANGEMASK) WRITE (ITOUT,FMT=6178)
            IF (ACCUMULATE) THEN
              IF (PRCFILM) THEN
                WRITE (ITOUT,FMT=6184)
              ELSE IF (PRBFILM) THEN
                WRITE (ITOUT,FMT=6182)
              ELSE
                IF (IMGP) THEN
                  WRITE (ITOUT,FMT=6183)
                ELSE
                  WRITE (ITOUT,FMT=6180)
                END IF
              END IF
            END IF
            IF (LPRINT(11)) WRITE (ITOUT,FMT=6186)
            IF (PRSAVE) WRITE (ITOUT,FMT=6188) PROFFNW(1:NCH4)
            IF (WEIGHT) WRITE (ITOUT,FMT=6192)
            IF (PKONLY) WRITE(ITOUT,FMT=6502)
            IF (WTPROFILE) WRITE(ITOUT,FMT=6199)
            IF (VARPRO) WRITE (ITOUT,FMT=6193)
            IF (PRPART) WRITE (ITOUT,FMT=6194)
            IF (DISCRIMINATE) WRITE(ITOUT,FMT=6201) DISCRIM
            WRITE(ITOUT,FMT=6203) PKWDLIM1,PKWDLIM2
            IF (PKWDOUTL.NE.0) WRITE(ITOUT,FMT=6205) IOUTL1,IOUTL2,
     +           PKWDOUTL
            IF (PREXICE) WRITE(ITOUT,FMT=6206) ICEWIDTH
            WRITE(ITOUT,FMT=6208) PRCUT
          END IF
C     
          WRITE(ITOUT,FMT=6210) NOVPIX,CUTOFF
          IF (USEOVRLD) WRITE(ITOUT,FMT=6195)
          WRITE(ITOUT,FMT=6220) BGFRAC,BGSIG
          IF (RECOVER) WRITE(ITOUT,FMT=6630) NINT(RECLEVEL*NBGMIN),
     +         NINT(RECLEVEL*NBGMIN)
          WRITE(ITOUT,FMT=6631) GRADMAX,NBGMIN,BGRAT,PKRAT
          IF (PKACCEPT) THEN
            WRITE(ITOUT,FMT=6223)
          ELSE
            WRITE(ITOUT,FMT=6221)
          END IF
          WRITE(ITOUT,FMT=6225)
          IF (USEDGE) WRITE(ITOUT,FMT=6197)
          IF (ADDPART) WRITE(ITOUT,6360)          
          IF (IANGLE.EQ.1) WRITE(ITOUT,FMT=6362) NWMAX
        END IF
C     
C---- End of IF (PROFILE) THEN block
C     
      END IF
C     
C---- Other options
C     
      IF (OTHERS) THEN       

        WRITE(IOUT,FMT=6500)
 6500   FORMAT(//,1X,'Other options',/,1X,'=============')
        
        IF (AVPR) THEN
          WRITE (IOUT,FMT=6014)
 6014     FORMAT (1X,
     $         'Average spot profile will be printed for all packs')
          IF (ONLINE) WRITE (ITOUT,FMT=6014)
        END IF
        IF (PRINTL) THEN
          IF (IMGP) THEN
            I = 2
          ELSE
            I = 1
          END IF
          WRITE (IOUT,FMT=6086) (PRINTOP(K)(1:NCHPR(K)),
     +         K=I,4)
 6086     FORMAT (1X,'Additional printout for the following:',/,1X,4A)
          IF (ONLINE) WRITE (ITOUT,FMT=6086) 
     +         (PRINTOP(K)(1:NCHPR(K)),K=I,4)
        END IF
        IF (MATCH.AND.RMOSAIC) THEN
          TEMP = MIN(RESOL1,RESOL2)
          IF (TEMP.EQ.0.0) TEMP = MAX(RESOL1,RESOL2)
          WRITE(IOUT,FMT=6133) TEMP,ETAMAX
          IF (ONLINE) WRITE(ITOUT,FMT=6133) TEMP,ETAMAX
 6133     FORMAT(/,1X,'Beam divergence refinement',/,1X,
     +         '==========================',/,1X,
     +         'The beam divergence will be refined using pattern',
     +         ' matching and data to',F5.1,'A'/,1X,
     +         'The divergence will be varied between',
     +         ' 0.0 and',F4.1,' degrees')
          IF (NBEAM.EQ.1) THEN
            WRITE(IOUT,FMT=6135)
            IF (ONLINE) WRITE(ITOUT,FMT=6135)
          ELSE
C     
C---- If anisotropic divergence is to be refined the initial values must
C     be non-zero. If there are zero, transfer mosaic spread to
C     diveregences
C     and if still zero, stop.
C     
            IF ((DIVH+DIVV).EQ.0.0) THEN
              IF (ETA.GT.0) THEN
                DIVH = ETA
                DIVV = DIVH
                ETA = 0
                ETAD = 0
                DIVHD = 2.0*DIVH/DTOR
                DIVVD = 2.0*DIVV/DTOR
                WRITE(IOUT,FMT=7250) DIVHD
                IF (ONLINE) WRITE(ITOUT,FMT=7250) DIVHD
 7250           FORMAT(1X,'If refining vertical and horizontal ',
     +               'beam divergences independently',/,1X,'the',
     +               'ir initial values must be non-zero.',/,1X,
     $               'Therefore horizontal and vertical divergen',
     +               'ce have been set to ',F5.2,' degrees and',/,
     +               1X,'the mosaic spread set to zero.')
              ELSE
                DIVHD = 0.1
                DIVH = 0.5*DIVHD*DTOR
                DIVVD = DIVHD
                DIVH = DIVV
                WRITE(IOUT,FMT=7252) DIVHD
                IF (ONLINE) WRITE(ITOUT,FMT=7252) DIVHD
 7252           FORMAT(1X,'If refining vertical and horizontal be',
     +               'am divergences independently',/,1X,'their i',
     +               'nitial values must be non-zero.',/,1X,'Ther',
     +               'efore horizontal and vertical divergence ha',
     $               've been set to ',F5.2,' degrees')
              END IF
            END IF
            WRITE(IOUT,FMT=6137)
            IF (ONLINE) WRITE(ITOUT,FMT=6137)
          END IF
 6135     FORMAT(1X,'An isotropic divergence will be refined ',
     +         '(this is the sum of mosaic',/,1X,'spread and ',
     $         'beam divergence).')
 6137     FORMAT(1X,'Horizontal and vertical divergence will be',
     +         ' refined separately (refined parameters',/,1X,
     +         'are the sum of mosaic spread and beam divergence',
     $         ').')
          IF (NOREFINE) THEN
            WRITE(IOUT,FMT=6145)
            IF (ONLINE) WRITE(ITOUT,FMT=6145)
 6145       FORMAT(1X,
     $           'No orientation refinement will be carried out.')
          ELSE
            WRITE(IOUT,FMT=6147)
            IF (ONLINE) WRITE(ITOUT,FMT=6147)
 6147       FORMAT(1X,'Mosaic spread refinement will be carried ou',
     +           't following orientation refinement')
          END IF
          WRITE(IOUT,FMT=6149)
          IF (ONLINE) WRITE(ITOUT,FMT=6149)
 6149     FORMAT(1X,'Intensities will NOT be measured')
        END IF
        IF (MATCH.AND.(.NOT.NOREFINE)) THEN
          WRITE (IOUT,FMT=6132) RCONV,OVRLAP,NSTEP,SECANGLE,NPASS,
     +         DAMP,TRUECCOM
 6132     FORMAT(/,1X,'Orientation refinement',/,1X,'============',
     +         '==========',/,1X,'Automatic pattern matching will ',
     +         'be performed (AUTO)',/,1X,'The radius of convergen',
     +         'ce has been set to ',F5.1,' degrees (RCONV)',/,1X,
     +         'Minimum overlap of calculated and observed',
     +         ' patterns is',F5.1,' Degrees (OVERLAP)',
     +         /,1X,'Number of steps used in the matching',I3,
     +         ' (NSTEP)'/,1X,'The acceptance semi-angle is',F6.1,
     +         ' Degrees (ANGLE)',
     +         /,1X,I2,' Passes (NPASS) will be made at each reso',
     +         'lution with a step damping factor',/,1X,'of',F5.2,
     +         ' (DAMP)'/,1X,'The camera constant ccomega will be ',
     +         'set to ',F6.2,' degrees (CCOMEGA) after the',/,1X,
     +         'central refinement and the missetting angle PSIX ',
     +         'adjusted accordingly',/,1X,
     +         '(This should be the nominally correct value of',
     +         ' this camera constant)')
          IF (ONLINE) WRITE (ITOUT,FMT=6132) RCONV,OVRLAP,NSTEP,
     +         SECANGLE,NPASS,DAMP,TRUECCOM
          IF (RESOL2.EQ.0) THEN
            IF (ONLINE) WRITE (ITOUT,FMT=6134) RESOL1
 6134       FORMAT (1X,'Resolution for pattern matching is ',F5.1,
     +           ' Angstroms (RESOL)')
            WRITE (IOUT,FMT=6134) RESOL1
          ELSE
            IF (ONLINE) WRITE (ITOUT,FMT=6136) RESOL1,RESOL2
 6136       FORMAT (1X,'Two passes will be made using data to ',
     $           F5.1,'Angstr',
     +           'oms in the first pass',/,1X,'and',F5.1,
     +           ' Angstroms in the second (RESOL)')
            WRITE (IOUT,FMT=6136) RESOL1,RESOL2
          END IF
          IF (NOCENT) THEN
            IF (ONLINE) WRITE (ITOUT,FMT=6138)
 6138       FORMAT (/,1X,'CENTRS Refinement will be skipped ',
     $           '(NOCENT)')
            WRITE (IOUT,FMT=6138)
          ELSE
            IF (RWEIGHT) THEN
              TEMP = 2.0
              WRITE(IOUT,FMT=6530) NCYCA,TEMP,AWRMSLIM
              IF (ONLINE) WRITE(ITOUT,FMT=6530) NCYCA,TEMP,
     +             AWRMSLIM
 6530         FORMAT(/,1X,'The initial refinement of the central ',
     +             'region will consist of',I2,' cycles (NCYC)',/,
     $             1X,'The refinement ',
     +             'will be repeated if the initial weighted ',
     +             'residual',/,1X,' exceeds',F4.1,
     +             /,1X,' and if the final weighted residual',
     +             ' exceeds',F4.1,/,1X,' (RESID) then ',
     +             'processing will be abandoned')
            ELSE
              WRITE(IOUT,FMT=6139) NCYCA,0.01*AELIMIT,0.01*ARMSLIM
              IF (ONLINE) WRITE(ITOUT,FMT=6139) NCYCA,0.01*AELIMIT,
     +             0.01*ARMSLIM
 6139         FORMAT(/,1X,'The initial refinement of the central',
     +             ' region will consist of',I2,' cycles (NCYC)',
     $             /,1X,'The refinement will be repeated ',
     +             'if the initial residual exceeds',F6.3,
     +             /,1X,'mm (ELIMIT) and if the final residual',
     +             ' exceeds',F6.3,' mm ',/,1X,' (RESID) then ',
     +             'processing will be abandoned')
            END IF
          END IF
          IF (NOMEAS) THEN
            IF (ONLINE) WRITE (ITOUT,FMT=6140)
 6140       FORMAT (/,1X,'*** The film will not be measured after ',
     +           'orientation refinement (NOMEAS)***')
            WRITE (IOUT,FMT=6140)
          END IF
        END IF
C     
C---- End of IF (OTHERS) THEN  block
C     
      END IF
C     
C     
C---- Check incompatible combinations of keywords
C     
C---- If doing AUTO, can ONLY measure one pack per generate file,
C     unless ALSO doing POSTREF, in which case AUTO will be run for
C     the first pack ONLY, although all will be POSTREF'd.
C     
      IF (MATCH) THEN
        IF (NPRUN.GT.1) THEN
          IF (POSTREF) THEN
            WRITE(IOUT,FMT=6143)
            IF (ONLINE) WRITE(ITOUT,FMT=6143)
 6143       FORMAT(/,1X,'Orientation refinement using AUTO will ',
     +           'ONLY be performed on the first image',/,1X,'Po',
     +           'st refinement will be carried out on all images')
          ELSE
            WRITE(IOUT,FMT=6141)
            IF (ONLINE) WRITE(ITOUT,FMT=6141)
 6141       FORMAT(//,1X,'***** fatal error *****',/,1X,'If AUTO',
     +           ' orientation refinement is to be performed, th',
     +           'en only ONE pack can',/,1X,'be processed UNLES',
     +           'S post refinement (POSTREF) is also specified.',
     $           /,1X,'If post refinement is requested',
     +           ' then only the FIRST image will be refined',
     +           /,1X,'using pattern matching')
            STOP
          END IF
        END IF
      END IF
C     
C---- Set default EFAC if not doing profile fitting to zero.
C     
      IF ((.NOT.PROFILE) .AND. (EFAC.LT.-900.0)) EFAC = 0.0
C     
      IF ((PROCES) .AND. (.NOT.PROFILE)) THEN
        IF (ONLINE) THEN
          WRITE (ITOUT,FMT=6168)
 6168     FORMAT (/,1X,'**** PROCESS Option has been requested bu',
     +         't not PROFILE ***',/,1X,'Please give PROFILE keyw',
     +         'ord or end this run')                     
          GO TO 50
        ELSE
          WRITE (IOUT,FMT=6170)
 6170     FORMAT (/,1X,'*** PROCESS Keyword has been given which ',
     +         'requires PROFILE fitting ***',/,1X,'PROFILE has t',
     +         'herefore been turned on')
          PROFILE = .TRUE.
        END IF
      END IF
C     
C     
      NCH = LENSTR(GENFILE)
C     AL        FWORK = ' '
C     AL        CALL UGTENV('HKLOUT',FWORK)
C---- Only reset MTZNAM if no HKLOUT keyword given and environment
C     variable
C     HKLOUT has been set
C     
C     AL        IF ((FWORK(1:1).NE.' ').AND.(IHKLOUT.EQ.0)) MTZNAM =
C     FWORK
C     
      IF (NPRUN.GT.1) THEN
        WRITE (IOUT,FMT=6164) IPACK1A(1),IPACK2A(NSERRUN),
     +       GENFILE(1:NCH),MTZNAM(1:LENSTR(MTZNAM))
 6164   FORMAT (//,1X,
     $       '*************************************************',
     +       '******************************',
     +       /,1X,'Processing images',I7,' to',I7,/,1X,'The output ',
     +       ' generate file is ',A,/,1X,'The output MTZ file is ',A)
        IF (ONLINE) WRITE(ITOUT,FMT=6164) IPACK1A(1),
     +       IPACK2A(NSERRUN),GENFILE(1:NCH),MTZNAM(1:LENSTR(MTZNAM))
        IF (NSERRUN.EQ.1) THEN
          IF (ISERADD.GT.0) THEN
            WRITE(IOUT,FMT=6165) ISERADD
            IF (ONLINE) WRITE(ITOUT,FMT=6165) ISERADD
 6165       FORMAT(1X,I6,' will be added to all image numbers to ',
     +           'generate the batch numbers',/,1X,'in the MTZ file')
          END IF
        ELSE
          DO 940 I = 1,NSERRUN
            WRITE(IOUT,FMT=7530) ISERAR(I),I
            IF (ONLINE) WRITE(ITOUT,FMT=7530) ISERAR(I),I
 940      CONTINUE
        END IF
 7530   FORMAT(1X,I6,' will be added to image numbers in part',I2,
     +       ' to generate the batch numbers',/,1X,
     +       'in the MTZ file')
        WRITE(IOUT,FMT=6167)
        IF (ONLINE) WRITE(ITOUT,FMT=6167)
 6167   FORMAT(1X,
     +       '***************************************************',  
     +       '****************************'//)
      ELSE
        WRITE (IOUT,FMT=6166) IPACK1A(1),GENFILE(1:NCH),
     +       MTZNAM(1:LENSTR(MTZNAM))
 6166   FORMAT(//,1X,
     $       '**************************************************',
     +       /,1X,'Processing image',I7,/,1X,'The output generate ',
     +       'file is ',A,/,1X,'The output MTZ file is ',A)  
        IF (ONLINE) WRITE (ITOUT,FMT=6166) IPACK1A(1),GENFILE(1:NCH),
     +       MTZNAM(1:LENSTR(MTZNAM))
        IF (ISERADD.GT.0) THEN
          WRITE(IOUT,FMT=6165) ISERAR(1)
          IF (ONLINE) WRITE(ITOUT,FMT=6165) ISERAR(1)
        END IF
        WRITE(IOUT,FMT=6167)
        IF (ONLINE) WRITE(ITOUT,FMT=6167)
      END IF
C     
C     
C---- Set up pass number if using accumulated profiles
C     
 745  FIRSTPASS = (ACCUMULATE .AND. PROFILE)
C     
C     
      DONERUN = .TRUE.
C     
C---- Convert spot separations (IXSEP,IYSEP) into "ideal detector"
C     coordinate
C     frame, as the spot coordinates (generate file coords) are in this
C     frame
C     
      MINDTX = NINT(ABS(IXSEP*COSOM0 + IYSEP*SINOM0))
      MINDTY = NINT(ABS(IYSEP*COSOM0 + IXSEP*SINOM0))
      RETURN
C     
C---- Eof on input
C     Check that a RUN card has been included
C     
 750  IF ((.NOT.DONERUN).AND.(.NOT.STRATEGY)) THEN
        WRITE (IOUT,FMT=6200)
 6200   FORMAT (//,1X,'*** WARNING **   NO "RUN" CARD GIVEN')
        IF (ONLINE) WRITE (ITOUT,FMT=6200)
      END IF
      CALL CTIMER(ENDTIME)
c      WRITE(IOUT,FMT=7800)ENDTIME-STARTTIME
c      IF(ONLINE)WRITE(ITOUT,FMT=7800)ENDTIME-STARTTIME
      WRITE (IOUT,FMT=6196)
      IF (ONLINE) WRITE (ITOUT,FMT=6196)
C     
C     ******************
c     Again, forgetting to set genopen
      IF (GENOPEN) then
C
C---- Print for program flow
C
        CALLED = 'QCLOSE'
        COMMENT = 'Close generate file when EOF on input'
        IF (DEBUG(99)) CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
        CALL QCLOSE(IUNIT)
        genopen = .false.
      end if
      IF (MTZOPEN) THEN
        MTZPRT = 1
C     *********************
        if (debug(52))then
          write(iout,fmt=7754) NRUN
          if (online) write(itout,fmt=7754) NRUN
 7754     format(1x,'**closing MTZ file in control 8 ***',
     +         ' NRUN=',I3)
        end if
        CALLED = 'LWCLOS'
        COMMENT = 'To close MTZ file '
        IF (DEBUG(99)) 
     +     CALL PRINTFLOW(THIS,CALLED,COMMENT,INDENT)
        CALL LWCLOS(MTZOUT,MTZPRT)
        MTZOPEN = .FALSE.
        IF(LBEST.AND.BESTOPEN)THEN
          CLOSE(BESTHKL)
          BESTOPEN = .FALSE.
        ENDIF
c---  add call mharvest here
c     when no end card given and EOF reached in input stream
c     
C     
C        IF(CCP4VERSION.LT.'4.0')HARVESTREADY = .FALSE.
C        IF(HARVESTREADY)
        IF(DOHARVEST)CALL MHARVEST(2)
C     *********************
      END IF
C     
c     socket    call close_socket(serverfd)


c     this could be bad

      call dna_end_output
      if(dnaout) close(dnafd)
      STOP
C
C---- temporary home for FORMAT statements which have been transferred
C     to s/r SECTION but are still referenced here.
 6650 FORMAT(76('*'),/,6(' FATAL ERROR'),/,6X,
     $     'You must change all occurrences',
     $     ' of the lines ',/,/,12X,'PARAMETER (IXWDTH=8192)',
     $     /,6X,'and',/,12X,'PARAMETER (IYLENGTH=4096)',/,6X,
     $     'to',/,12X,'PARAMETER (IXWDTH=12288)',/,
     $     6X,'and',/,12X,'PARAMETER (IYLENGTH=6144)',/,/,6X,
     $     'and rebuild the program in order to process ',
     $     'Quantum 315 unbinned images.',/,6X,'These lines ',
     $     'occur in the files:',/,/,12X,
     $     'mosflm_all_ip_inc.for',
     $     /,12X,'control.f',/,12X,'celref.f',/,8X,
     $     'and reek.f.',/,/,6X,'You may also need to increa',
     $     'se the swap space on your machine',/,6X,
     $     '(it must be more than ~150Mb)',/,/,76('*'))
 6713 FORMAT(4(1X,'***** WARNING *****',/),
     +     1X,'Input crystal to detector distance (',F7.2,
     +     'mm) does NOT agree with value in the image ',
     +     'header (',F7.2,'mm)',/,1X,'The input distance',
     +     ' will be used.')
 6715 FORMAT(4(1X,'***** WARNING *****',/),
     +     1X,'Input wavelength (',F6.4,
     +     'A) does NOT agree with value in the image ',
     +     'header (',F6.4,'A)',/,1X,'The input wavelength',
     +     ' will be used.')
 6717 FORMAT(/,1X,'Crystal to detector distance of',F8.2,
     +     'mm taken from image header')
 6719 FORMAT(/,1X,'Wavelength of',F8.5,'A taken',
     +     ' from image header')
 6725 FORMAT(4(/,1X,'***** WARNING *****',/),1X,
     +     'Oscillation angle derived from image ',
     +     'header is zero. If this is an oscillation',
     $     /,1X,' image, then the header information',
     $     ' is not correct. The phi values must be',
     $     /,1X,'given on the IMAGE',
     +     ' (or PROCESS) keyword.',/,/)
 6726 FORMAT(/,1X,'Pixel size of ',F6.4,'mm taken',
     +     ' from image header.')
 6727 FORMAT(/,1X,'Saturation value set to',I12,' as set',
     +      ' in image header.')
 6728 FORMAT(1X,'Nullpix value from header:',I5)
 6729 FORMAT(4(1X,'***** WARNING *****',/),
     +     1X,'Input pixel size (',F6.4,
     +     'mm) does NOT agree with value in the image ',
     +     'header (',F6.4,'mm)',/,1X,'The input value',
     +     ' will be used.')
 6730 FORMAT(1X,'Profile cutoff also set to',I12)
      
C
      END
