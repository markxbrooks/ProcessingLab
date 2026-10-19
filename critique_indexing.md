
 ***** CORRECT ***** (VERSION January 10, 2014)  12-Sep-2014

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=    3 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    62.046   186.898    72.370  90.000 100.183  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=   1048500     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    68.950     0.000
 NAME_TEMPLATE_OF_DATA_FRAMES=/data/results/LLS_Williams/MBP4-Rac1/76c1b_0002212/76c1b_0002212_1_0???.cbf CBF
 DATA_RANGE=       1     900
 ROTATION_AXIS=  0.999945  0.000085 -0.010492
 OSCILLATION_RANGE=  0.200000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       1
 X-RAY_WAVELENGTH=  0.979490
 INCIDENT_BEAM_DIRECTION=  0.009664  0.003608  1.020887
 FRACTION_OF_POLARIZATION=   0.990
 POLARIZATION_PLANE_NORMAL=  0.000000  1.000000  0.000000
 AIR=   0.00100
 DETECTOR=PILATUS
 SILICON=  3.706798 SENSOR_THICKNESS=  0.320000
 NUMBER OF DETECTOR SEGMENTS   1
 NX=  2463  NY=  2527    QX=  0.172000  QY=  0.172000
 ORGX=   1222.98  ORGY=   1201.88
 DETECTOR_DISTANCE=   384.986
 DIRECTION_OF_DETECTOR_X-AXIS=   1.00000   0.00000   0.00000
 DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000   1.00000   0.00000
 BEAM_DIVERGENCE_E.S.D.=   0.046
 REFLECTING_RANGE_E.S.D.=   0.427
 MINIMUM_ZETA=     0.050
 MAXIMUM_ERROR_OF_SPOT_POSITION=   3.0
 MAXIMUM_ERROR_OF_SPINDLE_POSITION=   2.0
 MINIMUM_I/SIGMA=   3.0
 REFLECTIONS/CORRECTION_FACTOR=      50
 STRICT_ABSORPTION_CORRECTION=FALSE 
 CORRECTIONS= DECAY MODULATION ABSORPTION                      
 REFERENCE_DATA_SET=                                                  


  324517 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       1 ...     900
       0 OVERLOADED REFLECTIONS (IGNORED)
     225 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
  324292 REFLECTIONS ACCEPTED



 ******************************************************************************
                      AUTOMATIC SPACE GROUP ASSIGNMENT
 ******************************************************************************

 XDS adopts the following approach.
 (1) it looks for possible symmetries of the crystal lattice
 (2) it computes a redundancy independent R-factor for all enantiomorphous
     point groups compatible with the observed lattice symmetry.
 (3) it selects the group which explains the intensity data at an acceptable,
     redundancy-independent R-factor (Rmeas, Rrim) using a minimum number of
     unique reflections.

 This approach does not test for the presence of screw axes. Consequently,
 orthorhombic cell axes will be specified in increasing length (following
 conventions), despite the possibility that different assignments for the
 cell axes could become necessary for space groups P222(1) and P2(1)2(1)2
 containing one or two screw axes, respectively.

 The user can always override the automatic decisions by specifying the
 correct space group number and unit cell constants in XDS.INP and repeating
 the CORRECT step of XDS. This provides a simple way to rename orthorhombic
 cell constants if screw axes are present.

 In addition, the user has the option to specify in XDS.INP
 (a) a reference data set or
 (b) a reindexing transformation or
 (c) the three basis vectors (if known from processing a previous data set
 taken at the same crystal orientation in a multi-wavelength experiment).
 These features of XDS are useful for resolving the issue of alternative
 settings of polar or rhombohedral cells (like P4, P6, R3).

 =====> Specifications provided by the user in XDS.INP:
 SPACE_GROUP_NUMBER=    3 (0: unknown symmetry and unit cell)
 UNIT_CELL_CONSTANTS=    61.80   185.60    72.00  90.000 100.200  90.000
 UNIT_CELL_A-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 UNIT_CELL_B-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 UNIT_CELL_C-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 REIDX=  0  0  0  0  0  0  0  0  0  0  0  0 (all 0 : not specified)
 TEST_RESOLUTION_RANGE=   10.00    5.00 (Angstrom) for space group determination
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas            21016
 NUMBER OF ACCEPTED UNIQUE REFLECTIONS FROM REFERENCE DATA SET       0


 *********** DETERMINATION OF LATTICE CHARACTER AND BRAVAIS LATTICE ***********

 The CHARACTER OF A LATTICE is defined by the metrical parameters of its
 reduced cell as described in the INTERNATIONAL TABLES FOR CRYSTALLOGRAPHY
 Volume A, p. 746 (KLUWER ACADEMIC PUBLISHERS, DORDRECHT/BOSTON/LONDON, 1989).
 Note that more than one lattice character may have the same BRAVAIS LATTICE.

 A lattice character is marked "*" to indicate a lattice consistent with the
 observed locations of the diffraction spots. These marked lattices must have
 low values for the QUALITY OF FIT and their implicated UNIT CELL CONSTANTS
 should not violate the ideal values by more than
 MAX_CELL_AXIS_ERROR=  0.03
 MAX_CELL_ANGLE_ERROR=           3.0 (Degrees)

 The REINDEXING TRANSFORMATION REIDX() consists of 12 integers that relate
 the original indices H,K,L from file INTEGRATE.HKL to the indices H',K',L'
 with respect to the new cell.
    H'=(REIDX(1)*H+REIDX( 2)*K+REIDX( 3)*L)/IDXV+REIDX( 4)
    K'=(REIDX(5)*H+REIDX( 6)*K+REIDX( 7)*L)/IDXV+REIDX( 8)
    L'=(REIDX(9)*H+REIDX(10)*K+REIDX(11)*L)/IDXV+REIDX(12)
 The value of the integer IDXV depends on the lattice type used for
 specifying reflections on file INTEGRATE.HKL; IDXV is
 1 for a primitive,
 2 for a face or body centred,
 3 for a rhombohedral,
 4 for a lattice centred on all faces.
 IDXV is set by XDS and cannot be input by the user.

  LATTICE-  BRAVAIS-   QUALITY  UNIT CELL CONSTANTS (ANGSTROEM & DEGREES)    REINDEXING TRANSFORMATION
 CHARACTER  LATTICE     OF FIT      a      b      c   alpha  beta gamma

 *  44        aP          0.0      62.1   72.4  186.7  90.0  90.0 100.2   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  31        aP          0.0      62.1   72.4  186.7  90.0  90.0  79.8    1  0  0  0  0  0 -1  0  0  1  0  0
 *  34        mP          0.0      62.1  186.7   72.4  90.0 100.2  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    35        mP        103.7      72.4   62.1  186.7  90.0  90.0 100.2    0  0  1  0  1  0  0  0  0  1  0  0
    33        mP        103.7      62.1   72.4  186.7  90.0  90.0 100.2   -1  0  0  0  0  0 -1  0  0 -1  0  0
    32        oP        103.7      62.1   72.4  186.7  90.0  90.0 100.2   -1  0  0  0  0  0 -1  0  0 -1  0  0
    39        mC        146.3     147.0   62.1  186.7  90.0  90.0  75.7    1  0  2  0  1  0  0  0  0  1  0  0
    38        oC        146.3      62.1  147.0  186.7  90.0  90.0 104.3    1  0  0  0 -1  0 -2  0  0  1  0  0
    29        mC        146.3      62.1  147.0  186.7  90.0  90.0  75.7    1  0  0  0  1  0  2  0  0 -1  0  0
    14        mC        180.6      86.6  103.4  186.7  90.0  90.0  81.1   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    13        oC        180.6      86.6  103.4  186.7  90.0  90.0  81.1   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    10        mC        180.6     103.4   86.6  186.7  90.0  90.0  98.9    1  0 -1  0  1  0  1  0  0 -1  0  0
    11        tP        284.2      62.1   72.4  186.7  90.0  90.0 100.2   -1  0  0  0  0  0 -1  0  0 -1  0  0
    28        mC        301.8      62.1  378.5   72.4  91.7 100.2  80.6   -1  0  0  0 -1  2  0  0  0  0 -1  0
    12        hP        326.9      62.1   72.4  186.7  90.0  90.0 100.2   -1  0  0  0  0  0 -1  0  0 -1  0  0
    37        mC        353.7     378.5   62.1   72.4 100.2  91.7  80.6   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    36        oC        353.7      62.1  378.5   72.4  88.3 100.2  99.4   -1  0  0  0  1  2  0  0  0  0 -1  0
    30        mC        392.1      72.4  380.3   62.1  91.9 100.2  79.0    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    41        mC        443.9     380.3   72.4   62.1 100.2  91.9  79.0    0  2  1  0  0  0  1  0  1  0  0  0
    40        oC        443.9      72.4  380.3   62.1  88.1 100.2 101.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    27        mC        500.0     147.0   62.1  200.2  86.3 109.2  75.7   -1  0 -2  0 -1  0  0  0  0  1  1  0
    26        oF        625.0      62.1  147.0  378.5  87.7  99.4 104.3    1  0  0  0 -1  0 -2  0 -1  2  0  0
    42        oI        693.9      62.1   72.4  383.3  99.2  97.4 100.2    1  0  0  0  0  0  1  0 -1 -2 -1  0
    15        tI        784.2      62.1   72.4  383.3  80.8  82.6 100.2   -1  0  0  0  0  0 -1  0 -1 -2 -1  0
     9        hR        930.6      62.1   86.6  569.5  91.6  97.6 124.6    1  0  0  0 -1  0 -1  0 -1  3  1  0
     1        cF        999.0     205.8  213.4  205.8 139.6  49.8 145.6    1  1  1  0  1 -1 -1  0 -1  1 -1  0
     2        hR        999.0      86.6  196.7  213.4 136.0  85.7 100.3   -1  0 -1  0  1  1  0  0  1 -1 -1  0
     3        cP        999.0      62.1   72.4  186.7  90.0  90.0 100.2   -1  0  0  0  0  0 -1  0  0 -1  0  0
     5        cI        999.0     196.7   86.6  200.2  75.1  30.2  79.7   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     4        hR        999.0     103.4  196.7  205.8 141.8  86.3 103.2   -1  0  1  0  1 -1  0  0  1  1  1  0
     6        tI        999.0     200.2  196.7   86.6  79.7  75.1  30.2    0 -1 -1  0 -1 -1  0  0 -1  0 -1  0
     7        tI        999.0     196.7   86.6  200.2  75.1  30.2  79.7   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     8        oI        999.0      86.6  196.7  200.2  30.2  75.1  79.7    1  0  1  0  1  1  0  0  0  1  1  0
    16        oF        999.0      86.6  103.4  383.3  92.0 103.1  81.1    1  0  1  0 -1  0  1  0 -1 -2 -1  0
    17        mC        999.0     103.4   86.6  186.7  90.0  90.0  81.1   -1  0  1  0  1  0  1  0  0  1  0  0
    18        tI        999.0     200.2  205.8   62.1  76.2  93.7 137.6    0  1  1  0  1 -1  1  0  1  0  0  0
    19        oI        999.0      62.1  200.2  205.8  42.4  76.2  86.3   -1  0  0  0  0  1  1  0 -1  1 -1  0
    20        mC        999.0     200.2  200.2   62.1  93.7  93.7 137.6    0  1  1  0  0 -1  1  0  1  0  0  0
    21        tP        999.0      72.4  186.7   62.1  90.0 100.2  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    22        hP        999.0      72.4  186.7   62.1  90.0 100.2  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    23        oC        999.0     200.2  200.2   62.1  86.3  93.7  42.4    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    24        hR        999.0     237.6  200.2   62.1  86.3  81.2  58.1   -1 -1 -2  0  0 -1  1  0 -1  0  0  0
    25        mC        999.0     200.2  200.2   62.1  86.3  93.7  42.4    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    43        mI        999.0      86.6  387.4   72.4 102.4 135.1  87.6   -1  0 -1  0  1 -2 -1  0  0  0  1  0

 ****** LATTICE SYMMETRY IMPLICATED BY SPACE GROUP SYMMETRY ******

BRAVAIS-            POSSIBLE SPACE-GROUPS FOR PROTEIN CRYSTALS
  TYPE                     [SPACE GROUP NUMBER,SYMBOL]
  aP      [1,P1]
  mP      [3,P2] [4,P2(1)]
 mC,mI    [5,C2]
  oP      [16,P222] [17,P222(1)] [18,P2(1)2(1)2] [19,P2(1)2(1)2(1)]
  oC      [21,C222] [20,C222(1)]
  oF      [22,F222]
  oI      [23,I222] [24,I2(1)2(1)2(1)]
  tP      [75,P4] [76,P4(1)] [77,P4(2)] [78,P4(3)] [89,P422] [90,P42(1)2]
          [91,P4(1)22] [92,P4(1)2(1)2] [93,P4(2)22] [94,P4(2)2(1)2]
          [95,P4(3)22] [96,P4(3)2(1)2]
  tI      [79,I4] [80,I4(1)] [97,I422] [98,I4(1)22]
  hP      [143,P3] [144,P3(1)] [145,P3(2)] [149,P312] [150,P321] [151,P3(1)12]
          [152,P3(1)21] [153,P3(2)12] [154,P3(2)21] [168,P6] [169,P6(1)]
          [170,P6(5)] [171,P6(2)] [172,P6(4)] [173,P6(3)] [177,P622]
          [178,P6(1)22] [179,P6(5)22] [180,P6(2)22] [181,P6(4)22] [182,P6(3)22]
  hR      [146,R3] [155,R32]
  cP      [195,P23] [198,P2(1)3] [207,P432] [208,P4(2)32] [212,P4(3)32]
          [213,P4(1)32]
  cF      [196,F23] [209,F432] [210,F4(1)32]
  cI      [197,I23] [199,I2(1)3] [211,I432] [214,I4(1)32]


 Several solutions exist for the given space group that explain the data. The
 selected solution (marked "*") correlates best with the given reference data set.
 The independent settings of the cell can be obtained by application of the
 REINDEXING TRANSFORMATION to the original indices H,K,L from file INTEGRATE.HKL

 CORRELATION = correlation factor with the reference data set
 NPAIR       = number of unique reflection pairs correlated 
 Rmeas=Rrim  = redundancy independent R-factor (intensities)
               For definition see:
               Rmeas : Diederichs & Karplus,
                       Nature Struct. Biol. 4, 269-275 (1997);
               Rrim  : Weiss & Hilgenfeld,
                       J.Appl.Cryst. 30,203-205 (1997).
 COMPARED    = number of reflections used for calculating Rmeas
 ESD         = Agreement with given unit cell geometry

 CORRELATION  NPAIR  Rmeas  COMPARED  ESD       REINDEX TRANSFORMATION

  *   0.00        0   20.7   14885    0.01  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=    3
 UNIT_CELL_CONSTANTS=    61.80   185.60    72.00  90.000 100.200  90.000
 UNIT_CELL_A-AXIS=    56.069    18.596    18.992
 UNIT_CELL_B-AXIS=   -14.487  -109.797   150.272
 UNIT_CELL_C-AXIS=    18.388   -57.368   -40.143
 REFLECTIONS H,K,L FROM "INTEGRATE.HKL" WILL BE REINDEXED BY
 H' = (  1*H +  0*K +  0*L)/ 1 +  0
 K' = (  0*H +  1*K +  0*L)/ 1 +  0
 L' = (  0*H +  0*K +  1*L)/ 1 +  0

 *****  2 EQUIVALENT POSITIONS IN SPACE GROUP #  3 *****

    If x',y',z' is an equivalent position to x,y,z, then
        x'=x*ML(1)+y*ML( 2)+z*ML( 3)+ML( 4)/12.0
        y'=x*ML(5)+y*ML( 6)+z*ML( 7)+ML( 8)/12.0
        z'=x*ML(9)+y*ML(10)+z*ML(11)+ML(12)/12.0

    #    1  2  3  4    5  6  7  8    9 10 11 12 
    1    1  0  0  0    0  1  0  0    0  0  1  0
    2   -1  0  0  0    0  1  0  0    0  0 -1  0

 ******************************************************************************
      MEAN DISCREPANCIES BETWEEN OBSERVED AND CALCULATED SPOT LOCATIONS
 ******************************************************************************

 The discrepancies in X- and Y-coordinates of the spots are depicted in the
 two images DX-CORRECTIONS.cbf and DY-CORRECTIONS.cbf for inspection with
 the XDS-Viewer.



 ******************************************************************************
  REFINEMENT OF DIFFRACTION PARAMETERS USING ALL IMAGES
 ******************************************************************************


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM  65516 INDEXED SPOTS
 REFINED PARAMETERS:   DISTANCE BEAM ORIENTATION CELL AXIS
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.27
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.21
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.040   186.911    72.363  90.000 100.183  90.000
 E.S.D. OF CELL PARAMETERS  1.3E-01 4.3E-01 1.7E-01 0.0E+00 4.0E-02 0.0E+00
 REC. CELL PARAMETERS   0.016377  0.005350  0.014040  90.000  79.817  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.090    18.591    18.901
 COORDINATES OF UNIT CELL B-AXIS   -14.250  -109.975   150.460
 COORDINATES OF UNIT CELL C-AXIS    18.382   -57.322   -40.158
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=    3
 UNIT_CELL_CONSTANTS=    62.04   186.91    72.36  90.000 100.183  90.000


 ******************************************************************************
        MEAN INTENSITY AS FUNCTION OF SPINDLE POSITION WITHIN DATA IMAGE
 ******************************************************************************

 This statistics could serve as a diagnostic tool for detecting
 shutter problems (suggested by Kay Diederichs).
 Data are corrected for this effect if PATCH_SHUTTER_PROBLEM=TRUE in XDS.INP.
 ===> Selected : PATCH_SHUTTER_PROBLEM=FALSE

 INTERVAL = Angular interval in units of the oscillation range
            covered by a data image.
 NUMBER   = Number of reflections in each interval
 INTENSITY= Mean intensity of the reflections in each interval
 FACTOR   = Correction factor applied to intensities;
            INTENSITY(overall)/INTENSITY(interval)

   INTERVAL     NUMBER      INTENSITY      FACTOR
 -0.05  0.05     32582         56.783       0.980
  0.05  0.15     32309         55.684       0.999
  0.15  0.25     32707         56.139       0.991
  0.25  0.35     32201         54.954       1.012
  0.35  0.45     32427         55.526       1.002
  0.45  0.55     32155         56.289       0.988
  0.55  0.65     32373         54.426       1.022
  0.65  0.75     32460         55.409       1.004
  0.75  0.85     32325         55.436       1.004
  0.85  0.95     32753         55.642       1.000
 
 -0.05  0.95    324292         55.630       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED      720
 DEGREES OF FREEDOM OF CHI^2 FIT             44545.9
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.005
 NUMBER OF CYCLES CARRIED OUT                      5

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     1.2 XMAX=   898.5 NXBIN=   36
 YMIN= 0.00023 YMAX= 0.45672 NYBIN=   20
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS      68509


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED     1369
 DEGREES OF FREEDOM OF CHI^2 FIT             44471.1
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.997
 NUMBER OF CYCLES CARRIED OUT                      4

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     5.8 XMAX=  2456.9 NXBIN=   37
 YMIN=    12.9 YMAX=  2515.5 NYBIN=   37
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS      68509


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED      468
 DEGREES OF FREEDOM OF CHI^2 FIT             44582.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.984
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     1.2 XMAX=   898.5 NXBIN=   36
 DETECTOR_SURFACE_POSITION=    1231    1264
 DETECTOR_SURFACE_POSITION=    1648    1689
 DETECTOR_SURFACE_POSITION=     815    1689
 DETECTOR_SURFACE_POSITION=     815     839
 DETECTOR_SURFACE_POSITION=    1648     839
 DETECTOR_SURFACE_POSITION=    2173    1663
 DETECTOR_SURFACE_POSITION=    1622    2226
 DETECTOR_SURFACE_POSITION=     841    2226
 DETECTOR_SURFACE_POSITION=     289    1663
 DETECTOR_SURFACE_POSITION=     289     866
 DETECTOR_SURFACE_POSITION=     841     302
 DETECTOR_SURFACE_POSITION=    1622     302
 DETECTOR_SURFACE_POSITION=    2173     866
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS      68509


 ******************************************************************************
    CORRECTION PARAMETERS FOR THE STANDARD ERROR OF REFLECTION INTENSITIES
 ******************************************************************************

 The variance v0(I) of the intensity I obtained from counting statistics is
 replaced by v(I)=a*(v0(I)+b*I^2). The model parameters a, b are chosen to
 minimize the discrepancies between v(I) and the variance estimated from
 sample statistics of symmetry related reflections. This model implicates
 an asymptotic limit ISa=1/SQRT(a*b) for the highest I/Sigma(I) that the
 experimental setup can produce (Diederichs (2010) Acta Cryst D66, 733-740).

     a        b          ISa
 1.503E+00  1.964E-02    5.82


 ******************************************************************************
      STANDARD ERROR OF REFLECTION INTENSITIES AS FUNCTION OF RESOLUTION
      FOR DATA SET  XDS_ASCII.HKL                                     
 ******************************************************************************

  I/Sigma  = mean intensity/Sigma of a reflection in shell
  Chi^2    = goodness of fit between sample variances of
             symmetry-related intensities and their errors
             (Chi^2 = 1 for perfect agreement)
  R-FACTOR
  observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
  expected = expected R-FACTOR derived from Sigma(I)

   NUMBER  = number of reflections in resolution shell
             used for calculation of R-FACTOR
  ACCEPTED = number of accepted reflections
  REJECTED = number of rejected reflections (MISFITS),
             recognized by comparison with symmetry-related
             reflections.

 RESOLUTION RANGE  I/Sigma  Chi^2  R-FACTOR  R-FACTOR  NUMBER ACCEPTED REJECTED
                                   observed  expected

   66.555  14.448     5.60   0.93     11.00     11.63     913     930      41
   14.448  10.339     5.62   0.93     11.80     11.64    1498    1524      95
   10.339   8.476     5.52   1.02     11.31     11.82    2252    2265     136
    8.476   7.355     5.26   1.04     11.66     11.92    2266    2283     128
    7.355   6.587     4.93   0.98     12.07     12.62    2875    2892     114
    6.587   6.018     4.68   1.05     12.57     13.11    3239    3259     121
    6.018   5.575     4.53   1.15     14.09     13.11    3162    3186     116
    5.575   5.217     4.56   1.14     14.79     13.38    3667    3680     160
    5.217   4.920     4.85   1.09     12.83     12.77    4053    4062     182
    4.920   4.669     4.99   1.11     12.95     12.54    4146    4165     196
    4.669   4.453     5.08   1.19     13.69     12.30    4038    4052     180
    4.453   4.264     4.90   1.20     13.75     12.66    4495    4508     185
    4.264   4.097     4.80   1.14     14.01     12.88    4743    4773     186
    4.097   3.949     4.68   1.13     14.18     13.22    4969    5015     158
    3.949   3.815     4.58   1.13     14.19     13.18    4680    4728     149
    3.815   3.694     4.38   1.11     14.22     13.61    4724    4793     166
    3.694   3.585     4.18   1.08     14.86     14.22    5018    5074     132
    3.585   3.484     3.98   1.04     15.36     14.95    5214    5278     128
    3.484   3.391     3.66   1.08     16.75     16.03    5334    5423     113
    3.391   3.305     3.37   1.02     18.00     17.51    5359    5434      91
    3.305   3.226     3.23   1.02     18.51     17.81    4788    4908      70
    3.226   3.152     2.96   1.00     19.41     19.09    5239    5384      58
    3.152   3.083     2.56   0.92     22.49     23.03    5171    5339      66
    3.083   3.018     2.31   0.88     24.76     25.53    5191    5428      45
    3.018   2.957     2.05   0.89     28.60     30.01    5060    5286      42
    2.957   2.900     1.91   0.87     30.91     33.21    5004    5266      54
    2.900   2.846     1.76   0.88     33.27     35.35    4506    4762      42
    2.846   2.795     1.59   0.83     37.11     40.36    4698    4993      24
    2.795   2.746     1.55   0.85     39.23     42.11    4731    5037      32
    2.746   2.700     1.40   0.82     43.09     46.89    4838    5103      29
    2.700   2.656     1.22   0.79     51.37     57.57    4833    5175      14
    2.656   2.614     1.16   0.80     53.98     60.48    4923    5211      17
    2.614   2.575     1.08   0.76     57.14     65.61    4566    4907       7
    2.575   2.536     1.01   0.78     60.63     68.75    4302    4623       8
    2.536   2.500     0.94   0.77     66.04     75.63    4578    4905      22
    2.500   2.465     0.89   0.79     71.03     79.93    4554    4883       7
    2.465   2.432     0.81   0.75     77.31     88.11    4626    4947      15
    2.432   2.399     0.79   0.77     80.16     91.19    4664    5001       8
    2.399   2.369     0.72   0.73     87.07    101.51    4636    4971      17
    2.369   2.339     0.70   0.76     91.15    104.75    4419    4785       3
    2.339   2.310     0.63   0.75     99.22    114.33    4151    4514       8
    2.310   2.282     0.59   0.71    105.36    125.38    4353    4702       7
    2.282   2.256     0.53   0.67    120.22    146.42    4464    4835       4
    2.256   2.230     0.50   0.69    123.59    148.64    4472    4841       4
    2.230   2.205     0.45   0.67    140.31    169.64    4338    4727       7
    2.205   2.181     0.45   0.66    142.25    175.36    4442    4872       2
    2.181   2.158     0.44   0.69    142.09    171.28    4428    4843       4
    2.158   2.135     0.39   0.67    159.61    193.84    4168    4650       2
    2.135   2.113     0.38   0.65    165.59    204.66    4015    4526       2
    2.113   2.092     0.39   0.67    163.92    198.09    4236    4704       3
    2.092   2.072     0.31   0.65    205.60    250.97    4277    4791       1
    2.072   2.052     0.28   0.64    227.73    283.47    4424    4896       4
    2.052   2.032     0.28   0.68    244.22    294.00    4405    4909       0
    2.032   2.013     0.26   0.64    239.40    298.48    4376    4883       1
    2.013   1.995     0.27   0.64    241.80    299.90    4047    4552       0
    1.995   1.977     0.22   0.63    285.65    358.20    3833    4338       1
    1.977   1.960     0.26   0.61    233.51    294.58    3677    4204       1
    1.960   1.943     0.25   0.62    234.88    297.30    3563    4101       2
    1.943   1.926     0.17   0.62    336.41    425.32    3446    3969       3
    1.926   1.910     0.21   0.64    272.07    336.61    3361    3899       2
    1.910   1.894     0.17   0.64    331.49    409.55    3176    3702       0
    1.894   1.879     0.17   0.63    326.78    407.55    3050    3592       0
    1.879   1.864     0.16   0.61    337.45    428.66    2561    3171       5
    1.864   1.849     0.16   0.58    319.73    418.09    2331    2911       0
    1.849   1.835     0.15   0.58    338.94    446.11    2205    2778       0
    1.835   1.821     0.17   0.59    313.98    405.90    2028    2600       0
    1.821   1.808     0.15   0.58    331.43    428.46    1882    2425       0
    1.808   1.794     0.17   0.63    323.38    401.77    1787    2329       0
    1.794   1.781     0.16   0.67    329.18    393.66    1651    2171       0
    1.781   1.768     0.14   0.59    406.98    526.08    1529    2021       0
    1.768   1.756     0.13   0.64    350.03    431.71    1460    1944       0
    1.756   1.744     0.13   0.59    370.79    470.06    1349    1825       1
    1.744   1.732     0.10   0.61    424.89    547.19    1215    1664       0
    1.732   1.720     0.14   0.62    331.44    416.97    1190    1638       0
    1.720   1.708     0.16   0.66    370.02    447.52    1079    1509       0
    1.708   1.697     0.11   0.58    521.10    676.40     971    1396       0
    1.697   1.686     0.15   0.58    312.17    400.39     922    1329       0
    1.686   1.675     0.14   0.64    375.04    465.99     877    1265       0
    1.675   1.665     0.11   0.58    468.75    615.01     825    1196       0
    1.665   1.654     0.14   0.54    371.98    515.97     737    1083       0
    1.654   1.644     0.12   0.59    397.76    512.25     643     999       0
    1.644   1.634     0.14   0.59    261.75    334.32     560     909       0
    1.634   1.624     0.13   0.63    285.23    350.77     532     866       0
    1.624   1.614     0.15   0.59    287.97    368.30     479     796       0
    1.614   1.605     0.12   0.57    328.78    439.82     401     723       0
    1.605   1.596     0.09   0.65    470.77    572.68     351     658       0
    1.596   1.586     0.18   0.48    154.75    218.47     268     565       0
    1.586   1.577     0.13   0.53    284.50    392.27     221     523       0
    1.577   1.568     0.08   0.74   1281.82   1503.47     186     465       0
    1.568   1.560     0.11   0.79    -99.90    -99.90     127     402       0
    1.560   1.551     0.23   0.44    165.93    240.28      88     355       0
    1.551   1.543     0.13   0.54    236.11    332.08      42     279       0
    1.543   1.534     0.18   0.34    -99.90    -99.90       6     227       0
    1.534   1.526     0.17   0.11     51.19    120.81       2     182       0
    1.526   1.518     0.15   0.00    -99.90    -99.90       0     141       0
    1.518   1.510     0.19   0.00    -99.90    -99.90       0      99       0
    1.510   1.502     0.09   0.00    -99.90    -99.90       0      69       0
    1.502   1.495     0.02   0.00    -99.90    -99.90       0      37       0
    1.495   1.487     0.04   0.00    -99.90    -99.90       0      24       0
    1.487   1.480    -0.06   0.00    -99.90    -99.90       0       9       0
   --------------------------------------------------------------------------
   66.555   1.480     1.71   0.83     23.72     25.36  291179  320871    3421
 


 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR VARIOUS SUBSETS OF INCLUDED DATA IMAGES
 ******************************************************************************

 Data set statistics is reported below several times, each with a different
 upper limit on the number of images included. This provides the user with
 the information for deciding which data images should be excluded from the
 final data set because of radiation damage or other defects. If the user
 decides for a subset of "good" images that differs from the specification
 as given by the input parameter (see DATA_RANGE= in XDS.INP) the CORRECT-
 step of XDS must be repeated with the new parameter values for DATA_RANGE=.


 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1     101
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     4.43        3673    3060     10156       30.1%       4.5%     10.5%     1226    5.31      6.3%    99.2*     0    0.000       0
     3.14        6792    6495     18271       35.5%       8.6%     12.9%      594    4.01     12.2%    97.1*     0    0.000       0
     2.56        8613    8556     23617       36.2%      18.6%     21.4%      114    1.76     26.4%    78.1*     0    0.000       0
     2.22        8866    8805     27876       31.6%      18.9%     22.8%      122    0.82     26.7%    84.0*     0    0.000       0
     1.98        8298    8241     31673       26.0%      46.3%     49.3%      114    0.39     65.5%     9.9      0    0.000       0
     1.81        6237    6221     34907       17.8%     127.1%    130.3%       32    0.16    179.8%   -20.8      0    0.000       0
     1.68        2724    2724     37953        7.2%     -99.9%    -99.9%        0  -99.00    -99.9%     0.0      0    0.000       0
     1.57        1236    1236     40832        3.0%     -99.9%    -99.9%        0  -99.00    -99.9%     0.0      0    0.000       0
     1.48         338     338     43372        0.8%     -99.9%    -99.9%        0  -99.00    -99.9%     0.0      0    0.000       0
    total       46777   45676    268657       17.0%       6.0%     11.5%     2202    1.52      8.4%    98.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES   47029
 NUMBER OF REJECTED MISFITS                             252
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                  0
 NUMBER OF ACCEPTED OBSERVATIONS                      46777
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                45676



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1     201
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     4.43        7404    5248     10156       51.7%       4.9%     10.6%     4312    5.70      6.9%    99.1*     0    0.000       0
     3.14       13709   10762     18271       58.9%       9.6%     13.4%     5894    4.24     13.5%    96.3*     0    0.000       0
     2.56       17613   14527     23617       61.5%      35.0%     41.3%     6172    1.77     49.5%    71.1*     0    0.000       0
     2.22       19393   16446     27876       59.0%      86.9%    103.3%     5894    0.82    122.9%    23.8*     0    0.000       0
     1.98       20091   17334     31673       54.7%     181.1%    221.3%     5514    0.41    256.2%     4.3      0    0.000       0
     1.81       15760   14171     34907       40.6%     264.6%    338.6%     3178    0.20    374.2%     1.2      0    0.000       0
     1.68        7934    7202     37953       19.0%     383.9%    471.5%     1464    0.15    543.0%     1.3      0    0.000       0
     1.57        3731    3523     40832        8.6%     650.3%    926.6%      416    0.11    919.7%     5.7      0    0.000       0
     1.48         902     902     43372        2.1%     -99.9%    -99.9%        0  -99.00    -99.9%     0.0      0    0.000       0
    total      106537   90115    268657       33.5%      16.7%     23.2%    32844    1.40     23.6%    97.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES  107112
 NUMBER OF REJECTED MISFITS                             575
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                  0
 NUMBER OF ACCEPTED OBSERVATIONS                     106537
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                90115



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1     301
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     4.43       11003    6232     10156       61.4%       8.6%     11.2%     8731    6.33     11.3%    97.1*   -78    0.368      58
     3.14       20562   12127     18271       66.4%      13.4%     14.3%    15242    4.81     17.9%    93.3*     0    0.000       0
     2.56       26586   15987     23617       67.7%      40.0%     44.7%    18987    2.02     53.5%    69.4*     0    0.000       0
     2.22       29879   18222     27876       65.4%      98.3%    117.5%    20912    0.96    131.5%    28.7*     0    0.000       0
     1.98       31901   19693     31673       62.2%     207.0%    258.3%    21921    0.49    277.1%     8.4      0    0.000       0
     1.81       25337   16685     34907       47.8%     325.7%    409.8%    15834    0.24    441.4%     1.8      0    0.000       0
     1.68       13238    9559     37953       25.2%     336.9%    421.9%     6682    0.18    454.8%     0.1      0    0.000       0
     1.57        6352    5329     40832       13.1%     372.8%    488.3%     1928    0.13    512.4%     4.2      0    0.000       0
     1.48        1492    1470     43372        3.4%     293.6%    327.6%       44    0.17    415.3%   -30.5      0    0.000       0
    total      166350  105304    268657       39.2%      28.1%     33.6%   110281    1.56     37.6%    95.5*   -78    0.368      58


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES  167372
 NUMBER OF REJECTED MISFITS                            1022
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                  0
 NUMBER OF ACCEPTED OBSERVATIONS                     166350
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS               105304



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1     401
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     4.43       14552    6356     10156       62.6%      12.9%     12.1%    12822    7.10     15.9%    95.1*   -51    0.541    1033
     3.14       27261   12322     18271       67.4%      16.9%     15.3%    23071    5.44     21.0%    91.4*   -36    0.682    1676
     2.56       35371   16257     23617       68.8%      42.1%     45.4%    29438    2.28     52.4%    76.6*   -16    0.662    1996
     2.22       40164   18596     27876       66.7%      91.5%    108.2%    33714    1.10    114.6%    48.6*    -4    0.594    1817
     1.98       42894   20190     31673       63.7%     173.2%    215.4%    35993    0.57    218.1%    27.5*    -5    0.527    1726
     1.81       33142   17385     34907       49.8%     271.9%    345.4%    26472    0.28    353.4%     7.6     -4    0.494     552
     1.68       17088   10323     37953       27.2%     334.9%    418.1%    11515    0.19    437.6%     0.5     -1    0.558     147
     1.57        7948    5843     40832       14.3%     332.7%    434.3%     3871    0.14    451.3%     5.2     -6    0.632      14
     1.48        1767    1679     43372        3.9%     441.1%    514.6%      176    0.17    623.9%    -2.4      0    0.000       0
    total      220187  108951    268657       40.6%      33.2%     36.6%   177072    1.74     41.8%    94.4*   -21    0.600    8961


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES  221830
 NUMBER OF REJECTED MISFITS                            1643
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                  0
 NUMBER OF ACCEPTED OBSERVATIONS                     220187
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS               108951



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1     500
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     4.43       18046    6591     10156       64.9%      14.0%     12.4%    16902    7.77     17.0%    94.2*   -54    0.516    2034
     3.14       33849   12790     18271       70.0%      17.3%     15.3%    31328    5.98     21.3%    91.3*   -36    0.658    3603
     2.56       43940   16988     23617       71.9%      39.9%     42.9%    40524    2.50     49.2%    80.2*   -15    0.673    4680
     2.22       48027   19224     27876       69.0%      85.6%    100.6%    43378    1.21    105.2%    58.3*    -5    0.608    4757
     1.98       49167   20434     31673       64.5%     179.0%    220.8%    42880    0.59    218.8%    31.1*    -4    0.538    4963
     1.81       36780   17451     34907       50.0%     288.6%    365.3%    30506    0.29    364.0%     8.5     -4    0.512    2156
     1.68       19033   10368     37953       27.3%     353.5%    442.8%    13898    0.19    449.1%    -0.4      2    0.527     791
     1.57        8773    5901     40832       14.5%     339.4%    444.1%     5188    0.15    456.9%     4.2      2    0.592      67
     1.48        1841    1704     43372        3.9%     475.1%    582.6%      274    0.16    671.9%     2.6      0    0.000       0
    total      259456  111451    268657       41.5%      32.2%     34.6%   224878    1.92     39.7%    94.2*   -18    0.594   23051


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES  261733
 NUMBER OF REJECTED MISFITS                            2277
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                  0
 NUMBER OF ACCEPTED OBSERVATIONS                     259456
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS               111451



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1     600
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     4.43       21653    7846     10156       77.3%      13.4%     12.4%    20375    7.83     16.3%    95.2*   -54    0.512    2433
     3.14       40365   14831     18271       81.2%      16.6%     15.2%    37554    6.14     20.3%    92.5*   -36    0.652    4374
     2.56       49557   18425     23617       78.0%      38.7%     41.8%    45856    2.62     47.2%    81.9*   -15    0.670    5514
     2.22       51343   19636     27876       70.4%      86.5%    100.6%    47105    1.23    105.6%    57.5*    -5    0.608    5297
     1.98       51436   20534     31673       64.8%     181.6%    222.7%    45794    0.59    221.2%    29.5*    -4    0.536    5329
     1.81       37693   17507     34907       50.2%     290.0%    367.2%    31670    0.28    364.6%     8.8     -4    0.513    2225
     1.68       19250   10369     37953       27.3%     354.3%    443.8%    14179    0.19    448.9%     0.6      2    0.527     795
     1.57        8795    5902     40832       14.5%     339.8%    444.1%     5219    0.15    457.1%     4.5      2    0.592      67
     1.48        1841    1704     43372        3.9%     475.1%    582.6%      274    0.16    671.9%     2.6      0    0.000       0
    total      281933  116754    268657       43.5%      29.8%     32.0%   248026    2.10     36.5%    95.1*   -19    0.594   26034


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES  284665
 NUMBER OF REJECTED MISFITS                            2732
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                  0
 NUMBER OF ACCEPTED OBSERVATIONS                     281933
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS               116754



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1     700
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     4.43       25365    9640     10156       94.9%      12.0%     12.1%    23884    7.73     14.7%    96.4*   -54    0.514    2449
     3.14       45101   16371     18271       89.6%      15.2%     14.5%    43161    6.33     18.7%    94.1*   -36    0.654    4412
     2.56       50730   18694     23617       79.2%      36.8%     39.8%    47300    2.68     45.0%    84.6*   -15    0.670    5515
     2.22       51447   19649     27876       70.5%      86.4%    100.3%    47225    1.23    105.4%    57.4*    -5    0.608    5297
     1.98       51437   20535     31673       64.8%     181.6%    222.7%    45794    0.59    221.2%    29.5*    -4    0.536    5329
     1.81       37693   17507     34907       50.2%     290.0%    367.2%    31670    0.28    364.6%     8.5     -4    0.513    2225
     1.68       19250   10369     37953       27.3%     354.3%    443.8%    14179    0.19    448.9%     0.4      2    0.527     795
     1.57        8795    5902     40832       14.5%     339.8%    444.1%     5219    0.15    457.1%     4.0      2    0.592      67
     1.48        1841    1704     43372        3.9%     475.1%    582.6%      274    0.16    671.9%     3.0      0    0.000       0
    total      291659  120371    268657       44.8%      26.2%     28.7%   258706    2.26     32.2%    96.2*   -19    0.595   26089


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES  294646
 NUMBER OF REJECTED MISFITS                            2987
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                  0
 NUMBER OF ACCEPTED OBSERVATIONS                     291659
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS               120371



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1     800
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     4.43       29193    9950     10156       98.0%      12.1%     12.2%    28618    8.26     14.7%    96.1*   -55    0.472    3146
     3.14       49605   17180     18271       94.0%      15.3%     14.3%    48060    6.58     18.8%    94.0*   -37    0.638    4776
     2.56       51920   18980     23617       80.4%      35.7%     38.4%    48625    2.73     43.6%    85.1*   -15    0.669    5524
     2.22       51448   19649     27876       70.5%      86.4%    100.3%    47227    1.23    105.4%    57.3*    -5    0.608    5297
     1.98       51437   20535     31673       64.8%     181.6%    222.7%    45794    0.59    221.2%    29.1*    -4    0.536    5329
     1.81       37693   17507     34907       50.2%     290.0%    367.2%    31670    0.28    364.6%     7.6     -4    0.513    2225
     1.68       19250   10369     37953       27.3%     354.3%    443.8%    14179    0.19    448.9%     0.2      2    0.527     795
     1.57        8795    5902     40832       14.5%     339.8%    444.1%     5219    0.15    457.1%     3.7      2    0.592      67
     1.48        1841    1704     43372        3.9%     475.1%    582.6%      274    0.16    671.9%     2.8      0    0.000       0
    total      301182  121776    268657       45.3%      24.4%     26.3%   269666    2.39     30.0%    96.1*   -20    0.586   27159


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES  304342
 NUMBER OF REJECTED MISFITS                            3160
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                  0
 NUMBER OF ACCEPTED OBSERVATIONS                     301182
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS               121776



 ******************************************************************************
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       1     900)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

    0    0    2    35.612   0.8779E+02  0.1527E+02         5.75           1 
    0    0    3    23.741   0.9191E+02  0.1141E+02         8.05           2 
    0    0    4    17.806   0.3564E+02  0.4903E+01         7.27           2 
    0    0    5    14.245   0.8404E+01  0.2482E+01         3.39           2 
    0    0    6    11.871   0.7112E+03  0.8688E+02         8.19           2 
    0    0    7    10.175   0.2589E+03  0.3196E+02         8.10           2 
    0    0    8     8.903   0.7823E+03  0.9557E+02         8.19           2 
    0    0    9     7.914   0.7515E+02  0.1023E+02         7.35           2 
    0    0   10     7.122   0.3420E+03  0.4226E+02         8.09           2 
    0    0   11     6.475   0.5823E+02  0.8758E+01         6.65           2 
    0    0   12     5.935   0.5030E+03  0.6194E+02         8.12           2 
    0    0   13     5.479   0.5640E+02  0.1297E+02         4.35           1 
    0    0   14     5.087   0.5498E+02  0.9349E+01         5.88           2 
    0    0   15     4.748   0.1851E+03  0.2432E+02         7.61           2 
    0    0   16     4.451   0.4912E+02  0.9653E+01         5.09           2 
    0    0   17     4.190   0.6085E+02  0.1548E+02         3.93           1 
    0    0   18     3.957   0.5469E+03  0.9662E+02         5.66           1 
    0    4    0    46.728   0.6234E+04  0.1071E+04         5.82           1 
    0    5    0    37.382   0.7855E+01  0.2537E+01         3.10           1 
    0    6    0    31.152   0.6660E+03  0.8102E+02         8.22           2 
    0    7    0    26.702   0.4000E+02  0.5414E+01         7.39           2 
    0    8    0    23.364   0.1365E+04  0.1660E+03         8.22           2 
    0    9    0    20.768   0.4093E+02  0.5826E+01         7.03           2 
    0   10    0    18.691   0.8190E+03  0.9973E+02         8.21           2 
    0   11    0    16.992   0.1091E+02  0.3589E+01         3.04           2 
    0   12    0    15.576   0.1048E+03  0.1345E+02         7.80           2 
    0   13    0    14.378   0.2970E+02  0.7596E+01         3.91           1 
    0   14    0    13.351   0.3109E+03  0.5407E+02         5.75           1 
    0   15    0    12.461   0.5840E+02  0.8547E+01         6.83           2 
    0   16    0    11.682   0.2833E+03  0.3505E+02         8.08           2 
    0   17    0    10.995   0.3293E+02  0.6573E+01         5.01           2 
    0   18    0    10.384   0.2484E+03  0.3098E+02         8.02           2 
    0   19    0     9.837   0.2108E+02  0.6109E+01         3.45           2 
    0   20    0     9.346   0.2340E+03  0.2942E+02         7.95           2 
    0   21    0     8.901   0.3029E+02  0.7103E+01         4.26           2 
    0   22    0     8.496   0.1339E+03  0.2524E+02         5.30           1 
    0   23    0     8.127   0.4699E+02  0.1257E+02         3.74           1 
    0   24    0     7.788   0.9741E+02  0.1386E+02         7.03           2 
    0   25    0     7.476   0.3740E+01  0.6727E+01         0.56           2 
    0   26    0     7.189   0.4742E+02  0.9325E+01         5.09           2 
    0   27    0     6.923   0.4207E+02  0.8983E+01         4.68           2 
    0   28    0     6.675   0.1946E+03  0.2535E+02         7.68           2 
    0   29    0     6.445   0.3358E+02  0.8900E+01         3.77           2 
    0   30    0     6.230   0.1783E+03  0.2356E+02         7.57           2 
    0   31    0     6.029   0.7520E+02  0.1772E+02         4.24           1 
    0   32    0     5.841   0.3279E+03  0.5839E+02         5.62           1 
    0   33    0     5.664   0.6432E+02  0.1197E+02         5.37           2 
    0   34    0     5.497   0.1095E+03  0.1653E+02         6.62           2 
    0   35    0     5.340   0.8097E+02  0.1400E+02         5.78           2 
    0   36    0     5.192   0.1135E+03  0.1699E+02         6.68           2 
    0   37    0     5.052   0.1465E+03  0.2093E+02         7.00           2 
    0   38    0     4.919   0.6927E+02  0.1282E+02         5.40           2 
    0   39    0     4.793   0.1263E+03  0.1881E+02         6.71           2 
    0   40    0     4.673   0.3701E+03  0.6572E+02         5.63           1 
    0   41    0     4.559   0.1985E+03  0.2778E+02         7.15           2 
    0   42    0     4.450   0.1570E+03  0.2229E+02         7.04           2 
    0   43    0     4.347   0.5791E+02  0.1293E+02         4.48           2 
    0   44    0     4.248   0.1077E+03  0.1666E+02         6.47           2 
    0   45    0     4.154   0.4925E+02  0.1173E+02         4.20           2 
    0   46    0     4.063   0.1035E+03  0.1674E+02         6.18           2 
    0   47    0     3.977   0.3158E+02  0.1077E+02         2.93           2 
    0   48    0     3.894   0.1837E+03  0.3576E+02         5.14           1 
    0   49    0     3.815   0.3144E+02  0.1547E+02         2.03           1 
    0   50    0     3.738   0.4118E+02  0.1162E+02         3.54           2 
    0   51    0     3.665   0.3749E+02  0.1178E+02         3.18           2 
    0   52    0     3.594   0.1161E+03  0.1833E+02         6.33           2 
    0   53    0     3.527   0.3939E+02  0.1206E+02         3.27           2 
    0   54    0     3.461   0.1212E+02  0.1062E+02         1.14           2 
    0   55    0     3.398   0.2551E+01  0.1068E+02         0.24           2 
    0   56    0     3.338   0.1967E+02  0.1588E+02         1.24           1 
    0   57    0     3.279   0.1666E+02  0.1656E+02         1.01           1 
    0   58    0     3.223   0.4366E+02  0.1266E+02         3.45           2 
    0   59    0     3.168   0.2879E+02  0.1198E+02         2.40           2 
    0   60    0     3.115   0.2154E+02  0.1171E+02         1.84           2 
    0   61    0     3.064  -0.1562E+02  0.1128E+02        -1.39           2 
    0   62    0     3.015   0.4348E+01  0.1124E+02         0.39           2 
    0   63    0     2.967   0.1140E+02  0.1177E+02         0.97           2 
    0   64    0     2.920   0.2153E+02  0.1281E+02         1.68           2 
    0   65    0     2.876   0.2273E+02  0.1930E+02         1.18           1 
    0   66    0     2.832   0.4373E+02  0.1415E+02         3.09           2 
    0   67    0     2.790   0.9199E+01  0.1275E+02         0.72           2 
    0   68    0     2.749   0.1373E+02  0.1263E+02         1.09           2 
    0   69    0     2.709   0.2004E+02  0.1278E+02         1.57           2 
    0   70    0     2.670  -0.9680E+01  0.1209E+02        -0.80           2 
    0   71    0     2.633   0.3242E+01  0.1176E+02         0.28           2 
    0   72    0     2.596   0.1762E+01  0.1618E+02         0.11           1 
    0   73    0     2.560   0.3619E+01  0.1597E+02         0.23           1 
    0   74    0     2.526   0.4376E+01  0.1208E+02         0.36           2 
    0   75    0     2.492  -0.5426E+01  0.1247E+02        -0.44           2 
    0   76    0     2.459   0.3496E+01  0.1278E+02         0.27           2 
    0   77    0     2.427  -0.2248E+02  0.1342E+02        -1.68           2 
    0   78    0     2.396  -0.4082E+01  0.1337E+02        -0.31           2 
    0   79    0     2.366   0.1869E+01  0.1381E+02         0.14           2 
    0   80    0     2.336   0.1329E+02  0.2019E+02         0.66           1 
    0   81    0     2.308   0.7309E+01  0.1953E+02         0.37           1 
    0   82    0     2.279  -0.3615E+01  0.1432E+02        -0.25           2 
    0   83    0     2.252   0.1197E+02  0.1453E+02         0.82           2 
    0   84    0     2.225   0.4139E+01  0.1429E+02         0.29           2 
    0   85    0     2.199   0.1626E+02  0.1442E+02         1.13           2 
    0   86    0     2.173   0.1379E+02  0.1414E+02         0.98           2 
    0   87    0     2.148  -0.1045E+02  0.1937E+02        -0.54           1 
    0   88    0     2.124  -0.1487E+02  0.1479E+02        -1.01           2 
    0   89    0     2.100   0.1550E+00  0.1455E+02         0.01           2 
    0   90    0     2.077  -0.1175E+02  0.1540E+02        -0.76           2 
    0   91    0     2.054  -0.6884E+01  0.1513E+02        -0.46           2 
    0   92    0     2.032   0.1874E+02  0.1553E+02         1.21           2 
    0   93    0     2.010   0.9327E+01  0.1541E+02         0.61           2 
    0   94    0     1.988  -0.6457E+01  0.2096E+02        -0.31           1 
    0   95    0     1.967  -0.1419E+02  0.2081E+02        -0.68           1 
    0   96    0     1.947   0.2461E+02  0.2165E+02         1.14           1 
    0   97    0     1.927  -0.1884E+02  0.2096E+02        -0.90           1 
    0   98    0     1.907  -0.1085E+02  0.2049E+02        -0.53           1 
    0   99    0     1.888   0.1713E+02  0.2127E+02         0.81           1 
    4    0    0    15.266   0.4179E+03  0.7193E+02         5.81           1 
    5    0    0    12.213   0.2517E+04  0.4325E+03         5.82           1 
    6    0    0    10.177   0.3494E+04  0.6005E+03         5.82           1 
    7    0    0     8.723   0.1183E+04  0.1441E+03         8.21           2 
    8    0    0     7.633   0.1361E+03  0.1690E+02         8.05           2 
    9    0    0     6.785   0.7805E+02  0.1011E+02         7.72           2 
   10    0    0     6.106   0.9049E+03  0.1104E+03         8.19           2 
   11    0    0     5.551   0.2150E+03  0.2667E+02         8.06           2 
   12    0    0     5.089   0.4116E+02  0.6199E+01         6.64           2 
   13    0    0     4.697   0.6753E+03  0.8272E+02         8.16           2 
   14    0    0     4.362   0.4978E+03  0.6092E+02         8.17           2 
   15    0    0     4.071   0.4175E+02  0.6680E+01         6.25           2 
   16    0    0     3.816   0.3873E+03  0.4849E+02         7.99           2 
   17    0    0     3.592   0.4988E+03  0.8642E+02         5.77           1 
   18    0    0     3.392   0.4749E+02  0.1151E+02         4.12           1 
   19    0    0     3.214   0.3130E+03  0.5478E+02         5.71           1 
   20    0    0     3.053   0.3761E+03  0.6545E+02         5.75           1 
   21    0    0     2.908   0.2974E+03  0.5205E+02         5.71           1 
   22    0    0     2.776   0.7093E+02  0.1505E+02         4.71           1 
   23    0    0     2.655   0.5861E+02  0.1308E+02         4.48           1 
   24    0    0     2.544   0.1559E+03  0.2834E+02         5.50           1 
   25    0    0     2.443   0.4714E+02  0.1217E+02         3.87           1 
   26    0    0     2.349   0.8643E+02  0.1783E+02         4.85           1 
   27    0    0     2.262   0.1474E+02  0.9873E+01         1.49           1 
   28    0    0     2.181   0.9377E+02  0.1916E+02         4.89           1 
   29    0    0     2.106   0.3088E+02  0.1089E+02         2.83           1 
   30    0    0     2.035   0.2442E+02  0.9966E+01         2.45           1 
   31    0    0     1.970   0.5239E+02  0.1266E+02         4.14           1 

          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     4.43       32808    9987     10156       98.3%      12.7%     12.4%    32619    8.78     15.2%    96.2*   -59    0.440    4363
     3.14       55988   17425     18271       95.4%      15.3%     14.5%    55195    7.00     18.3%    94.9*   -40    0.609    6118
     2.56       56883   19640     23617       83.2%      34.7%     37.5%    53824    2.85     41.9%    86.9*   -15    0.669    5690
     2.22       53529   19754     27876       70.9%      87.3%    101.5%    49684    1.24    106.0%    58.1*    -5    0.608    5297
     1.98       52662   20598     31673       65.0%     184.2%    225.9%    47345    0.59    223.9%    29.9*    -4    0.536    5329
     1.81       38552   17793     34907       51.0%     293.7%    371.9%    32439    0.28    368.6%     8.2     -4    0.513    2225
     1.68       19604   10493     37953       27.6%     357.9%    448.2%    14448    0.19    452.1%    -0.1      2    0.527     795
     1.57        8969    5988     40832       14.7%     337.7%    439.8%     5342    0.15    452.7%     3.6      3    0.592      67
     1.48        1876    1734     43372        4.0%     471.5%    578.2%      283    0.16    665.9%     2.0      0    0.000       0
    total      320871  123412    268657       45.9%      23.7%     25.4%   291179    2.51     28.7%    96.4*   -23    0.573   29884


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES  324292
 NUMBER OF REJECTED MISFITS                            3421
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                  0
 NUMBER OF ACCEPTED OBSERVATIONS                     320871
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS               123412



 ******************************************************************************
    WILSON STATISTICS OF DATA SET  "XDS_ASCII.HKL"
 ******************************************************************************

 Data is divided into resolution shells and a straight line 
 A - 2*B*SS is fitted to log<I>, where
   RES    = mean resolution (Angstrom) in shell
   SS     = mean of (sin(THETA)/LAMBDA)**2 in shell
   <I>    = mean reflection intensity in shell
   BO     = (A - log<I>)/(2*SS)
    #     = number of reflections in resolution shell

   WILSON LINE (using all data) : A=   5.198 B=  23.822 CORRELATION=  0.95
      #      RES      SS        <I>       log(<I>)       BO
    3300     8.275   0.004  2.4930E+02       5.519     -43.9
    5919     5.201   0.009  1.6924E+02       5.131       3.6
    7602     4.060   0.015  1.9259E+02       5.261      -2.1
    8584     3.439   0.021  1.1346E+02       4.731      11.0
    9219     3.037   0.027  5.2900E+01       3.968      22.7
    9319     2.747   0.033  2.8127E+01       3.337      28.1
    9329     2.528   0.039  1.8670E+01       2.927      29.0
    9362     2.354   0.045  1.3231E+01       2.583      29.0
    9545     2.211   0.051  9.3820E+00       2.239      28.9
    9778     2.092   0.057  6.9767E+00       1.943      28.5
    9723     1.990   0.063  4.9819E+00       1.606      28.5
    8853     1.903   0.069  3.3262E+00       1.202      28.9
    6702     1.826   0.075  3.0669E+00       1.121      27.2
    5183     1.757   0.081  2.2928E+00       0.830      27.0
    4093     1.695   0.087  2.6471E+00       0.973      24.3
    3120     1.639   0.093  2.5874E+00       0.951      22.8
    2258     1.589   0.099  2.3259E+00       0.844      22.0
    1254     1.545   0.105  3.3577E+00       1.211      19.0
     268     1.507   0.110  2.1030E+00       0.743      20.2


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

     142     8.275        1.876        3.164        4.881
     152     5.201        1.343        1.586        1.887
     131     4.060        1.340        1.027        0.698
      81     3.439        3.611        5.691        8.161
      59     3.037        3.287        4.065        4.276
      31     2.747        2.504        2.911        3.064
      20     2.528        3.376        3.877        3.769
      15     2.354       11.880       35.905       93.143
      16     2.211       10.873       25.042       47.880
      13     2.092       26.444      145.894      622.324
       9     1.990        7.214       11.070       14.398
       2     1.903        0.069        0.004        0.001
     671   overall        3.038        6.899       18.686


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

    3158     8.275        1.070        1.540        2.820
    5767     5.201        0.922        0.978        1.183
    7471     4.060        0.846        0.772        0.791
    8503     3.439        0.911        0.936        1.126
    9160     3.037        1.132        1.727        3.076
    9288     2.747        1.124        1.514        2.241
    9309     2.528        1.469        3.507       12.189
    9347     2.354        1.491        2.622        5.788
    9529     2.211        1.808        2.932        5.829
    9765     2.092        2.428        4.601       13.644
    9714     1.990        3.256        3.975        8.545
    8851     1.903        6.189        7.365       28.177
    6702     1.826        7.091        6.800       29.877
    5183     1.757       12.910       11.231      103.011
    4093     1.695       10.652       11.443       65.351
    3120     1.639       12.807       12.770       97.292
    2258     1.589       16.885       19.023      161.475
    1254     1.545        9.470        9.550       46.219
     268     1.507       28.999       48.791      425.540
  122740   overall        3.701        4.452       21.271

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

     142     8.275  0.20  0.51  0.73  0.86  0.89  0.92  1.03  1.00  1.02  1.00
     152     5.201  0.24  0.48  0.57  0.74  0.81  0.84  0.83  0.87  0.88  0.89
     131     4.060  0.06  0.20  0.35  0.44  0.50  0.54  0.58  0.57  0.59  0.60
      81     3.439  0.20  0.21  0.24  0.31  0.38  0.44  0.50  0.53  0.60  0.62
      59     3.037  0.07  0.10  0.08  0.07  0.23  0.30  0.34  0.40  0.44  0.47
      31     2.747  0.26  0.19  0.23  0.27  0.37  0.40  0.38  0.41  0.49  0.47
      20     2.528  0.00  0.00  0.00  0.00  0.00  0.18  0.17  0.16  0.23  0.29
      15     2.354  0.00  0.00  0.00  0.14  0.13  0.12  0.11  0.11  0.20  0.20
      16     2.211  0.00  0.00  0.00  0.13  0.12  0.11  0.10  0.10  0.19  0.27
      13     2.092  0.31  0.22  0.37  0.33  0.30  0.27  0.26  0.24  0.23  0.34
       9     1.990  0.45  0.32  0.27  0.23  0.21  0.20  0.19  0.18  0.17  0.16
       2     1.903  2.01  1.45  1.20  1.06  0.96  1.78  1.67  1.59  1.52  1.47
     671   overall  0.17  0.31  0.41  0.51  0.57  0.61  0.65  0.66  0.69  0.70


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

    3158     8.275  0.34  0.62  0.78  0.90  0.99  1.03  1.04  1.05  1.06  1.06
    5767     5.201  0.16  0.39  0.61  0.73  0.84  0.91  0.95  0.98  1.01  1.02
    7471     4.060  0.09  0.29  0.49  0.66  0.78  0.86  0.93  0.97  1.00  1.02
    8503     3.439  0.15  0.34  0.54  0.73  0.85  0.94  0.99  1.03  1.04  1.06
    9160     3.037  0.42  0.54  0.70  0.85  0.95  1.02  1.06  1.07  1.08  1.08
    9288     2.747  0.88  0.75  0.76  0.84  0.92  0.98  1.01  1.05  1.07  1.08
    9309     2.528  1.56  1.09  0.97  0.95  0.97  0.99  1.01  1.03  1.04  1.05
    9347     2.354  2.23  1.38  1.13  1.04  0.99  0.97  0.96  0.96  0.96  0.97
    9529     2.211  2.82  1.66  1.28  1.11  1.02  0.96  0.94  0.93  0.92  0.92
    9765     2.092  3.35  1.88  1.41  1.19  1.06  0.99  0.93  0.90  0.88  0.86
    9714     1.990  3.69  2.04  1.50  1.23  1.07  0.98  0.92  0.87  0.84  0.82
    8851     1.903  4.10  2.23  1.61  1.30  1.13  1.02  0.93  0.88  0.83  0.80
    6702     1.826  4.19  2.26  1.63  1.31  1.13  1.02  0.93  0.87  0.83  0.80
    5183     1.757  4.38  2.35  1.67  1.34  1.14  1.02  0.93  0.87  0.82  0.79
    4093     1.695  4.39  2.36  1.69  1.37  1.17  1.04  0.95  0.89  0.84  0.81
    3120     1.639  4.52  2.42  1.72  1.37  1.17  1.04  0.94  0.88  0.83  0.79
    2258     1.589  4.60  2.47  1.75  1.39  1.18  1.04  0.95  0.88  0.83  0.79
    1254     1.545  4.32  2.32  1.66  1.32  1.13  1.03  0.93  0.87  0.82  0.79
     268     1.507  5.21  2.76  1.96  1.55  1.34  1.17  1.05  0.97  0.91  0.86
  122740   overall  2.31  1.40  1.15  1.05  1.01  0.98  0.96  0.95  0.94  0.94


 List of    500 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

    9   92    5    1.92   28.74  0.9559E+02  0.1972E+02 "alien"
   29    2   -1    2.12   26.22  0.1830E+03  0.2402E+02 "alien"
   24    3   -5    2.58   26.00  0.4854E+03  0.8521E+02 "alien"
   24    3   -4    2.58   25.41  0.4743E+03  0.8315E+02 "alien"
   24    2   -6    2.58   24.86  0.4642E+03  0.8165E+02 "alien"
   28   23   -5    2.14   24.57  0.1714E+03  0.2426E+02 "alien"
   16   93   17    1.59   23.62  0.5494E+02  0.2369E+02 "alien"
   20   99   12    1.51   23.57  0.4957E+02  0.2215E+02 "alien"
   28   22   -7    2.14   23.26  0.1623E+03  0.2358E+02 "alien"
   24    2   -4    2.58   22.96  0.4286E+03  0.7555E+02 "alien"
    8   94   16    1.74   22.13  0.5074E+02  0.2342E+02 "alien"
   14   94   19    1.58   21.67  0.5039E+02  0.2280E+02 "alien"
   24   88   14    1.50   21.31  0.4481E+02  0.2007E+02 "alien"
   21   63   19    1.73   21.00  0.4816E+02  0.2091E+02 "alien"
   16   88   20    1.59   20.84  0.4846E+02  0.2423E+02 "alien"
   34   67  -12    1.52   20.82  0.4379E+02  0.2214E+02 "alien"
   29    1   -1    2.12   20.64  0.1440E+03  0.1946E+02 "alien"
   12   98    1    1.78   20.38  0.4672E+02  0.2059E+02 "alien"
   26    1   -1    2.36   20.33  0.2690E+03  0.4788E+02 "alien"
   33   33   -7    1.78   20.16  0.4622E+02  0.1936E+02 "alien"
   31   22   -3    1.94   20.13  0.6695E+02  0.1625E+02 "alien"
   31   63   -7    1.66   19.96  0.5164E+02  0.2158E+02 "alien"
   22   54   20    1.75   19.58  0.4489E+02  0.1932E+02 "alien"
   17   99    1    1.67   19.51  0.5165E+02  0.2105E+02 "alien"
   13   91   21    1.60   19.31  0.4491E+02  0.2304E+02 "alien"
   23   72   19    1.59   19.15  0.4455E+02  0.2263E+02 "alien"
   16   87   20    1.60   19.10  0.4442E+02  0.2315E+02 "alien"
   16   92   22    1.52   19.00  0.3996E+02  0.2460E+02 "alien"
   10   96   10    1.77   18.96  0.4348E+02  0.2039E+02 "alien"
   26   76   12    1.58   18.86  0.4387E+02  0.1865E+02 "alien"
    8   88   22    1.69   18.82  0.4981E+02  0.2401E+02 "alien"
   18   90   17    1.58   18.78  0.4369E+02  0.2104E+02 "alien"
   29   46   -4    1.89   18.64  0.6202E+02  0.2039E+02 "alien"
   15   71   27    1.62   18.58  0.4809E+02  0.2294E+02 "alien"
   32   63   -4    1.62   18.58  0.4808E+02  0.2086E+02 "alien"
   11   86   20    1.71   18.44  0.4880E+02  0.2268E+02 "alien"
   15   75   24    1.66   18.17  0.4701E+02  0.2285E+02 "alien"
   14   88   23    1.57   18.00  0.4187E+02  0.2357E+02 "alien"
   15   88   19    1.63   17.94  0.4643E+02  0.2351E+02 "alien"
   24   10   -5    2.56   17.94  0.3350E+03  0.5945E+02 "alien"
   13  100    9    1.67   17.88  0.4734E+02  0.2163E+02 "alien"
   26   89   -9    1.57   17.81  0.4143E+02  0.1767E+02 "alien"
   13   79   24    1.66   17.44  0.4513E+02  0.2322E+02 "alien"
   32   42   -6    1.78   17.35  0.3978E+02  0.1878E+02 "alien"
    8   96   13    1.76   17.32  0.3971E+02  0.2088E+02 "alien"
   33   74   -9    1.51   17.30  0.3639E+02  0.2157E+02 "alien"
   13   82   24    1.63   17.28  0.4472E+02  0.2294E+02 "alien"
   23   14  -10    2.60   17.26  0.3222E+03  0.5882E+02 "alien"
   29   20   -3    2.08   17.22  0.1202E+03  0.1817E+02 "alien"
   31   74   -8    1.57   17.21  0.4002E+02  0.2119E+02 "alien"
   24   73   17    1.59   17.05  0.3965E+02  0.1461E+02 "alien"
   24    3   -1    2.56   17.02  0.3178E+03  0.3984E+02 "alien"
   29   79   -3    1.58   17.02  0.3958E+02  0.1977E+02 "alien"
   32   60   -8    1.65   16.94  0.4382E+02  0.2105E+02 "alien"
   13   78   26    1.62   16.93  0.4381E+02  0.2363E+02 "alien"
   14   99   13    1.62   16.90  0.4373E+02  0.2157E+02 "alien"
   25   14    2    2.37   16.87  0.2232E+03  0.4010E+02 "alien"
    7   77   25    1.77   16.83  0.3860E+02  0.2334E+02 "alien"
   28   22   -5    2.14   16.82  0.1174E+03  0.1822E+02 "alien"
   33   61   -9    1.60   16.71  0.3887E+02  0.2091E+02 "alien"
   15   76   24    1.65   16.56  0.4284E+02  0.2292E+02 "alien"
   31   54   -8    1.73   16.53  0.3790E+02  0.2010E+02 "alien"
   16  100    6    1.64   16.46  0.4260E+02  0.2062E+02 "alien"
   29   16   -3    2.10   16.46  0.1148E+03  0.1737E+02 "alien"
   29    4    0    2.10   16.44  0.1147E+03  0.1580E+02 "alien"
   27   12   -3    2.27   16.39  0.1537E+03  0.2128E+02 "alien"
   28   11   -1    2.17   16.33  0.1532E+03  0.2067E+02 "alien"
   23   51   16    1.84   16.29  0.4995E+02  0.2032E+02 "alien"
    8   89   20    1.73   16.23  0.3722E+02  0.2247E+02 "alien"
   25   19   -6    2.41   16.21  0.2145E+03  0.2885E+02 "alien"
   27   84    9    1.52   16.20  0.3406E+02  0.1962E+02 "alien"
   11   85   23    1.66   16.18  0.4186E+02  0.2282E+02 "alien"
    8   91   17    1.76   16.15  0.3704E+02  0.2181E+02 "alien"
   18   99   -2    1.65   16.15  0.4179E+02  0.1977E+02 "alien"
   17   87  -15    1.77   16.14  0.3701E+02  0.1784E+02 "alien"
   21   73   21    1.59   16.09  0.3742E+02  0.2138E+02 "alien"
   34   58  -10    1.58   16.04  0.3731E+02  0.2118E+02 "alien"
   31   19   -1    1.94   15.98  0.5316E+02  0.1059E+02 "alien"
   15  100   12    1.60   15.96  0.3713E+02  0.2169E+02 "alien"
   10   97    9    1.77   15.96  0.3659E+02  0.2044E+02 "alien"
   33   34   -7    1.78   15.93  0.3652E+02  0.1813E+02 "alien"
   17   27   27    1.89   15.85  0.5273E+02  0.1809E+02 "alien"
   32    4   -3    1.93   15.83  0.5265E+02  0.1511E+02 "alien"
   23   16  -10    2.58   15.80  0.2950E+03  0.5420E+02 "alien"
   26   85   12    1.50   15.78  0.3319E+02  0.1888E+02 "alien"
   30   77   -4    1.57   15.77  0.3668E+02  0.1914E+02 "alien"
    8   87   22    1.71   15.76  0.4171E+02  0.2322E+02 "alien"
   28   32  -14    2.01   15.74  0.7841E+02  0.2071E+02 "alien"
   29   62  -19    1.66   15.65  0.4049E+02  0.1812E+02 "alien"
   19   64   27    1.58   15.64  0.3638E+02  0.1943E+02 "alien"
   29   85    3    1.51   15.59  0.3279E+02  0.2034E+02 "alien"
    1   87   20    1.83   15.59  0.4782E+02  0.2345E+02 "alien"
   32   77   -8    1.51   15.59  0.3278E+02  0.2035E+02 "alien"
   23   69   19    1.61   15.56  0.3620E+02  0.2143E+02 "alien"
   24   45   17    1.82   15.56  0.4773E+02  0.1982E+02 "alien"
   30   30  -12    1.94   15.56  0.5174E+02  0.1633E+02 "alien"
   26    8   -9    2.36   15.53  0.2055E+03  0.3927E+02 "alien"
   17   86   21    1.57   15.51  0.3608E+02  0.2273E+02 "alien"
   18    1   -4    3.45   15.51  0.1760E+04  0.1525E+03 "alien"
   14   51   31    1.69   15.51  0.4106E+02  0.1941E+02 "alien"
   25   85    5    1.60   15.50  0.3604E+02  0.1943E+02 "alien"
   31   27  -11    1.91   15.46  0.5142E+02  0.1639E+02 "alien"
    9   91   18    1.73   15.45  0.3542E+02  0.2171E+02 "alien"
   25   20   -4    2.40   15.44  0.2043E+03  0.3834E+02 "alien"
   22   89   12    1.57   15.43  0.3588E+02  0.1869E+02 "alien"
   31   66    0    1.62   15.32  0.3963E+02  0.1839E+02 "alien"
   19   11   -2    3.19   15.31  0.8101E+03  0.8192E+02 "alien"
   34   53  -10    1.62   15.27  0.3950E+02  0.2139E+02 "alien"
   32   77   -4    1.51   15.17  0.3191E+02  0.1921E+02 "alien"
    8   83   24    1.71   15.07  0.3988E+02  0.2239E+02 "alien"
   23   48   19    1.78   15.04  0.3448E+02  0.1882E+02 "alien"
   28   10    2    2.14   15.03  0.1049E+03  0.2068E+02 "alien"
   31   78   -8    1.54   15.01  0.5040E+02  0.2249E+02 "alien"
   33   33   -5    1.78   15.01  0.3440E+02  0.1740E+02 "alien"
   24   71  -18    1.75   15.00  0.3440E+02  0.1773E+02 "alien"
   31   54   -9    1.73   15.00  0.3439E+02  0.1293E+02 "alien"
   26   55  -22    1.78   14.95  0.3428E+02  0.1639E+02 "alien"
   24   32   17    1.91   14.95  0.4973E+02  0.1907E+02 "alien"
   34   70  -10    1.50   14.95  0.3143E+02  0.2115E+02 "alien"
   21   88   14    1.58   14.95  0.3476E+02  0.1918E+02 "alien"
   32   59  -19    1.59   14.91  0.3468E+02  0.1818E+02 "alien"
   22   52   21    1.74   14.91  0.3419E+02  0.1822E+02 "alien"
    8   85   21    1.76   14.89  0.3415E+02  0.2184E+02 "alien"
   26   62  -19    1.76   14.82  0.3397E+02  0.1658E+02 "alien"
   22   48   22    1.74   14.81  0.3397E+02  0.1750E+02 "alien"
   16   47   28    1.76   14.79  0.3391E+02  0.1719E+02 "alien"
   16   99    5    1.67   14.79  0.3826E+02  0.2073E+02 "alien"
   32    7   -3    1.92   14.77  0.4914E+02  0.1563E+02 "alien"
   24   62   13    1.76   14.77  0.3386E+02  0.1332E+02 "alien"
   18   98    1    1.66   14.77  0.3821E+02  0.1901E+02 "alien"
   28    9    1    2.16   14.76  0.1385E+03  0.2640E+02 "alien"
   15   98    5    1.70   14.74  0.3903E+02  0.1990E+02 "alien"
   27    1   -1    2.27   14.74  0.1383E+03  0.2646E+02 "alien"
   27   65  -16    1.74   14.74  0.3379E+02  0.1786E+02 "alien"
   29   56    6    1.72   14.72  0.3374E+02  0.1634E+02 "alien"
   26   87   12    1.49   14.70  0.3092E+02  0.1868E+02 "alien"
   29    8    0    2.10   14.70  0.1026E+03  0.1161E+02 "alien"
   28   25   -5    2.12   14.67  0.1023E+03  0.1681E+02 "alien"
    7    1   -1    8.84   14.66  0.3654E+04  0.3658E+03 "alien"
   28    8    1    2.16   14.61  0.1371E+03  0.2600E+02 "alien"
   16   48   34    1.57   14.60  0.3395E+02  0.1910E+02 "alien"
   28   55    5    1.79   14.56  0.3337E+02  0.1618E+02 "alien"
   33   23   -7    1.83   14.55  0.4462E+02  0.1754E+02 "alien"
   27   12   -4    2.27   14.54  0.1364E+03  0.2647E+02 "alien"
   12   46   31    1.77   14.52  0.3330E+02  0.1748E+02 "alien"
   29    3    0    2.10   14.51  0.1013E+03  0.1442E+02 "alien"
   11  100    4    1.75   14.51  0.3327E+02  0.1997E+02 "alien"
    8   92   20    1.69   14.48  0.3833E+02  0.2274E+02 "alien"
   29    6    0    2.10   14.46  0.1009E+03  0.1149E+02 "alien"
   30   52  -12    1.77   14.46  0.3316E+02  0.1616E+02 "alien"
   25   69  -15    1.77   14.44  0.3312E+02  0.1708E+02 "alien"
   17   97    7    1.65   14.42  0.3732E+02  0.2032E+02 "alien"
   16   65   28    1.62   14.41  0.3729E+02  0.1929E+02 "alien"
   29   53  -16    1.77   14.39  0.3300E+02  0.1680E+02 "alien"
   22   71   19    1.62   14.36  0.3716E+02  0.2100E+02 "alien"
   14   60   30    1.65   14.33  0.3707E+02  0.1899E+02 "alien"
   26   88   11    1.49   14.32  0.3012E+02  0.1882E+02 "alien"
   21   51   20    1.81   14.31  0.4389E+02  0.1826E+02 "alien"
   27   54  -23    1.73   14.27  0.3272E+02  0.1755E+02 "alien"
    9   88   23    1.66   14.23  0.3683E+02  0.2232E+02 "alien"
   25   12    2    2.38   14.22  0.1882E+03  0.3421E+02 "alien"
   20   60   21    1.74   14.20  0.3255E+02  0.1826E+02 "alien"
   24    8   -5    2.57   14.18  0.2648E+03  0.4788E+02 "alien"
   15   88   23    1.55   14.18  0.4761E+02  0.2335E+02 "alien"
   10   92   16    1.73   14.17  0.3250E+02  0.2245E+02 "alien"
    4   82   24    1.77   14.09  0.3230E+02  0.2182E+02 "alien"
   31   57  -20    1.62   14.08  0.3644E+02  0.1770E+02 "alien"
   31   43  -16    1.77   14.08  0.3227E+02  0.1557E+02 "alien"
   31    9   -3    1.98   14.08  0.7012E+02  0.1882E+02 "alien"
   28   67  -10    1.73   14.06  0.3223E+02  0.1645E+02 "alien"
    9   86   20    1.75   14.06  0.3223E+02  0.2202E+02 "alien"
   29   62  -14    1.71   14.03  0.3714E+02  0.1773E+02 "alien"
   21   98   -2    1.60   14.01  0.3260E+02  0.1878E+02 "alien"
   33   73  -11    1.51   14.01  0.2946E+02  0.2122E+02 "alien"
   23   63   16    1.73   14.00  0.3210E+02  0.1320E+02 "alien"
   22   73   21    1.57   14.00  0.3256E+02  0.2208E+02 "alien"
    3   84   24    1.76   13.98  0.3205E+02  0.2215E+02 "alien"
    7    1   -2    8.83   13.97  0.3484E+04  0.3001E+03 "alien"
   21   62   19    1.74   13.97  0.3204E+02  0.1944E+02 "alien"
   32   51   -8    1.71   13.97  0.3698E+02  0.1966E+02 "alien"
   13   83   24    1.62   13.97  0.3614E+02  0.2181E+02 "alien"
   10   93   18    1.68   13.97  0.3697E+02  0.2281E+02 "alien"
   10   93   14    1.76   13.96  0.3202E+02  0.2118E+02 "alien"
   24   14  -10    2.50   13.96  0.2607E+03  0.4892E+02 "alien"
   20   57   23    1.71   13.94  0.3690E+02  0.1922E+02 "alien"
   29   32  -10    2.00   13.92  0.6936E+02  0.1925E+02 "alien"
   10   97   -7    1.83   13.91  0.4267E+02  0.2134E+02 "alien"
   21   10   23    1.94   13.90  0.4623E+02  0.1525E+02 "alien"
   11   96   -6    1.83   13.89  0.4259E+02  0.2041E+02 "alien"
   33   66   -9    1.56   13.88  0.3229E+02  0.2097E+02 "alien"
   22   54   22    1.70   13.87  0.3671E+02  0.1805E+02 "alien"
   33   27   -7    1.81   13.83  0.4242E+02  0.1781E+02 "alien"
   28   65  -15    1.71   13.82  0.3658E+02  0.1748E+02 "alien"
   34   49  -10    1.64   13.81  0.3574E+02  0.2042E+02 "alien"
   32   20   -1    1.88   13.79  0.4586E+02  0.1504E+02 "alien"
   28    8    3    2.13   13.78  0.9615E+02  0.1801E+02 "alien"
   24    1   -4    2.58   13.77  0.2571E+03  0.4647E+02 "alien"
   31   13   -3    1.97   13.76  0.6856E+02  0.1247E+02 "alien"
   29   63  -17    1.68   13.75  0.3640E+02  0.1860E+02 "alien"
   29   18   -3    2.09   13.74  0.9589E+02  0.1550E+02 "alien"
   29    5    0    2.10   13.74  0.9587E+02  0.1376E+02 "alien"
   27   87    7    1.51   13.73  0.2887E+02  0.1927E+02 "alien"
   29    3   -1    2.11   13.71  0.9567E+02  0.1424E+02 "alien"
   10  100    3    1.78   13.70  0.3141E+02  0.2002E+02 "alien"
   12   92   18    1.66   13.68  0.3540E+02  0.2205E+02 "alien"
   30   69  -14    1.62   13.67  0.3538E+02  0.1710E+02 "alien"
   11   93   15    1.72   13.65  0.3612E+02  0.2049E+02 "alien"
   16   94   16    1.60   13.64  0.3173E+02  0.2106E+02 "alien"
   31   26    1    1.89   13.63  0.4533E+02  0.1516E+02 "alien"
   30   67    1    1.64   13.62  0.3524E+02  0.1844E+02 "alien"
   29   32   -2    2.00   13.62  0.6785E+02  0.1361E+02 "alien"
   28   68  -15    1.68   13.61  0.3603E+02  0.1759E+02 "alien"
   18   62   22    1.76   13.60  0.3118E+02  0.1981E+02 "alien"
   31    7   -3    1.99   13.60  0.6774E+02  0.1756E+02 "alien"
   33   61  -16    1.57   13.59  0.3162E+02  0.1832E+02 "alien"
   18   67   27    1.58   13.58  0.3158E+02  0.1937E+02 "alien"
   16   92   23    1.50   13.56  0.2851E+02  0.2375E+02 "alien"
   31   21    0    1.92   13.55  0.4508E+02  0.1550E+02 "alien"
   32   38   -1    1.79   13.55  0.3107E+02  0.1652E+02 "alien"
   13   92   20    1.61   13.54  0.3149E+02  0.2260E+02 "alien"
    5   94   14    1.81   13.54  0.4151E+02  0.2282E+02 "alien"
   25   78  -10    1.71   13.53  0.3582E+02  0.1368E+02 "alien"
   16   98   -3    1.71   13.53  0.3582E+02  0.1959E+02 "alien"
   15   73   24    1.68   13.52  0.3580E+02  0.2197E+02 "alien"
   29   54    7    1.73   13.51  0.3099E+02  0.1678E+02 "alien"
   21   59   19    1.76   13.51  0.3098E+02  0.1877E+02 "alien"
   15   94   -2    1.79   13.51  0.4143E+02  0.1911E+02 "alien"
   31   68  -17    1.57   13.51  0.3142E+02  0.1775E+02 "alien"
   12   99   -1    1.77   13.50  0.3094E+02  0.1962E+02 "alien"
    9   53   32    1.75   13.49  0.3092E+02  0.1803E+02 "alien"
   27    6   -8    2.28   13.48  0.1784E+03  0.3438E+02 "alien"
   18   46   29    1.68   13.47  0.3567E+02  0.1730E+02 "alien"
   32   50   -4    1.72   13.47  0.3567E+02  0.1828E+02 "alien"
   14   81   26    1.57   13.46  0.3131E+02  0.2296E+02 "alien"
   11   48   34    1.67   13.46  0.3562E+02  0.1897E+02 "alien"
   28   23   -6    2.14   13.44  0.9379E+02  0.1605E+02 "alien"
   11   86   22    1.67   13.44  0.3557E+02  0.2271E+02 "alien"
   19   67   27    1.56   13.43  0.4511E+02  0.2017E+02 "alien"
   15   66   22    1.80   13.42  0.4116E+02  0.1453E+02 "alien"
   21   82   17    1.59   13.41  0.3120E+02  0.1806E+02 "alien"
   31   68  -15    1.59   13.39  0.3115E+02  0.1728E+02 "alien"
   21   61   18    1.77   13.37  0.3066E+02  0.1901E+02 "alien"
   22   52   20    1.76   13.36  0.3062E+02  0.1823E+02 "alien"
   26   92   -2    1.54   13.35  0.4484E+02  0.2182E+02 "alien"
   17   40   29    1.74   13.34  0.3058E+02  0.1718E+02 "alien"
   30   67   -2    1.65   13.32  0.3447E+02  0.1252E+02 "alien"
   26   69    7    1.71   13.32  0.3525E+02  0.1343E+02 "alien"
   23   15  -10    2.59   13.31  0.2485E+03  0.4670E+02 "alien"
   26   17   -4    2.33   13.31  0.1761E+03  0.2405E+02 "alien"
   17   90   -9    1.79   13.31  0.3051E+02  0.1740E+02 "alien"
   13   98    8    1.71   13.30  0.3520E+02  0.2020E+02 "alien"
   28   10  -11    2.17   13.30  0.1247E+03  0.2661E+02 "alien"
   15   90   21    1.57   13.29  0.3092E+02  0.2244E+02 "alien"
   11   47   33    1.72   13.28  0.3516E+02  0.1838E+02 "alien"
   21   68   23    1.59   13.28  0.3089E+02  0.2030E+02 "alien"
    8   99    8    1.78   13.24  0.3037E+02  0.1987E+02 "alien"
   12   58   30    1.71   13.24  0.3506E+02  0.1867E+02 "alien"
   28   86    7    1.49   13.24  0.2784E+02  0.2004E+02 "alien"
   31   21   -1    1.93   13.24  0.4403E+02  0.1021E+02 "alien"
   32   57   -8    1.67   13.23  0.3502E+02  0.1981E+02 "alien"
   25   85   -8    1.64   13.22  0.3420E+02  0.1266E+02 "alien"
   32   46   -6    1.75   13.21  0.3028E+02  0.1794E+02 "alien"
   29    7    0    2.10   13.19  0.9200E+02  0.1069E+02 "alien"
   23   70   19    1.60   13.17  0.3064E+02  0.2106E+02 "alien"
    8   97    7    1.82   13.17  0.4038E+02  0.2060E+02 "alien"
   25   67   16    1.62   13.17  0.3406E+02  0.1733E+02 "alien"
   32   78   -8    1.51   13.16  0.2768E+02  0.2069E+02 "alien"
   13   83   26    1.57   13.15  0.3060E+02  0.2277E+02 "alien"
    5   96   14    1.78   13.15  0.3015E+02  0.2237E+02 "alien"
   27   25  -16    2.09   13.14  0.9164E+02  0.2258E+02 "alien"
   24    3   -6    2.58   13.13  0.2451E+03  0.4466E+02 "alien"
   31   11   -3    1.98   13.11  0.6532E+02  0.1215E+02 "alien"
   14   67   27    1.68   13.11  0.3470E+02  0.1903E+02 "alien"
   14   71   26    1.67   13.11  0.3470E+02  0.2197E+02 "alien"
   23   94   -4    1.60   13.10  0.3048E+02  0.1728E+02 "alien"
   21   80   17    1.61   13.10  0.3047E+02  0.1406E+02 "alien"
    5  100    3    1.84   13.10  0.4016E+02  0.2150E+02 "alien"
   27   14   -4    2.26   13.09  0.1229E+03  0.1805E+02 "alien"
   24   36   17    1.88   13.09  0.4352E+02  0.1893E+02 "alien"
   27   76    9    1.59   13.06  0.3037E+02  0.1418E+02 "alien"
   29   61  -20    1.66   13.04  0.3374E+02  0.1761E+02 "alien"
   31   80   -1    1.51   13.02  0.2737E+02  0.1854E+02 "alien"
   25   76   14    1.58   13.01  0.3026E+02  0.2363E+02 "alien"
   23   58   19    1.70   13.01  0.3443E+02  0.1965E+02 "alien"
    0   77   27    1.79   13.00  0.2981E+02  0.2083E+02 "alien"
   16   63   28    1.64   13.00  0.3363E+02  0.1806E+02 "alien"
   17   48   32    1.60   12.99  0.3022E+02  0.1870E+02 "alien"
   17   90   19    1.57   12.98  0.3019E+02  0.2159E+02 "alien"
   32   41   -4    1.78   12.98  0.2976E+02  0.1727E+02 "alien"
   28   14    2    2.13   12.98  0.9056E+02  0.1905E+02 "alien"
   17   66   22    1.75   12.96  0.2972E+02  0.1911E+02 "alien"
   14   29   31    1.82   12.96  0.3974E+02  0.1586E+02 "alien"
    7   79   24    1.77   12.94  0.2967E+02  0.2200E+02 "alien"
   27    5   -3    2.29   12.94  0.1712E+03  0.3167E+02 "alien"
   13   85   22    1.64   12.94  0.3347E+02  0.2246E+02 "alien"
   32   66   -5    1.60   12.93  0.3008E+02  0.1980E+02 "alien"
   15   52   33    1.60   12.93  0.3007E+02  0.1943E+02 "alien"
   22   61   17    1.76   12.88  0.2954E+02  0.1947E+02 "alien"
   29   87   -9    1.51   12.86  0.2705E+02  0.1655E+02 "alien"
   23   17    0    2.58   12.85  0.2399E+03  0.3037E+02 "alien"
   30   32    1    1.91   12.84  0.4272E+02  0.1539E+02 "alien"
   31   23   -5    1.94   12.84  0.4271E+02  0.1115E+02 "alien"
   23   55   20    1.70   12.83  0.3396E+02  0.1919E+02 "alien"
   30   45  -14    1.81   12.81  0.3929E+02  0.1620E+02 "alien"
   26   17   -8    2.32   12.81  0.1695E+03  0.3379E+02 "alien"
   34   68  -13    1.51   12.80  0.2692E+02  0.2172E+02 "alien"
   21   12  -14    2.70   12.80  0.3599E+03  0.6518E+02 "alien"
   25   63   16    1.66   12.79  0.3311E+02  0.1777E+02 "alien"
   19   90    1    1.74   12.79  0.2933E+02  0.1104E+02 "alien"
   18   88   18    1.58   12.79  0.2975E+02  0.2122E+02 "alien"
   20   47   26    1.70   12.79  0.3384E+02  0.1714E+02 "alien"
   16   53   32    1.60   12.78  0.2973E+02  0.1818E+02 "alien"
   31   20   -1    1.94   12.76  0.4244E+02  0.1011E+02 "alien"
   31   19   -3    1.95   12.76  0.6356E+02  0.1273E+02 "alien"
   15   95   -2    1.78   12.75  0.2924E+02  0.1784E+02 "alien"
   21   84   -7    1.77   12.75  0.2923E+02  0.1345E+02 "alien"
   24   87    5    1.61   12.73  0.2961E+02  0.1348E+02 "alien"
   11   93   17    1.69   12.73  0.3369E+02  0.2243E+02 "alien"
   28   52  -20    1.77   12.73  0.2918E+02  0.1628E+02 "alien"
   32   65  -12    1.60   12.72  0.2959E+02  0.2118E+02 "alien"
   16   83   22    1.60   12.72  0.2957E+02  0.2158E+02 "alien"
   24   19  -12    2.43   12.70  0.1680E+03  0.3413E+02 "alien"
   14   75   23    1.70   12.69  0.3360E+02  0.2200E+02 "alien"
   12   47   32    1.73   12.69  0.2910E+02  0.1781E+02 "alien"
   31   21   -3    1.94   12.67  0.4214E+02  0.1321E+02 "alien"
   10   89   20    1.69   12.66  0.3352E+02  0.2190E+02 "alien"
   13   89   22    1.60   12.66  0.2945E+02  0.2258E+02 "alien"
   26   58  -15    1.86   12.65  0.3880E+02  0.1782E+02 "alien"
   29    0   -1    2.12   12.65  0.2116E+03  0.3843E+02 "alien"
   16   94   20    1.53   12.62  0.4236E+02  0.2244E+02 "alien"
   25   88  -10    1.60   12.61  0.2932E+02  0.1680E+02 "alien"
   18   61   25    1.68   12.57  0.3329E+02  0.1868E+02 "alien"
   13   94   12    1.72   12.55  0.3323E+02  0.1352E+02 "alien"
   15   67   23    1.76   12.55  0.2877E+02  0.1379E+02 "alien"
   19   99   14    1.51   12.55  0.2639E+02  0.2078E+02 "alien"
   32   51   -4    1.71   12.55  0.3321E+02  0.1815E+02 "alien"
   13   82   19    1.74   12.53  0.2874E+02  0.1776E+02 "alien"
   23   61   18    1.70   12.53  0.3317E+02  0.1989E+02 "alien"
   19    4  -10    3.14   12.53  0.6627E+03  0.8300E+02 "alien"
   16   90   18    1.61   12.50  0.2908E+02  0.2094E+02 "alien"
    0    4    0   46.73   12.50  0.6234E+04  0.1071E+04 "alien"
   24    8   -7    2.56   12.50  0.2334E+03  0.4327E+02 "alien"
   24   74   16    1.60   12.48  0.2903E+02  0.1393E+02 "alien"
   31   57  -16    1.66   12.47  0.3227E+02  0.1745E+02 "alien"
   32   65   -7    1.61   12.47  0.2900E+02  0.1968E+02 "alien"
   12   81   23    1.68   12.47  0.3300E+02  0.2213E+02 "alien"
   30   74  -15    1.57   12.46  0.2899E+02  0.1733E+02 "alien"
   17   79   23    1.60   12.46  0.2898E+02  0.2179E+02 "alien"
   30   70  -13    1.62   12.45  0.3222E+02  0.1372E+02 "alien"
   21   11    5    2.72   12.43  0.3497E+03  0.6080E+02 "alien"
   21   72   21    1.60   12.42  0.2888E+02  0.2000E+02 "alien"
   25   73   -9    1.77   12.41  0.2845E+02  0.1336E+02 "alien"
   29   54  -20    1.71   12.41  0.3285E+02  0.1677E+02 "alien"
   20    2    0    3.05   12.40  0.6560E+03  0.8059E+02 "alien"
   28   10    1    2.15   12.40  0.1163E+03  0.2316E+02 "alien"
   32   28   -2    1.85   12.39  0.3800E+02  0.1613E+02 "alien"
    6   94   16    1.77   12.39  0.2840E+02  0.2119E+02 "alien"
   15   75   25    1.63   12.38  0.3203E+02  0.2282E+02 "alien"
   10   98   13    1.70   12.38  0.3276E+02  0.2191E+02 "alien"
   14   96    9    1.71   12.36  0.3272E+02  0.1985E+02 "alien"
   30   58  -11    1.73   12.36  0.2834E+02  0.1358E+02 "alien"
   21   54   23    1.70   12.36  0.3271E+02  0.1838E+02 "alien"
    1   89   20    1.80   12.36  0.3789E+02  0.2284E+02 "alien"
   23   12   -9    2.62   12.35  0.2305E+03  0.4329E+02 "alien"
   25   74   -6    1.77   12.33  0.2827E+02  0.8952E+01 "alien"
    7   96   15    1.74   12.32  0.2825E+02  0.2232E+02 "alien"
   19   74   22    1.62   12.31  0.3185E+02  0.2085E+02 "alien"
    6   87   23    1.71   12.30  0.3257E+02  0.2173E+02 "alien"
   15   49   34    1.58   12.29  0.2860E+02  0.1892E+02 "alien"
   11   96   17    1.65   12.29  0.3180E+02  0.2251E+02 "alien"
   20   69   19    1.70   12.29  0.3253E+02  0.2011E+02 "alien"
   13   87   23    1.60   12.28  0.2856E+02  0.2206E+02 "alien"
   33   54  -14    1.63   12.28  0.3177E+02  0.1727E+02 "alien"
   33   73  -10    1.51   12.27  0.2581E+02  0.2101E+02 "alien"
   11   95   12    1.74   12.26  0.2811E+02  0.1979E+02 "alien"
   27   10   -3    2.27   12.25  0.1149E+03  0.1651E+02 "alien"
   31    8   -3    1.99   12.24  0.6100E+02  0.1736E+02 "alien"
   32   52  -14    1.68   12.24  0.3239E+02  0.1703E+02 "alien"
   12   16   -1    4.71   12.22  0.2067E+04  0.2076E+03 "alien"
   15   96   10    1.68   12.21  0.3231E+02  0.2038E+02 "alien"
   25   71   12    1.66   12.21  0.3158E+02  0.1721E+02 "alien"
   23   11   -9    2.63   12.19  0.2277E+03  0.4280E+02 "alien"
   31   48   -5    1.78   12.19  0.2795E+02  0.1183E+02 "alien"
   20   67   21    1.67   12.19  0.3227E+02  0.2039E+02 "alien"
   28   70  -11    1.69   12.18  0.3223E+02  0.1348E+02 "alien"
   30   17   -4    2.03   12.17  0.6064E+02  0.1263E+02 "alien"
   19    6   -7    3.22   12.17  0.6437E+03  0.7962E+02 "alien"
   26   86   10    1.52   12.16  0.2557E+02  0.2008E+02 "alien"
   29   21   -9    2.07   12.15  0.8478E+02  0.1544E+02 "alien"
   30   21   -4    2.01   12.13  0.6044E+02  0.2005E+02 "alien"
   28   70    7    1.63   12.12  0.3137E+02  0.1809E+02 "alien"
    7   90   19    1.75   12.12  0.2779E+02  0.2170E+02 "alien"
    9   56   31    1.76   12.12  0.2779E+02  0.1761E+02 "alien"
   21   64   22    1.65   12.10  0.3132E+02  0.2033E+02 "alien"
   31   64    0    1.63   12.10  0.3131E+02  0.1756E+02 "alien"
   11   80   25    1.66   12.09  0.3129E+02  0.2206E+02 "alien"
   13   76   26    1.64   12.08  0.3126E+02  0.2141E+02 "alien"
   29   80   -2    1.57   12.08  0.2809E+02  0.1875E+02 "alien"
   18   61   23    1.74   12.07  0.2768E+02  0.1950E+02 "alien"
   32   72   -5    1.55   12.07  0.4053E+02  0.2132E+02 "alien"
   16   93   22    1.51   12.06  0.2536E+02  0.2330E+02 "alien"
   17   82   14    1.75   12.06  0.2765E+02  0.1116E+02 "alien"
   27   28  -13    2.12   12.06  0.8411E+02  0.2105E+02 "alien"
   30   63  -16    1.65   12.05  0.3118E+02  0.1762E+02 "alien"
    6   80   27    1.69   12.05  0.3190E+02  0.2274E+02 "alien"
   24   71  -16    1.77   12.05  0.2762E+02  0.1632E+02 "alien"
   32   57  -20    1.59   12.04  0.2801E+02  0.1783E+02 "alien"
   14   76   23    1.69   12.04  0.3186E+02  0.2152E+02 "alien"
   29   30  -10    2.01   12.04  0.5996E+02  0.1779E+02 "alien"
   31   27    1    1.89   12.03  0.4001E+02  0.1481E+02 "alien"
   24   81   11    1.61   12.03  0.2798E+02  0.1409E+02 "alien"
   19   75   22    1.61   12.01  0.2794E+02  0.2088E+02 "alien"
   30   34   -6    1.94   12.01  0.3994E+02  0.1532E+02 "alien"
   13   98    1    1.76   12.01  0.2753E+02  0.1966E+02 "alien"
   12   49   31    1.75   12.00  0.2752E+02  0.1726E+02 "alien"
   14   87   23    1.58   12.00  0.2791E+02  0.2263E+02 "alien"
   25   56  -19    1.86   12.00  0.3679E+02  0.1689E+02 "alien"
   29   34   -6    1.99   11.98  0.5969E+02  0.1376E+02 "alien"
   26   64   14    1.65   11.97  0.3098E+02  0.1702E+02 "alien"
   24   89   10    1.54   11.97  0.4018E+02  0.2064E+02 "alien"
   11   55   30    1.76   11.96  0.2743E+02  0.1808E+02 "alien"
   16   85   19    1.64   11.96  0.3095E+02  0.1793E+02 "alien"
   14   75   27    1.61   11.95  0.2779E+02  0.2143E+02 "alien"
   32   63  -11    1.61   11.95  0.3091E+02  0.2098E+02 "alien"
    0   87   20    1.84   11.94  0.3661E+02  0.2222E+02 "alien"
    8   68   24    1.89   11.93  0.3970E+02  0.1411E+02 "alien"
   29   68  -16    1.64   11.92  0.3083E+02  0.1737E+02 "alien"
   16   99    4    1.67   11.91  0.3152E+02  0.2008E+02 "alien"
    7   55   28    1.92   11.90  0.3960E+02  0.1845E+02 "alien"
   16   99   14    1.57   11.90  0.2769E+02  0.2058E+02 "alien"
   24   72  -16    1.76   11.90  0.2729E+02  0.1675E+02 "alien"
   32   60   -3    1.64   11.90  0.3079E+02  0.1886E+02 "alien"
   18   58   30    1.57   11.89  0.2767E+02  0.1852E+02 "alien"
    7   93   17    1.75   11.88  0.2725E+02  0.2207E+02 "alien"
   21   47   23    1.75   11.88  0.2724E+02  0.1664E+02 "alien"
   25   46  -21    1.92   11.88  0.3951E+02  0.1727E+02 "alien"
   13  101    5    1.70   11.87  0.3141E+02  0.1965E+02 "alien"
   31   43    1    1.79   11.86  0.2718E+02  0.1511E+02 "alien"
    7   98   -3    1.86   11.85  0.3940E+02  0.2170E+02 "alien"
   30   82    0    1.52   11.85  0.2491E+02  0.1976E+02 "alien"
   31   22   -4    1.94   11.84  0.3938E+02  0.1780E+02 "alien"
   13   96   -1    1.80   11.84  0.3631E+02  0.1957E+02 "alien"
   21  100    1    1.57   11.84  0.2753E+02  0.1930E+02 "alien"
   17   75   24    1.61   11.83  0.2750E+02  0.2187E+02 "alien"
   23   56   19    1.72   11.82  0.3130E+02  0.1953E+02 "alien"
   27    4    2    2.23   11.82  0.1109E+03  0.2057E+02 "alien"
   30   77   -8    1.57   11.82  0.2749E+02  0.2082E+02 "alien"
   31   63   -5    1.66   11.82  0.3057E+02  0.1903E+02 "alien"
    5   78  -22    1.93   11.81  0.3930E+02  0.1978E+02 "alien"
   27   64  -18    1.72   11.81  0.3127E+02  0.1825E+02 "alien"
   24   90   -9    1.61   11.80  0.2744E+02  0.1605E+02 "alien"
   13   82   23    1.65   11.79  0.3052E+02  0.2223E+02 "alien"
   22   29   21    1.90   11.78  0.3919E+02  0.1616E+02 "alien"
   28   25   -9    2.11   11.77  0.8214E+02  0.2023E+02 "alien"
   24   81   -6    1.72   11.75  0.3110E+02  0.1061E+02 "alien"
   14   93   18    1.61   11.74  0.3037E+02  0.2076E+02 "alien"
   16   61   24    1.77   11.73  0.2690E+02  0.1386E+02 "alien"
   28   74   -1    1.66   11.70  0.3027E+02  0.1324E+02 "alien"
   23   85   13    1.57   11.69  0.2718E+02  0.1844E+02 "alien"
   24   59   13    1.79   11.68  0.2678E+02  0.1285E+02 "alien"
   17   32   32    1.68   11.68  0.3091E+02  0.1640E+02 "alien"
   26   73   12    1.61   11.68  0.2716E+02  0.1414E+02 "alien"
   15   65   28    1.64   11.67  0.3020E+02  0.1868E+02 "alien"
   32    3   -3    1.93   11.67  0.3881E+02  0.1361E+02 "alien"
   31   70  -16    1.57   11.67  0.2714E+02  0.1729E+02 "alien"
   19   11   -6    3.19   11.67  0.6172E+03  0.5496E+02 "alien"
   20   30   27    1.77   11.67  0.2675E+02  0.1609E+02 "alien"
   27    7   -9    2.28   11.66  0.1094E+03  0.2342E+02 "alien"
   27   72    8    1.64   11.66  0.3017E+02  0.1355E+02 "alien"
   30   51    4    1.74   11.66  0.2673E+02  0.1611E+02 "alien"
   19   97   14    1.53   11.66  0.3914E+02  0.2089E+02 "alien"
    8   97   10    1.79   11.65  0.2671E+02  0.2027E+02 "alien"
   30   27    4    1.91   11.65  0.3874E+02  0.1350E+02 "alien"
   32   43  -15    1.73   11.64  0.2668E+02  0.1510E+02 "alien"
    5   88   21    1.76   11.62  0.2665E+02  0.2163E+02 "alien"
   16   71   25    1.65   11.62  0.3006E+02  0.2178E+02 "alien"
   32   30   -7    1.85   11.61  0.3561E+02  0.1382E+02 "alien"
   14   99    8    1.68   11.61  0.3073E+02  0.1994E+02 "alien"
    8   99    0    1.83   11.60  0.3559E+02  0.1990E+02 "alien"
   24   48   14    1.87   11.60  0.3858E+02  0.1940E+02 "alien"
   32   63   -9    1.62   11.60  0.3001E+02  0.2059E+02 "alien"
   22   56   15    1.85   11.60  0.3557E+02  0.1056E+02 "alien"
   24    1   -6    2.58   11.59  0.2164E+03  0.3964E+02 "alien"
   17   37   27    1.83   11.58  0.3551E+02  0.1642E+02 "alien"
   11   57   25    1.92   11.57  0.3849E+02  0.1888E+02 "alien"
   17   87   20    1.58   11.57  0.2691E+02  0.2081E+02 "alien"
   31   16   -3    1.96   11.57  0.5763E+02  0.1171E+02 "alien"
   31   82   -3    1.50   11.57  0.2432E+02  0.1942E+02 "alien"
   25   20   -8    2.39   11.56  0.1530E+03  0.2184E+02 "alien"
   28   19   -5    2.16   11.56  0.1085E+03  0.1713E+02 "alien"
   24   11   -5    2.56   11.55  0.2157E+03  0.3965E+02 "alien"
   25   22   -5    2.38   11.55  0.1528E+03  0.2128E+02 "alien"
   28   70  -16    1.66   11.54  0.2986E+02  0.1753E+02 "alien"
   16   56   26    1.75   11.53  0.2643E+02  0.1811E+02 "alien"
   23   61   20    1.66   11.52  0.2981E+02  0.1982E+02 "alien"
   25   72   12    1.65   11.51  0.2979E+02  0.1386E+02 "alien"
   19    6  -10    3.13   11.51  0.6089E+03  0.1071E+03 "alien"
   15   73   25    1.65   11.51  0.2977E+02  0.2198E+02 "alien"
   16   13  -17    3.04   11.50  0.6086E+03  0.6531E+02 "alien"
   32   54  -13    1.67   11.50  0.3045E+02  0.1645E+02 "alien"


 List of    500 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   26   87   12    1.49   14.70  0.3092E+02  0.1868E+02
   26   88   11    1.49   14.32  0.3012E+02  0.1882E+02
   28   86    7    1.49   13.24  0.2784E+02  0.2004E+02
   16   92   23    1.50   13.56  0.2851E+02  0.2375E+02
   31   82   -3    1.50   11.57  0.2432E+02  0.1942E+02
   34   70  -10    1.50   14.95  0.3143E+02  0.2115E+02
   24   88   14    1.50   21.31  0.4481E+02  0.2007E+02
   26   85   12    1.50   15.78  0.3319E+02  0.1888E+02
   29   85    3    1.51   15.59  0.3279E+02  0.2034E+02
   33   74   -9    1.51   17.30  0.3639E+02  0.2157E+02
   32   78   -8    1.51   13.16  0.2768E+02  0.2069E+02
   16   93   22    1.51   12.06  0.2536E+02  0.2330E+02
   34   68  -13    1.51   12.80  0.2692E+02  0.2172E+02
   33   73  -11    1.51   14.01  0.2946E+02  0.2122E+02
   31   80   -1    1.51   13.02  0.2737E+02  0.1854E+02
   19   99   14    1.51   12.55  0.2639E+02  0.2078E+02
   33   73  -10    1.51   12.27  0.2581E+02  0.2101E+02
   27   87    7    1.51   13.73  0.2887E+02  0.1927E+02
   32   77   -4    1.51   15.17  0.3191E+02  0.1921E+02
   29   87   -9    1.51   12.86  0.2705E+02  0.1655E+02
   32   77   -8    1.51   15.59  0.3278E+02  0.2035E+02
   20   99   12    1.51   23.57  0.4957E+02  0.2215E+02
   16   92   22    1.52   19.00  0.3996E+02  0.2460E+02
   26   86   10    1.52   12.16  0.2557E+02  0.2008E+02
   30   82    0    1.52   11.85  0.2491E+02  0.1976E+02
   34   67  -12    1.52   20.82  0.4379E+02  0.2214E+02
   27   84    9    1.52   16.20  0.3406E+02  0.1962E+02
   16   94   20    1.53   12.62  0.4236E+02  0.2244E+02
   19   97   14    1.53   11.66  0.3914E+02  0.2089E+02
   31   78   -8    1.54   15.01  0.5040E+02  0.2249E+02
   24   89   10    1.54   11.97  0.4018E+02  0.2064E+02
   26   92   -2    1.54   13.35  0.4484E+02  0.2182E+02
   32   72   -5    1.55   12.07  0.4053E+02  0.2132E+02
   15   88   23    1.55   14.18  0.4761E+02  0.2335E+02
   19   67   27    1.56   13.43  0.4511E+02  0.2017E+02
   33   66   -9    1.56   13.88  0.3229E+02  0.2097E+02
   31   70  -16    1.57   11.67  0.2714E+02  0.1729E+02
   31   74   -8    1.57   17.21  0.4002E+02  0.2119E+02
   22   73   21    1.57   14.00  0.3256E+02  0.2208E+02
   17   90   19    1.57   12.98  0.3019E+02  0.2159E+02
   16   48   34    1.57   14.60  0.3395E+02  0.1910E+02
   21  100    1    1.57   11.84  0.2753E+02  0.1930E+02
   22   89   12    1.57   15.43  0.3588E+02  0.1869E+02
   23   85   13    1.57   11.69  0.2718E+02  0.1844E+02
   33   61  -16    1.57   13.59  0.3162E+02  0.1832E+02
   16   99   14    1.57   11.90  0.2769E+02  0.2058E+02
   30   74  -15    1.57   12.46  0.2899E+02  0.1733E+02
   17   86   21    1.57   15.51  0.3608E+02  0.2273E+02
   26   89   -9    1.57   17.81  0.4143E+02  0.1767E+02
   29   80   -2    1.57   12.08  0.2809E+02  0.1875E+02
   15   90   21    1.57   13.29  0.3092E+02  0.2244E+02
   14   88   23    1.57   18.00  0.4187E+02  0.2357E+02
   18   58   30    1.57   11.89  0.2767E+02  0.1852E+02
   30   77   -4    1.57   15.77  0.3668E+02  0.1914E+02
   31   68  -17    1.57   13.51  0.3142E+02  0.1775E+02
   30   77   -8    1.57   11.82  0.2749E+02  0.2082E+02
   13   83   26    1.57   13.15  0.3060E+02  0.2277E+02
   14   81   26    1.57   13.46  0.3131E+02  0.2296E+02
   21   88   14    1.58   14.95  0.3476E+02  0.1918E+02
   17   87   20    1.58   11.57  0.2691E+02  0.2081E+02
   18   90   17    1.58   18.78  0.4369E+02  0.2104E+02
   26   76   12    1.58   18.86  0.4387E+02  0.1865E+02
   25   76   14    1.58   13.01  0.3026E+02  0.2363E+02
   14   87   23    1.58   12.00  0.2791E+02  0.2263E+02
   19   64   27    1.58   15.64  0.3638E+02  0.1943E+02
   29   79   -3    1.58   17.02  0.3958E+02  0.1977E+02
   18   88   18    1.58   12.79  0.2975E+02  0.2122E+02
   18   67   27    1.58   13.58  0.3158E+02  0.1937E+02
   15   49   34    1.58   12.29  0.2860E+02  0.1892E+02
   34   58  -10    1.58   16.04  0.3731E+02  0.2118E+02
   14   94   19    1.58   21.67  0.5039E+02  0.2280E+02
   24   73   17    1.59   17.05  0.3965E+02  0.1461E+02
   23   72   19    1.59   19.15  0.4455E+02  0.2263E+02
   21   82   17    1.59   13.41  0.3120E+02  0.1806E+02
   31   68  -15    1.59   13.39  0.3115E+02  0.1728E+02
   27   76    9    1.59   13.06  0.3037E+02  0.1418E+02
   32   59  -19    1.59   14.91  0.3468E+02  0.1818E+02
   16   93   17    1.59   23.62  0.5494E+02  0.2369E+02
   16   88   20    1.59   20.84  0.4846E+02  0.2423E+02
   32   57  -20    1.59   12.04  0.2801E+02  0.1783E+02
   21   68   23    1.59   13.28  0.3089E+02  0.2030E+02
   21   73   21    1.59   16.09  0.3742E+02  0.2138E+02
   24   74   16    1.60   12.48  0.2903E+02  0.1393E+02
   16   94   16    1.60   13.64  0.3173E+02  0.2106E+02
   32   65  -12    1.60   12.72  0.2959E+02  0.2118E+02
   15   52   33    1.60   12.93  0.3007E+02  0.1943E+02
   16   53   32    1.60   12.78  0.2973E+02  0.1818E+02
   17   79   23    1.60   12.46  0.2898E+02  0.2179E+02
   13   91   21    1.60   19.31  0.4491E+02  0.2304E+02
   32   66   -5    1.60   12.93  0.3008E+02  0.1980E+02
   13   89   22    1.60   12.66  0.2945E+02  0.2258E+02
   13   87   23    1.60   12.28  0.2856E+02  0.2206E+02
   21   98   -2    1.60   14.01  0.3260E+02  0.1878E+02
   23   94   -4    1.60   13.10  0.3048E+02  0.1728E+02
   33   61   -9    1.60   16.71  0.3887E+02  0.2091E+02
   16   87   20    1.60   19.10  0.4442E+02  0.2315E+02
   16   83   22    1.60   12.72  0.2957E+02  0.2158E+02
   15  100   12    1.60   15.96  0.3713E+02  0.2169E+02
   17   48   32    1.60   12.99  0.3022E+02  0.1870E+02
   25   85    5    1.60   15.50  0.3604E+02  0.1943E+02
   21   72   21    1.60   12.42  0.2888E+02  0.2000E+02
   23   70   19    1.60   13.17  0.3064E+02  0.2106E+02
   25   88  -10    1.60   12.61  0.2932E+02  0.1680E+02
   14   75   27    1.61   11.95  0.2779E+02  0.2143E+02
   16   90   18    1.61   12.50  0.2908E+02  0.2094E+02
   24   81   11    1.61   12.03  0.2798E+02  0.1409E+02
   13   92   20    1.61   13.54  0.3149E+02  0.2260E+02
   26   73   12    1.61   11.68  0.2716E+02  0.1414E+02
   21   80   17    1.61   13.10  0.3047E+02  0.1406E+02
   19   75   22    1.61   12.01  0.2794E+02  0.2088E+02
   32   65   -7    1.61   12.47  0.2900E+02  0.1968E+02
   24   87    5    1.61   12.73  0.2961E+02  0.1348E+02
   17   75   24    1.61   11.83  0.2750E+02  0.2187E+02
   23   69   19    1.61   15.56  0.3620E+02  0.2143E+02
   24   90   -9    1.61   11.80  0.2744E+02  0.1605E+02
   14   93   18    1.61   11.74  0.3037E+02  0.2076E+02
   32   63  -11    1.61   11.95  0.3091E+02  0.2098E+02
   30   70  -13    1.62   12.45  0.3222E+02  0.1372E+02
   19   74   22    1.62   12.31  0.3185E+02  0.2085E+02
   34   53  -10    1.62   15.27  0.3950E+02  0.2139E+02
   31   66    0    1.62   15.32  0.3963E+02  0.1839E+02
   30   69  -14    1.62   13.67  0.3538E+02  0.1710E+02
   13   83   24    1.62   13.97  0.3614E+02  0.2181E+02
   15   71   27    1.62   18.58  0.4809E+02  0.2294E+02
   14   99   13    1.62   16.90  0.4373E+02  0.2157E+02
   32   63   -4    1.62   18.58  0.4808E+02  0.2086E+02
   32   63   -9    1.62   11.60  0.3001E+02  0.2059E+02
   16   65   28    1.62   14.41  0.3729E+02  0.1929E+02
   13   78   26    1.62   16.93  0.4381E+02  0.2363E+02
   22   71   19    1.62   14.36  0.3716E+02  0.2100E+02
   31   57  -20    1.62   14.08  0.3644E+02  0.1770E+02
   25   67   16    1.62   13.17  0.3406E+02  0.1733E+02
   13   82   24    1.63   17.28  0.4472E+02  0.2294E+02
   33   54  -14    1.63   12.28  0.3177E+02  0.1727E+02
   15   88   19    1.63   17.94  0.4643E+02  0.2351E+02
   28   70    7    1.63   12.12  0.3137E+02  0.1809E+02
   15   75   25    1.63   12.38  0.3203E+02  0.2282E+02
   31   64    0    1.63   12.10  0.3131E+02  0.1756E+02
   27   72    8    1.64   11.66  0.3017E+02  0.1355E+02
   16   63   28    1.64   13.00  0.3363E+02  0.1806E+02
   30   67    1    1.64   13.62  0.3524E+02  0.1844E+02
   32   60   -3    1.64   11.90  0.3079E+02  0.1886E+02
   13   76   26    1.64   12.08  0.3126E+02  0.2141E+02
   16   85   19    1.64   11.96  0.3095E+02  0.1793E+02
   34   49  -10    1.64   13.81  0.3574E+02  0.2042E+02
   13   85   22    1.64   12.94  0.3347E+02  0.2246E+02
   25   85   -8    1.64   13.22  0.3420E+02  0.1266E+02
   29   68  -16    1.64   11.92  0.3083E+02  0.1737E+02
   15   65   28    1.64   11.67  0.3020E+02  0.1868E+02
   16  100    6    1.64   16.46  0.4260E+02  0.2062E+02
   32   60   -8    1.65   16.94  0.4382E+02  0.2105E+02
   25   72   12    1.65   11.51  0.2979E+02  0.1386E+02
   16   71   25    1.65   11.62  0.3006E+02  0.2178E+02
   15   76   24    1.65   16.56  0.4284E+02  0.2292E+02
   21   64   22    1.65   12.10  0.3132E+02  0.2033E+02
   11   96   17    1.65   12.29  0.3180E+02  0.2251E+02
   14   60   30    1.65   14.33  0.3707E+02  0.1899E+02
   26   64   14    1.65   11.97  0.3098E+02  0.1702E+02
   15   73   25    1.65   11.51  0.2977E+02  0.2198E+02
   30   63  -16    1.65   12.05  0.3118E+02  0.1762E+02
   13   82   23    1.65   11.79  0.3052E+02  0.2223E+02
   17   97    7    1.65   14.42  0.3732E+02  0.2032E+02
   30   67   -2    1.65   13.32  0.3447E+02  0.1252E+02
   18   99   -2    1.65   16.15  0.4179E+02  0.1977E+02
   28   74   -1    1.66   11.70  0.3027E+02  0.1324E+02
   25   71   12    1.66   12.21  0.3158E+02  0.1721E+02
   23   61   20    1.66   11.52  0.2981E+02  0.1982E+02
   28   70  -16    1.66   11.54  0.2986E+02  0.1753E+02
    9   88   23    1.66   14.23  0.3683E+02  0.2232E+02
   25   63   16    1.66   12.79  0.3311E+02  0.1777E+02
   15   75   24    1.66   18.17  0.4701E+02  0.2285E+02
   11   85   23    1.66   16.18  0.4186E+02  0.2282E+02
   29   61  -20    1.66   13.04  0.3374E+02  0.1761E+02
   31   63   -5    1.66   11.82  0.3057E+02  0.1903E+02
   18   98    1    1.66   14.77  0.3821E+02  0.1901E+02
   31   63   -7    1.66   19.96  0.5164E+02  0.2158E+02
   13   79   24    1.66   17.44  0.4513E+02  0.2322E+02
   12   92   18    1.66   13.68  0.3540E+02  0.2205E+02
   29   62  -19    1.66   15.65  0.4049E+02  0.1812E+02
   11   80   25    1.66   12.09  0.3129E+02  0.2206E+02
   31   57  -16    1.66   12.47  0.3227E+02  0.1745E+02
   16   99    5    1.67   14.79  0.3826E+02  0.2073E+02
   14   71   26    1.67   13.11  0.3470E+02  0.2197E+02
   17   99    1    1.67   19.51  0.5165E+02  0.2105E+02
   32   57   -8    1.67   13.23  0.3502E+02  0.1981E+02
   11   86   22    1.67   13.44  0.3557E+02  0.2271E+02
   32   54  -13    1.67   11.50  0.3045E+02  0.1645E+02
   16   99    4    1.67   11.91  0.3152E+02  0.2008E+02
   13  100    9    1.67   17.88  0.4734E+02  0.2163E+02
   20   67   21    1.67   12.19  0.3227E+02  0.2039E+02
   11   48   34    1.67   13.46  0.3562E+02  0.1897E+02
   15   96   10    1.68   12.21  0.3231E+02  0.2038E+02
   29   63  -17    1.68   13.75  0.3640E+02  0.1860E+02
   15   73   24    1.68   13.52  0.3580E+02  0.2197E+02
   18   46   29    1.68   13.47  0.3567E+02  0.1730E+02
   14   67   27    1.68   13.11  0.3470E+02  0.1903E+02
   14   99    8    1.68   11.61  0.3073E+02  0.1994E+02
   32   52  -14    1.68   12.24  0.3239E+02  0.1703E+02
   18   61   25    1.68   12.57  0.3329E+02  0.1868E+02
   12   81   23    1.68   12.47  0.3300E+02  0.2213E+02
   17   32   32    1.68   11.68  0.3091E+02  0.1640E+02
   10   93   18    1.68   13.97  0.3697E+02  0.2281E+02
   28   68  -15    1.68   13.61  0.3603E+02  0.1759E+02
   11   93   17    1.69   12.73  0.3369E+02  0.2243E+02
   14   51   31    1.69   15.51  0.4106E+02  0.1941E+02
    8   92   20    1.69   14.48  0.3833E+02  0.2274E+02
    6   80   27    1.69   12.05  0.3190E+02  0.2274E+02
   28   70  -11    1.69   12.18  0.3223E+02  0.1348E+02
   10   89   20    1.69   12.66  0.3352E+02  0.2190E+02
    8   88   22    1.69   18.82  0.4981E+02  0.2401E+02
   14   76   23    1.69   12.04  0.3186E+02  0.2152E+02
   22   54   22    1.70   13.87  0.3671E+02  0.1805E+02
   13  101    5    1.70   11.87  0.3141E+02  0.1965E+02
   15   98    5    1.70   14.74  0.3903E+02  0.1990E+02
   20   47   26    1.70   12.79  0.3384E+02  0.1714E+02
   23   61   18    1.70   12.53  0.3317E+02  0.1989E+02
   20   69   19    1.70   12.29  0.3253E+02  0.2011E+02
   10   98   13    1.70   12.38  0.3276E+02  0.2191E+02
   21   54   23    1.70   12.36  0.3271E+02  0.1838E+02
   23   58   19    1.70   13.01  0.3443E+02  0.1965E+02
   23   55   20    1.70   12.83  0.3396E+02  0.1919E+02
   14   75   23    1.70   12.69  0.3360E+02  0.2200E+02
    8   87   22    1.71   15.76  0.4171E+02  0.2322E+02
    8   83   24    1.71   15.07  0.3988E+02  0.2239E+02
   14   96    9    1.71   12.36  0.3272E+02  0.1985E+02
   26   69    7    1.71   13.32  0.3525E+02  0.1343E+02
   20   57   23    1.71   13.94  0.3690E+02  0.1922E+02
   32   51   -4    1.71   12.55  0.3321E+02  0.1815E+02
    6   87   23    1.71   12.30  0.3257E+02  0.2173E+02
   12   58   30    1.71   13.24  0.3506E+02  0.1867E+02
   16   98   -3    1.71   13.53  0.3582E+02  0.1959E+02
   28   65  -15    1.71   13.82  0.3658E+02  0.1748E+02
   13   98    8    1.71   13.30  0.3520E+02  0.2020E+02
   29   62  -14    1.71   14.03  0.3714E+02  0.1773E+02
   25   78  -10    1.71   13.53  0.3582E+02  0.1368E+02
   11   86   20    1.71   18.44  0.4880E+02  0.2268E+02
   32   51   -8    1.71   13.97  0.3698E+02  0.1966E+02
   29   54  -20    1.71   12.41  0.3285E+02  0.1677E+02
   11   47   33    1.72   13.28  0.3516E+02  0.1838E+02
   13   94   12    1.72   12.55  0.3323E+02  0.1352E+02
   32   50   -4    1.72   13.47  0.3567E+02  0.1828E+02
   23   56   19    1.72   11.82  0.3130E+02  0.1953E+02
   11   93   15    1.72   13.65  0.3612E+02  0.2049E+02
   24   81   -6    1.72   11.75  0.3110E+02  0.1061E+02
   27   64  -18    1.72   11.81  0.3127E+02  0.1825E+02
   29   56    6    1.72   14.72  0.3374E+02  0.1634E+02
   23   63   16    1.73   14.00  0.3210E+02  0.1320E+02
    9   91   18    1.73   15.45  0.3542E+02  0.2171E+02
    8   89   20    1.73   16.23  0.3722E+02  0.2247E+02
   28   67  -10    1.73   14.06  0.3223E+02  0.1645E+02
   21   63   19    1.73   21.00  0.4816E+02  0.2091E+02
   29   54    7    1.73   13.51  0.3099E+02  0.1678E+02
   12   47   32    1.73   12.69  0.2910E+02  0.1781E+02
   30   58  -11    1.73   12.36  0.2834E+02  0.1358E+02
   31   54   -9    1.73   15.00  0.3439E+02  0.1293E+02
   31   54   -8    1.73   16.53  0.3790E+02  0.2010E+02
   32   43  -15    1.73   11.64  0.2668E+02  0.1510E+02
   27   54  -23    1.73   14.27  0.3272E+02  0.1755E+02
   10   92   16    1.73   14.17  0.3250E+02  0.2245E+02
   21   62   19    1.74   13.97  0.3204E+02  0.1944E+02
   27   65  -16    1.74   14.74  0.3379E+02  0.1786E+02
   22   52   21    1.74   14.91  0.3419E+02  0.1822E+02
   20   60   21    1.74   14.20  0.3255E+02  0.1826E+02
   18   61   23    1.74   12.07  0.2768E+02  0.1950E+02
    8   94   16    1.74   22.13  0.5074E+02  0.2342E+02
   22   48   22    1.74   14.81  0.3397E+02  0.1750E+02
   19   90    1    1.74   12.79  0.2933E+02  0.1104E+02
   11   95   12    1.74   12.26  0.2811E+02  0.1979E+02
    7   96   15    1.74   12.32  0.2825E+02  0.2232E+02
   17   40   29    1.74   13.34  0.3058E+02  0.1718E+02
   13   82   19    1.74   12.53  0.2874E+02  0.1776E+02
   30   51    4    1.74   11.66  0.2673E+02  0.1611E+02
   17   66   22    1.75   12.96  0.2972E+02  0.1911E+02
   17   82   14    1.75   12.06  0.2765E+02  0.1116E+02
   22   54   20    1.75   19.58  0.4489E+02  0.1932E+02
   24   71  -18    1.75   15.00  0.3440E+02  0.1773E+02
    9   53   32    1.75   13.49  0.3092E+02  0.1803E+02
    7   93   17    1.75   11.88  0.2725E+02  0.2207E+02
   12   49   31    1.75   12.00  0.2752E+02  0.1726E+02
    7   90   19    1.75   12.12  0.2779E+02  0.2170E+02
    9   86   20    1.75   14.06  0.3223E+02  0.2202E+02
   32   46   -6    1.75   13.21  0.3028E+02  0.1794E+02
   16   56   26    1.75   11.53  0.2643E+02  0.1811E+02
   11  100    4    1.75   14.51  0.3327E+02  0.1997E+02
   21   47   23    1.75   11.88  0.2724E+02  0.1664E+02
    8   85   21    1.76   14.89  0.3415E+02  0.2184E+02
   18   62   22    1.76   13.60  0.3118E+02  0.1981E+02
   10   93   14    1.76   13.96  0.3202E+02  0.2118E+02
    5   88   21    1.76   11.62  0.2665E+02  0.2163E+02
    3   84   24    1.76   13.98  0.3205E+02  0.2215E+02
   22   61   17    1.76   12.88  0.2954E+02  0.1947E+02
    9   56   31    1.76   12.12  0.2779E+02  0.1761E+02
   11   55   30    1.76   11.96  0.2743E+02  0.1808E+02
   24   62   13    1.76   14.77  0.3386E+02  0.1332E+02
   16   47   28    1.76   14.79  0.3391E+02  0.1719E+02
    8   96   13    1.76   17.32  0.3971E+02  0.2088E+02
   26   62  -19    1.76   14.82  0.3397E+02  0.1658E+02
    8   91   17    1.76   16.15  0.3704E+02  0.2181E+02
   24   72  -16    1.76   11.90  0.2729E+02  0.1675E+02
   22   52   20    1.76   13.36  0.3062E+02  0.1823E+02
   21   59   19    1.76   13.51  0.3098E+02  0.1877E+02
   13   98    1    1.76   12.01  0.2753E+02  0.1966E+02
   15   67   23    1.76   12.55  0.2877E+02  0.1379E+02
   28   52  -20    1.77   12.73  0.2918E+02  0.1628E+02
    6   94   16    1.77   12.39  0.2840E+02  0.2119E+02
   31   43  -16    1.77   14.08  0.3227E+02  0.1557E+02
   20   30   27    1.77   11.67  0.2675E+02  0.1609E+02
    7   77   25    1.77   16.83  0.3860E+02  0.2334E+02
   16   61   24    1.77   11.73  0.2690E+02  0.1386E+02
   10   97    9    1.77   15.96  0.3659E+02  0.2044E+02
   21   61   18    1.77   13.37  0.3066E+02  0.1901E+02
   25   74   -6    1.77   12.33  0.2827E+02  0.8952E+01
    7   79   24    1.77   12.94  0.2967E+02  0.2200E+02
   12   46   31    1.77   14.52  0.3330E+02  0.1748E+02
   29   53  -16    1.77   14.39  0.3300E+02  0.1680E+02
   10   96   10    1.77   18.96  0.4348E+02  0.2039E+02
   12   99   -1    1.77   13.50  0.3094E+02  0.1962E+02
    4   82   24    1.77   14.09  0.3230E+02  0.2182E+02
   24   71  -16    1.77   12.05  0.2762E+02  0.1632E+02
   21   84   -7    1.77   12.75  0.2923E+02  0.1345E+02
   30   52  -12    1.77   14.46  0.3316E+02  0.1616E+02
   25   73   -9    1.77   12.41  0.2845E+02  0.1336E+02
   17   87  -15    1.77   16.14  0.3701E+02  0.1784E+02
   25   69  -15    1.77   14.44  0.3312E+02  0.1708E+02
   10  100    3    1.78   13.70  0.3141E+02  0.2002E+02
   15   95   -2    1.78   12.75  0.2924E+02  0.1784E+02
   26   55  -22    1.78   14.95  0.3428E+02  0.1639E+02
   32   42   -6    1.78   17.35  0.3978E+02  0.1878E+02
   33   34   -7    1.78   15.93  0.3652E+02  0.1813E+02
   31   48   -5    1.78   12.19  0.2795E+02  0.1183E+02
    8   99    8    1.78   13.24  0.3037E+02  0.1987E+02
   32   41   -4    1.78   12.98  0.2976E+02  0.1727E+02
    5   96   14    1.78   13.15  0.3015E+02  0.2237E+02
   33   33   -5    1.78   15.01  0.3440E+02  0.1740E+02
   12   98    1    1.78   20.38  0.4672E+02  0.2059E+02
   23   48   19    1.78   15.04  0.3448E+02  0.1882E+02
   33   33   -7    1.78   20.16  0.4622E+02  0.1936E+02
   32   38   -1    1.79   13.55  0.3107E+02  0.1652E+02
    0   77   27    1.79   13.00  0.2981E+02  0.2083E+02
   31   43    1    1.79   11.86  0.2718E+02  0.1511E+02
   28   55    5    1.79   14.56  0.3337E+02  0.1618E+02
   17   90   -9    1.79   13.31  0.3051E+02  0.1740E+02
   24   59   13    1.79   11.68  0.2678E+02  0.1285E+02
    8   97   10    1.79   11.65  0.2671E+02  0.2027E+02
   15   94   -2    1.79   13.51  0.4143E+02  0.1911E+02
   13   96   -1    1.80   11.84  0.3631E+02  0.1957E+02
    1   89   20    1.80   12.36  0.3789E+02  0.2284E+02
   15   66   22    1.80   13.42  0.4116E+02  0.1453E+02
   21   51   20    1.81   14.31  0.4389E+02  0.1826E+02
    5   94   14    1.81   13.54  0.4151E+02  0.2282E+02
   30   45  -14    1.81   12.81  0.3929E+02  0.1620E+02
   33   27   -7    1.81   13.83  0.4242E+02  0.1781E+02
   24   45   17    1.82   15.56  0.4773E+02  0.1982E+02
   14   29   31    1.82   12.96  0.3974E+02  0.1586E+02
    8   97    7    1.82   13.17  0.4038E+02  0.2060E+02
   10   97   -7    1.83   13.91  0.4267E+02  0.2134E+02
   33   23   -7    1.83   14.55  0.4462E+02  0.1754E+02
   11   96   -6    1.83   13.89  0.4259E+02  0.2041E+02
   17   37   27    1.83   11.58  0.3551E+02  0.1642E+02
    8   99    0    1.83   11.60  0.3559E+02  0.1990E+02
    1   87   20    1.83   15.59  0.4782E+02  0.2345E+02
   23   51   16    1.84   16.29  0.4995E+02  0.2032E+02
    5  100    3    1.84   13.10  0.4016E+02  0.2150E+02
    0   87   20    1.84   11.94  0.3661E+02  0.2222E+02
   32   28   -2    1.85   12.39  0.3800E+02  0.1613E+02
   32   30   -7    1.85   11.61  0.3561E+02  0.1382E+02
   22   56   15    1.85   11.60  0.3557E+02  0.1056E+02
   26   58  -15    1.86   12.65  0.3880E+02  0.1782E+02
   25   56  -19    1.86   12.00  0.3679E+02  0.1689E+02
    7   98   -3    1.86   11.85  0.3940E+02  0.2170E+02
   24   48   14    1.87   11.60  0.3858E+02  0.1940E+02
   32   20   -1    1.88   13.79  0.4586E+02  0.1504E+02
   24   36   17    1.88   13.09  0.4352E+02  0.1893E+02
   31   27    1    1.89   12.03  0.4001E+02  0.1481E+02
   29   46   -4    1.89   18.64  0.6202E+02  0.2039E+02
   31   26    1    1.89   13.63  0.4533E+02  0.1516E+02
   17   27   27    1.89   15.85  0.5273E+02  0.1809E+02
    8   68   24    1.89   11.93  0.3970E+02  0.1411E+02
   22   29   21    1.90   11.78  0.3919E+02  0.1616E+02
   30   27    4    1.91   11.65  0.3874E+02  0.1350E+02
   31   27  -11    1.91   15.46  0.5142E+02  0.1639E+02
   24   32   17    1.91   14.95  0.4973E+02  0.1907E+02
   30   32    1    1.91   12.84  0.4272E+02  0.1539E+02
    9   92    5    1.92   28.74  0.9559E+02  0.1972E+02
   25   46  -21    1.92   11.88  0.3951E+02  0.1727E+02
   11   57   25    1.92   11.57  0.3849E+02  0.1888E+02
   31   21    0    1.92   13.55  0.4508E+02  0.1550E+02
    7   55   28    1.92   11.90  0.3960E+02  0.1845E+02
   32    7   -3    1.92   14.77  0.4914E+02  0.1563E+02
   32    4   -3    1.93   15.83  0.5265E+02  0.1511E+02
   32    3   -3    1.93   11.67  0.3881E+02  0.1361E+02
   31   21   -1    1.93   13.24  0.4403E+02  0.1021E+02
    5   78  -22    1.93   11.81  0.3930E+02  0.1978E+02
   30   34   -6    1.94   12.01  0.3994E+02  0.1532E+02
   31   20   -1    1.94   12.76  0.4244E+02  0.1011E+02
   30   30  -12    1.94   15.56  0.5174E+02  0.1633E+02
   31   22   -3    1.94   20.13  0.6695E+02  0.1625E+02
   31   19   -1    1.94   15.98  0.5316E+02  0.1059E+02
   31   23   -5    1.94   12.84  0.4271E+02  0.1115E+02
   21   10   23    1.94   13.90  0.4623E+02  0.1525E+02
   31   22   -4    1.94   11.84  0.3938E+02  0.1780E+02
   31   21   -3    1.94   12.67  0.4214E+02  0.1321E+02
   31   19   -3    1.95   12.76  0.6356E+02  0.1273E+02
   31   16   -3    1.96   11.57  0.5763E+02  0.1171E+02
   31   13   -3    1.97   13.76  0.6856E+02  0.1247E+02
   31   11   -3    1.98   13.11  0.6532E+02  0.1215E+02
   31    9   -3    1.98   14.08  0.7012E+02  0.1882E+02
   31    8   -3    1.99   12.24  0.6100E+02  0.1736E+02
   31    7   -3    1.99   13.60  0.6774E+02  0.1756E+02
   29   34   -6    1.99   11.98  0.5969E+02  0.1376E+02
   29   32  -10    2.00   13.92  0.6936E+02  0.1925E+02
   29   32   -2    2.00   13.62  0.6785E+02  0.1361E+02
   29   30  -10    2.01   12.04  0.5996E+02  0.1779E+02
   30   21   -4    2.01   12.13  0.6044E+02  0.2005E+02
   28   32  -14    2.01   15.74  0.7841E+02  0.2071E+02
   30   17   -4    2.03   12.17  0.6064E+02  0.1263E+02
   29   21   -9    2.07   12.15  0.8478E+02  0.1544E+02
   29   20   -3    2.08   17.22  0.1202E+03  0.1817E+02
   29   18   -3    2.09   13.74  0.9589E+02  0.1550E+02
   27   25  -16    2.09   13.14  0.9164E+02  0.2258E+02
   29   16   -3    2.10   16.46  0.1148E+03  0.1737E+02
   29    8    0    2.10   14.70  0.1026E+03  0.1161E+02
   29    7    0    2.10   13.19  0.9200E+02  0.1069E+02
   29    6    0    2.10   14.46  0.1009E+03  0.1149E+02
   29    5    0    2.10   13.74  0.9587E+02  0.1376E+02
   29    4    0    2.10   16.44  0.1147E+03  0.1580E+02
   29    3    0    2.10   14.51  0.1013E+03  0.1442E+02
   29    3   -1    2.11   13.71  0.9567E+02  0.1424E+02
   28   25   -9    2.11   11.77  0.8214E+02  0.2023E+02
   29    2   -1    2.12   26.22  0.1830E+03  0.2402E+02
   29    1   -1    2.12   20.64  0.1440E+03  0.1946E+02
   29    0   -1    2.12   12.65  0.2116E+03  0.3843E+02
   27   28  -13    2.12   12.06  0.8411E+02  0.2105E+02
   28   25   -5    2.12   14.67  0.1023E+03  0.1681E+02
   28   14    2    2.13   12.98  0.9056E+02  0.1905E+02
   28    8    3    2.13   13.78  0.9615E+02  0.1801E+02
   28   23   -5    2.14   24.57  0.1714E+03  0.2426E+02
   28   23   -6    2.14   13.44  0.9379E+02  0.1605E+02
   28   10    2    2.14   15.03  0.1049E+03  0.2068E+02
   28   22   -7    2.14   23.26  0.1623E+03  0.2358E+02
   28   22   -5    2.14   16.82  0.1174E+03  0.1822E+02
   28   10    1    2.15   12.40  0.1163E+03  0.2316E+02
   28    9    1    2.16   14.76  0.1385E+03  0.2640E+02
   28    8    1    2.16   14.61  0.1371E+03  0.2600E+02
   28   19   -5    2.16   11.56  0.1085E+03  0.1713E+02
   28   10  -11    2.17   13.30  0.1247E+03  0.2661E+02
   28   11   -1    2.17   16.33  0.1532E+03  0.2067E+02
   27    4    2    2.23   11.82  0.1109E+03  0.2057E+02
   27   14   -4    2.26   13.09  0.1229E+03  0.1805E+02
   27   12   -3    2.27   16.39  0.1537E+03  0.2128E+02
   27   12   -4    2.27   14.54  0.1364E+03  0.2647E+02
   27   10   -3    2.27   12.25  0.1149E+03  0.1651E+02
   27    1   -1    2.27   14.74  0.1383E+03  0.2646E+02
   27    7   -9    2.28   11.66  0.1094E+03  0.2342E+02
   27    6   -8    2.28   13.48  0.1784E+03  0.3438E+02
   27    5   -3    2.29   12.94  0.1712E+03  0.3167E+02
   26   17   -8    2.32   12.81  0.1695E+03  0.3379E+02
   26   17   -4    2.33   13.31  0.1761E+03  0.2405E+02
   26    8   -9    2.36   15.53  0.2055E+03  0.3927E+02
   26    1   -1    2.36   20.33  0.2690E+03  0.4788E+02
   25   14    2    2.37   16.87  0.2232E+03  0.4010E+02
   25   12    2    2.38   14.22  0.1882E+03  0.3421E+02
   25   22   -5    2.38   11.55  0.1528E+03  0.2128E+02
   25   20   -8    2.39   11.56  0.1530E+03  0.2184E+02
   25   20   -4    2.40   15.44  0.2043E+03  0.3834E+02
   25   19   -6    2.41   16.21  0.2145E+03  0.2885E+02
   24   19  -12    2.43   12.70  0.1680E+03  0.3413E+02
   24   14  -10    2.50   13.96  0.2607E+03  0.4892E+02
   24   11   -5    2.56   11.55  0.2157E+03  0.3965E+02
   24    3   -1    2.56   17.02  0.3178E+03  0.3984E+02
   24   10   -5    2.56   17.94  0.3350E+03  0.5945E+02
   24    8   -7    2.56   12.50  0.2334E+03  0.4327E+02
   24    8   -5    2.57   14.18  0.2648E+03  0.4788E+02
   23   16  -10    2.58   15.80  0.2950E+03  0.5420E+02
   23   17    0    2.58   12.85  0.2399E+03  0.3037E+02
   24    3   -6    2.58   13.13  0.2451E+03  0.4466E+02
   24    3   -4    2.58   25.41  0.4743E+03  0.8315E+02
   24    2   -6    2.58   24.86  0.4642E+03  0.8165E+02
   24    2   -4    2.58   22.96  0.4286E+03  0.7555E+02
   24    3   -5    2.58   26.00  0.4854E+03  0.8521E+02
   24    1   -6    2.58   11.59  0.2164E+03  0.3964E+02
   24    1   -4    2.58   13.77  0.2571E+03  0.4647E+02
   23   15  -10    2.59   13.31  0.2485E+03  0.4670E+02
   23   14  -10    2.60   17.26  0.3222E+03  0.5882E+02
   23   12   -9    2.62   12.35  0.2305E+03  0.4329E+02
   23   11   -9    2.63   12.19  0.2277E+03  0.4280E+02
   21   12  -14    2.70   12.80  0.3599E+03  0.6518E+02
   21   11    5    2.72   12.43  0.3497E+03  0.6080E+02
   16   13  -17    3.04   11.50  0.6086E+03  0.6531E+02
   20    2    0    3.05   12.40  0.6560E+03  0.8059E+02
   19    6  -10    3.13   11.51  0.6089E+03  0.1071E+03
   19    4  -10    3.14   12.53  0.6627E+03  0.8300E+02
   19   11   -6    3.19   11.67  0.6172E+03  0.5496E+02
   19   11   -2    3.19   15.31  0.8101E+03  0.8192E+02
   19    6   -7    3.22   12.17  0.6437E+03  0.7962E+02
   18    1   -4    3.45   15.51  0.1760E+04  0.1525E+03
   12   16   -1    4.71   12.22  0.2067E+04  0.2076E+03
    7    1   -2    8.83   13.97  0.3484E+04  0.3001E+03
    7    1   -1    8.84   14.66  0.3654E+04  0.3658E+03
    0    4    0   46.73   12.50  0.6234E+04  0.1071E+04

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        21
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"      324292
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)  320845
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))     3447


 cpu time used                 55.8 sec
 elapsed wall-clock time       63.8 sec


 XDS prints this sort of information during integration. Please also show the current statistics in the logs. ***** INTEGRATE ***** (VERSION January 10, 2014)   2-Oct-2014

 NAME_TEMPLATE_OF_DATA_FRAMES=/data/results/LLS_Williams/MBP4-Rac1/76c1b_0002212/76c1b_0002212_1_0???.cbf CBF
 DATA_RANGE=       1     900
 MINIMUM_VALID_PIXEL_VALUE=     0    OVERLOAD=   1048500
 BACKGROUND_PIXEL=    6.00      SIGNAL_PIXEL=    3.00
 MAXIMUM_ERROR_OF_SPOT_POSITION=   3.0
 MAXIMUM_ERROR_OF_SPINDLE_POSITION=   2.0
 MINPK=  75.00000
 PROFILE_FITTING= TRUE
 SPACE_GROUP_NUMBER=    3
 UNIT_CELL_CONSTANTS=    62.04   186.91    72.36  90.000 100.183  90.000
 ROTATION_AXIS=  0.99994 -0.00029 -0.01089
 OSCILLATION_RANGE=  0.200000 DEGREES
 STARTING_ANGLE=   0.000 STARTING_FRAME=       1
 X-RAY_WAVELENGTH=   0.97949 ANGSTROM
 INCIDENT_BEAM_DIRECTION=  0.009950  0.003805  1.020884
 NUMBER OF DETECTOR SEGMENTS   1
 NX=  2463 NY=  2527 QX= 0.17200 QY= 0.17200
 ORGX=   1222.30  ORGY=   1201.42
 DETECTOR_DISTANCE=   384.979
 DIRECTION_OF_DETECTOR_X-AXIS=   1.00000   0.00000   0.00000
 DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000   1.00000   0.00000
 SEGMENT=    1 2463    1 2527  ! x1,x2,y1,y2 (pixel units)
    SEGMENT_ORGX=      0.00   SEGMENT_ORGY=      0.00
    SEGMENT_DISTANCE=     0.00
    DIRECTION_OF_SEGMENT_X-AXIS=  1.000000  0.000000  0.000000
    DIRECTION_OF_SEGMENT_Y-AXIS=  0.000000  1.000000  0.000000
 RESOLUTION RANGE RECORDED BY DETECTOR (ANGSTROM)    69.090     1.468
 NUMBER OF TRUSTED DETECTOR PIXELS    5516754
 MEAN CONTENTS OF TRUSTED PIXELS IN BACKGROUND TABLE    168.100
 MEAN VALUE OF NON-XRAY BACKGROUND (OFFSET)=    0.00
 NUMBER OF X-RAY COUNTS EQUIVALENT TO PIXEL CONTENTS IN A DATA IMAGE
 (PIXEL VALUE - OFFSET)/GAIN ,    WITH GAIN=    0.99



 ******************************************************************************
                     PROCESSING OF IMAGES        1 ...      26
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
     1   0  1.069  5481561    0  14648     464    26  0.0341  0.4626
     2   0  1.078  5482498    0  14617     457    33  0.0356  0.4730
     3   0  1.070  5482149    0  14577     452    31  0.0354  0.4662
     4   0  1.086  5482864    0  14613     454    32  0.0342  0.4975
     5   0  1.073  5482487    0  14588     459    33  0.0343  0.4839
     6   0  1.077  5483432    0  14617     433    22  0.0348  0.4732
     7   0  1.068  5484280    0  14603     428    25  0.0344  0.4531
     8   0  1.062  5482017    0  14651     444    24  0.0347  0.4581
     9   0  1.082  5481630    0  14641     434    22  0.0368  0.4391
    10   0  1.067  5481934    0  14610     442    18  0.0359  0.4258
    11   0  1.066  5481621    0  14576     438    17  0.0359  0.4650
    12   0  1.061  5481851    0  14617     433    20  0.0383  0.4640
    13   0  1.060  5481853    0  14586     442    23  0.0350  0.4694
    14   0  1.069  5482349    0  14564     438    23  0.0350  0.4898
    15   0  1.069  5482911    0  14637     438    25  0.0347  0.5090
    16   0  1.059  5482627    0  14621     457    26  0.0346  0.4970
    17   0  1.079  5482993    0  14593     437    31  0.0375  0.4671
    18   0  1.068  5482225    0  14597     467    28  0.0350  0.4662
    19   0  1.077  5482780    0  14598     441    31  0.0361  0.4638
    20   0  1.070  5481588    0  14560     460    26  0.0356  0.4637
    21   0  1.070  5481815    0  14579     463    25  0.0356  0.4544
    22   0  1.070  5482076    0  14591     449    19  0.0361  0.5256
    23   0  1.066  5483665    0  14574     432    22  0.0355  0.5115
    24   0  1.072  5484728    0  14536     432    25  0.0353  0.5163
    25   0  1.059  5484031    0  14622     407    22  0.0354  0.5011
    26   0  1.066  5484371    0  14597     418    25  0.0357  0.5186

  1747 OUT OF   1951 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.25
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.42
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.010   187.029    72.279  90.000 100.190  90.000
 REC. CELL PARAMETERS   0.016385  0.005347  0.014057  90.000  79.810  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.066    18.518    18.947
 COORDINATES OF UNIT CELL B-AXIS   -14.408  -110.252   150.389
 COORDINATES OF UNIT CELL C-AXIS    18.334   -57.212   -40.186
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.012     0.055    -0.055
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1747 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  515    141    124     95     96    159    233    200    175
 SIGMAB (degree)       0.066  0.068  0.070  0.071  0.070  0.067  0.066  0.067  0.068
 SIGMAR (degree)       0.433  0.379  0.361  0.390  0.415  0.407  0.402  0.436  0.426


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1747 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  0  0  0     0  0  0  1  1  1  0  0  0     0  0  0  2  2  2  1  0  0
     0  0  0  0  0  2  1  1  0     0  0  1  4  4  4  2  1  0     0  0  3  9  9  9  4  1  0
     0  0  0  0  0  2  1  0  0     0  0  2  3  4  7  3  1  0     0  1  5  8 12 17  7  1  0
     0  0  0  0  0  2  0  0  0     0  1  1  3  6  6  1  0  0     0  1  3  7 14 14  4  0  0
     0  0  0  0  0  0  0  0  0     0  0  1  2  2  0  0  0  0     0  0  2  4  5  2  0  0  0
     0  0  0  0  0  0  0  1  0     0  0  0  1  0  0  0  0  0     0  0  0  1  0 -1 -1 -1  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -1  0 -1 -1  0  0  0


     0  0 -1 -1 -1  0 -1  0  0     0  0 -1 -1  0  0 -1  0  0     0  0 -1 -1 -1  0 -1  0  0
     0 -1 -1  0  1  1  1  0  0     0 -1 -1  0  1  2  1  0  0     0 -1 -1  0  1  1  0  0  0
     0 -1  1  4  5  6  3  1  0    -1  0  1  2  6  7  3  1  0     0  0  1  3  6  4  3  0  0
     0  0  5 13 21 17  6  1 -1     0  1  4 13 33 20  6  1  0     0  1  5 15 25 15  6  1  0
     1  2  8 20 46 30  9  2  0     1  2  7 35100 35  8  2  0     1  2  8 27 45 26  9  2  0
     0  1  6 15 28 17  6  1  0     0  2  6 23 42 15  5  1  0     0  2  7 18 26 17  7  1 -1
     0  0  3  6  7  3  1 -1  0     0  1  3  6  6  2  1  0  0    -1  1  3  5  7  5  1 -1  0
     0  0  1  1  1 -1 -1 -1  0     0  0  1  1  1 -1 -1 -2  0     0  0  0  1  1 -1 -1 -1  0
     0  0 -1  0 -1 -1 -2  0  0     0  0 -1  0  0 -1 -2  0  0     0  0 -1 -1 -1 -1 -2  0  0


     0  0  0 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1 -1 -1  0  0  0  0  0     0  0  0  0  0  1  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  2  3  2  1  0  0     0  0  0  1  1  2  1  0  0     0  0  0  0  0  0  0  0  0
     0  1  4 13 13  7  3  0  0     0  1  2  5  5  3  2  1  0     0  0  0  3  0  0  0  0  0
     1  1  8 15 11 10  5  1  0     0  1  4  6  4  4  2  0  0     0  0  1  3  0  0  0  0  0
     0  1  5  9 10 12  3  0  0     0  1  2  3  4  5  1  0  0     0  0  1  2  0  0  0  0  0
    -1  0  1  3  4  3  0  0  0     0  0  1  2  2  1  0  0  0     0  0  1  1  1  0  0  0  0
     0 -1  0  0  1 -1 -1  0  0     0  0  0  0  1  0 -1  0  0     0  0  0  0  1  0  0  0  0
     0  0 -1  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30469 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES       27 ...      52
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
    27   0  1.067  5484952    0  14550     417    13  0.0346  0.5482
    28   0  1.062  5483956    0  14566     426    20  0.0361  0.4962
    29   0  1.054  5483611    0  14605     440    27  0.0350  0.4521
    30   0  1.059  5483862    0  14591     436    36  0.0361  0.4174
    31   0  1.055  5483288    0  14549     416    23  0.0375  0.4244
    32   0  1.060  5483507    0  14619     423    22  0.0367  0.4403
    33   0  1.059  5482522    0  14636     427    24  0.0366  0.4244
    34   0  1.062  5482368    0  14612     419    22  0.0369  0.4627
    35   0  1.058  5482547    0  14607     424    25  0.0375  0.4428
    36   0  1.061  5483635    0  14598     422    17  0.0357  0.4758
    37   0  1.063  5484164    0  14619     399    14  0.0361  0.5043
    38   0  1.069  5483925    0  14617     421    19  0.0366  0.4811
    39   0  1.073  5483495    0  14644     422    29  0.0360  0.4168
    40   0  1.064  5483820    0  14595     400    30  0.0374  0.4219
    41   0  1.054  5483001    0  14608     409    28  0.0367  0.4576
    42   0  1.045  5483755    0  14616     406    28  0.0369  0.4924
    43   0  1.058  5483862    0  14642     404    34  0.0376  0.4391
    44   0  1.063  5483464    0  14588     412    31  0.0367  0.4582
    45   0  1.069  5483764    0  14640     411    41  0.0363  0.4337
    46   0  1.058  5484949    0  14638     406    39  0.0379  0.4291
    47   0  1.055  5485087    0  14628     396    34  0.0385  0.4415
    48   0  1.064  5484951    0  14678     386    27  0.0364  0.4932
    49   0  1.056  5485177    0  14680     373    29  0.0395  0.4481
    50   0  1.047  5484724    0  14608     377    34  0.0386  0.4206
    51   0  1.063  5485829    0  14622     379    23  0.0372  0.4650
    52   0  1.052  5486491    0  14640     361    29  0.0387  0.3980

  1610 OUT OF   1794 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.30
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.40
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.022   186.907    72.256  90.000 100.200  90.000
 REC. CELL PARAMETERS   0.016382  0.005350  0.014062  90.000  79.800  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.069    18.520    18.973
 COORDINATES OF UNIT CELL B-AXIS   -14.464  -110.175   150.288
 COORDINATES OF UNIT CELL C-AXIS    18.331   -57.198   -40.167
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.010     0.061    -0.075
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1610 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  480    135    111     82     89    149    211    177    161
 SIGMAB (degree)       0.068  0.070  0.071  0.072  0.071  0.069  0.067  0.069  0.070
 SIGMAR (degree)       0.420  0.351  0.347  0.380  0.403  0.407  0.402  0.418  0.400


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1610 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -1 -1  0  0  0  0  0
     0  0  0  0  0  1  0 -1  0     0  0  0  0  1  0  0 -1  0     0  0 -1  0  0  0  0 -1  0
     0  0  0  0  1  1  1  0  0     0  0  0  1  1  1  1  0  0     0 -1  0  3  3  2  1  0  0
     0 -1  0  3  2  1  1  1  0     0  0  1  4  4  3  2  1  0     0  0  3  8  8  8  4  1 -1
     0  0  1  2  1  2  1  0  0     0  0  2  2  3  5  3  1  0     0  1  5  6 10 15  7  1  0
     0  0  1  3  3  2  0  0  0     0  0  2  3  6  5  1  0  0     0  1  4  6 13 14  5  0  0
     0  0  0  2  1  0  0  1  0     0  0  1  2  3  0  0  0  0    -1  0  2  4  5  3  0 -1  0
     0  0  0  0  0  0  0  1  0     0  0  0  1  1  0 -1  0  0     0 -1  0  1  0 -1 -1 -1  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -1  0  0 -1  0  0  0


     0  0 -1  0  0  0  0  0  0     0  0 -1 -1  0  0 -1  0  0     0  0 -2 -1 -1  0  0  0  0
     0 -1 -1  0  1  1  0 -1  0     0 -1  0  0  1  2  1 -1  0     0 -1 -1 -1  1  1  0 -1  0
    -1 -1  2  6  6  6  3  0  0    -1  0  2  3  7  7  3  1  0    -1  0  1  4  7  5  3  0  0
    -1  0  5 12 20 17  6  1  0    -1  1  5 15 35 20  6  2  0    -1  1  6 16 24 15  6  1  0
     0  2  7 19 43 30  9  2  0     0  2  7 36100 36  8  2  0     0  2  9 26 42 25  8  1  0
     0  1  6 15 29 18  7  1 -1     0  2  6 22 43 17  6  1 -1     0  2  7 16 24 17  7  1  0
    -1  0  3  7  8  5  1  0  0     0  1  3  6  7  4  2  0  0    -1  0  3  5  7  5  2 -1  0
     0 -1  0  1  0 -1 -1 -1  0     0  0  1  2  1  0 -1 -2  0     0 -1  0  1  1 -1 -1 -1  0
     0  0 -1  0 -1 -1 -2  0  0     0  0 -1 -1  0 -1 -2  0  0     0  0 -1  0 -1 -1 -2  0  0


     0  0 -1 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1 -1 -1  0  0  0  0  0     0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  3  4  3  1  0  0     0  1  0  1  2  1  1  0  0     0  0 -1  0  1  0  0  0  0
     1  1  5 12 11  6  3  1  0     0  1  2  5  4  3  2  0  0     0  1  0  2  2  0  0  0  0
     0  1  8 14  9  8  5  1  0     0  1  3  4  2  3  2  0  0     0  0  1  2  0  0  0  0  0
     0  1  5  7  9 11  4  0  0     0  0  2  3  4  4  2  0  0     0  0  1  1  2  0  0  0  0
    -1 -1  1  2  4  3  0 -1  0     0  0  1  2  2  1  0  0  0     0  0  1  1  1  0  0  0  0
     0 -1 -1  0  1  0 -1 -1  0     0 -1  0  0  1  0 -1  0  0     0  0  0  0  1  0  0  0  0
     0  0 -1  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30456 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES       53 ...      78
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
    53   0  1.051  5486366    0  14632     358    26  0.0397  0.3727
    54   0  1.059  5485906    0  14639     364    31  0.0376  0.3804
    55   0  1.050  5485309    0  14601     379    25  0.0370  0.4386
    56   0  1.070  5484775    0  14601     373    28  0.0384  0.4396
    57   0  1.062  5484859    0  14622     362    22  0.0397  0.4598
    58   0  1.067  5485249    0  14611     362    23  0.0380  0.4571
    59   0  1.053  5484162    0  14581     356    24  0.0396  0.4171
    60   0  1.064  5484587    0  14572     359    24  0.0377  0.4263
    61   0  1.066  5485187    0  14570     362    31  0.0385  0.3850
    62   0  1.056  5486516    0  14565     332    31  0.0390  0.4006
    63   0  1.052  5486173    0  14547     336    25  0.0383  0.4377
    64   0  1.060  5486590    0  14602     331    33  0.0394  0.3772
    65   0  1.054  5486233    0  14586     337    23  0.0417  0.4140
    66   0  1.061  5486130    0  14556     338    18  0.0397  0.4398
    67   0  1.041  5485890    0  14577     351    24  0.0390  0.4379
    68   0  1.057  5485313    0  14571     361    21  0.0393  0.4781
    69   0  1.066  5486182    0  14589     343    23  0.0399  0.4204
    70   0  1.053  5485712    0  14605     358    26  0.0380  0.3880
    71   0  1.046  5486322    0  14641     335    19  0.0399  0.3765
    72   0  1.045  5486155    0  14630     328    10  0.0404  0.3985
    73   0  1.057  5485987    0  14597     330    15  0.0403  0.3967
    74   0  1.043  5485091    0  14612     362    22  0.0391  0.4528
    75   0  1.045  5485713    0  14585     345    15  0.0402  0.4380
    76   0  1.051  5485450    0  14539     350    33  0.0388  0.3827
    77   0  1.041  5485261    0  14561     371    28  0.0405  0.4078
    78   0  1.067  5485516    0  14597     347    14  0.0401  0.4367

  1476 OUT OF   1626 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.24
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.39
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.039   186.645    72.217  90.000 100.208  90.000
 REC. CELL PARAMETERS   0.016378  0.005358  0.014070  90.000  79.792  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.075    18.526    19.005
 COORDINATES OF UNIT CELL B-AXIS   -14.518  -110.014   150.074
 COORDINATES OF UNIT CELL C-AXIS    18.331   -57.170   -40.136
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.014     0.065    -0.100
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1476 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  447    131    100     71     77    138    190    157    150
 SIGMAB (degree)       0.069  0.072  0.072  0.072  0.072  0.070  0.068  0.070  0.072
 SIGMAR (degree)       0.399  0.327  0.328  0.362  0.386  0.390  0.389  0.396  0.377


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1476 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  1 -1  0  0  0     0  0 -1  0  0  0  0  0  0
     0  0  0  1  0  0  0 -1  0     0  0  0  0  0  0  0  0  0     0 -1 -1  0  0  0  0 -1  0
     0  0 -1  1  1  1  1  0  0     0  0  0  1  1  1  0  0  0     0 -1  0  3  2  2  1  0  0
     0 -1 -1  2  2  1  1  1  0     0 -1  1  3  3  2  2  1  0     1  0  3  6  7  7  4  0 -1
     0 -1  0  0  0  1  1  0  0     0  0  1  1  1  4  3  1  0     1  1  5  5  9 14  7  1  0
     0  0  0  2  2  2  0  1  0     0  0  1  2  5  5  2  0  0     0  1  4  6 13 13  5  0  0
     0  0  1  1  1  0 -1  0  0     0  0  1  2  3  1  0  0  0    -1  0  2  3  5  3  0 -1  0
     0  1  0  1  0  0  0  0  0     0  0  0  1  1 -1 -1  0  0     0 -1  0  1  1 -1 -2 -1  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0 -1  0 -1  0  0  0


     0  0 -2 -2 -1  0  0  0  0     0  0 -2 -1 -1  0 -2  0  0     0  0 -2 -2  0  0  0  0  0
     0 -1 -1  0  1  1  0 -1  0     0 -1  0  0  1  2  1 -1  0     0 -1 -1 -1  1  1  0 -1  0
     0 -1  1  6  7  7  3  0  0    -1  0  2  4  8  8  4  1  0    -1  0  2  6  9  6  2  0  0
     0  0  4 11 20 17  7  1 -1    -1  1  5 16 38 22  7  2  0    -1  1  6 16 25 15  6  1  0
     0  1  7 19 41 31 10  2 -1     0  2  7 37100 38  9  2  0     0  2  8 24 39 23  8  2  0
     0  1  7 15 28 19  7  1 -1     0  2  6 24 43 18  7  2 -1    -1  2  7 15 22 16  6  0 -1
    -1  0  2  7  9  6  2  0  0     0  1  3  6  7  4  2  0  0    -1  0  3  5  7  6  1 -1  0
     0  0  0  0  0 -1 -1 -2  0     0  0  1  2  1  0 -1 -2  0     0 -1  0  1  1  0 -1 -1  0
     0  0 -1 -1 -1 -2 -2  0  0     0  0 -1 -1  0 -1 -2  0  0     0  0 -1 -1 -1 -1 -2  0  0


     0  0 -1 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1 -1 -1  0  0  0 -1  0     0 -1  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0
     0  0  0  3  5  3  2  0  0     0  1  0  1  3  2  1  0  0     0  0  0  0  1  1  0  0  0
     0  1  5 11 10  6  3  0  0     0  1  2  4  3  2  2  0  0     0  0  0  1  1  1  0  1  0
     0  1  7 11  8  7  4  0  0     0  0  3  3  1  2  1  0  0     0  0  1  0  0  1  1 -1  0
    -1  0  4  7  8  9  4  0  0     0  0  2  3  3  2  0  0  0     0  0  1  1  2  1  0  0  0
    -1 -1  1  2  4  3  0 -1  0     0  0  1  2  2  1  0  0  0     0  0  1  2  1  0  0  0  0
     0 -1  0  0  1  0 -1 -1  0     0 -1  0  1  1  0  0  0  0     0  0  0  1  0  0  0  1  0
     0  0 -1 -1  0 -1  0  0  0     0  0  0  0  1  0  0  0  0     0  0  0  0  0  0  0  0  0

  30448 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES       79 ...     104
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
    79   0  1.058  5486071    0  14554     330    19  0.0401  0.3942
    80   0  1.046  5484468    0  14581     347    22  0.0409  0.4026
    81   0  1.048  5485191    0  14608     350    21  0.0418  0.4100
    82   0  1.055  5485176    0  14611     334    18  0.0422  0.4410
    83   0  1.055  5484774    0  14596     333    29  0.0431  0.3684
    84   0  1.049  5485657    0  14599     323    23  0.0406  0.4142
    85   0  1.047  5486348    0  14580     321    32  0.0410  0.3485
    86   0  1.066  5486328    0  14560     316    30  0.0423  0.3558
    87   0  1.046  5485958    0  14558     315    23  0.0422  0.3807
    88   0  1.060  5486101    0  14579     321    24  0.0416  0.4025
    89   0  1.039  5486360    0  14572     317    17  0.0418  0.4353
    90   0  1.056  5486271    0  14543     328    22  0.0408  0.4215
    91   0  1.063  5486477    0  14586     300    18  0.0440  0.4609
    92   0  1.048  5485620    0  14602     324    22  0.0442  0.4476
    93   0  1.045  5485886    0  14647     321    31  0.0424  0.3834
    94   0  1.045  5485213    0  14629     317    23  0.0432  0.3745
    95   0  1.064  5485734    0  14629     308    23  0.0410  0.3556
    96   0  1.061  5485023    0  14636     327    34  0.0430  0.3110
    97   0  1.046  5485813    0  14589     322    34  0.0414  0.3139
    98   0  1.051  5485558    0  14599     307    30  0.0407  0.3153
    99   0  1.060  5485910    0  14674     298    15  0.0419  0.3721
   100   0  1.047  5486363    0  14700     309    27  0.0416  0.3534
   101   0  1.050  5486243    0  14671     301    27  0.0440  0.3651
   102   0  1.028  5485813    0  14624     313    26  0.0422  0.3668
   103   0  1.055  5486782    0  14632     309    23  0.0424  0.3679
   104   0  1.053  5486636    0  14583     303    19  0.0420  0.3947

  1384 OUT OF   1494 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.23
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.33
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.056   186.423    72.162  90.000 100.210  90.000
 REC. CELL PARAMETERS   0.016374  0.005364  0.014081  90.000  79.790  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.075    18.550    19.037
 COORDINATES OF UNIT CELL B-AXIS   -14.561  -109.844   149.920
 COORDINATES OF UNIT CELL C-AXIS    18.351   -57.134   -40.079
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.007     0.045    -0.116
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1384 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  433    121     93     67     78    133    176    140    135
 SIGMAB (degree)       0.070  0.073  0.073  0.072  0.072  0.071  0.070  0.071  0.073
 SIGMAR (degree)       0.385  0.313  0.316  0.360  0.381  0.369  0.370  0.381  0.358


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1384 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -1  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0 -1  0     0 -1 -1  0  0  0 -1 -1  0
     0  0  0  1  2  1  1  0  0     0  0  0  2  2  2  1  0  0     0 -1  0  3  4  2  0 -1  0
     0  0  0  1  1  1  2  1  0     0  0  1  3  3  3  3  1  0     1  0  4  7  8  8  4  0 -1
     0  0  0  0  0  1  1  1  0     0  1  2  0  1  3  3  1  0     0  1  5  5  8 12  6  1 -1
     0  0  0  1  1  1  1  0  0     0  0  1  2  3  4  2  1  0     0  0  3  6 11 11  5  1  0
     0  0  1  1  1  0  0  0  0     0  0  1  2  3  1  0  0  0    -1  0  1  4  6  4  0 -1  0
     0  0  0  1  1  0  0  0  0     0  0  0  1  1  0 -1 -1  0     0 -1  0  1  1 -1 -2 -1  0
     0  0  0  0  0  0  0  0  0     0  0  0  2  0 -1  0  0  0     0  0 -1  0  0 -2  0  0  0


     0  0 -2 -1  0  0  0  0  0     0  0 -2 -2 -1 -1 -2  0  0     0  0 -2 -1 -1 -1  0  0  0
     0 -1 -1  0  1  1  0 -2  0     0 -2  0  1  2  2  1 -1  0     0 -2 -1 -1  1  1  0 -1  0
    -1 -1  1  7  7  6  4  0  0    -1  0  3  5  8  8  4  0  0    -2  0  3  8  9  6  2  0  0
     0  0  5 13 21 17  7  1 -1    -1  1  5 18 41 23  7  1 -1     0  1  6 16 24 14  6  0 -1
     0  1  8 19 40 31  9  2 -1     0  2  7 38100 40  9  2 -1    -1  2  7 24 37 22  7  1  0
     0  1  6 15 27 19  7  1 -2     0  2  6 25 44 20  8  2 -1    -1  1  7 15 21 16  6  1  0
    -1  0  2  6  9  8  3 -1  0    -1  0  4  7  8  5  3  0  0    -2  0  3  6  7  6  2 -1  0
     0 -1  0  0  0 -1 -1 -2  0     0 -1  1  2  1  0  0 -1  0     0 -1  0  1  0  0 -1 -1  0
     0  0 -1 -1 -1 -2 -2  0  0     0  0 -1 -1 -1 -2 -3  0  0     0  0 -1 -1 -1 -1 -1  0  0


     0  0 -2 -1 -1  0  0  0  0     0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1 -1 -1  0  1  0 -1  0     0  0 -1  0  0  1  0  0  0     0  0  0 -1  1  1  0  0  0
     0  0  0  4  5  3  1  0  0     0  0  0  2  3  2  1  0  0     0  0  0  1  2  1  0  0  0
     0  1  5 10 10  6  3  1  0     0  1  3  5  3  2  2  0  0     0  1  0  2  1  0  0  0  0
     0  1  7 10  7  6  4  0  0     0  1  4  3  1  1  1  0  0     0  1  2  1 -1  1  0  0  0
     0  0  4  6  8  7  3  0  0     0  0  2  3  3  2  1  0  0     0  0  2  2  1  1  0  0  0
     0 -1  0  2  4  3  0 -1  0     0  0  0  1  2  1  0  0  0     0  0  1  1  1  0  0  0  0
     0 -1 -1  0  1  0 -1 -1  0     0 -1  0  0  0  0  0  0  0     0  0  0  1  0  1  1  1  0
     0  0  0  0  1 -1  0  0  0     0  0  0  0  1  0  0  0  0     0  0  0  0  0  0  0  0  0

  30365 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      105 ...     129
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   105   0  1.053  5487537    0  14608     297    29  0.0414  0.3200
   106   0  1.049  5487073    0  14654     301    38  0.0422  0.2815
   107   0  1.045  5487178    0  14597     303    32  0.0411  0.3110
   108   0  1.055  5487320    0  14585     291    26  0.0438  0.3222
   109   0  1.050  5487899    0  14584     290    34  0.0424  0.3199
   110   0  1.048  5488352    0  14610     272    17  0.0439  0.3848
   111   0  1.053  5487471    0  14585     282    19  0.0443  0.4311
   112   0  1.059  5488048    0  14610     276    15  0.0443  0.4476
   113   0  1.050  5487866    0  14631     272    16  0.0464  0.4443
   114   0  1.043  5486825    0  14572     286    33  0.0450  0.3514
   115   0  1.054  5486702    0  14583     296    23  0.0443  0.4082
   116   0  1.058  5487051    0  14630     282    23  0.0446  0.3719
   117   0  1.050  5486533    0  14591     283    18  0.0432  0.3775
   118   0  1.052  5487544    0  14618     265    17  0.0448  0.3656
   119   0  1.048  5487540    0  14530     269    18  0.0440  0.3610
   120   0  1.050  5488853    0  14488     274    29  0.0435  0.3268
   121   0  1.046  5488126    0  14541     282    39  0.0441  0.2959
   122   0  1.034  5487607    0  14547     287    31  0.0436  0.3289
   123   0  1.034  5487144    0  14558     289    24  0.0448  0.3540
   124   0  1.048  5486698    0  14660     283    39  0.0453  0.2957
   125   0  1.052  5487116    0  14653     273    27  0.0437  0.3371
   126   0  1.037  5487608    0  14473     276    19  0.0458  0.4018
   127   0  1.046  5487685    0  14526     272    23  0.0474  0.3698
   128   0  1.042  5488310    0  14531     264    17  0.0458  0.3844
   129   0  1.052  5487495    0  14514     269    19  0.0462  0.3605

  1197 OUT OF   1289 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.23
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.32
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.082   186.280    72.135  90.000 100.251  90.000
 REC. CELL PARAMETERS   0.016369  0.005368  0.014088  90.000  79.749  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.075    18.588    19.085
 COORDINATES OF UNIT CELL B-AXIS   -14.570  -109.813   149.763
 COORDINATES OF UNIT CELL C-AXIS    18.357   -57.097   -40.080
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.022     0.033    -0.134
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1197 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  378    105     81     59     72    122    141    113    113
 SIGMAB (degree)       0.070  0.075  0.073  0.072  0.072  0.072  0.071  0.072  0.074
 SIGMAR (degree)       0.374  0.297  0.310  0.357  0.369  0.360  0.365  0.373  0.343


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1197 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -2  0  0 -1  0  0  0
     0  0  0  0  0  0  0 -1  0     0  0 -1  0  0  0  0 -1  0     0 -1 -1  0  0  0 -1 -1  0
     0  0  0  1  1  1  1  0  0     0  0  0  2  2  2  0 -1  0     0  0  0  4  4  3  0  0  0
     0  0  0  0  1  1  1  1  0     0  0  1  1  2  2  2  0  0     0  1  4  5  7  7  4  0  0
     0 -1  0 -1 -1  0  1  0  0     0  1  1  0  0  2  2  1  0     0  1  4  4  6 10  6  1  0
     0  1  0  0  1  1  0  1  0     0  0  1  2  2  3  1  0  0     0  0  3  5 10 10  4  0  0
     0  0  0  1  1  1  0  0  0     0  0  1  2  4  2  0  0  0    -1  0  1  4  6  4  0 -1  0
     0  1  1  1  1  0  0  0  0     0 -1  0  1  1  0 -1 -1  0     0 -1  0  1  0 -1 -2 -1  0
     0  0  0  0  0  0  0  0  0     0  0  0  1  0  0  0  0  0     0  0 -1 -1 -1 -1  0  0  0


     0  0 -1  0 -1 -2  0  0  0     0  0 -2 -2 -1 -2  0  0  0     0  0 -3 -2  0 -1  0  0  0
     0 -1 -1  0  1  0  0 -2  0     0 -2 -1  1  2  2  0 -2  0     0 -2 -1  0  1  1  0 -1  0
    -1 -1  1  7  8  6  3 -1  0    -2  0  4  7  9  8  4  0  0    -1  0  4  9  9  6  2 -1  0
    -1  0  5 14 22 16  7  1 -2    -1  1  6 21 42 22  7  1 -1    -2  1  7 17 25 14  6  0 -1
    -1  1  8 20 38 30  9  2 -1    -1  1  7 39100 41  9  2 -1    -1  2  7 24 39 22  7  1 -1
     0  0  6 15 26 20  7  1 -2    -1  1  6 25 45 23  8  2 -2    -1  1  6 14 20 15  6  0  0
    -1  0  2  6  9  8  3 -1  0    -2  0  3  8  9  6  4  0  0     0  0  2  5  7  6  2 -1  0
     0 -1  0  0  0 -1 -2 -2  0     0 -1  0  2  2  1 -1 -2  0     0 -2  0  0  1  0 -1 -2  0
     0  0 -1 -1  0 -2  0  0  0     0  0 -2 -2 -1 -2 -2  0  0     0  0 -1 -1 -1 -1  0  0  0


     0  0 -2 -2 -1 -1  0  0  0     0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1 -1  0  1  0  0 -1  0     0  0 -1  0  0  1  0  0  0     0  0  0  0  0  1  1  0  0
     0  0  0  3  5  3  1 -1  0     0  0  0  2  3  3  1  0  0     0  0  0  1  2  1  1  1  0
     0  1  5 10  9  6  3  1  0     0  1  3  5  2  2  1  0  0     0  1  1  1  0  0  0  0  0
     0  1  7 10  7  5  4  0  1     0  1  4  3  0  0  1  0  0     0  1  2  0 -1  0 -1  0  0
     0  0  3  6  8  6  3  0 -1     0  0  2  3  3  2  1  0  0     0  0  2  2  1  0  0  0  0
     0 -1  0  1  4  3  0 -1  0     0  0  0  1  3  1  0  0  0     0  0  0  1  1  0  0  0  0
     0 -1 -1  0  0  0 -1 -1  0     0 -1 -1  0  1  0 -1  0  0     0 -1  0  0  1  0  0  0  0
     0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  29718 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      130 ...     155
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   130   0  1.033  5487277    0  14568     269    10  0.0461  0.4575
   131   0  1.035  5486914    0  14627     272    15  0.0447  0.4166
   132   0  1.039  5486498    0  14625     277    12  0.0446  0.4095
   133   0  1.046  5485532    0  14623     306    15  0.0439  0.4259
   134   0  1.042  5487435    0  14597     278    13  0.0443  0.4161
   135   0  1.059  5487935    0  14684     273    19  0.0445  0.3859
   136   0  1.043  5488583    0  14659     254    24  0.0458  0.3418
   137   0  1.029  5487552    0  14657     252    27  0.0460  0.3318
   138   0  1.041  5488021    0  14696     266    29  0.0450  0.3561
   139   0  1.037  5487904    0  14722     244    29  0.0433  0.3075
   140   0  1.054  5487567    0  14744     267    17  0.0439  0.4121
   141   0  1.041  5487060    0  14770     266    22  0.0444  0.4196
   142   0  1.054  5487143    0  14793     268    17  0.0468  0.4812
   143   0  1.034  5486229    0  14724     266    19  0.0472  0.4509
   144   0  1.042  5486661    0  14657     273    24  0.0458  0.3944
   145   0  1.040  5486420    0  14646     271    22  0.0459  0.3834
   146   0  1.041  5486235    0  14640     289    39  0.0459  0.3201
   147   0  1.049  5487233    0  14637     266    28  0.0454  0.3222
   148   0  1.047  5486821    0  14647     272    14  0.0456  0.4140
   149   0  1.040  5487250    0  14686     263    12  0.0442  0.4348
   150   0  1.027  5487766    0  14621     261    13  0.0464  0.4541
   151   0  1.035  5487664    0  14597     246    25  0.0456  0.3402
   152   0  1.055  5487424    0  14652     251    15  0.0472  0.4189
   153   0  1.039  5488107    0  14683     242    11  0.0459  0.4388
   154   0  1.033  5488217    0  14661     241    15  0.0459  0.4022
   155   0  1.059  5488629    0  14699     245    13  0.0461  0.4359

  1164 OUT OF   1255 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.25
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.29
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.075   186.512    72.290  90.000 100.258  90.000
 REC. CELL PARAMETERS   0.016371  0.005362  0.014058  90.000  79.742  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.073    18.599    19.058
 COORDINATES OF UNIT CELL B-AXIS   -14.358  -110.196   149.791
 COORDINATES OF UNIT CELL C-AXIS    18.392   -57.144   -40.276
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.100     0.001    -0.090
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1164 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  375    102     83     60     69    118    139    106    105
 SIGMAB (degree)       0.071  0.076  0.074  0.073  0.073  0.072  0.071  0.072  0.075
 SIGMAR (degree)       0.373  0.295  0.310  0.356  0.368  0.356  0.363  0.372  0.339


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1164 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -2 -1  0 -1  0  0  0
     0  0  0  0  1  0  0  0  0     0 -1  0  1  0  0  0 -1  0     0 -1 -1  0  1  0 -1 -1  0
     0  0  0  1  2  1  1  0  0     0  0  1  2  3  2  0  0  0     0  0  1  4  4  2  0 -1  0
     0  0  1  0  0  0  1  1  0     0  0  2  2  1  2  2  0  0     0  1  4  5  7  7  4  0  0
     0  0 -1 -2 -1  0  0  1  0     0  1  1 -1  0  2  2  1  0     1  1  4  3  6 11  6  1  0
     0  0  0  0  1  1  0  0  0     0  1  1  1  3  4  2  0  0     0  1  3  5 10 11  5  0  0
     0  0  0  1  2  2  0  0  0     0  0  1  3  5  2  0 -1  0    -1 -1  1  4  6  3 -1 -1  0
     0  0  0  1  2  1 -1 -1  0     0 -1  0  1  2  0 -2 -1  0     0 -2 -1  0  1 -1 -2 -2  0
     0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0  0  0     0  0 -1  0  0 -1  0  0  0


     0  0 -2 -2 -1 -2  0  0  0     0  0 -2 -2 -2 -2  0  0  0     0  0 -2 -1  0 -1  0  0  0
     0 -1 -1  0  1  1 -1 -3  0     0 -2 -1  1  2  2  0 -2  0     0 -2 -1  0  1  1  0 -2  0
    -1 -1  2  8  8  7  4 -1  0    -3  0  4  8 10  9  6  0  0    -1 -1  4  9 10  6  2 -1  0
    -1  1  6 13 21 17  8  1 -1    -1  2  6 20 42 24  9  2 -2    -2  1  6 16 25 15  7  0 -1
     0  1  8 19 39 31  9  1 -2    -1  2  7 35100 43 10  2 -2    -1  2  7 22 40 24  8  1  0
     0  0  6 13 25 19  8  1 -2    -2  1  6 25 48 25  9  1 -2     0  1  7 14 23 17  6  0 -1
    -2 -1  2  5  9  8  3 -1  0    -2  0  5 10 11  9  5 -1  0     0 -1  3  6  8  6  2 -1  0
     0 -1  0  0  0 -1 -1 -2  0     0 -2  0  2  2  1 -1 -2  0     0 -1 -1  1  1 -1 -1 -2  0
     0  0 -2  0 -1 -2  0  0  0     0  0  0 -2 -2 -2 -2  0  0     0  0  0 -2 -2 -1  0  0  0


     0  0 -2 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1 -1 -1  0  0 -1 -1  0     0 -1 -1  0  0  1  0 -1  0     0  0  0  0  1  1  0  1  0
     0  0  1  4  5  4  1  0  0     0  0  0  2  3  3  2  0  0     0  0  0  2  3  1  1  1  0
     0  1  5  9  9  5  4  1  0     0  1  3  4  3  1  1  1  0     0  0  1  1  0 -1 -1  0  0
     0  0  7 10  7  5  4  1  1     0  1  3  3  0  0  1  0  0     0  1  2  0 -1 -1  0  0  0
     0 -1  4  6  8  7  3  0  0     0  0  1  3  3  1  1  0  0     0  1  2  1  1  0  0  0  0
     0 -1  0  2  4  3  0 -1  0     0  0  0  2  3  2  0 -1  0     0 -1  0  1  2  1  0  0  0
     0 -1  0  0  1  0 -1 -2  0     0 -1  0  0  1  0  0 -1  0     0 -1  0  1  1  1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30517 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      156 ...     181
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   156   0  1.047  5488267    0  14657     254    24  0.0460  0.3552
   157   0  1.047  5488938    0  14636     263    18  0.0452  0.3994
   158   0  1.029  5489621    0  14602     245    19  0.0464  0.3851
   159   0  1.034  5489027    0  14599     256    24  0.0462  0.3450
   160   0  1.050  5488438    0  14519     258    28  0.0469  0.3265
   161   0  1.032  5488275    0  14478     264    25  0.0467  0.3665
   162   0  1.030  5487875    0  14596     264    23  0.0455  0.3771
   163   0  1.028  5488353    0  14645     255    24  0.0453  0.3588
   164   0  1.040  5488650    0  14691     243    36  0.0456  0.2857
   165   0  1.050  5489430    0  14636     224    29  0.0479  0.3085
   166   0  1.023  5489253    0  14663     217    22  0.0467  0.3654
   167   0  1.047  5489244    0  14612     234    24  0.0467  0.3585
   168   0  1.043  5490196    0  14582     229    31  0.0462  0.3006
   169   0  1.035  5489560    0  14558     238    25  0.0462  0.3398
   170   0  1.033  5489333    0  14462     249    40  0.0457  0.2496
   171   0  1.021  5488819    0  14474     248    34  0.0456  0.2888
   172   0  1.049  5489562    0  14534     244    32  0.0443  0.3053
   173   0  1.032  5489482    0  14645     225    20  0.0469  0.3416
   174   0  1.029  5489843    0  14723     219    14  0.0454  0.3848
   175   0  1.023  5489305    0  14698     230    16  0.0486  0.4115
   176   0  1.041  5489613    0  14686     232    26  0.0467  0.3273
   177   0  1.043  5489758    0  14623     231    28  0.0466  0.3311
   178   0  1.035  5489284    0  14574     247    43  0.0454  0.2669
   179   0  1.025  5489972    0  14486     238    21  0.0455  0.3543
   180   0  1.036  5489485    0  14495     246    35  0.0457  0.3039
   181   0  1.042  5490265    0  14486     230    30  0.0442  0.3305

  1080 OUT OF   1168 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.27
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.29
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.089   186.294    72.217  90.000 100.278  90.000
 REC. CELL PARAMETERS   0.016368  0.005368  0.014073  90.000  79.722  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.069    18.630    19.086
 COORDINATES OF UNIT CELL B-AXIS   -14.370  -110.044   149.631
 COORDINATES OF UNIT CELL C-AXIS    18.391   -57.091   -40.221
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.096    -0.009    -0.103
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1080 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  356     86     67     54     73    115    127     97     94
 SIGMAB (degree)       0.071  0.076  0.074  0.072  0.073  0.073  0.072  0.072  0.075
 SIGMAR (degree)       0.368  0.295  0.317  0.358  0.364  0.348  0.360  0.368  0.332


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1080 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0 -1  0  0  0  0     0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  1  0 -1  0     0 -1  0  0  0  0 -1 -2  0     0 -2 -1  0  1  0 -1 -3  0
     0  0  0  1  1  2  2 -1  0     0  0  0  2  3  3  1 -1  0    -1 -1  1  4  5  3  1 -1  0
     0 -1  0 -1  0  1  2  1  0     0  0  1  0  1  3  3  0  0     0  1  3  4  7  8  5  0  0
     0  0 -1 -1 -1 -1  1  2  0     0  1  0 -1  0  2  3  2  0     1  1  3  3  6 11  7  0  0
     0  1  0  0  1  1  1  0  0     0  1  1  1  3  3  2  1  0    -1  1  2  5 10  9  5  0  0
     0  0  1  1  2  2  0 -1  0     0  0  1  2  3  2 -1  0  0    -1  0  1  4  5  4  0 -1  0
     0  0  0  1  2  0 -1 -1  0     0 -1  0  1  1  0 -1 -1  0     0 -1  0  1  0 -1 -1 -2  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -2  0  0  0     0  0 -1 -1 -1 -2  0  0  0


     0  0 -2 -1 -1 -2  0  0  0     0  0 -1 -3 -2 -3  0  0  0     0  0 -3 -2 -2 -1  0  0  0
     0 -2 -1  0  1  0 -2 -3  0     0 -2 -1  1  2  1 -1 -3  0     0 -2 -2  0  1  1  0 -2  0
    -1 -1  3  8  9  6  3 -2  0    -2 -1  5 10 11  8  5 -1  0     0 -2  5  9 10  6  2  0  0
    -1  1  5 13 23 17  8  1  0    -1  1  6 21 43 22  8  2 -2    -2  0  7 18 26 14  6  0 -2
    -1  1  7 19 38 30  9  2 -2    -1  2  7 37100 43 10  2 -2    -2  1  7 23 41 23  7  1 -1
    -2  0  5 13 25 19  8  1 -2    -2  1  7 27 48 28 10  2 -3    -2  1  6 14 21 15  6  0 -2
    -1  0  2  5  9  9  4 -2  0    -3  0  5 10 12 11  6 -1  0     0 -1  3  6  6  6  2  0  0
     0 -2  0  1  1  0 -2 -3  0     0 -2 -1  1  2  1 -2 -3  0     0 -2 -1  1  0 -1 -1 -2  0
     0  0 -2 -1 -1 -3  0  0  0     0  0  0 -2 -2 -2  0  0  0     0  0  0 -2 -1 -2  0  0  0


     0  0 -2 -2 -1 -1  0  0  0     0  0  0 -1  0  1  0  0  0     0  0  0  0  0  0  0  0  0
     0 -2 -1 -1  0  0  0 -1  0     0 -1 -2 -1  1  1  0  0  0     0 -1 -1  0  0  1  1  0  0
     0  0  1  5  6  4  1  0  0     0  0  0  4  4  2  1  0  0     0  0  1  3  2  1  1  1  0
     0  0  6  9 10  6  3  0  0     0  1  4  4  2  1  1  1 -1     0  1  2  1  0  0  0  0  0
     0  0  7  9  7  5  4  1  0     0  0  4  2  0  0  0  0  0     0  2  3  0 -1 -1  0  0  0
     0  0  3  6  9  6  3  0  0     0  0  2  4  3  1  1  0  0     0  0  2  2  0 -1  0  0  0
     0 -1  0  2  4  3  0 -1  0     0 -1  0  2  3  2  0  0  0     0 -1  1  2  2  1  0  0  0
     0 -2  0  0  1  0 -1 -1  0     0 -1 -1  0  1  0 -1 -1  0     0 -1 -1  0  2  0  0 -1  0
     0  0  0 -1  0 -2  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30418 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      182 ...     207
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   182   0  1.030  5489955    0  14462     231    36  0.0470  0.2832
   183   0  1.031  5490189    0  14513     232    25  0.0460  0.3384
   184   0  1.027  5489185    0  14526     234    20  0.0479  0.3538
   185   0  1.041  5490524    0  14605     225    27  0.0478  0.3351
   186   0  1.032  5491150    0  14631     229    27  0.0462  0.3591
   187   0  1.037  5491691    0  14727     223    21  0.0459  0.3634
   188   0  1.031  5491850    0  14703     204    20  0.0475  0.3474
   189   0  1.040  5491957    0  14618     216    22  0.0477  0.3317
   190   0  1.038  5491785    0  14565     223    27  0.0474  0.3318
   191   0  1.029  5491141    0  14506     231    15  0.0461  0.4630
   192   0  1.031  5491284    0  14493     224    33  0.0466  0.2972
   193   0  1.029  5490771    0  14494     232    26  0.0466  0.3689
   194   0  1.028  5490509    0  14508     222    26  0.0464  0.3576
   195   0  1.022  5489672    0  14488     241    25  0.0451  0.3518
   196   0  1.030  5490710    0  14520     223    23  0.0476  0.3554
   197   0  1.031  5490310    0  14545     229    22  0.0462  0.3801
   198   0  1.016  5489751    0  14609     243    30  0.0459  0.3424
   199   0  1.017  5489905    0  14655     238    23  0.0480  0.3358
   200   0  1.038  5490536    0  14716     238    29  0.0462  0.2869
   201   0  1.020  5490221    0  14707     241    22  0.0447  0.3164
   202   0  1.038  5491226    0  14684     226    19  0.0460  0.3363
   203   0  1.015  5490934    0  14707     228    25  0.0447  0.2964
   204   0  1.025  5490981    0  14657     223    21  0.0454  0.3036
   205   0  1.015  5490419    0  14664     220    17  0.0452  0.3442
   206   0  1.012  5490112    0  14609     227    16  0.0479  0.3702
   207   0  1.011  5488868    0  14594     241    20  0.0469  0.3239

  1049 OUT OF   1114 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.35
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.26
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.102   186.267    72.220  90.000 100.281  90.000
 REC. CELL PARAMETERS   0.016365  0.005369  0.014073  90.000  79.719  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.081    18.643    19.081
 COORDINATES OF UNIT CELL B-AXIS   -14.304  -110.069   149.585
 COORDINATES OF UNIT CELL C-AXIS    18.392   -57.079   -40.242
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.111    -0.015    -0.091
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1049 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  353     88     67     57     72    105    113     93     98
 SIGMAB (degree)       0.072  0.077  0.075  0.073  0.073  0.073  0.073  0.073  0.076
 SIGMAR (degree)       0.360  0.300  0.325  0.362  0.365  0.338  0.352  0.356  0.325


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1049 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0 -1 -1 -1  0  0  0     0  0  0 -2 -1 -1  0  0  0
     0  0  0  1  1  0  0 -2  0     0 -2 -1  0  1  0 -1 -2  0     0 -1 -1  0  1  0 -2 -2  0
     0  0  0  1  2  2  2  0  0     0  0  1  3  4  2  1 -1  0     0  0  2  5  6  3  0 -1  0
     0  0  0 -1  0  0  2  1  0     0  1  1  1  2  2  3  1  0     0  2  4  5  8  7  4  0  0
     0  0 -1 -2 -2 -1  1  1  0     0  1  0 -1  0  1  3  1  0     0  2  4  4  5  8  6  1  0
     0  1  0  0  0  0  1  0  0     0  1  1  1  2  3  2  0  0    -2  0  3  5 10  8  4  1  0
     0  1  1  1  2  2  0 -1  0     0  0  1  2  3  2  0  0  0    -3 -1  1  4  6  3  1 -1  0
     0  0  0  1  1  0 -1 -1  0     0 -1  0  1  1 -1 -1 -1  0     0 -2 -1  1  0 -1 -2 -2  0
     0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0  0  0     0  0  0 -1  0 -2  0  0  0


     0  0 -2 -2 -3 -3  0  0  0     0  0 -2 -2 -2 -3  0  0  0     0  0 -2 -2 -2 -2  0  0  0
     0 -2 -1 -1  0  0 -2 -3  0     0 -3 -2  1  3  1 -2 -3  0     0 -3 -2  0  1  1  0 -3  0
     0  1  4  7  8  6  3 -2  0    -3 -2  7 12 12  8  5 -2  0     0 -2  4 10 10  6  2 -1  0
    -1  1  8 16 24 17  8  0  0    -1  1  8 23 42 22  9  1 -2    -2  0  6 17 27 14  6  0 -2
    -2  2  9 21 38 28 10  2 -2    -2  2  8 38100 44 11  2 -2    -3  2  8 25 43 22  8  2 -1
    -2  0  6 14 26 19  8  1  0    -2  1  8 29 49 30 11  1 -2     0  0  7 16 23 15  7  1 -2
    -3 -1  2  5  9  9  5 -2  0    -2 -2  5 12 13 13  7 -2  0    -2 -1  3  6  6  5  2  0  0
     0 -3 -1  0  0 -1 -2 -3  0     0 -2 -1  1  2  0 -2 -2  0     0 -1 -1  0  0 -1 -2 -2  0
     0  0  0 -1 -2 -2  0  0  0     0  0  0 -1 -2 -3  0  0  0     0  0  0 -3 -2 -2  0  0  0


     0  0  0 -3 -1 -1  0  0  0     0  0  0 -2 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0
     0 -3 -2 -1  1  0 -1 -2  0     0 -2 -2  0  1  1 -1 -1  0     0 -1 -2  0  1  1  1  0  0
     0 -2  1  5  6  4  2 -1  0     0  0  1  4  4  3  1  0  0     0  0  1  3  2  0  0  0  0
     0  1  6 10  9  6  4  1 -1     0  1  5  4  2  1  1  1  0     0  1  3  1 -1 -1 -1  0  0
     0  0  7  9  7  4  4  1  0     0  1  4  2  0  0  0  0  0     0  1  1 -1 -2 -1 -2 -1  0
     0  0  4  7  9  6  3  1 -1     0  0  2  3  3  2  1  0  0     0  0  1  1  0 -1  0  0  0
     0 -1  0  3  4  3  0 -1  0     0 -1  0  2  4  2  0 -1  0     0 -1  1  2  2  1  1  1  0
     0 -2 -1  0  0 -1 -2 -2  0     0 -1 -1  0  1  0 -1 -1  0     0 -1  0  0  2  1  0  0  0
     0  0  0 -2 -2 -2  0  0  0     0  0  0  0 -2 -1  0  0  0     0  0  0  0  0  0  0  0  0

  30336 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      208 ...     233
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   208   0  1.018  5489214    0  14571     238    25  0.0457  0.3059
   209   0  1.019  5489489    0  14583     227    38  0.0458  0.2446
   210   0  1.014  5488591    0  14525     253    27  0.0468  0.3608
   211   0  1.021  5488457    0  14565     261    23  0.0462  0.3648
   212   0  1.021  5488492    0  14584     263    27  0.0478  0.3795
   213   0  1.026  5488556    0  14613     260    27  0.0459  0.3743
   214   0  1.020  5488861    0  14661     255    38  0.0479  0.3115
   215   0  1.022  5488633    0  14698     254    36  0.0455  0.3488
   216   0  1.027  5490136    0  14695     229    39  0.0470  0.3239
   217   0  1.023  5490298    0  14715     235    39  0.0449  0.3425
   218   0  1.010  5489290    0  14709     238    40  0.0460  0.2926
   219   0  1.008  5488862    0  14690     250    36  0.0471  0.3282
   220   0  1.016  5489314    0  14677     243    35  0.0472  0.3077
   221   0  1.017  5489487    0  14673     243    34  0.0486  0.3178
   222   0  1.018  5490458    0  14684     243    32  0.0462  0.2976
   223   0  1.004  5489621    0  14703     234    25  0.0485  0.3483
   224   0  1.014  5490070    0  14710     239    31  0.0471  0.3527
   225   0  1.011  5489835    0  14674     231    27  0.0496  0.3070
   226   0  1.011  5489849    0  14679     222    22  0.0484  0.3409
   227   0  1.004  5488850    0  14636     249    22  0.0479  0.3449
   228   0  1.014  5489528    0  14620     237    24  0.0481  0.3211
   229   0  1.012  5488930    0  14558     251    36  0.0472  0.3189
   230   0  1.011  5489445    0  14518     242    33  0.0464  0.3241
   231   0  1.014  5489077    0  14515     257    40  0.0481  0.2992
   232   0  1.015  5489874    0  14548     248    44  0.0461  0.2955
   233   0  1.021  5489787    0  14554     253    45  0.0446  0.3154

  1063 OUT OF   1143 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.38
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.30
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.104   186.330    72.270  90.000 100.288  90.000
 REC. CELL PARAMETERS   0.016365  0.005367  0.014063  90.000  79.712  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.086    18.646    19.069
 COORDINATES OF UNIT CELL B-AXIS   -14.234  -110.170   149.596
 COORDINATES OF UNIT CELL C-AXIS    18.393   -57.101   -40.302
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.134    -0.016    -0.077
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1063 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  361     85     69     59     77    111    110     91     93
 SIGMAB (degree)       0.072  0.076  0.074  0.073  0.073  0.074  0.073  0.073  0.076
 SIGMAR (degree)       0.352  0.299  0.326  0.365  0.364  0.345  0.357  0.360  0.328


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1063 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0 -2 -2  0  0  0  0
     0  0  0  1  1  1  0 -1  0     0 -2 -1  1  2  0 -1 -3  0     0 -3 -2  0  1  0 -2 -3  0
     0  1  1  0  1  1  1 -1  0     0  0  1  3  3  2  1 -2  0     0 -1  2  5  7  4  0 -2  0
     0  0  0 -1 -1  0  2  1  0     0  1  0  0  1  2  2  0  0     0  1  3  5  8  7  5  0  0
     0  0 -1 -2 -2 -1  1  1  0     0  0  0 -1 -1  1  3  1  0     0  2  3  4  6  9  6  1  0
     0  1  0 -1  0  0  2  0  0     0  1  0  1  3  3  3  0  0     0  1  3  6 10  8  4  0  0
     0  1  1  1  2  2  1  0  0     0  0  1  2  4  2  0 -1  0     0 -1  1  4  6  3  1 -3  0
     0  0  1  1  1  0 -1 -2  0     0 -1  0  1  2  0 -2 -2  0     0 -2 -1  1  1 -1 -3 -3  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0 -1 -2  0  0  0


     0  0 -2 -3 -4  0  0  0  0     0  0 -2 -3 -3 -4  0  0  0     0  0 -3 -3 -2 -2  0  0  0
     0 -3 -2 -1  1  0 -2 -4  0     0 -3 -2  0  3  0 -3 -3  0     0 -3 -2 -1  1  1 -2 -4  0
     0 -1  4  7 10  6  3 -3  0     0 -2  6 13 14  9  5 -3  0     0 -2  5  8  9  6  3 -2  0
     0  1  7 16 25 17  7  0  0    -2  1  7 22 44 22  8  0 -3     0 -1  7 18 28 16  6  1 -3
    -2  2  7 22 40 29 10  2 -2    -3  2  8 38100 45 11  3 -3     0  1  8 26 45 24  9  2 -2
    -2  0  5 13 26 19  9 -1 -3    -2  0  7 30 49 31 11  1 -3     0  0  6 18 24 15  7  1 -2
     0 -2  2  5  7  8  5 -3  0    -2 -2  5 13 15 14  8 -2  0     0 -1  3  7  5  5  3 -1  0
     0 -3 -1  1  0 -1 -2 -3  0     0 -2 -2  0  3 -1 -2 -2  0     0 -2 -2  0  0 -1 -2 -3  0
     0  0  0 -2 -2 -3  0  0  0     0  0  0 -3 -3 -2  0  0  0     0  0  0 -3 -2 -3  0  0  0


     0  0  0 -3 -1  0  0  0  0     0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -3 -2 -1  1  0 -1 -3  0     0 -2 -2 -1  1  0  0 -2  0     0 -2 -1  1  1  1  1 -1  0
     0 -2  2  5  6  4  1 -1  0     0 -1  2  4  4  2  1  0  0     0  0  2  3  2  1  0  0  0
     0  0  6 10  9  7  4  1  0     0  0  5  3  2  1  1  1  0     0  1  2  0  0 -1 -1 -1  0
     0  1  5 10  8  5  4  2  0     0  0  3  2  0  0  0  0  0     0  0  1 -1 -1 -2 -1 -1  0
     0  0  3  7 10  6  3  0  0     0 -1  2  3  4  2  1  0  0     0  0  1  1  1 -1 -1  0  0
     0 -1  1  3  5  3  1 -1  0     0 -1  0  3  4  3  1 -1  0     0 -1  0  2  2  1  1  0  0
     0 -1 -1  0  0 -1 -2 -2  0     0 -2 -1  0  1  0 -1 -2  0     0 -2 -1  1  1  1  0 -1  0
     0  0  0 -2 -2  0  0  0  0     0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0

  30345 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      234 ...     258
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   234   0  1.015  5490215    0  14587     234    34  0.0456  0.3400
   235   0  1.006  5490783    0  14573     236    36  0.0462  0.3167
   236   0  1.012  5489792    0  14626     243    21  0.0487  0.3820
   237   0  1.027  5490568    0  14579     233    23  0.0497  0.3207
   238   0  1.000  5490268    0  14627     245    20  0.0507  0.3333
   239   0  1.010  5489629    0  14592     242    32  0.0497  0.2699
   240   0  0.999  5490638    0  14627     235    31  0.0493  0.2639
   241   0  1.022  5492291    0  14630     219    26  0.0498  0.3177
   242   0  1.012  5492164    0  14623     211    18  0.0491  0.3886
   243   0  1.004  5493235    0  14583     213    11  0.0477  0.4962
   244   0  1.002  5493557    0  14591     191     8  0.0481  0.5656
   245   0  1.005  5493182    0  14597     196    11  0.0474  0.5298
   246   0  1.013  5492937    0  14577     197    11  0.0490  0.5262
   247   0  0.992  5491818    0  14627     199    18  0.0479  0.4408
   248   0  0.999  5491372    0  14584     204    18  0.0475  0.4556
   249   0  1.013  5491284    0  14593     206    26  0.0468  0.3647
   250   0  1.015  5491067    0  14584     188    26  0.0484  0.3086
   251   0  1.016  5490587    0  14601     216    23  0.0454  0.2821
   252   0  0.992  5490523    0  14615     205    20  0.0476  0.3256
   253   0  0.995  5490402    0  14597     222    20  0.0462  0.3469
   254   0  1.007  5490869    0  14603     212    18  0.0478  0.3506
   255   0  0.991  5490841    0  14622     222    26  0.0478  0.3206
   256   0  1.003  5490614    0  14621     215    24  0.0488  0.3627
   257   0  1.000  5491365    0  14625     209    37  0.0472  0.2907
   258   0  0.993  5490730    0  14592     209    31  0.0476  0.3804

   957 OUT OF   1034 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.45
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.29
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.070   186.425    72.256  90.000 100.216  90.000
 REC. CELL PARAMETERS   0.016370  0.005364  0.014063  90.000  79.784  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.083    18.608    19.003
 COORDINATES OF UNIT CELL B-AXIS   -14.160  -110.200   149.699
 COORDINATES OF UNIT CELL C-AXIS    18.409   -57.090   -40.285
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.125    -0.031    -0.042
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     957 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  323     78     61     52     70     96     99     85     89
 SIGMAB (degree)       0.073  0.077  0.075  0.074  0.074  0.074  0.073  0.073  0.076
 SIGMAR (degree)       0.378  0.299  0.333  0.370  0.375  0.355  0.362  0.369  0.333


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    957 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0 -2  0  0  0  0
     0  0  0  1  2  1 -1 -2  0     0 -2 -1  1  2  0 -2 -3  0     0 -3 -2  1  1 -1 -3 -4  0
     0  0 -1  0  1  1  0 -1  0     0 -1  1  2  4  2  0 -2  0     0  0  3  6  7  3  1 -3  0
     0 -1 -2 -1 -1  0  1  1  0     0  1  0  0  2  3  3  1  0     0  2  4  6  8  8  5 -1  0
     0  0 -1 -2 -1 -1  0  1  0     0  1  0 -1  0  2  3  2  0     0  2  4  5  8 10  7  2  0
     0  0  0  0  1  1  1  0  0     0  1  1  2  4  4  3  0  0     0  0  3  7 12  9  5  0  0
     0  0  1  2  3  2  1 -1  0     0  0  2  3  5  2  0 -2  0     0 -1  1  5  6  3  0 -3  0
     0  0  1  2  2  1 -2 -3  0     0 -2  0  1  2 -1 -3 -3  0     0 -3 -1  1  1 -1 -3 -3  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0 -2 -3  0  0  0


     0  0 -3 -3 -3 -3  0  0  0     0  0 -3 -2 -3 -4  0  0  0     0  0 -3  0 -2 -3  0  0  0
     0 -3 -2  0  1 -1 -3 -4  0     0 -3 -2  0  3 -1 -3 -4  0     0 -4 -2 -1  2  1 -3 -4  0
     0 -1  5  8  8  5  3 -3  0     0 -2  6 15 14  9  5 -3  0     0 -3  5  9 10  6  3 -2  0
    -3  1  7 15 24 17  7 -1  0    -2  0  7 22 43 21  8  0 -4     0 -1  8 19 29 17  7  1 -4
    -2  2  8 22 43 31 11  1 -3    -3  2  8 39100 47 11  3 -3     0  1  8 26 45 25  9  2 -3
     0  1  5 15 28 21  9 -1  0    -1  0  9 32 49 36 13  0 -3     0 -1  6 18 23 17  8  1  0
     0 -2  3  5  7  8  5 -3  0     0 -2  7 13 15 18 10 -2  0     0 -2  4  7  5  5  4 -1  0
     0 -3 -2  0  0 -2 -3 -3  0     0 -2 -2  0  2 -1 -3 -3  0     0 -2 -2 -1  0 -2 -2 -3  0
     0  0  0 -2 -3 -4  0  0  0     0  0  0 -2 -3 -2  0  0  0     0  0  0 -2 -3 -3  0  0  0


     0  0  0  0 -2 -3  0  0  0     0  0  0  0 -2  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -3 -3 -1  2  1 -1 -4  0     0 -3 -2  1  1  1  0 -2  0     0 -2 -1  2  1  1  1  0  0
     0 -2  3  5  7  5  2 -1  0     0 -2  3  4  4  3  2  0  0     0  0  3  2  2  1  1  1  0
     0 -1  6 10 10  7  5  1  0     0  0  4  3  2  1  1  2  0     0  3  2 -1  0 -1 -1  0  0
     0  0  6 10  8  4  4  2  0     0  1  4  3  0  0  0  0  0     0  3  2 -1 -2 -2 -1 -1  0
     0 -1  4  8 10  5  2  1  0     0 -1  3  5  4  2  0 -1  0     0  1  2  1  0  0 -1 -1  0
     0 -2  1  3  5  3  1 -2  0     0 -2  0  4  4  2  0 -2  0     0 -2  1  3  2  1  0  0  0
     0 -3 -2 -1  0 -1 -2 -2  0     0 -3 -2  1  1  0 -1 -2  0     0 -3 -1  2  2  1  0 -1  0
     0  0  0 -2  0  0  0  0  0     0  0  0  0 -2  0  0  0  0     0  0  0  0  0  0  0  0  0

  29816 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      259 ...     284
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   259   0  1.006  5489837    0  14653     213    28  0.0498  0.3823
   260   0  1.012  5489828    0  14635     213    23  0.0474  0.4521
   261   0  1.014  5489301    0  14650     208    15  0.0505  0.4771
   262   0  1.024  5489630    0  14679     213    17  0.0490  0.4784
   263   0  0.999  5489856    0  14646     207    12  0.0508  0.5231
   264   0  1.007  5488960    0  14674     220    21  0.0499  0.4233
   265   0  0.997  5488128    0  14694     241    24  0.0500  0.4398
   266   0  1.017  5488025    0  14671     229    17  0.0508  0.4678
   267   0  1.005  5487195    0  14670     247    25  0.0497  0.4203
   268   0  0.994  5487679    0  14670     247    23  0.0508  0.4249
   269   0  1.003  5486990    0  14660     252    20  0.0493  0.4153
   270   0  1.005  5487395    0  14663     239    20  0.0507  0.4128
   271   0  1.010  5488289    0  14625     240    20  0.0498  0.3780
   272   0  1.013  5487536    0  14688     241    22  0.0509  0.4012
   273   0  0.990  5486786    0  14677     244    23  0.0511  0.4065
   274   0  1.008  5487082    0  14662     242    26  0.0508  0.3925
   275   0  0.991  5485374    0  14688     269    21  0.0510  0.5066
   276   0  1.010  5484854    0  14642     255    16  0.0501  0.5049
   277   0  0.991  5484604    0  14646     258    17  0.0507  0.4686
   278   0  1.019  5485588    0  14677     255    20  0.0512  0.4867
   279   0  1.006  5484818    0  14630     245    24  0.0499  0.4171
   280   0  1.019  5485332    0  14662     244    24  0.0483  0.4433
   281   0  0.972  5485254    0  14606     241    17  0.0493  0.4671
   282   0  0.998  5485304    0  14624     242    19  0.0484  0.4389
   283   0  0.984  5486221    0  14617     238    19  0.0482  0.4504
   284   0  1.010  5486304    0  14641     224    24  0.0487  0.3902

  1016 OUT OF   1088 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.47
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.33
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.055   186.730    72.395  90.000 100.214  90.000
 REC. CELL PARAMETERS   0.016374  0.005355  0.014036  90.000  79.786  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.098    18.562    18.956
 COORDINATES OF UNIT CELL B-AXIS   -14.145  -110.382   149.946
 COORDINATES OF UNIT CELL C-AXIS    18.374   -57.208   -40.380
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.130     0.016    -0.028
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1016 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  348     84     63     55     77     97     98     90     97
 SIGMAB (degree)       0.073  0.078  0.076  0.074  0.074  0.074  0.073  0.074  0.077
 SIGMAR (degree)       0.381  0.307  0.345  0.373  0.372  0.356  0.367  0.363  0.329


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1016 STRONG REFLECTIONS *****


     1  1  1  1  1  1  1  1  1     1  1  1  1  0  1  1  1  1     1  1  1 -1 -2  1  1  1  1
     1  1  2  2  3  3  1 -1  1     1 -2  1  3  4  3 -1 -2  1     1 -2  0  3  4  1 -2 -3  1
     1  2  2  2  3  4  4 -1  1     1  2  3  4  5  5  3 -1  1     1  0  4  6  9  6  2 -3  1
     1  1  0  0  1  2  4  3  1     1  3  2  2  4  6  5  2  1     1  4  7  8 10 11  7  0  1
     1  2  0 -1 -1  0  2  4  1     1  3  2  0  2  3  6  3  1     1  4  6  7  9 12  9  1  1
     1  2  1  0  2  2  3  3  1     1  2  3  3  5  5  5  1  1     1  2  4  7 12 10  7  1  1
     1  3  2  2  3  3  4  1  1     1  2  2  4  5  4  4 -1  1     1  0  2  5  7  5  4 -1  1
     1  2  2  3  3  2  0 -1  1     1 -1  0  3  3  1 -1 -2  1     1 -2 -1  1  2  0 -1 -2  1
     1  1  1  1  1  1  1  1  1     1  1  1  1  1  1  1  1  1     1  1  1  1 -2  1  1  1  1


     1  1  1 -1 -3 -3  1  1  1     1  1 -1 -1 -2 -2  1  1  1     1  1 -2 -1 -2  1  1  1  1
     1 -1 -1  2  4  0 -2 -3  1     1 -1 -1  1  5  0 -2 -3  1     1 -2 -1  0  3  2 -2 -3  1
     1 -1  8 11 12  6  3 -2  1     1 -1 10 19 16  9  6 -2  1     1 -1  7 11 12  9  5 -2  1
     1  2  9 17 26 18  7 -1 -3     0  1 10 24 44 23  9  1 -2     1  0  9 21 32 20  9  2 -2
    -1  4  9 22 46 35 12  1 -3    -1  3  9 35100 53 13  3 -2     0  2  8 27 49 29 12  5 -2
    -2  2  7 14 29 26 13  0 -1     0  1  9 30 50 42 16  1 -2     1  0  8 18 28 21 10  3 -2
     1 -1  4  6  9 11  9 -1  1    -1 -1  7 14 15 20 12 -1  1     1  0  5  8  8  8  6  0  1
     1 -1 -1  1  2  0 -1 -2  1     1 -1 -1  0  3  0 -2 -2  1     1 -1 -1  0  1  0 -2 -3  1
     1  1  1 -1 -2 -2  1  1  1     1  1  1 -1  0 -2  1  1  1     1  1  1  1 -1 -2  1  1  1


     1  1  1 -2 -1  1  1  1  1     1  1  1  1  1  1  1  1  1     1  1  1  1  1  1  1  1  1
     1 -3 -2  0  4  3 -1 -3  1     1 -2 -1  2  4  4  1 -1  1     1 -2  0  3  4  3  2  2  1
     1 -1  4  7 10  8  4  0  1     1 -1  3  5  6  5  3  2  1     1 -1  3  3  3  2  1  3  1
     1  0  7 13 12  9  7  3  1     1  1  6  6  4  2  2  2  1     1  3  4  2  0  0  0  1  1
     1  2  6 12  9  6  6  5  1     1  3  6  4  1  0  1  2  1     1  5  4  1 -1 -1  0  0  1
     1  1  5  7 12  8  4  3  1     1  2  3  5  5  2  1  1  1     1  2  3  2  1  0  0  0  1
     1  0  3  3  7  6  3  0  1     1 -1  1  3  6  4  3  1  1     1 -1  2  4  3  2  1  2  1
     1 -1  0  0  2  2 -1 -2  1     1 -1 -1  1  3  3  1 -1  1     1 -1  0  2  3  3  2  1  1
     1  1  1  1 -2  1  1  1  1     1  1  1  1  0  1  1  1  1     1  1  1  1  1  1  1  1  1

  30530 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      285 ...     310
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   285   0  0.972  5486147    0  14642     213    37  0.0493  0.2772
   286   0  0.994  5485803    0  14635     231    36  0.0490  0.3171
   287   0  0.991  5486772    0  14651     223    37  0.0482  0.3327
   288   0  1.000  5486604    0  14673     209    27  0.0494  0.3965
   289   0  0.978  5486237    0  14665     211    17  0.0486  0.5098
   290   0  0.989  5487311    0  14687     195    16  0.0504  0.5350
   291   0  0.979  5486175    0  14686     214    14  0.0503  0.5651
   292   0  0.975  5487315    0  14703     214    18  0.0488  0.5224
   293   0  0.984  5487492    0  14710     217    18  0.0469  0.5180
   294   0  0.976  5486668    0  14692     217    20  0.0486  0.4648
   295   0  0.968  5486509    0  14660     225    31  0.0499  0.4462
   296   0  0.981  5486578    0  14652     242    29  0.0470  0.4178
   297   0  0.974  5486670    0  14666     226    19  0.0493  0.4226
   298   0  0.971  5487933    0  14675     221    15  0.0496  0.4570
   299   0  0.969  5488496    0  14710     216    14  0.0488  0.4910
   300   0  0.968  5489315    0  14706     221    21  0.0504  0.4678
   301   0  0.978  5489162    0  14741     221    26  0.0492  0.4301
   302   0  0.973  5489153    0  14715     222    27  0.0492  0.4825
   303   0  0.972  5489684    0  14700     215    11  0.0504  0.6707
   304   0  0.965  5488654    0  14700     221    13  0.0489  0.6386
   305   0  0.962  5488594    0  14673     212    15  0.0481  0.5054
   306   0  0.967  5487945    0  14710     220    20  0.0482  0.4305
   307   0  0.979  5487749    0  14636     228    28  0.0475  0.3776
   308   0  0.977  5486429    0  14663     240    28  0.0491  0.3760
   309   0  0.984  5486824    0  14629     246    29  0.0497  0.4010
   310   0  0.965  5485852    0  14613     247    20  0.0503  0.5030

   928 OUT OF   1012 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.42
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.39
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.003   186.694    72.341  90.000 100.078  90.000
 REC. CELL PARAMETERS   0.016381  0.005356  0.014040  90.000  79.922  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.092    18.472    18.889
 COORDINATES OF UNIT CELL B-AXIS   -14.123  -110.395   149.894
 COORDINATES OF UNIT CELL C-AXIS    18.415   -57.147   -40.353
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.132    -0.004    -0.006
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     928 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  323     78     61     53     76     88     78     77     85
 SIGMAB (degree)       0.073  0.078  0.075  0.073  0.074  0.074  0.072  0.074  0.077
 SIGMAR (degree)       0.397  0.315  0.355  0.389  0.383  0.367  0.380  0.371  0.337


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    928 STRONG REFLECTIONS *****


     1  1  1  1  1  1  1  1  1     1  1  1  1 -1 -1  1  1  1     1  1  1 -1 -1 -2  1  1  1
     1  0  1  3  3  2  0 -1  1     1 -1  1  3  3  1 -1 -2  1     1 -2  1  3  3  0 -1 -2  1
     1  1  2  3  4  4  2  0  1     1  1  3  4  6  5  2 -2  1     1  1  4  6  7  4  2 -2  1
     1  1  1  2  3  3  4  2  1     1  2  2  4  6  6  6  1  1     1  3  6  9 13 11  6  0  1
     1  1  0  0  0  0  3  3  1     1  3  3  1  2  4  6  2  1     1  4  7  8 12 14  9  2  1
     1  1  1  1  3  2  3  1  1     1  2  3  4  7  5  4  1  1     1  3  5  9 15 12  9  0  1
     1  2  4  3  5  4  1 -1  1     1  0  3  6  6  4  2 -1  1     1  0  3  6  7  6  4 -1  1
     1  0  2  4  3  2 -1 -1  1     1 -1  0  3  3  2  0 -1  1     1 -2  0  1  2  1 -1 -1  1
     1  1  1  1  1  1  1  1  1     1  1 -2  1  1  1  1  1  1     1  1  1 -2  1  1  1  1  1


     1  1  1 -1 -2 -1 -2  1  1     1  1 -2 -1  0 -2 -3  1  1     1  1  1 -2 -1 -3  1  1  1
     1 -2 -1  4  4  0 -2 -2  1     1 -1 -1  2  6  0 -2 -3  1     1 -1 -1  1  3  1 -1 -3  1
     1  0  7 10 11  6  3 -2  1     1 -1  8 18 18 10  6 -2  1     1 -1  7 12 10  7  4 -1  1
     1  2  9 17 26 16  7  0 -3    -1  1 10 23 44 21  8  1 -3     1  0 10 21 28 18  9  2 -3
     1  4 10 26 52 35 12  3  1    -1  4 10 38100 49 12  3 -1     1  3  9 27 48 28 12  4 -2
     1  3  8 18 31 26 14  1  1     1  2 10 30 50 38 14  1 -1     1  0  7 19 31 21 10  2 -2
     1  0  5  6  8 11  8 -1  1     0 -1  9 13 16 18 11 -1  1     0 -1  5  9 10 10  7 -1  1
     1 -2 -1  0  1  1 -1 -1  1     1 -1 -1  0  3  0 -1 -2  1     1 -1 -1  0  3  2 -1 -3  1
     1  1 -1 -2 -2 -1  1  1  1     1  1 -1 -2 -2 -2  1  1  1     1  1 -2 -1  0 -2  1  1  1


     1  1  1 -2 -1  1  1  1  1     1  1  1  1  1  1  1  1  1     1  1  1  1  1  1  1  1  1
     1 -2 -1  1  3  1 -1 -3  1     1 -2  0  3  4  3  1  0  1     1 -2  1  4  3  3  3  2  1
     1 -1  4  7  8  6  4  0  1     1 -1  3  5  6  5  4  2  1     1  0  3  3  3  2  2  3  1
     1  1  8 12 11  9  6  3  1     1  2  6  5  3  3  2  3  1     1  2  4  2  1  1  1  1  1
     1  2  7 12 10  6  5  4  1     1  3  6  5  2  0  1  2  1     1  3  4  1  0 -1  0  0  1
     1  0  5 10 14  8  5  3  1     1  1  4  6  6  3  1  2  1     1  2  4  3  2  1  0  0  1
     1 -1  2  5  8  6  3  0  1     1  0  1  4  7  5  2  0  1     1  0  2  4  5  3  1  1  1
     1 -1 -1  0  2  3  0 -3  1     1 -1 -1  0  4  4  1 -1  1     1 -1  0  2  4  4  2  1  1
     1  1 -2  1 -1  1  1  1  1     1  1  1  1  0 -1  1  1  1     1  1  1  1  1  1  1  1  1

  30547 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      311 ...     336
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   311   0  0.961  5485608    0  14622     236    28  0.0508  0.4068
   312   0  0.961  5485016    0  14597     259    27  0.0488  0.4690
   313   0  0.956  5484788    0  14634     266    42  0.0504  0.3282
   314   0  0.970  5485216    0  14636     258    41  0.0521  0.3669
   315   0  0.952  5484107    0  14616     258    39  0.0500  0.3657
   316   0  0.973  5483170    0  14622     262    32  0.0496  0.3817
   317   0  0.965  5482624    0  14640     261    33  0.0481  0.4026
   318   0  0.968  5481720    0  14638     254    36  0.0488  0.3612
   319   0  0.967  5482148    0  14669     255    30  0.0490  0.3632
   320   0  0.952  5482110    0  14726     254    20  0.0486  0.4029
   321   0  0.973  5483576    0  14718     251    31  0.0476  0.3760
   322   0  0.952  5484469    0  14748     244    12  0.0464  0.4821
   323   0  0.974  5486005    0  14720     233    13  0.0457  0.4726
   324   0  0.956  5485301    0  14759     240    13  0.0464  0.5113
   325   0  0.967  5486407    0  14706     235    18  0.0455  0.4680
   326   0  0.967  5486051    0  14679     239    21  0.0468  0.4629
   327   0  0.971  5486260    0  14682     216    32  0.0475  0.3552
   328   0  0.973  5485821    0  14697     219    27  0.0493  0.3939
   329   0  0.953  5486257    0  14756     216    24  0.0491  0.4710
   330   0  0.948  5485487    0  14776     221    19  0.0467  0.4752
   331   0  0.959  5485689    0  14749     217    22  0.0484  0.4589
   332   0  0.948  5485191    0  14701     216     8  0.0489  0.5967
   333   0  0.956  5484681    0  14639     224     7  0.0496  0.5824
   334   0  0.951  5485530    0  14623     223    20  0.0498  0.4461
   335   0  0.949  5484612    0  14584     226    18  0.0484  0.4862
   336   0  0.953  5484784    0  14540     229    21  0.0473  0.4998

   986 OUT OF   1065 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.38
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.41
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     61.967   186.756    72.461  90.000 100.023  90.000
 REC. CELL PARAMETERS   0.016388  0.005355  0.014014  90.000  79.977  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.092    18.385    18.859
 COORDINATES OF UNIT CELL B-AXIS   -14.243  -110.382   149.969
 COORDINATES OF UNIT CELL C-AXIS    18.420   -57.264   -40.399
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.128     0.030    -0.021
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM     986 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  348     80     59     57     81     99     87     81     89
 SIGMAB (degree)       0.071  0.076  0.074  0.072  0.073  0.073  0.072  0.073  0.076
 SIGMAR (degree)       0.411  0.331  0.376  0.408  0.395  0.382  0.394  0.389  0.351


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF    986 STRONG REFLECTIONS *****


     0  0  0 -1 -1  0  0  0  0     0  0  0 -2 -1 -1  0  0  0     0  0  0 -2 -1 -2  0  0  0
     0 -1 -1  1  1  0 -2 -3  0     0 -1 -1  0  1  0 -2 -4  0     0 -3  0  2  1 -1 -3 -3  0
     0  1  1  2  3  3  0  0  0     0  0  1  3  5  4  0 -1  0     0  0  2  4  6  3  1 -2  0
     0  1  0  0  2  3  4  1  0     0  1  2  3  6  7  6  0  0     0  1  6  9 13 11  5  0  0
     0  0  0 -2 -2  0  4  1  0     0  2  3  1  2  5  7  1  0    -1  2  6  9 13 16  9  1  0
     0  0  0  0  1  3  2  0  0     0  0  3  3  6  6  5  0  0     0  0  4  9 14 12  8  0  0
     0 -1  1  2  4  4  1 -1  0     0 -2  2  4  5  4  1 -2  0     0 -2  1  5  6  5  1 -2  0
     0 -2  0  1  2  1 -1 -4  0     0 -2 -1  0  1  0 -1 -3  0     0 -3 -2  0  1  1 -3 -4  0
     0  0  0  0  0  0  0  0  0     0  0  0 -1 -1  0  0  0  0     0  0  0 -2 -3 -3  0  0  0


     0  0  0 -3 -2 -3  0  0  0     0  0  0 -2 -3 -3  0  0  0     0  0  0 -3 -3 -3  0  0  0
     0 -3 -1  4  3 -1 -3 -4  0     0 -3 -2  3  5 -1 -3 -4  0     0 -3 -2  1  1 -1 -2 -4  0
     0 -1  4  9  9  5  2 -3  0     0 -2  7 17 16  8  3 -3  0     0 -2  5 11 10  5  2 -2  0
    -2  1  7 16 25 15  5 -1 -3    -2  1  7 21 42 20  6  0 -5     0  1  7 19 28 16  7  0 -4
    -1  2  8 25 55 32  9  1 -3    -3  2  8 37100 45  9  2 -3    -2  1  8 26 48 25 10  2 -4
    -2  1  6 16 29 23 10  0  0    -3  0  8 29 49 34 10  0 -3    -2  0  6 19 30 19  9  1 -4
     0 -1  2  5  8 10  6 -2  0    -2 -2  5 12 15 15  8 -2  0    -2 -2  2  7  8  8  5 -2  0
     0 -3 -2 -1  1  1 -3 -3  0     0 -2 -2  0  3  0 -2 -3  0     0 -3 -2 -1  1  2 -1 -3  0
     0  0 -2 -3 -4 -3  0  0  0     0  0 -2 -2 -3 -3  0  0  0     0  0 -4 -2 -3 -3  0  0  0


     0  0  0 -4 -2  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -3 -2  0  1  0 -2 -4  0     0 -3  0  0  2  1 -1 -2  0     0 -3  0  1  2  1  1  1  0
     0 -1  2  5  7  5  2 -1  0     0  0  1  3  5  3  3  0  0     0  1  1  3  3  2  1  2  0
     0  1  6 11 12  8  5  1 -3     0  1  4  5  3  3  3  2  0     0  1  2  1  0  0  0  1  0
     0  1  6 11 11  6  6  3 -2     0  0  4  4  1  1  2  2  0     0  0  3  0 -1 -2 -1  0  0
     0  0  4  9 14  8  5  1  0     0  0  3  5  5  2  1  0  0     0 -1  2  2  1  0 -1 -1  0
     0 -2  0  4  7  5  2 -1  0     0 -2 -1  4  6  4  2 -1  0     0 -1  0  4  4  2  1  0  0
     0 -3 -3 -1  1  2 -1 -4  0     0 -2 -3 -1  2  2  0 -1  0     0 -2 -2  0  3  2  1 -1  0
     0  0 -3 -2 -2 -2  0  0  0     0  0  0 -2 -1  0  0  0  0     0  0  0  0  0  0  0  0  0

  30417 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      337 ...     362
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   337   0  0.942  5484537    0  14573     241    17  0.0484  0.5145
   338   0  0.957  5484513    0  14549     235    15  0.0511  0.5236
   339   0  0.946  5483600    0  14575     259    12  0.0478  0.5307
   340   0  0.957  5483116    0  14655     259    12  0.0462  0.5823
   341   0  0.947  5482684    0  14706     271    20  0.0470  0.5423
   342   0  0.943  5482871    0  14733     256    18  0.0488  0.4871
   343   0  0.949  5483307    0  14725     250    21  0.0485  0.4382
   344   0  0.941  5483608    0  14771     249    17  0.0492  0.4860
   345   0  0.942  5483607    0  14713     271    17  0.0486  0.5347
   346   0  0.934  5484528    0  14684     241    21  0.0493  0.4747
   347   0  0.934  5484233    0  14690     246    16  0.0492  0.5220
   348   0  0.955  5484520    0  14711     242    18  0.0486  0.5279
   349   0  0.941  5484965    0  14651     242    13  0.0494  0.6115
   350   0  0.932  5483866    0  14699     244    13  0.0489  0.5707
   351   0  0.945  5484574    0  14678     238    13  0.0487  0.5958
   352   0  0.949  5484944    0  14607     242    19  0.0470  0.5569
   353   0  0.952  5484155    0  14567     218    23  0.0472  0.4942
   354   0  0.936  5482613    0  14592     250    20  0.0473  0.4984
   355   0  0.940  5482572    0  14581     252    21  0.0476  0.5007
   356   0  0.954  5481267    0  14638     254    27  0.0481  0.4566
   357   0  0.943  5480621    0  14698     262    22  0.0471  0.4775
   358   0  0.944  5480803    0  14775     261    25  0.0484  0.4683
   359   0  0.941  5481308    0  14812     248    24  0.0480  0.4760
   360   0  0.945  5482019    0  14791     258    20  0.0446  0.5120
   361   0  0.949  5481501    0  14868     256    22  0.0486  0.4900
   362   0  0.936  5481976    0  14835     243    14  0.0483  0.5674

  1020 OUT OF   1126 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.38
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.46
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     61.936   186.798    72.615  90.000  99.990  90.000
 REC. CELL PARAMETERS   0.016394  0.005353  0.013983  90.000  80.010  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.087    18.305    18.846
 COORDINATES OF UNIT CELL B-AXIS   -14.467  -110.224   150.116
 COORDINATES OF UNIT CELL C-AXIS    18.418   -57.452   -40.410
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.096     0.078    -0.055
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1020 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  359     76     57     59     87     99     92     87     92
 SIGMAB (degree)       0.070  0.075  0.073  0.070  0.072  0.072  0.072  0.073  0.074
 SIGMAR (degree)       0.422  0.345  0.388  0.419  0.415  0.409  0.411  0.395  0.355


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1020 STRONG REFLECTIONS *****


     0  0  0  0  0 -1  0  0  0     0  0  0  1  1 -2  0  0  0     0  0  0  1 -1 -3  0  0  0
     0 -1  0  1  2  0 -2 -2  0     0 -1  0  1  2  0 -2 -3  0     0 -2  0  2  1 -1 -1 -3  0
     0  1  2  3  4  3  0 -1  0     0  0  1  3  6  3  0 -1  0     0  0  2  4  6  3  1 -2  0
     0  1  1  1  2  4  3  0  0     0  1  3  4  7  8  4  0  0     0  1  7  9 14 11  4  0  0
     0  0  0 -2 -1  2  4  1  0     0  2  4  2  3  8  7  1  0     0  2  7 10 15 18  8  1  0
     0  0  1  1  3  5  3  1  0     0  1  4  5  7  7  4  0  0     0  1  5  9 15 14  8  1  0
     0 -1  2  3  6  5  1 -1  0     0 -1  1  4  6  4  0 -2  0     0 -1  0  4  7  6  2 -2  0
     0 -1  0  1  3  0 -1 -3  0     0 -2 -1  0  1  1 -1 -3  0     0 -2 -2  0  1  1 -1 -3  0
     0  0  0  0  0  0  0  0  0     0  0  0 -1  0 -1  0  0  0     0  0  0 -2 -2  0  0  0  0


     0  0  0 -1 -2 -2 -2  0  0     0  0  0 -2 -3 -2 -3  0  0     0  0  0 -4 -3 -2  0  0  0
     0 -2  1  3  2 -1 -2 -3  0     0 -3 -2  4  4  0 -2 -3  0     0 -3 -2  1  2  0 -2 -3  0
     0  0  4  9  8  4  3 -2  0     0 -2  7 14 13  6  4 -2  0     0 -2  5 11 10  8  3 -2  0
     0  0  7 15 25 13  5  0 -3     0  1  7 21 41 19  7  1 -3     0  0  7 20 28 18  7  1 -3
    -2  2  9 26 60 30  8  2 -2    -2  2  8 38100 42  8  1 -3    -2  3  8 28 49 25 10  1 -2
     0  2  7 16 29 23  9  0  0    -2  1  9 27 45 30  9  0 -3     0  1  7 18 30 17  8  0 -3
    -2  0  2  5  9 11  6 -2  0    -3 -2  6  8 11 13  7 -2  0     0 -1  3  7  8  8  4 -1  0
     0 -2 -2 -1  1  0 -2 -3  0     0 -2 -2  0  2  2 -2 -3  0     0 -2 -2 -1  1  1  0 -3  0
     0  0 -4 -2 -3 -3  0  0  0     0  0 -2 -3 -2 -2  0  0  0     0  0 -2 -3 -2 -3  0  0  0


     0  0  0  0 -1 -1  0  0  0     0  0  0  0  1  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -3 -1  1  1  0 -2 -3  0     0 -3 -1  1  2  1 -1 -2  0     0 -2 -1  1  2  1  0 -1  0
     0 -1  2  6  8  6  2 -1  0     0 -1  1  5  6  4  2 -1  0     0  1  2  4  4  3  1  0  0
     0  1  6 13 13 10  6  1 -1     0  1  5  6  5  4  3  1 -1     0  2  3  2  1  1  1  1  0
     0  2  7 13 11  8  6  2  0     0  2  5  4  2  1  2  1  0     0  2  2  0 -1 -1  0  0  0
     0  1  4  9 14  8  4  0 -1     0  0  3  5  5  2  1  0  0     0  0  2  2  2  1  1  0  0
     0 -1  0  4  8  5  1 -2  0     0 -1 -1  3  6  4  1 -1  0     0 -1 -1  4  4  2  0  0  0
     0 -3 -2 -1  1  2  0 -3  0     0 -2 -2 -1  3  2  0 -2  0     0 -2 -2  0  3  3  0 -1  0
     0  0 -3 -2 -1  0  0  0  0     0  0  0 -1  0  1  0  0  0     0  0  0  0  0  0  0  0  0

  30593 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      363 ...     387
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   363   0  0.946  5483373    0  14788     240    12  0.0476  0.5973
   364   0  0.941  5482890    0  14719     248     9  0.0460  0.6204
   365   0  0.948  5482657    0  14661     262    11  0.0452  0.5939
   366   0  0.941  5482688    0  14651     255    13  0.0474  0.6033
   367   0  0.923  5482814    0  14645     256    14  0.0441  0.5236
   368   0  0.943  5482961    0  14643     254    22  0.0434  0.5104
   369   0  0.931  5482144    0  14693     259    23  0.0435  0.4695
   370   0  0.945  5482197    0  14678     261    13  0.0442  0.5711
   371   0  0.935  5482323    0  14737     256    18  0.0466  0.5181
   372   0  0.936  5481781    0  14763     266    18  0.0429  0.5295
   373   0  0.946  5481277    0  14783     271    14  0.0434  0.5422
   374   0  0.935  5480207    0  14759     280    13  0.0451  0.6101
   375   0  0.936  5480996    0  14734     265    16  0.0445  0.5708
   376   0  0.928  5481252    0  14666     273    31  0.0456  0.4986
   377   0  0.943  5480624    0  14617     267    27  0.0478  0.4802
   378   0  0.941  5480649    0  14600     275    20  0.0471  0.5217
   379   0  0.932  5480105    0  14654     264    19  0.0463  0.5691
   380   0  0.929  5479868    0  14691     268    18  0.0453  0.5394
   381   0  0.942  5479960    0  14715     274    11  0.0446  0.6039
   382   0  0.943  5479636    0  14682     268    11  0.0446  0.6109
   383   0  0.950  5479260    0  14724     290    20  0.0440  0.5506
   384   0  0.935  5478536    0  14719     283    14  0.0464  0.5919
   385   0  0.939  5479328    0  14758     270    17  0.0444  0.5220
   386   0  0.939  5478879    0  14746     282    22  0.0443  0.5278
   387   0  0.936  5479459    0  14704     278    15  0.0453  0.5553

  1005 OUT OF   1134 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.37
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.46
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     61.917   186.929    72.633  90.000  99.962  90.000
 REC. CELL PARAMETERS   0.016398  0.005350  0.013979  90.000  80.038  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.079    18.260    18.854
 COORDINATES OF UNIT CELL B-AXIS   -14.595  -110.286   150.221
 COORDINATES OF UNIT CELL C-AXIS    18.426   -57.475   -40.406
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.079     0.092    -0.077
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1005 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  350     81     57     61     91     94     80     85     98
 SIGMAB (degree)       0.069  0.075  0.071  0.069  0.070  0.071  0.070  0.072  0.075
 SIGMAR (degree)       0.433  0.356  0.395  0.432  0.427  0.417  0.411  0.403  0.364


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1005 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  2  2 -2  0  0  0     0  0  0  2  1 -1  0  0  0
     0  0  1  2  3  0 -2 -2  0     0 -1  0  2  3 -1 -2 -2  0     0 -1  1  2  2 -1 -2 -3  0
     0  1  1  3  5  3 -1 -1  0     0  0  1  4  6  2 -1 -1  0     0 -1  1  5  7  2 -1 -1  0
     0  1  2  1  4  5  1  0  0     0  1  3  4  8  8  2  0  0     0  1  6 11 15  9  3  0  0
     0  0  2 -1  0  4  3  1  0     0  1  5  3  3  9  4  0  0     0  1  8 12 18 18  6  0  0
     0  0  3  1  3  4  2  0  0     0  0  5  6  7  7  3  0  0     0  1  5 10 16 16  6  0  0
     0  0  2  4  4  3  0 -1  0     0 -1  1  5  6  4  1  0  0    -1 -1  0  4  7  6  1 -1  0
     0 -2  0  1  2  2  0 -2  0     0 -2 -1  1  3  2  0 -2  0     0 -2 -2  0  2  2  0 -2  0
     0  0  0  0  0  0  0  0  0     0  0  0 -1  1  0  0  0  0     0  0 -3 -2  1 -1  0  0  0


     0  0  0 -1 -1 -2 -2  0  0     0  0  0 -1 -3 -2 -3  0  0     0  0  0 -2 -3 -2  0  0  0
     0 -2  1  3  2  0 -1 -2  0     0 -3 -1  7  4  0 -2 -3  0     0 -3 -2  3  2  0 -2 -4  0
     0  0  3 10  9  5  2 -2  0     0 -1  6 14 13  6  3 -2  0     0 -1  5 12 10  6  1 -2  0
    -1  1  7 18 27 12  4  1 -3    -2  1  7 22 41 19  6  1 -3     0  2  7 20 28 17  7  0 -3
    -1  1  9 28 62 27  6  2 -2    -2  2  8 40100 40  8  1 -3    -2  2  8 28 50 24 10  1 -3
    -2  1  6 16 29 23  8  1  0    -2  1  8 23 43 29  9  0  0    -2  1  6 18 31 17  9  0 -3
    -2  0  1  5 10 12  5 -1  0    -2 -1  4  7 10 13  7 -1  0    -3 -2  2  6  9  8  4  0  0
     0 -2 -1  0  1  2 -1 -3  0     0 -3 -1  1  3  4 -1 -3  0     0 -3 -2 -1  2  2  0 -2  0
     0  0 -3 -2 -2 -2  0  0  0     0  0 -1 -2 -3 -3  0  0  0     0  0 -1 -1 -1 -1  0  0  0


     0  0  0  0  0 -2  0  0  0     0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0
     0 -2 -1  1  2  0 -2 -3  0     0 -1 -1  0  1  0 -1 -2  0     0  0  0  1  1  0  0 -1  0
     0  1  2  4  6  6  1 -2  0     0  1  1  4  6  4  1  0  0     0  1  2  3  5  3  1  0  0
     0  1  6 13 13 10  6  1 -1     0  2  4  5  5  4  4  1  0     0  1  3  2  2  1  1  0  0
     0  1  7 13 12  8  6  1 -2     0  1  5  5  2  2  2  1  0     0  2  2  1 -1  0  2  1  0
     0  1  4 10 15  7  5  0 -1     0  0  3  7  6  3  2  0  0     0  0  1  4  2  1  1  0  0
    -1 -1 -1  4  8  5  1 -1  0     0 -1 -1  4  6  4  1 -1  0     0 -2 -1  5  5  3  1  0  0
     0 -2 -2 -1  2  2  0 -2  0     0 -2 -2  0  3  2  0 -1  0     0 -2 -2  1  4  2  0 -1  0
     0  0 -3 -1  1  1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  2  0  0  0  0

  30008 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      388 ...     413
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   388   0  0.930  5480043    0  14681     267    12  0.0449  0.5918
   389   0  0.932  5480349    0  14655     259    11  0.0455  0.6356
   390   0  0.947  5479975    0  14676     272    17  0.0438  0.5408
   391   0  0.933  5479764    0  14707     258    19  0.0445  0.5575
   392   0  0.940  5479734    0  14708     256    23  0.0453  0.5234
   393   0  0.931  5480062    0  14722     267    17  0.0433  0.5704
   394   0  0.932  5479840    0  14734     282    22  0.0415  0.5411
   395   0  0.931  5479430    0  14750     268    18  0.0429  0.5284
   396   0  0.931  5479530    0  14748     254    13  0.0438  0.5392
   397   0  0.929  5479184    0  14702     271    15  0.0435  0.5403
   398   0  0.925  5478624    0  14736     260    17  0.0460  0.5156
   399   0  0.928  5478462    0  14764     275    13  0.0462  0.5629
   400   0  0.933  5478493    0  14672     291    14  0.0442  0.5885
   401   0  0.947  5478278    0  14676     282    20  0.0428  0.5203
   402   0  0.929  5478399    0  14745     290    20  0.0446  0.5421
   403   0  0.937  5479131    0  14699     278    18  0.0423  0.5272
   404   0  0.940  5479010    0  14713     279    22  0.0441  0.5382
   405   0  0.938  5480713    0  14745     274    13  0.0438  0.6063
   406   0  0.926  5480915    0  14756     244    12  0.0440  0.5725
   407   0  0.934  5481056    0  14783     250     8  0.0436  0.6070
   408   0  0.937  5480795    0  14720     256    12  0.0450  0.5617
   409   0  0.941  5479639    0  14695     259    20  0.0442  0.4977
   410   0  0.934  5478441    0  14710     278    24  0.0420  0.5006
   411   0  0.928  5478622    0  14725     272    16  0.0429  0.5741
   412   0  0.931  5478718    0  14641     275    16  0.0436  0.5338
   413   0  0.943  5479499    0  14686     259    16  0.0434  0.5532

  1063 OUT OF   1197 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.30
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.49
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     61.962   187.102    72.497  90.000  99.998  90.000
 REC. CELL PARAMETERS   0.016388  0.005345  0.014006  90.000  80.002  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.086    18.356    18.888
 COORDINATES OF UNIT CELL B-AXIS   -14.528  -110.362   150.388
 COORDINATES OF UNIT CELL C-AXIS    18.445   -57.362   -40.314
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.048     0.040    -0.058
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1063 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  363     83     61     67    103    109     84     85    100
 SIGMAB (degree)       0.068  0.074  0.071  0.068  0.069  0.070  0.069  0.071  0.073
 SIGMAR (degree)       0.442  0.355  0.412  0.440  0.437  0.427  0.423  0.408  0.369


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1063 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0     0  0 -1 -1  0 -1  0  0  0
     0  0  0  2  3  1 -1 -1  0     0  0  0  1  2  0 -2 -1  0     0 -1  0  1  1 -1 -1 -1  0
     0  1  2  4  6  3  0  0  0     0  1  2  5  6  3  0  0  0    -1  0  2  5  6  2  1  0  0
     0  2  2  3  5  6  1  1  0     0  2  4  6  9  8  3  0  0     0  2  7 12 16  9  3  0 -2
     0  1  2  1  0  4  3  1  0     0  1  5  4  5  9  4  0  0     0  1  8 12 18 16  5  0  0
     0  1  2  3  3  4  2  0  0     2  1  4  5  7  7  2 -1  0     0  0  4  9 14 13  5  0  0
     0  0  1  2  3  3  1 -1  0     0 -1  1  3  4  3  1  0  0     0 -1 -1  4  7  6  2  0  0
     0 -1  0  1  2  2  0 -1  0     0 -2 -1  0  1  2  0 -1  0     0 -2 -2  0  2  2  1 -1  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -2  0  0  0  0  0  0


     0  0  0  0  0 -2 -2  0  0     0  0  0 -1 -2 -2 -2  0  0     0  0  0 -2  0 -2 -3  0  0
     0 -2  1  2  1 -1 -1 -2  0     0 -3  0  6  4  0 -2 -3  0     0 -3  0  2  2  0 -2 -3  0
    -2  0  3  9  8  5  3 -1  0     0  0  5 11 12  6  3 -1  0     0 -1  4 10  9  5  0 -1  0
    -1  1  7 18 27 12  5  0 -3    -2  2  6 19 39 17  5  1 -2    -3  1  6 17 25 14  5  0 -1
     0  2  8 27 60 26  6  1 -3    -2  2  8 38100 38  7  2 -2    -2  2  8 29 51 26  8  1 -2
    -1  0  6 16 32 21  7  1 -3    -2  1  8 23 42 27  7  1 -2    -2  1  7 18 29 19  9  0 -1
    -2 -1  1  6 11 11  5 -1  0    -2 -1  4  6 10 12  6  0  0    -3 -1  3  6  9  9  3  0  0
     0 -2 -2  0  2  3  0 -2  0     0 -2 -1  1  3  4  0 -2  0     0 -2 -1 -1  2  2  0 -2  0
     0  0 -3 -1  0 -2  0  0  0     0  0 -2 -2 -2 -2  0  0  0     0  0 -2 -1  0  1  0  0  0


     0  0  0  0  1  0  0  0  0     0  0  0  0  0  1  0  0  0     0  0  0  0  0  0  0  0  0
     0 -2  0  0  1  0 -1 -2  0     0 -2 -1  1  1  1 -1 -2  0     0  0  0  1  1  1  1 -1  0
     0  0  1  4  6  4  0 -2  0     0  0  1  4  6  5  1 -1  0     0  0  1  3  5  2  1  1  0
     0  0  5 12 13 10  5  0 -2     0  1  4  7  7  5  4  1  0     0  1  2  3  3  2  0  1  0
    -1  1  6 14 13 10  7  1  0     0  1  4  6  3  3  4  1  0     0  2  2  2 -1 -1  2  0  0
    -1  1  4 11 14  9  6  1  0     0  1  2  6  6  3  3  1  0     0  1  0  4  2  1  1  1  0
    -1  0  0  3  8  5  1  0  0     0  0 -1  3  7  4  2  0  0     0  0 -1  3  5  3  1  1  0
     0 -2 -2 -1  2  1  0 -2  0     0 -1 -1  0  3  2  0 -1  0     0 -2 -1  1  4  3  1  0  0
     0  0 -2 -2  1 -1  0  0  0     0  0 -2 -1  1  0  0  0  0     0  0  0  0  0  0  0  0  0

  30718 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      414 ...     439
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   414   0  0.923  5479250    0  14738     267    21  0.0430  0.4810
   415   0  0.924  5479217    0  14723     271    24  0.0428  0.4970
   416   0  0.924  5479450    0  14760     259    16  0.0399  0.5321
   417   0  0.939  5479802    0  14755     261    22  0.0413  0.4620
   418   0  0.935  5480018    0  14731     258    15  0.0406  0.5563
   419   0  0.923  5479609    0  14715     264    16  0.0412  0.5541
   420   0  0.936  5479813    0  14739     264    20  0.0427  0.4800
   421   0  0.930  5479922    0  14722     267    15  0.0432  0.5654
   422   0  0.938  5480574    0  14712     272    21  0.0426  0.5409
   423   0  0.927  5479582    0  14679     291    12  0.0436  0.5834
   424   0  0.933  5479723    0  14741     271    12  0.0425  0.6721
   425   0  0.944  5480462    0  14794     273    23  0.0435  0.5568
   426   0  0.935  5480228    0  14766     274    18  0.0432  0.5716
   427   0  0.943  5480894    0  14781     277    26  0.0441  0.4808
   428   0  0.943  5481576    0  14769     257    16  0.0432  0.4743
   429   0  0.949  5481501    0  14747     277    15  0.0415  0.5648
   430   0  0.943  5482354    0  14729     266    12  0.0407  0.5574
   431   0  0.937  5480713    0  14698     296    12  0.0401  0.5873
   432   0  0.929  5480703    0  14716     279    13  0.0418  0.4906
   433   0  0.940  5480975    0  14752     274    15  0.0414  0.4983
   434   0  0.939  5480164    0  14723     295    14  0.0427  0.5649
   435   0  0.941  5480831    0  14713     294    10  0.0405  0.6141
   436   0  0.929  5480132    0  14682     285    11  0.0403  0.5668
   437   0  0.940  5479815    0  14724     288    24  0.0411  0.4745
   438   0  0.929  5478880    0  14696     297    29  0.0405  0.4482
   439   0  0.946  5479669    0  14690     293    27  0.0420  0.4433

  1096 OUT OF   1231 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.33
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.48
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     61.963   187.218    72.643  90.000 100.023  90.000
 REC. CELL PARAMETERS   0.016389  0.005341  0.013979  90.000  79.977  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.076    18.357    18.918
 COORDINATES OF UNIT CELL B-AXIS   -14.723  -110.231   150.608
 COORDINATES OF UNIT CELL C-AXIS    18.465   -57.543   -40.311
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION    -0.010     0.050    -0.089
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1096 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  372     79     59     70    113    115     85     91    103
 SIGMAB (degree)       0.067  0.072  0.069  0.067  0.068  0.069  0.069  0.070  0.072
 SIGMAR (degree)       0.449  0.374  0.425  0.447  0.451  0.444  0.441  0.418  0.380


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1096 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  1  0 -1  0  0  0     0  0  1  1 -1 -2 -2  0  0
     0  0  1  2  3  0 -2 -1  0     0 -1  0  2  1 -1 -2 -1  0     0 -1  1  2  1 -1 -2 -1  0
     0  0  2  4  5  2 -1  0  0     0  0  2  4  6  2  0  0  0    -1  0  2  6  6  2  0  0  0
     0  1  2  3  6  5  1  1  0     0  1  3  6 11  8  2  1 -1     0  1  5 11 16 10  4  0 -1
     0  1  1  2  1  5  3  1  0     0  1  4  5  6 10  4  1  0     0  1  7 13 20 17  5  0 -1
     0  0  0  4  4  5  3  1  0     0  0  3  6  7  7  3  0  0     0  0  4 10 14 13  4 -1  0
     0 -1  0  3  3  2  1  0  0     0 -1  1  4  5  3  0 -1  0    -1 -1  0  3  6  6  2  0  0
     0 -1  0  1  3  1  0 -1  0     0 -2 -1  0  3  2  0 -1  0     0 -3 -2 -1  2  2  0 -1  0
     0  0  0  0  0  0  0  0  0     0  0  0 -1  0  1  0  0  0     0  0 -4 -2  0  1  0  0  0


     0  0  0  1  1 -2 -2  0  0     0  0  0 -3 -2 -2 -3  0  0     0  0  0 -1  0 -1 -3  0  0
     0 -1  1  2  2 -1 -1 -2  0     0 -2  2  6  4  0 -2 -3  0     0 -2  0  3  2  0 -2 -3  0
    -1  0  3  8  8  4  2  0  0     0  0  6 11 11  4  3 -1  0     0  0  4 10  9  5  0 -2  0
     0  0  7 18 28 13  5  1 -3    -1  2  6 22 42 18  5  1 -2    -2  1  6 19 28 15  5 -1 -1
     0  1  7 27 61 26  6  1 -2    -1  2  7 39100 37  6  1 -2    -1  2  8 29 56 28  8  0 -1
    -1  0  5 15 29 20  6  1 -2    -1  2  7 18 35 24  7  2 -1    -1  2  6 15 28 18  9  0  0
    -2 -1  1  5  9  9  4  0  0    -2  0  3  4  8 10  5  1  0    -2  0  2  4  7  8  3  0  0
     0 -2 -1 -1  1  2  0 -2  0     0 -2 -1  0  3  4  1 -2  0     0 -2 -1 -1  1  2  1 -1  0
     0  0 -3 -2  0 -1  0  0  0     0  0 -2 -2 -2 -3  0  0  0     0  0 -2  0  0  0  0  0  0


     0  0  0  1 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  1  1  0 -2 -3  0     0  0  0  0  0  0 -1 -2  0     0  0  0  1  1  1  0 -1  0
     0 -1  2  4  6  4  0 -2  0     0  0  1  4  5  3  1 -1  0     0  1  1  3  4  2  1  0  0
     0  0  5 13 15 10  5  0 -1     0  1  4  7  8  6  3  0  0     0  2  4  4  3  2  0  0  0
    -1  1  7 15 16 10  7  1  0     0  2  4  6  4  3  3  0  0     0  2  1  1  0  1  1  0  0
     0  1  4 10 15  9  5  1  0     0  1  1  6  6  4  2  1  0     0  1 -1  2  3  2  1  0  0
    -1 -1 -1  2  7  5  2  0  0     0 -1 -1  2  6  4  2  1  0     0  0 -1  2  4  3  2  1  0
     0 -2 -2 -1  1  1  0 -1  0     0 -2 -2  0  3  2  0  0  0     0 -1 -2  0  2  2  1  1  0
     0  0 -2 -1  0  0  0  0  0     0  0 -3 -2  2  1  0  0  0     0  0  0  0  0  0  0  0  0

  30692 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      440 ...     465
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   440   0  0.939  5480243    0  14685     286    32  0.0400  0.3781
   441   0  0.935  5479793    0  14705     283    24  0.0402  0.4365
   442   0  0.957  5480607    0  14709     290    24  0.0381  0.4938
   443   0  0.929  5479621    0  14691     291    17  0.0401  0.5210
   444   0  0.933  5480671    0  14712     283    19  0.0388  0.5016
   445   0  0.945  5480601    0  14711     286    24  0.0399  0.5024
   446   0  0.942  5480138    0  14708     275    22  0.0405  0.5071
   447   0  0.952  5480197    0  14670     291    22  0.0380  0.5056
   448   0  0.942  5480907    0  14677     287    21  0.0379  0.5130
   449   0  0.947  5480932    0  14693     267    20  0.0391  0.4541
   450   0  0.948  5480554    0  14700     288    17  0.0394  0.4403
   451   0  0.938  5480287    0  14738     293    23  0.0411  0.4669
   452   0  0.951  5479454    0  14739     308    14  0.0386  0.5536
   453   0  0.945  5479174    0  14696     313    16  0.0398  0.5522
   454   0  0.944  5479691    0  14680     303    16  0.0409  0.5371
   455   0  0.953  5479561    0  14731     310    26  0.0402  0.5185
   456   0  0.944  5479898    0  14712     295    27  0.0395  0.4892
   457   0  0.951  5480477    0  14691     279    18  0.0398  0.5219
   458   0  0.936  5480638    0  14679     287    21  0.0391  0.4493
   459   0  0.944  5479008    0  14690     283    30  0.0384  0.4295
   460   0  0.939  5478939    0  14711     290    28  0.0375  0.4636
   461   0  0.943  5480251    0  14710     296    22  0.0382  0.4958
   462   0  0.945  5480104    0  14743     313    14  0.0400  0.5825
   463   0  0.961  5479987    0  14772     309    18  0.0393  0.5566
   464   0  0.952  5479638    0  14712     302    19  0.0401  0.5173
   465   0  0.950  5479203    0  14719     315    24  0.0390  0.4756

  1191 OUT OF   1313 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.28
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.49
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     61.969   187.464    72.625  90.000 100.033  90.000
 REC. CELL PARAMETERS   0.016388  0.005334  0.013983  90.000  79.967  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.077    18.361    18.931
 COORDINATES OF UNIT CELL B-AXIS   -14.839  -110.253   150.886
 COORDINATES OF UNIT CELL C-AXIS    18.454   -57.566   -40.249
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.040     0.063    -0.101
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1191 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  400     89     74     84    126    115     83     95    110
 SIGMAB (degree)       0.066  0.071  0.069  0.067  0.066  0.067  0.067  0.069  0.070
 SIGMAR (degree)       0.456  0.401  0.435  0.451  0.451  0.455  0.453  0.429  0.400


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1191 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  1 -1  0  0  0     0  0 -1  0 -1 -1  0  0  0
     0  1  1  2  1 -1 -1  0  0     0  0  1  2  1  0 -1 -1  0     0 -1  0  1  1 -1 -2 -1  0
     0  0  2  4  4  1  0  1  0     0  0  2  5  5  1  0  0  0    -1  0  2  6  6  2  1  0  0
     0  1  3  5  7  3  1  1  0     0  1  4  8 11  5  2  1  0     0  1  6 13 17  9  4  1 -2
     0  1  1  3  2  4  3  2  0     0  2  4  6  8  9  4  1  0     0  2  7 14 20 15  6  1 -1
     0 -1  0  4  4  4  2  0  0     0  0  2  6  8  7  3  0  0     0  1  4  9 14 11  4  0 -1
     0 -1  0  2  3  2  1  0  0     0 -1  0  3  5  4  1  0  0    -1 -1  0  3  6  6  2  0  0
     0 -1  0  1  2  1  0 -1  0     0 -2 -1  0  2  2  0 -1  0     0 -2 -2  0  2  2  1 -1  0
     0  0  0  0  0  0  0  0  0     0  0  0 -1  0  1  0  0  0     0  0 -2 -1  0  1  0  0  0


     0  0 -1  1 -1 -2 -2  0  0     0  0 -3 -1 -1 -2 -3  0  0     0  0  0  0  0 -1 -2  0  0
     0  0  1  2  1 -1 -1 -1  0     0 -1  2  4  3  0 -1 -2  0     0 -1  1  2  1  0 -2 -2  0
    -1  0  3  8  8  4  3  0  0    -1  1  5  9  9  4  2  0  0     0  0  4  7  7  4  1 -1  0
     0  1  7 18 27 13  6  2 -2     0  2  6 20 38 16  5  1 -1    -1  1  5 16 26 14  5  0 -1
     0  2  8 26 59 25  6  2 -1     0  3  8 40100 37  7  2 -1     0  3  8 29 58 29  9  1 -1
     0  1  5 13 26 18  6  2 -1     0  2  7 17 34 23  7  2 -1     0  2  7 18 31 19  8  1 -1
    -1 -1  1  4  7  8  4  1  0    -2  0  2  3  6  8  5  0  0    -2  0  2  5  7  7  3  0  0
     0 -1 -1 -1  1  2  1 -2  0     0 -2 -1  0  2  3  2  0  0     0 -2 -1  0  1  1  1 -1  0
     0  0 -2 -1 -1 -1  0  0  0     0  0 -2 -2 -1 -1  0  0  0     0  0 -2 -1 -1  0  0  0  0


     0  0  0  1  0 -1  0  0  0     0  0  0  0  0 -2  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1  0  1  1  0 -2 -2  0     0 -1  0  1  2  0 -2 -2  0     0  0  0  1  1  1  0 -1  0
     0 -1  1  4  5  4  0 -1  0     0  0  1  3  5  3  0 -1  0     0  1  1  2  3  2  0  0  0
     0  1  5 11 15 11  5  0  0     0  1  3  6  8  7  3  1  0     0  1  3  3  4  3  0  0  0
     0  1  6 14 17 12  8  1  0     0  2  4  6  5  5  3  1  0     0  2  2  2  2  2  1  1  0
     1  1  5 12 17 12  6  1 -1     0  1  1  6  8  6  3  0  0     0  1 -1  2  5  4  2  0  0
     0  0  0  3  8  6  1 -1  0     0  0 -1  2  7  5  2  0  0     0  0 -1  1  5  4  1  0  0
     0 -1 -2 -1  2  2  0 -1  0     0 -1 -1  0  3  2  1 -1  0     0 -1 -1  0  4  3  1  0  0
     0  0 -2 -1  0  0  0  0  0     0  0 -1 -1  1  1  0  0  0     0  0  0  0  0  0  0  0  0

  30664 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      466 ...     491
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   466   0  0.942  5479531    0  14785     308    37  0.0389  0.3805
   467   0  0.963  5481084    0  14736     287    33  0.0390  0.3909
   468   0  0.948  5479969    0  14729     288    21  0.0410  0.5003
   469   0  0.952  5481107    0  14784     278    21  0.0392  0.4808
   470   0  0.957  5481572    0  14746     279    27  0.0375  0.4767
   471   0  0.939  5481003    0  14715     295    28  0.0402  0.4590
   472   0  0.953  5480660    0  14743     311    37  0.0393  0.4349
   473   0  0.952  5480010    0  14745     310    26  0.0371  0.4511
   474   0  0.958  5479830    0  14755     293    29  0.0391  0.4577
   475   0  0.964  5479488    0  14784     294    24  0.0386  0.5244
   476   0  0.954  5479668    0  14769     299    19  0.0363  0.5270
   477   0  0.971  5480002    0  14757     301    22  0.0381  0.5151
   478   0  0.971  5480642    0  14782     284    19  0.0380  0.5288
   479   0  0.949  5481117    0  14755     281    18  0.0371  0.5103
   480   0  0.948  5480055    0  14769     280    16  0.0378  0.4801
   481   0  0.950  5480278    0  14792     282    22  0.0387  0.4422
   482   0  0.959  5479479    0  14774     306    28  0.0385  0.3778
   483   0  0.948  5481137    0  14756     309    31  0.0376  0.4094
   484   0  0.954  5481276    0  14745     302    28  0.0371  0.3790
   485   0  0.959  5481108    0  14776     292    22  0.0379  0.4078
   486   0  0.953  5480738    0  14749     288    21  0.0382  0.4651
   487   0  0.960  5480577    0  14754     288    22  0.0379  0.4610
   488   0  0.957  5479238    0  14767     300    12  0.0406  0.5335
   489   0  0.962  5478424    0  14821     309    25  0.0391  0.4317
   490   0  0.981  5478893    0  14790     308    19  0.0377  0.4988
   491   0  0.967  5479066    0  14760     315    32  0.0364  0.4419

  1237 OUT OF   1373 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.27
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.45
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     61.992   187.669    72.575  90.000 100.055  90.000
 REC. CELL PARAMETERS   0.016383  0.005329  0.013994  90.000  79.945  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.084    18.406    18.941
 COORDINATES OF UNIT CELL B-AXIS   -14.838  -110.295   151.111
 COORDINATES OF UNIT CELL C-AXIS    18.452   -57.545   -40.190
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.075     0.049    -0.091
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1237 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  402     95     76     86    141    128     83     99    119
 SIGMAB (degree)       0.065  0.070  0.069  0.066  0.065  0.066  0.066  0.067  0.068
 SIGMAR (degree)       0.448  0.401  0.436  0.445  0.444  0.450  0.444  0.426  0.396


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1237 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  1  2 -2  0  0  0     0  0  1  1  0  0 -1  0  0
     0  0  1  2  2  0  0  0  0     0 -1  1  2  2 -1 -1 -1  0     0  0  1  2  1 -1 -1 -1  0
     0  1  2  3  3  0  0  0  0     0  0  2  4  4  1  0  0  0    -1  0  2  5  4  2  1  0  0
     0  2  4  6  6  1  0  0  0     0  1  4  8  8  4  2  0  0     0  1  5 12 15  9  4  0  0
     0  1  2  4  3  4  2  1  0     0  1  4  7  7  7  4  0  0     0  2  6 13 18 13  5  1  0
     0  0 -1  3  4  3  2  1  0     0  0  1  5  6  5  3  1  0     0  0  3  8 12 10  4  1  0
     0 -1 -1  1  3  2  2  0  0     0 -1  0  2  4  3  1  1  0    -1 -1  0  3  6  5  2  0  0
     0 -1 -1  0  2  1  0  0  0     0 -2 -1  0  2  2  0 -1  0     0 -1 -1  0  2  2  0 -1  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0     0  0 -3 -1  0  0  0  0  0


     0  0  0  1  0 -1 -2  0  0     0  0 -1  1  0 -2 -3  0  0     0  0 -1  0 -1 -1 -2  0  0
     0  0  1  2  2 -1 -1 -1  0     0  0  3  5  3  0 -1 -1  0     0 -1  1  2  1 -1 -2 -2  0
    -1 -1  3  8  8  3  2  0  0    -1  1  4  8  8  3  2  0  0    -1  1  3  7  8  4  1 -1  0
     0  1  6 18 27 12  5  1 -1    -1  2  5 19 38 16  4  1 -1     0  2  5 15 27 14  4  1  0
     0  2  7 25 55 23  6  2 -1     0  2  8 38100 38  6  1 -1     0  3  7 27 58 29  8  1  0
     0  1  5 12 24 17  5  1 -1    -1  2  6 16 34 23  6  2 -1    -1  2  6 16 29 19  7  1 -1
    -1  0  1  3  7  8  4  1  0    -1  0  2  3  6  8  4  0  0    -1  0  2  4  7  7  3 -1  0
     0 -1 -1  0  2  3  2 -1  0     0 -1 -1  0  2  4  2  0  0     0 -2 -1  0  1  2  1 -1  0
     0  0 -2 -2  0  1  0  0  0     0  0 -2 -1  0  1  0  0  0     0  0 -1 -1  0  0  0  0  0


     0  0 -2  0  0 -1  0  0  0     0  0  0  0  1  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1  0  0  1  0 -1 -2  0     0  0  0  0  1  0 -1 -1  0     0  0  1  1  1  1  0  0  0
     0  0  1  3  5  3  0 -1  0     0  0  1  3  4  2  0 -1  0     0  0  1  3  3  1 -1 -1  0
     1  1  4 10 15 10  3  1  0     0  2  3  6  8  5  1  0  0     0  1  3  4  5  3  0  0  0
     0  1  5 14 18 13  6  1  1     0  1  2  6  6  5  3  1  0     0  1  1  3  3  3  1  0  0
     0  1  4 11 16 12  5  1 -1     0  1  1  4  8  7  3  1  0     0  0  0  1  4  4  3  0  0
     0  0  0  2  6  5  1 -1  0     0  0  0  1  6  5  2  0  0     0  0 -1  1  4  3  2  1  0
     0 -1 -1 -1  1  2  0 -1  0     0 -1 -1 -1  2  2  1 -1  0     0  0 -1  0  3  2  1  0  0
     0  0 -2 -1  0  1  0  0  0     0  0 -2 -1  2  0  0  0  0     0  0  0  0  0  0  0  0  0

  30770 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      492 ...     516
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   492   0  0.982  5479566    0  14745     311    27  0.0374  0.4600
   493   0  0.967  5478303    0  14776     332    26  0.0380  0.4850
   494   0  0.964  5479426    0  14760     321    22  0.0380  0.4806
   495   0  0.964  5479808    0  14810     293    28  0.0392  0.4448
   496   0  0.966  5479506    0  14790     303    21  0.0392  0.4451
   497   0  0.967  5478913    0  14797     317    19  0.0383  0.3974
   498   0  0.970  5478984    0  14782     310    16  0.0390  0.4155
   499   0  0.969  5479455    0  14796     304    20  0.0396  0.4755
   500   0  0.968  5479193    0  14835     312    16  0.0380  0.4652
   501   0  0.964  5479712    0  14810     302    18  0.0398  0.4561
   502   0  0.979  5479106    0  14768     320    23  0.0389  0.4658
   503   0  0.974  5478705    0  14832     319    21  0.0373  0.4921
   504   0  0.967  5478220    0  14826     335    25  0.0370  0.4847
   505   0  0.960  5476695    0  14841     357    31  0.0381  0.4170
   506   0  0.984  5478058    0  14848     330    35  0.0377  0.3921
   507   0  0.977  5478785    0  14803     328    29  0.0376  0.4792
   508   0  0.972  5479429    0  14788     318    29  0.0371  0.4465
   509   0  0.986  5480319    0  14825     306    26  0.0374  0.4484
   510   0  0.985  5480383    0  14845     308    27  0.0359  0.3977
   511   0  0.985  5480991    0  14807     310    21  0.0366  0.4813
   512   0  0.979  5481371    0  14841     309    24  0.0368  0.4697
   513   0  0.970  5480771    0  14864     292    25  0.0384  0.4658
   514   0  0.972  5478927    0  14837     293    24  0.0387  0.4804
   515   0  0.977  5478289    0  14872     305    13  0.0386  0.5086
   516   0  0.973  5477962    0  14878     318    17  0.0381  0.5056

  1235 OUT OF   1382 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.25
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.47
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.011   187.745    72.570  90.000 100.099  90.000
 REC. CELL PARAMETERS   0.016380  0.005326  0.013997  90.000  79.901  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.082    18.473    18.945
 COORDINATES OF UNIT CELL B-AXIS   -14.755  -110.290   151.217
 COORDINATES OF UNIT CELL C-AXIS    18.457   -57.550   -40.173
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.084     0.023    -0.068
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1235 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  392     92     71     84    141    128     86    104    123
 SIGMAB (degree)       0.063  0.068  0.067  0.065  0.064  0.065  0.066  0.066  0.068
 SIGMAR (degree)       0.452  0.416  0.449  0.455  0.428  0.438  0.451  0.433  0.411


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1235 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0     0  0  0  0 -1 -1 -1  0  0
     0  0  1  1  1  0 -1 -1  0     0 -1  1  2  1 -1 -1  0  0     0 -1  1  1  1 -1 -1 -1  0
     0  1  1  3  2  0 -1  0  0     0  0  2  4  3  0  0  0  0    -1  0  2  4  4  1  0  0  0
     0  0  0  0  4  2  1  1  0     0  1  3  6  8  3  1  0  0     0  1  4 10 13  7  3  0 -1
     0  0  0  0  4  3  1  1  0     0  1  3  6  8  7  2  0  0    -1  1  5 13 20 13  5  1 -1
     0  0  0  0  5  3  1  0  0     0  0  0  5  7  6  2  0  0     0  0  3  9 13 10  4  1 -1
     0  0  0  1  3  1  0 -1  0     0  0  0  1  4  2  1 -1  0     0 -1  0  1  4  4  1 -1  0
     0  0  0  0  1  1  0 -1  0     0 -1 -1  0  1  1  0 -1  0     0 -2 -2 -1  1  2  0 -1  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -2 -2  0  0  0  0  0


     0  0  0  0 -1 -2 -2  0  0     0  0  0  1  0 -1 -2  0  0     0  0 -1 -1  0 -2 -3  0  0
     0  0  0  2  1 -1 -1 -1  0     0  0  2  4  2 -1 -1 -2  0     0  0  1  1  1 -1 -2 -2  0
     0  0  3  7  7  2  1  0 -1    -1  1  4  7  7  2  1  0 -2    -1  1  3  6  7  4  0 -1  0
    -1  1  5 17 24 11  4  1 -1     0  1  5 20 37 14  3  1 -1     0  1  5 16 27 14  4  0  0
    -1  1  6 25 53 21  6  2  0     0  2  7 37100 36  6  2  0     0  2  7 28 60 28  7  1  0
    -1  1  5 13 23 16  5  2  0     0  1  5 13 33 23  6  2  0     0  1  5 13 26 18  6  1 -1
    -1 -1  1  3  6  6  3  0  0    -1  0  1  2  5  6  3  1  0    -2  0  1  2  5  6  3  0  0
     0 -2 -2 -1  1  2  1  0  0     0 -1 -1  0  2  3  2 -1  0     0 -2 -1  0  1  1  1 -1  0
     0  0 -3 -1 -1  0  0  0  0     0  0 -2 -1  0  1 -1  0  0     0  0 -2 -1  0  1  0  0  0


     0  0 -1  0 -1 -2  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1  0  1  0  0 -2 -2  0     0  0  0  1  1  0 -1 -1  0     0  1  0  0  0  0  1  0  0
    -1  0  2  3  4  3 -1 -1  0     0  0  1  2  3  2  0 -1  0     0  0  0  1  1  0 -1  0  0
     0  1  4 10 16 11  3  0  0     0  0  2  5  7  6  1  0  0     0  0  1  3  5  3  0  0  0
     0  1  5 14 20 14  6  1  1     0  0  2  5  6  7  3  2  0     0  0  1  3  4  4  0  0  0
    -1  1  3 10 15 11  4  1  0     0  0  1  3  7  7  3  1  0     0  0  0  1  4  5  0  1  0
    -1  0  0  1  4  5  1 -1  0     0  0  0  0  3  4  2  0  0     0  0  0  0  2  2  1  0  0
     0 -1 -1 -1  1  1  0 -1  0     0  0 -1 -1  1  1  0 -1  0     0  0 -1  0  1  1  1  0  0
     0  0 -2 -1  0  0  0  0  0     0  0  0  0  0  1  0  0  0     0  0  0  0  0  0  0  0  0

  30181 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      517 ...     542
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   517   0  0.977  5478240    0  14809     313     8  0.0387  0.5571
   518   0  0.984  5479239    0  14774     310    24  0.0385  0.4141
   519   0  1.000  5479129    0  14796     320    27  0.0358  0.4091
   520   0  0.985  5478099    0  14804     326    22  0.0374  0.4831
   521   0  0.975  5476980    0  14797     319    20  0.0373  0.4686
   522   0  0.981  5478405    0  14747     320    18  0.0401  0.4743
   523   0  0.989  5478270    0  14805     317    21  0.0387  0.4974
   524   0  0.990  5477840    0  14780     337    21  0.0383  0.4767
   525   0  0.984  5478272    0  14717     341    24  0.0377  0.4725
   526   0  0.987  5478053    0  14727     325    20  0.0376  0.4947
   527   0  0.980  5478808    0  14733     307    20  0.0372  0.4673
   528   0  0.984  5477880    0  14768     338    11  0.0374  0.5095
   529   0  0.988  5478674    0  14731     328    13  0.0381  0.5370
   530   0  0.981  5479443    0  14759     304    11  0.0371  0.5272
   531   0  0.991  5478565    0  14727     297    12  0.0357  0.5232
   532   0  1.012  5479419    0  14743     316    20  0.0352  0.4775
   533   0  1.002  5478917    0  14772     314    21  0.0373  0.4331
   534   0  1.001  5478702    0  14783     321    13  0.0379  0.4799
   535   0  0.989  5478504    0  14667     316    19  0.0367  0.5080
   536   0  0.997  5478314    0  14737     323    15  0.0352  0.5280
   537   0  1.008  5478575    0  14684     347    24  0.0349  0.5326
   538   0  0.981  5477352    0  14661     346    23  0.0372  0.5083
   539   0  1.005  5477783    0  14725     343    25  0.0357  0.4695
   540   0  0.999  5476888    0  14717     364    28  0.0358  0.4400
   541   0  1.004  5477384    0  14683     358    33  0.0365  0.4182
   542   0  0.991  5476788    0  14690     369    31  0.0350  0.4653

  1373 OUT OF   1540 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.30
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.49
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.036   187.294    72.441  90.000 100.169  90.000
 REC. CELL PARAMETERS   0.016377  0.005339  0.014025  90.000  79.831  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.084    18.573    18.923
 COORDINATES OF UNIT CELL B-AXIS   -14.405  -110.143   150.798
 COORDINATES OF UNIT CELL C-AXIS    18.416   -57.403   -40.168
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.031    -0.005    -0.017
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1373 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  434    108     84    103    155    129     88    115    146
 SIGMAB (degree)       0.063  0.066  0.065  0.066  0.064  0.064  0.066  0.067  0.067
 SIGMAR (degree)       0.463  0.434  0.468  0.468  0.436  0.448  0.460  0.434  0.409


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1373 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1 -1  0  0
     0  0  0  1  0 -1 -1  0  0     0  0  1  1  1  0 -1 -1  0     0  0  0  1  1 -1 -1 -1  0
     0  0  0  0  0  0 -1  0  0     0  1  3  3  2  0  0 -1  0     0  0  2  3  3  1  0 -1  0
     0  0  0  0  0  2  0  0  0     0  1  4  7  7  4  1  0  0     0  0  4 11 14  8  3  0 -1
     0  0  0  0  0  6  2  0  0     0  1  4  8 10  8  3  0  0     1  1  6 16 22 14  6  1 -1
     0  0  0  0  0  4  2  1  0     0  0  1  4  7  6  2  0  0     0  0  3 10 12 10  5  1  0
     0  0  0  0  2  2  0  0  0     0 -1 -1  1  3  3  1  0  0    -1 -1 -1  1  3  4  1  0  0
     0  0  0  0  1  1  0  0  0     0 -1 -1 -1  2  1  0  0  0     0 -2 -1 -1  1  1  0 -1  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -2 -2  0  1  0  0  0


     0  0  0  0  0 -1 -3  0  0     0  0  0  0  0 -1 -3  0  0     0  0  0  0 -1 -2 -3  0  0
     0  0  1  2  1 -1 -2 -2  0     0  0  2  2  1 -1 -1 -2  0     0  0  1  1  1 -1 -2 -2  0
     0  1  3  6  5  2  1  0  0    -1  1  4  7  6  1  1  0 -3    -1  1  4  6  6  2  0 -2  0
     0  1  6 17 22  9  3  0 -1    -1  2  7 25 34 10  2  0 -1     0  2  8 21 27 11  3  0 -1
     0  1  8 30 52 21  5  0 -1    -1  2  9 40100 32  5  1 -1     0  2  9 33 60 24  6  1 -1
    -1  1  5 14 23 15  6  2  0    -1  1  5 12 35 23  6  2  0     0  1  5 12 25 18  6  1  0
    -1 -1  0  2  6  6  4  1  0    -2  0  0  2  6  7  4  2  0    -1 -1  0  2  5  7  3  1  0
     0 -2 -2 -1  1  1  2  1  0     0 -2 -1  0  1  2  2  1  0     0 -2 -1 -1  1  2  1  1  0
     0  0 -2 -2  0  0  0  0  0     0  0 -1 -1  0  1 -1  0  0     0  0 -2 -1 -1  0  0  0  0


     0  0  0 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1  0  0  0 -1 -2 -2  0     0 -1 -1  0  0  0 -1 -1  0     0 -1  0  0  0  0  0  0  0
    -2  0  1  4  4  2 -1 -1  0     0  0  1  2  2  1  0 -1  0     0  0  0  1  1  0  0  0  0
    -1  1  5 11 13  8  2  0 -1     0  1  3  5  7  4  0  0  0     0  1  2  4  0  0  0  0  0
     0  1  6 19 22 13  5  1  0     0  1  2  8  9  6  3  1  0     0  0  2  5  0  0  0  0  0
    -1  0  4 11 16 10  4  1  0     0  0  1  4  8  6  4  1  0     0  0  0  1  0  0  0  0  0
    -1 -1  0  1  4  5  2  0  0     0  0  0  0  3  4  3  0  0     0  0 -1 -1  1  0  0  0  0
     0 -1 -2 -1  1  2  1  0  0     0 -1 -1 -1  0  2  1  0  0     0  0 -1 -1  1  2  1  0  0
     0  0 -2 -2  0  1  0  0  0     0  0  0  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0

  30702 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      543 ...     568
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   543   0  1.000  5477892    0  14766     359    23  0.0383  0.4583
   544   0  1.007  5478432    0  14740     346    25  0.0366  0.4651
   545   0  0.994  5478209    0  14723     330    22  0.0379  0.4643
   546   0  1.015  5479078    0  14711     321    20  0.0374  0.4862
   547   0  0.992  5478391    0  14709     317    21  0.0364  0.4582
   548   0  1.006  5479043    0  14732     310    19  0.0365  0.4252
   549   0  1.012  5478921    0  14743     324    22  0.0377  0.3791
   550   0  1.008  5478606    0  14681     353    19  0.0364  0.4787
   551   0  1.011  5477600    0  14711     357    25  0.0362  0.4052
   552   0  0.998  5477352    0  14728     348    25  0.0374  0.4669
   553   0  0.999  5477232    0  14734     352    22  0.0347  0.4925
   554   0  1.018  5477745    0  14701     359    24  0.0340  0.4713
   555   0  1.004  5477024    0  14723     349    22  0.0374  0.4909
   556   0  1.016  5476626    0  14727     363    23  0.0367  0.4794
   557   0  1.005  5477004    0  14764     352    20  0.0382  0.4237
   558   0  1.020  5476224    0  14794     358    19  0.0375  0.4789
   559   0  1.026  5476583    0  14782     345    23  0.0365  0.4743
   560   0  1.008  5476829    0  14816     340    22  0.0373  0.4081
   561   0  1.016  5477255    0  14791     325    20  0.0368  0.4918
   562   0  1.027  5478125    0  14758     336    17  0.0360  0.5173
   563   0  1.030  5477804    0  14781     340    22  0.0351  0.5437
   564   0  1.015  5476915    0  14734     363    23  0.0356  0.5104
   565   0  1.006  5475986    0  14732     346    24  0.0375  0.4741
   566   0  1.025  5476376    0  14781     342    31  0.0372  0.3775
   567   0  1.027  5476367    0  14762     353    30  0.0376  0.4292
   568   0  1.009  5475443    0  14734     350    28  0.0375  0.4509

  1441 OUT OF   1680 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.27
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.47
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.038   187.737    72.441  90.000 100.163  90.000
 REC. CELL PARAMETERS   0.016376  0.005327  0.014024  90.000  79.837  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.083    18.575    18.931
 COORDINATES OF UNIT CELL B-AXIS   -14.465  -110.384   151.167
 COORDINATES OF UNIT CELL C-AXIS    18.429   -57.407   -40.156
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.064    -0.013    -0.014
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1441 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  449    119     93    105    150    129    102    128    158
 SIGMAB (degree)       0.062  0.064  0.064  0.064  0.063  0.063  0.066  0.066  0.065
 SIGMAR (degree)       0.465  0.444  0.476  0.474  0.407  0.441  0.465  0.443  0.412


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1441 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -1  0  0 -1  0  0  0
     0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0 -1  0  0     0 -1  0  0  0 -1 -1 -1  0
     0  0  0  0  0  0 -1 -1  0     0  0  2  2  1  0  0  0  0     0  1  2  3  2  1  0  0  0
     0  0  0  0  0  0  1  1  0     0  1  3  6  6  4  2  1  0     0  1  4 10 12  8  3  1  0
     0  0  0  0  0  0  3  1  0     0  1  4  9 12  9  3  0  0     0  1  5 17 24 14  4  1  0
     0  0  0  0  0  0  1  0  0     0  0  1  4  7  6  2  0  0     0  0  2 10 13  9  4  0  0
     0  0  0  0  0  2  0 -1  0     0 -1  0  1  2  3  1 -1  0     0 -1  0  1  4  4  2  0  0
     0  0  0  0  0  1  1  0  0     0 -1 -1 -1  1  1  1 -1  0     0 -2 -2  0  0  1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -2 -1  0 -1  0  0  0


     0  0  0  0  0 -1 -2  0  0     0  0  0  0  0 -2 -2  0  0     0  0 -1 -1 -1 -2 -3  0  0
     0  0  0  1  1 -1 -1 -1  0     0  0  1  2  1 -1 -1 -1  0     0  0  1  1  0 -1 -2 -2  0
    -1  1  3  5  5  1  0  0  0     0  1  3  6  6  2  1  0  0     0  1  4  5  5  3  0 -1  0
    -1  1  5 16 22 10  3  1  0     0  1  6 22 37 12  3  1  0     0  1  7 20 28 13  4  1  0
     0  1  8 30 53 21  5  1  0     0  2  8 38100 34  5  1  0     0  2  8 32 64 26  6  2  0
     0  1  6 14 22 14  5  1 -1    -1  1  5 12 31 22  5  1  0     0  1  5 11 23 17  6  2  0
    -1 -1  1  2  5  6  3  0  0    -1  0  1  1  5  6  3  1  0    -1  0  0  1  4  6  3  0  0
     0 -2 -2 -1  1  2  1  0  0     0 -1 -1 -1  1  2  1  1  0     0 -1 -1 -1  0  2  1  0  0
     0  0 -2 -2  0  0 -1  0  0     0  0 -1 -1  0  0 -2  0  0     0  0 -2 -1 -1 -1  0  0  0


     0  0  0 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1  0  0  0 -1 -2 -2  0     0 -1  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0
    -1  0  2  4  4  1 -1 -1  0     0  0  1  2  2  0 -1 -1  0     0  0  0  2  0  0  0  0  0
    -1  0  4 11 14  9  2  0  0     0  0  2  5  6  3  0  1  0     0  0  1  3  0  0  0  0  0
     0  1  6 19 26 15  5  1  0     0  1  3  8 11  8  4  1  0     0  1  2  0  0  0  0  0  0
     0  1  4 10 15 10  3  1  0     0  1  2  4  6  6  3  1  0     0  0  0  0  0  0  0  0  0
    -1  0  0  0  2  4  2  0  0     0  0  0  0  1  3  3  0  0     0 -1 -1  0  0  0  0  0  0
     0 -1 -1 -1  0  1  1  0  0     0 -1 -1 -1  0  1  1  1  0     0 -1 -1 -1  0  0  0  0  0
     0  0 -2 -1 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30717 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      569 ...     594
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   569   0  1.021  5475408    0  14758     361    28  0.0363  0.4490
   570   0  1.010  5475550    0  14770     352    22  0.0348  0.4849
   571   0  1.032  5475777    0  14805     338    13  0.0368  0.5255
   572   0  1.016  5476007    0  14813     340    14  0.0378  0.4814
   573   0  1.014  5476361    0  14772     358    15  0.0354  0.4404
   574   0  1.018  5475466    0  14777     385    23  0.0363  0.4535
   575   0  1.029  5476391    0  14773     364    24  0.0368  0.4045
   576   0  1.030  5477088    0  14817     356    21  0.0364  0.4779
   577   0  1.014  5476137    0  14814     371    12  0.0351  0.5473
   578   0  1.019  5476490    0  14769     358    12  0.0370  0.5295
   579   0  1.032  5476997    0  14798     365    17  0.0340  0.5077
   580   0  1.024  5476800    0  14763     354    13  0.0358  0.5249
   581   0  1.029  5477486    0  14766     348    15  0.0378  0.5173
   582   0  1.027  5476957    0  14774     339    14  0.0353  0.4992
   583   0  1.030  5476895    0  14772     346    11  0.0361  0.5138
   584   0  1.037  5475772    0  14755     365    16  0.0369  0.4768
   585   0  1.020  5475013    0  14753     357    16  0.0374  0.4707
   586   0  1.034  5475470    0  14720     358    18  0.0387  0.4420
   587   0  1.035  5475401    0  14737     341    13  0.0405  0.4750
   588   0  1.033  5475709    0  14748     351    14  0.0392  0.4737
   589   0  1.028  5474937    0  14799     371    16  0.0382  0.4640
   590   0  1.028  5474876    0  14734     364    15  0.0363  0.4745
   591   0  1.046  5474729    0  14789     368    16  0.0358  0.5056
   592   0  1.044  5475061    0  14774     385    23  0.0363  0.5014
   593   0  1.041  5475041    0  14756     377    22  0.0343  0.4998
   594   0  1.030  5474920    0  14779     375    19  0.0370  0.5157

  1382 OUT OF   1711 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.25
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.49
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.039   188.062    72.423  90.000 100.149  90.000
 REC. CELL PARAMETERS   0.016375  0.005317  0.014027  90.000  79.851  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.088    18.554    18.939
 COORDINATES OF UNIT CELL B-AXIS   -14.561  -110.555   151.436
 COORDINATES OF UNIT CELL C-AXIS    18.425   -57.401   -40.134
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.096    -0.001    -0.024
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1382 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  436    114     87     98    127    112    102    137    163
 SIGMAB (degree)       0.061  0.062  0.061  0.062  0.061  0.062  0.064  0.065  0.064
 SIGMAR (degree)       0.480  0.462  0.496  0.492  0.436  0.449  0.475  0.457  0.417


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1382 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -1  0 -1 -1 -1  0  0
     0  0  0  0  0  0  0  0  0     0 -1  0  0  0 -1 -1 -1  0     0  0  0  0  0 -1 -1 -1  0
     0  0  0  0  0  0  0  0  0     0  0  1  1  1  0  0  0  0     0  1  1  2  1  1  0 -1  0
     0  0  0  0  0  0  0  1  0     0  1  3  6  6  4  2  1  0     0  1  4 10 11  7  3  0  0
     0  0  0  0  0  0  0  1  0     0  2  4 10 14  9  3  1  0     1  2  6 21 30 17  5  1  0
     0  0  0  0  0  0  3  1  0     0  0  1  3  7  7  3  1  0     0  1  2  9 14 11  4  1  1
     0  0  0  0  0  0  1  0  0     0 -1  0  0  2  3  1 -1  0     0  0 -1  1  3  4  1  0  0
     0  0  0  0  0  0  1  0  0     0 -1 -1  0  0  1  0  0  0     0 -1 -1  0  0  1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -2 -1 -1 -1  0  0  0


     0  0  0 -1 -1 -1 -1  0  0     0  0  0 -1 -1 -2 -2  0  0     0  0 -1  0 -1 -2 -2  0  0
     0  0  0  0  1 -1 -1 -1  0     0  0  1  1  1 -1 -1 -1  0     0  0  0  0  1 -1 -2 -2  0
    -1  1  2  4  3  0  0  0 -1     0  1  3  5  5  1  0 -1 -1     0  1  3  4  5  2  0 -1  0
     0  2  6 16 19  8  2  0 -1     0  2  7 21 34 11  3  0  0     1  2  7 19 28 12  3  0 -1
     0  2  9 33 59 22  5  1  0     0  2  8 37100 34  6  1 -1     1  2  7 30 67 29  8  2  0
     0  1  5 12 21 14  5  2  1    -1  1  4 10 29 20  5  2  1     0  1  3  9 21 17  7  2  0
    -1 -1  0  1  4  5  2  1  0    -2  0  0  1  4  6  3  1  0    -1 -1  0  1  3  5  3  1  0
     0 -2 -2 -1  1  2  1  1  0     0 -2 -1 -1  1  2  2  1  0     0 -1 -1 -1  1  1  1  0  0
     0  0 -2 -1  0  0 -2  0  0     0  0 -2 -1 -1  0  0  0  0     0  0 -2 -1 -1 -1 -1  0  0


     0  0 -1 -1 -1 -2  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0 -1 -1 -2  0     0 -1  1  0  0 -1 -1 -1  0     0 -1  0  0  0  0  0  0  0
    -1  0  2  3  3  1  0 -1  0     0  0  2  2  1  0 -1  0  0     0  1  1  0  0  0  0  0  0
     0  1  5 11 13  8  2  0  0     0  0  3  7  7  3  1  0  0     0  0  2  0  0  0  0  0  0
     0  1  6 22 32 18  6  1  0     0  0  3 10 14 10  4  1  0     0  0  2  0  0  0  0  0  0
     0  1  3  9 14 10  5  1  0    -1  0  1  4  6  8  4  1  0     0  0  0  0  0  0  0  0  0
    -1 -1  0  0  2  4  2  1  0     0  0  0  0  1  3  2  1  0     0  0  0  0  0  0  0  0  0
     0 -1 -1 -1  0  1  0 -1  0     0 -1 -1  0  1  1  1  0  0     0 -1  0  0  0  0  0  0  0
     0  0 -1 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30758 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      595 ...     619
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   595   0  1.041  5475751    0  14790     376    27  0.0372  0.4745
   596   0  1.055  5474988    0  14781     366    30  0.0366  0.4516
   597   0  1.042  5475068    0  14826     356    23  0.0390  0.4423
   598   0  1.035  5475067    0  14810     343    21  0.0380  0.4620
   599   0  1.042  5474627    0  14811     358    21  0.0365  0.4486
   600   0  1.043  5475342    0  14763     353    24  0.0380  0.4619
   601   0  1.045  5475152    0  14761     380    19  0.0372  0.4680
   602   0  1.042  5475654    0  14779     377    22  0.0366  0.4332
   603   0  1.046  5474672    0  14794     397    19  0.0373  0.4845
   604   0  1.044  5474852    0  14806     375    17  0.0355  0.4854
   605   0  1.047  5475491    0  14809     384    10  0.0354  0.5227
   606   0  1.065  5474488    0  14817     373    14  0.0375  0.5211
   607   0  1.049  5474036    0  14817     389    15  0.0359  0.5039
   608   0  1.053  5474879    0  14801     375    14  0.0387  0.4929
   609   0  1.046  5474992    0  14834     394    17  0.0386  0.4816
   610   0  1.057  5475581    0  14842     373    14  0.0366  0.4801
   611   0  1.041  5475266    0  14793     370    17  0.0369  0.4768
   612   0  1.054  5475708    0  14795     368    23  0.0346  0.4473
   613   0  1.058  5475554    0  14803     395    20  0.0360  0.4373
   614   0  1.052  5475967    0  14796     376    18  0.0363  0.4939
   615   0  1.041  5476053    0  14782     374    25  0.0359  0.4974
   616   0  1.045  5474871    0  14854     385    17  0.0365  0.4964
   617   0  1.051  5474955    0  14827     366    17  0.0372  0.4968
   618   0  1.057  5474557    0  14787     392    26  0.0354  0.4787
   619   0  1.046  5473838    0  14761     394    17  0.0359  0.5268

  1220 OUT OF   1701 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.28
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.51
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.044   188.431    72.403  90.000 100.155  90.000
 REC. CELL PARAMETERS   0.016374  0.005307  0.014031  90.000  79.845  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.093    18.555    18.942
 COORDINATES OF UNIT CELL B-AXIS   -14.578  -110.806   151.710
 COORDINATES OF UNIT CELL C-AXIS    18.413   -57.377   -40.138
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.114     0.007    -0.019
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1220 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  384    110     71     71     87     89     92    137    168
 SIGMAB (degree)       0.059  0.061  0.060  0.060  0.059  0.060  0.063  0.064  0.063
 SIGMAR (degree)       0.488  0.443  0.508  0.504  0.476  0.452  0.482  0.446  0.421


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1220 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0 -1 -1  0  0  0
     0  0  0  0  0  0  0 -1  0     0  0  0  0  0  0 -1  0  0     0  0  0  0  0 -1 -1 -1  0
     0  0  0  0  0  0  0  0  0     0  0  1  1  0  0  0  0  0     0  0  1  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  1  2  5  5  3  1  0 -1     0  1  5 10 11  6  1  0 -1
     0  0  0  0  0  0  0  1  0     0  1  3 12 15  8  2  0 -1     0  1  6 25 32 13  3  0  0
     0  0  0  0  0  0  2  0  0     0  0  1  4  6  6  2  1 -1     0  0  2  9 14 10  3  1  0
     0  0  0  0  0  0  2  1  0     0  0  0  0  1  3  3  1  0     0  0 -1  0  3  4  3  1  0
     0  0  0  0  0  0  1  1  0     0 -1 -1 -1  0  1  1  1  0     0 -1 -2 -1  0  1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -1 -1 -1 -1  0  0  0


     0  0 -1 -1  0 -1 -2  0  0     0  0 -1  0  0 -1 -2  0  0     0  0 -1  0  0 -1 -2  0  0
     0  0  0  1  1 -1 -1 -1  0     0  0  0  1  1  0 -1 -1  0     0  1  1  1  1 -1 -2 -1  0
     0  1  2  4  3  1  0 -1  0     1  1  3  4  5  1  0 -1 -1     1  1  2  4  4  2  0 -1  0
     0  1  7 17 19  7  1  0 -1     0  1  6 20 33  9  2  0 -1     0  1  5 16 28 12  3  0 -1
     0  1  8 35 59 18  3  1  0     0  1  6 34100 31  4  1 -1     0  1  5 25 68 30  6  1  0
     0  1  4 11 20 13  4  1  0    -1  1  3  8 28 21  5  1 -1     0  1  3  7 19 18  7  1  0
    -2 -1  0  1  4  6  3  1  0    -1  0  0  1  4  6  4  1  0     0  0  0  1  3  5  4  1  0
     0 -2 -2 -1  1  2  1  1  0     0 -1 -1 -1  1  2  2  1  0     0 -1 -1 -1  1  1  1  0  0
     0  0 -2 -1  0  0 -1  0  0     0  0 -2 -1  0  0  1  0  0     0  0 -1 -1  0 -1  0  0  0


     0  0 -1  0 -1 -2  0  0  0     0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0 -1 -2 -1  0     0  0  1  0  0 -1 -1 -1  0     0 -1  0  0  0  0  0  0  0
     0  1  2  2  2  1  0  0  0     1  0  2  2  1  0  0  0  0     0  0  1  0  0  0  0  0  0
     0  1  3 10 13  8  3  1  0     0  1  2  6  7  3  1  0  0     0  1  0  0  0  0  0  0  0
     0  1  5 19 35 21  6  1  0    -1  0  3 10 17 11  4  1  0     0  1  0  0  0  0  0  0  0
     0  0  3  7 12 11  4  0  0    -1  0  1  4  6  6  3  1  0     0  0  0  0  0  0  0  0  0
    -1 -1  0  0  2  3  1  0  0     0 -1  0  0  1  2  1  1  0     0 -1  0  0  0  0  0  0  0
     0 -1 -1 -1  0  1  0  0  0     0 -1 -1 -1  0  1  0  0  0     0 -1  0  0  0  0  0  0  0
     0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30167 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      620 ...     644
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   620   0  1.051  5474418    0  14759     405    21  0.0347  0.5126
   621   0  1.061  5474755    0  14808     400    21  0.0364  0.4841
   622   0  1.048  5474280    0  14794     408    17  0.0352  0.5085
   623   0  1.063  5474412    0  14844     392    18  0.0376  0.4792
   624   0  1.050  5473286    0  14792     388    15  0.0388  0.5129
   625   0  1.048  5474592    0  14858     386    12  0.0372  0.5101
   626   0  1.073  5473726    0  14859     399    14  0.0365  0.4817
   627   0  1.054  5473678    0  14851     390    18  0.0382  0.4815
   628   0  1.055  5473355    0  14838     396    21  0.0344  0.4936
   629   0  1.058  5474016    0  14848     381    21  0.0386  0.4791
   630   0  1.069  5473862    0  14882     406    23  0.0362  0.4473
   631   0  1.065  5473864    0  14852     389    26  0.0359  0.4325
   632   0  1.059  5473449    0  14867     413    24  0.0357  0.4540
   633   0  1.050  5473886    0  14834     395    22  0.0365  0.4378
   634   0  1.040  5474047    0  14801     383    19  0.0374  0.4395
   635   0  1.058  5473804    0  14820     397    16  0.0379  0.4795
   636   0  1.072  5473409    0  14833     416    25  0.0344  0.4678
   637   0  1.051  5472946    0  14803     422    30  0.0369  0.4763
   638   0  1.068  5473707    0  14809     404    37  0.0347  0.4218
   639   0  1.067  5472865    0  14768     411    26  0.0367  0.4640
   640   0  1.065  5473416    0  14811     414    21  0.0356  0.5082
   641   0  1.057  5473467    0  14824     419    19  0.0354  0.5200
   642   0  1.052  5473032    0  14842     382    25  0.0351  0.4601
   643   0  1.062  5472445    0  14794     404    19  0.0355  0.4926
   644   0  1.074  5471505    0  14812     430    21  0.0343  0.4697

  1237 OUT OF   1925 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.28
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.51
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.049   188.329    72.428  90.000 100.170  90.000
 REC. CELL PARAMETERS   0.016374  0.005310  0.014027  90.000  79.830  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.090    18.620    18.902
 COORDINATES OF UNIT CELL B-AXIS   -14.402  -110.640   151.720
 COORDINATES OF UNIT CELL C-AXIS    18.432   -57.415   -40.119
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.113    -0.036     0.029
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1237 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  376    112     61     55     79    102    116    149    178
 SIGMAB (degree)       0.057  0.059  0.057  0.058  0.057  0.060  0.061  0.062  0.062
 SIGMAR (degree)       0.503  0.411  0.515  0.512  0.476  0.437  0.478  0.464  0.408


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1237 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0 -1 -1 -1 -2  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0 -1 -1  0     0  0  0  0 -1 -1 -1  0  0
     0  0  0  0  0  0  0  1  0     0  0  0  1  0  0  1  0  0     0  0  1  2  1  0  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  6  4  3  2  1  0     0  0  4 10 10  6  2  1 -1
     0  0  0  0  0  0  0  1  0     0  0  4 17 19  8  3  1  0     0  1  8 32 37 12  3  0  0
     0  0  0  0  0  0  0  0  0     0  0  1  3  6  6  2  1 -1     0  1  2  9 14  8  2  0  0
     0  0  0  0  0  0  0  1  0     0  0 -1  0  0  2  2  1  0     0  0 -1  0  2  3  2  0  0
     0  0  0  0  0  0  0 -1  0     0  0 -1  0  0  0  0  0  0     0 -1 -1 -1  0  1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -1 -1 -1 -1  0  0  0


     0  0 -1 -1  0 -1 -1  0  0     0  0 -1  0  0 -1 -1  0  0     0  0 -1  0  0 -1 -2  0  0
     0 -1  0  0  0 -1 -1  0  0     0 -1  0  1  1  0 -1 -1  0     0  0  0  1  1  0 -2 -1  0
     0  0  1  4  3  1  1  0  0     1  0  2  4  5  2  1  0 -1     0  0  1  3  5  2  0  0  0
     0  1  6 16 19  7  2  0  0     0  1  5 17 34 11  3  0  0     0  1  4 13 28 13  3  0  0
     0  2  9 35 59 16  3  1  0     0  1  6 30100 30  4  1  0     0  1  4 21 70 32  6  1  0
     0  1  4 11 19 10  3  0  0     0  1  3  8 28 17  4  0  0     0  1  2  6 18 16  5  1  0
    -1 -1  0  1  4  4  2  0  0    -1  0  0  2  5  4  2  0  0    -1  0  0  1  3  4  2  0  0
     0 -2 -1 -1  1  1  0  0  0     0 -1 -1  0  1  1  0  0  0     0 -1 -1  0  1  1  0  0  0
     0  0 -1 -1  0  0 -2  0  0     0  0 -1 -1  0  1 -1  0  0     0  0 -1 -1 -1 -1 -1  0  0


     0  0 -1 -1 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1  0  0  0  0 -2 -1  0     0 -1  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0
     0  0  1  2  3  1  0  0  0     0  1  1  1  1  0 -1  0  0     0  0  0  0  0  0  0  0  0
     0  1  4  9 14  9  2  0  0     0  0  2  6  8  3  0  0  0     0  1  0  0  0  0  0  0  0
     0  1  4 16 39 24  6  1  0    -1  1  3 11 22 12  3  0  0     0  1  0  0  0  0  0  0  0
    -1  1  2  5 10 11  4  1  0     0  0  2  3  5  6  3  0  0     0  0  0  0  0  0  0  0  0
    -1  0  0  0  1  2  2  0  0     0  0  0  0  0  1  2  0  0     0  0  0  0  0  0  0  0  0
     0 -1 -1 -1  0  0  1  0  0     0  0  0 -1  0  0  1  0  0     0  0  0  0  0  0  0  0  0
     0  0 -1 -1 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30176 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      645 ...     670
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   645   0  1.060  5472095    0  14766     424    21  0.0361  0.4969
   646   0  1.058  5471978    0  14779     426    28  0.0377  0.4793
   647   0  1.055  5472404    0  14756     413    23  0.0351  0.4969
   648   0  1.078  5471917    0  14711     401    20  0.0386  0.5035
   649   0  1.063  5471271    0  14735     398    21  0.0390  0.4991
   650   0  1.047  5471618    0  14735     411    12  0.0374  0.5565
   651   0  1.061  5472215    0  14704     396    10  0.0367  0.5796
   652   0  1.063  5471681    0  14753     410    18  0.0354  0.5207
   653   0  1.073  5473120    0  14757     410    21  0.0356  0.5121
   654   0  1.056  5473932    0  14749     412    21  0.0348  0.5119
   655   0  1.062  5473828    0  14767     388    25  0.0358  0.4643
   656   0  1.071  5473305    0  14781     405    25  0.0360  0.4667
   657   0  1.069  5473258    0  14778     403    23  0.0358  0.4722
   658   0  1.074  5472736    0  14721     400    22  0.0364  0.4592
   659   0  1.051  5470817    0  14665     423    14  0.0355  0.5053
   660   0  1.067  5471077    0  14712     423    20  0.0367  0.4921
   661   0  1.059  5471067    0  14757     432    14  0.0343  0.4825
   662   0  1.065  5470392    0  14764     442    16  0.0376  0.5046
   663   0  1.065  5469988    0  14746     420    15  0.0367  0.4849
   664   0  1.067  5471400    0  14732     411    20  0.0364  0.4809
   665   0  1.069  5471541    0  14723     415    20  0.0366  0.4797
   666   0  1.062  5470221    0  14720     437    22  0.0361  0.4935
   667   0  1.059  5470216    0  14711     441    20  0.0360  0.4872
   668   0  1.079  5472155    0  14738     449    31  0.0351  0.4651
   669   0  1.078  5471194    0  14744     446    23  0.0373  0.4776
   670   0  1.081  5471839    0  14725     448    19  0.0356  0.4845

  1219 OUT OF   2054 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.29
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.53
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.054   187.885    72.371  90.000 100.179  90.000
 REC. CELL PARAMETERS   0.016373  0.005322  0.014039  90.000  79.821  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.098    18.599    18.914
 COORDINATES OF UNIT CELL B-AXIS   -14.323  -110.582   151.219
 COORDINATES OF UNIT CELL C-AXIS    18.399   -57.317   -40.171
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.051    -0.012     0.013
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1219 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  364    119     59     47     76    100    112    151    186
 SIGMAB (degree)       0.055  0.057  0.055  0.056  0.056  0.059  0.060  0.061  0.060
 SIGMAR (degree)       0.529  0.408  0.541  0.539  0.499  0.431  0.491  0.483  0.411


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1219 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0 -1 -1 -1  0  0
     0  0  0  0  0  0  0  0  0     0  0  1  1  0 -1  0 -1  0     0  0  0  0  0  0 -1 -1  0
     0  0  0  0  0  0  0  0  0     0  0  1  1  0  0  0  0  0     0  0  1  2  1  1  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  2  6  6  3  1  0 -1     0  1  4 12 11  4  2  1  0
     0  0  0  0  0  0  0  0  0     0  0  4 18 26 10  3  1 -1     0  1  8 34 41 11  3  1  0
     0  0  0  0  0  0  0  0  0     0  0  2  3  7  5  2  0 -1     0  1  3  8 13  8  2  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  1  0  0     0 -1 -1  0  2  4  2  0  0
     0  0  0  0  0  0  0  0  0     0  0 -1  0  0  0  0  0  0     0 -1 -1 -1  0  1  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -1  0 -1  0 -1  0  0


     0  0 -1  0 -1 -2 -1  0  0     0  0  0  0  1  0 -1  0  0     0  0  0  0  0 -1 -1  0  0
     0 -1  0  1  0 -1 -1  0  0     0  0  0  1  1  0 -1  0  0     0 -1  0  1  1  0 -1 -1  0
    -1  0  2  3  3  0  1  0  0     0  1  1  4  5  1  1  0  0     0  0  1  3  4  2  0 -1  0
     0  2  6 17 20  5  2  0  0     0  1  4 15 33  9  2  0 -1     0  1  3 10 23 10  3  0  0
     0  2  8 36 63 14  3  1  0     0  1  5 29100 27  4  1  0     1  1  3 18 65 28  6  1  0
     0  1  4 11 18  9  2  1  0    -1  1  3  9 27 15  3  1  0     0  1  2  6 17 15  5  1  0
    -2 -1  0  1  4  4  2  0  0     0  0  0  2  5  5  2  1  0     0  0  0  1  3  4  2  0  0
     0 -1 -1  0  1  1  1  0  0     0 -1 -1  0  1  1  1  0  0     0 -1 -1  0  1  1  1  0  0
     0  0 -2  0  0  0  0  0  0     0  0 -2  0  0  0 -1  0  0     0  0 -1 -1  0  0 -1  0  0


     0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  1  2  2  1  0  0  0     0  0  0  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  1  3  7 14  7  2  1  0     0  1  3  6  8  3  0  0  0     0  0  0  0  0  0  0  0  0
     0  1  3 15 43 23  6  1  0     0  1  3 13 31 14  0  0  0     0  0  0  0  0  0  0  0  0
     0  1  2  5 11 11  4  1  0     0  1  2  3  6  6  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  1  3  2  0  0     0  0  0  0  1  2  0  0  0     0  0  0  0  0  0  0  0  0
     0  0 -1  0  0  0  0 -1  0     0 -1 -1 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0 -1 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30719 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      671 ...     696
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   671   0  1.060  5471720    0  14735     447    16  0.0359  0.5071
   672   0  1.058  5471438    0  14718     426    24  0.0351  0.4557
   673   0  1.085  5471614    0  14794     439    14  0.0386  0.4771
   674   0  1.071  5472654    0  14753     422    18  0.0387  0.4798
   675   0  1.074  5471209    0  14756     449    20  0.0374  0.4666
   676   0  1.067  5471553    0  14759     430    21  0.0397  0.4668
   677   0  1.072  5471144    0  14743     431    24  0.0387  0.4836
   678   0  1.075  5470092    0  14734     435    14  0.0386  0.5048
   679   0  1.073  5471227    0  14789     436    17  0.0365  0.4940
   680   0  1.074  5470948    0  14734     452    23  0.0385  0.4505
   681   0  1.083  5472696    0  14742     440    20  0.0364  0.4734
   682   0  1.074  5472677    0  14785     431    27  0.0375  0.4616
   683   0  1.076  5472172    0  14754     413    18  0.0384  0.4621
   684   0  1.075  5471537    0  14737     451    25  0.0383  0.4581
   685   0  1.071  5471129    0  14775     442    15  0.0369  0.4837
   686   0  1.077  5471630    0  14811     445    13  0.0382  0.5058
   687   0  1.074  5471133    0  14769     439    19  0.0362  0.4622
   688   0  1.071  5472380    0  14779     432    21  0.0376  0.4706
   689   0  1.066  5471404    0  14792     427    19  0.0347  0.4846
   690   0  1.084  5471976    0  14772     431    18  0.0359  0.5123
   691   0  1.073  5472091    0  14782     443    13  0.0359  0.5090
   692   0  1.081  5472879    0  14779     425    17  0.0368  0.4665
   693   0  1.074  5472397    0  14737     441    27  0.0360  0.4426
   694   0  1.086  5473147    0  14735     440    24  0.0354  0.4552
   695   0  1.072  5472591    0  14794     443    18  0.0359  0.4566
   696   0  1.086  5472366    0  14782     438    22  0.0345  0.4227

  1095 OUT OF   2090 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.29
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.55
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.054   188.108    72.381  90.000 100.177  90.000
 REC. CELL PARAMETERS   0.016373  0.005316  0.014037  90.000  79.823  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.100    18.604    18.905
 COORDINATES OF UNIT CELL B-AXIS   -14.409  -110.528   151.528
 COORDINATES OF UNIT CELL C-AXIS    18.396   -57.378   -40.104
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.106    -0.006     0.015
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1095 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  323    125     58     43     76     97     88    108    171
 SIGMAB (degree)       0.053  0.056  0.054  0.055  0.058  0.058  0.057  0.057  0.058
 SIGMAR (degree)       0.530  0.393  0.535  0.534  0.423  0.419  0.507  0.490  0.388


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1095 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  1  1  0  1  0  0     0  0  0  1  1  1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  3  8  8  3  1  0  0     0  1  3 11 13  4  1  1  1
     0  0  0  0  0  0  0  0  0     0  0  3 14 28  7  1  0  0     0  1  5 26 41  9  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  2  2  6  4  1  0  0     0  1  3  6 11  6  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  1  2  1  0  0     0  0  0  1  3  3  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  1  0  0     0  0 -1  0  0  1  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0


     0  0 -1  0  0 -1  0  0  0     0  0  0  1  1  1 -1  0  0     0  0 -1  0  0  0 -1  0  0
     0 -1  0  1  1  0  0  0  0     0 -1  0  2  2  1  0  0  0     0  0  0  1  1  0  0  0  0
     0  0  1  4  4  1  1  0  0    -1 -1  1  4  7  2  1  0 -2    -1  0  0  2  5  2  1  0  0
     0  1  3 14 21  6  1  0  0    -1  0  2 14 37  9  2  0  0     0  0  2  8 24 10  2  0 -1
     0  1  5 27 57 12  1  0  0     0  0  4 24100 23  2  0 -1     0  0  3 15 66 23  2  0 -1
     0  1  4  9 16  7  2  0  0    -1  1  3  8 26 12  2  0 -1     0  0  2  6 17 12  2  0  0
    -1  0  1  2  4  4  2  0  0    -1  0  1  2  6  5  2  0  0     0  0  1  1  3  4  2 -1  0
     0 -1 -1  0  1  1  1  0  0     0  0  0  1  2  1  1  0  0     0  0  0  0  1  1  0  0  0
     0  0 -2  0  0  1 -1  0  0     0  0 -1  0  1  1 -1  0  0     0  0 -1  0  1  1  0  0  0


     0  0 -1  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  2  2  1  0  0  0     0  0  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  2  6 13  6  2  0  0     0  0  2  5  9  3  1  0  0     0  0  0  0  0  0  0  0  0
     1  1  3 12 44 19  2  0  0     0  1  2 10 33 12  2  0  0     0  0  0  0  0  0  0  0  0
     0  1  2  4 10  9  2  0  0     0  1  1  2  5  5  2  0  0     0  0  0  0  0  0  0  0  0
    -1  0  0  0  1  2  1  0  0     0  0  0  0  0  1  1  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0 -1 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0 -1  0 -1  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30788 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      697 ...     722
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   697   0  1.082  5472492    0  14729     437    17  0.0360  0.4224
   698   0  1.082  5472770    0  14744     457    15  0.0374  0.4642
   699   0  1.082  5472249    0  14765     428    12  0.0367  0.4422
   700   0  1.085  5470195    0  14740     472    15  0.0382  0.4598
   701   0  1.088  5471219    0  14712     464    15  0.0364  0.4556
   702   0  1.075  5470924    0  14776     464    12  0.0347  0.4701
   703   0  1.072  5472255    0  14753     441    10  0.0347  0.4867
   704   0  1.094  5472951    0  14739     460    17  0.0347  0.4726
   705   0  1.081  5472285    0  14745     466    11  0.0361  0.4809
   706   0  1.074  5472630    0  14755     445    13  0.0357  0.4761
   707   0  1.079  5471704    0  14810     441    18  0.0351  0.4740
   708   0  1.082  5471203    0  14760     465    16  0.0360  0.4747
   709   0  1.096  5471998    0  14771     450    13  0.0358  0.5027
   710   0  1.086  5472936    0  14737     441    17  0.0365  0.4974
   711   0  1.081  5472601    0  14732     458    20  0.0385  0.5016
   712   0  1.079  5473231    0  14746     447    17  0.0353  0.5079
   713   0  1.104  5473942    0  14755     432    17  0.0340  0.5026
   714   0  1.077  5472565    0  14738     445    26  0.0347  0.4622
   715   0  1.072  5472422    0  14732     444    29  0.0362  0.4646
   716   0  1.074  5471614    0  14750     437    26  0.0347  0.4708
   717   0  1.094  5471455    0  14716     436    25  0.0351  0.4534
   718   0  1.076  5471050    0  14729     463    18  0.0384  0.4752
   719   0  1.074  5471559    0  14792     441    18  0.0345  0.4568
   720   0  1.081  5471167    0  14723     455    23  0.0345  0.4703
   721   0  1.071  5471545    0  14720     455    18  0.0349  0.4748
   722   0  1.074  5471821    0  14734     443    17  0.0333  0.4745

  1133 OUT OF   2202 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.30
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.49
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.059   187.792    72.385  90.000 100.196  90.000
 REC. CELL PARAMETERS   0.016372  0.005325  0.014037  90.000  79.804  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.107    18.618    18.886
 COORDINATES OF UNIT CELL B-AXIS   -14.399  -110.177   151.392
 COORDINATES OF UNIT CELL C-AXIS    18.365   -57.431   -40.050
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.122     0.015     0.016
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1133 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  332    131     69     52     99    109     83     93    161
 SIGMAB (degree)       0.051  0.055  0.056  0.055  0.058  0.057  0.055  0.056  0.056
 SIGMAR (degree)       0.529  0.398  0.490  0.518  0.375  0.392  0.524  0.521  0.392


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1133 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1 -1 -1  0     0  0  0  0  0 -1 -1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  1  1  0  0  0  0     0  0  1  2  1  1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  7  7  2  1  1  0     0  0  4  9 12  4  2  1  0
     0  0  0  0  0  0  0  0  0     0  0  0 13 29  6  2  1  0     0  1  4 20 43  8  2  1  0
     0  0  0  0  0  0  0  0  0     0  0  0  2  6  4  2  0  0     0  1  2  5 10  4  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0     0  0  0  0  1  1  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -1 -1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0


     0  0  0  0 -1 -1  0  0  0     0  0 -1  0  1  0 -1  0  0     0  0  0  0  0  0 -1  0  0
     0 -1  0  1  1  0  0  0  0     0 -1  0  1  2  0  0  0  0     0  0  0  1  2  0 -1  0  0
    -1  0  1  3  3  1  1  0  0    -1 -1  1  4  6  2  1  0  0     0  0  0  3  5  2  0  0  0
    -1  0  3 10 18  5  2  1  0    -1  0  2 10 33  9  2  0 -1     0  0  2  6 22  9  2  0  0
     0  0  4 22 57 12  2  0 -1    -1  0  3 20100 23  2  0 -1     0  0  2 12 64 23  3  0 -1
     0  1  3  8 16  6  1  0  0     0  1  2  8 28 11  1  0 -1     0  1  2  5 17 11  1  0  0
     0  0  1  1  3  2  0  0  0     0  0  1  2  6  3  0  0  0     0  0  1  2  3  3  1  0  0
     0  0 -1  0  1  1  0  0  0     0  0  0  1  2  1  0 -1  0     0  0  0  0  1  1  0 -1  0
     0  0 -1  0  1  0  0  0  0     0  0 -1  0  2  1  0  0  0     0  0 -1  0  0  0 -1  0  0


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0 -1 -1  0  0     0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  2  2  0  0  0  0     0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0
    -1  0  2  5 12  6  2  1  0     0  0  2  4  8  3  1  0  0     0  0  0  0  0  0  0  0  0
     0  0  2  9 42 18  2  1  0     0  0  2  7 32 13  1  0  0     0  0  0  0  0  0  0  0  0
     0  1  2  4 10  8  1  0  0     0  1  2  2  6  5  1  0  0     0  0  0  0  0  0  0  0  0
     0  0  1  1  1  1  1  0  0     0  0  0  0  0  1  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0     0 -1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0 -1 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30683 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      723 ...     747
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   723   0  1.070  5471278    0  14682     454    18  0.0364  0.4531
   724   0  1.092  5473314    0  14760     433     9  0.0342  0.4772
   725   0  1.101  5473139    0  14762     454    15  0.0336  0.4574
   726   0  1.087  5473379    0  14733     456    15  0.0335  0.4550
   727   0  1.098  5473238    0  14749     465    10  0.0331  0.4888
   728   0  1.081  5473057    0  14763     459    11  0.0346  0.4622
   729   0  1.102  5473570    0  14772     471     9  0.0356  0.4769
   730   0  1.088  5471727    0  14774     452    10  0.0358  0.4761
   731   0  1.072  5472368    0  14765     462    17  0.0357  0.4724
   732   0  1.076  5472452    0  14772     483    13  0.0368  0.4796
   733   0  1.097  5472951    0  14730     467    13  0.0361  0.4487
   734   0  1.091  5472969    0  14782     464    15  0.0355  0.4816
   735   0  1.085  5473629    0  14743     463    17  0.0361  0.4611
   736   0  1.078  5472806    0  14763     456    15  0.0350  0.4525
   737   0  1.091  5472291    0  14726     448    16  0.0351  0.4348
   738   0  1.090  5470496    0  14736     472    18  0.0368  0.4454
   739   0  1.095  5470906    0  14772     457    21  0.0359  0.4442
   740   0  1.086  5470159    0  14730     472    15  0.0358  0.4604
   741   0  1.077  5470226    0  14700     481    21  0.0354  0.4650
   742   0  1.102  5471492    0  14681     457    14  0.0354  0.4642
   743   0  1.085  5471054    0  14725     460    11  0.0366  0.4728
   744   0  1.089  5471688    0  14753     480    12  0.0331  0.5035
   745   0  1.088  5471905    0  14761     448    12  0.0344  0.4836
   746   0  1.078  5472513    0  14765     473    16  0.0353  0.4673
   747   0  1.084  5471631    0  14763     455    15  0.0346  0.4894

  1155 OUT OF   2113 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.29
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.48
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.062   187.535    72.341  90.000 100.193  90.000
 REC. CELL PARAMETERS   0.016371  0.005332  0.014045  90.000  79.807  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.112    18.611    18.888
 COORDINATES OF UNIT CELL B-AXIS   -14.427  -109.976   151.217
 COORDINATES OF UNIT CELL C-AXIS    18.349   -57.413   -40.004
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.127     0.029     0.003
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1155 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  345    131     81     72    118    113     76     76    137
 SIGMAB (degree)       0.052  0.055  0.055  0.055  0.057  0.055  0.054  0.055  0.055
 SIGMAR (degree)       0.535  0.407  0.524  0.522  0.386  0.385  0.528  0.529  0.396


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1155 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0 -1  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  1  2  1  0  0  0  0     0  0  1  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  4  4  7  2  0  0  0     0  2  4  7 10  3  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0 13 34  6  1  0  0     0  2  3 20 42  9  2  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  4  9  3  2  1  0     0  1  2  7 13  5  2  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0     0  0  0  0  1  2  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0  0  0     0 -1 -1  0  0  1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0 -1 -1  0  0  0  0


     0  0  0  0  0 -1  0  0  0     0  0  0  1  2  0 -1  0  0     0  0 -1  1  1  0 -1  0  0
     0  0  1  2  1  0  0  0  0     0 -1  0  2  3  1  0 -1  0     0  0  1  1  2  0 -1  0  0
    -1  0  2  3  2  1  1  0  0     0  0  1  4  6  2  1  0 -1     0  0  1  3  4  2  1  0  0
     0  0  3 10 15  5  2  0  0     0  0  2 11 33  9  2  0  0     0  0  2  6 22 10  2  0 -1
     0  1  3 23 58 12  2  0  0     0  0  3 21100 22  2  0 -1     0  0  3 12 63 24  3  0 -1
     0  0  2  9 18  6  2  0  0     0  0  2  8 28 10  2  0 -1     0  1  2  5 16 10  2  0 -1
     0  0  1  2  3  2  1  0  0     0  0  1  2  5  4  1  0  0     0  0  1  1  3  3  1  0  0
     0  0  0  0  1  1  0  0  0     0  0  0  1  2  2  0 -1  0     0  0  0  0  1  1  0  0  0
     0  0 -1 -1 -1  1 -1  0  0     0  0 -1  0  1  1 -1  0  0     0  0 -1  0  1  0  0  0  0


     0  0  0  0 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  1  0  0  0 -1  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  1  2  1  1  0  0  0     0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  2  4 13  8  2  1  0     0  1  2  4  8  3  0  0  0     0  0  0  0  0  0  0  0  0
     0  1  3 10 43 22  3  1  0     0  1  2  7 32 14  0  0  0     0  0  0  0  0  0  0  0  0
     0  1  2  4 10  7  2  1  0     0  1  2  2  5  3  1  0  0     0  0  0  0  0  0  0  0  0
     0  0  1  1  1  1  1  1  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0 -1 -1  0 -1  0  1  0  0     0  0  0  0  0  0  0  0  0
     0  0 -1 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30026 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      748 ...     772
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   748   0  1.075  5472106    0  14705     461    20  0.0355  0.4467
   749   0  1.069  5471522    0  14734     466    19  0.0367  0.4587
   750   0  1.094  5472454    0  14726     472    23  0.0362  0.4547
   751   0  1.092  5472392    0  14704     463    14  0.0351  0.4694
   752   0  1.097  5473207    0  14702     470    18  0.0352  0.4667
   753   0  1.079  5472381    0  14687     465    18  0.0338  0.4456
   754   0  1.085  5473297    0  14670     466    24  0.0346  0.4406
   755   0  1.089  5473679    0  14710     467    18  0.0350  0.4417
   756   0  1.069  5472499    0  14739     480    26  0.0353  0.4367
   757   0  1.083  5472346    0  14704     483    24  0.0354  0.4436
   758   0  1.079  5473539    0  14723     472    14  0.0337  0.4755
   759   0  1.082  5473760    0  14744     478    13  0.0355  0.4794
   760   0  1.076  5473348    0  14740     483    12  0.0338  0.4881
   761   0  1.078  5473342    0  14710     475    14  0.0341  0.4786
   762   0  1.097  5474501    0  14730     457    16  0.0341  0.4644
   763   0  1.081  5474444    0  14711     456    20  0.0355  0.4590
   764   0  1.095  5473941    0  14682     451    27  0.0347  0.4377
   765   0  1.092  5473676    0  14727     466    23  0.0366  0.4419
   766   0  1.082  5472882    0  14731     475    16  0.0333  0.4816
   767   0  1.092  5473069    0  14727     461    19  0.0354  0.4684
   768   0  1.088  5474760    0  14726     467    23  0.0337  0.4540
   769   0  1.084  5473802    0  14717     483    17  0.0328  0.4871
   770   0  1.076  5474097    0  14735     469    17  0.0331  0.4928
   771   0  1.081  5474201    0  14789     487    28  0.0352  0.4304
   772   0  1.085  5473496    0  14786     481    20  0.0335  0.4524

  1225 OUT OF   2075 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.21
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.47
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.066   187.388    72.300  90.000 100.191  90.000
 REC. CELL PARAMETERS   0.016370  0.005337  0.014053  90.000  79.809  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.116    18.613    18.887
 COORDINATES OF UNIT CELL B-AXIS   -14.386  -109.925   151.076
 COORDINATES OF UNIT CELL C-AXIS    18.341   -57.369   -39.996
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.116     0.028     0.005
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1225 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  377    127     97    100    134    114     72     66    128
 SIGMAB (degree)       0.052  0.054  0.057  0.058  0.061  0.057  0.052  0.051  0.054
 SIGMAR (degree)       0.530  0.414  0.479  0.478  0.388  0.382  0.529  0.536  0.393


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1225 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  1  2  0  0  0  0  0     0  1  1  2  1  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  5  9  2  1  1  0     0  1  2  8 10  2  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0 15 30  7  1  0  0     0  0  3 23 38  7  2  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  5  8  5  2  1  0     0  1  2  9 16  6  2  1  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0     0  0  0  0  3  2  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0     0 -1 -1  0  0  1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0 -1  0  0  0  0  0


     0  0  0  0  0 -1 -1  0  0     0  0  0  0  1 -1 -1  0  0     0  0  0  0  1 -1 -1  0  0
     0  0  0  1  0  0  0  0  0     0  0  0  1  2  0  0  0  0     0  0  0  1  1  0 -1 -1  0
     0  0  1  3  2  1  1  0  0     0  0  1  4  6  2  1  0 -1     0  0  1  3  4  2  1  0  0
     0  0  3 11 16  4  2  0  0     0  0  2 12 34  9  2  0 -1     0  1  2  8 24 12  3  0  0
     0  0  3 27 59 11  2  0  0     0  0  3 24100 22  2  0 -1     0  0  2 13 61 28  3  0  0
     0  1  2 11 21  7  2  0  0     0  0  2  8 27 11  2  0 -1     0  0  2  5 15 12  3  0 -1
     0  0  0  2  4  3  1  0  0     0  0  1  2  6  5  1  0  0     0  0  1  1  3  4  2  0  0
     0 -1  0  0  1  1  0  0  0     0 -1  0  1  2  1  0  0  0     0  0  0  0  1  1  0  0  0
     0  0 -1 -1  0  0  0  0  0     0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  1  0  0 -1  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  1  2  2  0  0  0  0     0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  1  5 16 10  3  0  0     1  0  2  4 10  4  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  1  8 37 23  4  1  0     0  0  1  7 28 13  0  0  0     0  0  0  0  0  0  0  0  0
     0  1  1  3  8  7  3  1  0     0  0  0  2  6  3  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0     0  0  0  0  1  1  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  0  0  0     0  0 -1 -1  0  1  1  0  0     0  0  0  0  0  0  0  0  0
     0  0  0 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30077 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      773 ...     798
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   773   0  1.082  5473250    0  14776     505    28  0.0338  0.4400
   774   0  1.078  5473527    0  14733     507    19  0.0360  0.4631
   775   0  1.080  5473751    0  14737     483    21  0.0338  0.4580
   776   0  1.085  5474398    0  14731     494    23  0.0344  0.4701
   777   0  1.081  5474064    0  14720     494    22  0.0339  0.4657
   778   0  1.069  5474101    0  14718     478    22  0.0352  0.4591
   779   0  1.070  5474589    0  14726     491    18  0.0334  0.4814
   780   0  1.080  5474689    0  14703     493    23  0.0330  0.4686
   781   0  1.087  5474518    0  14713     497    16  0.0361  0.4795
   782   0  1.077  5474864    0  14712     508    25  0.0349  0.4781
   783   0  1.066  5473602    0  14750     505    33  0.0325  0.4585
   784   0  1.069  5472732    0  14704     515    29  0.0309  0.4646
   785   0  1.079  5473164    0  14705     511    31  0.0342  0.4670
   786   0  1.082  5474496    0  14714     500    24  0.0326  0.4938
   787   0  1.078  5473769    0  14694     517    19  0.0343  0.5120
   788   0  1.073  5474096    0  14725     494    19  0.0345  0.4845
   789   0  1.079  5474478    0  14697     501    24  0.0333  0.4586
   790   0  1.072  5474863    0  14687     487    25  0.0352  0.4488
   791   0  1.074  5475197    0  14696     479    23  0.0358  0.4563
   792   0  1.062  5474420    0  14695     511    25  0.0342  0.4678
   793   0  1.076  5475287    0  14691     503    25  0.0341  0.4786
   794   0  1.081  5476295    0  14654     487    28  0.0336  0.4652
   795   0  1.072  5476103    0  14652     481    24  0.0341  0.4756
   796   0  1.066  5476100    0  14673     484    24  0.0346  0.4678
   797   0  1.074  5476748    0  14676     473    28  0.0344  0.4020
   798   0  1.094  5476328    0  14684     477    17  0.0352  0.4914

  1427 OUT OF   2254 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.19
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.48
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.070   187.080    72.329  90.000 100.198  90.000
 REC. CELL PARAMETERS   0.016370  0.005345  0.014048  90.000  79.802  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.119    18.622    18.883
 COORDINATES OF UNIT CELL B-AXIS   -14.241  -109.911   150.717
 COORDINATES OF UNIT CELL C-AXIS    18.351   -57.342   -40.083
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.049     0.014     0.017
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1427 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  441    149    116    123    151    125     85     84    148
 SIGMAB (degree)       0.055  0.057  0.060  0.061  0.062  0.058  0.054  0.052  0.053
 SIGMAR (degree)       0.520  0.411  0.461  0.462  0.391  0.396  0.525  0.534  0.410


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1427 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1 -1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  1  0  0  0  0  0  0     0  0  1  1  1  1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  4  9  3  1  0  0     0  1  3  8 11  4  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0 13 22  8  1  0  0     0  1  4 22 32  9  2  1  0
     0  0  0  0  0  0  0  0  0     0  0  0  5  8  4  1  0  0     0  0  2 12 15  6  2  1  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  1  0  0     0  0  0  1  2  2  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0 -1  0  0  0  0     0 -1 -1 -1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0 -1 -1  0  0  0  0


     0  0  0  0  0 -1 -1  0  0     0  0  0  0  1  0 -2  0  0     0  0 -1 -1  0 -1 -2  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  1  2  0 -1 -1  0     0  0  0  1  1  0 -1 -1  0
     0  0  1  3  2  1  0  0  0     0  0  2  4  6  2  0  0  0     0  1  1  3  5  2  0 -1  0
     1  1  4 12 17  5  2  0  0     0  1  3 13 35 10  2  0  0     0  1  2  8 25 14  3  0  0
     0  1  4 28 59 13  2  1  0     0  1  3 24100 26  3  0  0     0  1  3 14 56 30  5  0 -1
     0  1  2 12 22  8  2  1  0     0  1  2  8 26 12  3  0  0     0  0  1  5 14 13  3  1  0
     0  0  0  1  3  3  1  0  0     0  0  1  2  5  4  1  0  0     0  0  0  1  2  3  1  0  0
     0 -1 -1  0  1  1  0  0  0     0 -1  0  0  2  2  1  0  0     0  0  0  0  1  1  0  0  0
     0  0 -2 -1  0  0  0  0  0     0  0 -1  0  1  1  0  0  0     0  0 -1  0  0  1 -1  0  0


     0  0 -1 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0 -1 -1  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  1  2  2  1  0 -1  0     0  0  1  1  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  1  2  6 15 12  3  0  0     1  1  1  4 10  7  0  0  0     0  0  0  0  0  0  0  0  0
     0  1  2  9 30 25  5  1  0     0  1  1  9 23 16  3  0  0     0  0  0  0  0  0  0  0  0
     0  0  1  4  9  7  4  1  0     0  0  0  3  6  4  3  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  1  1  0  0     0  0  0  0  0  1  1  0  0     0  0  0  0  0  0  0  0  0
     0  0  0 -1  0  0  0  0  0     0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0 -1 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30651 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      799 ...     824
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   799   0  1.081  5475953    0  14645     481    19  0.0349  0.4986
   800   0  1.079  5476452    0  14694     487    19  0.0338  0.4990
   801   0  1.073  5476234    0  14736     481    26  0.0340  0.4690
   802   0  1.079  5476988    0  14695     493    24  0.0344  0.4917
   803   0  1.064  5475733    0  14673     501    23  0.0335  0.5053
   804   0  1.066  5474794    0  14685     497    23  0.0339  0.4986
   805   0  1.061  5474674    0  14686     507    24  0.0334  0.4782
   806   0  1.074  5474656    0  14689     505    27  0.0338  0.4835
   807   0  1.068  5476089    0  14746     475    27  0.0339  0.4710
   808   0  1.078  5476644    0  14724     457    22  0.0335  0.4761
   809   0  1.079  5476330    0  14642     493    31  0.0338  0.4838
   810   0  1.072  5476567    0  14679     461    24  0.0356  0.4766
   811   0  1.087  5476511    0  14702     490    22  0.0339  0.4672
   812   0  1.070  5476107    0  14690     474    41  0.0352  0.4166
   813   0  1.079  5477143    0  14731     465    38  0.0353  0.4251
   814   0  1.082  5477146    0  14699     466    29  0.0347  0.4464
   815   0  1.084  5476492    0  14684     482    19  0.0365  0.4661
   816   0  1.080  5475954    0  14716     506    24  0.0336  0.4957
   817   0  1.063  5475623    0  14696     522    27  0.0349  0.4959
   818   0  1.080  5476202    0  14707     509    33  0.0350  0.4643
   819   0  1.085  5477728    0  14725     487    36  0.0339  0.4600
   820   0  1.071  5478171    0  14692     478    33  0.0332  0.4779
   821   0  1.068  5476841    0  14732     497    29  0.0356  0.4873
   822   0  1.060  5476730    0  14691     503    36  0.0355  0.4858
   823   0  1.075  5477176    0  14684     481    31  0.0347  0.4892
   824   0  1.068  5476348    0  14643     481    23  0.0364  0.4851

  1565 OUT OF   2143 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.20
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.46
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.071   187.157    72.301  90.000 100.203  90.000
 REC. CELL PARAMETERS   0.016369  0.005343  0.014053  90.000  79.797  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.120    18.616    18.888
 COORDINATES OF UNIT CELL B-AXIS   -14.260  -109.977   150.763
 COORDINATES OF UNIT CELL C-AXIS    18.336   -57.317   -40.076
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.055     0.026     0.011
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1565 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  491    155    127    131    147    124    105    115    162
 SIGMAB (degree)       0.058  0.059  0.063  0.064  0.064  0.060  0.057  0.056  0.055
 SIGMAR (degree)       0.504  0.438  0.465  0.440  0.406  0.402  0.512  0.521  0.478


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1565 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1 -1  0  0
     0  0  0  0  0  0  0  0  0     0  0  1  1  0  0  0  0  0     0  1  1  1  1  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  2  5  9  4  1  0  0     0  1  3  8 11  5  2  0  0
     0  0  0  0  0  0  0  0  0     0  0  2 12 17  9  2  1  1     0  1  6 20 26 13  3  1  1
     0  0  0  0  0  0  0  0  0     0  0  0  7  7  3  1  1  1     0  0  2 14 14  7  3  1  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  0  0  0     0  0  0  1  2  2  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -1 -1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0 -1 -1  0  0  0  0


     0  0 -1  0  0 -1 -1  0  0     0  0  0  0  0  0 -1  0  0     0  0  0  0  0 -1 -1  0  0
     0  0  0  1  1  0 -1  0  0     0  0  0  1  1  0  0  0  0     0  0  0  1  1 -1 -1 -1  0
    -1  1  2  3  2  1  0  0  0     0  1  2  3  5  2  0  0 -1     0  0  1  3  5  2  0 -1  0
     0  1  5 12 17  6  2  0  0     0  1  4 14 34 11  2  1  0     1  1  3 11 25 15  4  0 -1
     0  1  7 29 56 17  3  1  1     1  1  5 29100 28  4  1  0     1  1  4 18 56 31  6  1  0
     0  1  4 15 23 10  3  1  0     0  1  3  9 26 14  4  1  0     0  1  2  7 18 14  5  1  0
    -1 -1  0  2  3  3  2  1  0     0  0  1  2  4  4  2  0  0    -1  0  1  1  3  4  2  0  0
     0 -2 -1  0  0  1  1  0  0     0  0  0  0  1  1  1  0  0     0  0  0  0  1  1  1  0  0
     0  0 -2 -1 -1  0 -1  0  0     0  0 -1  0  0  1 -1  0  0     0  0 -1  0  0  0 -1  0  0


     0  0 -1 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0 -1  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     1  0  1  2  3  1 -1  0  0     0  0  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  1  2  7 15 12  2  0  0     0  0  2  5  9  8  1  0  0     0  0  0  0  0  0  0  0  0
     0  1  3 13 26 23  5  1  0     0  1  2 10 16 14  3  0  0     0  0  0  0  0  0  0  0  0
     1  0  2  6 12  9  3  1  0     0  0  1  4  7  5  2  0  0     0  0  0  0  0  0  0  0  0
    -1  0  0  0  1  2  1  0  0     0  0  0  0  0  1  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30513 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      825 ...     850
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   825   0  1.054  5475657    0  14694     484    20  0.0350  0.4800
   826   0  1.065  5477299    0  14701     475    22  0.0353  0.4887
   827   0  1.081  5477867    0  14726     467    25  0.0352  0.4532
   828   0  1.083  5477965    0  14711     469    21  0.0343  0.4662
   829   0  1.083  5478315    0  14671     476    23  0.0337  0.4809
   830   0  1.055  5478326    0  14666     462    30  0.0343  0.4746
   831   0  1.058  5479106    0  14740     452    23  0.0336  0.4898
   832   0  1.060  5479262    0  14753     457    23  0.0345  0.5095
   833   0  1.048  5478295    0  14729     474    22  0.0347  0.5163
   834   0  1.066  5478543    0  14691     453    16  0.0355  0.4958
   835   0  1.057  5478249    0  14697     468    28  0.0339  0.4708
   836   0  1.071  5477701    0  14702     477    28  0.0358  0.4518
   837   0  1.077  5477926    0  14716     470    25  0.0344  0.4651
   838   0  1.070  5478025    0  14753     469    26  0.0345  0.4700
   839   0  1.064  5478079    0  14711     464    23  0.0346  0.4786
   840   0  1.051  5477283    0  14683     483    25  0.0348  0.4694
   841   0  1.085  5478399    0  14725     457    21  0.0342  0.5008
   842   0  1.061  5477955    0  14751     482    21  0.0321  0.5092
   843   0  1.077  5478616    0  14688     475    25  0.0337  0.4931
   844   0  1.075  5478789    0  14675     457    28  0.0332  0.5004
   845   0  1.059  5479112    0  14743     463    22  0.0327  0.5125
   846   0  1.053  5478628    0  14770     466    27  0.0344  0.4995
   847   0  1.056  5478287    0  14727     486    30  0.0342  0.5006
   848   0  1.055  5478413    0  14726     487    26  0.0329  0.5106
   849   0  1.057  5478061    0  14735     489    29  0.0342  0.4833
   850   0  1.062  5478603    0  14730     474    24  0.0340  0.5098

  1638 OUT OF   2043 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.17
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.43
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.071   187.424    72.280  90.000 100.190  90.000
 REC. CELL PARAMETERS   0.016369  0.005335  0.014057  90.000  79.810  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.124    18.612    18.881
 COORDINATES OF UNIT CELL B-AXIS   -14.323  -110.036   151.045
 COORDINATES OF UNIT CELL C-AXIS    18.334   -57.327   -40.024
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.098     0.029     0.014
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1638 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  506    160    119    120    137    122    120    151    193
 SIGMAB (degree)       0.060  0.060  0.065  0.067  0.065  0.061  0.059  0.059  0.058
 SIGMAR (degree)       0.479  0.433  0.448  0.422  0.399  0.407  0.483  0.490  0.428


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1638 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  1  0 -1  0  0  0     0  0  0  0  0 -1 -1  0  0
     0  0  0  0  0  0  0  0  0     0  0  1  1  1  0  0  0  0     0  0  1  2  2  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  2  4  7  4  1  0  0     0  1  3  8 13  8  2  0  0
     0  0  0  0  0  0  0  0  0     0  1  2  8 11  8  1  0  0     0  1  6 16 22 16  3  1  0
     0  0  0  0  0  0  0  0  0     0  0  1  6  6  3  1  1  1     0  0  3 13 12  8  3  1  1
     0  0  0  0  0  0  0  0  0     0  0  0  0  1  1  0  0  0     0  0  0  1  2  2  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  1  0  0  0     0  0 -1 -1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0 -1 -1 -1  0  0  0


     0  0  0  0  0 -1 -1  0  0     0  0  0  0  0 -1 -1  0  0     0  0  0  0  0 -2 -2  0  0
     0  0  1  1  0  0 -1 -1  0     0  0  1  1  1  0 -1 -1  0     0  0  0  1  1 -1 -1 -1  0
     0  1  2  3  3  1  0  0  0     0  1  2  4  5  2  1  0  0     0  0  2  3  5  3  0 -1  0
     0  1  6 13 20  8  2  0  0     0  1  4 15 36 11  3  1  0     0  1  4 12 24 15  4  0 -1
     0  1  7 26 54 22  4  1  1     0  1  5 29100 29  4  1  0     1  1  4 21 52 28  6  1  0
    -1  1  4 14 22 12  3  1  0     0  1  3  9 27 16  4  1  0     0  1  3  8 19 14  5  1  0
    -1  0  0  2  3  4  2  1  0     0  0  1  2  4  5  2  1  0     0  0  0  1  3  4  2  0  0
     0 -1 -1  0  0  1  1  0  0     0 -1  0  0  1  1  1  0  0     0  0 -1  0  0  1  1  0  0
     0  0 -1 -1 -1  0 -1  0  0     0  0 -1  0  0  0  0  0  0     0  0 -1  0  0  0 -1  0  0


     0  0 -1  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1  0  0  0 -1 -1  0  0     0 -1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  1  2  3  1  0  0  0     0  0  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  3  8 13 12  2  0  0     0  0  2  4  8  8  1  0  0     0  0  0  0  0  0  0  0  0
     0  1  4 14 20 17  5  1  0     0  1  3  8 10 10  3  0  0     0  0  0  0  0  0  0  0  0
     0  1  2  8 12  7  3  1  0     0  0  1  3  6  4  2  1  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  1  1  1  0  0     0  0 -1 -1  0  1  0  0  0     0  0  0  0  0  0  0  0  0
     0  0 -1 -1  0  0  0  0  0     0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0 -1  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  30652 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      851 ...     875
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   851   0  1.062  5478920    0  14764     471    27  0.0340  0.4761
   852   0  1.059  5478493    0  14664     466    28  0.0351  0.4869
   853   0  1.055  5478015    0  14691     467    31  0.0357  0.4619
   854   0  1.079  5479381    0  14710     447    33  0.0351  0.4110
   855   0  1.052  5479072    0  14660     451    26  0.0350  0.4599
   856   0  1.058  5480467    0  14701     431    17  0.0354  0.4625
   857   0  1.054  5479973    0  14705     447    15  0.0348  0.4943
   858   0  1.062  5480659    0  14650     435    23  0.0347  0.4600
   859   0  1.060  5481888    0  14677     422    19  0.0351  0.4653
   860   0  1.047  5480548    0  14708     465    27  0.0334  0.4879
   861   0  1.083  5481200    0  14654     451    22  0.0353  0.4648
   862   0  1.045  5479699    0  14674     444    37  0.0329  0.4002
   863   0  1.043  5479294    0  14674     466    24  0.0326  0.4586
   864   0  1.055  5480286    0  14685     452    24  0.0352  0.4649
   865   0  1.058  5480352    0  14690     467    32  0.0345  0.4491
   866   0  1.066  5480795    0  14668     449    28  0.0343  0.4818
   867   0  1.061  5481033    0  14640     440    17  0.0338  0.5149
   868   0  1.056  5481081    0  14618     452    14  0.0351  0.5020
   869   0  1.042  5481747    0  14683     435    21  0.0348  0.4848
   870   0  1.046  5481785    0  14670     423    21  0.0360  0.4858
   871   0  1.057  5481997    0  14725     434    19  0.0351  0.5044
   872   0  1.049  5481846    0  14710     440    26  0.0331  0.4877
   873   0  1.043  5481451    0  14673     440    26  0.0345  0.4895
   874   0  1.045  5481219    0  14682     430    35  0.0351  0.4280
   875   0  1.040  5481905    0  14683     441    28  0.0347  0.4937

  1663 OUT OF   1881 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.16
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.42
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.071   187.271    72.292  90.000 100.201  90.000
 REC. CELL PARAMETERS   0.016369  0.005340  0.014055  90.000  79.799  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.125    18.609    18.884
 COORDINATES OF UNIT CELL B-AXIS   -14.301  -109.992   150.890
 COORDINATES OF UNIT CELL C-AXIS    18.323   -57.326   -40.051
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.077     0.039     0.009
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1663 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  495    164    108    105    125    122    139    175    219
 SIGMAB (degree)       0.062  0.064  0.067  0.068  0.067  0.064  0.063  0.062  0.062
 SIGMAR (degree)       0.461  0.414  0.440  0.415  0.372  0.384  0.463  0.469  0.406


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1663 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0     0  0  0  0  0 -1 -1  0  0
     0  0  0  0  0  0  0  0  0     0  0  1  1  1  0  0  0  0     0  0  1  2  3  1  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  2  5  7  5  1  1  0     0  1  3  7 13 10  2  0  0
     0  0  0  0  0  0  0  0  0     0  0  3  8  7  7  2  0  0     0  1  5 12 17 17  4  1  1
     0  0  0  0  0  0  0  0  0     0  0  0  6  5  3  1  1  0     0  0  3 11 10  8  3  1  1
     0  0  0  0  0  0  0  0  0     0  0 -1  0  1  0  0  0  0     0  0  0  1  2  2  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0 -1 -1 -1  0  0  0 -1  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0


     0  0 -1  0  0 -1 -1  0  0     0  0 -1  0  0 -1 -1  0  0     0  0 -1 -1 -1 -1 -2  0  0
     0  0  0  1  1 -1 -1 -1  0     0  0  1  1  1  0 -1 -1  0     0  0  0  0  1 -1 -1 -1  0
     0  1  2  4  5  2  0  0  0     0  1  2  4  5  2  1  0 -1     0  1  2  4  5  3  1  0  0
     0  2  5 13 23 11  3  0  0     1  2  5 17 38 12  3  1  0     0  1  4 13 23 14  4  0  0
     1  2  8 22 49 25  5  1  0     0  2  6 31100 31  5  1  0     1  1  6 24 48 25  7  1  0
     0  1  5 13 22 15  4  1  1     0  1  3 10 29 17  5  1  0     0  1  3 11 21 14  5  1  0
    -1 -1  0  2  3  4  2  0  0    -1  0  1  1  5  5  2  0  0     0  0  0  1  3  4  2  0  0
     0 -2 -1 -1  0  1  1  0  0     0 -1 -1  0  1  2  1  0  0     0 -1 -1 -1  0  0  0  0  0
     0  0 -2 -1 -1  0 -1  0  0     0  0 -1 -1  0  0 -1  0  0     0  0 -1 -1  0  0 -1  0  0


     0  0 -1 -1  0 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1  0  0  0 -1 -1  0  0     0  0  0  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  1  2  3  1  0  0  0     0  0  1  1  1  0  0  0  0     0  0  1  0  0  0  0  0  0
     0  1  3  8 11 10  2  0  0     0  1  2  4  6  5  1  0  0     0 -1  0  0  0  0  0  0  0
     1  1  5 15 15 13  5  1  0     0  1  2  7  7  7  3  0  0     0  0  0  0  0  0  0  0  0
     0  1  2 10 12  7  3  1  0     0  0  1  5  6  4  2  1  0     0  0  0  0  0  0  0  0  0
     0  0  0  1  2  2  1  0  0     0  0 -1  0  1  1  1  0  0     0  0  0  0  0  0  0  0  0
     0 -1 -1  0  0  0  0  0  0     0  0  0  0  0  0  1  0  0     0  0  0  0  0  0  0  0  0
     0  0  0 -1 -1 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  29986 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
                     PROCESSING OF IMAGES      876 ...     900
 ******************************************************************************

 USING   1 PROCESSORS

 *** DEFINITION OF SYMBOLS ***
 IER     = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT
 SCALE   = SCALING FACTOR FOR THIS DATA IMAGE
 NBKG    = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NOVL    = NUMBER OF OVERLOADED REFLECTIONS ON DATA IMAGE
 NEWALD  = NUMBER OF REFLECTIONS CLOSE TO THE EWALD SPHERE
 NSTRONG = NUMBER OF STRONG REFLECTIONS ON DATA IMAGE
 NREJ    = NUMBER OF UNEXPECTED REFLECTIONS
 SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB
 SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)

 IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR
   876   0  1.050  5482179    0  14642     421    26  0.0346  0.4949
   877   0  1.043  5482629    0  14637     421    21  0.0343  0.5169
   878   0  1.041  5482644    0  14634     448    23  0.0322  0.5262
   879   0  1.043  5482163    0  14668     429    28  0.0350  0.4870
   880   0  1.046  5482275    0  14648     417    30  0.0351  0.4250
   881   0  1.053  5481135    0  14670     444    23  0.0351  0.4674
   882   0  1.043  5481681    0  14632     421    29  0.0368  0.4437
   883   0  1.044  5481500    0  14633     420    22  0.0352  0.4771
   884   0  1.039  5482311    0  14644     413    15  0.0346  0.4752
   885   0  1.038  5482912    0  14656     392    18  0.0357  0.4691
   886   0  1.045  5482322    0  14638     399    17  0.0362  0.4573
   887   0  1.041  5484074    0  14668     385    14  0.0340  0.4970
   888   0  1.049  5484453    0  14658     386    18  0.0352  0.4861
   889   0  1.033  5484998    0  14705     397    15  0.0360  0.5071
   890   0  1.042  5485622    0  14714     381    22  0.0354  0.4771
   891   0  1.033  5484943    0  14718     393    22  0.0353  0.4956
   892   0  1.046  5485506    0  14691     385    14  0.0350  0.5008
   893   0  1.046  5484897    0  14728     391    11  0.0366  0.5075
   894   0  1.030  5484329    0  14760     382    22  0.0366  0.4823
   895   0  1.043  5483782    0  14754     387    26  0.0346  0.4444
   896   0  1.043  5483321    0  14734     405    25  0.0351  0.4108
   897   0  1.039  5483495    0  14755     414    23  0.0358  0.4696
   898   0  1.028  5482972    0  14685     404    20  0.0352  0.4986
   899   0  1.022  5483134    0  14695     405    22  0.0361  0.4833
   900   0  1.034  5483096    0  14726     398    27  0.0363  0.4975

  1552 OUT OF   1725 REFLECTIONS ACCEPTED FOR REFINEMENT
 REFINED PARAMETERS:   ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     1.20
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.40
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     62.075   187.121    72.298  90.000 100.201  90.000
 REC. CELL PARAMETERS   0.016368  0.005344  0.014054  90.000  79.799  90.000
 COORDINATES OF UNIT CELL A-AXIS    56.133    18.596    18.884
 COORDINATES OF UNIT CELL B-AXIS   -14.275  -109.966   150.724
 COORDINATES OF UNIT CELL C-AXIS    18.313   -57.315   -40.082
 CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION     0.053     0.051     0.004
 shown as x,y,z components of rotation axis X angle (degrees)
 CRYSTAL MOSAICITY (DEGREES)     0.427
 LAB COORDINATES OF ROTATION AXIS  0.999941 -0.000287 -0.010885
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.009950  0.003805  1.020884
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1244.12   1209.77
 DETECTOR ORIGIN (PIXELS) AT                     1222.30   1201.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       384.98
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 STANDARD DEVIATIONS OF BEAM DIVERGENCE AND REFLECTING RANGE OBTAINED
 FROM    1552 REFLECTION PROFILES AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0
 NUMBER                  465    159     96     91    106    110    130    168    213
 SIGMAB (degree)       0.063  0.064  0.068  0.070  0.069  0.066  0.065  0.065  0.063
 SIGMAR (degree)       0.441  0.412  0.422  0.402  0.353  0.393  0.442  0.451  0.415


 ***** AVERAGE THREE-DIMENSIONAL PROFILE OF   1552 STRONG REFLECTIONS *****


     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0 -1  0  0  0
     0  0  0  0  0 -1  0  1  0     0  0  0  0  1 -1 -1  0  0     0  0  0  0  0 -1 -1 -1  0
     0  0  0  0  0  0  0  0  0     0  0  1  2  1  0  0  0  0     0  0  1  3  3  1  0  0  0
     0  0  0  0  0  0  0  0  0     0  1  2  4  6  5  1  0  0     1  1  3  7 13 11  3  1  0
     0  0  0  0  0  0  1  0  0     0  0  2  5  5  6  3  0  0     1  1  5  9 14 17  6  1  0
     0  0  0  0  0  0  0  0  0     0  0  1  5  4  3  1  0  0     0  0  2  9  9  8  4  1  0
     0  0  0  0  0  0  0  0  0     0  0 -1  0  1  1  1  0  0     0  0  0  1  2  2  1  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -1 -1  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0 -1  0  0  0  0  0  0


     0  0 -1  0  0 -1 -2  0  0     0  0 -1  0 -1 -1 -2  0  0     0  0 -1 -1 -1 -1 -2  0  0
     0  0  0  1  0 -1 -1 -1  0     0  0  1  1  1  0 -1 -1  0     0  0  0  0  0 -1 -1 -2  0
     0  0  2  4  5  2  0  0  0     1  1  2  4  5  2  1  0  0     0  1  2  4  5  3  1 -1  0
     0  1  5 13 23 12  3  1 -1     1  2  5 17 38 13  3  0  0     0  2  5 13 22 14  4  0 -1
     1  2  7 21 45 27  6  1  0     1  2  6 33100 32  6  1  0     0  2  6 25 45 24  7  1  0
     0  1  5 13 21 16  5  1  0     0  1  3 11 30 19  5  1  1    -1  1  4 13 22 14  6  1  1
    -1 -1  0  2  4  5  3  1  0    -1  0  1  1  5  6  3  1  0     0  0  1  2  4  4  2  0  0
     0 -2 -1 -1  0  1  1  0  0     0 -1 -1  0  1  2  1  1  0     0 -1 -1 -1  0  0  1  0  0
     0  0 -2 -2  0  0 -1  0  0     0  0 -2 -1  0  0 -1  0  0     0  0 -2 -1 -1 -1 -1  0  0


     0  0 -1  0  0 -1  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0 -1 -1  0  0 -1 -1  0  0     0 -1  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  1  2  3  2  0  0  0     0  0  1  1  1  0  0  0  0     0  0  0  0  0  0  0  0  0
    -1  1  3  7  9  9  2  0  0     0  0  1  3  4  4  1  0  0     0  0  0  0  0  0  0  0  0
     0  1  6 15 12 11  4  0  0     0  1  3  6  4  4  2  0  0     0  0  1  0  0  0  0  0  0
     0  1  4 12 12  6  3  1  1     0  0  1  4  5  3  2  1  0     0  0 -1  0  0  0  0  0  0
     0  0  0  1  2  2  1  0  0     0  0 -1  0  1  1  0  0  0     0  0 -1  0  0  0  0  0  0
     0 -1 -1 -1  0  0  0  0  0     0  0 -1 -1  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0 -1 -1 -1  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0

  29926 REFLECTION INTENSITIES INTEGRATED


 ******************************************************************************
  REFERENCE PROFILES AS FUNCTION OF DETECTOR POSITION
 ******************************************************************************

 PROFILES ARE DETERMINED AT 9 POSITIONS ON THE DETECTOR SURFACE.
 POSITION NUMBER        1      2      3      4      5      6      7      8      9
 X-COORDINATE (pixel) 1231.0 2148.9 1880.1 1231.0  582.0  313.1  581.9 1231.0 1880.0
 Y-COORDINATE (pixel) 1263.0 1263.0 1928.9 2204.8 1928.9 1263.0  597.1  321.2  597.0




 ***** RUN-AVERAGE OF PROFILE #  1 *****
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  5  4  3  0  0  0
     0  0  0  0  5  0  0  0  0     0  0  0  6  7  5  0  0  0     0  0  5  9 12  8  4  0  0
     0  0  0  7 12  7  0  0  0     0  0  4 10 15  7  0  0  0     0  0  6 14 23 14  6  0  0
     0  0  0  0  6  0  0  0  0     0  0  0  5  6  5  0  0  0     0  0  3  8 12  9  5  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  3  4  4  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  4  6  5  4  3  0  0     0  0  5  6  6  5  4  0  0     0  0  5  5  5  4  3  0  0
     0  0  6 14 21 11  5  0  0     0  0  6 17 32 14  5  0  0     0  0  6 14 23 13  5  0  0
     0  0  7 26 56 25  7  3  0     0  3  7 33100 35  7  3  0     0  0  7 24 58 27  7  3  0
     0  0  5 12 22 15  6  0  0     0  0  6 15 31 19  7  0  0     0  0  5 12 22 16  6  0  0
     0  0  0  4  5  6  5  0  0     0  0  4  6  5  7  6  0  0     0  0  0  4  5  5  3  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  4  4  3  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  5  9 13  8  4  0  0     0  0  0  5  7  5  0  0  0     0  0  0  4  7  0  0  0  0
     0  0  6 14 23 14  5  0  0     0  0  3  7 15 10  3  0  0     0  0  0  7 11 10  0  0  0
     0  0  4  8 12  9  4  0  0     0  0  0  4  6  5  0  0  0     0  0  0  0  4  0  0  0  0
     0  0  0  3  5  4  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  2 *****
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  4  4  3  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  7  6  4  0  0  0     0  0  7  8 10  7  3  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  6  7  7  0  0  0     0  0  7 12 18  9  3  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  4  4  3  0  0  0     0  0  4  9 11  6  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  3  4  3  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  3  3  0  0  0  0     0  0  0  4  7  3  0  0  0     0  0  0  0  4  0  0  0  0
     0  0  5  7  8  6  4  0  0     0  0  4  9 14  8  4  0  0     0  0  4  7 10  6  3  0  0
     0  0  6 15 23 10  4  0  0     0  0  4 16 46 18  4  0  0     0  0  3  9 27 17  6  0  0
     0  0  7 25 43 14  3  0  0     0  0  4 29100 32  4  0  0     0  0  3 13 46 26  7  0  0
     0  0  7 16 23 11  4  0  0     0  0  6 20 41 20  6  0  0     0  0  4 11 23 15  7  0  0
     0  0  4  6  9  8  4  0  0     0  0  6  9 12 10  6  0  0     0  0  4  7  8  6  4  0  0
     0  0  0  2  4  4  0  0  0     0  0  0  4  4  3  0  0  0     0  0  0  3  3  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  5  4  3  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  6 11 10  4  0  0     0  0  0  3  4  3  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  4 11 18 12  7  0  0     0  0  0  9  8  6  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  4  7  9  8  6  0  0     0  0  0  4  5  6  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  3  4  3  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  3 *****
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  9  9  6  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  9 10  8  0  0  0     0  0  7 12 15 10  5  0  0
     0  0  0  5  0  0  0  0  0     0  0  5  8 16  9  4  0  0     0  0  9 15 21 13  7  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  3  9  5  0  0  0     0  0  6 12 13  8  4  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  4  4  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  4  5  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  7 10 10  6  0  0     0  0  3  5  7  9  9  0  0     0  0  0  3  5  7  6  0  0
     0  0  7 19 29 18  8  0  0     0  0  5 14 36 24 10  4  0     0  0  4  8 19 18 10  0  0
     0  5 13 34 51 22  7  6  0     0  4 10 42100 39  9  4  0     0  5  8 24 57 33 11  5  0
     0  6 11 17 17  9  5  0  0     0  5 12 30 36 15  5  0  0     0  6  8 22 33 17  6  0  0
     0  0  5  7  4  3  0  0  0     0  0 10 12  7  5  3  0  0     0  0  7 14 11  5  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  5  5  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  4  8 16 14  6  0  0     0  0  0  6 12  3  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  6 14 23 16  8  0  0     0  0  4 10 21  9  5  0  0     0  0  0  0  4  0  0  0  0
     0  0  6 12 17 13  8  0  0     0  0  0  8  9  8  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  7 11  8  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  4 *****
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0 10 11  8  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  4  7  9  9  0  0  0     0  0  8 14 18 13  7  0  0
     0  0  0  0  4  0  0  0  0     0  0  4  9 16 11  0  0  0     0  0  8 15 21 15  8  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  6 11  6  0  0  0     0  0  6 14 15 12  7  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  5  4  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  7  8  9  7  4  0  0     0  0  7  8  5  6  5  0  0     0  0  6  6  4  3  3  0  0
     0  0 10 22 29 14  8  0  0     0  4 10 23 34 15  7  4  0     0  5  8 16 19 11  8  3  0
     0  5 11 30 53 28 10  4  0     0  5 10 40100 42 11  5  0     0  3  9 29 64 34 11  4  0
     0  0  7 10 18 15  9  3  0     0  5  7 16 33 26 11  4  0     0  0  7 16 34 26 11  3  0
     0  0  0  2  4  6  6  0  0     0  0  5  5  4  7  7  0  0     0  0  4  6  7  8  5  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  4  5  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  6  9 13 13  5  0  0     0  0  3  6 11  6  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  7 16 26 17  8  0  0     0  0  4  9 18 10  0  0  0     0  0  0  0  3  0  0  0  0
     0  0  6 14 21 15  7  0  0     0  0  0  9 10  8  5  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  7 11  9  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  5 *****
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  6 10  0  0  0  0     0  0  0  8  9  7  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  6  7 11  0  0  0     0  0  5 10 13 11  0  0  0
     0  0  0  0  2  0  0  0  0     0  0  0  6  9  7  0  0  0     0  0  8 15 16 13  8  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  5  6  0  0  0     0  0  0 10 10  8  6  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  4  5  5  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  8  6  0  0  0  0     0  0  0  7  5  0  0  0  0     0  0  0  3  0  0  0  0  0
     0  0  6 11 11  4  0  0  0     0  0  7 13 11  3  0  0  0     0  0  6  9  8  3  0  0  0
     0  0  7 18 26 12  4  0  0     0  4  7 23 43 10  3  0  0     0  0  7 17 23  8  3  0  0
     0  0  6 19 41 21  8  0  0     0  0  4 27100 33  7  3  0     0  0  6 21 54 24  7  0  0
     0  0  3  7 19 15  8  4  0     0  0  2 10 41 29  9  5  0     0  0  3 10 29 25  9  3  0
     0  0  0  3  7  8  6  0  0     0  0  0  4 11 13  8  0  0     0  0  0  3  9 12  7  0  0
     0  0  0  0  0  3  0  0  0     0  0  0  0  4  6  0  0  0     0  0  0  0  5  6  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  4  4  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  7  9  9 10  0  0  0     0  0  0  5  5  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  7 16 20 13  8  0  0     0  0  0  8  8  4  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  5 11 15 12  7  0  0     0  0  0  8  9  6  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  5  9  8  3  0  0     0  0  0  0 10  4  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  6 *****
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  5  0  0  0  0     0  0  0  4  6  4  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  4  5  6  0  0  0     0  0  2  8 12  9  4  0  0
     0  0  0  0  3  0  0  0  0     0  0  0  4  7  8  0  0  0     0  0  6 17 17 12  5  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  3  5  7  0  0  0     0  0  4  8 10  9  4  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  4  0  0  0  0     0  0  0  3  4  4  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  6  5  0  0  0  0     0  0  0  5  6  4  0  0  0     0  0  0  2  4  3  0  0  0
     0  0  3  7  9  6  0  0  0     0  0  4  9 13  6  3  0  0     0  0  4  8  9  5  0  0  0
     0  0  4 11 23 12  4  0  0     0  0  3 16 45 14  3  0  0     0  0  4 13 24 11  4  0  0
     0  0  5 18 40 17  4  0  0     0  0  3 28100 25  3  0  0     0  0  4 21 46 18  4  0  0
     0  0  5 12 20 11  5  0  0     0  0  4 19 46 18  5  0  0     0  0  4 13 24 13  5  0  0
     0  0  0  5  8  7  4  0  0     0  0  4  8 14  9  5  0  0     0  0  4  6  9  8  3  0  0
     0  0  0  3  3  2  0  0  0     0  0  0  3  6  4  0  0  0     0  0  0  0  5  4  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  5  4  3  0  0  0     0  0  0  0  6  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  6 11 10  8  4  0  0     0  0  0  8  5  4  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  6 13 17 13  5  0  0     0  0  0  8  7  4  0  0  0     0  0  0  0  3  0  0  0  0
     0  0  4  8 12  7  3  0  0     0  0  0  6  5  4  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  3  6  4  0  0  0     0  0  0  0  6  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  7 *****
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  5  4  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  6 10  6  0  0  0     0  0  4  9 12 10  7  0  0
     0  0  0  0  4  7  0  0  0     0  0  3  9 15  9  5  0  0     0  0  7 20 24 17  7  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  6  9  9  0  0  0     0  0  6 11 16 11  5  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  7 10  6  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  3  6  0  0  0     0  0  0  0  4  0  0  0  0
     0  0  0  3  5  8  7  0  0     0  0  0  4  8  9  8  0  0     0  0  0  6 10 11  5  0  0
     0  0  4  9 20 17  8  0  0     0  0  4 13 37 22  8  4  0     0  0  5 16 29 18  7  0  0
     0  4  9 27 58 29  8  3  0     0  4  8 40100 32  7  4  0     0  4  9 32 52 22  7  0  0
     0  5  8 20 29 13  4  0  0     0  5  9 27 38 12  4  0  0     0  4 10 18 19  8  5  0  0
     0  0  4 11 10  5  0  0  0     0  0  8 13  9  4  0  0  0     0  0  7  8  5  3  0  0  0
     0  0  0  5  3  0  0  0  0     0  0  0  5  3  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  7 10  7  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  6 14 16 11  5  0  0     0  0  0  9  8  7  0  0  0     0  0  0  0  0  0  0  0  0
     0  0 10 17 21 15  6  0  0     0  0  6  8 15 10  0  0  0     0  0  0  5  4  0  0  0  0
     0  0  7 10 11  8  4  0  0     0  0  3  7 10  5  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  4  3  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  8 *****
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  5  5  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  7 13  8  0  0  0     0  0  6 11 13  8  5  0  0
     0  0  0  0  0  0  0  0  0     0  0  3  9 17  9  3  0  0     0  0  8 17 26 16  7  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  7  9  8  0  0  0     0  0  6 13 19 16  7  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  9 11  8  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  4  5  4  4  4  0  0     0  0  6  7  5  6  5  0  0     0  0  6  7  8  7  5  0  0
     0  0  8 14 16 10  8  0  0     0  4  9 20 33 16  8  5  0     0  0  9 21 33 18  9  0  0
     0  4 10 30 63 32 10  4  0     0  4 10 39100 42 11  4  0     0  5 10 29 55 29 10  4  0
     0  0  8 16 33 24  9  0  0     0  4  7 16 31 24 10  4  0     0  0  8 10 15 13  9  4  0
     0  0  4  7  7  8  6  0  0     0  0  5  6  4  7  7  0  0     0  0  4  3  3  4  4  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  8 11  9  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  7 15 19 13  8  0  0     0  0  0  7  8  7  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  8 15 21 14  8  0  0     0  0  4 11 17  8  7  0  0     0  0  0  0  2  0  0  0  0
     0  0  6  9 15 10  5  0  0     0  0  0  6 11  5  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  4  5  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0


 ***** RUN-AVERAGE OF PROFILE #  9 *****
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  4  4  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  4  7  9  0  0  0     0  0  6 10 10  5  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  6  8  4  0  0  0     0  0  8 13 20  9  5  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  5  7  5  0  0  0     0  0  5 10 15 11  7  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  6  9  7  3  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  6  6  0  0  0  0     0  0  0  6  6  0  0  0  0
     0  0  5  9  7  2  0  0  0     0  0  7 11 12  3  0  0  0     0  0  6 11 12  6  0  0  0
     0  0  8 16 17  5  3  0  0     0  4  7 22 43 12  4  0  0     0  0  6 17 33 16  6  0  0
     0  0  8 23 49 18  6  0  0     0  3  6 28100 37  8  3  0     0  0  5 15 45 27 10  5  0
     0  0  5 12 31 22  8  0  0     0  0  4 10 38 30 10  4  0     0  0  3  5 15 16 10  6  0
     0  0  0  4 11 13  7  0  0     0  0  0  4 10 13  8  0  0     0  0  0  2  4  7  6  0  0
     0  0  0  0  4  6  0  0  0     0  0  0  0  3  4  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
 
 
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  6  8  8  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  7 10 14 10  7  0  0     0  0  0  6  7  6  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  5 10 16 11  8  0  0     0  0  0  5  8  6  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  2  6 11  9  6  0  0     0  0  0  7  7  4  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  5  3  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0
     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0


 STANDARD DEVIATION (DEGREES) OF BEAM DIVERGENCE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.049  0.050  0.053  0.055  0.055  0.049  0.050  0.051  0.050
     2    0.051  0.052  0.054  0.057  0.056  0.052  0.051  0.054  0.052
     3    0.052  0.054  0.055  0.057  0.056  0.052  0.052  0.055  0.054
     4    0.053  0.055  0.055  0.056  0.056  0.054  0.054  0.055  0.055
     5    0.053  0.057  0.056  0.056  0.056  0.055  0.055  0.056  0.056
     6    0.055  0.058  0.057  0.057  0.057  0.056  0.056  0.057  0.058
     7    0.055  0.058  0.056  0.056  0.057  0.057  0.056  0.057  0.058
     8    0.055  0.060  0.057  0.057  0.057  0.056  0.056  0.059  0.060
     9    0.055  0.059  0.057  0.057  0.057  0.057  0.057  0.059  0.060
    10    0.057  0.060  0.058  0.058  0.058  0.057  0.056  0.058  0.059
    11    0.057  0.060  0.060  0.058  0.060  0.058  0.056  0.059  0.061
    12    0.056  0.061  0.059  0.057  0.060  0.059  0.056  0.058  0.060
    13    0.055  0.060  0.058  0.057  0.060  0.058  0.055  0.057  0.060
    14    0.054  0.059  0.057  0.054  0.058  0.058  0.056  0.057  0.059
    15    0.053  0.060  0.055  0.053  0.057  0.057  0.054  0.056  0.060
    16    0.052  0.058  0.056  0.052  0.056  0.055  0.053  0.055  0.058
    17    0.050  0.056  0.053  0.051  0.055  0.054  0.053  0.054  0.057
    18    0.049  0.054  0.056  0.051  0.052  0.050  0.052  0.053  0.055
    19    0.048  0.052  0.055  0.051  0.051  0.050  0.050  0.051  0.053
    20    0.047  0.050  0.051  0.050  0.048  0.048  0.049  0.050  0.052
    21    0.046  0.048  0.049  0.052  0.050  0.046  0.051  0.051  0.050
    22    0.045  0.047  0.049  0.049  0.048  0.045  0.051  0.050  0.049
    23    0.044  0.045  0.046  0.047  0.046  0.045  0.049  0.049  0.048
    24    0.042  0.043  0.044  0.045  0.044  0.043  0.046  0.049  0.048
    25    0.040  0.041  0.042  0.042  0.041  0.044  0.046  0.045  0.045
    26    0.038  0.039  0.040  0.041  0.040  0.042  0.044  0.044  0.045
    27    0.036  0.039  0.038  0.040  0.042  0.041  0.041  0.041  0.041
    28    0.034  0.037  0.040  0.038  0.041  0.039  0.040  0.040  0.039
    29    0.035  0.037  0.039  0.039  0.040  0.038  0.039  0.040  0.039
    30    0.036  0.037  0.041  0.042  0.043  0.040  0.036  0.035  0.039
    31    0.039  0.040  0.043  0.044  0.045  0.041  0.038  0.037  0.038
    32    0.042  0.042  0.046  0.047  0.048  0.042  0.040  0.041  0.039
    33    0.043  0.043  0.048  0.050  0.049  0.043  0.043  0.043  0.042
    34    0.045  0.046  0.051  0.052  0.051  0.045  0.046  0.045  0.045
    35    0.046  0.046  0.052  0.054  0.052  0.047  0.047  0.048  0.048


 STANDARD DEVIATION (DEGREES) OF REFLECTING RANGE
 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1    0.489  0.362  0.339  0.405  0.454  0.435  0.411  0.491  0.486
     2    0.477  0.330  0.328  0.392  0.430  0.435  0.397  0.462  0.446
     3    0.441  0.310  0.310  0.366  0.393  0.389  0.377  0.425  0.408
     4    0.421  0.297  0.297  0.359  0.393  0.363  0.373  0.404  0.377
     5    0.379  0.281  0.288  0.352  0.380  0.365  0.354  0.394  0.346
     6    0.377  0.271  0.291  0.356  0.378  0.301  0.352  0.387  0.331
     7    0.363  0.261  0.294  0.360  0.379  0.291  0.349  0.387  0.317
     8    0.350  0.266  0.304  0.381  0.374  0.283  0.326  0.358  0.299
     9    0.348  0.267  0.311  0.384  0.364  0.289  0.338  0.368  0.299
    10    0.408  0.258  0.321  0.388  0.400  0.307  0.371  0.411  0.303
    11    0.402  0.283  0.347  0.389  0.356  0.331  0.368  0.385  0.317
    12    0.440  0.288  0.342  0.423  0.409  0.384  0.412  0.402  0.333
    13    0.470  0.306  0.379  0.441  0.414  0.362  0.448  0.432  0.348
    14    0.489  0.333  0.435  0.474  0.456  0.471  0.483  0.438  0.361
    15    0.506  0.360  0.438  0.495  0.487  0.449  0.463  0.449  0.363
    16    0.502  0.358  0.451  0.506  0.478  0.473  0.467  0.451  0.369
    17    0.508  0.380  0.465  0.516  0.473  0.462  0.506  0.469  0.390
    18    0.515  0.429  0.494  0.522  0.475  0.483  0.514  0.489  0.416
    19    0.499  0.434  0.492  0.500  0.463  0.446  0.499  0.469  0.406
    20    0.496  0.448  0.517  0.504  0.456  0.469  0.514  0.471  0.418
    21    0.517  0.470  0.556  0.537  0.445  0.501  0.509  0.465  0.407
    22    0.513  0.471  0.566  0.540  0.409  0.474  0.530  0.482  0.409
    23    0.551  0.470  0.613  0.587  0.462  0.504  0.552  0.506  0.429
    24    0.571  0.475  0.625  0.610  0.504  0.493  0.557  0.482  0.429
    25    0.610  0.411  0.613  0.600  0.504  0.457  0.539  0.514  0.409
    26    0.691  0.405  0.682  0.623  0.578  0.443  0.554  0.542  0.415
    27    0.733  0.384  0.638  0.604  0.417  0.402  0.585  0.532  0.378
    28    0.719  0.392  0.554  0.549  0.346  0.370  0.627  0.639  0.369
    29    0.708  0.400  0.589  0.604  0.369  0.360  0.605  0.659  0.367
    30    0.644  0.410  0.538  0.515  0.383  0.361  0.593  0.662  0.365
    31    0.672  0.419  0.512  0.498  0.385  0.364  0.626  0.620  0.383
    32    0.583  0.450  0.510  0.477  0.393  0.373  0.574  0.628  0.552
    33    0.539  0.442  0.496  0.456  0.383  0.412  0.545  0.545  0.468
    34    0.519  0.419  0.467  0.420  0.354  0.359  0.508  0.522  0.407
    35    0.495  0.419  0.440  0.410  0.339  0.372  0.498  0.505  0.423


 ***** SUGGESTED VALUES FOR INPUT PARAMETERS *****
 BEAM_DIVERGENCE=   0.502  BEAM_DIVERGENCE_E.S.D.=   0.050
 REFLECTING_RANGE=  3.091  REFLECTING_RANGE_E.S.D.=  0.442



 NUMBER OF STRONG REFLECTIONS USED FOR LEARNING REFERENCE PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1      515    141    124     95     96    159    233    200    175
     2      480    135    111     82     89    149    211    177    161
     3      447    131    100     71     77    138    190    157    150
     4      433    121     93     67     78    133    176    140    135
     5      378    105     81     59     72    122    141    113    113
     6      375    102     83     60     69    118    139    106    105
     7      356     86     67     54     73    115    127     97     94
     8      353     88     67     57     72    105    113     93     98
     9      361     85     69     59     77    111    110     91     93
    10      323     78     61     52     70     96     99     85     89
    11      348     84     63     55     77     97     98     90     97
    12      323     78     61     53     76     88     78     77     85
    13      348     80     59     57     81     99     87     81     89
    14      359     76     57     59     87     99     92     87     92
    15      350     81     57     61     91     94     80     85     98
    16      363     83     61     67    103    109     84     85    100
    17      372     79     59     70    113    115     85     91    103
    18      400     89     74     84    126    115     83     95    110
    19      402     95     76     86    141    128     83     99    119
    20      392     92     71     84    141    128     86    104    123
    21      434    108     84    103    155    129     88    115    146
    22      449    119     93    105    150    129    102    128    158
    23      436    114     87     98    127    112    102    137    163
    24      384    110     71     71     87     89     92    137    168
    25      376    112     61     55     79    102    116    149    178
    26      364    119     59     47     76    100    112    151    186
    27      323    125     58     43     76     97     88    108    171
    28      332    131     69     52     99    109     83     93    161
    29      345    131     81     72    118    113     76     76    137
    30      377    127     97    100    134    114     72     66    128
    31      441    149    116    123    151    125     85     84    148
    32      491    155    127    131    147    124    105    115    162
    33      506    160    119    120    137    122    120    151    193
    34      495    164    108    105    125    122    139    175    219
    35      465    159     96     91    106    110    130    168    213
 margin   13899   3891   2820   2648   3577   4013   3904   4005   4757



 NUMBER OF GRID POINTS IN REFERENCE PROFILES WITH  SIGNAL >   2.% OF MAXIMUM

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1      116     98     85    104    106    102     87    101     99
     2      125     94     84    110    107    110     95    106     97
     3      126     96     86    107    104    106     89    101     94
     4      126     95     87    103    102    100     94    105     89
     5      121     94     88    107    100    100     95    107     89
     6      125     98     89    106    104    102     97    110     93
     7      124     96     93    100     99    100     98    109     97
     8      119    106     94    104    102     99     99    107     96
     9      114    103     92    107    102     97    102    109     94
    10      130    100     93    106    102     98     99    108     93
    11      123     98    105    103    102    104     96    105     94
    12      130    107    101    110    103    105    103    105    100
    13      129    109    105    110    109    104    110    112    101
    14      129    114    107    105    108    117    109    113    100
    15      128    115    104    106    103    115    107    111    104
    16      129    107    112    105    112    112    111    111    101
    17      124    111    111    102    112    115    114    110    101
    18      124    118    119    107    101    110    115    110    104
    19      121    113    113    108     93    106    113    111     99
    20      116    115    109    100     89    101    110    105    102
    21      117    113    106    106     89     98    116    108     93
    22      111    109    108    106     83     95    110    110     95
    23      109    106    100    101     90    104    106    102     94
    24      108     97    101    102     92     95    101     99     89
    25      103     71    100     94     82     87     92     90     79
    26       97     65     90     96     80     81     88     90     81
    27       97     64     88     95     71     72     87     89     62
    28       85     60     86     88     58     61     88     94     62
    29       86     59     85     85     56     57     85     91     60
    30       84     57     86     86     66     57     75     74     61
    31       96     68     86     87     73     64     84     81     60
    32      101     81     94     88     79     69     91     90     74
    33      106     87     94     96     79     79     95     91     75
    34      110     98    104     97     79     81     96     92     79
    35      110     92    101    105     83     95     93     97     88
 margin     151    147    149    153    144    154    151    149    141

 MINIMUM INTEGRATION VOLUME IS   56 GRID POINTS



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH CENTRAL PROFILE #1

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     1.00   0.89   0.86   0.87   0.82   0.90   0.89   0.87   0.82
     2     1.00   0.88   0.87   0.88   0.83   0.90   0.89   0.88   0.84
     3     1.00   0.88   0.88   0.89   0.85   0.90   0.90   0.90   0.85
     4     1.00   0.88   0.88   0.91   0.86   0.91   0.90   0.91   0.85
     5     1.00   0.87   0.87   0.92   0.87   0.89   0.90   0.91   0.86
     6     1.00   0.86   0.86   0.91   0.88   0.90   0.89   0.90   0.85
     7     1.00   0.88   0.86   0.92   0.88   0.90   0.89   0.91   0.87
     8     1.00   0.87   0.86   0.91   0.89   0.90   0.87   0.90   0.85
     9     1.00   0.90   0.85   0.90   0.88   0.88   0.87   0.89   0.87
    10     1.00   0.88   0.83   0.90   0.86   0.89   0.86   0.90   0.86
    11     1.00   0.87   0.83   0.91   0.87   0.90   0.84   0.89   0.84
    12     1.00   0.85   0.82   0.91   0.86   0.88   0.87   0.87   0.84
    13     1.00   0.86   0.86   0.91   0.87   0.87   0.86   0.88   0.82
    14     1.00   0.86   0.85   0.91   0.87   0.88   0.84   0.88   0.82
    15     1.00   0.83   0.87   0.91   0.88   0.88   0.87   0.87   0.79
    16     1.00   0.86   0.88   0.92   0.90   0.88   0.87   0.89   0.82
    17     1.00   0.87   0.88   0.92   0.90   0.89   0.87   0.90   0.83
    18     1.00   0.89   0.88   0.93   0.91   0.91   0.88   0.89   0.83
    19     1.00   0.90   0.89   0.94   0.91   0.91   0.91   0.90   0.85
    20     1.00   0.91   0.89   0.94   0.92   0.91   0.90   0.92   0.87
    21     1.00   0.91   0.91   0.93   0.92   0.92   0.89   0.90   0.87
    22     1.00   0.93   0.91   0.94   0.93   0.93   0.88   0.92   0.88
    23     1.00   0.93   0.93   0.94   0.94   0.94   0.90   0.91   0.89
    24     1.00   0.94   0.94   0.95   0.95   0.95   0.93   0.93   0.91
    25     1.00   0.93   0.95   0.96   0.95   0.96   0.93   0.95   0.92
    26     1.00   0.95   0.96   0.97   0.95   0.96   0.94   0.95   0.93
    27     1.00   0.96   0.98   0.98   0.95   0.96   0.95   0.97   0.95
    28     1.00   0.97   0.97   0.98   0.90   0.94   0.96   0.98   0.96
    29     1.00   0.96   0.97   0.97   0.92   0.94   0.96   0.97   0.95
    30     1.00   0.96   0.96   0.95   0.91   0.93   0.95   0.95   0.93
    31     1.00   0.96   0.94   0.93   0.91   0.92   0.92   0.93   0.94
    32     1.00   0.94   0.91   0.91   0.90   0.92   0.89   0.89   0.91
    33     1.00   0.93   0.89   0.90   0.89   0.92   0.87   0.87   0.92
    34     1.00   0.93   0.87   0.89   0.88   0.92   0.87   0.88   0.91
    35     1.00   0.92   0.87   0.88   0.87   0.92   0.87   0.88   0.91
 margin    1.00   0.96   0.96   0.97   0.95   0.95   0.96   0.97   0.95



 CORRELATIONS BETWEEN REFERENCE PROFILES WITH RUN-AVERAGED PROFILES

 BATCH-                           PROFILE NUMBER
 NUMBER      1      2      3      4      5      6      7      8      9
     1     0.97   0.94   0.90   0.88   0.86   0.94   0.91   0.87   0.85
     2     0.97   0.94   0.90   0.90   0.88   0.94   0.91   0.88   0.87
     3     0.96   0.93   0.91   0.91   0.91   0.94   0.91   0.89   0.89
     4     0.96   0.92   0.91   0.92   0.90   0.94   0.91   0.90   0.90
     5     0.96   0.91   0.90   0.93   0.91   0.93   0.91   0.91   0.90
     6     0.95   0.91   0.91   0.92   0.91   0.93   0.90   0.91   0.90
     7     0.95   0.91   0.91   0.92   0.90   0.92   0.91   0.92   0.92
     8     0.95   0.90   0.91   0.92   0.91   0.93   0.90   0.92   0.90
     9     0.95   0.91   0.90   0.91   0.91   0.91   0.90   0.93   0.92
    10     0.95   0.89   0.89   0.91   0.88   0.91   0.90   0.92   0.90
    11     0.94   0.90   0.90   0.90   0.87   0.91   0.90   0.91   0.88
    12     0.95   0.89   0.90   0.92   0.88   0.91   0.92   0.92   0.90
    13     0.96   0.90   0.92   0.91   0.90   0.92   0.92   0.93   0.90
    14     0.97   0.92   0.93   0.91   0.90   0.91   0.91   0.94   0.91
    15     0.97   0.91   0.93   0.90   0.89   0.91   0.92   0.93   0.90
    16     0.98   0.92   0.92   0.92   0.91   0.92   0.93   0.94   0.92
    17     0.98   0.94   0.92   0.93   0.91   0.93   0.93   0.95   0.93
    18     0.98   0.92   0.90   0.93   0.93   0.95   0.94   0.94   0.93
    19     0.98   0.92   0.90   0.94   0.95   0.97   0.95   0.95   0.94
    20     0.99   0.94   0.90   0.94   0.95   0.97   0.94   0.95   0.94
    21     0.98   0.95   0.88   0.89   0.94   0.97   0.92   0.93   0.96
    22     0.98   0.94   0.88   0.91   0.95   0.96   0.92   0.94   0.97
    23     0.97   0.95   0.86   0.89   0.93   0.95   0.90   0.92   0.96
    24     0.97   0.95   0.88   0.90   0.91   0.93   0.89   0.91   0.96
    25     0.95   0.93   0.88   0.89   0.91   0.91   0.86   0.89   0.94
    26     0.92   0.94   0.85   0.85   0.87   0.91   0.85   0.87   0.93
    27     0.92   0.93   0.85   0.86   0.92   0.92   0.87   0.89   0.93
    28     0.91   0.92   0.85   0.87   0.92   0.93   0.86   0.88   0.92
    29     0.90   0.90   0.84   0.86   0.91   0.93   0.85   0.85   0.91
    30     0.93   0.90   0.86   0.87   0.91   0.93   0.85   0.84   0.90
    31     0.95   0.92   0.86   0.86   0.90   0.93   0.86   0.87   0.91
    32     0.97   0.92   0.84   0.84   0.89   0.94   0.86   0.87   0.90
    33     0.98   0.93   0.83   0.85   0.89   0.95   0.86   0.86   0.90
    34     0.99   0.93   0.83   0.84   0.90   0.95   0.83   0.84   0.90
    35     0.98   0.94   0.85   0.84   0.89   0.94   0.82   0.83   0.89
 margin    1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00



      CORRELATION BETWEEN OBSERVED AND EXPECTED PROFILES

  I/Sigma  = mean of intensity/Sigma, where
  Sigma^2  =  4.0*[variance(I; from counting statistics)
                   + 0.00010*I^2]
  CORR     = mean correlation factor between observed and
             expected reflection profiles.
  E.S.D.   = estimated standard deviation of CORR.
  <I>      = mean LP-corrected reflection intensity,
             assuming unpolarized incident beam.
  NUMBER   = number of accepted reflections used to calculate
             I/Sigma, CORR, E.S.D., and <I>.

       I/Sigma        CORR      E.S.D.      <I>      NUMBER

       -3...-2      -0.115      0.021       -44.          5
       -2...-1      -0.082      0.015       -38.       1734
       -1... 0      -0.067      0.016       -10.      61926
        0... 1       0.073      0.018        13.     115136
        1... 2       0.104      0.023        38.      36085
        2... 3       0.152      0.030        63.      17280
        3... 6       0.238      0.056       107.      30486
        6... 9       0.350      0.065       178.      17080
        9...12       0.416      0.077       252.      10426
       12...15       0.438      0.088       323.       6131
       15...18       0.441      0.092       394.       3831
       18...21       0.428      0.093       454.       2495
       21...24       0.410      0.091       509.       1519
       24...27       0.384      0.090       565.       1017
       27...30       0.365      0.084       642.        700
       30...33       0.342      0.076       740.        508
       33...36       0.340      0.078       894.        353
       36...39       0.303      0.064       982.        240
       39...42       0.302      0.074      1311.        117
       42...45       0.279      0.068      1454.         63
       45...48       0.285      0.077      2443.         22
       48...51       0.112     -0.000      6061.          1

  307155 REFLECTIONS SAVED ON FILE "INTEGRATE.HKL"


 Time spent in the multi-tasking/processing part of INTEGRATE

     job #  cpu-time  wall-clock time (sec)
        1      321.7     332.3
        2      312.0     324.1
        3      311.2     321.7
        4      318.3     328.3
        5      318.9     330.7
        6      314.2     325.7
        7      312.1     313.7


 Time spent in the final, sequential part of INTEGRATE
            cpu-time  wall-clock time (sec)
                 2.8      14.3

 total cpu time used               2211.3 sec
 total elapsed wall-clock time      346.6 sec

