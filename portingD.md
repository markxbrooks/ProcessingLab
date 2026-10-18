Phase D: Retire Fortran (program plan)

What complete means

Per porting.md Phase D: retire Fortran only when all callers migrated.

Definition of Done:





A C++-owned path runs known-indexed insulin (geometry, predict, measure, integrate, MTZ) without linking mosflm_all_ip_inc_noX.o / REEK / Fortran PROCESS.



Differentials vs [mosflm-cpp/testdata/oracle/known_indexed/](mosflm-cpp/testdata/oracle/known_indexed/) meet documented HKL / XDET / intensity gates.



Default release artifact does not include the amalgamation; Autotools ipmosflm remains an optional oracle build.



porting.md marks Phase D Done.

This is a multi-milestone program (months). After approval, implement D0 only, then review before D1.

flowchart LR
  D0[D0_Close_Phase_C] --> D1[D1_Cpp_PROCESS]
  D1 --> D2[D2_MEAS]
  D2 --> D3[D3_INTEG3]
  D3 --> D4[D4_PostRefinement]
  D4 --> D5[D5_CONTROL]
  D5 --> D6[D6_Drop_amalgamation]

Baseline already landed





Vertical slice under [mosflm-cpp/](mosflm-cpp/): geometry, miniCBF, predict, summation, Gaussian profile, MTZ.



Lab frame / XDET: [LabFrame.hpp](mosflm-cpp/include/mosflm/geometry/LabFrame.hpp); CTests insulin_hkl_overlap, insulin_xdet.



Phase C splice: MOSFLM_CPP_PREDICT=1 on GENERATE/REFLFIND; insulin PROCESS ~89% exact-HKL recall ([known_indexed_cpp_predict/](mosflm-cpp/testdata/oracle/known_indexed_cpp_predict/), [PHASE_C_SPLICE.md](mosflm-cpp/docs/PHASE_C_SPLICE.md)).



Still Fortran-owned: CONTROL, PROCESS, REEK, MEAS, full INTEG3, post-refinement, SMV I/O.

Essential Fortran OBJS today ([mosflm/Makefile](mosflm/Makefile)): control.o, unoptimizable.o, mosflm_all_ip_inc_noX.o, torefix_noX.o.

Fortran call chain (inventory)

flowchart LR
  CONTROL --> GENERATE
  GENERATE --> REEK
  GENERATE --> MEAS
  MEAS --> PROCESS
  PROCESS --> INTEG
  PROCESS --> POSTREFL







Capability



Location



Role





CONTROL



mosflm/control.f:4



Keyword loop; PROCESS/GO/POSTREF





PROGRAM MOSFLM



unoptimizable.f:5



Binary entry





GENERATE



amalgamation ~66626



Per-image predict + genfile





REEK



unoptimizable.f:10656



Osc prediction (mosaicity/divergence/ice)





PREDICT / REFLFIND



amalgamation ~171811 / ~194120



Interactive; Phase C can splice C++





MEAS



amalgamation ~119126



Stripe to shoebox; feeds PROCESS





PROCESS



unoptimizable.f:25854



Profiles, integrate, write intensities





INTEG / INTEG2



unoptimizable.f:22256 / 20799



Summation / profile-path integrate





INTEG3 / INTEG4



amalgamation ~97463 / ~98162



Full profile fitting





POSTREFL



amalgamation ~168264



Post-refinement reflection list

I/O/runtime OBJS (CBF, HDF5, sockets, JPEG) are not diffraction kernels; replace or keep as C later—science amalgamation is the Phase D blocker.



D0 — Close Phase C (1–2 weeks)





CTest comparing known_indexed_cpp_predict vs known_indexed MTZ (exact-HKL recall gate >= 0.85).



Predictor mosaicity / Ewald width so spot counts approach REEK (η, divergence still the main gap vs Fortran REEK).



Detector TILT/TWIST/YSCALE; tighten insulin_xdet toward ~0.25 mm refine residual.



Sync porting.md / README (Phase C Done; step 5 MTZ todo closed).

Exit: Phase C Done; known-indexed CTests green with published gates.

D1 — C++ PROCESS + GENERATE path (2–4 weeks)

Own the known-indexed pipeline without Fortran GENERATE/REEK/PROCESS orchestration:





ProcessRequest / ProcessResult in [DOMAIN_API.md](mosflm-cpp/docs/DOMAIN_API.md).



SMV image reader (insulin ADSC) under mosflm-cpp/src/image/.



In-memory spot table (genfile optional later); predict -> summation integrate -> MTZ.



CLI: geometry from [geometry.txt](mosflm-cpp/testdata/images/known_indexed/geometry.txt).



CTest without amalgamation / mosflm_predict_hook; reuse D0 HKL gate.



Optional transitional: splice C++ INTEG/INTEG2 behind PROCESS via ABI while D1 CLI lands (does not replace D1 exit).

Exit: Pure C++ insulin PROCESS CTest green.

D2 — MEAS measurement boxes (3–6 weeks)





Port measurement-box / background / overlap contract from Fortran MEAS (~amalgamation 119126).



Soft intensity/sigma differential vs oracle I/SIGI on matched HKL.

Exit: C++ PROCESS path no longer needs Fortran MEAS for insulin.

D3 — Full INTEG3 parity (1–2 months)





Replace Gaussian ProfileFitter with standard-profile behavior from INTEG3/INTEG4 (~97463 / ~98162).



Differentials on profile-related MTZ columns / flags.

Exit: porting.md INTEG3 row = parity, not stand-in.

D4 — Post-refinement Stage 2 (1–2 months)





New mosflm-cpp/src/refinement/ (cell, missetting, mosaic from partials); port POSTREFL-class behavior (~168264).



Capture/extend oracles; differentials before enabling by default.

Exit: Stage 2 checklist Done in porting.md.

D5 — CONTROL keyword subset





Batch keywords for known-indexed script: CELL, SYMMETRY, BEAM, DISTANCE, WAVELENGTH, PIXEL, MATRIX/AMAT, PROCESS, GO, EXIT.



Run equivalent of [ipmosflm.in](mosflm-cpp/testdata/images/known_indexed/ipmosflm.in) on Cdriver (thin CONTROL or pure C parser).

Exit: Keyword-driven C++ PROCESS without Fortran CONTROL for this corpus.

D6 — Drop amalgamation (Phase D complete)





Default build/install is C++ core (+ thin ipmosflm name for CCP4 if needed).



Amalgamation behind MOSFLM_FORTRAN_ORACLE=1 (or similar), CI optional/nightly.



Default binary links without mosflm_all_ip_inc_noX.o / science unoptimizable.o.



Mark Phase D Done in porting.md.

Out of Phase D DoD





iMosflm/Tcl sockets, full STRATEGY/multi-lattice, every detector, REFIX (already stubbed).

Execution rule

Start D0 only after plan approval; do not begin D1 until D0 exit criteria pass.