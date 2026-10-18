! Read geometry from amalgamation COMMON, then call PREDICT_CPP.
!
! Caller only supplies GENDATA / IOSP buffers. CELL, XTOFD, AMAT, WAVE,
! RAST*, detector size, PHIBEG/PHIEND, and DSTMAX are read here.
! COMMON layouts must match mosflm/inc/{cell,xy,scn,misc,reeke}.f.
! See docs/PHASE_C_SPLICE.md.

subroutine predict_cpp_from_common(iflagpr, nmax, ihg, ikg, ilg, xg, yg, &
                                   phig, phiwg, img, nspot)
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int), intent(out) :: iflagpr
  integer(c_int), intent(in) :: nmax
  integer(c_short), intent(out) :: ihg(*)
  integer(c_short), intent(out) :: ikg(*)
  integer(c_short), intent(out) :: ilg(*)
  real(c_float), intent(out) :: xg(*)
  real(c_float), intent(out) :: yg(*)
  real(c_float), intent(out) :: phig(*)
  real(c_float), intent(out) :: phiwg(*)
  integer(c_short), intent(out) :: img(*)
  integer(c_int), intent(out) :: nspot

  ! /CELLCOM/
  real(c_float) :: umat(3, 3), bmat(3, 3), gmat(3, 3)
  real(c_float) :: cell(6), rcell(6), umatcell(6), softcell, volume
  integer(c_int) :: lcell(6), icryst, numspg, nlaue
  common /cellcom/ umat, bmat, gmat, cell, rcell, umatcell, softcell, &
                   volume, lcell, icryst, numspg, nlaue

  ! /XY/
  real(c_float) :: xtofd, sinv, cosv, tanv, twotheta
  integer(c_int) :: icass
  common /xy/ xtofd, sinv, cosv, tanv, twotheta, icass

  ! /SCN/
  real(c_float) :: fact, scnsz, rast, rasty, rowreadt, rotspeed
  integer(c_int) :: nbyte, nrec, nword, iylen, icurr, nextra, nhbyte
  logical :: bytswap, lactive, reallyforceread
  common /scn/ fact, scnsz, rast, rasty, rowreadt, rotspeed, nbyte, nrec, &
               nword, iylen, icurr, nextra, nhbyte, bytswap, lactive, &
               reallyforceread

  ! /MISC/
  real(c_float) :: delphi(3), resani(3)
  real(c_double) :: starttime, endtime
  real(c_float) :: phibeg, phiend, phirng, rmin, rmax, wave
  integer(c_int) :: iax(3), ipackid, minint, ierrflg
  logical :: anites
  common /misc/ delphi, resani, starttime, endtime, phibeg, phiend, &
                phirng, rmin, rmax, wave, iax, ipackid, minint, ierrflg, &
                anites

  ! /REEKE/ — storage order matches amalgamation COMMON statement.
  real(c_float) :: rmc(3, 3), amat(3, 3), scaxis(3), s0(3), s0unit(3)
  real(c_float) :: camcon(8), beammu(2), x1, x2, y1, y2, z1, z2, xys
  real(c_float) :: dstar2, dstmax, dstpl, dstpl2, divh, divv, delamb
  real(c_float) :: eta, delcor, dstmaxs, tor, dstmin, wmax, etafrac, etasm
  real(c_float) :: gon_chi, gon_phis, gon_kappa, gon_omegas, gon_2theta
  real(c_float) :: gon_phie, gon_omegae, gon_chi_set, gon_kappa_set
  real(c_float) :: blksize
  integer(c_int) :: isyn, imono, nwmax, ipad, nsmooth
  integer(c_int) :: igon_chi, igon_phi, igon_kappa, igon_omega, igon_2theta
  integer(c_int) :: igon_chi_set, igon_kappa_set, ihoscaxis
  logical :: logeta, nureek, lcamcon, gengon, fourcircle, threecircle
  logical :: vertical
  common /reeke/ rmc, amat, scaxis, s0, s0unit, camcon, beammu, x1, x2, &
                 y1, y2, z1, z2, xys, dstar2, dstmax, dstpl, dstpl2, divh, &
                 divv, delamb, eta, delcor, dstmaxs, tor, dstmin, wmax, &
                 etafrac, etasm, gon_chi, gon_phis, gon_kappa, gon_omegas, &
                 gon_2theta, gon_phie, gon_omegae, gon_chi_set, &
                 gon_kappa_set, isyn, imono, nwmax, ipad, nsmooth, &
                 igon_chi, igon_phi, igon_kappa, igon_omega, igon_2theta, &
                 igon_chi_set, igon_kappa_set, ihoscaxis, logeta, nureek, &
                 lcamcon, gengon, fourcircle, threecircle, blksize, &
                 vertical

  real(c_float) :: psi(3), dmin_a, etol
  integer(c_int) :: nfast, nslow

  external :: predict_cpp

  psi = 0.0_c_float
  dmin_a = 0.0_c_float
  if (dstmax > 1.0e-12_c_float) dmin_a = wave / dstmax
  etol = 0.0025_c_float
  nfast = iylen
  nslow = nrec

  call predict_cpp(psi, iflagpr, cell, xtofd, wave, rast, rasty, nfast, &
                   nslow, phibeg, phiend, dmin_a, etol, amat, nmax, ihg, &
                   ikg, ilg, xg, yg, phig, phiwg, img, nspot)
end subroutine predict_cpp_from_common
