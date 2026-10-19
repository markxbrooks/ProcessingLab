! Fill amalgamation COMMONs, call PREDICT_CPP_FROM_COMMON (AMAT + osc).

program test_predict_cpp_from_common
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int), parameter :: nmax = 8192_c_int

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

  ! /REEKE/ (minimal fill; AMAT + DSTMAX matter)
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

  integer(c_int) :: nspot, iflagpr
  integer(c_short) :: ihg(nmax), ikg(nmax), ilg(nmax), img(nmax)
  real(c_float) :: xg(nmax), yg(nmax), phig(nmax), phiwg(nmax)
  real(c_float) :: a_star, c_star

  external :: predict_cpp_from_common

  ! Zero commons then fill analytical stills fixture geometry.
  cell = [58.0_c_float, 58.0_c_float, 150.0_c_float, 90.0_c_float, &
          90.0_c_float, 90.0_c_float]
  xtofd = 25000.0_c_float
  rast = 0.172_c_float
  rasty = 0.172_c_float
  iylen = 1475_c_int
  nrec = 1679_c_int
  wave = 1.0_c_float
  phibeg = 0.0_c_float
  phiend = 0.0_c_float
  dstmax = wave / 3.5_c_float

  ! AMAT such that mosflm_lab_to_cpp * AMAT = diag(1/a,1/b,1/c) in C++ frame.
  ! R is an involution, so AMAT = R * A*_cpp.
  a_star = 1.0_c_float / 58.0_c_float
  c_star = 1.0_c_float / 150.0_c_float
  amat = 0.0_c_float
  amat(1, 3) = c_star
  amat(2, 2) = -a_star
  amat(3, 1) = a_star

  call predict_cpp_from_common(iflagpr, nmax, ihg, ikg, ilg, xg, yg, phig, &
                               phiwg, img, nspot)

  if (iflagpr /= 0 .or. nspot <= 0) then
    write (*, '(A,I0,A,I0)') 'FAIL: iflagpr=', iflagpr, ' nspot=', nspot
    stop 1
  end if

  write (*, '(A,I0)') 'test_predict_cpp_from_common: OK n=', nspot
  stop 0
end program test_predict_cpp_from_common
