! Fill /CELLCOM/ + /XY/ fixture values, call PREDICT_CPP_FROM_COMMON.

program test_predict_cpp_from_common
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int), parameter :: nmax = 8192_c_int
  integer(c_int), parameter :: expect = 656_c_int

  real(c_float) :: umat(3, 3), bmat(3, 3), gmat(3, 3)
  real(c_float) :: cell(6), rcell(6), umatcell(6), softcell, volume
  integer(c_int) :: lcell(6), icryst, numspg, nlaue
  common /cellcom/ umat, bmat, gmat, cell, rcell, umatcell, softcell, &
                   volume, lcell, icryst, numspg, nlaue

  real(c_float) :: xtofd, sinv, cosv, tanv, twotheta
  integer(c_int) :: icass
  common /xy/ xtofd, sinv, cosv, tanv, twotheta, icass

  real(c_float) :: psi(3), wave, rast, rasty, phi, dmin, etol
  integer(c_int) :: nf, ns, nspot, iflagpr
  integer(c_short) :: ihg(nmax), ikg(nmax), ilg(nmax), img(nmax)
  real(c_float) :: xg(nmax), yg(nmax), phig(nmax), phiwg(nmax)

  external :: predict_cpp_from_common

  umat = 0.0_c_float
  bmat = 0.0_c_float
  gmat = 0.0_c_float
  cell = [58.0_c_float, 58.0_c_float, 150.0_c_float, 90.0_c_float, &
          90.0_c_float, 90.0_c_float]
  rcell = 0.0_c_float
  umatcell = 0.0_c_float
  softcell = 0.0_c_float
  volume = 0.0_c_float
  lcell = 0_c_int
  icryst = 0_c_int
  numspg = 0_c_int
  nlaue = 0_c_int

  xtofd = 25000.0_c_float
  sinv = 0.0_c_float
  cosv = 0.0_c_float
  tanv = 0.0_c_float
  twotheta = 0.0_c_float
  icass = 0_c_int

  wave = 1.0_c_float
  rast = 0.172_c_float
  rasty = 0.172_c_float
  nf = 1475_c_int
  ns = 1679_c_int
  phi = 0.0_c_float
  dmin = 3.5_c_float
  etol = 0.0025_c_float
  psi = [0.0_c_float, 0.0_c_float, 0.0_c_float]

  call predict_cpp_from_common(psi, iflagpr, wave, rast, rasty, nf, ns, &
                               phi, dmin, etol, nmax, ihg, ikg, ilg, xg, &
                               yg, phig, phiwg, img, nspot)

  if (iflagpr /= 0 .or. nspot /= expect) then
    write (*, '(A,I0,A,I0)') 'FAIL: iflagpr=', iflagpr, ' nspot=', nspot
    stop 1
  end if

  write (*, '(A,I0)') 'test_predict_cpp_from_common: OK n=', nspot
  stop 0
end program test_predict_cpp_from_common
