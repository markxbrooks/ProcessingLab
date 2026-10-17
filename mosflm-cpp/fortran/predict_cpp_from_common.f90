! Read CELL / XTOFD from amalgamation COMMON, then call PREDICT_CPP.
!
! REFLFIND does not include /CELLCOM/ or /XY/; this helper supplies them.
! Caller passes geometry already in scope (WAVE, RAST, GENDATA arrays, …).
!
! See docs/PHASE_C_SPLICE.md.

subroutine predict_cpp_from_common(psi, iflagpr, wave, rast, rasty, &
                                   nfast, nslow, phi, dmin_A, ewald_tol, &
                                   nmax, ihg, ikg, ilg, xg, yg, phig, &
                                   phiwg, img, nspot)
  use, intrinsic :: iso_c_binding
  implicit none

  real(c_float), intent(inout) :: psi(3)
  integer(c_int), intent(out) :: iflagpr
  real(c_float), intent(in) :: wave
  real(c_float), intent(in) :: rast
  real(c_float), intent(in) :: rasty
  integer(c_int), intent(in) :: nfast
  integer(c_int), intent(in) :: nslow
  real(c_float), intent(in) :: phi
  real(c_float), intent(in) :: dmin_A
  real(c_float), intent(in) :: ewald_tol
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

  ! /CELLCOM/ layout must match mosflm/inc/cell.f
  real(c_float) :: umat(3, 3), bmat(3, 3), gmat(3, 3)
  real(c_float) :: cell(6), rcell(6), umatcell(6), softcell, volume
  integer(c_int) :: lcell(6), icryst, numspg, nlaue
  common /cellcom/ umat, bmat, gmat, cell, rcell, umatcell, softcell, &
                   volume, lcell, icryst, numspg, nlaue

  ! /XY/ layout must match mosflm/inc/xy.f
  real(c_float) :: xtofd, sinv, cosv, tanv, twotheta
  integer(c_int) :: icass
  common /xy/ xtofd, sinv, cosv, tanv, twotheta, icass

  external :: predict_cpp

  call predict_cpp(psi, iflagpr, cell, xtofd, wave, rast, rasty, nfast, &
                   nslow, phi, dmin_A, ewald_tol, nmax, ihg, ikg, ilg, &
                   xg, yg, phig, phiwg, img, nspot)
end subroutine predict_cpp_from_common
