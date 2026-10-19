! Amalgamation-facing PREDICT_CPP hook (Phase C).
!
! Drop-in near ``CALL PREDICT`` / ``CALL REEK`` when ``MOSFLM_CPP_PREDICT`` is
! defined. Uses MOSFLM AMAT (lab-frame mapped) + oscillation φ range.
! Writes GENDATA-shaped ``INTEGER*2`` miller / ``IMG`` arrays and ``NSPOT``.
!
! See docs/PHASE_C_SPLICE.md.

subroutine predict_cpp(psi, iflagpr, cell, xtofd, wave, rast, rasty, &
                       nfast, nslow, phi_beg, phi_end, dmin_A, ewald_tol, &
                       amat, nmax, ihg, ikg, ilg, xg, yg, phig, phiwg, &
                       img, nspot)
  use, intrinsic :: iso_c_binding
  use predict_gendata_pack
  implicit none

  real(c_float), intent(inout) :: psi(3)
  integer(c_int), intent(out) :: iflagpr
  real(c_float), intent(inout) :: cell(6)
  real(c_float), intent(in) :: xtofd          ! 10 µm units (/XY/)
  real(c_float), intent(in) :: wave           ! Å
  real(c_float), intent(in) :: rast           ! pixel slow, mm
  real(c_float), intent(in) :: rasty          ! pixel fast, mm (RASTY)
  integer(c_int), intent(in) :: nfast
  integer(c_int), intent(in) :: nslow
  real(c_float), intent(in) :: phi_beg        ! oscillation start, degrees
  real(c_float), intent(in) :: phi_end        ! oscillation end, degrees
  real(c_float), intent(in) :: dmin_A
  real(c_float), intent(in) :: ewald_tol
  real(c_float), intent(in) :: amat(3, 3)     ! MOSFLM Orientation Matrix [A]
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

  real(c_double) :: distance_mm, wavelength_A, pixel_x_mm, pixel_y_mm
  real(c_double) :: beam_x_mm, beam_y_mm, phi0, phi1, dmin, etol, step
  integer(c_int) :: nf, ns, ierr
  real(c_float) :: cell_work(6), amat_work(3, 3)

  iflagpr = 0_c_int
  nspot = 0_c_int
  ! psi retained for call-site ABI; orientation comes from amat.
  psi(1) = psi(1)

  cell_work = cell
  amat_work = amat
  distance_mm = mosflm_xtofd_to_mm(xtofd)
  wavelength_A = real(wave, kind=c_double)
  ! C++ pixel_x = fast, pixel_y = slow; MOSFLM RAST=slow, RASTY=fast.
  pixel_x_mm = real(rasty, kind=c_double)
  pixel_y_mm = real(rast, kind=c_double)
  nf = nfast
  ns = nslow
  beam_x_mm = 0.0_c_double
  beam_y_mm = 0.0_c_double
  phi0 = real(phi_beg, kind=c_double)
  phi1 = real(phi_end, kind=c_double)
  dmin = real(dmin_A, kind=c_double)
  etol = real(ewald_tol, kind=c_double)
  step = 0.05_c_double

  call predict_cpp_osc_amat_to_gendata_i2(cell_work, distance_mm, &
                                          wavelength_A, pixel_x_mm, &
                                          pixel_y_mm, nf, ns, beam_x_mm, &
                                          beam_y_mm, phi0, phi1, dmin, etol, &
                                          step, amat_work, nmax, ihg, ikg, &
                                          ilg, xg, yg, phig, phiwg, img, &
                                          nspot, iflagpr, ierr)
  if (ierr /= 0_c_int .and. iflagpr == 0_c_int) iflagpr = -1_c_int
end subroutine predict_cpp
