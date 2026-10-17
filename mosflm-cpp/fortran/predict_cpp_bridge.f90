! Amalgamation-oriented stills PREDICT bridge (Phase C).
!
! Name ``predict_cpp_still`` is the intended drop-in near
! ``CALL PREDICT(PSI, IFLAGPR)`` when building with MOSFLM_CPP_PREDICT.
! It forwards to the C ABI ``predict_still_cpp_`` (gfortran underscore).
!
! See docs/PHASE_C_SPLICE.md for the ifdef recipe and COMMON mapping notes.

subroutine predict_cpp_still(cell, distance_mm, wavelength_A, pixel_x_mm, &
                             pixel_y_mm, nfast, nslow, beam_x_mm, beam_y_mm, &
                             phi_deg, dmin_A, ewald_tol, h, k, l, x_mm, y_mm, &
                             phi_beg, phi_end, ifull, nmax, nspot, ierr)
  use, intrinsic :: iso_c_binding
  implicit none

  real(c_float), intent(inout) :: cell(6)
  real(c_double), intent(inout) :: distance_mm
  real(c_double), intent(inout) :: wavelength_A
  real(c_double), intent(inout) :: pixel_x_mm
  real(c_double), intent(inout) :: pixel_y_mm
  integer(c_int), intent(inout) :: nfast
  integer(c_int), intent(inout) :: nslow
  real(c_double), intent(inout) :: beam_x_mm
  real(c_double), intent(inout) :: beam_y_mm
  real(c_double), intent(inout) :: phi_deg
  real(c_double), intent(inout) :: dmin_A
  real(c_double), intent(inout) :: ewald_tol
  integer(c_int), intent(out) :: h(*)
  integer(c_int), intent(out) :: k(*)
  integer(c_int), intent(out) :: l(*)
  real(c_double), intent(out) :: x_mm(*)
  real(c_double), intent(out) :: y_mm(*)
  real(c_double), intent(out) :: phi_beg(*)
  real(c_double), intent(out) :: phi_end(*)
  integer(c_int), intent(out) :: ifull(*)
  integer(c_int), intent(in) :: nmax
  integer(c_int), intent(out) :: nspot
  integer(c_int), intent(out) :: ierr

  interface
    subroutine predict_still_cpp(cell, distance_mm, wavelength_A, pixel_x_mm, &
                                 pixel_y_mm, nfast, nslow, beam_x_mm, &
                                 beam_y_mm, phi_deg, dmin_A, ewald_tol, h, k, &
                                 l, x_mm, y_mm, phi_beg, phi_end, ifull, nmax, &
                                 nspot, ierr) &
      bind(C, name='predict_still_cpp_')
      import :: c_float, c_double, c_int
      real(c_float), intent(inout) :: cell(6)
      real(c_double), intent(inout) :: distance_mm
      real(c_double), intent(inout) :: wavelength_A
      real(c_double), intent(inout) :: pixel_x_mm
      real(c_double), intent(inout) :: pixel_y_mm
      integer(c_int), intent(inout) :: nfast
      integer(c_int), intent(inout) :: nslow
      real(c_double), intent(inout) :: beam_x_mm
      real(c_double), intent(inout) :: beam_y_mm
      real(c_double), intent(inout) :: phi_deg
      real(c_double), intent(inout) :: dmin_A
      real(c_double), intent(inout) :: ewald_tol
      integer(c_int), intent(out) :: h(*)
      integer(c_int), intent(out) :: k(*)
      integer(c_int), intent(out) :: l(*)
      real(c_double), intent(out) :: x_mm(*)
      real(c_double), intent(out) :: y_mm(*)
      real(c_double), intent(out) :: phi_beg(*)
      real(c_double), intent(out) :: phi_end(*)
      integer(c_int), intent(out) :: ifull(*)
      integer(c_int), intent(in) :: nmax
      integer(c_int), intent(out) :: nspot
      integer(c_int), intent(out) :: ierr
    end subroutine predict_still_cpp
  end interface

  call predict_still_cpp(cell, distance_mm, wavelength_A, pixel_x_mm, &
                         pixel_y_mm, nfast, nslow, beam_x_mm, beam_y_mm, &
                         phi_deg, dmin_A, ewald_tol, h, k, l, x_mm, y_mm, &
                         phi_beg, phi_end, ifull, nmax, nspot, ierr)
end subroutine predict_cpp_still

! Stills PREDICT with missetting PSI(3) degrees → mosflm_ub_from_psi.
subroutine predict_cpp_still_psi(cell, distance_mm, wavelength_A, &
                                 pixel_x_mm, pixel_y_mm, nfast, nslow, &
                                 beam_x_mm, beam_y_mm, phi_deg, dmin_A, &
                                 ewald_tol, psi_deg, h, k, l, x_mm, y_mm, &
                                 phi_beg, phi_end, ifull, nmax, nspot, ierr)
  use, intrinsic :: iso_c_binding
  implicit none

  real(c_float), intent(inout) :: cell(6)
  real(c_double), intent(inout) :: distance_mm
  real(c_double), intent(inout) :: wavelength_A
  real(c_double), intent(inout) :: pixel_x_mm
  real(c_double), intent(inout) :: pixel_y_mm
  integer(c_int), intent(inout) :: nfast
  integer(c_int), intent(inout) :: nslow
  real(c_double), intent(inout) :: beam_x_mm
  real(c_double), intent(inout) :: beam_y_mm
  real(c_double), intent(inout) :: phi_deg
  real(c_double), intent(inout) :: dmin_A
  real(c_double), intent(inout) :: ewald_tol
  real(c_double), intent(inout) :: psi_deg(3)
  integer(c_int), intent(out) :: h(*)
  integer(c_int), intent(out) :: k(*)
  integer(c_int), intent(out) :: l(*)
  real(c_double), intent(out) :: x_mm(*)
  real(c_double), intent(out) :: y_mm(*)
  real(c_double), intent(out) :: phi_beg(*)
  real(c_double), intent(out) :: phi_end(*)
  integer(c_int), intent(out) :: ifull(*)
  integer(c_int), intent(in) :: nmax
  integer(c_int), intent(out) :: nspot
  integer(c_int), intent(out) :: ierr

  interface
    subroutine predict_still_psi_cpp(cell, distance_mm, wavelength_A, &
                                     pixel_x_mm, pixel_y_mm, nfast, nslow, &
                                     beam_x_mm, beam_y_mm, phi_deg, dmin_A, &
                                     ewald_tol, psi_deg, h, k, l, x_mm, &
                                     y_mm, phi_beg, phi_end, ifull, nmax, &
                                     nspot, ierr) &
      bind(C, name='predict_still_psi_cpp_')
      import :: c_float, c_double, c_int
      real(c_float), intent(inout) :: cell(6)
      real(c_double), intent(inout) :: distance_mm
      real(c_double), intent(inout) :: wavelength_A
      real(c_double), intent(inout) :: pixel_x_mm
      real(c_double), intent(inout) :: pixel_y_mm
      integer(c_int), intent(inout) :: nfast
      integer(c_int), intent(inout) :: nslow
      real(c_double), intent(inout) :: beam_x_mm
      real(c_double), intent(inout) :: beam_y_mm
      real(c_double), intent(inout) :: phi_deg
      real(c_double), intent(inout) :: dmin_A
      real(c_double), intent(inout) :: ewald_tol
      real(c_double), intent(inout) :: psi_deg(3)
      integer(c_int), intent(out) :: h(*)
      integer(c_int), intent(out) :: k(*)
      integer(c_int), intent(out) :: l(*)
      real(c_double), intent(out) :: x_mm(*)
      real(c_double), intent(out) :: y_mm(*)
      real(c_double), intent(out) :: phi_beg(*)
      real(c_double), intent(out) :: phi_end(*)
      integer(c_int), intent(out) :: ifull(*)
      integer(c_int), intent(in) :: nmax
      integer(c_int), intent(out) :: nspot
      integer(c_int), intent(out) :: ierr
    end subroutine predict_still_psi_cpp
  end interface

  call predict_still_psi_cpp(cell, distance_mm, wavelength_A, pixel_x_mm, &
                             pixel_y_mm, nfast, nslow, beam_x_mm, beam_y_mm, &
                             phi_deg, dmin_A, ewald_tol, psi_deg, h, k, l, &
                             x_mm, y_mm, phi_beg, phi_end, ifull, nmax, &
                             nspot, ierr)
end subroutine predict_cpp_still_psi
