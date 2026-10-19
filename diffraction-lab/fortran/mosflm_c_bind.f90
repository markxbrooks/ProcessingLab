! Phase C: iso_c_binding mirror of include/mosflm/legacy/MosflmCApi.h
!
! Layout must stay in lockstep with the C header (padding-sensitive).
! Rebuild and run test_fortran_c_api after any MosflmCApi.h change.

module mosflm_c_bind
  use, intrinsic :: iso_c_binding
  implicit none
  private

  integer(c_int), parameter, public :: MOSFLM_C_MAX_SPOTS = 8192_c_int
  integer(c_int), parameter, public :: MOSFLM_SUCCESS = 0_c_int
  integer(c_int), parameter, public :: MOSFLM_ERROR_INVALID_ARGUMENT = 1_c_int
  integer(c_int), parameter, public :: MOSFLM_ERROR_BUFFER_TOO_SMALL = 2_c_int
  integer(c_int), parameter, public :: MOSFLM_ERROR_INTERNAL = 3_c_int

  type, bind(C), public :: mosflm_cell6
    real(c_double) :: a
    real(c_double) :: b
    real(c_double) :: c
    real(c_double) :: alpha_deg
    real(c_double) :: beta_deg
    real(c_double) :: gamma_deg
  end type mosflm_cell6

  type, bind(C), public :: mosflm_detector_c
    real(c_double) :: distance_mm
    real(c_double) :: beam_x_mm
    real(c_double) :: beam_y_mm
    real(c_double) :: pixel_size_x_mm
    real(c_double) :: pixel_size_y_mm
    integer(c_int) :: nfast
    integer(c_int) :: nslow
    real(c_double) :: wavelength_A
  end type mosflm_detector_c

  type, bind(C), public :: mosflm_ub9
    real(c_double) :: m(9)
  end type mosflm_ub9

  type, bind(C), public :: mosflm_predicted_c
    integer(c_int) :: h
    integer(c_int) :: k
    integer(c_int) :: l
    real(c_double) :: x_mm
    real(c_double) :: y_mm
    real(c_double) :: phi_begin_deg
    real(c_double) :: phi_end_deg
    integer(c_int) :: full
  end type mosflm_predicted_c

  type, bind(C), public :: mosflm_integrated_c
    integer(c_int) :: h
    integer(c_int) :: k
    integer(c_int) :: l
    real(c_double) :: x_mm
    real(c_double) :: y_mm
    real(c_double) :: intensity
    real(c_double) :: sigma
    integer(c_int) :: ok
  end type mosflm_integrated_c

  public :: mosflm_predict_still
  public :: mosflm_predict_osc
  public :: mosflm_ub_from_psi
  public :: mosflm_ub_from_matrix9
  public :: mosflm_analytical_still_fixture
  public :: mosflm_analytical_osc_fixture
  public :: mosflm_lp_factor
  public :: mosflm_integrate_summation
  public :: predict_still_cpp
  public :: predict_osc_cpp

  interface
    function mosflm_ub_from_psi(cell, psi_deg, out) result(rc) &
      bind(C, name='mosflm_ub_from_psi')
      import :: c_int, c_double, mosflm_cell6, mosflm_ub9
      type(mosflm_cell6), intent(in) :: cell
      real(c_double), intent(in) :: psi_deg(3)
      type(mosflm_ub9), intent(out) :: out
      integer(c_int) :: rc
    end function mosflm_ub_from_psi

    function mosflm_ub_from_matrix9(m, out) result(rc) &
      bind(C, name='mosflm_ub_from_matrix9')
      import :: c_int, c_double, mosflm_ub9
      real(c_double), intent(in) :: m(9)
      type(mosflm_ub9), intent(out) :: out
      integer(c_int) :: rc
    end function mosflm_ub_from_matrix9

    function mosflm_predict_still(cell, det, ub, spacegroup, phi_deg, &
                                  dmin_A, ewald_tol, out, out_cap, out_n) &
      result(rc) bind(C, name='mosflm_predict_still')
      import :: c_int, c_double, c_ptr, mosflm_cell6, mosflm_detector_c, &
                mosflm_predicted_c
      type(mosflm_cell6), intent(in) :: cell
      type(mosflm_detector_c), intent(in) :: det
      type(c_ptr), value :: ub
      type(c_ptr), value :: spacegroup
      real(c_double), value :: phi_deg
      real(c_double), value :: dmin_A
      real(c_double), value :: ewald_tol
      type(mosflm_predicted_c), intent(out) :: out(*)
      integer(c_int), value :: out_cap
      integer(c_int), intent(out) :: out_n
      integer(c_int) :: rc
    end function mosflm_predict_still

    function mosflm_predict_osc(cell, det, ub, spacegroup, phi_start, &
                                phi_end, dmin_A, ewald_tol, phi_step, out, &
                                out_cap, out_n) result(rc) &
      bind(C, name='mosflm_predict_osc')
      import :: c_int, c_double, c_ptr, mosflm_cell6, mosflm_detector_c, &
                mosflm_predicted_c
      type(mosflm_cell6), intent(in) :: cell
      type(mosflm_detector_c), intent(in) :: det
      type(c_ptr), value :: ub
      type(c_ptr), value :: spacegroup
      real(c_double), value :: phi_start
      real(c_double), value :: phi_end
      real(c_double), value :: dmin_A
      real(c_double), value :: ewald_tol
      real(c_double), value :: phi_step
      type(mosflm_predicted_c), intent(out) :: out(*)
      integer(c_int), value :: out_cap
      integer(c_int), intent(out) :: out_n
      integer(c_int) :: rc
    end function mosflm_predict_osc

    subroutine mosflm_analytical_still_fixture(cell, det, phi_deg, dmin_A, &
                                               ewald_tol) &
      bind(C, name='mosflm_analytical_still_fixture')
      import :: c_double, mosflm_cell6, mosflm_detector_c
      type(mosflm_cell6), intent(out) :: cell
      type(mosflm_detector_c), intent(out) :: det
      real(c_double), intent(out) :: phi_deg
      real(c_double), intent(out) :: dmin_A
      real(c_double), intent(out) :: ewald_tol
    end subroutine mosflm_analytical_still_fixture

    subroutine mosflm_analytical_osc_fixture(cell, det, phi_start, phi_end, &
                                             dmin_A, ewald_tol, phi_step) &
      bind(C, name='mosflm_analytical_osc_fixture')
      import :: c_double, mosflm_cell6, mosflm_detector_c
      type(mosflm_cell6), intent(out) :: cell
      type(mosflm_detector_c), intent(out) :: det
      real(c_double), intent(out) :: phi_start
      real(c_double), intent(out) :: phi_end
      real(c_double), intent(out) :: dmin_A
      real(c_double), intent(out) :: ewald_tol
      real(c_double), intent(out) :: phi_step
    end subroutine mosflm_analytical_osc_fixture

    function mosflm_lp_factor(det, x_mm, y_mm) result(lp) &
      bind(C, name='mosflm_lp_factor')
      import :: c_double, mosflm_detector_c
      type(mosflm_detector_c), intent(in) :: det
      real(c_double), value :: x_mm
      real(c_double), value :: y_mm
      real(c_double) :: lp
    end function mosflm_lp_factor

    function mosflm_integrate_summation(det, pixels, nfast, nslow, preds, &
                                        n_preds, box_half_px, bg_frame_px, &
                                        apply_lp, out, out_cap, out_n) &
      result(rc) bind(C, name='mosflm_integrate_summation')
      import :: c_int, mosflm_detector_c, mosflm_predicted_c, &
                mosflm_integrated_c
      type(mosflm_detector_c), intent(in) :: det
      integer(c_int), intent(in) :: pixels(*)
      integer(c_int), value :: nfast
      integer(c_int), value :: nslow
      type(mosflm_predicted_c), intent(in) :: preds(*)
      integer(c_int), value :: n_preds
      integer(c_int), value :: box_half_px
      integer(c_int), value :: bg_frame_px
      integer(c_int), value :: apply_lp
      type(mosflm_integrated_c), intent(out) :: out(*)
      integer(c_int), value :: out_cap
      integer(c_int), intent(out) :: out_n
      integer(c_int) :: rc
    end function mosflm_integrate_summation
  end interface

contains

  !> Convenience stills PREDICT (optional analytical fixture fill).
  subroutine predict_still_cpp(use_fixture, cell, det, phi_deg, dmin_A, &
                               ewald_tol, spots, nspots, rc)
    logical, intent(in) :: use_fixture
    type(mosflm_cell6), intent(inout) :: cell
    type(mosflm_detector_c), intent(inout) :: det
    real(c_double), intent(inout) :: phi_deg
    real(c_double), intent(inout) :: dmin_A
    real(c_double), intent(inout) :: ewald_tol
    type(mosflm_predicted_c), intent(out) :: spots(MOSFLM_C_MAX_SPOTS)
    integer(c_int), intent(out) :: nspots
    integer(c_int), intent(out) :: rc

    if (use_fixture) then
      call mosflm_analytical_still_fixture(cell, det, phi_deg, dmin_A, &
                                           ewald_tol)
    end if
    rc = mosflm_predict_still(cell, det, c_null_ptr, c_null_ptr, phi_deg, &
                              dmin_A, ewald_tol, spots, MOSFLM_C_MAX_SPOTS, &
                              nspots)
  end subroutine predict_still_cpp

  !> Convenience oscillation PREDICT (optional analytical fixture fill).
  subroutine predict_osc_cpp(use_fixture, cell, det, phi_start, phi_end, &
                             dmin_A, ewald_tol, phi_step, spots, nspots, rc)
    logical, intent(in) :: use_fixture
    type(mosflm_cell6), intent(inout) :: cell
    type(mosflm_detector_c), intent(inout) :: det
    real(c_double), intent(inout) :: phi_start
    real(c_double), intent(inout) :: phi_end
    real(c_double), intent(inout) :: dmin_A
    real(c_double), intent(inout) :: ewald_tol
    real(c_double), intent(inout) :: phi_step
    type(mosflm_predicted_c), intent(out) :: spots(MOSFLM_C_MAX_SPOTS)
    integer(c_int), intent(out) :: nspots
    integer(c_int), intent(out) :: rc

    if (use_fixture) then
      call mosflm_analytical_osc_fixture(cell, det, phi_start, phi_end, &
                                         dmin_A, ewald_tol, phi_step)
    end if
    rc = mosflm_predict_osc(cell, det, c_null_ptr, c_null_ptr, phi_start, &
                            phi_end, dmin_A, ewald_tol, phi_step, spots, &
                            MOSFLM_C_MAX_SPOTS, nspots)
  end subroutine predict_osc_cpp

end module mosflm_c_bind
