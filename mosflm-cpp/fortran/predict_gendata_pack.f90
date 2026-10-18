! Phase C: pack/unpack between C++ stills PREDICT and MOSFLM /GENDATA/ fields.
!
! Offline helper (no amalgamation COMMON linkage). Arrays are shaped like
! gendata.f but sized by the caller (``nmax``). See docs/PHASE_C_SPLICE.md.
!
! Unit conventions (from mosflm/inc/gendata.f comments):
!   XG, YG  — virtual detector coords in **10 micron units**, beam-relative
!             (mm * 100). C++ bridge uses millimetres.
!   Axis map (insulin XDET differential / docs/PHASE_C_SPLICE.md):
!             XG ↔ C++ y_mm,  YG ↔ C++ x_mm.
!   IMG     — 0 = full; non-zero = partial (1..100 in legacy REEK).
!   PHIG    — phi begin (degrees); PHIWG — phi width (end - begin).

module predict_gendata_pack
  use, intrinsic :: iso_c_binding
  implicit none
  private

  public :: mosflm_mm_to_10um
  public :: mosflm_10um_to_mm
  public :: mosflm_xtofd_to_mm
  public :: predict_unpack_gendata
  public :: predict_unpack_gendata_i2
  public :: predict_cpp_still_to_gendata
  public :: predict_cpp_still_to_gendata_i2

contains

  !> Crystal-to-detector ``XTOFD`` (10 µm units) → millimetres.
  elemental function mosflm_xtofd_to_mm(xtofd) result(mm)
    real(c_float), intent(in) :: xtofd
    real(c_double) :: mm
    mm = real(xtofd, kind=c_double) * 0.01_c_double
  end function mosflm_xtofd_to_mm

  !> Convert beam-relative millimetres to MOSFLM 10-micron units.
  elemental function mosflm_mm_to_10um(mm) result(u)
    real(c_double), intent(in) :: mm
    real(c_float) :: u
    u = real(mm * 100.0_c_double, kind=c_float)
  end function mosflm_mm_to_10um

  !> Convert MOSFLM 10-micron units to millimetres.
  elemental function mosflm_10um_to_mm(u) result(mm)
    real(c_float), intent(in) :: u
    real(c_double) :: mm
    mm = real(u, kind=c_double) * 0.01_c_double
  end function mosflm_10um_to_mm

  !> Scatter bridge outputs into GENDATA-shaped arrays; set ``nspot_out`` /
  !! ``iflagpr`` (−1 on truncate / error).
  subroutine predict_unpack_gendata(h, k, l, x_mm, y_mm, phi_beg, phi_end, &
                                    ifull, nspot, nmax, ierr, ihg, ikg, ilg, &
                                    xg, yg, phig, phiwg, img, nspot_out, &
                                    iflagpr)
    integer(c_int), intent(in) :: h(*)
    integer(c_int), intent(in) :: k(*)
    integer(c_int), intent(in) :: l(*)
    real(c_double), intent(in) :: x_mm(*)
    real(c_double), intent(in) :: y_mm(*)
    real(c_double), intent(in) :: phi_beg(*)
    real(c_double), intent(in) :: phi_end(*)
    integer(c_int), intent(in) :: ifull(*)
    integer(c_int), intent(in) :: nspot
    integer(c_int), intent(in) :: nmax
    integer(c_int), intent(in) :: ierr
    integer(c_int), intent(out) :: ihg(*)
    integer(c_int), intent(out) :: ikg(*)
    integer(c_int), intent(out) :: ilg(*)
    real(c_float), intent(out) :: xg(*)
    real(c_float), intent(out) :: yg(*)
    real(c_float), intent(out) :: phig(*)
    real(c_float), intent(out) :: phiwg(*)
    integer(c_int), intent(out) :: img(*)
    integer(c_int), intent(out) :: nspot_out
    integer(c_int), intent(out) :: iflagpr

    integer(c_int) :: i, nwrite

    iflagpr = 0_c_int
    nspot_out = 0_c_int
    if (nmax <= 0) then
      iflagpr = -1_c_int
      return
    end if

    if (ierr /= 0_c_int .and. ierr /= 2_c_int) then
      ! MOSFLM_ERROR_BUFFER_TOO_SMALL == 2: still write what we have.
      iflagpr = -1_c_int
      return
    end if

    nwrite = nspot
    if (nwrite > nmax) nwrite = nmax
    if (ierr == 2_c_int) iflagpr = -1_c_int

    do i = 1, nwrite
      ihg(i) = h(i)
      ikg(i) = k(i)
      ilg(i) = l(i)
      ! XG ∥ lab Y ↔ C++ y_mm; YG ∥ lab Z ↔ C++ x_mm.
      xg(i) = mosflm_mm_to_10um(y_mm(i))
      yg(i) = mosflm_mm_to_10um(x_mm(i))
      phig(i) = real(phi_beg(i), kind=c_float)
      phiwg(i) = real(phi_end(i) - phi_beg(i), kind=c_float)
      if (ifull(i) /= 0_c_int) then
        img(i) = 0_c_int
      else
        img(i) = 50_c_int
      end if
    end do
    nspot_out = nwrite
  end subroutine predict_unpack_gendata

  !> Run ``predict_cpp_still_psi`` and unpack into GENDATA-shaped arrays.
  subroutine predict_cpp_still_to_gendata(cell, distance_mm, wavelength_A, &
                                          pixel_x_mm, pixel_y_mm, nfast, &
                                          nslow, beam_x_mm, beam_y_mm, &
                                          phi_deg, dmin_A, ewald_tol, &
                                          psi_deg, nmax, ihg, ikg, ilg, xg, &
                                          yg, phig, phiwg, img, nspot, &
                                          iflagpr, ierr)
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
    integer(c_int), intent(in) :: nmax
    integer(c_int), intent(out) :: ihg(*)
    integer(c_int), intent(out) :: ikg(*)
    integer(c_int), intent(out) :: ilg(*)
    real(c_float), intent(out) :: xg(*)
    real(c_float), intent(out) :: yg(*)
    real(c_float), intent(out) :: phig(*)
    real(c_float), intent(out) :: phiwg(*)
    integer(c_int), intent(out) :: img(*)
    integer(c_int), intent(out) :: nspot
    integer(c_int), intent(out) :: iflagpr
    integer(c_int), intent(out) :: ierr

    integer(c_int), allocatable :: h(:), k(:), l(:), ifull(:)
    real(c_double), allocatable :: x_mm(:), y_mm(:), phi_beg(:), phi_end(:)
    integer(c_int) :: n_bridge

    external :: predict_cpp_still_psi

    iflagpr = 0_c_int
    nspot = 0_c_int
    ierr = 1_c_int
    if (nmax <= 0) then
      iflagpr = -1_c_int
      return
    end if

    allocate (h(nmax), k(nmax), l(nmax), ifull(nmax))
    allocate (x_mm(nmax), y_mm(nmax), phi_beg(nmax), phi_end(nmax))

    call predict_cpp_still_psi(cell, distance_mm, wavelength_A, pixel_x_mm, &
                               pixel_y_mm, nfast, nslow, beam_x_mm, &
                               beam_y_mm, phi_deg, dmin_A, ewald_tol, &
                               psi_deg, h, k, l, x_mm, y_mm, phi_beg, &
                               phi_end, ifull, nmax, n_bridge, ierr)

    call predict_unpack_gendata(h, k, l, x_mm, y_mm, phi_beg, phi_end, &
                                ifull, n_bridge, nmax, ierr, ihg, ikg, ilg, &
                                xg, yg, phig, phiwg, img, nspot, iflagpr)
  end subroutine predict_cpp_still_to_gendata

  !> Like ``predict_unpack_gendata`` but miller / ``IMG`` are ``INTEGER*2``
  !! (``c_short``), matching amalgamation ``/GENDATA/``.
  subroutine predict_unpack_gendata_i2(h, k, l, x_mm, y_mm, phi_beg, &
                                       phi_end, ifull, nspot, nmax, ierr, &
                                       ihg, ikg, ilg, xg, yg, phig, phiwg, &
                                       img, nspot_out, iflagpr)
    integer(c_int), intent(in) :: h(*)
    integer(c_int), intent(in) :: k(*)
    integer(c_int), intent(in) :: l(*)
    real(c_double), intent(in) :: x_mm(*)
    real(c_double), intent(in) :: y_mm(*)
    real(c_double), intent(in) :: phi_beg(*)
    real(c_double), intent(in) :: phi_end(*)
    integer(c_int), intent(in) :: ifull(*)
    integer(c_int), intent(in) :: nspot
    integer(c_int), intent(in) :: nmax
    integer(c_int), intent(in) :: ierr
    integer(c_short), intent(out) :: ihg(*)
    integer(c_short), intent(out) :: ikg(*)
    integer(c_short), intent(out) :: ilg(*)
    real(c_float), intent(out) :: xg(*)
    real(c_float), intent(out) :: yg(*)
    real(c_float), intent(out) :: phig(*)
    real(c_float), intent(out) :: phiwg(*)
    integer(c_short), intent(out) :: img(*)
    integer(c_int), intent(out) :: nspot_out
    integer(c_int), intent(out) :: iflagpr

    integer(c_int), allocatable :: ih32(:), ik32(:), il32(:), im32(:)
    integer(c_int) :: i

    allocate (ih32(max(nmax, 1)), ik32(max(nmax, 1)), il32(max(nmax, 1)), &
              im32(max(nmax, 1)))
    call predict_unpack_gendata(h, k, l, x_mm, y_mm, phi_beg, phi_end, &
                                ifull, nspot, nmax, ierr, ih32, ik32, &
                                il32, xg, yg, phig, phiwg, im32, nspot_out, &
                                iflagpr)
    do i = 1, nspot_out
      ihg(i) = int(ih32(i), kind=c_short)
      ikg(i) = int(ik32(i), kind=c_short)
      ilg(i) = int(il32(i), kind=c_short)
      img(i) = int(im32(i), kind=c_short)
    end do
  end subroutine predict_unpack_gendata_i2

  !> Stills PREDICT → GENDATA with ``INTEGER*2`` miller / ``IMG``.
  subroutine predict_cpp_still_to_gendata_i2(cell, distance_mm, &
                                             wavelength_A, pixel_x_mm, &
                                             pixel_y_mm, nfast, nslow, &
                                             beam_x_mm, beam_y_mm, phi_deg, &
                                             dmin_A, ewald_tol, psi_deg, &
                                             nmax, ihg, ikg, ilg, xg, yg, &
                                             phig, phiwg, img, nspot, &
                                             iflagpr, ierr)
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
    integer(c_int), intent(out) :: iflagpr
    integer(c_int), intent(out) :: ierr

    integer(c_int), allocatable :: h(:), k(:), l(:), ifull(:)
    real(c_double), allocatable :: x_mm(:), y_mm(:), phi_beg(:), phi_end(:)
    integer(c_int) :: n_bridge

    external :: predict_cpp_still_psi

    iflagpr = 0_c_int
    nspot = 0_c_int
    ierr = 1_c_int
    if (nmax <= 0) then
      iflagpr = -1_c_int
      return
    end if

    allocate (h(nmax), k(nmax), l(nmax), ifull(nmax))
    allocate (x_mm(nmax), y_mm(nmax), phi_beg(nmax), phi_end(nmax))

    call predict_cpp_still_psi(cell, distance_mm, wavelength_A, pixel_x_mm, &
                               pixel_y_mm, nfast, nslow, beam_x_mm, &
                               beam_y_mm, phi_deg, dmin_A, ewald_tol, &
                               psi_deg, h, k, l, x_mm, y_mm, phi_beg, &
                               phi_end, ifull, nmax, n_bridge, ierr)

    call predict_unpack_gendata_i2(h, k, l, x_mm, y_mm, phi_beg, phi_end, &
                                   ifull, n_bridge, nmax, ierr, ihg, ikg, &
                                   ilg, xg, yg, phig, phiwg, img, nspot, &
                                   iflagpr)
  end subroutine predict_cpp_still_to_gendata_i2

end module predict_gendata_pack
