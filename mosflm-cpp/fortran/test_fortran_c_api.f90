! Differential smoke: Fortran → MosflmCApi stills PREDICT.
! Expects the analytical fixture spot count to match the frozen oracle (656).

program test_fortran_c_api
  use, intrinsic :: iso_c_binding
  use mosflm_c_bind
  implicit none

  type(mosflm_cell6) :: cell
  type(mosflm_detector_c) :: det
  type(mosflm_predicted_c) :: spots(MOSFLM_C_MAX_SPOTS)
  real(c_double) :: phi, dmin, etol, lp0, lp_off
  integer(c_int) :: nspots, rc
  integer :: i
  logical :: found_nonzero

  ! Expected count from testdata/oracle/predict_still.dump
  integer, parameter :: expect_n = 656

  call predict_still_cpp(.true., cell, det, phi, dmin, etol, spots, nspots, rc)

  if (rc /= 0) then
    write (*, '(A,I0)') 'FAIL: predict_still_cpp rc=', rc
    stop 1
  end if
  if (nspots /= expect_n) then
    write (*, '(A,I0,A,I0)') 'FAIL: nspots=', nspots, ' expected ', expect_n
    stop 1
  end if
  if (abs(cell%a - 58.0_c_double) > 1.0e-9_c_double) then
    write (*, '(A)') 'FAIL: fixture cell a'
    stop 1
  end if
  if (det%nfast /= 1475_c_int) then
    write (*, '(A)') 'FAIL: fixture nfast'
    stop 1
  end if

  found_nonzero = .false.
  do i = 1, int(nspots)
    if (spots(i)%h /= 0 .or. spots(i)%k /= 0 .or. spots(i)%l /= 0) then
      found_nonzero = .true.
      exit
    end if
  end do
  if (.not. found_nonzero) then
    write (*, '(A)') 'FAIL: all hkl zero'
    stop 1
  end if

  lp0 = mosflm_lp_factor(det, 0.0_c_double, 0.0_c_double)
  if (abs(lp0 - 1.0_c_double) > 1.0e-6_c_double) then
    write (*, '(A,ES12.5)') 'FAIL: LP at beam=', lp0
    stop 1
  end if
  lp_off = mosflm_lp_factor(det, 50.0_c_double, 0.0_c_double)
  if (lp_off <= 1.0_c_double) then
    write (*, '(A)') 'FAIL: LP off-axis not > 1'
    stop 1
  end if

  ! Oscillation path (predict_osc.dump n=368).
  block
    type(mosflm_cell6) :: ocell
    type(mosflm_detector_c) :: odet
    type(mosflm_predicted_c) :: ospots(MOSFLM_C_MAX_SPOTS)
    real(c_double) :: p0, p1, odmin, oetol, ostep
    integer(c_int) :: on, orc
    integer, parameter :: expect_osc = 368
    call predict_osc_cpp(.true., ocell, odet, p0, p1, odmin, oetol, ostep, &
                         ospots, on, orc)
    if (orc /= 0 .or. on /= expect_osc) then
      write (*, '(A,I0,A,I0)') 'FAIL: osc n=', on, ' rc=', orc
      stop 1
    end if
  end block

  write (*, '(A,I0,A,ES12.5)') &
    'test_fortran_c_api: OK (n=', nspots, ' lp_off=', lp_off, ')'
  stop 0
end program test_fortran_c_api
