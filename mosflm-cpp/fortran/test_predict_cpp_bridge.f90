! Smoke: amalgamation-named PREDICT_CPP_STILL → C++ stills PREDICT (n=656).

program test_predict_cpp_bridge
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int), parameter :: nmax = 8192_c_int
  integer(c_int), parameter :: expect = 656_c_int

  real(c_float) :: cell(6)
  real(c_double) :: dist, wave, px, py, bx, by, phi, dmin, etol
  integer(c_int) :: nf, ns, nspot, ierr
  integer(c_int) :: h(nmax), k(nmax), l(nmax), ifull(nmax)
  real(c_double) :: xmm(nmax), ymm(nmax), phib(nmax), phie(nmax)

  external :: predict_cpp_still

  cell = [58.0_c_float, 58.0_c_float, 150.0_c_float, 90.0_c_float, &
          90.0_c_float, 90.0_c_float]
  dist = 250.0_c_double
  wave = 1.0_c_double
  px = 0.172_c_double
  py = 0.172_c_double
  nf = 1475_c_int
  ns = 1679_c_int
  bx = 0.0_c_double
  by = 0.0_c_double
  phi = 0.0_c_double
  dmin = 3.5_c_double
  etol = 0.0025_c_double

  call predict_cpp_still(cell, dist, wave, px, py, nf, ns, bx, by, phi, &
                         dmin, etol, h, k, l, xmm, ymm, phib, phie, ifull, &
                         nmax, nspot, ierr)

  if (ierr /= 0) then
    write (*, '(A,I0)') 'FAIL: ierr=', ierr
    stop 1
  end if
  if (nspot /= expect) then
    write (*, '(A,I0,A,I0)') 'FAIL: nspot=', nspot, ' expected ', expect
    stop 1
  end if

  write (*, '(A,I0)') 'test_predict_cpp_bridge: OK n=', nspot
  stop 0
end program test_predict_cpp_bridge
