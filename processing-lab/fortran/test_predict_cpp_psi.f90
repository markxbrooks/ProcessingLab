! Zero PSI matches identity bridge; non-zero PSI changes spot count.

program test_predict_cpp_psi
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int), parameter :: nmax = 8192_c_int
  integer(c_int), parameter :: expect_id = 656_c_int

  real(c_float) :: cell(6)
  real(c_double) :: dist, wave, px, py, bx, by, phi, dmin, etol
  real(c_double) :: psi(3)
  integer(c_int) :: nf, ns, n_id, n_rot, ierr
  integer(c_int) :: h(nmax), k(nmax), l(nmax), ifull(nmax)
  real(c_double) :: xmm(nmax), ymm(nmax), phib(nmax), phie(nmax)

  external :: predict_cpp_still_psi

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

  psi = [0.0_c_double, 0.0_c_double, 0.0_c_double]
  call predict_cpp_still_psi(cell, dist, wave, px, py, nf, ns, bx, by, phi, &
                             dmin, etol, psi, h, k, l, xmm, ymm, phib, phie, &
                             ifull, nmax, n_id, ierr)
  if (ierr /= 0) then
    write (*, '(A,I0)') 'FAIL: psi0 ierr=', ierr
    stop 1
  end if
  if (n_id /= expect_id) then
    write (*, '(A,I0)') 'FAIL: psi0 nspot=', n_id
    stop 1
  end if

  psi = [5.0_c_double, -2.0_c_double, 1.0_c_double]
  call predict_cpp_still_psi(cell, dist, wave, px, py, nf, ns, bx, by, phi, &
                             dmin, etol, psi, h, k, l, xmm, ymm, phib, phie, &
                             ifull, nmax, n_rot, ierr)
  if (ierr /= 0) then
    write (*, '(A,I0)') 'FAIL: psi_rot ierr=', ierr
    stop 1
  end if
  if (n_rot == n_id .or. n_rot <= 0) then
    write (*, '(A,I0,A,I0)') 'FAIL: rotated n=', n_rot, ' identity=', n_id
    stop 1
  end if

  write (*, '(A,I0,A,I0)') 'test_predict_cpp_psi: OK id=', n_id, &
                          ' rotated=', n_rot
  stop 0
end program test_predict_cpp_psi
