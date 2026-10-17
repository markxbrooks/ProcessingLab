! Amalgamation-style PREDICT_CPP: XTOFD in 10 µm, I*2 GENDATA, n=656.

program test_predict_cpp_hook
  use, intrinsic :: iso_c_binding
  use predict_gendata_pack, only: mosflm_xtofd_to_mm
  implicit none

  integer(c_int), parameter :: nmax = 8192_c_int
  integer(c_int), parameter :: expect = 656_c_int

  real(c_float) :: psi(3), cell(6)
  real(c_float) :: xtofd, wave, rast, rasty, phi, dmin, etol
  integer(c_int) :: nf, ns, nspot, iflagpr
  integer(c_short) :: ihg(nmax), ikg(nmax), ilg(nmax), img(nmax)
  real(c_float) :: xg(nmax), yg(nmax), phig(nmax), phiwg(nmax)

  external :: predict_cpp

  if (abs(mosflm_xtofd_to_mm(25000.0_c_float) - 250.0_c_double) > &
      1.0e-9_c_double) then
    write (*, '(A)') 'FAIL: XTOFD 25000 → 250 mm'
    stop 1
  end if

  cell = [58.0_c_float, 58.0_c_float, 150.0_c_float, 90.0_c_float, &
          90.0_c_float, 90.0_c_float]
  xtofd = 25000.0_c_float
  wave = 1.0_c_float
  rast = 0.172_c_float
  rasty = 0.172_c_float
  nf = 1475_c_int
  ns = 1679_c_int
  phi = 0.0_c_float
  dmin = 3.5_c_float
  etol = 0.0025_c_float
  psi = [0.0_c_float, 0.0_c_float, 0.0_c_float]

  call predict_cpp(psi, iflagpr, cell, xtofd, wave, rast, rasty, nf, ns, &
                   phi, dmin, etol, nmax, ihg, ikg, ilg, xg, yg, phig, &
                   phiwg, img, nspot)

  if (iflagpr /= 0) then
    write (*, '(A,I0)') 'FAIL: iflagpr=', iflagpr
    stop 1
  end if
  if (nspot /= expect) then
    write (*, '(A,I0)') 'FAIL: nspot=', nspot
    stop 1
  end if
  if (img(1) /= 0_c_short .and. img(1) /= 50_c_short) then
    write (*, '(A,I0)') 'FAIL: img kind/value=', int(img(1))
    stop 1
  end if

  psi = [5.0_c_float, -2.0_c_float, 1.0_c_float]
  block
    integer(c_int) :: n_rot, flg
    call predict_cpp(psi, flg, cell, xtofd, wave, rast, rasty, nf, ns, &
                     phi, dmin, etol, nmax, ihg, ikg, ilg, xg, yg, phig, &
                     phiwg, img, n_rot)
    if (flg /= 0 .or. n_rot == nspot .or. n_rot <= 0) then
      write (*, '(A,I0,A,I0)') 'FAIL: rotated n=', n_rot, ' id=', nspot
      stop 1
    end if
  end block

  write (*, '(A,I0)') 'test_predict_cpp_hook: OK n=', nspot
  stop 0
end program test_predict_cpp_hook
