! Offline: C++ stills PREDICT → GENDATA-shaped arrays (n=656, XG=mm*100).

program test_predict_gendata_pack
  use, intrinsic :: iso_c_binding
  use predict_gendata_pack
  implicit none

  integer(c_int), parameter :: nmax = 8192_c_int
  integer(c_int), parameter :: expect = 656_c_int
  real(c_float) :: cell(6)
  real(c_double) :: dist, wave, px, py, bx, by, phi, dmin, etol
  real(c_double) :: psi(3)
  integer(c_int) :: nf, ns, nspot, iflagpr, ierr
  integer(c_int) :: ihg(nmax), ikg(nmax), ilg(nmax), img(nmax)
  real(c_float) :: xg(nmax), yg(nmax), phig(nmax), phiwg(nmax)

  ! Round-trip unit helpers (value exact in float32 ×100).
  if (abs(mosflm_10um_to_mm(mosflm_mm_to_10um(1.25_c_double)) - &
          1.25_c_double) > 1.0e-6_c_double) then
    write (*, '(A)') 'FAIL: mm ↔ 10um round-trip'
    stop 1
  end if

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

  call predict_cpp_still_to_gendata(cell, dist, wave, px, py, nf, ns, bx, &
                                    by, phi, dmin, etol, psi, nmax, ihg, &
                                    ikg, ilg, xg, yg, phig, phiwg, img, &
                                    nspot, iflagpr, ierr)

  if (ierr /= 0) then
    write (*, '(A,I0)') 'FAIL: ierr=', ierr
    stop 1
  end if
  if (iflagpr /= 0) then
    write (*, '(A,I0)') 'FAIL: iflagpr=', iflagpr
    stop 1
  end if
  if (nspot /= expect) then
    write (*, '(A,I0)') 'FAIL: nspot=', nspot
    stop 1
  end if

  ! Direct unpack smoke: XG↔y_mm, YG↔x_mm (10 µm = mm*100).
  block
    integer(c_int) :: h1(1), k1(1), l1(1), if1(1), n_out, flg
    real(c_double) :: xm(1), ym(1), pb(1), pe(1)
    integer(c_int) :: ih1(1), ik1(1), il1(1), im1(1)
    real(c_float) :: xg1(1), yg1(1), pg1(1), pw1(1)
    h1 = 1; k1 = 2; l1 = 3; if1 = 1
    xm = 1.5_c_double; ym = -2.0_c_double
    pb = 0.0_c_double; pe = 0.0_c_double
    call predict_unpack_gendata(h1, k1, l1, xm, ym, pb, pe, if1, 1_c_int, &
                                1_c_int, 0_c_int, ih1, ik1, il1, xg1, yg1, &
                                pg1, pw1, im1, n_out, flg)
    if (n_out /= 1 .or. flg /= 0 .or. ih1(1) /= 1 .or. ik1(1) /= 2 .or. &
        il1(1) /= 3) then
      write (*, '(A)') 'FAIL: unpack indices'
      stop 1
    end if
    if (abs(real(xg1(1), kind=c_double) + 200.0_c_double) > 1.0e-3_c_double &
        .or. abs(real(yg1(1), kind=c_double) - 150.0_c_double) > &
        1.0e-3_c_double) then
      write (*, '(A)') 'FAIL: unpack XG/YG (expect XG=y*100, YG=x*100)'
      stop 1
    end if
    if (im1(1) /= 0_c_int) then
      write (*, '(A)') 'FAIL: full → IMG=0'
      stop 1
    end if
  end block

  if (img(1) /= 0_c_int .and. img(1) /= 50_c_int) then
    write (*, '(A,I0)') 'FAIL: img(1)=', img(1)
    stop 1
  end if
  if (abs(phiwg(1)) > 1.0e-3_c_float .and. img(1) == 0_c_int) then
    write (*, '(A)') 'FAIL: full still should have ~0 PHIWG'
    stop 1
  end if

  ! Non-zero PSI must change count (same as predict_cpp_psi).
  block
    integer(c_int) :: n_rot, iflag2, ierr2
    psi = [5.0_c_double, -2.0_c_double, 1.0_c_double]
    call predict_cpp_still_to_gendata(cell, dist, wave, px, py, nf, ns, bx, &
                                      by, phi, dmin, etol, psi, nmax, ihg, &
                                      ikg, ilg, xg, yg, phig, phiwg, img, &
                                      n_rot, iflag2, ierr2)
    if (ierr2 /= 0 .or. n_rot == nspot .or. n_rot <= 0) then
      write (*, '(A,I0,A,I0)') 'FAIL: rotated n=', n_rot, ' id=', nspot
      stop 1
    end if
  end block

  write (*, '(A,I0)') 'test_predict_gendata_pack: OK n=', nspot
  stop 0
end program test_predict_gendata_pack
