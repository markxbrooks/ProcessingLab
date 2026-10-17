      subroutine invert( a, b )
C *************************************************************** INVERT
 
C TEMPORARY --- REPLACE with QR factorization

      implicit none
 
      real a(3,3), b(3,3), c(3,3), w_1, w_2, w_3
      real work( 3, 6 ), xmult
 
      integer irow, icol, irow_1, irow_2, irow_3
 
      do 1010 irow = 1, 3
        do 1000 icol = 1, 3
          work(irow,icol) = a(irow,icol)
          work(irow,icol+3) = 0.
          c(irow,icol) = 0.
 1000   continue
        work(irow,irow+3) = 1.
 1010 continue
 
C Find pivot, first column
      irow_1 = 1
      w_1    = work(1,1)
      if( abs( w_1 ) .lt. abs( work(2,1) ) ) then
        irow_1 = 2
        w_1    = work(2,1)
      end if
      if( abs( w_1 ) .lt. abs( work(3,1) ) ) then
        irow_1 = 3
        w_1    = work(3,1)
      end if
 
C Make pivot unity
      do 1020 icol = 1, 6
        work(irow_1,icol) = work(irow_1,icol) / w_1
 1020 continue
 
C Eliminate first pivot
      do 1040 irow = 1, 3
        if( irow .ne. irow_1 ) then
          xmult = work(irow,1)
          do 1030 icol = 1, 6
            work(irow,icol) = work(irow,icol)
     +                      - xmult * work(irow_1,icol)
 1030     continue
        end if
 1040 continue
 
C Find second pivot
       irow_2 = 0
       w_2    = 0.
       do 1050 irow = 1, 3
         if( irow .ne. irow_1 ) then
           if( abs( w_2 ) .lt. abs( work(irow,2) ) ) then
             irow_2 = irow
             w_2    = work(irow,2)
           end if
         end if
 1050 continue

C Make pivot unity
      do 1060 icol = 1, 6
        work(irow_2,icol) = work(irow_2,icol) / w_2
 1060 continue
 
C Eliminate second pivot
      do 1080 irow = 1, 3
        if( irow .ne. irow_2 ) then
          xmult = work(irow,2)
          do 1070 icol = 1, 6
            work(irow,icol) = work(irow,icol)
     +                      - xmult * work(irow_2,icol)
 1070     continue
        end if
 1080 continue
 
C Find third pivot
      do 1090 irow = 1, 3
        if( ( irow .ne. irow_1 ) .and. (irow .ne. irow_2 ) ) then
          irow_3 = irow
          w_3    = work(irow,3)
        end if
 1090 continue
 
C Make pivot unity
      do 1100 icol = 1, 6
        work(irow_3,icol) = work(irow_3,icol) / w_3
 1100 continue
 
C Eliminate third pivot
      do 1120 irow = 1, 3
        if( irow .ne. irow_3 ) then
          xmult = work(irow,3)
          do 1110 icol = 1, 6
            work(irow,icol) = work(irow,icol)
     +                      - xmult * work(irow_3,icol)
 1110     continue
        end if
 1120 continue
 
C Transfer to b
      do 1140 icol = 1, 3
        b(1,icol) = work(irow_1,icol+3)
        b(2,icol) = work(irow_2,icol+3)
        b(3,icol) = work(irow_3,icol+3)
 1140 continue

      return
      end
