      program readit

C Arrays for x, y and z coordinates of the reflections
      parameter (n=5000)
      real a(n),b(n),c(n)
      real ma(3,3)

C Number of reflections read.
      integer num

C Maximum expected cell edge length
      dmax=300.

C Mode of indexing (only one exiting at the moment
      mode=0

C Read the reflections (see below)
      call read_peaks (a,b,c,n,num)

C Call the indexing
      call index_fortran(num,a,b,c,dmax,mode,ma)

      write(*,999) ((ma(i,j),j=1,3),i=1,3)
999   format(3(3f10.7,/))
      end

      subroutine read_peaks (a,b,c,n,i)

      real a(n),b(n),c(n)
      integer i

      open (1,file='input',form='formatted',err=200)

      i = 0

1     read (1, *,end=100,err=201) a(i+1),b(i+1),c(i+1)

      i = i + 1

      goto 1

100   type *, '# of lines in file: ', i

      goto 2
200   stop 'error openning input file'
201   stop 'error reading input file'

2     return

      end

