      subroutine error_exit( prgnm, sbnm, message )
C *********************************************************** ERROR_EXIT
      implicit none

      character*(*) prgnm, sbnm, message

      print 1000, prgnm, sbnm, message
 1000 format(/' >>> ERROR in program ',a/
     +        '           subprogram ',a/'    ', a )

C      stop
      end
