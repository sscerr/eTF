      SUBROUTINE scatter_real(SEND,RECV,SIZE,root)
      IMPLICIT NONE
      INTEGER SIZE, ROOT
      REAL*8 SEND(*), RECV(*)
      INTEGER I
      do i = 1, size
         recv( i ) = send( i )
      end do
      RETURN
      END

