      SUBROUTINE gather_real(SEND,RECV,SIZE,root)
      IMPLICIT NONE
      INTEGER SIZE, ROOT
      REAL*8 SEND(*), RECV(*)
      INTEGER I
      do i = 1, size
         recv( i ) = send( i )
      end do 
      RETURN
      END

      SUBROUTINE gatherv_real( SEND, ssize, RECV, rsize, rdisp, root )
      IMPLICIT NONE
      INTEGER ssize, ROOT
      INTEGER rsize(*), rdisp(*)
      REAL*8 SEND(*), RECV(*)
      INTEGER I
      do i = 1, ssize
         recv( i ) = send( i )
      end do 
      RETURN
      END

