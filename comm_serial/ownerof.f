
       integer function  OWNEROF(ig,nx,npx)
!
!   INPUT :
!           ig  global index of the x dimension of array element
!          
!           nx  x dimension of the global array
!
!           npx  number of processor in the x dimension of the processors grid
!
!   OUTPUT :
!


       IMPLICIT NONE

       INTEGER ig,nx,npx,r,q,pex

       q = int(nx/npx)
       r = mod(nx,npx)
       if(r.gt.0) then
         if(ig.le.((q+1)*r)) then
            pex = int((ig-1)/(q+1))+1
         else
            pex = int((ig-1-r*(q+1))/q)+1+r
         end if
       else
         pex = int((ig-1)/q)+1
       end if

       ownerof = pex
 
       RETURN
       END
