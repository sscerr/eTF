       integer function LOCALINDEX(ig,nx,npx,pex)
!
!   INPUT :
!           ig  global index of the x dimension of array element
!          
!           nx  x dimension of the global array
!
!           npx  number of processor in the x dimension of the processors grid
!
!           pex  x index of the local processor in the processor grid
!                (starting from zero)
!
!   OUTPUT :
!
!           localindex  index of the element in the local block
!   


       IMPLICIT NONE

       INTEGER il,ig,nx,npx,pex,r,q

       Q = INT(nx/npx)
       R = MOD(nx,npx) 
       IF((pex+1).LE.R) THEN
         localindex = ig - (Q+1)*pex
       ELSE
         localindex = ig - Q*pex + R
       end if

       return 
       end
