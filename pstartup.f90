subroutine parallel_init

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!**************************************************

use parameter_mod   !, only: nx, ny, space_dim, nxl, nyl
use parallel_mod    !, only: npcol, nprow, mycol, myrow, &
                    !   mpime,nproc,group,root, parallel_build

use dom_distr_mod

IMPLICIT NONE

integer :: globalindex
integer :: ip, localdim

!-----------------------------------------------------------------
!   nx =number of grid points in the open x-direction
!   ny+1=number of grid points in the periodic y-direction 
!----------------------------------------------------------------


if( nproc .gt. nprocx ) then
    write(*,*) '*** too many processors ***'
    write(*,*) '***   increase nprocx   ***'
    call hangup
    stop
end if

! ---------- Distribuzione dati 1D (taglio solo lungo x) -----------

nprow = nproc
npcol = 1

! -----------

myrow = mod(mpime,nprow)
mycol = mpime / nprow

nxg = nx
nyg = ny
!nzg va ancora def nei moduli....per ora pongo solo nzl = nz

nxl = localdim(nxg,nprow,myrow)
nyl = localdim(nyg,npcol,mycol)
nzl = nz

if (nz.GT.1) then
 n2z = nz/2
else
 n2z = 1
endif

nyt = localdim(nyg,nprow,myrow)  ! distribuzione trasposta


!al solito mancano ancora tutti #z
ixlg = globalindex(1,nxg,nprow,myrow)
iylg = globalindex(1,nyg,npcol,mycol)

iytg = globalindex(1,nyg,nprow,myrow)  ! distribuzione trasposta

nxlm1 = nxl - 1
nxlm2 = nxl - 2
nxlm3 = nxl - 3
nxlm4 = nxl - 4

nylm1 = nyl - 1
nylm2 = nyl - 2
nylm3 = nyl - 3
nylm4 = nyl - 4

!manca ancora #z vari.
!if( space_dim == '3D' ) then
  !!!   write(*, 1000) mpime, nproc, myrow, mycol, nprow, npcol
  !!!   write(*, 1001) mpime, nxl, nyl, nzl, ixlg, iytg, nxg
!else
!     write(*, 1000) mpime, nproc, myrow, mycol, nprow, npcol
!     write(*, 1002) mpime, nxl, ixlg, nxg
!     if( nyl /= ny ) then
!          write(*,*) ' nyl different from ny in 1D? '
!          call hangup()
!          stop 
!     end if
!end if

1000   FORMAT(' mpime, nproc, myrow, mycol, nprow, npcol',6I4) 
1001   FORMAT(' mpime, nxl, nyl, ixlg, iytg, nxg',7I4) 
1002   FORMAT(' mpime, nxl, ixlg, nxg',4I4) 


!manca warning in z
if( nxl .lt. 5 ) then
    write(*,*) ' nxl less than 5, use less processors '
    call hangup()
    stop 
end if
if( ( space_dim == '3D' ) .and. ( nyl .lt. 5 ) ) then
    write(*,*) ' nyl less than 5, use less processors '
    call hangup()
    stop 
end if


!manca #z
do ip = 1, nproc
    iprow(ip) = mod((ip-1),nprow)
    ipcol(ip) = (ip-1) / nprow
    ixlgp(ip) = globalindex(1,nxg,nprow,iprow(ip))
    iylgp(ip) = globalindex(1,nyg,npcol,ipcol(ip))
    iytgp(ip) = globalindex(1,nyg,nprow,iprow(ip))
    nxlp(ip) = localdim(nxg,nprow,iprow(ip))
    nylp(ip) = localdim(nyg,npcol,ipcol(ip))
    nytp(ip) = localdim(nyg,nprow,iprow(ip))
!   IF(mpime .EQ. root) THEN
!       write(*,'(10X,7I5)') ip, iprow(ip), ipcol(ip), &
!            ixlgp(ip), iylgp(ip), nxp(ip), nyp(ip)
!   END IF
end do
   
    nxlmax = maxval(nxlp)
    nytmax = maxval(nytp)

end subroutine 
