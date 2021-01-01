! ------------------------------------------------------ !
!  F. Califano, Marzo 2006, M. Faganello, Marzo 2007     !
! Solver equazione di Poisson nel piano x,y.             !
! (1 - de^2 Nabla^2) E = Ohm(right hand side)            !
! Componenti Ex, Ey                                      !
! ------------------------------------------------------ !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!   3D PARALLEL VERSION: FAGANELLO 2010
!**************************************************

subroutine poisson( RHS )

use box_mod
use deriv_mod
use parameter_mod
use dom_distr_mod
use poisson_open
use fields_UJ_mod
use fields_EB_mod
use parallel_mod

IMPLICIT NONE

integer :: ix, iy, iz, ily, ilz

REAL(dp), DIMENSION ( nxl, ny, nz ) :: RHS

REAL(dp),ALLOCATABLE ::  sy(:), sz(:)
REAL(dp),ALLOCATABLE ::  syy1(:,:), syy2(:,:)
REAL(dp),ALLOCATABLE ::  g1(:)  
REAL(dp),ALLOCATABLE ::  Rt(:,:,:)
!real*8 :: kkyy, kkzz

allocate(g1(nx))
allocate(sy(ny),sz(nz))
allocate(syy1(ny,nz))
allocate(syy2(ny,nz))
allocate(Rt(nx, nyt, nz))


!test
!kkyy=2.*pgreco/yl
!kkzz=2.*pgreco/zl
!do iz = 1, nz
! do iy = 1, ny
!  do ix = 1, nxl
!   Ex(ix,iy,iz)=tanh(5*(x(ix+ixlg-1)-xl/2.))*cos(kkyy*y(iy))*sin(kkzz*z(iz)) 
!  enddo
! enddo
!enddo

!call outx

!do iz = 1, nz
! do iy = 1, ny
!  do ix = 1, nxl
!   RHS(ix,iy,iz)=tanh(5.*(x(ix+ixlg-1)-xl/2.))*cos(kkyy*y(iy))*sin(kkzz*z(iz)) - de2 *&
!(-2.*25.*tanh(5.*(x(ix+ixlg-1)-xl/2.))/(cosh(5.*(x(ix+ixlg-1)-xl/2.)))**2 *cos(kkyy*y(iy))*sin(kkzz*z(iz))&
!- kkyy**2 * tanh(5.*(x(ix+ixlg-1)-xl/2.))*cos(kkyy*y(iy))*sin(kkzz*z(iz)) &
!- kkzz**2 * tanh(5.*(x(ix+ixlg-1)-xl/2.))*cos(kkyy*y(iy))*sin(kkzz*z(iz)) ) 
!  enddo
! enddo
!enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ALTERNATIVA PER BOUNDARY CONDITIONS - DA RIVEDERE --------------------------------------------------

do iz = 1, nz
 do ix = 1, nxl
     sy = RHS(ix, :, iz)
     call drfftf(ny, sy, savey)
     RHS(ix, :, iz) = sy
  enddo
enddo

IF (nz.GT.1) THEN
do iy = 1, ny
 do ix = 1, nxl
     sz = RHS(ix, iy, :)
     call drfftf(nz, sz, savez)
     RHS(ix, iy, :) = sz
  enddo
enddo
ENDIF

 CALL traspdist( RHS, Rt, 1 )

do iz = 1, nz
 do iy=1,nyt
  syy1(iy+iylg-1,iz)=Rt(1,iy,iz)
  syy2(iy+iylg-1,iz)=Rt(nx,iy,iz)
 enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do iz=1,nz
 do iy=1,nyt

ily=int(0.5*(iy+iytg-1))
IF(ily==ny/2) ily=0

ilz=int(0.5*iz)
IF(ilz==nz/2) ilz=0

  do ix = 2, nx - 1
    g1(ix) = pss_rhs_1(ily,ilz) * (Rt(ix-1,iy,iz) + Rt(ix+1,iy,iz))  &
           + pss_rhs_0(ily,ilz) *  Rt(ix,iy,iz)
  enddo

! Condizioni al contorno
! E = cost (vedi subroutine inizializza, riga 190)
! Imponiamo E = RHS_bordi

 g1(1)	= syy1(iy,iz)
 g1(nx) = syy2(iy,iz)

! Condizioni al contorno
! dE/dx = 0.0 (vedi subroutine inizializza, riga 190)

! g1(1) = 0.0d0
! g1(nx)= 0.0d0


 CALL  DGTTRS(TRNAS,nx,irhs,poiss_m(:,ily,ilz),poiss_d(:,ily,ilz),poiss_p(:,ily,ilz),poiss_w,ipv_p(:,ily,ilz),g1,idb,ifno)

 if (ifno > 0 .or. ifno < 0) then
    write(*,*) 'Problemi soluzione, ifno:', ifno
    stop
 endif

 Rt(:,iy,iz)=g1

 enddo
enddo


 CALL traspdist( RHS, Rt, -1 )
 RHS(:,ny,:)=0.d0
IF (nz.GT.1) THEN
 RHS(:,:,nz)=0.d0
ENDIF
 
do iz = 1, nz
 do ix = 1, nxl
   sy = RHS(ix, :,iz)
  call drfftb( ny, sy, savey )
   RHS( ix, :, iz) = sy * Tot_inv_ny
 enddo
enddo 

IF (nz.GT.1) THEN
do iy = 1, ny
 do ix = 1, nxl
   sz = RHS(ix, iy, :)
  call drfftb( nz, sz, savez )
   RHS( ix, iy, :) = sz * Tot_inv_nz
 enddo
enddo 
ENDIF

!test
!Ex=RHS
!call outx




deallocate(g1)
deallocate(sy)
deallocate(sz)
deallocate(syy1)
deallocate(syy2)
deallocate(Rt)





  end subroutine  
