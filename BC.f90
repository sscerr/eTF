       subroutine  BC

!----------------------------------- !
!       M. Faganello, 2007           !
!                                    !
! Boundary Condition basate sulle    !
! onde caratteristiche MHD, lungo    !
!     direzione disomogenea x        !                           
!                                    !
! Calcolo L_vari, vedi Landi2005,    !
!     Del_Zanna2001                  !
!                                    !
!----------------------------------- !

!**************************************************
!  MPI PARALLEL VERSION: VALENTINI-FAGANELLO 2009
!   3D PARALLEL VERSION: FAGANELLO 2010
!**************************************************

use parameter_mod
use box_mod
use deriv_mod
use fields_UJ_mod
use fields_DP_mod
use fields_EB_mod
use dom_distr_mod

IMPLICIT NONE

integer :: ix,iy,iz

REAL(dp), ALLOCATABLE :: grad(:,:,:), T_tot(:,:,:)
REAL(dp), ALLOCATABLE :: ax_2(:,:,:), ay_2(:,:,:), az_2(:,:,:), c_s_2(:,:,:)

allocate( grad(2,ny,nz) )
allocate( ax_2(2,ny,nz) )
allocate( ay_2(2,ny,nz) )
allocate( az_2(2,ny,nz) )
allocate( c_s_2(2,ny,nz) )
allocate( T_tot(nxl,ny,nz) )

!**************NON ISOTERMI*****************************
!
!non-reflecting boundary condition...di conseguenza L_vari_+ = 0 in ix=1
!e L_vari_- = 0 in ix=nx.
!gli altri vanno calcolati utilizzando i punti interni per la derivata.



! *********** calcolo velocita' caratteristiche necessarie ai vari L *****
! in realta' usero' poi solo 1 del primo processore e 2(nxl) del'ultimo

T_tot = (Ti_para + Te_para + 2.0d0*(Ti_perp + Te_perp)) * terzo

!***************        ix=1       ***************************************

ax(1,:,:) = dsqrt( Dinv(1,:,:)  ) * Bx(1,:,:)
ay(1,:,:) = dsqrt( Dinv(1,:,:)  ) * By(1,:,:)
az(1,:,:) = dsqrt( Dinv(1,:,:)  ) * Bz(1,:,:)
segno(1,:,:) = sign(dble(1.0),Bx(1,:,:))
 c_s(1,:,:) = dsqrt(gam * T_tot(1,:,:))

!***************        ix=nxl       ***************************************

ax(2,:,:) = dsqrt( Dinv(nxl,:,:)  ) * Bx(nxl,:,:)
ay(2,:,:) = dsqrt( Dinv(nxl,:,:)  ) * By(nxl,:,:)
az(2,:,:) = dsqrt( Dinv(nxl,:,:)  ) * Bz(nxl,:,:)
segno(2,:,:) = sign(dble(1.0),Bx(nxl,:,:))
 c_s(2,:,:) = dsqrt(gam * T_tot(nxl,:,:))

! ********** ix = 1, nxl ****************************************************

ax_2 = ax**2
ay_2 = ay**2
az_2 = az**2

a_inv = 1.0d0 / dsqrt( ay_2 + az_2)

c_s_2 = c_s**2
f = c_s_2 + ax_2 + ay_2 + az_2

! parte modificata
s = 0.5d0 * f * ( 1.0d0 - dsqrt( 1.0d0 - 4.0d0 * c_s_2 * ax_2 / f**2 ) )
do iz = 1, nz
 do iy = 1, ny
  do ix = 1, 2
   IF ( s(ix,iy,iz) .LE. 0.0d0 ) THEN
    s(ix,iy,iz) = 0
   ENDIF
  enddo
 enddo
enddo 
f = - s + f
! fine parte modificata


! parte originale
!s = ( f )**2.0 - 4.0d0 * c_s_2 * ax_2
!s = 0.5d0 * ( f - dsqrt( s ) )
!f = - s + f
! fine parte originale

alpha_1 = dsqrt( ( f - ax_2 ) / ( f - s ) )
alpha_2 = dsqrt( ( f - c_s_2 ) / ( f - s ) )

!Tracciante(1,:,:) = alpha_1(1,:,:)
!Tracciante(2,:,:) = alpha_1(2,:,:)
!Tracciante(3,:,:) = alpha_2(1,:,:)
!Tracciante(4,:,:) = alpha_2(2,:,:)
!Tracciante(5,:,:) = f(1,:,:)
!Tracciante(6,:,:) = f(2,:,:)
!Tracciante(7,:,:) = s(1,:,:)
!Tracciante(8,:,:) = s(2,:,:)
!call outx
!call syncronize

f = dsqrt( f )
s = dsqrt( s )

f_inv = 1.0 / f
 c_s_inv =1.0d0 / c_s


deallocate(ax_2)
deallocate(ay_2)
deallocate(az_2)
deallocate(c_s_2)

!write(*,*), ax(1,1,1), ay(1,1,1), az(1,1,1), f(1,1,1), s(1,1,1)
!stop


! ***********       DEFINIRE L_VARI ********************************** 


 call derx_BC(By,grad)

L_s_m = ay * grad
L_a_m = - az * grad


call derx_BC(Bz,grad)

L_s_m = L_s_m + az * grad
L_f_m = L_s_m

L_s_m(1,:,:) = L_s_m(1,:,:) * dsqrt(Dinv(1,:,:)) * f_inv(1,:,:) * c_s(1,:,:)
L_s_m(2,:,:) = L_s_m(2,:,:) * dsqrt(Dinv(nxl,:,:)) * f_inv(2,:,:) * c_s(2,:,:)

L_s_p = L_s_m

L_a_m = L_a_m + ay * grad

L_a_m(1,:,:)= L_a_m(1,:,:) * dsqrt(Dinv(1,:,:)) * segno(1,:,:)
L_a_m(2,:,:)= L_a_m(2,:,:) * dsqrt(Dinv(nxl,:,:)) * segno(2,:,:)

L_a_p = - L_a_m

L_f_m(1,:,:) = L_f_m(1,:,:) * dsqrt(Dinv(1,:,:))
L_f_m(2,:,:) = L_f_m(2,:,:) * dsqrt(Dinv(nxl,:,:))

L_f_p = L_f_m


call derx_BC(Uy,grad)

L_s_m = L_s_m + segno * ay * grad
L_a_m = L_a_m - az * grad
L_f_m = L_f_m + ax * ay * f_inv * grad


call derx_BC(Uz,grad)

L_s_m = L_s_m + segno * az * grad
L_s_p = - L_s_m + 2.0d0 * L_s_p

L_s_m = L_s_m * alpha_1 * a_inv
L_s_p = L_s_p * alpha_1 * a_inv


L_a_m = L_a_m + ay * grad
L_a_p = L_a_m + 2.0d0 * L_a_p

L_a_m(1,:,:) = L_a_m(1,:,:) * a_inv(1,:,:) * (Ux(1,:,:) - segno(1,:,:) * ax(1,:,:))
L_a_m(2,:,:) = L_a_m(2,:,:) * a_inv(2,:,:) * (Ux(nxl,:,:) - segno(2,:,:) * ax(2,:,:))

L_a_p(1,:,:) = L_a_p(1,:,:) * a_inv(1,:,:) * (Ux(1,:,:) + segno(1,:,:) * ax(1,:,:))
L_a_p(2,:,:) = L_a_p(2,:,:) * a_inv(2,:,:) * (Ux(nxl,:,:) + segno(2,:,:) * ax(2,:,:))

L_f_m = L_f_m + ax * az * f_inv * grad
L_f_p = - L_f_m + 2.0d0 * L_f_p

L_f_m = L_f_m * alpha_2 * a_inv
L_f_p = L_f_p * alpha_2 * a_inv


call derx_BC(Ux,grad)

L_s_m = L_s_m + alpha_2 * segno * ax * f_inv * grad
L_s_p = L_s_p - alpha_2 * segno * ax * f_inv * grad

L_f_m = L_f_m - alpha_1 * grad
L_f_p = L_f_p + alpha_1 * grad


 call derx_BC(T_tot,grad)

L_s_m = L_s_m - alpha_2 * c_s_inv * grad
L_s_p = L_s_p - alpha_2 * c_s_inv * grad

L_f_m = L_f_m + alpha_1 * f_inv * grad
L_f_p = L_f_p + alpha_1 * f_inv * grad

L_0(1,:,:) = c_s(1,:,:)/(T_tot(1,:,:) * gam) * grad(1,:,:)
L_0(2,:,:) = c_s(2,:,:)/(T_tot(nxl,:,:) * gam) * grad(2,:,:)


 call derx_BC(Den,grad)

L_s_m(1,:,:) = L_s_m(1,:,:) - alpha_2(1,:,:) * c_s_inv(1,:,:) & 
                  * T_tot(1,:,:) * Dinv(1,:,:) * grad(1,:,:)
L_s_m(1,:,:) = L_s_m(1,:,:) * (Ux(1,:,:) - s(1,:,:))

L_s_p(1,:,:) = L_s_p(1,:,:) - alpha_2(1,:,:) * c_s_inv(1,:,:) & 
                  * T_tot(1,:,:) * Dinv(1,:,:) * grad(1,:,:)
L_s_p(1,:,:) = L_s_p(1,:,:) * (Ux(1,:,:) + s(1,:,:))




L_s_m(2,:,:) = L_s_m(2,:,:) - alpha_2(2,:,:) * c_s_inv(2,:,:) & 
                  * T_tot(nxl,:,:) * Dinv(nxl,:,:) * grad(2,:,:)
L_s_m(2,:,:) = L_s_m(2,:,:) * (Ux(nxl,:,:) - s(2,:,:))

L_s_p(2,:,:) = L_s_p(2,:,:) - alpha_2(2,:,:) * c_s_inv(2,:,:) & 
                  * T_tot(nxl,:,:) * Dinv(nxl,:,:) * grad(2,:,:)
L_s_p(2,:,:) = L_s_p(2,:,:) * (Ux(nxl,:,:) + s(2,:,:))




L_f_m(1,:,:) = L_f_m(1,:,:) + alpha_1(1,:,:) * f_inv(1,:,:) & 
        * T_tot(1,:,:) * Dinv(1,:,:) * grad(1,:,:)
L_f_m(1,:,:) = L_f_m(1,:,:) * (Ux(1,:,:) - f(1,:,:))

L_f_p(1,:,:) = L_f_p(1,:,:) + alpha_1(1,:,:) * f_inv(1,:,:) & 
        * T_tot(1,:,:) * Dinv(1,:,:) * grad(1,:,:)
L_f_p(1,:,:) = L_f_p(1,:,:) * (Ux(1,:,:) + f(1,:,:))




L_f_m(2,:,:) = L_f_m(2,:,:) + alpha_1(2,:,:) * f_inv(2,:,:) & 
        * T_tot(nxl,:,:) * Dinv(nxl,:,:) * grad(2,:,:)
L_f_m(2,:,:) = L_f_m(2,:,:) * (Ux(nxl,:,:) - f(2,:,:))

L_f_p(2,:,:) = L_f_p(2,:,:) + alpha_1(2,:,:) * f_inv(2,:,:) & 
        * T_tot(nxl,:,:) * Dinv(nxl,:,:) * grad(2,:,:)
L_f_p(2,:,:) = L_f_p(2,:,:) * (Ux(nxl,:,:) + f(2,:,:))




L_0(1,:,:) = L_0(1,:,:) + c_s(1,:,:) * gm1 / gam * Dinv(1,:,:) * grad(1,:,:)
L_0(1,:,:) = L_0(1,:,:) * Ux(1,:,:)


L_0(2,:,:) = L_0(2,:,:) + c_s(2,:,:) * gm1 / gam * Dinv(nxl,:,:) * grad(2,:,:)
L_0(2,:,:) = L_0(2,:,:) * Ux(nxl,:,:)



deallocate(grad)
deallocate(T_tot)


! ***********   CONDIZIONI AL BORDO    ********************************** 
! in questo modo tutto e' calcolato con i punti interni.
! selezioniamo ora le onde uscenti, date da (Ux +/- v_varie) uscente dal 
! dominio...imponiamo quindi il loro valore compatibilemente con
! 1) equilibrio (parte indipendente dal tempo)
! 2) perturbazione entrante (parte dipendente dal tempo).

! nel nostro caso (equilibrio senza inflow e bordi perf. non-riflettenti)
!  tutte caratteristiche entranti = 0.0 !!!!! 

do iz = 1, nz 
 do iy = 1, ny 

! ************* Bordo sx ********************************************

IF ((Ux(1,iy,iz) + f(1,iy,iz)) .GT. 0.0)  then
L_f_p(1,iy,iz) = 0.0
ENDIF

IF ((Ux(1,iy,iz) + segno(1,iy,iz) * ax(1,iy,iz)) .GT. 0.0)  then
L_a_p(1,iy,iz) = 0.0
ENDIF

IF ((Ux(1,iy,iz) + s(1,iy,iz)) .GT. 0.0)  then
L_s_p(1,iy,iz) = 0.0
ENDIF

IF ((Ux(1,iy,iz) - f(1,iy,iz)) .GT. 0.0)  then
L_f_m(1,iy,iz) = 0.0
ENDIF

IF ((Ux(1,iy,iz) - segno(1,iy,iz) * ax(1,iy,iz)) .GT. 0.0)  then
L_a_m(1,iy,iz) = 0.0
ENDIF

IF ((Ux(1,iy,iz) - s(1,iy,iz)) .GT. 0.0)  then
L_s_m(1,iy,iz) = 0.0
ENDIF

IF ((Ux(1,iy,iz)) .GT. 0.0d0)  then
L_0(1,iy,iz) = 0.0d0
ENDIF


! ************ Bordo dx *************************************************

IF ((Ux(nxl,iy,iz) + f(2,iy,iz)) .LT. 0.0)  then
L_f_p(2,iy,iz) = 0.0
ENDIF

IF ((Ux(nxl,iy,iz) + segno(2,iy,iz) * ax(2,iy,iz)) .LT. 0.0)  then
L_a_p(2,iy,iz) = 0.0
ENDIF

IF ((Ux(nxl,iy,iz) + s(2,iy,iz)) .LT. 0.0)  then
L_s_p(2,iy,iz) = 0.0
ENDIF

IF ((Ux(nxl,iy,iz) - f(2,iy,iz)) .LT. 0.0)  then
L_f_m(2,iy,iz) = 0.0
ENDIF

IF ((Ux(nxl,iy,iz) - segno(2,iy,iz) * ax(2,iy,iz)) .LT. 0.0)  then
L_a_m(2,iy,iz) = 0.0
ENDIF

IF ((Ux(nxl,iy,iz) - s(2,iy,iz)) .LT. 0.0)  then
L_s_m(2,iy,iz) = 0.0
ENDIF

IF ((Ux(nxl,iy,iz)) .LT. 0.0d0)  then
L_0(2,iy,iz) = 0.0d0
ENDIF


 enddo
enddo

! N.B. non abbiamo calcolato le caratteristiche uscenti!
! sono poste = 0. ...nulla entra nel dominio, allo stesso modo abbiamo posto
! = 0. anche le ''uscenti'' ma che presentino autovalote (u_x +/- vel) che 
! in realta' sia entrante.
! Questo e' corretto per un equilibrio che non contempli nulla di entrante,
! e condizioni al bordo perfettamente riflettenti.
! altrimenti: 1) parte costante nel tempo per garantire l'equilibrio (entrata)
! 2) parte dipendente dal tempo della perturbazione entrante. 




! TEST

!************** Onda  entrante, a ***************************

!IF (tempo .LT. 30.0)  then
!L_a_m(2,:,:) = 0.005 * dcos(0.20944 * tempo)
!else
!L_a_m(2,:,:) = 0.0
!ENDIF


!************** Onda  entrante, f ***************************

!IF (tempo .LT. 19.1083)  then
!L_f_m(2,:,:) = 0.05 * dcos(0.32882 * tempo)
!else
!L_f_m(2,:,:) = 0.0
!ENDIF


!************** Onda  entrante, s ***************************

!IF (tempo .LT. 53.5714)  then
!L_s_m(2,:,:) = 0.005 * dcos(0.117286 * tempo)
!else
!L_s_m(2,:,:) = 0.0
!ENDIF





end subroutine
