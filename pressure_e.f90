subroutine  pressure_e(RHS_x,RHS_y)

!------------------------------------------------------!
!       M.Faganello, 2011                              !
!                                                      !
!   Eq. chiusura per gli ioni: calcolo di RHS per      !
!                                                      !
!          p_para e p_perp                             !
!                                                      !
! Cerri 2011: banale riadattamento per gli elettroni   !
!                                                      !
!             pe_perp & pe_para                        !
!------------------------------------------------------!



use parameter_mod
use box_mod
use deriv_mod
use fields_UJ_mod
use fields_DP_mod
use fields_EB_mod
use dom_distr_mod
use parallel_mod

IMPLICIT NONE

integer :: ix, iy, iz

REAL(dp),DIMENSION (nxl,ny,nz) :: RHS_x, RHS_y

REAL(dp),ALLOCATABLE :: zx(:), zy(:), zz(:)
REAL(dp),ALLOCATABLE :: AA(:,:,:), bbx(:,:,:), bby(:,:,:), bbz(:,:,:) 
REAL(dp),ALLOCATABLE :: At(:,:,:)

allocate(zx(nx))
allocate(zy(ny))
allocate(zz(nz))
allocate(AA(nxl, ny, nz))
allocate(bbx(nxl, ny, nz))
allocate(bby(nxl, ny, nz))
allocate(bbz(nxl, ny, nz))
allocate(At(nx, nyt, nz))


AA = dsqrt( Bx**2 + By**2 + Bz**2 )
AA = 1.0d0 / AA
bbx = Bx * AA
bby = By * AA
bbz = Bz * AA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calcolo div b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL traspdist( bbx, At, 1 )

do iz = 1, nz 
 do iy = 1, nyt
  call  derx_1( At(:,iy, iz), zx )
  At(:,iy, iz) = zx
 enddo
enddo

 CALL traspdist( RHS_x, At, -1 )

do iz= 1, nz
 do ix = 1, nxl
  call  dery_1(bby(ix,:,iz), zy)
  RHS_x(ix,:,iz) = RHS_x(ix,:,iz) + zy
 enddo
enddo


RHS_x = RHS_x * qe_perp
!RHS_y = - RHS_x                     lo faccio dopo !!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calcolo ( bb div Ue )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL traspdist( Uex, At, 1 )

do iz = 1, nz 
 do iy = 1, nyt
  call  derx_1( At(:,iy, iz), zx )
  At(:,iy, iz) = zx
 enddo
enddo

 CALL traspdist( RHS_y, At, -1 )

do iz= 1, nz
 do ix = 1, nxl
  call  dery_1(Uex(ix,:,iz), zy)
  RHS_y(ix,:,iz) = bbx(ix,:,iz) * RHS_y(ix,:,iz) + bby(ix,:,iz) * zy
 enddo
enddo



AA = RHS_y * bbx


  CALL traspdist( Uey, At, 1 )

do iz = 1, nz 
 do iy = 1, nyt
  call  derx_1( At(:,iy, iz), zx )
  At(:,iy, iz) = zx
 enddo
enddo

 CALL traspdist( RHS_y, At, -1 )

do iz= 1, nz
 do ix = 1, nxl
  call  dery_1(Uey(ix,:,iz), zy)
  RHS_y(ix,:,iz) = bbx(ix,:,iz) * RHS_y(ix,:,iz) + bby(ix,:,iz) * zy
 enddo
enddo


AA = AA + bby * RHS_y

  CALL traspdist( Uez, At, 1 )

do iz = 1, nz 
 do iy = 1, nyt
  call  derx_1( At(:,iy, iz), zx )
  At(:,iy, iz) = zx
 enddo
enddo

 CALL traspdist( RHS_y, At, -1 )

do iz= 1, nz
 do ix = 1, nxl
  call  dery_1(Uez(ix,:,iz), zy)
  RHS_y(ix,:,iz) = bbx(ix,:,iz) * RHS_y(ix,:,iz) + bby(ix,:,iz) * zy
 enddo
enddo



AA = AA + bbz * RHS_y

RHS_y = - RHS_x
RHS_x = 2.0d0 * ( RHS_x - AA * pe_para ) 
RHS_y = RHS_y + AA * pe_perp  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                calcolo div Ue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  CALL traspdist( Uex, At, 1 )

do iz = 1, nz 
 do iy = 1, nyt
  call  derx_1( At(:,iy, iz), zx )
  At(:,iy, iz) = zx
 enddo
enddo

 CALL traspdist( AA, At, -1 )


do iz= 1, nz
 do ix = 1, nxl
  call  dery_1(Uey(ix,:,iz), zy)
  AA(ix,:,iz) = AA(ix,:,iz) + zy
 enddo
enddo


RHS_y = RHS_y - AA * pe_perp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    calcolo div [ U*pe_para + b*qe_para ]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



AA = Uex * pe_para + bbx * qe_para

  CALL traspdist( AA, At, 1 )

do iz = 1, nz 
 do iy = 1, nyt
  call  derx_1( At(:,iy, iz), zx )
  At(:,iy, iz) = zx
 enddo
enddo

 CALL traspdist( AA, At, -1 )

RHS_x = RHS_x - AA


AA =  Uey * pe_para + bby * qe_para

do iz= 1, nz
 do ix = 1, nxl
  call  dery_1(AA(ix,:,iz), zy)
  RHS_x(ix,:,iz) = RHS_x(ix,:,iz) - zy
 enddo
enddo


AA = Uez * pe_para + bbz * qe_para


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    calcolo div [ Ue*pe_perp + b*qe_perp ]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

AA = Uex * pe_perp + bbx * qe_perp

  CALL traspdist( AA, At, 1 )

do iz = 1, nz 
 do iy = 1, nyt
  call  derx_1( At(:,iy, iz), zx )
  At(:,iy, iz) = zx
 enddo
enddo

 CALL traspdist( AA, At, -1 )

RHS_y = RHS_y - AA


AA =  Uey * pe_perp + bby * qe_perp

do iz= 1, nz
 do ix = 1, nxl
  call  dery_1(AA(ix,:,iz), zy)
  RHS_y(ix,:,iz) = RHS_y(ix,:,iz) - zy
 enddo
enddo


AA = Uez * pe_perp + bbz * qe_perp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 BC MHD con caratteristiche lungo x
!
!
!  *************************************************************
!  NOTA:
!  riadattata per elettroni: ho semplicemente usato le quantita'
!  elettroniche al posto di quelle ioniche. Comunque chiedere
!  conferma a matteo che non ci sia altro da cambiare.
!  *************************************************************
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (ibc .EQ. 4) then
!N.B. va sommato  con il meno
!RHS_x=n*dT/dt+T*n*dn/dt
!dove le derivate sono da intendersi parziali

 if ( mpime == 0 ) then 

RHS_x(1,:,:)= +(Den(1,:,:)*(-Te_para(1,:,:)*c_s_inv(1,:,:)* &
             L_0(1,:,:) - Te_para(1,:,:)*gm1* 0.5 * & 
             (alpha_2(1,:,:) * c_s_inv(1,:,:) * &
        (L_s_p(1,:,:) + L_s_m(1,:,:)) - alpha_1(1,:,:) * f_inv(1,:,:)* &
           (L_f_p(1,:,:) + L_f_m(1,:,:)))+Te_para(1,:,:)* &
             ( c_s_inv(1,:,:) * L_0(1,:,:) + 0.5 *& 
             (alpha_2(1,:,:) * c_s_inv(1,:,:)* &
        (L_s_p(1,:,:) + L_s_m(1,:,:)) - alpha_1(1,:,:) * f_inv(1,:,:)* &
             (L_f_p(1,:,:) + L_f_m(1,:,:))))))

AA(1,:,:) = - pe_para(1,:,:) * Uey(1,:,:) * Den(1,:,:)**gm1 

do iz = 1, nz
  call  dery_1(AA(1,:,iz), zy)
  RHS_x(1,:,iz) = RHS_x(1,:,iz) + Den(1,:,iz)**(-gm1) * zy
enddo

AA(1,:,:) = Den(1,:,:) * Uey(1,:,:) 

do iz = 1, nz
  call  dery_1(AA(1,:,iz), zy)
  RHS_x(1,:,iz) = RHS_x(1,:,iz) + gm1 * Te_para(1,:,iz) * zy
enddo


RHS_y(1,:,:) = pe_perp(1,:,:) / pe_para(1,:,:) * RHS_x(1,:,:)



 endif
 if ( mpime == nprow -1 ) then 

RHS_x(nxl,:,:)= +(Den(nxl,:,:)*(-Te_para(nxl,:,:)*c_s_inv(2,:,:)* &
            L_0(2,:,:) - Te_para(nxl,:,:)*gm1* 0.5 * & 
             (alpha_2(2,:,:) * c_s_inv(2,:,:) * &
        (L_s_p(2,:,:) + L_s_m(2,:,:)) - alpha_1(2,:,:) * f_inv(2,:,:)* &
             (L_f_p(2,:,:) + L_f_m(2,:,:)))+ Te_para(nxl,:,:)* &
             ( c_s_inv(2,:,:) * L_0(2,:,:) + 0.5 * &
             (alpha_2(2,:,:) * c_s_inv(2,:,:) * &
        (L_s_p(2,:,:) + L_s_m(2,:,:)) - alpha_1(2,:,:) * f_inv(2,:,:) * &
             (L_f_p(2,:,:) + L_f_m(2,:,:))))))


AA(nxl,:,:) = - pe_para(nxl,:,:) * Uey(nxl,:,:) * Den(nxl,:,:)**gm1 

do iz = 1, nz
  call  dery_1(AA(nxl,:,iz), zy)
  RHS_x(nxl,:,iz) = RHS_x(nxl,:,iz) + Den(nxl,:,iz)**(-gm1) * zy
enddo

AA(nxl,:,:) = Den(nxl,:,:) * Uey(nxl,:,:) 

do iz = 1, nz
  call  dery_1(AA(nxl,:,iz), zy)
  RHS_x(nxl,:,iz) = RHS_x(nxl,:,iz) + gm1 * Te_para(nxl,:,iz) * zy
enddo


RHS_y(nxl,:,:) = pe_perp(nxl,:,:) / pe_para(nxl,:,:) * RHS_x(nxl,:,:)

 endif
ENDIF

deallocate(zx)
deallocate(zy)
deallocate(zz)
deallocate(AA)
deallocate(At)
deallocate(bbx)
deallocate(bby)
deallocate(bbz)


end subroutine
