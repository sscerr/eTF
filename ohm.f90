subroutine ohm

!****************************************************************!
! COMPUTE: generalized Ohm's law                                 !
!          [2D version, no electron-inertia terms]               !
! -------------------------------------------------------------- ! 
! original version:             F. Califano, 2006                !
! MPI parallel version:         M. Faganello/F. Valentini, 2008  !                  
! Anisotropic/FLR-MHD version:  S. S. Cerri, 2011                !
!****************************************************************!

use parameter_mod
use box_mod
use deriv_mod
use fields_UJ_mod
use fields_DP_mod
use fields_EB_mod
use dom_distr_mod


IMPLICIT NONE

integer :: ix, iy, iz

REAL(dp),ALLOCATABLE :: zy(:), zx(:), zz(:)
REAL(dp),ALLOCATABLE :: AA(:,:,:), AB(:,:,:)
REAL(dp),ALLOCATABLE :: At(:,:,:)

REAL(dp),ALLOCATABLE :: bbx(:,:,:), bby(:,:,:), bbz(:,:,:)
REAL(dp),ALLOCATABLE :: Pdiag(:,:,:), Pmix(:,:,:)



  !*******************************************!
  !     Computing:  E = - u_s x B + eta*J     !
  !                                           !
  !  [ ideal-MHD/or/Hall-MHD + resistivity ]  !
  !*******************************************!
  
  if (hall_on .eq. 1) then
    !case: Hall-MHD + resistivity
    ! E = - u_e x B + eta*J 
    Ex = Uez * By - Uey * Bz + eta*Jx
    Ey = Uex * Bz - Uez * Bx + eta*Jy
    Ez = Uey * Bx - Uex * By + eta*Jz
  else
    !case: ideal-MHD + resistivity
    ! E = - u_i x B + eta*J
    Ex = Uz * By - Uy * Bz + eta*Jx  !Ux = Uix when m_e = 0
    Ey = Ux * Bz - Uz * Bx + eta*Jy  !Uy = Uiy when m_e = 0
    Ez = Uy * Bx - Ux * By + eta*Jz  !Uz = Uiz when m_e = 0
  endif


  allocate(AA(nxl, ny, nz)) 


  if (divPe_on .eq. 1) then

    !***********************************************!
    !      Computing: - div(Pe)/n contribution      !
    !                                               !
    !  [ Pe = pe_perp*I + (pe_para - pe_perp)*bb ]  !
    !***********************************************!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                                           !!!
    !!!  [ NOTICE: 3D version not yet released ]  !!!
    !!!                                           !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(zy(ny))
    allocate(zx(nx))
    allocate(zz(nz))
    allocate(AB(nxl, ny, nz))
    allocate(At(nx, nyt, nz))
 
    allocate(bbx(nxl, ny, nz))
    allocate(bby(nxl, ny, nz))
    allocate(bbz(nxl, ny, nz))
    allocate(Pdiag(nxl, ny, nz))
    allocate(Pmix(nxl, ny, nz))
   
    !define useful quantities
    AA = dsqrt( Bx**2 + By**2 + Bz**2 )
    AA = 1.0d0 / AA
    bbx = Bx * AA
    bby = By * AA
    bbz = Bz * AA
   
    ! further finite-m_e contributions: 
    ! not available in this release 
    Pdiag = pe_perp ! + ... 
    Pmix = pe_para - pe_perp ! + ... 
 
    !-------------------------------------------!
    ! x-component:                              !
    !                                           !
    !   - {div[Pe]}_x  =  - d[Pe_xx]/dx         !
    !                     - d[Pe_yx]/dy         !
    !                     - d[Pe_zx]/dz         !
    !                                           !
    ! (derivatives w.r.t z only in 3D version)  !
    !-------------------------------------------!
   
    AA = Pdiag + Pmix*bbx*bbx 
   
    call traspdist( AA, At, 1 )
   
    do iz = 1, nz
     do iy = 1, nyt
      call  derx_1(At(:,iy,iz), zx)
      At(:,iy,iz) = -zx
     enddo
    enddo
   
    call traspdist( AB, At, -1 )
   
    AA = Pmix*bbx*bby 
   
    do iz = 1, nz
     do ix = 1, nxl
      call  dery_1(AA(ix,:,iz), zy)
      AB(ix,:,iz) = AB(ix,:,iz) - zy
     enddo
    enddo
   
   
    Ex = Ex + Dinv * AB
   
    !-------------------------------------------!
    ! y-component:                              !
    !                                           !
    !   - {div[Pe]}_y  =  - d[Pe_xy]/dx         !
    !                     - d[Pe_yy]/dy         !
    !                     - d[Pe_zy]/dz         !
    !                                           !
    ! (derivatives w.r.t z only in 3D version)  !
    !-------------------------------------------!
  
    ! AA = Pmix*bby*bbx  -> same as above  
   
    call traspdist( AA, At, 1 )
   
    do iz = 1, nz
      do iy = 1, nyt
        call  derx_1(At(:,iy,iz), zx)
        At(:,iy,iz) = -zx
      enddo
    enddo
   
    call traspdist( AB, At, -1 )
  
    AA = Pdiag + Pmix*bby*bby 
  
    do iz = 1, nz
      do ix = 1, nxl
        call  dery_1(AA(ix,:,iz), zy)
        AB(ix,:,iz) = AB(ix,:,iz) - zy
      enddo
    enddo
    
    Ey = Ey + Dinv * AB
   
    !-------------------------------------------!
    ! z-component:                              !
    !                                           !
    !   - {div[Pe]}_x  =  - d[Pe_xz]/dx         !
    !                     - d[Pe_yz]/dy         !
    !                     - d[Pe_zz]/dz         !
    !                                           !
    ! (derivatives w.r.t z only in 3D version)  !
    !-------------------------------------------!

    AA = Pmix*bbx*bbz 
   
    call traspdist( AA, At, 1 )
   
    do iz = 1, nz
      do iy = 1, nyt
        call  derx_1(At(:,iy,iz), zx)
        At(:,iy,iz) = -zx
      enddo
    enddo
   
    call traspdist( AB, At, -1 )
    
    AA = Pmix*bby*bbz
   
    do iz = 1, nz
      do ix = 1, nxl
        call  dery_1(AA(ix,:,iz), zy)
        AB(ix,:,iz) = AB(ix,:,iz) - zy
      enddo
    enddo
   
    Ez = Ez + Dinv * AB
   
   
    deallocate(bbx)
    deallocate(bby)
    deallocate(bbz)
    
    deallocate(Pdiag)
    deallocate(Pmix)
   
    deallocate(AB)
    deallocate(At)
    deallocate(zy)
    deallocate(zx)
    deallocate(zz)
 
  endif

  if ( (hall_on .eq. 1) .or. (divPe_on .eq. 1)) then

    !---------------------------------------------!
    ! smoothing -> ideal MHD at the x-boundaries  !
    !---------------------------------------------!

    AA =  Uz * By - Uy * Bz
    do ix = 1, nxl
      Ex(ix,:,:) = lambda(ixlg + ix - 1) * Ex(ix,:,:) + &
                   ( 1.0d0 - lambda(ixlg + ix - 1) ) * AA(ix,:,:)
    enddo

    AA =  Ux * Bz - Uz * Bx
    do ix = 1, nxl
      Ey(ix,:,:) = lambda(ixlg + ix - 1) * Ey(ix,:,:) + &
                   ( 1.0d0 - lambda(ixlg + ix - 1) ) * AA(ix,:,:)
    enddo

    AA =  Uy * Bx - Ux * By
    do ix = 1, nxl
      Ez(ix,:,:) = lambda(ixlg + ix - 1) * Ez(ix,:,:) + &
                   ( 1.0d0 - lambda(ixlg + ix - 1) ) * AA(ix,:,:)
    enddo

  endif

  deallocate(AA)

end subroutine

