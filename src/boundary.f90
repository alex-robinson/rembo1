
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  b o u n d a r y
  ! Purpose    :  Definition of the surface temperature (must be less
  !               than 0 deg C!) and of the accumulation-ablation
  !               function, updating mask.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine boundary(n_step,temp_s,accum_now, as_perp, abl_perp,   &
                      abl_o, m_shallow, zs, zb, zb0, z_sl,  &
                      maske, m_lakes, z_lakes, maske_migr)

    use emb_global
    
    implicit none
    
    integer, parameter :: nx = nxs, ny = nys
    real (8), parameter :: dx = dxs
      
    integer :: n_step, i, j
    
    integer, dimension(ny,nx) :: m_shallow, maske, m_lakes, maske_migr
    integer, dimension(nys,nxs) :: maske_neu
    
    real (8) :: z_sl, z_lakes, rho_ratio
    real (8), dimension(ny,nx) :: temp_s, accum_now, abl_perp, as_perp,  &
                                  abl_o, zs, zb, zb0
    real (8), dimension(ny,nx) :: m2
    
    ! Only perform boundary calcs if it is the right timestep
    ! (normally on a timestep when sicopolis was called)
    if(nstep.eq.nstep/dtime_bndry*dtime_bndry) then
    
      ! Temporary, put maske into double var
      m2 = dble(maske)
      
      ! Set emb global nstep to driver's n_step
      nstep = n_step
      
      !  ------ Update of the mask according to the sea level
      !         *** Currently some parts removed implemented!

      ! update the mask      
      call mask_updater(z_sl,zs,zb,zb0,m2)
      
      ! melting of inland ice through ocean heat exchange
      call ocean_melt(m2,zb0,abl_o,m_shallow) 
      
      !write(*,*) "abl_o min/max: ",minval(abl_o),maxval(abl_o)
      !write(*,*) "abl_perp min/max: ",minval(abl_perp),maxval(abl_perp)
      !write(*,*) "as_perp min/max: ",minval(as_perp),maxval(as_perp)
      !stop
      
! Exchange from Climber-2 to Sicopolis

!#if CLIMBER_ON==1
!     call climber_to_sico(maske_migr, temp_s, as_perp, abl_perp, abl_o)
!#endif

      ! Temporary, update maske with m2
      maske = int(m2)
      
    end if
     
    return
    
  end subroutine boundary
  
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Function :  m a s k _ u p d a t e r
!   Purpose  :  Update the topography mask due to the sea level (new).
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine mask_updater(z_sl,zs,zb,zb0,m2)
    
    use emb_global
    
    implicit none
    
    integer, parameter :: nx = nxs, ny = nys
    real (8), dimension(ny,nx) :: m2, zs, zb, zb0
    
    integer :: i, j
    real (8) :: z_sl, dz_thr
    real (8) :: rho_ratio
        
    dz_thr = 10.d0
    rho_ratio = rho_sw / rho_i
    
    do i = 1, nx
      do j = 1, ny
        
        ! -------- Previously ice-free land point or sea point
        if ( (m2(j,i) .eq. 1.d0) .or. (m2(j,i) .eq. 2.d0) ) then
          
          if (zb(j,i) .gt. z_sl) then
            m2(j,i) = 1.d0           ! now ice-free land point
          else
            m2(j,i) = 2.d0           ! now sea point
          end if
        
        ! Previously glaciated point
        else if (m2(j,i) .eq. 0.d0) then
          
          ! Will stay glaciated point unless:
          if ( zs(j,i) .lt. zb(j,i) + rho_ratio*(z_sl-zb(j,i)) ) then                           
              m2(j,i) = margin_value
          end if
         
        ! Previously floating-ice point
        else if (m2(j,i) .eq. 3.d0) then

          if (zb(j,i) .gt. z_sl) then
            m2(j,i) = 0.d0                 ! now glaciated land point
          else
            
            if ( zs(j,i) .lt. dz_thr + zb(j,i) + rho_ratio*(z_sl-zb(j,i))) then                           
              m2(j,i) = margin_value    ! floating-ice stays floating ice
            else
              m2(j,i) = 0.d0               ! returns to being grounded ice
            end if
          
          end if
        
        end if
      
      end do
    end do

    return
  
  end subroutine mask_updater

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Function :  o c e a n _ m e l t
! Author   :  Reinhard Calov
! Purpose  :  Calculation of inland ice melt through heat exchange with
!             the ocean
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocean_melt(m2, zb0, abl_o, m_shallow)
  
    use emb_global
    
    implicit none
    
    integer, parameter :: nx = nxs, ny = nys
    
    real (8), dimension(ny,nx) :: m2, zb0, abl_o
    integer, dimension(ny,nx) :: m_shallow
    
    integer i, j
    real (8) :: abl_o_500, abl_o_0, abl_err, abl_diff
    
    double precision, parameter :: depth = -500.d0
    
    ! Calculation of melting through the heat of the ocean

    abl_o_0   = ablo_min / sec_year        ! m/a --> m/s !
    abl_o_500 = ablo_max / sec_year        ! m/a --> m/s !
    abl_err   = 50.d0  / sec_year        ! m/a --> m/s !
    abl_diff  = abl_o_500 - abl_o_0
    
    if (abl_diff .lt. 0.d0) then
      write(*,*) "ocean_melt: deeper melt value should be greater"
      write(*,*) "            than shallow value. Stopping."
      stop
    end if
     
    ! If grid point is in the ocean, set heating to desired value
    ! *higher melting value as ocean gets deeper, so there are no 
    !  large discontinuities in melting
    where(m2 .ge. 2.d0)
      abl_o = abl_o_0 + abl_diff * max(depth,min(0.d0,zb0)) / depth
    elsewhere
      abl_o = 0.d0
    end where
    
    ! Make sure border points have very higher ablation so no ice goes there!
    abl_o(1:2,:)     = abl_err
    abl_o(ny-1:ny,:) = abl_err
    abl_o(:,1:2)     = abl_err
    abl_o(:,nx-1:nx) = abl_err

    return
  
  end subroutine ocean_melt
  
  
  
