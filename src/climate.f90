module rembo_sclimate

  implicit none 

  type rembo_param_class
    double precision :: dtime_emb
    double precision :: dtime_smb
    integer :: nx
    integer :: ny 
    real(8) :: dx
  end type

  type rembo_class  ! ij indices, times
    type(rembo_param_class) :: par 
    double precision, allocatable :: pr(:,:) 
    double precision, allocatable :: sf(:,:) 
    double precision, allocatable :: smb(:,:) 
    double precision, allocatable :: melt(:,:) 
    double precision, allocatable :: runoff(:,:)
    double precision, allocatable :: refrz(:,:) 
    double precision, allocatable :: H_snow(:,:)
    double precision, allocatable :: T_srf(:,:)  
    double precision, allocatable :: T_ann(:,:) 
    double precision, allocatable :: T_pann(:,:) 
    double precision, allocatable :: T_jja(:,:) 
    double precision, allocatable :: T_djf(:,:) 
    double precision, allocatable :: pdds(:,:)

    integer, allocatable :: mask_relax(:,:)

    double precision :: time_emb, time_smb
    double precision :: dt_emb,   dt_smb
    
    double precision :: time_dto_clim
    double precision :: time_dto_clim2D
  end type 

  type(rembo_class) :: rembo_ann 

  type precip_correction_type
    real(8), allocatable :: pp_rembo(:,:,:)
    real(8), allocatable :: pp_ref(:,:,:)
    real(8), allocatable :: dpp_corr(:,:,:)
  end type

  type(precip_correction_type) :: ppcorr0

  private
  public :: rembo_init 
  public :: rembo_equilibrate
  public :: rembo_update 
  public :: rembo_class
  public :: rembo_ann
  public :: rembo_restart_write

contains 

subroutine rembo_equilibrate(time,z_srf,H_ice,z_sl,mask_relax)

    implicit none

    real(8), intent(IN) :: time                         ! Current external driver time
    real(8), intent(IN), optional :: z_srf(:,:)         ! Input from ice-sheet model or data
    real(8), intent(IN), optional :: H_ice(:,:)         ! Input from ice-sheet model or data
    real(8), intent(IN), optional :: z_sl(:,:)          ! Input from ice-sheet model or data
    integer, intent(IN), optional :: mask_relax(:,:)    ! Relaxation mask

    ! Local variables
    integer :: nx, ny, n, ntot, k, nm
    real(8) :: dx
    real(8) :: time_now
    real(8) :: time_init
    real(8) :: time_end
    real(8), parameter :: time_ins  = 0.0   ! years BP
    real(8), parameter :: dT_summer = 0.0

    nx = rembo_ann%par%nx
    ny = rembo_ann%par%ny
    dx = rembo_ann%par%dx
    
    nm = 12 

    ntot = 100          ! Run for 100 years
    time_init = time
    time_end  = time_init + real(ntot,8)     

    do n = 1, ntot
      time_now = time_init + real(n-1,8)
      call rembo_update(time_now,time_ins,dT_summer,z_srf,H_ice,z_sl,mask_relax)
    end do

    ! Reset rembo time back to initial time
    rembo_ann%time_emb = time 
    rembo_ann%time_smb = time

    ! Load reference precipitation dataset

    ! TO DO 

    ! Calculate precipitation correction factor

    do k = 1, nm
      call rembo_calc_precip_corr(ppcorr0%dpp_corr(:,:,k),ppcorr0%pp_rembo(:,:,k), &
                                        ppcorr0%pp_ref(:,:,k),dx,sigma=100d3,max_corr=0.5d0)
    end do

    return

end subroutine rembo_equilibrate

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Subroutine :  s c l i m a t e
! Author   :  Alex Robinson
! Purpose  :  Wrapper subroutine to choose between emb module
!             and other methods of generating temperature and
!             accumulation distributions for forcing of sicopolis
!             Whatever method is chosen, it will return daily arrays 
!             over 1 year of accum and temperature
!             *All data sets should be on sicopolis grid size...
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  subroutine rembo_update(time,time_ins,dT_summer,z_srf,H_ice,z_sl,mask_relax,dT_mon,co2)
  
    use emb_global
    use emb_functions
    use rembo_main
  
    implicit none
    
    integer, parameter :: nx=nxs, ny=nys
    real(8), parameter :: dx=dxs
    
    real(8), intent(IN) :: time                         ! Current external driver time
    real(8), intent(IN) :: time_ins                     ! Time to use for insolation calculations (usually years BP)
    real(8), intent(IN) :: dT_summer                    ! Summer temp anomaly [K]
    real(8), intent(IN), optional :: z_srf(nx,ny)       ! Input from ice-sheet model or data
    real(8), intent(IN), optional :: H_ice(nx,ny)       ! Input from ice-sheet model or data
    real(8), intent(IN), optional :: z_sl(nx,ny)        ! Input from ice-sheet model or data
    integer, intent(IN), optional :: mask_relax(nx,ny)  ! Relaxation mask
    real(8), intent(IN), optional :: dT_mon(12)         ! Monthly temperature anomalies
    real(8), intent(IN), optional :: co2                ! Global atmospheric CO2 concentation

    real (8), dimension(ny,nx) :: m2, zs, lats, lons, aco2
    real (8), dimension(ny,nx) :: ZZ, tma, tmj, ampl
    real (8), dimension(ny,nx) :: restart_tt, restart_pp
    real (8), dimension(ny,nx) :: restart_snowh
    real (8), dimension(ny,nx) :: restart_dh, restart_ap
    real (8), dimension(ny,nx) :: mask
    integer,  dimension(ny,nx) :: mrelax
    real (8) :: wt, pscalar, sfac
    real (8) :: co2_now 

    character(len=10) :: c_dxs, dattype
    character(len=256) :: fnm, restart_folder
    character(len=15) :: c_yr
    
    integer :: i, j, k, qq, n, n_now
    real (8) :: dtime1, dtime2
    real (8) :: yearnow, yearnow1, get_year
    real (8) :: tmp_noise
    
    real (8) :: bndtmp
    
    real(8) :: T_warming        ! Summer temperature anomaly
    real(8) :: T_anomaly(nk)    ! Daily temperature anomalies

    type(choices) :: now

    logical :: init_summary2_file

    ! Safety check
    if ( present(z_srf) .and. (.not. present(H_ice)) .or. &
         (.not. present(z_srf)) .and. present(H_ice) ) then 
      write(*,*) "Error: rembo_update:: Both z_srf and H_ice must be provided as arguments."
      stop
    end if

    ! Update boundary climate forcing in program (stored in global variables Tanomaly and T_warming)
    if (present(dT_mon)) then
      ! Get daily temperature anomalies from input monthly anomalies
      call rembo_update_boundary_forcing_monthly(T_anomaly,T_warming,dT_mon,day_month,day_year)
    else
      ! Get daily temperature anomalies from input T_summer and T_wintfac values
      call rembo_update_boundary_forcing_sin(T_anomaly,T_warming,dT_summer,dT_summer*T_wintfac,day_year)
    end if

    if (present(co2)) then
      co2_now = co2
    else
      co2_now = calc_co2(T_warming)
    end if

if (.FALSE.) then
  ! ajr: check anomalies...
    write(*,*) "co2_now:   ", co2_now
    write(*,*) "T_warming: ", T_warming
    write(*,*) "T_anomaly: "
    do k = 1, nk
      write(*,*) T_anomaly(k)
    end do

    stop 
end if

    ! Kill program as needed (for fracs and other simulations)
    if ( kill .eq. 1 ) then
      write(*,*) "Kill switch activated, program finished."
      stop
    end if
  
    call cpu_time(timer%climate)           ! get current time in seconds   
    
    ! ##### Get time information #####

    ! yearnow needed for boundary forcing and sinsol2d
    yearnow = time_ins

    ! Get the current year, as relevant to output
    yearnow1 = time

    ! Adjust time step for calling SMB module
    ! HACK for long transient future simulations (ajr: 15.05.2012)
    !if ( yearnow .gt. 1e3 ) dtime_smb = 10
    
    ! Determine which functions should be called now
    ! Check which operations should be performed for this timestep
    now%clim = .FALSE.
    now%smb  = .FALSE.
    !if( nstep.eq.nstep/dtime_emb*dtime_emb ) now%clim = .TRUE.
    !if( nstep.eq.nstep/dtime_smb*dtime_smb ) now%smb  = .TRUE.

    if (time - rembo_ann%time_emb .ge. rembo_ann%par%dtime_emb) then 
      ! Updated emb 
      now%clim = .TRUE. 
      rembo_ann%dt_emb   = time - rembo_ann%time_emb
      rembo_ann%time_emb = time 
      write(*,"(a,4f12.3)") "emb: ", time, rembo_ann%time_emb, rembo_ann%dt_emb, rembo_ann%par%dtime_emb
    end if 

    if (time - rembo_ann%time_smb .ge. rembo_ann%par%dtime_smb) then 
      ! Updated smb 
      now%smb = .TRUE. 
      rembo_ann%dt_smb   = time - rembo_ann%time_smb
      rembo_ann%time_smb = time 
      write(*,"(a,4f12.3)") "smb: ", time, rembo_ann%time_smb, rembo_ann%dt_smb, rembo_ann%par%dtime_smb
    end if 
    
    ! If clim or smb is running, update forcing and topo
    if ( now%clim .or. now%smb ) then
      
      ! Update fields from exchange (if necessary)
      !call vars2clim(zs,m2,transT%dVdt)
      
      ! ================================================
      ! remboyelmo 

      if (present(z_srf) .and. present(H_ice)) then 
        ! Update fields from external model (zs and m2), if available
        
        call rembo_get_topo(zs,m2,z_srf,H_ice,z_sl)

      else 
        ! Define fields from initial setup 

        zs = fields0%zs 
        m2 = fields0%m2 
        
      end if

      ! Define relaxation mask
      if (present(mask_relax)) then

        ! Transpose the input relaxation mask
        do j = 1, ny
        do i = 1, nx
          mrelax(j,i) = mask_relax(i,j)
        end do
        end do

      else

        ! Define relaxation mask based on m2
        mrelax = 0
        where (m2 .eq. 2) mrelax = 1
        mrelax(1,:)  = 1
        mrelax(ny,:) = 1
        mrelax(:,1)  = 1
        mrelax(:,nx) = 1
        
      end if

      ! ================================================
      
      ! ! Get current paleo forcing 
      ! !write(*,*) "climate:: forcing:",size(forcing)
      ! if ( boundary_forcing .gt. 0 ) call get_forcing_now(yearnow,m2)
      
      ! ! Determine amount of atmospheric co2 based on warming
      ! ! T_warming = T_global_mean_warming = T_greenland_summer_warming
      ! T_warming = deltaT(0)
      ! aco2 = calc_co2(T_warming)     
      
      aco2 = co2_now 

    end if
    
    ! #### Call REMBO module ####
    if ( climchoice .eq. 1 ) then
    
      if ( equili .eq. 1 ) then ! only called when no restart file is used
        
        ! First store the boundary settings temporarily, equilibrate model
        ! with average over whole period
        bndtmp = bnd_ave
        
        if ( bnd_trans .eq. 1 ) bnd_ave   = 100
        
        ! Call the energy-moisture balance module (with fixed topography)
        ! until it approaches equilibrium
        now%clim = .TRUE.; now%smb = .TRUE.
        do qq = 1, n_equili
          call rembo(yearnow,yearnow1,now,m2,zs,lats,lons,mrelax,aco2,T_warming,T_anomaly)
          write(*,*) "rembo equili, ",qq
        end do
        
        write(*,*) "Finished with rembo equilibration."

        equili = 0
        
        ! Reset boundary variables
        bnd_ave   = bndtmp
        
      end if
                
      ! Call the energy-moisture balance module (output sent to variable "saved")  
      call rembo(yearnow,yearnow1,now,m2,zs,lats,lons,mrelax,aco2,T_warming,T_anomaly)      

    else   ! climchoice=0 => bilinear interp + data 
      
      ! Call conventional approach (output sent to variable "saved") 
      if ( now%smb ) call conventional(m2,zs,lats,lons,aco2,T_warming,T_anomaly,time)   
      
    end if
    
    ! #### OUTPUT SECTION ####
    
    ! Determine the current output index
    n_now = match(yearnow1,time_out)
    
    if ( n_now .ne. 0 ) then
      
      !! Output cntro file
      call rembo_nc(trim(outfldr)//"clima.cntro.nc","years",yearnow1, &
                lats*todegs,lons*todegs,sico_grid%x,sico_grid%y,mask_hydro=fields0%mask_hydro, &
                mask_relax=mrelax)
      call rembo_nc_write(trim(outfldr)//"clima.cntro.nc",ann,1,yearnow1,mask=m2,zs=zs, &
                      tjja=saved%tjja, tdjf=saved%tdjf,ttp=saved%ttp, &
!                       pkappa=saved%pkappa,qqfac=fields0%qqfac, &
                      pdds=saved%pdds, write_mon = .TRUE. )
  
      !! Output annual file
!       call rembo_nc_write(trim(outfldr)//"clima.nc",ann,n_now,yearnow1,mask=m2,zs=zs, &
!                           tjja=saved%tjja, tdjf=saved%tdjf,ttp=saved%ttp, &
! !                           pkappa=saved%pkappa,qqfac=fields0%qqfac, &
!                           pdds=saved%pdds, write_mon = .TRUE. )
      ! call rembo_nc_write_small(trim(outfldr)//"clima.nc",ann,mon,n_now,yearnow1,mask=m2,zs=zs, &
      !                           tjja=saved%tjja,tdjf=saved%tdjf,ttp=saved%ttp,pdds=saved%pdds)

      !! Output monthly file
      if ( write_rembo_m .eq. 1 ) then
        
        call rembo_nc(trim(outfldr)//"climm.cntro.nc","months",time_out_mon(1),lats*todegs,lons*todegs, &
                      sico_grid%x,sico_grid%y,mask=m2,zs=zs, &
                mask_relax=mrelax)
          
        do k = 1, nm
          call rembo_nc_write(trim(outfldr)//"climm.cntro.nc",mon(k),k,time_out_mon(k))
        end do
      end if   
      
      !! Output daily file
      if ( write_rembo_d .eq. 1 ) then      ! Write daily output
        ! The same bug applies as above in the monthly file!!
        call rembo_nc(trim(outfldr)//"climd.cntro.nc","days",time_out_day(1),lats*todegs,lons*todegs, &
                      sico_grid%x,sico_grid%y,mask=m2,zs=zs, &
                mask_relax=mrelax)
          
        do k = 1, nk
          call rembo_nc_write(trim(outfldr)//"climd.cntro.nc",day(k),k,time_out_day(k))
        end do
      end if

    end if
      
    ! Determine the current output index for restart file
    n_now = match(yearnow1,time_out_r)

    ! Obsolete, restart writing happens externally now.
    ! if ( write_rembo_r .eq. 1 .and. yearnow1 .ne. year0 .and. n_now .ne. 0 ) then
    !   call rembo_restart(yearnow1,m2,zs,day(nk),saved%pdds)  ! Write restart file  
    ! end if
      
    ! ! Now for fractional restart file
    ! if ( write_rembo_r .eq. 2 ) then
    
    !   ! (Initialize tmp file, so rembo won't write restart file,
    !   !  but has something to read)
    !   if ( yearnow1 .ne. year0 ) then
    !     open(99,file=trim(outfldr)//"tmp123456789",status="unknown")
    !     write(99,*) "na"
    !     close(99)
    !   end if
    
    !   ! Get the filename from the tmp file written in sicopolis
    !   open(99,file=trim(outfldr)//"tmp123456789",status="old")
    !   read(99,"(a)") fnm; fnm = adjustl(fnm)
    !   close(99)
      
    !   ! Check the filename, if it's good,
    !   ! write to the restart file
    !   if ( trim(fnm) .ne. "na" ) then
        
    !     write(*,*) "Writing fractional restart file: "//trim(fnm)
        
    !     call rembo_restart(yearnow1,m2,zs,day(nk),saved%pdds,file=trim(fnm))  ! Write restart file
        
    !     ! Wipe out filename in tmp file
    !     open(99,file=trim(outfldr)//"tmp123456789",status="unknown")
    !     write(99,*) "na"
    !     close(99)
        
    !   end if

    ! end if
    
    ! Output various 1D files, by region
    !if ((time-year0) .eq. (time-year0)/dto_clim*dto_clim ) then
    if (time - rembo_ann%time_dto_clim .ge. dto_clim ) then 

      init_summary2_file = init_summary2 .eq. 0

      ! Update current saved time 
      rembo_ann%time_dto_clim = time 

      ! First by sector (five sectors total, hard coded for now)
      if ( clim_coupled .lt. 0 .and. .FALSE. ) then
        do qq = 1, 5
          
          mask = 0.d0
          where( m2 .eq. 0.d0 .and. fields0%mask_hydro .eq. dble(qq) ) mask = 1.d0
          
          write(fnm,"(a11,i1,a3)") "rembo.gis.S",qq,".nc"
          call climchecker_new(trim(outfldr)//trim(fnm), &
                              day,mon,ann,mask,zs,lats,lons,T_warming,T_anomaly,co2_now,time,init_summary2_file)
        
        end do
      end if 
      
      ! Now by total gis area
      mask = 0.d0
      where( m2 .eq. 0.d0 ) mask = 1.d0
      call climchecker_new(trim(outfldr)//"rembo.gis.nc", &
                          day,mon,ann,mask,zs,lats,lons,T_warming,T_anomaly,co2_now,time,init_summary2_file)
      
!!!   ! Kill condition for transient simulations (enough ice points and enough time passed)
!!!      if ( sum(mask) .lt. 150.d0 .and. time .gt. 20d3 ) kill = 1
      
      ! And total greenland area
      mask = 0.d0
      where( m2 .le. 1.d0 ) mask = 1.d0
      call climchecker_new(trim(outfldr)//"rembo.grl.nc", &
                          day,mon,ann,mask,zs,lats,lons,T_warming,T_anomaly,co2_now,time,init_summary2_file)
      
      ! Write/calc comparison to observational fields
      !if ( trim(domain) .eq. "GRL" ) call climchecker2(m2,zs)
      
      ! Ensure init_summary2 is not zero
      init_summary2 = 1

    end if
    
    ! #### END OUTPUT SECTION ####
        
    if (now%clim) then
      write(*,"(a1,5x,a10,f12.2)") "e","yearnow = ",yearnow
      write(*,"(a1,5x,a10,f12.2)") "e","co2 = ",co2_now
    end if
    
    call cpu_time(timer%now)                   ! get current time in seconds
    timer%climate = (timer%now-timer%climate)  ! Get elapsed time in seconds
    
!     if ( clim_coupled .eq. 1 .and. (now%smb .or. now%clim) ) then
!       call vars2ice(saved%tt,saved%tdjf,saved%tjja,saved%ttp,saved%tts, &
!                     saved%tjan, saved%tjul,saved%precip,saved%snow, &
!                     saved%runoff_snow,saved%runoff_rain,saved%melted_ice, &
!                     saved%refrozen,saved%evap,saved%smb,saved%h_snow, "save")
!     end if
    
    ! ================================================
    ! remboyelmo 
if (.TRUE.) then 
    if (now%smb .or. now%clim) then 
      ! Update annual rembo variables ([nx,ny] indices) for external use

      do j = 1, ny 
      do i = 1, nx 

        ! [mm we / a]
        rembo_ann%pr(i,j)     = saved(j,i)%precip 
        rembo_ann%sf(i,j)     = saved(j,i)%snow 
        rembo_ann%smb(i,j)    = saved(j,i)%smb 
        rembo_ann%melt(i,j)   = saved(j,i)%melted_ice 
        rembo_ann%runoff(i,j) = saved(j,i)%runoff_snow + saved(j,i)%runoff_rain
        rembo_ann%refrz(i,j)  = saved(j,i)%refrozen 
        rembo_ann%H_snow(i,j) = saved(j,i)%h_snow

        ! [K] 
        rembo_ann%T_srf(i,j)  = saved(j,i)%tts  + 273.15   ! [degC] => [K]
        rembo_ann%T_ann(i,j)  = saved(j,i)%tt   + 273.15   ! [degC] => [K]
        rembo_ann%T_pann(i,j) = saved(j,i)%ttp  + 273.15   ! [degC] => [K]
        rembo_ann%T_jja(i,j)  = saved(j,i)%tjja + 273.15   ! [degC] => [K]
        rembo_ann%T_djf(i,j)  = saved(j,i)%tdjf + 273.15   ! [degC] => [K]
        
        rembo_ann%pdds(i,j)   = saved(j,i)%pdds 

        ! Relaxation mask
        rembo_ann%mask_relax(i,j) = mrelax(j,i)
      end do 
      end do 

    end if 
end if 
    ! ================================================

    ! Output variables for climber
    !if ( boundary_forcing .eq. 3 .and. (now%smb .or. now%clim) ) then
    !  call write_for_climber(trim(outfldr),nstep,lats*todegs,m2,zs,transT%dVdt)
    !end if 

    return
  
  end subroutine rembo_update
  
  subroutine rembo_init(time)
  
    use emb_global
    use emb_functions
    use rembo_main
  
    implicit none
    
    real(8), intent(IN) :: time 

    integer, parameter :: nx=nxs, ny=nys
    real (8), parameter :: dx=dxs
    
    real (8), dimension(ny,nx) :: m2, zs, lats, lons, aco2
    real (8), dimension(ny,nx) :: ZZ, tma, tmj, ampl
    real (8), dimension(ny,nx) :: restart_tt, restart_pp
    real (8), dimension(ny,nx) :: restart_snowh
    real (8), dimension(ny,nx) :: restart_dh, restart_ap
    real (8), dimension(ny,nx) :: mask
    real (8), dimension(ny,nx) :: Hi 
    real (8) :: wt, pscalar, sfac

    character(len=10) :: c_dxs, dattype
    character(len=256) :: fnm, restart_folder
    character(len=15) :: c_yr
    
    integer :: i, j, k, qq, n, n_now

    ! Set emb global nstep to zero
    nstep = 0
    
      ! ================================================
      ! remboyelmo 
      
      allocate(rembo_ann%pr(nx,ny))
      allocate(rembo_ann%sf(nx,ny))
      allocate(rembo_ann%smb(nx,ny))
      allocate(rembo_ann%melt(nx,ny))
      allocate(rembo_ann%runoff(nx,ny))
      allocate(rembo_ann%refrz(nx,ny))
      allocate(rembo_ann%H_snow(nx,ny))
      allocate(rembo_ann%T_srf(nx,ny))
      allocate(rembo_ann%T_ann(nx,ny))
      allocate(rembo_ann%T_pann(nx,ny))
      allocate(rembo_ann%T_jja(nx,ny))
      allocate(rembo_ann%T_djf(nx,ny))
      allocate(rembo_ann%pdds(nx,ny))
      allocate(rembo_ann%mask_relax(nx,ny))

      rembo_ann%pr     = 0.0 
      rembo_ann%sf     = 0.0 
      rembo_ann%smb    = 0.0 
      rembo_ann%melt   = 0.0 
      rembo_ann%runoff = 0.0 
      rembo_ann%refrz  = 0.0 
      rembo_ann%H_snow = 0.0 
      rembo_ann%T_srf  = 0.0 
      rembo_ann%T_ann  = 0.0 
      rembo_ann%T_pann = 0.0 
      rembo_ann%T_jja  = 0.0 
      rembo_ann%T_djf  = 0.0 
      rembo_ann%pdds   = 0.0 
      rembo_ann%mask_relax = 0 

      allocate(ppcorr0%pp_rembo(nx,ny,nm))
      allocate(ppcorr0%pp_ref(nx,ny,nm))
      allocate(ppcorr0%dpp_corr(nx,ny,nm))
      
      ppcorr0%pp_rembo = 0.0
      ppcorr0%pp_ref   = 0.0
      ppcorr0%dpp_corr = 0.0

      ! ================================================

      call emb_globinit(1)
      
      if (init_rembo .eq. 0) then
        
        ! Allocate day, mon and year structures
        !allocate( day(nk), mon(nm) )
        
        ! Initialize stuff needed for boundary routine called in ini_sico!
        m2        = fields0%m2
        zs        = fields0%zs
        lats      = fields0%lats
        lons      = fields0%lons
        
        ! Calculate the horizontal gradient of the initial topography
        call hgradient(fields0%zs,dx,fields0%dzs)
        
        write(*,"(a1,5x,f12.3,5x,a)") "e",time,"Initialized topography for rembo_update."
        
        ! Initialize output nc file
        call rembo_nc(trim(outfldr)//"clima.nc","years",time_out(1)*1d3, &
                      lats*todegs,lons*todegs,sico_grid%x,sico_grid%y,mask_hydro=fields0%mask_hydro)
                          
        init_rembo = 1

        rembo_ann%time_emb  = time - dtime_emb 
        rembo_ann%time_smb  = time - dtime_smb 
        rembo_ann%dt_emb    = dtime_emb
        rembo_ann%dt_smb    = dtime_smb

        rembo_ann%time_dto_clim   = time - dto_clim 
        rembo_ann%time_dto_clim2D = time - dto_clim2D 
        
        ! Store parameter values too 
        rembo_ann%par%dtime_emb = dtime_emb 
        rembo_ann%par%dtime_smb = dtime_smb 
        rembo_ann%par%nx = nx
        rembo_ann%par%ny = ny 
        rembo_ann%par%dx = dx 

        write(*,*) "rembo_init:: summary"
        write(*,*) "range m2  : ", minval(m2), maxval(m2) 
        write(*,*) "range zs  : ", minval(zs), maxval(zs) 
        
      end if

      call rembo_init_2(m2,time)

    return 

  end subroutine rembo_init

  subroutine rembo_get_topo(zs,m2,z_srf,H_ice,z_sl)

    implicit none 

    real(8), intent(OUT) :: zs(:,:)      ! [ny,nx]
    real(8), intent(OUT) :: m2(:,:)      ! [ny,nx]
    real(8), intent(IN)  :: z_srf(:,:)   ! [nx,ny]
    real(8), intent(IN)  :: H_ice(:,:)   ! [nx,ny] 
    real(8), intent(IN), optional  :: z_sl(:,:)    ! [nx,ny] 

    ! Local variables 
    integer :: i, j, nx, ny 
    real(8), allocatable :: Hi(:,:)  
    real(8), allocatable :: zsl(:,:)  

    nx = size(zs,2)
    ny = size(zs,1) 

    allocate(Hi(ny,nx)) 
    allocate(zsl(ny,nx))

    ! Update fields from external model (zs and m2), if available
    ! (transposed i,j => j,i)

    do j = 1, ny 
    do i = 1, nx 

      ! Update surface elevation
      zs(j,i) = z_srf(i,j)

      ! Update ice thickness 
      Hi(j,i) = H_ice(i,j) 

      ! Update sea level
      if (present(z_sl)) then
        zsl(j,i) = z_sl(i,j)
      else
        zsl(j,i) = 0.0
      end if
      
    end do 
    end do 

    ! Set ice(0)/land(1)/ocean(2) mask
    m2 = 0.0
    where( (zs-zsl) .gt. 0.0 .and. Hi .eq. 0.0) m2 = 1.0 
    where( (zs-zsl) .le. 0.0 .and. Hi .eq. 0.0) m2 = 2.0  
      
    return 

  end subroutine rembo_get_topo

  subroutine rembo_set_time(time)

    implicit none

    double precision, intent(IN) :: time 

    

    return 

  end subroutine rembo_set_time

  subroutine rembo_restart_write(filename,time,z_srf,H_ice,z_sl)

    use emb_global
    use rembo_main

    implicit none

    character(len=*), intent(IN) :: filename
    double precision, intent(IN) :: time
    double precision, intent(IN), optional :: z_srf(:,:)
    double precision, intent(IN), optional :: H_ice(:,:)
    double precision, intent(IN), optional :: z_sl(:,:)
    
    ! Local variables
    double precision :: zs(nys,nxs)
    double precision :: m2(nys,nxs)

    if (present(z_srf) .and. present(H_ice)) then
      ! Get transposed topo information
      call rembo_get_topo(zs,m2,z_srf,H_ice,z_sl)
    else
      ! Define fields from initial setup 
      zs = fields0%zs 
      m2 = fields0%m2 
    end if 

    ! Write restart file
    call rembo_restart(filename,day(nk),m2,zs,time)

    return

  end subroutine rembo_restart_write
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : c o n v e n t i o n a l
  ! Author     : Alex Robinson, 27. Oct 2009
  ! Purpose    : Determine the climate and melt based on conventional 
  !              approaches
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine conventional(mask,zs,lats,lons,aco2,T_warming,T_anomaly,time)
    
    use emb_global
    use emb_functions
    use rembo_functions
    use rembo_main

    implicit none
    
    integer, parameter :: nx = nxs, ny = nys

    double precision, dimension(ny,nx) :: mask, zs, lats, lons, aco2
    double precision, intent(IN) :: time 
    double precision, dimension(ny,nx) :: tmp
    integer :: i, j, k, m, km, km2
    
    double precision, intent(IN) :: T_warming       ! Summer warming average
    double precision, intent(IN) :: T_anomaly(:)    ! Daily temp anomaly values
    
    character(len=10)  :: c_dxs, dattype
    character(len=256) :: fnm
    
    call cpu_time(timer%pddold)           ! get current time in seconds
    
    ! Get the precipitation and temperature fields from data
    ! (they will be stored in the global day/mon/ann structures)
    call precip0(mask,zs,lats,lons,T_warming)
    call temper0(mask,zs,lats,lons,T_anomaly)

    ! #### Begin melt potential calculation ####
    ! Get annual pdd distribution [PDD/a]... pdd: Reinhard's formula; ppb: my formula
    call pdd(mask,ann%tte,ann%pp,ann%snow,ann%rf,ann%runoff, &
             ann%runoff_snow, ann%runoff_rain,ann%melt,ann%melted_ice, &
             ann%ice,ann%refrozen, ann%dh ) 
          
    !  ------ Ice-surface temperature (10-m firn temperature) temp_s,
    !         including empirical firn-warming correction due to
    !         refreezing meltwater when superimposed ice is formed
    ! Firn-warming correction factor, in (d*deg C)/(mm WE)
    saved%tts = min(-1d-3,ann%tt + firn_factor * max(0.d0,-ann%dh) )  ! Minus, bc dh is negative for excess refreezing
    where (mask .ge. 2.d0) saved%tts = -2.d0                ! Set surface T over ocean to -2degC
    
    saved%tt    = ann%tt
    saved%tjan  = mon(1)%tt
    saved%tjul  = mon(7)%tt
    saved%tdjf  = (mon(12)%tt+mon(1)%tt+mon(2)%tt) / 3.d0
    saved%tjja  = (mon(6)%tt+mon(7)%tt+mon(8)%tt) / 3.d0
    
    ! Get the precipitation-weighted annual temperature
    saved%ttp = 0.d0; tmp = 0.d0
    do m = 1, nm
      saved%ttp = saved%ttp + ( mon(m)%tt * mon(m)%pp )
      tmp       = tmp + mon(m)%pp
    end do
    saved%ttp = ( saved%ttp / tmp )
    
    ! ajr: this may be broken, added reset to zero (ajr, 2019-03-31)
    saved%pdds = 0.0
    do m = 1, nm
      saved%pdds = max( mon(m)%tt, 0.d0 ) * day_month
    end do
    
!     saved%precip      = ann%pp          *1d-3  / sec_year /rho_i_w     ! mm.w.e. / a => m.i.e. / s
!     saved%snow        = ann%snow        *1d-3  / sec_year /rho_i_w     ! mm.w.e. / a => m.i.e. / s                                                  
!     saved%melted_ice  = ann%melted_ice  *1d-3  / sec_year /rho_i_w     ! mm.w.e. / a => m.i.e. / s
!     saved%runoff_snow = ann%runoff_snow *1d-3  / sec_year /rho_i_w     ! mm.w.e. / a => m.i.e. / s
!     saved%runoff_rain = ann%runoff_rain *1d-3  / sec_year /rho_i_w     ! mm.w.e. / a => m.i.e. / s
!     saved%smb         = (ann%pp-ann%runoff) *1d-3  / sec_year /rho_i_w     ! mm.w.e. / a => m.i.e. / s
!     saved%evap        = 0.d0
                                                    ! *Note: divide by rho_i_w, bc m3 in denominator!
    
    ! ================================================
    ! remboyelmo 
if (.TRUE.) then 
    ! Make standard output in [mm.w.e. / a]
    saved%precip      = ann%pp                     
    saved%snow        = ann%snow        
    saved%melted_ice  = ann%melted_ice  
    saved%runoff_snow = ann%runoff_snow 
    saved%runoff_rain = ann%runoff_rain 
    saved%refrozen    = ann%refrozen    
    saved%smb         = (ann%ice-ann%melted_ice)
    saved%evap        = 0.d0
    saved%h_snow      = ann%snowh  
    saved%pdds        = ann%tte
end if 
    ! ================================================

    call cpu_time(timer%now)           ! get current time in seconds
    timer%pddold = (timer%now-timer%pddold)  ! Get elapsed time in seconds
    
    write(*,"(a1,5x,i10,5x,a)") "e",time, "Called pdd (annual)."
    
    return
  
  end subroutine conventional
  
  subroutine rembo_calc_precip_corr(dpp_corr,pp_rembo,pp_ref,dx,sigma,max_corr)
    ! This routine takes a reference rembo precipitation array and
    ! a reference target precipitation array and produces a
    ! correction factor that can be used to correct the precipitation 
    ! calculation online. 

    implicit none

    real(8), intent(OUT) :: dpp_corr(:,:)
    real(8), intent(IN)  :: pp_rembo(:,:)
    real(8), intent(IN)  :: pp_ref(:,:)
    real(8), intent(IN)  :: dx
    real(8), intent(IN)  :: sigma
    real(8), intent(IN)  :: max_corr

    ! Local variables
    integer :: i, j, k, nx, ny, nm
    real(8), allocatable :: pp_rembo_now(:,:)
    real(8), allocatable :: pp_ref_now(:,:)

    real(8), parameter :: eps = 1e-8

    nx = size(dpp_corr,1)
    ny = size(dpp_corr,2)

    ! Initially set correction factor to 1 everywhere for safety
    dpp_corr = 1.0

    ! Store precip fields locally to allow changes
    pp_rembo_now = pp_rembo
    pp_ref_now   = pp_ref

    ! Apply Gaussian smoothing to the precip fields
    
    if (sigma .gt. 0.0) then 
      call smooth_gauss_2D(pp_rembo,dx, sigma / dx)
    end if
    
    if (sigma .gt. 0.0) then 
      call smooth_gauss_2D(pp_ref,dx, sigma / dx)
    end if

    ! Calculate the correction factor, point by point

    do j = 1, ny
    do i = 1, nx

      ! Calculate correction factor
      dpp_corr(i,j) = pp_ref(i,j) / (pp_rembo(i,j)+eps)

      ! Limit to desired range (e.g., 1.0±0.5)
      if (dpp_corr(i,j) .gt. 1.0+max_corr) then
        dpp_corr(i,j) = 1.0+max_corr
      else if (dpp_corr(i,j) .lt. 1.0-max_corr) then
        dpp_corr(i,j) = 1.0-max_corr
      end if

    end do
    end do

    return

  end subroutine rembo_calc_precip_corr

end module rembo_sclimate
