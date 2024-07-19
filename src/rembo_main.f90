module rembo_main
  
  use emb_global
  use emb_functions
  use emb_pdesolver
  use rembo_functions

  ! Decide which surface mass balance module to use
  use smb_itm
  
  implicit none
  
  type emb_vars
    double precision, dimension(nye,nxe) ::   &
        tt, pp, qq, snow, rain, S, ap, lats,  &
        rhum, dh, rrco2, dS, &
        dpp_corr
  end type
  
  type(emb_vars) :: low
  type(emb_out) :: emb_d
  
  ! Global data (initialized once for emb - or updated as needed)
  type rembo_initvars
    integer :: niter
    double precision, allocatable, dimension(:,:) :: wtst
    double precision, allocatable, dimension(:,:) :: wtsp
    double precision, dimension(nye,nxe) :: zs, lats, lons, xx, yy
    double precision, dimension(nk,nye,nxe) :: tt,pp,rhum, S0
    integer, dimension(nk) :: m0, m1
    double precision, dimension(nk) :: wt0, wt1
    double precision, dimension(:,:,:), allocatable :: t2m0,rhum0
    double precision, dimension(:), allocatable :: time0
    integer :: bndstep

    double precision, dimension(nk,nye,nxe) :: dpp_corr
  end type
  
  type(rembo_initvars) :: rembo0

!   type(rembo_type), allocatable :: day(:), mon(:)
  type(rembo_type) :: day(nk), mon(nm)
  type(rembo_type) :: ann
    
  type choices
    logical :: clim, smb
    double precision :: time_now, time_prev, dt_now 
  end type   
  
  type tuning_vars
    double precision, dimension(nye,nxe) :: pp,pp00
    double precision, dimension(nye,nxe) :: qqfac
    double precision, dimension(nye,nxe) :: pkappa, p_tau
  end type
  
  type(tuning_vars) :: paramfields
  
  ! Counters, switches
  integer :: nout_rembo, nout2, ap_stored
  
contains
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  r e m b o
  ! Author     :  Alex Robinson
  ! Purpose    :  Energy moisture balance combined with surface interface
  !               Model should be run for an entire year, including
  !               a few days in the beginning for equilibration
  !               Daily data will be written to a file.
  !               Annual data will be output as arguments.
  !               Precip and temperature data for a yearly cycle
  !               will be output per day (or per timestep required by
  !               sicopolis).
  ! dt [days]
  ! dxl [m]
  ! P  [kg/m2-s]
  ! T  [°C]   - SURFACE TEMPERATURE !!
  ! rhum [0:1]
  !
  ! === NCEP DATA ===
  ! t0 [°C]      - surface temperature
  ! p0 [kg/m2-s] - precip rate
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine rembo(yearnow,yearnow1,now,m2,zs,lats,lons,mrelax,aco2,T_warming,T_anomaly)

    implicit none
    
    double precision, intent(IN) :: yearnow, yearnow1
    double precision, intent(IN) :: T_warming           ! Summer anomaly
    double precision, intent(IN) :: T_anomaly(:)        ! Daily anomalies
    integer, parameter :: nxl=nxe, nyl=nye
    double precision, parameter :: dxl = dxe
    double precision :: dt
    
    integer, parameter :: nx=nxs, ny=nys
    double precision, parameter :: dx = dxs
    
    integer :: i, j, k, iter, kd, qd, m, km, km2
    integer :: qq, qqmax
    
    double precision :: get_year
    
    ! Parameters
    double precision :: t_kappa, lapse
    double precision :: dt_sec
    double precision :: sigma
    
    ! Lo-res
    double precision, dimension(nyl,nxl) :: p_kappa, kappa_p, kappa_t, pkfac, tkfac
    double precision, dimension(nyl,nxl) :: aco2l  , tmp                       
    double precision, dimension(nyl,nxl) :: relaxt, relaxp
    double precision, dimension(nyl,nxl) :: m2l, m2lp, zsl, zslp, dzslp
    double precision, dimension(nk,nyl,nxl) :: S
    ! Hi-res
    double precision, dimension(ny,nx) :: m2, zs, zs_0, lats, lons, aco2, mmfac, pp_tmp, dzsp
    integer, dimension(ny,nx) :: mrelax

    type(choices) :: now
    
    ! Local boundary variables (can be modified for warming, etc)
    ! Also, array for storing daily value until total iterations met
    type(emb_vars) :: low0(nk), lowd
    
    double precision :: lonmin, lonmax, klon
    double precision, dimension(nyl,nxl) :: lonfac
    
    !double precision, dimension(nyl,nxl) :: latsl

    ! Set local variables equal to global values
    dt = dte
    dt_sec = dt*sec_day
    
    t_kappa   = tkappa
    lapse     = tlfac
    p_kappa   = pkappa     ! m2 / K,    diffusion coefficient
    pkfac     = 1.d0 
    kappa_p   = p_kappa*pkfac 

    sigma = Teff_sigma
    
    lonmin = minval(emb_grid%x)
    lonmax = maxval(emb_grid%x)
    klon   = 0.4d0
    
    ! ! #### Initialize variables to be used in rembo model,
    ! !      only if it not already initialized
    ! if ( init_rembo .lt. 2 ) then
      
    !   nout_rembo = 0
      
    !   call ini_rembo(rembo0,latsl)    
    !   write(stdout,"(f12.3,5x,a)") yearnow1, "Called ini_rembo."

    !   ! Set snow to limit for the first timestep
    !   ! on Greenland (land and ice)
    !   ! Only if starting with ice present, otherwise init. snow level is 0.d0
    !   if (anf_dat .eq. 1) then
    !     where (m2 .lt. 2.d0) day(nk)%snowh = h_snow_max
    !     write(*,"(a1,5x,f12.3,5x,a)") "e",yearnow1,"Snow height set to maximum on land!"
    !   end if
      
    !   ! Initialize the PDDs representing vegetation to zero
    !   veg_pdds = 0.d0 
      
    !   ! Set saved flag to zero (indicating that ap field has not been saved yet)
    !   ap_stored = 0
      
    !   if (anf_dat .eq. 3) then  ! Get initial fields from previously saved run
      
    !     day(nk)%tt    = fields0%tt
    !     day(nk)%snowh = fields0%snowh
    !     day(nk)%dh    = fields0%dh
    !     day(nk)%ap    = fields0%ap
        
    !     veg_pdds      = fields0%pdds
        
    !     write(*,"(a1,5x,f12.3,5x,a)") "e",yearnow1,"Loaded REMBO restart variables!"
    !   end if
      
    !   ! Determine what the melt rate factor should be
    !   ! (currently constant for the entire grid
    !   melt_ice = mm_teff_ice
      
    !   write(stdout,"(f12.3,5x,a)") yearnow1, "Further smb initializations."
    !   init_rembo = 2

    ! end if
    
    ! #### Initial calculations, depending on module to run ####
    
    if ( now%clim .or. now%smb ) then
      call cpu_time(timer%emb) ! Get start time in seconds  
    end if
    
    if ( now%clim ) then
      
      ! Update the boundary data fields
      call update_boundaries(rembo0,yearnow)
      
      ! Get the daily insolation for a year at low-resolution
      ! ajr: 1950 represents present day => change for transient runs!   
      call sinsol2d(S,rembo0%lats,yearnow) 

      ! Ensure elevation is greater than 0.0 for aggregation
      zs_0 = max(zs,0.d0)
      
      ! Rescale large scale fields to low-res diffusion grid
      call tolores(zs_0,  zsl, rembo0%wtst,nrt,ratio)
      call tolores(zs_0, zslp, rembo0%wtsp,nrp,ratio)
      call tolores(m2,    m2l, rembo0%wtst,nrt,ratio,m="m2")
      call tolores(m2,   m2lp, rembo0%wtsp,nrp,ratio,m="m2")
      
      ! TUNING (and finished equilibrating)
      if ( tuning .eq. 1 ) then
        
        !! Calculate p_tau matrix
        paramfields%p_tau = p_tau
        
!         lonfac = ( 1.d0 + klon*(2.d0*(rembo0%xx-lonmin)/(lonmax-lonmin) - 1.d0) )
!         where ( m2l .eq. 2.d0 ) lonfac = 1.d0
!         paramfields%p_tau = p_tau * lonfac     
!         where ( paramfields%p_tau .lt. 3.d0 ) paramfields%p_tau = p_tau
        
       !fields0%qqfac = 1.d0
        call tune_qqfac(ann%pp,fields0%precip,m2,lons,fields0%qqfac)
        
        ! Aggregate to lo-res field.
        call tolores(fields0%qqfac, paramfields%qqfac, rembo0%wtst, nrt, ratio)
        call tolores(ann%pp, paramfields%pp, rembo0%wtst, nrt, ratio)
        
        ! Tune the precipitation kappa and qqfac
        call tune_pkappa(paramfields%pp,paramfields%pp00,m2lp,paramfields%pkappa)
        
        
        p_kappa = paramfields%pkappa
      else
        
        fields0%qqfac = 1.d0
        paramfields%qqfac = 1.d0
        paramfields%p_tau = p_tau
        paramfields%pkappa = pkappa
        p_kappa = paramfields%pkappa
        
      end if
      
      ! Get scaled horizontal gradient of elevation (for precip diffusion)
      call precip_grad(zslp,dxl,rembo0%lats,dzslp,klat=p_k_lat,kfac=p_k_eastfrac)
      
!       call precip_grad(zs_0,dx,fields0%lats,dzsp,klat=p_k_lat,kfac=p_k_eastfrac)
!       call tolores(dzsp,dzslp,rembo0%wtst,nrt,ratio)
      
      ! ! Set relaxation term on grid according to m2
      ! relaxt = 0.d0
      ! where (m2l .eq. 2.d0)  relaxt = 1.d0       
      ! relaxp = 0.d0
      ! where (m2lp .eq. 2.d0) relaxp = 1.d0 
      
      ! Interpolate relaxation mask from high resolution input
      call tolores(relaxt,dble(mrelax),rembo0%wtst,nrt,ratio)
      where(relaxt .lt. 0.5d0) relaxt = 0.0d0
      where(relaxt .ge. 0.5d0) relaxt = 1.0d0 

      relaxp = relaxt 

      low%rrco2 = 0.d0
      if (co2_rad .ne. 0) then
        ! Get co2 concentration on lo-res grid, then
        ! determine radiative forcing
        call tolores(aco2,aco2l,rembo0%wtst,nrt,ratio)
        do j=1,nyl
          do i=1,nxl
            low%rrco2(j,i) = Rco2(aco2l(j,i))
          end do
        end do
      end if
      
      ! Load up boundary variables for the entire year
      do k = 1, nk  
        low0(k)%rhum = rembo0%rhum(k,:,:)
      end do
      
      ! Apply latitudinal scaling of temp anomalies
      do k = 1, nk 
        low0(k)%tt   = rembo0%tt(k,:,:)    + T_anomaly(k)   ! Apply boundary warming
        ! low0(k)%tt = rembo0%tt(k,:,:)  &
        !             + T_anomaly(k)*( 1.d0+lat_grad*(latsl-72.6d0*torads)/((76.d0-72.6d0)*torads) )  ! NEEMup - GISP2
      end do

      ! Get the precipitation fields from data if desired
      ! (they will be stored in the global day/mon/ann structures)
      if ( precip .eq. 0 ) then
        
        call precip0(m2,zs,lats,lons,T_warming)
        
        ! Scale precip to lo-resolution
        do k = 1, nk
          call tolores(day(k)%pp  /sec_day,low0(k)%pp,  rembo0%wtst,nrt,ratio)
          call tolores(day(k)%snow/sec_day,low0(k)%snow,rembo0%wtst,nrt,ratio)
          low0(k)%rain = max(low0(k)%pp - low0(k)%snow, 0.d0 )
        end do
      end if
      
      ! Assume precip flux-correction is zero by default
      do k = 1, nk
        low0(k)%dpp_corr = 0.d0 
      end do 

      ! Get the temperature fields from data if desired
      ! (they will be stored in the global day/mon/ann structures)
      if ( temper .eq. 0 )then
      
        call temper0(m2,zs,lats,lons,T_anomaly)
        
        ! Scale temper to lo-resolution (and convert to sea-level T)
        do k = 1, nk
          call tolores(day(k)%tt+lapse*zs_0,low0(k)%tt,rembo0%wtst,nrt,ratio)
        end do
        
      end if
      
! !       if ( temper .eq. 0 .or. precip .eq. 0 ) then
! !         ! Get monthly values
! !         do m = 1, nm
! !           km = (m-1)*day_month+1
! !           km2 = km+day_month-1
! !           call rembo_ave( day(km:km2),mon(m), clim = .TRUE., smb = .FALSE. )
! !         end do
! !         
! !         ! Get annual values
! !         call rembo_ave( mon,ann, clim = .TRUE., smb = .FALSE. )
! !       end if
      
      ! Store initial values to start diffusion process, using boundary data
      ! (should start with values from last day of previous year)
      low%tt   = low0(nk)%tt
      
      ! If hi-res fields have already been calculated at least once,
      ! or starting from a restart file, use these as the initial values
      ! pp not needed yet, because it is calculated directly from tt !
      if ( yearnow1 .gt. year0 .or. anf_dat .eq. 3 ) then
        
        ! Obtain initial T and P fields from global hi-res values
        call tolores(day(nk)%tt+lapse*zs_0, low%tt,rembo0%wtst,nrt,ratio)
        !!call tolores(day(nk)%pp/sec_day,    low%pp,rembo0%wtst,nrt,ratio)                              
        
        write(*,*) "REMBO: Aggregated hi-res field (tt) for initial values.",yearnow1
      end if
      
      ! Save initial planetary albedo (if using fixed albedo option)
      if ( ap_fixed .eq. 1 .and. yearnow1 .gt. year0 ) then
        
        if ( ap_stored .eq. 0 ) then ! Check to see if array has been stored
          
          ! Store the planetary albedo field in the saved array
          do k = 1, nk
            day(k)%ap0 = day(k)%ap
          end do
          
          ! Indicate that ap field has been stored
          ap_stored = 1

        end if
      
      end if
        
      
    end if
    
    if ( now%smb ) then
    
      if (climchoice .eq. 0) day(nk)%snowh = 0.d0
      
    end if
    
    timer%boundary = 0.d0

    ! #############################################################
    ! #### Main loop over all the days of the year (dt = days) ####
    do kd = 1, nk
      
      ! Get index of previous day
      qd = kd - 1
      if ( kd .eq. 1 ) qd = nk
      
      ! #### Perform climate calculations ####
      if ( now%clim ) then

        ! Get today's insolation onto low-res grid
        ! (as well as difference wrt present day, for pdd)
        low%S  = S(kd,:,:)
        low%dS = low%S - rembo0%S0(kd,:,:)

        ! Get yesterday's albedo, melt height
        call tolores(day(qd)%dh/(1d3*sec_day),low%dh, rembo0%wtst, nrt,ratio)  ! mm/d => m/s     

        if ( ap_fixed .eq. 1 .and. ap_stored .eq. 1 ) day(qd)%ap = day(qd)%ap0
        call tolores(day(qd)%ap,low%ap, rembo0%wtst, nrt,ratio)   
                
        ! Get more data from boundary values
        low%rhum  = low0(kd)%rhum
        low%rrco2 = low%rrco2  ! Same value throughout the year !!
        
        ! Get kappa adjustment (for today)
        call kappa(tkfac,kd,rembo0%lats,rembo0%lons, zsl, m2l, &
                   klat=kappalat,klon=0.d0,kamp=kappaamp,kzs=kappazs,kdir=0.d0)
        !call kappa(pkfac,kd,rembo0%lats,rembo0%lons,zslp, m2l, &
        !           klat=kappalat,klon=kappalon,kamp=kappaamp,kzs=0.0004d0,kdir=4.d0)
        call kappa(pkfac,kd,rembo0%lats,rembo0%lons,zslp, m2l, &
                   klat=kappalat,klon=0.d0,kamp=kappaamp,kzs=0.d0,kdir=0.d0)
                   
        ! Reset daily values to zero before iterations
        lowd%tt   = 0.d0
        lowd%pp   = 0.d0
        lowd%snow = 0.d0
        lowd%qq   = 0.d0
        
        ! Adjust kappa as desired
        kappa_p = p_kappa * pkfac
        
        ! #### solve energy-moisture balances ####
        do iter = 1, rembo0%niter  ! iterate to complete 1 day
          
          ! Get current precipitation
          if (precip .eq. 1) then
          
            ! Perform moisture balance
            call emb_precip(low%pp,low%tt,low%qq,lapse, &
                        low%rhum,zslp,dzslp,relaxp,dxl,dt_sec,kappa_p,&
                        paramfields%qqfac,paramfields%p_tau)
            
            ! scale the precip (for uncertainty studies)
            low%pp = low%pp * pp_scalar(T_anomaly(kd))
            
            ! ---------------------------------------------
            ! ajr: to do
            ! Add flux correction application here

            low%pp = low%pp + low0(kd)%dpp_corr

            ! ---------------------------------------------
            

            ! Calculate other quantities
            low%snow = low%pp * snowfrac(low%tt - lapse*zsl)
            low%rain = low%pp - low%snow
            
          else  ! Current precip is based on boundaries (ie, data)
            
            low%pp   = low0(kd)%pp
            low%snow = low0(kd)%snow
            low%rain = low0(kd)%rain
            
          end if
          
          ! Get current temperature
          if (temper .eq. 1) then
            
            ! Adjust kappa as desired
            kappa_t = t_kappa * tkfac
            
            ! Perform energy balance
            call emb_temper(low%tt,low0(kd)%tt,low%rain,low%snow,low%S,low%ap, &
                    low%dh,low%rrco2,zsl,relaxt,dxl,dt_sec,kappa_t)       
          
          else ! Current temperature is based on boundaries (ie, data)
            
            low%tt = low0(kd)%tt
            
          end if
          
          ! Add this iteration's calculations into the day's total
          lowd%tt   = lowd%tt   + low%tt       ! tt in sea-level T
          lowd%pp   = lowd%pp   + low%pp       ! pp in kg/m2-sec
          lowd%snow = lowd%snow + low%snow     ! snow in kg/m2-sec
          lowd%qq   = lowd%qq   + low%qq       ! qq in kg/m2 
          
        end do  ! End moisture-balance iterations to make one day
      
        ! Get the average value, is value for today
        lowd%tt   = lowd%tt   / dble(rembo0%niter)
        lowd%pp   = lowd%pp   / dble(rembo0%niter) * sec_day    ! kg/m2-sec => kg/m2-day => mm/day
        lowd%snow = lowd%snow / dble(rembo0%niter) * sec_day    ! kg/m2-sec => kg/m2-day => mm/day
        lowd%qq   = lowd%qq   / dble(rembo0%niter)
        
        ! Determine the current output index, and write to output file (if desired)
        ! Get the current year, as relevant to output
        !yearnow1 = time 

        if ( write_emb_d .eq. 1 .and. match(yearnow1,time_out*1d3) .ne. 0 ) then
          if ( kd .eq. 1 ) then
            call emb_nc(trim(outfldr)//"embd.nc",t0=dble(kd),units="days", &
                        lat=rembo0%lats*todegs,lon=rembo0%lons*todegs, &
                        xx=emb_grid%x,yy=emb_grid%y,mask=m2l,&
                        zs=zsl,zsp=zslp,dzsp=dzslp,pkap=kappa_p,ptau=paramfields%p_tau)
            
          end if
          
          write(*,*) kd,"S ave: ", sum(low%S) / (nxl*nyl)
          
          ! Get all of the necessary forcing quantities
          emb_d%tt   = lowd%tt-lapse*zsl
          emb_d%pp   = lowd%pp
          emb_d%qq   = lowd%qq
          emb_d%snow = lowd%snow
          emb_d%rhum = low%rhum
          emb_d%S    = low%S
          emb_d%co2  = low%rrco2
          emb_d%ap   = low%ap
          emb_d%apS  = (1.d0-low%ap) * low%S
          emb_d%ABT  = ta + tb*emb_d%tt
          call d2u_dx2(lowd%tt,dxl,tmp)        ! Use sea-level T
          emb_d%d2T  = kappa_t * tmp
          emb_d%Lr   = Lw*(lowd%pp-lowd%snow) / sec_day   ! kg/m2-s * J/kg = W/m2
          emb_d%Ls   = Ls*lowd%snow / sec_day             ! kg/m2-s * J/kg = W/m2
          emb_d%Ldh  = low%dh*rho_w*Lm
          
          ! Write the daily emb file for current day
          call emb_nc_write(trim(outfldr)//"embd.nc",emb_d,kd)
          
        end if 
        
        ! Interpolate back to high res grid       
        if (temper .ge. 1) then
            call tohires(lowd%tt,day(kd)%tt,ratio,0)
            day(kd)%tt = day(kd)%tt-lapse*zs_0    ! Convert sea-level T to surface T 
            
            ! smbtables
            day(kd)%tt = day(kd)%tt + smb_dT
        end if
          
        if (precip .ge. 1) then
            call tohires(lowd%pp,day(kd)%pp,ratio,0)    ! mm/day
            
            ! smbtables
            day(kd)%pp = day(kd)%pp*(1.d0+smb_ppfac)
            
            ! Determine the amount of precip falling as snow [mm]
            day(kd)%snow = day(kd)%pp * snowfrac(day(kd)%tt)
        end if

        ! Determine effective temperature at each point (+degC)
        day(kd)%tte = effectiveT(day(kd)%tt,sigma) * sec_frac
        
        ! Interpolate insolation to hi-resolution
        call tohires(low%S,day(kd)%S,ratio,0)

      end if
      
      ! #### Perform melt calculations for today ####
      if ( now%smb ) then

        ! Store the snow height and net melt from yesterday
        ! as initial conditions [mm]
        day(kd)%snowh = day(qd)%snowh
        day(kd)%dh    = day(qd)%dh
        
        mmfac = 1.d0   ! ratio between snow and ice melt is 1 (only if melt_choice==1)
        day(kd)%melt_insol = 0.d0 ! additional insolation melt is zero
        if (melt_choice .le. 0) then
          mmfac   = melt_ice / mm_teff_snow    ! ratio to convert snow melt to ice melt 

          ! Calculate change in melt due to insolation differences (mm/day)
          day(kd)%melt_insol = pdd_a * sec_day * (day(kd)%S - day(kd)%S0) *1e3

        else if (melt_choice .eq. 1 .and. itm_S0 .gt. 0) then 

          call get_pdd_corr(day(kd)%pdd_corr,day(kd)%tt,kd)
          !call get_pdd_corr_melt(day(kd)%pdd_corr,day(qd)%melt)

          !day(kd)%melt_insol = day(kd)%pdd_corr*(sec_day/sec_year)*(day(kd)%S-day(kd)%S0)
          day(kd)%melt_insol = day(kd)%pdd_corr*(day(kd)%S-day(kd)%S0)

        end if
      
        call cpu_time(timer%tmp)           ! get current time in seconds
        
        ! Get daily melt
        call melt_budget( m2, zs, day(kd)%S, day(kd)%tt, day(kd)%tte, veg_pdds,  &
                          day(kd)%pp, day(kd)%snow, mmfac, day(kd)%as, day(kd)%ap,  &
                          day(kd)%snowh, day(kd)%rf, day(kd)%runoff,   &
                          day(kd)%runoff_snow, day(kd)%runoff_rain, day(kd)%melt,  &
                          day(kd)%melted_ice, day(kd)%ice, day(kd)%refrozen, &
                          day(kd)%dh, day(kd)%melt_insol, day(kd)%melt_S0, day(kd)%S0  )                          
        
        ! Output the diagnostic pdd correction term
        if (melt_choice .eq. 1 .and. itm_S0 .eq. 1) then 
          day(kd)%pdd_corr = 0.0 
          where ( day(kd)%S-day(kd)%S0 .ne. 0.d0) &
              day(kd)%pdd_corr = (day(kd)%melt-day(kd)%melt_S0) / (day(kd)%S-day(kd)%S0)
        end if

        call cpu_time(timer%now)           ! get current time in seconds
        timer%boundary = timer%boundary + (timer%now-timer%tmp)
      
      end if
      
    end do  ! End daily loop
    
    ! #### Sum up totals over the year ####
    if ( now%clim .or. now%smb ) then
      
      ! Get monthly values
      do m = 1, nm
        km = (m-1)*day_month+1
        km2 = km+day_month-1
        call rembo_ave( day(km:km2),mon(m), now%clim )
      end do
      
      ! Get annual values
      call rembo_ave( mon,ann, now%clim )
      
      veg_pdds = 0.d0
      do k = 1, nk
          veg_pdds = veg_pdds + max( day(k)%tt, 0.d0 )
      end do
      
      ! #### Re-calculate melt for the year using annual PDD approach ####
      if ( melt_choice .eq. -1 ) then
        
        ! Set some annual variables to zero to match 'annual approach'
        ann%snowh = 0.d0; ann%dh = 0.d0
        
        ! Get annual melt (like standard annual PDD approach)
        call melt_budget( m2, zs, day(kd)%S, ann%tt, ann%tte*day_year, veg_pdds,  &
                          ann%pp, ann%snow, mmfac, ann%as, ann%ap,  &
                          ann%snowh, ann%rf, ann%runoff,   &
                          ann%runoff_snow, ann%runoff_rain, ann%melt,  &
                          ann%melted_ice, ann%ice, ann%refrozen, &
                          ann%dh, ann%melt_insol, ann%melt_S0, ann%S0, &
                          annual=.TRUE. ) 
                          
        ! Get annual pdd distribution [PDD/a]... pdd: Reinhard's formula; ppb: my formula
!         call pdd(m2,ann%tte,ann%pp,ann%snow,ann%rf,ann%runoff, &
!               ann%runoff_snow, ann%runoff_rain,ann%melt,ann%melted_ice, &
!               ann%ice,ann%refrozen, ann%dh ) 
        
      end if
      
      ! ## Get the Janssens/Huybrechts refreezing values
      !call calc_rf(m2,ann%tt,ann%pp,ann%snow,ann%melt,saved%jh_rf,saved%jh_refrozen)
      
      ! #### Save annual output to global fields ####
      ! (Store structures in standard arrays for output, in SI units)
      
      !  ------ Ice-surface temperature (10-m firn temperature) temp_s,
      !         including empirical firn-warming correction due to
      !         refreezing meltwater when superimposed ice is formed
      ! Firn-warming correction factor, in (d*deg C)/(mm WE)
      saved%tts = min(-1d-3,ann%tt + firn_factor * max(0.d0,-ann%dh) )  ! Minus, bc dh is negative for excess refreezing
      where (m2 .ge. 2.d0) saved%tts = -2.d0                ! Set surface T over ocean to -2degC
      
      saved%tt    = ann%tt
      saved%tjan  = mon(1)%tt
      saved%tjul  = mon(7)%tt
      saved%tdjf  = (mon(12)%tt+mon(1)%tt+mon(2)%tt) / 3.d0
      saved%tjja  = ( mon(6)%tt+mon(7)%tt+mon(8)%tt) / 3.d0

      call rembo_array(savedm%S,mon,varname="S")
      call rembo_array(savedm%S0,mon,varname="S0")
      call rembo_array(savedm%tt,mon,varname="tt")
      call rembo_array(savedm%melt,mon,varname="melt")
      call rembo_array(savedm%as,mon,varname="as")
      call rembo_array(savedm%pdd_corr,mon,varname="pdd_corr")
      
      ! Interpolate pkappa to hi-resolution
      call tohires(kappa_p,saved%pkappa,ratio,0)
      call tohires(paramfields%pkappa,fields0%pkappa,ratio,0)
      call tohires(paramfields%qqfac,fields0%qqfac,ratio,0)
      
      ! Get the precipitation-weighted annual temperature
      saved%ttp = 0.d0; pp_tmp = 0.d0
      do k = 1, nk
        saved%ttp = saved%ttp + ( day(k)%tt * day(k)%pp )
        pp_tmp       = pp_tmp + day(k)%pp
      end do
      saved%ttp = ( saved%ttp / pp_tmp )
      
      ! Save the vegetation PDDs field
      saved%pdds = veg_pdds
      
!       ! Convert saved quantities for sicopolis: mm.w.e. / a => m.i.e. / s
!       ! *Note: divide by rho_i_w, bc m3 in denominator!
!       saved%precip      = ann%pp          *1d-3  / sec_year /rho_i_w           
!       saved%snow        = ann%snow        *1d-3  / sec_year /rho_i_w
!       saved%melted_ice  = ann%melted_ice  *1d-3  / sec_year /rho_i_w
!       saved%runoff_snow = ann%runoff_snow *1d-3  / sec_year /rho_i_w
!       saved%runoff_rain = ann%runoff_rain *1d-3  / sec_year /rho_i_w
!       saved%refrozen    = ann%refrozen    *1d-3  / sec_year /rho_i_w
!       saved%smb         = (ann%ice-ann%melted_ice) *1d-3  / sec_year /rho_i_w
!       saved%evap        = 0.d0
!       saved%h_snow      = ann%snowh       *1d-3  / rho_i_w
      
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
      saved%pdds        = veg_pdds 
end if 
      ! ================================================

      call cpu_time(timer%now)           ! get current time in seconds
      timer%emb = (timer%now-timer%emb)  ! Get elapsed time in seconds
      
      write(stdout,*)
      write(stdout,"(f12.3,a30,f10.3,1x,f10.3)") yearnow1,"time EMB [sec]", &
                            timer%emb, timer%boundary

    end if
    
    return
    
  end subroutine rembo
  
  subroutine rembo_init_2(m2,yearnow1)

    implicit none 

    !double precision :: latsl(nxe,nye)
    double precision, dimension(nys,nxs) :: m2
    double precision :: yearnow1 

    ! #### Initialize variables to be used in rembo model,

      nout_rembo = 0
      
      call ini_rembo(rembo0) !,latsl)    
      write(stdout,"(f12.3,5x,a)") yearnow1, "Called ini_rembo."

      ! Set snow to limit for the first timestep
      ! on Greenland (land and ice)
      ! Only if starting with ice present, otherwise init. snow level is 0.d0
      if (anf_dat .eq. 1) then
        where (m2 .lt. 2.d0) day(nk)%snowh = h_snow_max
        write(*,"(a1,5x,f12.3,5x,a)") "e",yearnow1,"Snow height set to maximum on land!"
      end if
      
      ! Initialize the PDDs representing vegetation to zero
      veg_pdds = 0.d0 
      
      ! Set saved flag to zero (indicating that ap field has not been saved yet)
      ap_stored = 0
      
      if (anf_dat .eq. 3) then  ! Get initial fields from previously saved run
      
        day(nk)%tt    = fields0%tt
        day(nk)%snowh = fields0%snowh
        day(nk)%dh    = fields0%dh
        day(nk)%ap    = fields0%ap
        
        veg_pdds      = fields0%pdds
        
        write(*,"(a1,5x,f12.3,5x,a)") "e",yearnow1,"Loaded REMBO restart variables!"
      end if
      
      ! Determine what the melt rate factor should be
      ! (currently constant for the entire grid
      melt_ice = mm_teff_ice
      
      write(stdout,"(f12.3,5x,a)") yearnow1, "Further smb initializations."
      init_rembo = 2

    return 

  end subroutine rembo_init_2

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i n i _ r e m b o
  !   Author     :  Alex Robinson
  !   Purpose    :  Initialize energy moisture balance model
  !   dt [days]
  !   dx [m]
  !   P  [kg/m2-s]
  !   T  [°C]   - SURFACE TEMPERATURE !!
  !   rhum [0:1]
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ini_rembo(init) !,latsl)
            
    implicit none

    integer,  parameter :: nxl = nxe, nyl = nye
    double precision, parameter :: dxl = dxe
    
    double precision :: dt
    double precision :: chck,chck2, t_tot
    double precision :: wt0, wt1
    integer :: i, j, k, n_iter, nt, q, m
    character(len=512) :: fnm
    
    type(rembo_initvars) :: init
    !double precision :: latsl(nxl,nyl)

    dt = dte
    
    ! Determine the number of iterations for 1 day of calculations
    init%niter = ceiling( 1.d0 / dt )
    
    ! Check the chosen timesteps for numerical stability
    chck  = 4.d0*tkappa*(dt*sec_day) / &
            (dxl*dxl*tce)
    chck2 = 4.d0*(pkappa+10d4)*(dt*sec_day) /  &
            (dxl*dxl*pce)
    write(stdout,"(a30,f10.3)")  &
          "Max temperature dt possible:", dt/chck
    write(stdout,"(a30,f10.3)")  &
          "Max precip dt possible:", dt/chck2
     
    if (solver .eq. 0 .and. max(chck,chck2) .gt. 1.d0) then
      write(stdout,*) "An emb dt is too large for stability!"
      write(stdout,*) "Try again!"
      call embout(1)
    end if
    
! #### MAKE CALCULATIONS FOR RESOLUTION SCALING ####     
           
    ! Get ratio of lowres to highres grid points to know weighting
    if (mod(1.d0/ratio,1.d0) .ne. 0.d0) then
      write(stdout,*) "Low-res grid is not a multiple of high-res grid."
      write(stdout,*) "Interpolation algorithm won't work."
      write(stdout,*) "Try again!"
      call embout(1)
    end if
    
    ! Allocate the interpolation weighting arrays
    if (allocated(init%wtst)) deallocate(init%wtst)
    if (allocated(init%wtsp)) deallocate(init%wtsp)
    allocate(init%wtst(2*nrt+1,2*nrt+1), init%wtsp(2*nrp+1,2*nrp+1))
    
    ! Determine weights of neighbors for different smoothing radii, 
    ! one for precip (wtstp), one for everything else (wtst)...          
    call tolowts(init%wtst,nrt)
    call tolowts(init%wtsp,nrp)
    
    ! Calculate indices and weights for temporal interpolation (months=>days)
    call monthly2daily(init%m0,init%wt0,init%m1,init%wt1)
    
! #### END CALCULATIONS ####
    
    ! Allocate the day and month objects 
!     allocate(day(nk),mon(nm))

    ! #### LOAD / GENERATE BOUNDARY DATA ####
    ! Get filename of boundary data
    fnm = trim(emb_clim_file)
    
    ! Determine length of boundary fields (in time - how many months of data?)
    nt = nc_size(fnm,"time")*12 
    
    if (allocated(init%rhum0)) deallocate(init%rhum0)
    if (allocated(init%t2m0) ) deallocate(init%t2m0)
    if (allocated(init%time0)) deallocate(init%time0)
    allocate(init%time0(nt),init%rhum0(nt,nyl,nxl),init%t2m0(nt,nyl,nxl))
    
    ! Specify the initial year of the total boundary data set being loaded
    ! 1957.09 for ERA40, etc..
    init%time0(1) = bnd_yr0
    
    do k = 2, nt
      init%time0(k) = init%time0(k-1)+0.01d0
      if ( get_month(init%time0(k)) .eq. 13 ) &
               init%time0(k) = dble(floor(init%time0(k)))+1.01d0
    end do
    
    ! Set boundary nsteps to zero (since it is just starting)
    init%bndstep = 0
    
    ! Get main topo variables
    call nc_read_t(fnm,"lat2D",init%lats)
    call nc_read_t(fnm,"lon2D",init%lons)
    call nc_read_t(fnm,"x2D",init%xx)
    call nc_read_t(fnm,"y2D",init%yy)
    call nc_read_t(fnm,"zs",init%zs)

    ! read climate variables first (monthly for many years)
    !call nc_read_t(fnm,"t2m", init%t2m0)
    !call nc_read_t(fnm,"rhum",init%rhum0)
    
    q = 0 
    do k = 1, nt/12
      do m = 1, 12 
        q = q+1
        call nc_read_t(fnm,"t2m", init%t2m0(q,:,:), start=[1,1,m,k],count=[nxl,nyl,1,1])
        call nc_read_t(fnm,"rhum",init%rhum0(q,:,:),start=[1,1,m,k],count=[nxl,nyl,1,1])

      end do 
    end do 
    
    ! Convert from Kelvin to Celcius 
    init%t2m0 = init%t2m0 - 273.15 

    ! Convert lat/lon to rads
    init%lats = init%lats * torads
    init%lons = init%lons * torads
    
    ! Scale rhum to a fraction from percent
    init%rhum0 = init%rhum0 / 100.d0 
    
    ! Get data precip field to lo-res
    call tolores(fields0%precip, paramfields%pp00, rembo0%wtst, nrt, ratio)
    
    ! Get the daily insolation for a year at low-resolution for present day,
    ! Then store as hi-res daily and monthly fields
    call sinsol2d(init%S0,init%lats,0.d0)
    do k = 1,nk
      call tohires(init%S0(k,:,:),day(k)%S0,ratio,0)
    end do

    ! Get the latitude on lo-resolution
    !call tolores(init%lats,latsl,rembo0%wtst,nrt,ratio)

    write(*,*) "ini_rembo summary"
    write(*,*) "range zs    : ", minval(init%zs),    maxval(init%zs) 
    write(*,*) "range t2m0  : ", minval(init%t2m0),  maxval(init%t2m0)
    write(*,*) "range rhum0 : ", minval(init%rhum0), maxval(init%rhum0)
    write(*,*) "range S0    : ", minval(init%S0),    maxval(init%S0)
    
    return
          
  end subroutine ini_rembo
  
  ! Get the integer corresponding to the current month
  function get_month(year)
    
    implicit none
    
    integer :: get_month, k
    double precision :: year
    
    
    k = idnint( (year - floor(year))*100 )
    
    get_month = k

    return
    
  end function get_month
  
  
  
  subroutine update_boundaries(init,time)
    
    implicit none
    
    double precision, intent(IN) :: time 

    integer,  parameter :: nxl = nxe, nyl = nye
    double precision, parameter :: dxl = dxe
    
    double precision :: year_now
    double precision, dimension(nm,nyl,nxl) :: t2m,rhum
    double precision, dimension(nyl,nxl) :: hgt0
    double precision :: lapse, wt0, wt1, mon(12)
    double precision :: bnd_start1, bnd_end1, bnd_offset
    integer :: k, m0, m1, nt
    character(len=50) :: fnm
    
    type(rembo_initvars) :: init

    lapse = tlfac
    
    nt = size(init%time0)
    
    ! Get the current year
    year_now = time
    
    ! Calculate how much to offset the start date, if averaging is used
    ! (so that the current year lies in the middle of the averaging period)
    bnd_offset = floor(bnd_ave / 2.d0)
    
    ! Figure out bnd_start and bnd_end
    select case(bnd_trans)
      case(1,2)      ! transient boundary fields
        
        ! Adjust start year of boundary averaging with timesteps and the offset
        bnd_start1 = year_now - bnd_offset
        
        ! Ensure the start year is withing bounds
        if ( bnd_start1 .gt. init%time0(nt) ) then
          
          if ( bnd_trans .eq. 1 ) then   ! Keep at last year range
            bnd_start1 = init%time0(nt) - bnd_offset
          else if ( bnd_trans .eq. 2 ) then    ! Loop to initial conditions
            bnd_start1 = bnd_start  ! ajr: doesnt work yet!! Need to fix nstep problem...
          end if
        end if
        
      case DEFAULT ! 0 - intransient boundary fields
        
        bnd_start1 = bnd_yr0
!         if ( timetype .eq. 1 ) bnd_start1 = year0
        
    end select
    
    bnd_end1 = bnd_start1+bnd_ave
    
    write(*,"(a,f12.1,2f10.2)") "year_now,bnd_start,bnd_end",year_now,bnd_start1,bnd_end1  
    !write(*,"(a,i6,2f10.2)")  "nt, time(1),time(nt)", nt,init%time0(1),init%time0(nt)
    
    ! #### GENERATE CURRENT BOUNDARY DATA ####

    select case(trim(domain))
      case("calov-island")  ! (reinhard, calov-island)
        init%lats = 60.d0
        init%lons = -10.d0
        init%zs   = 0.d0
        init%rhum = 90.d0

        ! Generate daily temperature field for boundaries...
        do k = 1, nk
          init%tt(k,:,:) = -tempamp*dcos((k-15)*2.d0*pi/day_year)
        end do
        ! End daily temperature field
      
      case DEFAULT  ! GRL = greenland

        !! Average the correct months and years to get
        !! the monthly boundary data
        mon = 0.d0; t2m = 0.d0; rhum = 0.d0
        do k = 1, nt
          
          ! If the current month/year fits in range, then add it to the
          ! boundary field
          if ( init%time0(k) .ge. bnd_start1 .and. init%time0(k) .lt. bnd_end1 ) then
            
            m1 = get_month(init%time0(k))
            mon(m1) = mon(m1) + 1.d0
            t2m(m1,:,:)  =  t2m(m1,:,:) +  init%t2m0(k,:,:)
            rhum(m1,:,:) = rhum(m1,:,:) + init%rhum0(k,:,:)
            
          end if
        
        end do
        
        ! Make sure enough months have been used...
        !write(*,"(a,12f8.1)") "mon = ",mon
        if ( minval(mon) .eq. 0.d0 ) then
          write(*,*) "Not enough months for averaging boundary data. Try again!"
          stop
        end if

        ! Now make average by dividing by number of months
        ! that contributed to the sum
        do m1 = 1, nm
          t2m(m1,:,:)  =  t2m(m1,:,:) / mon(m1)
          rhum(m1,:,:) = rhum(m1,:,:) / mon(m1)
        end do
        
        !! Finished averaging months
        
        
        ! Interpolate to daily values
        do k = 1, nk
          
          ! Load proper interpolation values for today
          m0 = init%m0(k); wt0 = init%wt0(k)
          m1 = init%m1(k); wt1 = init%wt1(k)
          
          ! Get daily fields
          init%tt(k,:,:)   = wt0*t2m(m0,:,:)  + wt1*t2m(m1,:,:)
          init%rhum(k,:,:) = wt0*rhum(m0,:,:) + wt1*rhum(m1,:,:)
        
        end do
        
    end select  

    ! Make sure ocean is set to 0.0 elevation, for lapse rate scaling [m]
    hgt0 = max(0.d0,init%zs)
    
    ! Convert air temperature to sea-level T, and any offset required
    do k = 1, nk
      init%tt(k,:,:) = init%tt(k,:,:) +lapse*hgt0 + T_offset
    end do   
    
! #### END LOAD BOUNDARY DATA ####
  
    return
    
  end subroutine update_boundaries
  
  function deltaT_2D(dT,lat,lat_fac) result(dT2D)

    implicit none 

    integer, parameter :: nx = nxe, ny = nye
    
    double precision, intent(IN) :: dT, lat(ny,nx), lat_fac 
    double precision :: dT2D(ny,nx) 

    dT2D = dT + lat_fac*(lat-72.6d0)/(76.d0-72.6d0)     ! NEEMup - GISP2 
    return 

  end function deltaT_2D

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  p p _ s c a l a r
  ! Author     :  Alex Robinson
  ! Purpose    :  Function to apply uncertainty scaling to precip
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function pp_scalar(dT)
  
    implicit none
    
    double precision :: pp_scalar, dT
    double precision :: ppfac1
    
    ppfac1 = ppfac
!     if (ppfac1 .eq. 0.d0) ppfac1 = smb_ppfac
    
    pp_scalar = 1.d0 + ppfac * dT
  
    return
    
  end function pp_scalar
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c l i m c h e c k e r 2
  ! Author     :  Alex Robinson
  ! Purpose    :  Output values to check climate with observations
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine climchecker2(m20,zs0)

    use projector 
    
    implicit none
    
    character (len=30) :: units(15)
    
    type obsout
      real (8) :: tmean, pmean, smean, as, dp0, dp1, dp10
      real (8) :: dsnowf, dsnowc, snowhmean 
    end type
    
    type(obsout) :: out(13)
    
    integer,  parameter :: nx = nxs, ny = nys
    real (8), parameter :: dx = dxs
    
    real (8), dimension(ny,nx) :: m20, m2, zs0, zs, mobs, id
    character (len=256) :: fnm
    
    integer :: i, j, k, q, k0, k1
    
    real (8) :: tmp(int(day_month))
    
    ! Save mask to local grid
    m2 = m20
    zs = max(zs0,-0.1)    ! Make sure elevation doesn't include water [m]
    
    ! Store names of each variable being output
!!    units = (/ "id ", "lat", "lon", "zs", "month", &
!!               "tmean", "pmean", "smean", "as",  "dp0.1", "dp1.0", "dp10.0", &
!!               "dsnowf", "dsnowc", "snowhmean" /)
                
    if (init_obs .eq. 0) then
    
      ! Open the station file
      open(10,file="input/observations/stations.txt",status="old")
      
      read(10,*)    ! Skip the first line
      
      ! Read in all station information
      do q = 1, n_stations
        read(10,"(2x,a4,i9,a33,2f11.4,2f11.2)")  &
                                     stations(q)%src,  stations(q)%id,  &
                                     stations(q)%name, stations(q)%lat, &
                                     stations(q)%lon,  stations(q)%zs,  &
                                     stations(q)%StartDate
        
        ! From lon/lat, get x/y
        call plane_ellipsoid(stations(q)%lon*torads, stations(q)%lat*torads, &
                             stations(q)%x,   stations(q)%y, 1)
        
        stations(q)%x = stations(q)%x *1e-3
        stations(q)%y = stations(q)%y *1e-3
        
        ! Determine the indices on the sico grid
        call xy2index(sico_grid, stations(q)%x, stations(q)%y,  &
                                 stations(q)%i, stations(q)%j )
       
      end do
      
      ! Close stations file
      close(10)
      
      ! Reopen the file, and read the id, but into a string!!
      open(10,file="input/observations/stations.txt",status="old")
      read(10,*)    ! Skip the first line
      do q = 1, n_stations
        read(10,"(2x,a4,a9)")  stations(q)%src,  stations(q)%idc
        stations(q)%idc = adjustl(stations(q)%idc)
      end do         
      close(10)
      
      init_obs = 1
      
    end if
    ! end initializing stations
    
    
    ! Go through stations and write out all observations to each one
    do q = 1, n_stations
      
      i = stations(q)%i
      j = stations(q)%j
      
      out%tmean     = (/ mon%tt(j,i), ann%tt(j,i) /)
      out%pmean     = (/ mon%pp(j,i), ann%pp(j,i) /)
      out%smean     = (/ mon%snow(j,i), ann%snow(j,i) /)
      out%as        = (/ mon%as(j,i), ann%as(j,i) /)
      out%snowhmean = (/ mon%snowh(j,i), ann%snowh(j,i) /)
      
      ! Make monthly calculations where not available already
      do k = 1, 12
        k0 = (k-1)*day_month+1
        k1 = k0 + day_month-1
        
        tmp = 0.d0
        where (day(k0:k1)%pp(j,i) .ge. 0.1d0) tmp = 1.d0
        out(k)%dp0 = sum(tmp)
        
        tmp = 0.d0
        where (day(k0:k1)%pp(j,i) .ge. 1.0d0) tmp = 1.d0
        out(k)%dp1 = sum(tmp)
        
        tmp = 0.d0
        where (day(k0:k1)%pp(j,i) .ge. 10.0d0) tmp = 1.d0
        out(k)%dp10 = sum(tmp)
        
        tmp = 0.d0
        where (day(k0:k1)%snowh(j,i) .gt. 0.0d0) tmp = 1.d0
        out(k)%dsnowf = sum(tmp)
        
        tmp = 0.d0
        where (day(k0:k1)%snowh(j,i) .gt. 0.0d0) tmp = 1.d0
        out(k)%dsnowc = sum(tmp)
          
      end do
      
      ! Make the annual value (sums of days)
      k = 13
      out(k)%dp0       = sum(out(1:12)%dp0)  
      out(k)%dp1       = sum(out(1:12)%dp1)    
      out(k)%dp10      = sum(out(1:12)%dp10)  
      out(k)%dsnowf    = sum(out(1:12)%dsnowf)
      out(k)%dsnowc    = sum(out(1:12)%dsnowc)

      
      ! ## Now, output the results to files ##
      !    (one file per station)

      fnm = trim(stations(q)%idc)//".ltm"

      open(unit=10,file=trim(outfldr)//"obs/"//trim(fnm),status="unknown")
      
      ! Write the header
      write(10,"(a,1x,a1,1x,a)")     "ID","=",trim(stations(q)%idc)
      write(10,"(a,1x,a1,1x,a)")     "Name","=",trim(stations(q)%name)
      write(10,"(a,1x,a1,1x,f11.4)") "Latitude", "=",stations(q)%lat
      write(10,"(a,1x,a1,1x,f11.4)") "Longitude","=",stations(q)%lon
      write(10,"(a,1x,a1,1x,f11.4)") "Elevation","=",stations(q)%zs
      write(10,"(a,1x,a1,1x,f11.4)") "StartDate","=",stations(q)%StartDate

      ! Write the header of the output
      !!write(10,"(15a12)") (trim(units(k)), k=1,size(units))
      write(10,"(15a12)") "id ", "lat", "lon", "zs", "month", &
          "tmean", "pmean", "smean", "as",  "dp0.1", "dp1.0", &
          "dp10.0","dsnowf", "dsnowc", "snowhmean"
      
      ! Write the montly output of each variable
      do k = 1, 13
        write(10,"(a12,3f12.3,i12,11f12.3)") &
            trim(stations(q)%idc), stations(q)%lat,      &
            stations(q)%lon, zs(j,i), k,                 &
            out(k)%tmean,  out(k)%pmean,  out(k)%smean,  &
            out(k)%as,     out(k)%dp0,    out(k)%dp1,    &
            out(k)%dp10,   out(k)%dsnowf, out(k)%dsnowc, &
            out(k)%snowhmean 
      end do

      close(10)
      
    end do
    
    
    return
    
  end subroutine climchecker2
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  e m b _ t e m p e r
  ! Author     :  Alex Robinson
  ! Purpose    :  Determine the temperature over land through 
  !               diffusion, relaxing to data values over ocean
  ! dt [days]
  ! dxl [m]
  ! pp  [kg/m2-s]
  ! tt, t0  [°C]   - Sea-level T
  ! rhum [0:1]
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine emb_temper(tt,t0,rain,snow,S,ap,dh,rrco2,zs,mr,dxl,dt_sec,kappa)
    
    implicit none
    
    integer, parameter :: nx = nxe, ny = nye
    
    real (8) :: gama, beta, dxl, dt_sec, chck
    !real (8), parameter :: ce      = tce
    !real (8), parameter :: A       = ta      ! W / m2,    reflective rad. parameter, intercept   
    !real (8), parameter :: B       = tb      ! W / m2-K,  reflective rad. parameter, slope 
    !real (8), parameter :: L       = Lw       ! J / kg,    Specific latent heat flux of water in atmos... 
    !real (8), parameter :: lapse   = tlfac      ! °C / m,  atmospheric lapse rate
    !real (8), parameter :: rfac    = trfac    ! W / m2-K,  relaxation coeff.
    real (8), dimension(ny,nx) :: alpha, F, relax
    real (8), dimension(:,:) :: t0, tt, S, ap, dh, rrco2, zs, mr, kappa
    real (8), dimension(:,:) :: rain, snow
  
    integer :: i, j
    
    ! Constants, in time units of sec
    alpha = kappa * dt_sec / (dxl*dxl*tce)
    gama = dt_sec / tce 
    beta = tb

    ! Scale the relaxation mask
    relax = mr * trfac
        
    ! Get temperature forcing 
    F =  S*(1-ap) + rrco2 + (Lw*rain + Ls*snow) &
         - (ta - tb*tlfac*zs) - dh*rho_w*Lm

    if (solver .eq. 1) then
      call adi(t0,tt,F,alpha,gama,beta,relax)
    else      
      call explicit(t0,tt,F,alpha,gama,beta,relax)
    end if

    return
    
  end subroutine emb_temper

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  e m b _ p r e c i p
  ! Author     :  Alex Robinson
  ! Purpose    :  Determine an amount of precipitation based
  !               on the water content in the atmosphere for current
  !               temperature
  ! dt [days]
  ! dxl [m]
  ! P  [kg/m2-s]
  ! T  [°C]   - 'sea-level' temperature, to be corrected by lapse
  ! rhum [0:1]
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine emb_precip(pp,tt,qq,lapse,rhum,zs,dzs,mr,dxl,dt_sec,kappa,qqfac,ptau)
    
    implicit none
    
    integer,  parameter :: nx = nxe, ny = nye
    real (8), parameter :: dx = dxe

    integer :: i, j
    real (8) :: dxl, dt_sec, lapse, gama, chck
    !real (8), parameter :: ce    = pce
    !real (8), parameter :: rfac = prfac
    
    !real (8), parameter :: h_e = 2000.d0 ! meters elevation
    !real (8), parameter :: tau = p_tau * sec_day ! usually 5 days
    !real (8), parameter :: k = p_k
    
    real (8) :: h_e
    
    real (8), dimension(:,:) :: pp, tt, qq, rhum, zs, dzs, mr, kappa, qqfac, ptau
    real (8), dimension(ny,nx) :: p0, q0, q_sat, relax, alpha
    real (8), dimension(ny,nx) :: tmp, fr, tau
    
    ! Constants, in time units of sec
    alpha = kappa*dt_sec/(dx*dx*pce)
    gama = dt_sec / pce   
    
    tau = ptau * sec_day ! usually 5 days
    
    ! Characteristic elevation (usually 2000m)
    h_e = p_he
    
    ! Scale the relaxation mask
    relax = mr * prfac
    
    ! scale the temperature to surface T, 
    ! and increase by scaling factor where no relaxation occurs
    tmp = tt - lapse*zs + p_scaleT*(1.d0-mr)       
         
    do i = 1, nx
      do j = 1, ny
        q_sat(j,i) = emb_FQSAT(tmp(j,i))  &
        * h_e * airdens(zs(j,i))     ! rel. humidity is 1 for satuation
      end do
    end do
            
    ! Get current Q based on relative humidity
    q0 = q_sat*rhum
    
    ! Adjust by qfactor (only increases precip outside Greenland)
    q0 = q0 * qqfac
    
    ! Use diffusion to determine the new Q, forced by current P 
    if (solver .eq. 1) then     
      call adi(q0,qq,-pp,alpha,gama,0.d0,relax)
    else
      call explicit(q0,qq,-pp,alpha,gama,0.d0,relax)
    end if
       
    ! Hack in case q comes out negative (which shouldn't happen, but does in some cases!)
    ! ajr: find out why it goes negative, emb_PRECIP=1,emb_TEMPER=2
    ! ...or if that happens anymore?? ajr, 10.01.2009
    qq = max(qq,0.d0)
    
    ! Using the new value of Q, determine P
    ! Precip is either the total Q spread over eg tau=5days or
    ! the amount of Q over the saturation level
    !pp = max(qq/tau,(qq - q_sat)/dt_sec) 
    
    ! Only p1
    !pp = qq/tau
    
    ! Only p1b
    !n = 4.d0
    !fr = (qq/q_sat)**n
    !pp = (qq/tau)* fr
    
    ! Only p1c
    !fr = 1.d0 + p_k*dzs   !! now this is calculated outside in subroutine precip_grad...
    pp = (qq/tau) * (1.d0 + p_k*dzs)
            
    return
    
  end subroutine emb_precip    
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  p r e c i p _ g r a d
  ! Author     :  Alex Robinson
  ! Purpose    :  Horizontal gradient of u0, scaled for precip
  !               Returns gradient in km / km or m / m, depends on input
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
  subroutine precip_grad(u0,dx,lats,uu,klat,kfac)
  
    implicit none
    
    !integer, parameter :: nx = nxe, ny = nye
    
    integer :: i, j, nx, ny
    real (8), dimension(:,:), intent(IN) :: u0, lats
    real (8), dimension(:,:) :: uu
    real (8), dimension(:,:), allocatable :: uu_dir, uu_lat
    real (8) :: dx, inv_2dx, dudx, dudy
    real (8) :: klat, kfac, lmin, lmax
    
    nx = size(u0,2)
    ny = size(u0,1)
    
    allocate(uu_dir(ny,nx),uu_lat(ny,nx))

    lmin = minval(lats); lmax = maxval(lats)
    
    inv_2dx = 1.d0 / (2.d0 * dx)
    
    ! Assign boundary values
    uu(1:ny,1)  = 0.d0
    uu(1:ny,nx) = 0.d0
    uu(1,1:nx)  = 0.d0
    uu(ny,1:nx) = 0.d0
    
    ! Assign fraction as 1, assuming all points are East-facing
    uu_dir = 1.d0
    
    do i = 2, nx-1
      do j = 2, ny-1
      
        !uu(j,i) = dsqrt( inv_4dx2 * (u0(j,i+1) - u0(j,i-1))**2 + &
        !                 (u0(j+1,i) - u0(j-1,i))**2 )
        
        dudx = (u0(j,i+1) - u0(j,i-1)) * inv_2dx
        dudy = (u0(j+1,i) - u0(j-1,i)) * inv_2dx
        
        uu(j,i) = dsqrt( dudx**2 + dudy**2 )
        
        ! Correct the fraction for western facing points
        if (dudx .gt. 0) uu_dir(j,i) = kfac
        
      end do
    end do
    
    ! Apply an additional factor for latitude (same as in kappalat)
    ! Make a scale that goes from 1 at low lat to the fraction desired at top of grid
    ! kappalat (-1:1), negative means lower latitudes have higher kappa,
    ! positive kappa lat means higher latitudes have higher kappa
    uu_lat = 1.d0 + klat*(lats-lmin) 
    
    ! Now calculate the precipitation effect due to the gradient of topography
    ! scaled by the direction of the gradient and a latitudal factor
    uu = uu * uu_dir * uu_lat
    
    return
  
  end subroutine precip_grad
  
  
  subroutine kappa(kapfac,day,lats,lons,zs,mask,klat,klon,kamp,kzs,kdir)
    
    implicit none
      
      integer :: nx, ny, i, j
      integer  :: day,f_lat,f_season,f_zs, f_lon
      real (8) :: kapfac(:,:), lats(:,:), lons(:,:), zs(:,:), mask(:,:)
      real (8) :: latmin, latmax, lonmin, lonmax
      real (8) :: klat,klon,kamp,kzs,kdir
      
      real (8), dimension(:,:), allocatable :: uu_dir
      real (8) :: dudx, dudy, inv_2dx, dx
      
      nx = size(kapfac,2)
      ny = size(kapfac,1)
      
      allocate(uu_dir(ny,nx))
      
      dx = 20e3
      inv_2dx = 1.d0 / (2.d0 * dx)
      
!       klat      = kappalat
!       klon      = kappalon
!       kamp      = kappaamp   ! dimensionless
!       kzs1      = kappazs    ! dimensionless
!       if (present(kzs)) kzs1 = kzs
      
      latmin = minval(lats); latmax = maxval(lats)
      lonmin = minval(lons); lonmax = maxval(lons)
      
      kapfac = 1.d0

      ! Adjust kappa with latitude
      ! Make a scale that goes from 1 at low lat to the fraction desired at top of grid
      ! kappalat (-1:1), negative means lower latitudes have higher kappa,
      ! positive kappa lat means higher latitudes have higher kappa
      kapfac = 1.d0 + klat*(2.d0*(lats-latmin)/(latmax-latmin) - 1.d0) 
      
      ! Adjust kappa with latitude
      ! Make a scale that goes from 1 at low lat to the fraction desired at top of grid
      ! kappalat (-1:1), negative means lower latitudes have higher kappa,
      ! positive kappa lat means higher latitudes have higher kappa
      kapfac = kapfac * ( 1.d0 + klon*(2.d0*(lons-lonmin)/(lonmax-lonmin) - 1.d0) )
      
      
      ! Adjust kappa seasonally
      kapfac = kapfac * (1.d0 + kamp*dcos((day-15)*2.d0*pi/day_year))

      ! Adjust kappa with elevation
      kapfac = kapfac* (1.d0 + kzs*zs)
      
      ! Adjust kappa with direction of elevation
      ! Assign fraction as 0, assuming all points are East-facing
      uu_dir = 0.d0
      
      do i = 2, nx-1
        do j = 2, ny-1
        
          dudx = (zs(j,i+1) - zs(j,i-1)) * inv_2dx
          dudy = (zs(j+1,i) - zs(j-1,i)) * inv_2dx
          
          ! Correct the fraction for western facing points
          if (dudx .gt. 0) uu_dir(j,i) = kdir
          
        end do
      end do
      
      where ( mask .ne. 2.d0 ) kapfac = kapfac* (1.d0 + uu_dir)
  
    return
    
  end subroutine kappa
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : e f f e c t i v e T
  ! Author     : Reinhard Calov
  ! Purpose    : Computation of the positive degree days (PDD) with
  !              statistical temperature fluctuations;
  !              based on semi-analytical solution by Reinhard Calov.
  !              This subroutine uses days as time unit, each day
  !              is added individually
  !              (the same sigma as for pdd monthly can be used)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function effectiveT(temp, sigma)
    
    implicit none
    
    integer, parameter :: ny = nys, nx = nxs
    
    real (8), intent(in), dimension(:,:) :: temp
    real (8), dimension(ny,nx) :: effectiveT
    real (8) :: sigma

    real (8) :: inv_sigma

    real (8), parameter :: inv_sqrt2   = 1.d0/dsqrt(2.d0)
    real (8), parameter :: inv_sqrt2pi = 1.d0/dsqrt(2.d0*pi)
    
    inv_sigma   = 1.d0/sigma
    
    effectiveT = sigma*inv_sqrt2pi*dexp(-0.5d0*(temp*inv_sigma)**2)  &
                 + temp*0.5d0*erfcc(-temp*inv_sigma*inv_sqrt2)
    
    ! Result is the assumed/felt/effective positive degrees, 
    ! given the actual temperature (accounting for fluctuations in day/month/etc, 
    ! based on the sigma chosen)
                 
    return
  
  contains
  
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Function :  e r f c c
    ! Author   :  Reinhard Calov and Ralf Greve
    ! Purpose  :  Returns the complementary error function erfc(x) with 
    !             fractional error everywhere less than 1.2 x 10^(-7).
    !             Credit: Press et al., 'Numerical recipes in Fortran 77'.
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    function erfcc(x)

      implicit none
            
      real (8), intent(in), dimension(:,:) :: x
      real (8), dimension(ny,nx) :: erfcc

      real (8), dimension(ny,nx) :: t, z

      z = abs(x)
      t = 1.d0/(1.d0+0.5d0*z)

      erfcc = t*exp(-z*z-1.26551223d0+t*(1.00002368d0+t*(0.37409196d0+  &
      t*(0.09678418d0+t*(-0.18628806d0+t*(0.27886807d0+t*  &
      (-1.13520398d0+t*(1.48851587d0+t*(-0.82215223d0+  &
      t*0.17087277d0)))))))))

      where (x .lt. 0.d0) erfcc = 2.d0-erfcc

      return
      
    end function erfcc
    
  end function effectiveT
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s n o w f r a c
  ! Author   :  Reinhard Calov and Ralf Greve
  ! Purpose  :  Returns the fraction of precip that would fall as snow
  !             depending on the temperature (from 0 to 1)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function snowfrac(T)
    
    implicit none
    
    !integer, parameter :: nx = nxs, ny = nys
    integer :: ny, nx
    
    real (8), dimension(:,:) :: T
    real (8), dimension(:,:), allocatable :: snowfrac
    real (8) :: Tmin, Tmax
    
    ny = size(T,1)
    nx = size(T,2)
    allocate(snowfrac(ny,nx))
    
    Tmin = snow_Tmin                  ! min °C for snow 
    Tmax = snow_Tmax                  ! max °C for snow

    where (T .lt. Tmin)
      ! All precip is snow
      snowfrac = 1.d0
    else where (T .ge. Tmax)
      ! All precip is water - will not accumulate on ice sheet
      snowfrac = 0.d0
    elsewhere
      ! A fraction of precip is snow
      snowfrac = 0.5d0 * (cos(pi*(T-Tmin)/(Tmax-Tmin)) + 1.d0)
    end where
    
    return
    
  end function snowfrac 
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  p r e c i p 0
  ! Author   :  Alex Robinson
  ! Purpose  :  Calculate the precipitation based on a default field
  !             and anomalies
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine precip0(mask,zs,lats,lons,T_warming)

    implicit none
    
    integer, parameter :: nx = nxs, ny = nys
    double precision, parameter :: dx = dxs
    
    double precision, dimension(ny,nx) :: mask, zs, lats, lons, aco2 
    double precision, dimension(ny,nx) :: gamma, eta, dzs0, dzs
    double precision, dimension(ny,nx) :: pscalar
    double precision, dimension(nye,nxe) :: dzsl
    double precision :: sigma, sfac
    double precision :: T_warming          ! Summer warming anomaly

    integer :: i, j, k, m, km, km2

    character(len=256) :: fnm
    
    ! Assign global
    sigma = Teff_sigma
    
    ! Scale precipitation based on the boundary temperature anomaly
      
    ! Letreguilly, Huybrechts, 1991
!!    sfac   = 0.0533d0        
!!    pscalar = (1+sfac)**min(0.d0,T_warming_now)
      
    ! Huybrechts, DeWolde, 1999
! !     if (T_warming_now .le. -10.d0) then
! !       sfac = 0.10d0
! !     else if (T_warming_now .le. 0.d0) then
! !       sfac = 0.05d0 - 0.005d0 * T_warming_now
! !     else
! !       sfac = 0.05d0
! !     end if
! !     pscalar = (1.d0+sfac)**T_warming_now
    
    ! #### Huybrechts, 2002, QSR ####
    ! ## TO BE IMPLEMENTED... needs d18O record as input... ##
    
    ! #### Tarasov & Peltier, 2003, JGR ####
    ! comment by Reinhard:
    ! #### eta_xy is missing in implemented parameterization
    ! #### of precip sensitivity to temperature change.
    ! #### Total precip ddoes not chnages strong enough with
    ! #### temperature. Work arount use value from Ritz et al (1997)
    ! #### for eta. Well. This valus is too small. Check own one ...
    
    ! Get elevation gradients
    dzs0 = fields0%dzs
    call hgrad(zs,dx,dzs) 
    
    ! Smooth them!
    call tolores(dzs0,  dzsl, rembo0%wtst,nrt,ratio)
    call tohires(dzsl,dzs0,ratio,0)
    
    call tolores(dzs,  dzsl, rembo0%wtst,nrt,ratio)
    call tohires(dzsl,dzs,ratio,0)
    
    eta   = 0.1d0
    gamma = (1.d0-(zs/4d3)) &
            *min(1.5d0, (dzs+1d-3)/(dzs0+1d-3)) &
            + zs/4d3
            
    pscalar = gamma * dexp(eta*T_warming)
    
    ! Apply the scaling factor to the data loaded using function emb_load_input
    ann%pp   = fields0%precip * pscalar
    
    ! Get daily pp values.
    do k = 1, nk
      day(k)%pp   = ann%pp   / day_year    ! mm/a => mm/d (kg/m2-d)  water equivalent
      day(k)%snow = day(k)%pp * snowfrac(day(k)%tt)
    end do
        
    write(*,"(a1,5x,a,f8.2,f8.3)") "e", &
                      "Scaled accumulation (dT, pscalar):", T_warming, pscalar
  
    return
    
  end subroutine precip0
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t e m p e r 0
  ! Author   :  Alex Robinson
  ! Purpose  :  Calculate the temperature based on a default field
  !             and anomalies
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine temper0(mask,zs,lats,lons,T_anomaly)

    implicit none
    
    integer, parameter :: nx = nxs, ny = nys
    
    double precision, dimension(:,:) :: mask, zs, lats, lons
    double precision, dimension(ny,nx) :: ZZ, tma, tmj, tampl
    double precision :: gma, gmj, cma, cmj, dma, dmj
    character (len=100) :: method
    
    double precision, intent(IN) :: T_anomaly(:)    ! Daily prescribed anomaly values 

    integer :: i, j, k, m, km, km1, km2
    double precision :: sigma
    
    ! Assign global
    sigma = Teff_sigma
    
    select case(temper)
      
      case(0)
        ! Temperature  = EISMINT parameterization
        ! *lats are in degrees
        method = "EISMINT"
        
        ! First get height of inversion, then
        ! set ZZ properly
        ZZ = 20.d0*(lats*todegs - 65.d0)    
        ZZ = max(zs,ZZ); ZZ = max(0.d0,ZZ)
        
        tmj = 30.38d0 - 0.006277d0 * ZZ - 0.3262d0 * (lats*todegs)
        tma = 49.13d0 - 0.007992d0 * ZZ - 0.7576d0 * (lats*todegs) 
        
        !dma = 49.13d0; gma = -7.992d-3; cma = -7.576d-1
        
        
        tampl = tmj - tma
        
      case(-1)
        ! Temperature  = Reinhard's parameterization
        ! * lats should be in radians
        method = "REINHARD"
        
        ZZ = max(0.d0,zs)
        tma=48.38d0 - 43.d0*lats - 0.008d0*ZZ
        tampl = -23.d0 + 31.d0*lats
      
      case(-2)
        ! Temperature  = Fausto parameterization
        ! *lats/lons should be in degrees
        method = "FAUSTO"

        ZZ = max(0.d0,zs)
        
        tmj = 14.70d0 - 0.005426d0 * ZZ - 0.1585d0 * (lats*todegs) &
              - 0.0518d0 * (lons*todegs)
              
        tma = 41.83d0 - 0.006309d0 * ZZ - 0.7189d0 * (lats*todegs) &
              - 0.0672d0 * (lons*todegs)
        
        tampl = tmj - tma
      
      case(-3)
        ! Reinhard's change for circular setup
        ! Temperature  = schematic parameterization
        method = "SETUP_FAST"

        ZZ = max(0.d0,zs)

        tma = - 0.008d0 * ZZ
        tampl = tempamp
        
      case DEFAULT
        
        write(*,*) "conventional: temper setting should be less than zero!"
        write(*,*) "temper = ",temper
        stop
        
    end select
    
!     ! Get the temperature for each month
!     ann%tte = 0.d0
!     do m = 1, nm
!       k = int(m*day_month-15)  ! Make day equal to the middle of the month
!       mon(m)%tt  = tma - tampl*dcos(2.d0*pi*(k-15)/day_year) + deltaT(k) + T_offset
!       mon(m)%tte = effectiveT(mon(m)%tt,sigma) * sec_frac
!       ann%tte    = ann%tte + mon(m)%tte
!     end do                 
!     
!     ! Store the annual temperature
!     ann%tt  = tma + deltaT(0)
!     ann%tte = ann%tte / nm
    
    ! Get daily T and Teff values.
    do k = 1, nk
      day(k)%tt  = tma - tampl*dcos(2.d0*pi*(k-15)/day_year) + T_anomaly(k) + T_offset
      day(k)%tte = effectiveT(day(k)%tt,sigma) * sec_frac
    end do   
    
    ! Get the monthly mean temp and effective temp
    do m = 1, 12
        km1 = (m-1)*ndm+1
        km2 = km1+ndm
        mon(m)%tt  = 0.d0
        mon(m)%tte = 0.d0 
        do k = km1, km2
          mon(m)%tt  = mon(m)%tt  + day(k)%tt*(1.d0/dble(ndm))
          mon(m)%tte = mon(m)%tte + day(k)%tte
        end do 
    end do 

    ! Get the annual mean temp and effective temp 
    ann%tt  = 0.d0
    ann%tte = 0.d0 
    do k = 1, nk
      ann%tt  = ann%tt  + day(k)%tt*(1.d0/dble(nk)) 
      ann%tte = ann%tte + day(k)%tte 
    end do 

    write(*,"(a1,5x,i10,5x,a)") "e",nstep,"Calculated "//trim(method)//" temperature."
    
    return
    
  end subroutine temper0
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine :  p d d
    ! Author     : Alex Robinson
    ! Purpose    :  Computation of the positive degree days (PDD) with
    !               statistical temperature fluctuations;
    !               based on semi-analytical solution by Reinhard Calov.
    !               This subroutine uses year as time unit, each day
    !               is added individually
    !               (the same sigma as for pdd monthly can be used)
    !  * essentially identical to melt_budget, but only for use 
    !    as an annual approach
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine pdd(m2, tte, pp, snow, rf, runoff, runoff_snow,   &
                   runoff_rain, melt, melted_ice, new_ice, refrozen, dh )
      
      implicit none
      
      integer, parameter :: nx = nxs, ny = nys
      double precision, dimension(ny,nx) :: &
                m2, tte, pp, snow, rf, refrozen, runoff, melt, smb, &
                new_ice, smbi, dh
      double precision, dimension(ny,nx) :: &
                rain, pddsum, melted_snow, melted_ice, mmfac, mm_snow, &
                refrozen_snow, refrozen_rain, runoff_rain, runoff_snow, &
                refreezing_max, refrozen_snow_max, snow_to_ice
      
      ! Convert precip: [mm/a]
      rain   = pp - snow
      
      pddsum = tte * day_year    !   [a * deg C] => [D * deg C]
      
      mm_snow = pddsum * mm_teff_snow
      mmfac   = mm_teff_ice / mm_teff_snow
      
      where (mm_snow .gt. snow)
        
        ! All snow is melted, the rest of energy converted to melt some ice
        melted_snow = snow
        
        ! The rest of energy will go into melting ice
        melted_ice = (mm_snow - snow) * mmfac

      elsewhere
        
        ! Snow melt will use all energy, none left for ice melt
        melted_snow = mm_snow
        melted_ice  = 0.d0
        
      end where
      
      ! Determine how much snow will turn into ice (any snow that wasn't originally melted)
      snow_to_ice = snow - melted_snow
      
      ! Determine the refreezing factor (potential for refreezing)
      ! (like in JanssensHuybrechts2000)
      rf = 0.d0
      where (m2 .eq. 0.d0) rf = Pmaxfrac * snow / max(1d-5, (snow + rain))
      
      ! Determine the actual maximum amount of refreezing
      refreezing_max    = snow                                  ! Total refreezing depends on amount of snow left!
      refrozen_rain     = min(rain*rf,refreezing_max)           ! First rain takes up refreezing capacity
      refrozen_snow_max = refreezing_max - refrozen_rain        ! Subtract rain from refreezing capacity to determine what's left for snow
      refrozen_snow     = min(melted_snow*rf,refrozen_snow_max) ! melted_snow uses remaining capacity if it needs it
      refrozen          = refrozen_snow + refrozen_rain         ! Amount of ice created from refreezing
      
      runoff_snow = melted_snow - refrozen_snow                 ! Net snow melt
      runoff_rain = rain - refrozen_rain                        ! Net rain
      runoff      = runoff_snow + runoff_rain + melted_ice      ! Total runoff (includes melted ice)
      
      ! Get the total ice accumulated for the day
      new_ice = snow_to_ice + refrozen
               
      ! Get values of total ablation and smb (mm/a)!!
      melt = melted_snow + melted_ice
      smb  = snow + rain - runoff
      
      ! Other quantities (to be consistent with daily melt budget)
      smbi    = new_ice - melted_ice

      where ( m2 .eq. 0 ) 
        dh = melt-refrozen       ! Technically refrozen=0 here, but include it anyway
      elsewhere
        dh = melted_snow-refrozen
      end where
      
      return
    
    end subroutine pdd
    
    subroutine calc_rf(mask,tt,pp,snow,melt,rf,refrozen)
    
      implicit none
      
      double precision, parameter :: rho_wet = 960.d0, rho_dry = 300.d0   ! km/m3
      double precision, parameter :: rho_wet_dry = rho_wet/rho_dry
      double precision, parameter :: c_ice = 2.05d3  ! J/kg-K
      
      double precision, dimension(:,:) :: mask, tt, pp, snow, melt
      double precision, dimension(:,:) :: rf, refrozen
      double precision, dimension(nys,nxs) :: melted_snow
      
      melted_snow = min(snow,melt)
            
      rf = -min(0.d0,tt)*(c_ice/Lm)*(snow/pp)                    &
          + (rho_wet_dry-1.d0)*(snow-melted_snow)/pp
      
      rf = min( 1.d0, rf )
      
      where ( mask .ne. 0.d0 ) rf = 0.d0
      
      refrozen = rf*(melted_snow + (pp-snow))
      
      write(*,*) "min/max melted_snow: ",minval(melted_snow), maxval(melted_snow)
      write(*,*) "min/max melt: ",minval(melt), maxval(melt)
      write(*,*) "min/max snow-melted_snow: ",minval(snow-melted_snow),maxval(snow-melted_snow)
      write(*,*) "min/max rf: ",minval(rf), maxval(rf)
      
      if ( minval(rf) .lt. 0.d0 ) stop
      
      return
      
    end subroutine calc_rf
  
end module rembo_main
