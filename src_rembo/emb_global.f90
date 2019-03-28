! emb globals

module emb_global
  
  use parameters
  use exchange
  use ncio
  use ncio_transpose 

  implicit none

  ! Math constants
  real (8), parameter :: pi = 3.141592653589793d0
  real (8), parameter :: torads = pi/180.d0, todegs = 180.d0/pi
  double precision, parameter :: ebs = 1d-8
  
  ! Physical constants
  real (8), parameter :: rho_i   = 9.1d2       ! Density of ice
  real (8), parameter :: rho_w   = 1.d3        ! Density of pure water
  real (8), parameter :: rho_sw  = 1.28d3      ! Density of sea-water (phys_para_grl.dat)
  real (8), parameter :: rho_i_w = rho_i / rho_w
  
  real (8), parameter :: R_earth = 6394800.d0
  real (8), parameter :: a_earth = 6378137.d0, b_earth = 6356752.3142d0
  ! Grid: reference for stereographic projection
  real (8), parameter :: lambda0 = -39.d0 * torads, &
                         phi0    =  71.d0 * torads
                         
  ! Lookup table
  real (8), parameter :: Teff_min = -80.d0
  real (8), parameter :: Teff_max = 40.d0
  real (8), parameter :: dTeff = 0.5d0
  integer, parameter  :: nTeff = (Teff_max - Teff_min) / dTeff
  integer, parameter  :: nTeffa = nTeff / (Teff_max - Teff_min)
  real (8), dimension(2,nTeff) :: Teffs
  
  ! binary adjustment
  integer, parameter :: bin_fac = 4
  
  
  ! Dimensional variables

  ! **Note: later these values will come from the sico_global file, not here!
  !integer,  parameter :: nxs = 301, nys = 561     ! FROM SICOPOLIS
  !real (8), parameter :: dxs =  5d3               ! FROM SICOPOLIS
  !integer,  parameter :: nxs = 151, nys = 281     ! FROM SICOPOLIS
  !real (8), parameter :: dxs = 10d3               ! FROM SICOPOLIS
  !integer,  parameter :: nxs = 76, nys = 141       ! FROM SICOPOLIS
  !real (8), parameter :: dxs = 20d3                ! FROM SICOPOLIS
  integer,  parameter :: nxs = 85, nys = 145       ! FROM SICOPOLIS
  real (8), parameter :: dxs = 20d3                ! FROM SICOPOLIS

  ! (reinhard, for calov-island)
!  integer,  parameter :: nxs = 161, nys = 161      ! FROM SICOPOLIS
!  real (8), parameter :: dxs = 10d3                ! FROM SICOPOLIS

  ! Time step choices
  real (8) :: dte
  integer :: dtime_emb, dtime_bndry, dtime_smb
  integer :: dto_clim, dto_clim2d
  integer :: dto_timer
  
  real (8) :: year0, yearf, year_offset
  integer,  parameter :: nk  = 360
  integer,  parameter :: nm  = 12
  integer,  parameter :: ndm = 30
  
  ! Time conversions
  real (8), parameter :: day_year   = dble(nk)
  real (8), parameter :: day_month  = dble(ndm)
  real (8), parameter :: month_year = dble(nm)
  real (8), parameter :: sec_year   = 31556926.d0
  real (8), parameter :: sec_day    = sec_year / day_year   ! 8.765813d4
  real (8), parameter :: sec_day0   = 8.64d4
  real (8), parameter :: sec_frac   = sec_day / sec_day0
  
  ! ### Constants for resolution
  !real (8), parameter :: dxe = 50d3
  !integer,  parameter :: nxe = 31, nye = 57
  !real (8), parameter :: dxe = 100d3
  !integer,  parameter :: nxe = 16, nye = 29
  real (8), parameter :: dxe = 80d3
  integer,  parameter :: nxe = 22, nye = 37
  real (8), parameter :: x0  = -720000.d0
  real (8), parameter :: y0  = -3450000.d0
  
  ! (reinhard, for calov island)
!  real (8), parameter :: dxe = 100d3
!  integer,  parameter :: nxe = 17, nye = 17
!  real (8), parameter :: x0  = -800000.d0
!  real (8), parameter :: y0  = -800000.d0
  
  type location
    double precision :: x, y, lat, lon
  end type
  
  type(location), dimension(nys,nxs) :: sico_grid
  type(location), dimension(nye,nxe) :: emb_grid
  
  ! Smoothing radii & neighborhoods
  real (8) :: prad, trad
  integer :: nrt, nrp
  
  real (8), parameter :: ratio = dxs / dxe
  
  
  ! ### Constants for temperature diffusion
  real (8) :: tce, tkappa, tb, ta, trfac, tlfac, s0
  real (8) :: T_offset, T_warming, T_wintfac, Teff_sigma
  real (8) :: T_noise, T_noise_period, clim_sens, firn_factor, dT_factor, dT_min, dT_width, lat_grad
  real (8) :: T_warming_delay, T_trans_max, T_diff, dT_rate
  real (8) :: f_eem, f_hol, f_seas, f_glac 

  ! (reinhard)
  real (8) :: tempamp
  
  ! ### Constants for precipitation diffusion
  real (8) :: pce, pkappa, prfac, p_scaleT, p_k, p_k_eastfrac, p_k_lat
  real (8) :: p_tau, p_he
  
  ! ### Constant for both
  real (8) :: kappaamp, kappalat, kappalon, kappazs
  real (8) :: ap0_intercept, ap0_slope
  real (8) :: ap_snow, ap_ice, ap_land
  real (8) :: as_snow0, as_snow1, as_snow2, as_snow_forest, as_land, as_land_forest
  real (8) :: as_ice, as_ocean
  
  ! ### Constants for melt budget calculation
  real (8) :: Cp, Lw, Lm, Ls, T_melt, h_snow_max, melt_crit
  real (8) :: hsnow_crit0, hsnow_crit1
  real (8) :: mm_teff_snow, mm_teff_ice, pdd_factor, pdd_a, pdd_b 
  real (8) :: pdd_Tmax, pdd_Tmin, pdd_spread, pdd_amax
  real (8) :: Pmaxfrac
  real (8) :: itm_c, itm_b, itm_t
  real (8), dimension(:,:), allocatable :: itm_cc
  real (8) :: at_intercept, at_slope
  
  ! ### Snow fraction variables
  real (8) :: snow_Tmin, snow_Tmax

  ! Program switches
  integer :: solver, equili, precip, temper, climchoice
  integer :: melt_choice, refreezing, clim_coupled
  integer :: co2_rad, ap_fixed, anf_dat, transient, boundary_forcing
  integer :: kill, slow_hyst, timetype, tuning, itm_S0
  
  ! Related to observation calculations
  integer, allocatable :: elem_yes(:)
  
  ! Input / output locations, switches
  character (len=256) :: domain, outfldr
  character (len=256) :: topo_file, restart_file, forcing_file
  character (len=256) :: emb_clim_file, clim_file
  integer :: stdout
  
  integer :: write_rembo_r, write_emb_d, write_rembo_m, write_rembo_d
  
  ! Boundary data options
  double precision :: bnd_yr0, bnd_start, bnd_ave
  integer :: bnd_trans, bnd_trans_delay
  double precision :: forcing_yr0, forcing_yrf 

  ! Equilibration time steps
  integer :: n_equili 
  
  type output
    real (8), dimension(13) :: tt, pp, tte, as, ap, snow, snowh, ice
    real (8), dimension(13) :: melt, smb, melt2, smb2, ma, aa, frozen
    real (8), dimension(13) :: S, S65
  end type
  
  ! EXTRA : smb case
  double precision :: smb_ppfac, smb_dT
  double precision :: ppfac
  
  ! EXTRA : slow hyst
  double precision :: h_dVdt_max, h_dTdt_min, h_dTdt_max, h_fac
  integer :: h_dVdt_window
  
  ! EXTRA : paleo fractions
  double precision :: paleo_frac_dT
  
  ! Counting variables
  integer :: nout_pdd, nout_pddday, nout_pddmonth, nout_day
  integer :: nout_pdd_nc
  integer :: nstep, nstep_max, nstep_sclim_max
  
  integer :: n_out, n_out_r
  double precision, allocatable :: time_out_mon(:), time_out_day(:)
  double precision, pointer :: time_out(:), time_out_r(:)
  
  ! Initialization switches
  integer :: init_p, init_t, init_pddnew
  integer :: init_rembo, init_summary, init_summary2, init_obs, init_paleo
  integer :: init_sinsol
  
  ! Temporary globals for testing!!!
  real (8), dimension(:,:), allocatable :: melt_ice
  real (8), dimension(:,:), allocatable :: accum_in0, pp_in0
  
  type saved_type
    double precision :: precip, snow, melted_ice, runoff_snow, runoff_rain
    double precision :: refrozen, evap, smb, tts, h_snow
    double precision :: jh_rf, jh_refrozen, tt, tjja, tdjf, tjan, tjul, ttp, pdds, melt_days
    double precision :: pkappa, ptau
  end type
  
  type(saved_type), dimension(:,:), allocatable :: saved

  type saved_type_mon
    double precision, dimension(:,:,:), allocatable :: tt, pdd_corr
    double precision, dimension(:,:,:), allocatable :: melt, as, S, S0
  end type 

  type(saved_type_mon) :: savedm 

  type input_fields
    double precision :: m2, zs, lats, lons
    double precision :: mask_hydro
    double precision :: accum, precip
    double precision :: tt, snowh, dh, ap, pdds
    double precision :: pkappa, qqfac
    double precision :: dzs
  end type
  
  type(input_fields), dimension(:,:), allocatable :: fields0
  
  type transient_warming_type
    double precision :: T_warming, dTdt, dVdt
    double precision, dimension(:), allocatable :: dVdtm
    double precision, dimension(:,:), allocatable :: zs
    integer :: nstep, mstep
  end type
  
  type(transient_warming_type) :: transT
  
  ! Timing  
  type timers
    real (8) :: boundary, climate, emb, pddnew
    real (8) :: pddold, now, step, tmp
  end type
  
  type(timers) :: timer, timerave, timern
  
  ! Observation comparison variables
  type station
    character (len=100) :: idc, src, name
    real (8) :: lat, lon, zs, StartDate, x, y 
    integer  :: id, i, j
  end type
  
  integer, parameter :: n_stations = 53
  type(station) :: stations(n_stations)
  
  ! Cost optimization output vector
  ! length = n_stations * 13months * (2:tmean, pmean) + (2:accum%a,melt%a)
  real (8), dimension(n_stations*13*2+2) :: fit
  
  type boundary_forcing_type
    double precision :: time, co2e, co2, rsl, dT(12), dTbnd(12), dTamp
  end type
  
  type(boundary_forcing_type), allocatable, dimension(:) :: forcing
  type(boundary_forcing_type) :: forcing_now

  double precision :: Tanomaly(nk)
  
contains  
    
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  m a k e _ g r i d
  ! Author     :  Alex Robinson
  ! Purpose    :  Fill in a grid with x, y values based on 
  !               starting lower-left corner and dx, dy
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine make_grid(grid,dx,x0,y0,nx,ny)
  
    implicit none
    
    integer :: i, j, nx, ny
    
    real (8) :: x0, y0, dx
    
    type(location) :: grid(:,:)
    
    do i = 1, nx
      grid(:,i)%x = x0 + (i-1)*dx
    end do
    
    do j = 1, ny
      grid(j,:)%y = y0 + (j-1)*dx    
    end do
    
    ! Convert to kilometers
    grid%x = grid%x *1e-3
    grid%y = grid%y *1e-3
    
    return
  
  end subroutine make_grid

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  e m b _ g l o b i n i t
  ! Author     :  Alex Robinson
  ! Purpose    :  Assign values to global variables 
  !               located in emb_global.f90
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine emb_globinit(stat)
    
    implicit none
    
    integer :: stat
    
    integer,  parameter :: nx = nxs, ny = nys
    character (len=10) :: c_dx
    character (len=100) :: folder, tmp
    integer :: i, k, k0, k1
    real (8), parameter :: dx = dxs
    real (8) :: latmid
    
    double precision, allocatable, dimension(:) :: tmpt
    
    if (stat .eq. 1) then
      
      ! Get resolution suffix (10km, 20km, etc)
      if (int(dx/1d3) .lt. 100) then
        write(c_dx,"(i2)") int(dx/1d3)
      else
        write(c_dx,"(i3)") int(dx/1d3)
      end if

      ! Get the command line arguments defining output folder
      call args(outfldr)      
      
      ! Initialize emb txt output file (needs outfldr)
      ! Set the emb output file unit number
      stdout = 334
      call embout(0,trim(outfldr))
      
      ! Get the parameter choices from the options file,
      ! and save these values to the global parameter array
      call get_params(trim(outfldr)//"options_rembo",glob=1)
        
      ! Write the parameter list to a file
      call out_params(trim(outfldr)//"rembo.params")
      
      domain              = cparam("domain")
      restart_file        = cparam("restart_file")
      topo_file           = cparam("topo_file")
      forcing_file        = cparam("forcing_file")
      emb_clim_file       = cparam("emb_clim_file")
      clim_file           = cparam("clim_file")
      
      !! Program switches
      climchoice                  =  int(param("climchoice"))
      precip                      =  int(param("precip"))
      temper                      =  int(param("temper"))
      equili                      =  int(param("equili"))
      solver                      =  int(param("solver"))
      melt_choice                 =  int(param("melt_choice"))
      refreezing                  =  int(param("refreezing"))
      co2_rad                     =  int(param("co2_rad"))
      ap_fixed                    =  int(param("ap_fixed"))
      boundary_forcing            =  int(param("boundary_forcing"))
      anf_dat                     =  int(param("anf_dat"))
      clim_coupled                =  int(param("clim_coupled"))
      transient                   =  int(param("transient"))
      kill                        =  int(param("kill"))
      slow_hyst                   =  int(param("slow_hyst"))
      timetype                    =  int(param("timetype"))
      tuning                      =  int(param("tuning"))
      itm_S0                      =  int(param("itm_S0"))

      !! Time variables
      dte                         =  param("dte")
      dtime_emb                   =  int(param("dtime_emb"))
      dtime_smb                   =  int(param("dtime_smb"))
      dto_clim                    =  int(param("dto_clim"))
      dto_clim2d                  =  int(param("dto_clim2d"))
      dto_timer                   =  int(param("dto_timer"))
      n_equili                    =  int(param("n_equili"))
      
      !! Writing switches
      write_emb_d                 =  int(param("write_emb_d"))
      write_rembo_m               =  int(param("write_rembo_m"))
      write_rembo_d               =  int(param("write_rembo_d"))
      write_rembo_r               =  int(param("write_rembo_r"))
      
      !! Time options
      year0                       =  param("year0")
      yearf                       =  param("yearf")
      year_offset                 =  param("year_offset")
      
      !! Boundary data options
      bnd_yr0                     =  param("bnd_yr0")
      bnd_start                   =  param("bnd_start")
      bnd_ave                     =  param("bnd_ave")
      bnd_trans                   =  int(param("bnd_trans"))
      bnd_trans_delay             =  int(param("bnd_trans_delay"))
      forcing_yr0                 =  param("forcing_yr0")
      forcing_yrf                 =  param("forcing_yrf")

      !! EMB diffusion variables
      kappaamp                    =  param("kappaamp")
      kappalat                    =  param("kappalat")
      kappalon                    =  param("kappalon")
      kappazs                     =  param("kappazs")
       
      !! EMB diffusion variables, T
      tce                         =  param("tce")
      tkappa                      =  param("tkappa")
      tb                          =  param("tb")
      ta                          =  param("ta")
      trfac                       =  param("trfac")
      tlfac                       =  param("tlfac")
      s0                          =  param("s0")
        
      T_offset                    =  param("T_offset")
      T_warming                   =  param("T_warming")
      T_wintfac                   =  param("T_wintfac")
      T_noise                     =  param("T_noise")
      T_noise_period              =  param("T_noise_period")
      clim_sens                   =  param("clim_sens")
      firn_factor                 =  param("firn_factor")
      dT_factor                   =  param("dT_factor")
      dT_min                      =  param("dT_min")
      dT_width                    =  param("dT_width")
      lat_grad                    =  param("lat_grad")
      f_eem                       =  param("f_eem")
      f_hol                       =  param("f_hol")
      f_seas                      =  param("f_seas")
      f_glac                      =  param("f_glac")
      T_warming_delay             =  param("T_warming_delay")
      T_trans_max                 =  param("T_trans_max")
      T_diff                      =  param("T_diff")
      dT_rate                     =  param("dT_rate")

      ! (reinhard)
      tempamp                     =  param("tempamp")
      
      !! EMB diffusion variables, P
      pce                         =  param("pce")
      prfac                       =  param("prfac")
      p_scaleT                    =  param("p_scaleT")
      pkappa                      =  param("pkappa")
      p_k                         =  param("p_k")
      p_k_eastfrac                =  param("p_k_eastfrac")
      p_k_lat                     =  param("p_k_lat")
      p_tau                       =  param("p_tau")
      p_he                        =  param("p_he")
      ppfac                       =  param("ppfac")
      
      !! Snow fraction 
      snow_Tmin                   =  param("snow_Tmin")
      snow_Tmax                   =  param("snow_Tmax")

      !! Surface albedo
      as_snow0                    =  param("as_snow0")
      as_snow1                    =  param("as_snow1")
      as_snow2                    =  param("as_snow2")
      as_snow_forest              =  param("as_snow_forest")
      as_ice                      =  param("as_ice")
      as_land                     =  param("as_land")
      as_land_forest              =  param("as_land_forest")
      as_ocean                    =  param("as_ocean")
      hsnow_crit0                 =  param("hsnow_crit0")
      hsnow_crit1                 =  param("hsnow_crit1")
      melt_crit                   =  param("melt_crit")
      
      !! Planetary albedo, parameterization 1
      ap0_intercept               =  param("ap0_intercept")
      ap0_slope                   =  param("ap0_slope")
           
      !! Melt variables
      Teff_sigma                  =  param("Teff_sigma")
      mm_teff_snow                =  param("mm_teff_snow")
      mm_teff_ice                 =  param("mm_teff_ice")
      pdd_factor                  =  param("pdd_factor")
      pdd_Tmax                    =  param("pdd_Tmax")
      pdd_Tmin                    =  param("pdd_Tmin")
      pdd_spread                  =  param("pdd_spread")
      pdd_amax                    =  param("pdd_amax")
      Pmaxfrac                    =  param("Pmaxfrac")

      !! Oerleman's melt scheme
      itm_c                       =  param("itm_c")
      itm_b                       =  param("itm_b")
      itm_t                       =  param("itm_t")
      at_intercept                =  param("at_intercept")
      at_slope                    =  param("at_slope")
      
      !! Refreezing parameters (superimposed ice)
      Cp                          =  param("Cp")
      T_melt                      =  param("T_melt")
      h_snow_max                  =  param("h_snow_max")
        
      !! Smoothing radii and topography
      prad                        =  param("prad")
      trad                        =  param("trad")

      !! Physics
      Lw                          =  param("Lw")
      Lm                          =  param("Lm")
      Ls                          =  param("Ls")
      
      !!! EXTRA cases
      smb_ppfac                   =  param("smb_ppfac")
      smb_dT                      =  param("smb_dT")
      
      h_dVdt_max                  =  param("h_dVdt_max")
      h_dTdt_min                  =  param("h_dTdt_min")
      h_dTdt_max                  =  param("h_dTdt_max")
      h_dVdt_window               =  int(param("h_dVdt_window"))
      h_fac                       =  param("h_fac")
      
      paleo_frac_dT               =  param("paleo_frac_dT")
      
      !!! End extra cases
      
      
      ! * Once I know prad and trad, I can calculate nrt and nrp
      !   (the neighborhood radii)
      nrt = int(trad / dxs)
      nrp = int(prad / dxs)
        
      ! ## SPECIAL CASES ##
      
      ! If running 'calov-island', ensure parameters are correct
      if ( trim(domain) .eq. "calov-island" ) then
      
        kappalat = 0.d0
        kappalon = 0.d0
        
      end if
      
      ! If REMBO is in stand-alone mode, modify some variables
      if ( clim_coupled .lt. 0 ) write(*,*) "** REMBO only mode **"
      if ( tuning .eq. 1 ) write(*,*) "** TUNING **"
      if ( clim_coupled .eq. 0 ) write(*,*) "** SICOPOLIS only mode (no rembo) **"
      
      ! For consistency
      if ( clim_coupled .eq. 1 .and. climchoice .eq. 0 ) then
      
        ! When using the conventional forcing approach, make sure it
        ! updates every year that rembo is called...
        dtime_smb = 5
        write(*,"(a1,5x,a,a)") "e","climchoice=0 ==> dtime_smb = dt_ice (sicopolis timestep)"
        
      end if
      
      ! For consistency
      if (climchoice .eq. 0) then
        
        ! Make sure switches make sense
        melt_choice=0; equili=0
        if ( temper .gt. 0 ) temper = 0
        if ( precip .gt. 0 ) precip = 0
        
        write(*,"(a1,5x,a,a)") "e","climchoice=0 ==> temper,precip,melt_choice,equili=0"
        
      end if
      
      ! ## END SPECIAL CASES ##
      
      ! Adjust by the pdd factor ( result is still mmwe / degC )
      mm_teff_snow = mm_teff_snow * pdd_factor 
      mm_teff_ice  = mm_teff_ice  * pdd_factor
      
      ! Adjust h_snow_max relative to smb timestep
      !h_snow_max = h_snow_max * dt_smb / 10.d0
      
      ! observations
      allocate(elem_yes(3))

      allocate(saved(ny,nx) )
      allocate(savedm%tt(12,ny,nx),  savedm%pdd_corr(12,ny,nx))
      allocate(savedm%melt(12,ny,nx),savedm%as(12,ny,nx))
      allocate(savedm%S(12,ny,nx),   savedm%S0(12,ny,nx))

      allocate( itm_cc(ny,nx) )
      allocate( accum_in0(ny,nx), pp_in0(ny,nx), melt_ice(ny,nx) )
      
      allocate( fields0(ny,nx) )
      allocate( transT%zs(ny,nx) )
      allocate( transT%dVdtm(h_dVdt_window/dtime_smb) )
      
      !! Send pertinent information to the exchange
      call init_exchange(coupled=clim_coupled,ny=ny,nx=nx,     &
                         t0=year0,t1=yearf,dt_clim=dble(dtime_emb),  &
                         dt_surf=dble(dtime_smb))
      
      ! Get input fields
      call emb_load_input()
      
      ! Initialize the sicopolis resolution grid
      call make_grid(sico_grid,dxs,x0,y0,nx,ny)
      sico_grid%lat = fields0%lats
      sico_grid%lon = fields0%lons
      
      ! Initialize the emb grid
      call make_grid(emb_grid,dxe,x0,y0,nxe,nye)
      
      !! Modify the itm_cc field according to latitude
      latmid = (maxval(sico_grid%lat)+minval(sico_grid%lat)) / 2.d0
      itm_cc = itm_c + itm_b*(sico_grid%lat-latmid)
          
      ! Load 2d and restart output years
      call get_output_times("out2d",time_out,n_out,year0,yearf,dble(dto_clim2d))
      call get_output_times("restart",time_out_r,n_out_r,year0,yearf)
      
      ! Generate daily and monthly values for the hi-res output
      allocate(time_out_mon(nm),time_out_day(nk))
      do i = 1, nm; time_out_mon(i) = dble(i); end do
      do i = 1, nk; time_out_day(i) = dble(i); end do
      
      ! Make sure all intializations will occur (set flags to zero)
      init_rembo    = 0
      init_summary  = 0
      init_summary2 = 0
      init_p        = 0
      init_t        = 0
      init_obs      = 0
      init_paleo    = 0
      init_sinsol   = 0
      
      write(*,*)
      write(*,"(a1,5x,a,a)") "e","outfolder  = ",trim(outfldr)
      write(*,*)
      
    else if (stat .eq. -1) then
            
      ! Close emb txt output file
      call embout(-1)
    else
      
      write(stdout,*) "Incorrect allocation choice for emb_globinit. Try again!"
      call embout(1)
    
    end if
    
    
    return
  
  end subroutine emb_globinit
  
  subroutine emb_load_input()
    
    implicit none
    
    integer, parameter :: nx = nxs, ny = nys
    real (8), dimension(ny,nx) :: accum, precip
    real (8), dimension(ny,nx) :: m2, zs, zb, lats, lons, zb0
    real (8), dimension(ny,nx) :: mask_hydro
    real (8), dimension(ny,nx) :: H_ice 
    character (len=256) :: fnm
    real (8), parameter :: dx = dxs
    
    select case(trim(domain))  ! (reinhard, calov-island)
      case("calov-island")
      
        precip = 300.d0
        accum  = 300.d0
      
      case DEFAULT   ! GRL = greenland
        
        fnm = trim(clim_file)
        write(*,*) "fnm = ", trim(fnm)
        write(*,*) "size(precip): ", size(precip,1), size(precip,2)
        call nc_read_t(fnm,"pr",  precip)  ! mmwe/a
        !call nc_read_t(fnm,"sf",accum)     ! mmwe/a
        accum = precip
        
    end select
    
    ! Make sure accum is never higher than actual prepitation (for consistency)
    where ( accum .gt. precip ) accum = precip
    
    write(*,"(a1,5x,i10,5x,a16,a)") "e",nstep,"Loaded precip: ",trim(clim_file)
    write(*,"(a,2f12.2)") "        precip min/max: ",minval(precip),maxval(precip)
    write(*,"(a,2f12.2)") "  snow (accum) min/max: ",minval(accum),maxval(accum)
          
    fnm = trim(topo_file)
    
    ! Get default topo variables (present-day)
    !call nc_read_t(fnm,"mask",  m2)
    call nc_read_t(fnm,"lat2D", lats)
    call nc_read_t(fnm,"lon2D", lons)
    call nc_read_t(fnm,"z_srf", zs)
    call nc_read_t(fnm,"z_bed", zb)
    !call nc_read_t(fnm,"zb0", zb0)
    zb0 = zb 
    
    ! Set ice(0)/land(1)/ocean(2) mask
    call nc_read_t(fnm,"H_ice",H_ice)
    m2 = 0.0
    where(zs .gt. 0.0 .and. H_ice .eq. 0.0) m2 = 1.0 
    where(zs .le. 0.0 .and. H_ice .eq. 0.0) m2 = 2.0  

    ! Load the hydrological basin
    !call nc_read_t(fnm,"mask_hydro",mask_hydro)
    mask_hydro = 1.0 

    if ( anf_dat .eq. 3 ) then !! Load variables from restart file
      
      ! netcdf restart filename
      fnm = trim(outfldr)//trim(restart_file)      
      
      ! Get main topo variables
      call nc_read_t(fnm,"mask", m2)
      call nc_read_t(fnm,"zs",   zs)
      call nc_read_t(fnm,"zb",   zb)
      
      ! read climate variables first
      call nc_read_t(fnm,"tt",   fields0%tt)
      call nc_read_t(fnm,"snowh",fields0%snowh)
      call nc_read_t(fnm,"dh",   fields0%dh)
      call nc_read_t(fnm,"ap",   fields0%ap)
      call nc_read_t(fnm,"pdds", fields0%pdds)
      
      ! Turn off equilibration (since we are restarting!)
      equili = 0
      
    end if

    ! Convert lats/lons to rads
    lats = lats * torads
    lons = lons * torads
    
    if (anf_dat .eq. 2) then         ! ice-free, uplifted (relaxed) topography
      
      ! Set surface elevation equal to zb0
      zs = zb0; zb = zb0
      
      ! Set mask equal to ice-free mask
      where(m2 .eq. 0) m2 = 1      
    
    else if (anf_dat .eq. -2) then   ! ice-free, loaded (un-relaxed) topography
      
      ! Set surface elevation equal to zb
      zs = zb
      
      ! Set mask equal to ice-free mask
      where(m2 .eq. 0) m2 = 1
      
    else if (anf_dat > 10 .or. anf_dat < 10) then ! surface elevation adjusted by anf_dat meters
      
      ! Offset ice surface elevation
      where ( zs > zb ) zs = zs + anf_dat
      
      ! Correct the mask and surface elevation for consistency
      where ( zs < zb ) 
        zs = zb0
        zb = zb0
        m2 = 1
      end where

      ! Adjust the bedrock elevation to reflect new ice thickness
      zb = zb0 + 910.d0/3300.d0*(zs-zb)

    end if
          
    ! Store the loaded variables in the global input fields
    fields0%m2     = m2
    fields0%zs     = zs
    fields0%accum  = accum
    fields0%precip = precip
    fields0%lats   = lats
    fields0%lons   = lons
    
    fields0%mask_hydro = mask_hydro
        
    return
    
  end subroutine emb_load_input
  
    
end module emb_global

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  t i m i n g
  ! Author     :  Alex Robinson
  ! Purpose    :  Use this routine to keep track of different timers
  !               in the program and output the results to file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine timing(n_step,start,tot)
    
    use emb_global 

    implicit none
    
    integer :: n_step, check
    
    real (8) :: start, tot, now, prev
    real (8) :: step, boundary, climate, emb, pddnew, pddold
    
    real (8), parameter :: convert = 1.d0/60.d0     ! [sec] => [min]
    character (len=100) :: fnm
    
    call cpu_time(now)
    
    ! Open file for writing timing output
    if (n_step .eq. 0) then
      
      fnm=trim(outfldr)//"out.timing"
      
      open(329,iostat=check,file=trim(fnm),status="replace")
      if (check .ne. 1) then
        write(stdout,*) "error opening file: ",trim(fnm)
      end if

      write(329,"(a1,a14,7a15)") "#", "years","sec","sec","sec", "sec", "sec","sec","sec"
      write(329,"(8a15)") "time","tot","step","boundary", "climate", "emb","pddnew","pddold"
      
      ! Initialize the averages and counts to zero
      call timer_reset(timerave)
      call timer_reset(timern)
      call timer_reset(timer)
      
      tot = now - start
    
    else
        
      prev = tot                     ! Previous equals the total we had
      tot  = now - start             ! New total based on time now
      timer%step = tot - prev        ! Time between now and last step [sec]
      
      ! Add timer values to the global averages
      ! (global averages will be flushed every timerstep iterations)
      call get_t_ave(timerave%boundary,timer%boundary,timern%boundary)
      call get_t_ave(timerave%climate,timer%climate,timern%climate)
      call get_t_ave(timerave%emb,timer%emb,timern%emb)
      call get_t_ave(timerave%pddnew,timer%pddnew,timern%pddnew)
      call get_t_ave(timerave%pddold,timer%pddold,timern%pddold)
      call get_t_ave(timerave%step,timer%step,timern%step)
    
    end if
    
    ! If iterations equals timerstep, get the timer averages, output
    ! to file, and reset timers
    if(n_step.eq.n_step/dto_timer*dto_timer .or. n_step .eq. 1) then
      boundary = timerave%boundary / max(timern%boundary,1.d0)
      climate  = timerave%climate  / max(timern%climate,1.d0)
      emb      = timerave%emb      / max(timern%emb,1.d0)
      pddnew   = timerave%pddnew   / max(timern%pddnew,1.d0)
      pddold   = timerave%pddold   / max(timern%pddold,1.d0)
      step     = timerave%step     / max(timern%step,1.d0)

      write(329,"(f15.2,7f15.5)")     n_step/day_year, tot*convert,  &
                                      step, boundary, climate, emb, &
                                      pddnew, pddold
      
      ! Reset the averages and counts to zero
      call timer_reset(timerave)
      call timer_reset(timern)
      
    end if
    
    ! Reset the step timer
    call timer_reset(timer)

    
    return  
  
  end subroutine timing
  
  subroutine get_t_ave(tave,tstep,n)
   
    implicit none
    
    real (8) :: tstep, tave, n
    real (8), parameter :: convert = 1.d0/60.d0     ! [sec] => [min]
    
    if (tstep .ne. 0.d0) then
      n = n + 1.d0
      tave = tave+tstep*convert
    end if
    
    return
    
  end subroutine get_t_ave
  
  subroutine timer_reset(tm)
    
    use emb_global
    
    implicit none
    
    type(timers) :: tm
    
    tm%boundary = 0.d0
    tm%climate  = 0.d0
    tm%emb      = 0.d0
    tm%pddnew   = 0.d0
    tm%pddold   = 0.d0
    tm%step     = 0.d0
    
    return
  
  end subroutine timer_reset



