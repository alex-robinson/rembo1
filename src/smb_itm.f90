module smb_itm
  ! This module contains all subroutines related to the calculation of
  ! surface mass balance. The external interface to this module is via
  ! the melt_budget() subroutine, which takes arrays from the main
  ! program as input and outputs arrays of the calculated variables.
  
  use emb_global
  !use emb_functions
  ! Note: depends on emb_global for global parameter
  ! values and for input/output

  implicit none
  
contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : m e l t _ b u d g e t
  ! Author     : Alex Robinson
  ! Purpose    : Determine the total melt, accumulation and 
  !              surface mass balance at a given point
  !              * input in mm water equivalent
  ! * As in Ettema et al (2009) supplementary information,
  !   surface mass balance is defined as:
  !     SMB  = snow + rain - runoff   [kg m2 / d ] == [mm / d]
  !   runoff = rain + melt - refrozen
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine melt_budget( m2, zs, S, tt, tte, pdds, pp, snow, mmfac, as, ap, &
                          h_snow, rf, runoff, runoff_snow,  &
                          runoff_rain, melt, melted_ice, new_ice, refrozen, dh, &
                          melt_insol, melt_S0, S0, annual )

    implicit none
    
    integer, parameter :: ny = nys, nx = nxs

    double precision, dimension(ny,nx) :: m2, zs, tt, tte, pdds, mmfac
    double precision, dimension(ny,nx) :: rf, as, ap, S
    
    double precision, dimension(ny,nx) :: pp, snow, rain, h_snow, h_snow0
    double precision, dimension(ny,nx) :: refreezing_max, refrozen_snow_max
    double precision, dimension(ny,nx) :: smb, smbi, dh
    double precision, dimension(ny,nx) :: mm_snow
    
    double precision, dimension(ny,nx) ::              refrozen_rain, runoff_rain
    double precision, dimension(ny,nx) :: melted_snow, refrozen_snow, runoff_snow
    double precision, dimension(ny,nx) :: melted_firn, refrozen_firn, runoff_firn
    double precision, dimension(ny,nx) :: melted_ice, snow_to_ice, new_ice
    double precision, dimension(ny,nx) :: melt, refrozen, runoff 
    double precision, dimension(ny,nx) :: melt_insol, melt_S0, S0
    double precision, dimension(ny,nx) :: S_melt, S_diagnose

    double precision, dimension(ny,nx) :: mm_snow_2, melted_snow_2, melted_ice_2

    logical :: is_annual
    logical, optional :: annual
    
    is_annual = .FALSE.
    if ( present(annual)) is_annual = annual
    
    ! Decide which insolation values to use for actual melt calculations
    S_melt = S
    S_diagnose = S0
    if (itm_S0 .eq. 2 .or. itm_S0 .eq. 3) then
      S_melt = S0
      S_diagnose = S 
    end if

    ! Determine initial values   
    rain = pp - snow

    ! Determine preliminary surface and planetary albedo
    call get_albedo(as, ap, m2, h_snow, pdds)
    
    ! Add additional snowfall
    h_snow0 = h_snow 
    h_snow = h_snow + snow
    
    ! Make sure all snow disappears over the ocean
    ! (important for paleo runs with varying sea-level)
    !where ( m2 .eq. 2.d0 ) h_snow = 0.d0
    
    if (melt_choice .eq. 1) then
    
      ! Get amount of melt from oerlemans scheme
      call itm(mm_snow,zs,tt,as,ap,S_melt)

      ! Add additional melt from correction term
      if (itm_S0 .eq. 3) then 
        mm_snow = mm_snow + melt_insol
        where(mm_snow .lt. 0.d0) mm_snow = 0.d0
      end if

    else

      ! Melt of snow and ice is based on pdd factors for snow and ice
      !mm_snow = tte * mm_teff_snow
      
      ! Melt of snow and ice is based on pdd factors for snow and ice
      ! and additional melt due to insolation anomaly
      mm_snow = tte * mm_teff_snow + melt_insol
      where(mm_snow .lt. 0.d0) mm_snow = 0.d0


    end if
    
    ! Determine how much snow and ice would be melted today
    where (mm_snow .gt. h_snow)
      
      ! All snow is melted, the rest of energy converted to melt some ice
      ! The rest of energy will go into melting ice
      melted_snow = h_snow
      melted_ice  = (mm_snow - h_snow) * mmfac

    elsewhere
      
      ! Snow melt will use all energy, none left for ice melt
      melted_snow = mm_snow
      melted_ice  = 0.d0
      
    end where    
    
    ! Now calculate the melt (total ablation)
    melt   = melted_snow + melted_ice
    
    ! If desired, diagnose melt if insolation were equal to present day (or vice versa)
    if (melt_choice .eq. 1 .and. itm_S0 .gt. 0) then

      call itm(mm_snow_2,zs,tt,as,ap,S_diagnose)

      where (mm_snow .gt. h_snow)
        ! All snow is melted, the rest of energy converted to melt some ice
        ! The rest of energy will go into melting ice
        melted_snow_2 = h_snow
        melted_ice_2  = (mm_snow_2 - h_snow) * mmfac
      elsewhere
        ! Snow melt will use all energy, none left for ice melt
        melted_snow_2 = mm_snow_2
        melted_ice_2  = 0.d0
      end where    
      
      ! Total melt for S0
      melt_S0 = melted_snow_2 + melted_ice_2 

    end if 

    ! Remove melted snow, if any, from the snow height budget
    h_snow = h_snow - melted_snow
    
    ! Adjust the albedo (accounting for actual amount of melt)
    call get_albedo(as, ap, m2, h_snow, pdds, melt)
    
    ! Determine how much new ice is made from compression of remaining snow
    snow_to_ice = 0.d0
    where (h_snow .gt. h_snow_max)
      ! Assume excess contributes to new ice
      snow_to_ice = (h_snow - h_snow_max)
      
      ! Reset snow height down to maximum height
      h_snow = h_snow_max
    end where
    
    ! If using annual PDD approach...
    if ( is_annual ) then
      snow_to_ice = snow - melted_snow
      h_snow = snow
    end if
    
    ! Determine what fraction of the melted snow and rain will refreeze, 
    ! (Note: rf is zero if not on ice sheet or there is no snow cover)
    rf = 0.d0
    where ( m2 .eq. 0.d0 .and. h_snow .gt. 0.d0 )
    
      rf = Pmaxfrac * snow / max(1d-5, (snow + rain))    ! max() here ensures no division by zero, if (snow+rain)==0, then rf is zero anyway
      
      ! Modify refreezing factor based on height of snow
      where ( h_snow .gt. 2e3 )
        rf = 1.d0                                         ! refreezing factor is 1 for large snow heights.      
      else where ( h_snow .gt. 1e3 )
        rf = rf + ((h_snow-1e3)/(2e3-1e3) ) * (1.d0 - rf) ! linear function increasing to rf=1 as h_snow increases to 2m.
      end where
      
    end where
   
    ! Determine the actual maximum amount of refreezing
    refreezing_max    = h_snow                                ! Total refreezing depends on amount of snow left!
    refrozen_rain     = min(rain*rf,refreezing_max)           ! First rain takes up refreezing capacity
    refrozen_snow_max = refreezing_max - refrozen_rain        ! Subtract rain from refreezing capacity to determine what's left for snow
    refrozen_snow     = min(melted_snow*rf,refrozen_snow_max) ! melted_snow uses remaining capacity if it needs it
    refrozen          = refrozen_snow + refrozen_rain         ! Amount of ice created from refreezing

    ! Determine how much water will runoff for each component
    runoff_snow = melted_snow - refrozen_snow                 ! Net snow melt
    runoff_rain = rain - refrozen_rain                        ! Net rain
    runoff      = runoff_snow + runoff_rain + melted_ice      ! Total runoff
    
    ! Get the total ice accumulated for the day
    new_ice = snow_to_ice + refrozen
    
    ! Get the global surface mass balance
    smb = snow + rain - runoff
    
    ! Get the internal surface mass balance (what ice sheet model needs)
    ! These should apply at separate places, so
    ! smbi = -melted_ice, for negative mass balance
    ! smbi =     new_ice, for positive mass balance
    smbi = new_ice - melted_ice
    
    ! Determine the actual melted amount of snow or snow+ice,
    ! depending on where melt is taking place (for energy-balance)
    ! Note: no refreezing off of ice sheet, so there total melt can
    ! only equal the amount of snow melted
    ! Note: this value can be negative too (negative contribution to energy flux)
    where ( m2 .eq. 0 ) 
      dh = melt-refrozen       ! Technically refrozen=0 here, but include it anyway
    elsewhere
      dh = melted_snow-refrozen
    end where
    
    ! Output variables for mass balance
    !smb, smbi, snow, rain, runoff, runoff_rain, new_ice (accum), melt, melted_ice, refrozen
    
    return
  
  end subroutine melt_budget
    
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : g e t _ a l b e d o
  ! Author     : Alex Robinson
  ! Purpose    : Determine the current surface and planetary albedos
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine get_albedo(as, ap, mask, hsnow, PDDs, melt_in)
    
    implicit none
    
    integer, parameter :: ny = nys, nx = nxs
    
    real (8), dimension(:,:) :: as, ap, mask, hsnow, PDDs
    real (8), dimension(:,:), optional :: melt_in
    real (8), dimension(ny,nx) :: as_snow, as_ground, depth, as_snow_planet
    real (8), dimension(ny,nx) :: hsnow_critical, melt
    integer :: n, k
    
    ! Determine the critical snow height based on number of PDDs,
    ! which correspond to vegetation growth
    ! 0-100 PDDs    = desert, 10mm
    ! 100-1000 PDDs = tundra (grass), linearly increasing between 10mm => 100mm
    ! 1000 PDDs +   = forest, 100mm
    where ( PDDs .le. 1d2 ) 
      hsnow_critical = hsnow_crit0 ! 10 mm
    else where ( PDDs .ge. 1d3 )
      hsnow_critical = hsnow_crit1 ! 100 mm
    elsewhere
      hsnow_critical = hsnow_crit0 + (hsnow_crit1-hsnow_crit0) * (( PDDs - 1d2 ) / ( 1d3 - 1d2))
    end where
    
    ! Determine the scaled snowdepth
    depth = hsnow / hsnow_critical
    
    ! Determine the amount of melt to affect abledo calculation.
    ! Default is high melt everywhere, so that jump in albedo is avoided
    ! at onset of melt season.
    ! After calculating actual melt, this is included as an argument and
    ! albedo is recalculated.
    melt = melt_crit+1.d0
    if ( present(melt_in) ) melt = melt_in
    
    ! Figure out what the surface albedo would be with no snow, 
    ! based on the type of ground underneath ( ice or land )
    as_ground = as_ice
    where (mask .eq. 2.d0) as_ground = as_ocean
    where (mask .eq. 1.d0) as_ground = as_land*(1d3-min(PDDs,1d3))/(1d3-0d0) &
                               + as_land_forest*min(PDDs,1d3)/(1d3-0d0)

    ! Determine current snow albedo: if melt gt eg, 1mm/day, then change albedo to melting snow!
    as_snow = as_snow0
    where (melt .gt. melt_crit) as_snow = as_snow1   
    
    ! Determine snow albedo for use with REMBO (to get planetary albedo)
    ! Where PDDs > 1d3, forest snow
    as_snow_planet = as_snow
    where ( hsnow_critical .ge. hsnow_crit1 ) as_snow_planet = as_snow_forest  

    ! First, calculate surface albedo to use for REMBO, 
    ! which is not necessarily identical to snow albedo
    ! minimum albedo is that of bare ground
    as = min(as_ground + depth*(as_snow_planet-as_ground), as_snow_planet)
    
    ! Get the current planetary albedo
    ap = ap0_intercept + ap0_slope * as
    
    ! Get current surface albedo to be used for ITM melt scheme
    ! It will either be the maximum albedo (that of dry snow: as_snow)
    ! or the wet snow albedo plus a fraction depending on height of snow
    ! minimum albedo now should be that of wet snow minimum
    as = min(as_snow2 + depth*(as_snow-as_snow2), as_snow)
    
    return
    
  end subroutine get_albedo
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : i t m (insolation temperature melt)
  ! Author     : Alex Robinson
  ! Purpose    : Determine the total melt, accumulation and 
  !              surface mass balance at a given point
  !              * output in mm water equivalent
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine itm(melt,zs,tt,as,ap,S)
    
    implicit none
    
    integer, parameter :: ny = nys, nx = nxs
    
    real (8), dimension(:,:) :: melt, zs, tt, as, ap, S
    real (8), dimension(ny,nx) :: at
    
    integer :: i, j, k, q
    real (8) :: dt

    
    ! Find out how many seconds in given time period (always daily!)
    dt = sec_year / nk
    
    at = at_intercept + at_slope* max(zs,0.d0)
    
    ! Calculate potential melt
    melt = dt*(at*(1.d0 - as)*S + itm_cc + itm_t*tt) / (rho_w*Lm)
    
    melt = max( melt, 0.d0 ) * 1d3     ! [m/day] => [mm/day], only positive melt

    return
    
  end subroutine itm

end module smb_itm 

subroutine get_pdd_corr_melt(corr,melt)
  
  use emb_global 
  
  implicit none 

  double precision :: cmax 
  double precision, dimension(nys,nxs) :: corr, melt

  cmax = 1.d0-1.d0/(60.d0**pdd_spread)

  corr = 0.d0 
  where (melt .gt. 0.d0) corr = 1.d0-1.d0/((melt+1.d0)**pdd_spread)

  corr = corr / cmax*pdd_amax 
  where (corr .gt. pdd_amax) corr = pdd_amax

  return

end subroutine get_pdd_corr_melt 

subroutine get_pdd_corr(corr,tt,k)

  use emb_global 

  implicit none

  integer, parameter :: kmin = 1, kmid = 180, kmax = 360
  integer :: k
  double precision :: Tmax_now, Tmin_now, wt
  double precision, dimension(nys,nxs) :: corr, tt

  ! Specify the current Tmax value (from parameter choice)
  Tmax_now = pdd_Tmax
  Tmin_now = sin_nonlin(k,pdd_Tmax,pdd_Tmin,pdd_spread)

  ! Calculate the correction value as a function of temperature
  ! Limit the values to be between 0 and specified max value
  corr = (tt - Tmin_now) / (Tmax_now - Tmin_now) * pdd_amax
  where (corr .lt. 0.d0)     corr = 0.d0 
  where (corr .gt. pdd_amax) corr = pdd_amax
  
  return

contains 

  function sin_nonlin(x,Tmax,Tmin,spread)

    implicit none 

    integer :: x
    double precision :: Tmax, Tmin, spread 
    double precision :: amp, offset, ymax
    double precision :: y 
    double precision :: sin_nonlin

    ! Get sin parameters based on Tmax and Tmin
    amp    = (Tmax-Tmin)/2.d0
    offset = Tmin + amp

    ! Calculate the maximum y value (at day 180) for normalization
    ! and calculate today's value
    ymax = (sin((pi/180)*(180+90)-pi)+1.d0)**spread
    y = (sin((pi/180)*(x+90)-pi)+1.d0)**spread

    ! Normalize and rescale today's value
    y = (-y/ymax*2.d0+1.d0)*amp + offset

    sin_nonlin = y 
    
    return

  end function sin_nonlin


end subroutine get_pdd_corr


subroutine get_pdd_corr_old(corr,tt)
  use emb_global 

  implicit none
  double precision, dimension(:,:) :: corr, tt

  corr = pdd_a+pdd_b*tt
  where(corr .lt. 0.d0)     corr = 0.d0 
  where(corr .gt. pdd_amax) corr = pdd_amax

  return

end subroutine get_pdd_corr_old
