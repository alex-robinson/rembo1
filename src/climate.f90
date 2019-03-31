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
  subroutine sclimate(n_step,z_srf,H_ice)
  
    use emb_global
    use emb_functions
    use rembo_main
  
    implicit none
    
    integer, parameter :: nx=nxs, ny=nys
    real (8), parameter :: dx=dxs
    
    integer, intent(IN) :: n_step
    real(8), intent(IN), optional :: z_srf(nx,ny)    ! Input from ice-sheet model
    real(8), intent(IN), optional :: H_ice(nx,ny)    ! Input from ice-sheet model 
    
    real (8), dimension(ny,nx) :: m2, zs, lats, lons, aco2
    real (8), dimension(ny,nx) :: ZZ, tma, tmj, ampl
    real (8), dimension(ny,nx) :: restart_tt, restart_pp
    real (8), dimension(ny,nx) :: restart_snowh
    real (8), dimension(ny,nx) :: restart_dh, restart_ap
    real (8), dimension(ny,nx) :: mask
    real (8) :: wt, pscalar, sfac

    character(len=10) :: c_dxs, dattype
    character(len=256) :: fnm, restart_folder
    character(len=15) :: c_yr
    
    integer :: i, j, k, qq, n, n_now
    real (8) :: dtime1, dtime2
    real (8) :: yearnow, yearnow1, get_year
    real (8) :: T_warming_now, tmp_noise
    
    real (8) :: bndtmp
    
    type(choices) :: now

    ! Kill program as needed (for fracs and other simulations)
    if ( kill .eq. 1 ) then
      write(*,*) "Kill switch activated, program finished."
      stop
    end if
  
    call cpu_time(timer%climate)           ! get current time in seconds   
    
    ! Set emb global nstep to driver's n_step 
    nstep = n_step
    
    ! Initialize information output location
    if (nstep .eq. 0) then
      
      call emb_globinit(1)
      
      if(init_rembo .eq. 0) then
        
        ! Allocate day, mon and year structures
        !allocate( day(nk), mon(nm) )
        
        ! Initialize stuff needed for boundary routine called in ini_sico!
        m2        = fields0%m2
        zs        = fields0%zs
        lats      = fields0%lats
        lons      = fields0%lons
        
        ! Calculate the horizontal gradient of the initial topography
        call hgrad(fields0%zs,dx,fields0%dzs)
        
        write(*,"(a1,5x,i10,5x,a)") "e",nstep,"Initialized topography for sclimate."
        
        ! Initialize output nc file
        call rembo_nc(trim(outfldr)//"clima.nc","years",time_out(1)*1d3, &
                      lats*todegs,lons*todegs,sico_grid%x,sico_grid%y,mask_hydro=fields0%mask_hydro)  
                          
        init_rembo = 1
      end if

      if ( init_paleo .eq. 0 .and. boundary_forcing .gt. 0) then

        call load_forcing(trim(forcing_file),year0+year_offset,yearf+year_offset)
        write(*,"(a1,5x,i10,5x,a)") "e",nstep,"Initialized boundary forcing from file: "//trim(forcing_file)
        init_paleo = 1
        
      end if
      
    end if
    
    ! ##### Adjust current distribution for temporal shift #####
    ! yearnow needed for boundary forcing and sinsol2d
    if (transient .ne. 0) then           ! Transient runs
      
      yearnow = get_year(nstep)
      
      ! ## ADJUSTMENT FOR PALEO RUNS (after year 00) ##
      if ( transient .eq. 2 .and. yearnow .gt. 0.d0 ) yearnow = 0.d0  

    else                                 ! Equilibrium runs
      
      yearnow = year0 + year_offset
    
    end if

    ! Adjust time step for calling SMB module
    ! HACK for long transient future simulations (ajr: 15.05.2012)
    !if ( yearnow .gt. 1e3 ) dtime_smb = 10
    
    ! Determine which functions should be called now
    ! Check which operations should be performed for this timestep
    now%clim = .FALSE.
    now%smb  = .FALSE.
    if( nstep.eq.nstep/dtime_emb*dtime_emb ) now%clim = .TRUE.
    if( nstep.eq.nstep/dtime_smb*dtime_smb ) now%smb  = .TRUE.

    ! If clim or smb is running, update forcing and topo
    if ( now%clim .or. now%smb ) then
      
      ! Update fields from exchange (if necessary)
      call vars2clim(zs,m2,transT%dVdt)
      
      ! Get current paleo forcing 
      !write(*,*) "climate:: forcing:",size(forcing)
      if ( boundary_forcing .gt. 0 ) call get_forcing_now(yearnow,m2)
      
      ! Determine amount of atmospheric co2 based on warming
      ! T_warming = T_global_mean_warming = T_greenland_summer_warming
      T_warming_now = deltaT(0)
      aco2 = co2(T_warming_now)     
      
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
          call rembo(yearnow,now,m2,zs,lats,lons,aco2)
          write(*,*) "rembo equili, ",qq
        end do
        
        write(*,*) "Finished with rembo equilibration."

        equili = 0
        
        ! Reset boundary variables
        bnd_ave   = bndtmp
        
      end if
                
      ! Call the energy-moisture balance module (output sent to variable "saved")  
      call rembo(yearnow,now,m2,zs,lats,lons,aco2)      

    else   ! climchoice=0 => bilinear interp + data 
      
      ! Call conventional approach (output sent to variable "saved") 
      if ( now%smb ) call conventional(m2,zs,lats,lons,aco2)   
      
    end if
    
    ! #### OUTPUT SECTION ####
        
    ! Get the current year, as relevant to output
    yearnow1 = get_year(n_step)

    ! Determine the current output index
    n_now = match(yearnow1,time_out*1d3)
    
    if ( n_now .ne. 0 ) then
      
      !! Output cntro file
      call rembo_nc(trim(outfldr)//"clima.cntro.nc","years",yearnow1, &
                lats*todegs,lons*todegs,sico_grid%x,sico_grid%y,mask_hydro=fields0%mask_hydro)
      call rembo_nc_write(trim(outfldr)//"clima.cntro.nc",ann,1,yearnow1,mask=m2,zs=zs, &
                      tjja=saved%tjja, tdjf=saved%tdjf,ttp=saved%ttp, &
!                       pkappa=saved%pkappa,qqfac=fields0%qqfac, &
                      pdds=saved%pdds, write_mon = .TRUE. )
  
      !! Output annual file
!       call rembo_nc_write(trim(outfldr)//"clima.nc",ann,n_now,yearnow1,mask=m2,zs=zs, &
!                           tjja=saved%tjja, tdjf=saved%tdjf,ttp=saved%ttp, &
! !                           pkappa=saved%pkappa,qqfac=fields0%qqfac, &
!                           pdds=saved%pdds, write_mon = .TRUE. )
      call rembo_nc_write_small(trim(outfldr)//"clima.nc",ann,mon,n_now,yearnow1,mask=m2,zs=zs, &
                                tjja=saved%tjja,tdjf=saved%tdjf,ttp=saved%ttp,pdds=saved%pdds)

      !! Output monthly file
      if ( write_rembo_m .eq. 1 ) then
        
        call rembo_nc(trim(outfldr)//"climm.cntro.nc","months",time_out_mon(1),lats*todegs,lons*todegs, &
                      sico_grid%x,sico_grid%y,mask=m2,zs=zs)
          
        do k = 1, nm
          call rembo_nc_write(trim(outfldr)//"climm.cntro.nc",mon(k),k,time_out_mon(k))
        end do
      end if   
      
      !! Output daily file
      if ( write_rembo_d .eq. 1 ) then      ! Write daily output
        ! The same bug applies as above in the monthly file!!
        call rembo_nc(trim(outfldr)//"climd.cntro.nc","days",time_out_day(1),lats*todegs,lons*todegs, &
                      sico_grid%x,sico_grid%y,mask=m2,zs=zs)
          
        do k = 1, nk
          call rembo_nc_write(trim(outfldr)//"climd.cntro.nc",day(k),k,time_out_day(k))
        end do
      end if

    end if
      
    ! Determine the current output index for restart file
    n_now = match(yearnow1,time_out_r*1d3)

    if ( write_rembo_r .eq. 1 .and. nstep .ne. 0 .and. n_now .ne. 0 ) then
      call rembo_restart(yearnow1,m2,zs,day(nk),saved%pdds)  ! Write restart file  
    end if
      
    ! Now for fractional restart file
    if ( write_rembo_r .eq. 2 ) then
    
      ! (Initialize tmp file, so rembo won't write restart file,
      !  but has something to read)
      if ( nstep .eq. 0 ) then
        open(99,file=trim(outfldr)//"tmp123456789",status="unknown")
        write(99,*) "na"
        close(99)
      end if
    
      ! Get the filename from the tmp file written in sicopolis
      open(99,file=trim(outfldr)//"tmp123456789",status="old")
      read(99,"(a)") fnm; fnm = adjustl(fnm)
      close(99)
      
      ! Check the filename, if it's good,
      ! write to the restart file
      if ( trim(fnm) .ne. "na" ) then
        
        write(*,*) "Writing fractional restart file: "//trim(fnm)
        
        call rembo_restart(yearnow1,m2,zs,day(nk),saved%pdds,file=trim(fnm))  ! Write restart file
        
        ! Wipe out filename in tmp file
        open(99,file=trim(outfldr)//"tmp123456789",status="unknown")
        write(99,*) "na"
        close(99)
        
      end if

    end if
      
    ! Output various 1D files, by region
    if (nstep .eq. nstep/dto_clim*dto_clim ) then
       
      ! First by sector (five sectors total, hard coded for now)
      if ( clim_coupled .lt. 0 .and. .FALSE. ) then
        do qq = 1, 5
          
          mask = 0.d0
          where( m2 .eq. 0.d0 .and. fields0%mask_hydro .eq. dble(qq) ) mask = 1.d0
          
          write(fnm,"(a11,i1,a3)") "rembo.gis.S",qq,".nc"
          call climchecker_new(trim(outfldr)//trim(fnm), &
                              day,mon,ann,mask,zs,lats,lons,nstep,yearnow1)
          if ( nstep .eq. 0 ) init_summary2 = 0
        
        end do
      end if 
      
      ! Now by total gis area
      mask = 0.d0
      where( m2 .eq. 0.d0 ) mask = 1.d0
      call climchecker_new(trim(outfldr)//"rembo.gis.nc", &
                          day,mon,ann,mask,zs,lats,lons,nstep,yearnow1)
      if ( nstep .eq. 0 ) init_summary2 = 0
      
!!!   ! Kill condition for transient simulations (enough ice points and enough time passed)
!!!      if ( sum(mask) .lt. 150.d0 .and. yearnow1 .gt. 20d3 ) kill = 1
      
      ! And total greenland area
      mask = 0.d0
      where( m2 .le. 1.d0 ) mask = 1.d0
      call climchecker_new(trim(outfldr)//"rembo.grl.nc", &
                          day,mon,ann,mask,zs,lats,lons,nstep,yearnow1)

      ! Write/calc comparison to observational fields
      if ( trim(domain) .eq. "GRL" ) call climchecker2(m2,zs)
                                
    end if
    
    ! #### END OUTPUT SECTION ####
        
    if (now%clim) then
      write(*,"(a1,5x,a10,f12.2)") "e","yearnow = ",yearnow
      write(*,"(a1,5x,a10,f12.2)") "e","co2 = ",co2(T_warming_now)
    end if
    
    call cpu_time(timer%now)                   ! get current time in seconds
    timer%climate = (timer%now-timer%climate)  ! Get elapsed time in seconds
    
    if ( clim_coupled .eq. 1 .and. (now%smb .or. now%clim) ) then
      call vars2ice(saved%tt,saved%tdjf,saved%tjja,saved%ttp,saved%tts, &
                    saved%tjan, saved%tjul,saved%precip,saved%snow, &
                    saved%runoff_snow,saved%runoff_rain,saved%melted_ice, &
                    saved%refrozen,saved%evap,saved%smb,saved%h_snow, "save")
    end if
    
    ! Output variables for climber
    !if ( boundary_forcing .eq. 3 .and. (now%smb .or. now%clim) ) then
    !  call write_for_climber(trim(outfldr),nstep,lats*todegs,m2,zs,transT%dVdt)
    !end if 

    return
  
  end subroutine sclimate
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : c o n v e n t i o n a l
  ! Author     : Alex Robinson, 27. Oct 2009
  ! Purpose    : Determine the climate and melt based on conventional 
  !              approaches
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine conventional(mask,zs,lats,lons,aco2)
    
    use emb_global
    use emb_functions
    use rembo_functions
    use rembo_main

    implicit none
    
    integer, parameter :: nx = nxs, ny = nys

    double precision, dimension(ny,nx) :: mask, zs, lats, lons, aco2
    double precision, dimension(ny,nx) :: tmp
    integer :: i, j, k, m, km, km2
    
    character(len=10)  :: c_dxs, dattype
    character(len=256) :: fnm
    
    call cpu_time(timer%pddold)           ! get current time in seconds
    
    ! Get the precipitation and temperature fields from data
    ! (they will be stored in the global day/mon/ann structures)
    call precip0(mask,zs,lats,lons)
    call temper0(mask,zs,lats,lons)

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
    
    ! ajr: this appears to be broken... shouldn't it add onto the saved%pdds each month?
    ! check in detail!
    do m = 1, nm
      saved%pdds = max( mon(m)%tt, 0.d0 ) * day_month
    end do
    
    saved%precip      = ann%pp          *1d-3  / sec_year /rho_i_w     ! mm.w.e. / a => m.i.e. / s
    saved%snow        = ann%snow        *1d-3  / sec_year /rho_i_w     ! mm.w.e. / a => m.i.e. / s                                                  
    saved%melted_ice  = ann%melted_ice  *1d-3  / sec_year /rho_i_w     ! mm.w.e. / a => m.i.e. / s
    saved%runoff_snow = ann%runoff_snow *1d-3  / sec_year /rho_i_w     ! mm.w.e. / a => m.i.e. / s
    saved%runoff_rain = ann%runoff_rain *1d-3  / sec_year /rho_i_w     ! mm.w.e. / a => m.i.e. / s
    saved%smb         = (ann%pp-ann%runoff) *1d-3  / sec_year /rho_i_w     ! mm.w.e. / a => m.i.e. / s
    saved%evap        = 0.d0
                                                    ! *Note: divide by rho_i_w, bc m3 in denominator!
        
    call cpu_time(timer%now)           ! get current time in seconds
    timer%pddold = (timer%now-timer%pddold)  ! Get elapsed time in seconds
    
    write(*,"(a1,5x,i10,5x,a)") "e",nstep, "Called pdd (annual)."
    
    return
  
  end subroutine conventional
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : c l i m _ t i m e
  ! Author     : Alex Robinson, 27. Oct 2008
  ! Purpose    : Determine the current year given the timestep in years
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function get_year(n_step)
    
    use emb_global
    
    implicit none
    
    integer :: n_step
    double precision :: get_year
    
    get_year = year0 + n_step
  
    return
  
  end function get_year


