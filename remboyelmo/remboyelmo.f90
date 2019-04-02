program remboyelmo_driver

    !use emb_global
    use rembo_sclimate 
    use yelmo 

    ! External libraries
    use sealevel 
    use isostasy  
    
    use marine_shelf   
    use basal_hydrology 
    use sediments 
    use geothermal
    
    use hyster 
    use ncio 

    implicit none

    ! Control variables 
    real(prec) :: time_init, time_end, time_equil, time
    real(prec) :: dtt, dt1D_out, dt2D_out    
    logical    :: calc_transient_climate

    logical, parameter :: use_hyster = .FALSE. 

    ! REMBO variables
    integer :: n_step, nstep_max 
    double precision :: timer_start, timer_tot

    ! Yelmo variables
    character(len=512) :: path_out, path_par, path_const 
    character(len=512) :: path_file1D, path_file2D 
    real(prec) :: conv_km3_Gt
    real(prec) :: enh_out, cf_out 

    type(yelmo_class)      :: yelmo1 
    type(sealevel_class)   :: sealev
    type(hyster_class)     :: hyst1   
    type(isos_class)       :: isos1
    type(marshelf_class)   :: mshlf1
    type(hydro_class)      :: hyd1 
    type(sediments_class)  :: sed1 
    type(geothermal_class) :: gthrm1
    
    write(*,*)
    write(*,*) "                         ===== Greenland simulation ====="
    write(*,*) "Atmosphere (temperature and precipitation): 2D energy-moisture balance model (REMBO)"
    write(*,*) "Snowpack: One-layer snowpack model with ITM or PDD melt scheme (inside REMBO)"
    write(*,*) "Ice sheet: 3D thermomechanical hybrid ice sheet model (Yelmo)" 
    write(*,*)

    ! Initialize the ice sheet model Yelmo 
    path_out    = "./"
    path_par    = trim(path_out)//"yelmo_Greenland.nml"
    path_const  = trim(path_out)//"yelmo_const_Earth.nml"
    path_file1D = trim(path_out)//"yelmo1D.nc"
    path_file2D = trim(path_out)//"yelmo2D.nc"
  
    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)
    conv_km3_Gt = rho_ice *1e-3

    ! Initialize data objects
    call yelmo_init(yelmo1,filename=path_par)
    
    ! Define parameter values for printing 
    enh_out = yelmo1%mat%par%enh_stream 
    cf_out  = yelmo1%dyn%par%cf_stream 
    if (yelmo1%dyn%par%mix_method .eq. -2) cf_out = 0.0 

    ! Load control parameters (timing, etc)
    call nml_read(path_par,"control","time_init",    time_init)                 ! [yr] Starting time
    call nml_read(path_par,"control","time_end",     time_end)                  ! [yr] Ending time
    call nml_read(path_par,"control","time_equil",   time_equil)                ! [yr] Years to equilibrate first
    call nml_read(path_par,"control","dtt",          dtt)                       ! [yr] Main loop time step 
    call nml_read(path_par,"control","dt1D_out",     dt1D_out)                  ! [yr] Frequency of 1D output 
    call nml_read(path_par,"control","dt2D_out",     dt2D_out)                  ! [yr] Frequency of 2D output 
    call nml_read(path_par,"control","transient",    calc_transient_climate)    ! Calculate transient climate? 

    ! Get start time in seconds
    call cpu_time(timer_start) 

    ! Initialize the climate model REMBO, including loading parameters from options_rembo 
    call rembo_init()
    call timing(0,timer_start,timer_tot)

    ! Initialize state variables (topo only)
    call yelmo_init_state_1(yelmo1,path_par,time=time_init)
    
    ! Initialize all Yelmo helping modules 
    call initialize_yelmo_externals(sealev,hyst1,isos1,mshlf1,hyd1,sed1,gthrm1,yelmo1,path_par)

    ! === Update initial boundary conditions for current time and yelmo state =====
    ! ybound: z_bed, z_sl, H_sed, H_w, smb, T_srf, bmb_shlf , Q_geo

    call isos_init_state(isos1,z_bed=yelmo1%bnd%z_bed,z_bed_ref=yelmo1%bnd%z_bed, &
                               H_ice_ref=yelmo1%tpo%now%H_ice,z_sl=yelmo1%bnd%z_sl,time=time_init)

    call sealevel_update(sealev,year_bp=time_init)

    yelmo1%bnd%z_sl  = sealev%z_sl 
    yelmo1%bnd%H_sed = sed1%now%H 
    yelmo1%bnd%H_w   = 0.0   ! use hydro_init_state later

    ! Update REMBO, with ice sheet topography    
    call rembo_update(0,real(yelmo1%tpo%now%z_srf,8),real(yelmo1%tpo%now%H_ice,8))
            
    ! Update surface mass balance and surface temperature from REMBO
    yelmo1%bnd%smb   = rembo_ann%smb    *conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
    yelmo1%bnd%T_srf = rembo_ann%T_srf
    
    ! Set ocean temperature and temp anomaly
    call marshelf_set_Tshlf(mshlf1,to_ann=273.15,dto_ann=0.0)

    call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%z_sl,regions=yelmo1%bnd%regions,dx=real(yelmo1%grid%G%dx,prec))

    yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
    yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

    yelmo1%bnd%Q_geo    = gthrm1%now%ghf 
    
    call yelmo_print_bound(yelmo1%bnd)

    time = time_init 

    ! Initialize state variables (dyn,therm,mat)
    ! (initialize temps with robin method with a cold base)
    call yelmo_init_state_2(yelmo1,path_par,time=time_init,thrm_method="robin-cold")

    ! Run yelmo for several years with constant boundary conditions and topo
    ! to equilibrate thermodynamics and dynamics
!     call yelmo_update_equil(yelmo1,time,time_tot=50.0,topo_fixed=.FALSE.,dt=1.0,ssa_vel_max=500.0)
    
    ! 2D file 
    call yelmo_write_init(yelmo1,path_file2D,time_init=time,units="years")
    call write_step_2D_combined(yelmo1,isos1,mshlf1,hyd1,path_file2D,time=time)
    
    ! 1D file 
    call write_yreg_init(yelmo1,path_file1D,time_init=time_init,units="years",mask=yelmo1%bnd%ice_allowed)
    call write_yreg_step(yelmo1%reg,path_file1D,time=time) 
    
    ! Determine total iterations [yr]
    nstep_max = time_end - time_init 
  
    ! ### Run iterations ###
    do n_step = 1, nstep_max    ! in years

        ! Get current driver time [yr]
        time = time_init + n_step 

        ! == SEA LEVEL ==========================================================
        call sealevel_update(sealev,year_bp=time)
        yelmo1%bnd%z_sl  = sealev%z_sl 

        ! == Yelmo ice sheet ===================================================
        call yelmo_update(yelmo1,time)

        ! == ISOSTASY ==========================================================
        call isos_update(isos1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_sl,time)
        yelmo1%bnd%z_bed = isos1%now%z_bed


if (calc_transient_climate) then
        ! == CLIMATE (ATMOSPHERE AND OCEAN) ====================================
        if (mod(time,10.0)==0) then
            
            ! call REMBO1     
            call rembo_update(n_step,real(yelmo1%tpo%now%z_srf,8),real(yelmo1%tpo%now%H_ice,8))
            
            ! Update surface mass balance and surface temperature from REMBO
            yelmo1%bnd%smb   = rembo_ann%smb    *conv_we_ie*1e-3       ! [mm we/a] => [m ie/a]
            yelmo1%bnd%T_srf = rembo_ann%T_srf
                
        end if 

!         ! == MARINE AND TOTAL BASAL MASS BALANCE ===============================
        
        call marshelf_update(mshlf1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%f_grnd, &
                         yelmo1%bnd%z_sl,regions=yelmo1%bnd%regions,dx=real(yelmo1%grid%G%dx,prec))

end if 
        
        yelmo1%bnd%bmb_shlf = mshlf1%now%bmb_shlf  
        yelmo1%bnd%T_shlf   = mshlf1%now%T_shlf  

        ! == BASAL HYDROLOGY ====================================================
        call hydro_update(hyd1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_bed,yelmo1%tpo%now%z_srf,yelmo1%bnd%z_sl, &
                          yelmo1%tpo%now%bmb,yelmo1%tpo%now%f_grnd,yelmo1%tpo%par%dx,yelmo1%tpo%par%dx,time)
        yelmo1%bnd%H_w   = hyd1%now%H_water 

        if (mod(time,dt2D_out)==0) then 
            call write_step_2D_combined(yelmo1,isos1,mshlf1,hyd1,path_file2D,time=time)
        end if 

        if (mod(time,dt1D_out)==0) then 
            call write_yreg_step(yelmo1%reg,path_file1D,time=time) 
        end if 

        !if (mod(time,100.0)==0) then
        !    write(*,"(a,f14.4)") "yelmo::       time = ", time
        !end if 

        if (mod(time,1000.0)==0) then
            ! Write hysteresis summary 
            write(*,"(a5,2f8.1,f15.3,2f5.1)") "tble", time, yelmo1%reg%V_ice_g*1e-6, yelmo1%reg%dVicedt*conv_km3_Gt, enh_out, cf_out 
            call flush() 
        end if 

        ! Update the timers for each timestep and output
        call timing(n_step,timer_start,timer_tot)
  
    end do 


contains 

    subroutine initialize_yelmo_externals(sealev,hyst,isos,mshlf,hyd,sed,gthrm,ylmo,path_par)

        implicit none 

        type(sealevel_class)   :: sealev
        type(hyster_class)     :: hyst   
        type(isos_class)       :: isos
        type(marshelf_class)   :: mshlf
        type(hydro_class)      :: hyd 
        type(sediments_class)  :: sed 
        type(geothermal_class) :: gthrm
        
        type(yelmo_class)      :: ylmo
        character(len=*)       :: path_par  

        ! Local variables 
        type(grid_class)       :: grd 
        character(len=256)     :: domain 

        ! === Initialize external models (forcing for ice sheet) ======

        ! Store yelmo grid in grd and domain name as a shortcut 
        grd  = ylmo%grid 
        domain = ylmo%par%domain 

        ! Initialize global sea level model (bnd%z_sl)
        call sealevel_init(sealev,path_par)

        ! Initialize bedrock model 
        call isos_init(isos,path_par,grd%G%nx,grd%G%ny,real(grd%G%dx*grd%xy_conv))

        ! Initialize hysteresis module for transient forcing experiments 
        call hyster_init(hyst,path_par,time_init) 
        
        ! Initialize marine melt model (bnd%bmb_shlf)
        call marshelf_init(mshlf,path_par,grd%G%nx,grd%G%ny,domain,grd%name,ylmo%bnd%basins)
        
        ! Initialize basal hydrology (bnd%H_w)
        call hydro_init(hyd,path_par,grd%G%nx,grd%G%ny)

        ! Load other constant boundary variables (bnd%H_sed, bnd%Q_geo)
        call sediments_init(sed,path_par,grd%G%nx,grd%G%ny,domain,grd%name)
        call geothermal_init(gthrm,path_par,grd%G%nx,grd%G%ny,domain,grd%name)

        return 

    end subroutine initialize_yelmo_externals


    subroutine write_step_2D_combined(ylmo,isos,mshlf,hyd,filename,time)

        implicit none 
        
        type(yelmo_class),      intent(IN) :: ylmo
        type(isos_class),       intent(IN) :: isos 
        !type(snapclim_class),   intent(IN) :: snp 
        type(marshelf_class),   intent(IN) :: mshlf  
        type(hydro_class),      intent(IN) :: hyd  
        !type(sediments_class),  intent(IN) :: sed 
        !type(geothermal_class), intent(IN) :: gthrm
        !type(isos_class),       intent(IN) :: isos
        
        character(len=*),  intent(IN) :: filename
        real(prec), intent(IN) :: time

        ! Local variables
        integer    :: ncid, n
        real(prec) :: time_prev 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write model speed 
        call nc_write(filename,"speed",ylmo%par%model_speed,units="kyr/hr",long_name="Model speed (Yelmo only)", &
                      dim1="time",start=[n],count=[1],ncid=ncid)
        
        ! == yelmo_topography ==
        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"mb_applied",ylmo%tpo%now%mb_applied,units="m",long_name="Applied net mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"H_grnd",ylmo%tpo%now%H_grnd,units="m",long_name="Ice thickness overburden", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"N_eff",ylmo%tpo%now%N_eff,units="bar",long_name="Effective pressure", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"calv",ylmo%tpo%now%calv,units="m/a ice equiv.",long_name="Calving rate", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded fraction", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Ice fraction in grid cell", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="km", &
!                       long_name="Distance to nearest grounding-line point", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"C_bed",ylmo%dyn%now%C_bed,units="m a^-1 Pa^-2",long_name="Dragging constant", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m-1",long_name="Dragging coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff",ylmo%dyn%now%visc_eff,units="Pa a^-1",long_name="Vertically averaged viscosity", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"uxy_i_bar",ylmo%dyn%now%uxy_i_bar,units="m/a",long_name="Internal shear velocity magnitude", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal sliding velocity magnitude", &
                     dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically integrated velocity magnitude", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
!                       dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"Q_strn",ylmo%thrm%now%Q_strn/(rho_ice*ylmo%thrm%now%cp),units="K a-1",long_name="Strain heating", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="J a-1 m-2",long_name="Basal frictional heating", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"bmb_grnd",ylmo%thrm%now%bmb_grnd,units="m/a ice equiv.",long_name="Basal mass balance (grounded)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"ATT",ylmo%mat%now%ATT,units="a^-1 Pa^-3",long_name="Rate factor", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)

        call nc_write(filename,"enh_bar",ylmo%mat%now%enh_bar,units="1",long_name="Vertically averaged enhancement factor", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"visc_bar",ylmo%mat%now%visc_bar,units="Pa a^-1",long_name="Vertically averaged viscosity", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! Boundaries
        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid) 
        call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_w",ylmo%bnd%H_w,units="m",long_name="Basal water layer thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb",ylmo%bnd%smb,units="m/a ice equiv.",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb_errpd",ylmo%bnd%smb-ylmo%dta%pd%smb_ann,units="m/a ice equiv.",long_name="Surface mass balance error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"T_srf",ylmo%bnd%T_srf,units="K",long_name="Surface temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a ice equiv.",long_name="Basal mass balance (shelf)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_sl",ylmo%bnd%z_sl,units="m",long_name="Sea level rel. to present", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! External data
!         call nc_write(filename,"dzbdt",isos%now%dzbdt,units="m/a",long_name="Bedrock uplift rate", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
 
!         call nc_write(filename,"dT_shlf",mshlf%now%dT_shlf,units="K",long_name="Shelf temperature anomaly", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"Ta_ann",snp%now%ta_ann,units="K",long_name="Near-surface air temperature (ann)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"Ta_sum",snp%now%ta_sum,units="K",long_name="Near-surface air temperature (sum)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"Pr_ann",snp%now%pr_ann*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"dTa_ann",snp%now%ta_ann-snp%clim0%ta_ann,units="K",long_name="Near-surface air temperature anomaly (ann)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"dTa_sum",snp%now%ta_sum-snp%clim0%ta_sum,units="K",long_name="Near-surface air temperature anomaly (sum)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"dPr_ann",(snp%now%pr_ann-snp%clim0%pr_ann)*1e-3,units="m/a water equiv.",long_name="Precipitation anomaly (ann)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Climate and surface 

        call nc_write(filename,"Ta_ann",rembo_ann%T_ann,units="K",long_name="Near-surface air temperature (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Ta_jja",rembo_ann%T_jja,units="K",long_name="Near-surface air temperature (jja)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Ta_djf",rembo_ann%T_djf,units="K",long_name="Near-surface air temperature (djf)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Pr_ann",rembo_ann%pr*1e-3,units="m/a water equiv.",long_name="Precipitation (ann)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"smb_ann",rembo_ann%smb*1e-3,units="m/a water equiv.",long_name="Surface mass balance (ann)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        

        ! Ice thickness comparison with present-day 
        call nc_write(filename,"H_ice_errpd",ylmo%tpo%now%H_ice-ylmo%dta%pd%H_ice,units="m",long_name="Ice thickness error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_s_errpd",ylmo%dyn%now%uxy_s-ylmo%dta%pd%uxy_s,units="m",long_name="Ice thickness error wrt present day", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"ssa_mask_acx",ylmo%dyn%now%ssa_mask_acx,units="1",long_name="SSA mask (acx)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"ssa_mask_acy",ylmo%dyn%now%ssa_mask_acy,units="1",long_name="SSA mask (acy)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"dzsdx",ylmo%tpo%now%dzsdx,units="m/m",long_name="Surface slope", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"dzsdy",ylmo%tpo%now%dzsdy,units="m/m",long_name="Surface slope", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"f_grnd_acx",ylmo%tpo%now%f_grnd_acx,units="1",long_name="Grounded fraction (acx)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"f_grnd_acy",ylmo%tpo%now%f_grnd_acy,units="1",long_name="Grounded fraction (acy)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"taub_acx",ylmo%dyn%now%taub_acx,units="Pa",long_name="Basal stress (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"taub_acy",ylmo%dyn%now%taub_acy,units="Pa",long_name="Basal stress (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"taud_acx",ylmo%dyn%now%taud_acx,units="Pa",long_name="Driving stress (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"taud_acy",ylmo%dyn%now%taud_acy,units="Pa",long_name="Driving stress (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
!         call nc_write(filename,"ux_b",ylmo%dyn%now%ux_b,units="m/a",long_name="Basal sliding velocity (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_b",ylmo%dyn%now%uy_b,units="m/a",long_name="Basal sliding velocity (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"ux_i_bar",ylmo%dyn%now%ux_i_bar,units="m/a",long_name="Internal shear velocity (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_i_bar",ylmo%dyn%now%uy_i_bar,units="m/a",long_name="Internal shear velocity (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"ux_bar",ylmo%dyn%now%ux_bar,units="m/a",long_name="Vertically averaged velocity (x)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"uy_bar",ylmo%dyn%now%uy_bar,units="m/a",long_name="Vertically averaged velocity (y)", &
!                        dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

!         call nc_write(filename,"beta_acx",ylmo%dyn%now%beta_acx,units="Pa a m-1",long_name="Dragging coefficient (x)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
!         call nc_write(filename,"beta_acy",ylmo%dyn%now%beta_acy,units="Pa a m-1",long_name="Dragging coefficient (y)", &
!                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine write_step_2D_combined

end program remboyelmo_driver
      

