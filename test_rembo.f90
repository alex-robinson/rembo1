program emb_driver
    
    use nml 
    use ncio
    use rembo_sclimate 
    
    implicit none

    integer :: n_step, n_step_max 
    double precision :: timer_start, timer_tot

    real(8) :: time_init
    real(8) :: time_end 
    real(8) :: time_ins
    real(8) :: T_summer 
    real(8) :: T_mon(12)
    real(8) :: co2 

    integer :: nx, ny 
    real(8), allocatable :: z_srf(:,:)
    real(8), allocatable :: H_ice(:,:)
    real(8), allocatable :: z_sl(:,:)
    real(8), allocatable :: reg_mask(:,:)

    ! Get start time in seconds
    call cpu_time(timer_start) 
    
    write(*,*)
    write(*,*) "  *** 2D Diffusion: Temperature and precipitation ***"
    write(*,*) "                   over Greenland" 
    write(*,*)
    
    ! Read in parameters
    call nml_read("rembo_Greenland.nml","ctrl","time_init",  time_init)
    call nml_read("rembo_Greenland.nml","ctrl","time_end",   time_end)
    call nml_read("rembo_Greenland.nml","ctrl","time_ins",   time_ins)
    call nml_read("rembo_Greenland.nml","ctrl","T_summer",   T_summer)
    
    ! Initialize the climate model REMBO, including loading parameters from options_rembo 
    call rembo_init(time_init)
    call timing(0,timer_start,timer_tot)
    
    
    !time_init = 0.0 
    !time_end  = 10.0 
    !T_summer  = 5.0 
    !T_mon      = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,5.0,4.0,3.0,2.0,1.0]
    !T_summer  = 0.0
    T_mon     = 0.0
    co2       = 350.0 

    ! Define a topography to work with
    nx = size(rembo_ann%smb,1)
    ny = size(rembo_ann%smb,2)

    allocate(z_srf(nx,ny))
    allocate(H_ice(nx,ny))
    allocate(z_sl(nx,ny))
    allocate(reg_mask(nx,ny))

    ! Load the topography 
    call load_topo_rtopo2(z_srf,H_ice,reg_mask,path="ice_data/Greenland/GRL-16KM", &
                                            domain="Greenland",grid_name="GRL-16KM")
    z_sl = 0.0 

    ! Testing modifications to topography and sea level
    !z_srf = z_srf*0.5
    !H_ice = H_ice*0.5
    !z_sl = 10.0 

    ! Update REMBO, with ice sheet topography    
    !call rembo_update(time_init,T_summer)
    call rembo_update(time_init,time_ins,T_summer,z_srf,H_ice,z_sl,dT_mon=T_mon,co2=co2)
    
    n_step_max = time_end - time_init
  
    ! ### Run iterations ###
    do n_step = 1, n_step_max    ! in years

        ! call REMBO1
        call rembo_update(time_init+n_step,time_ins,T_summer,z_srf,H_ice,z_sl,dT_mon=T_mon,co2=co2) 
              
        ! Update the timers for each timestep and output
        call timing(n_step,timer_start,timer_tot)

    end do 


    ! Write a restart file
    call rembo_restart_write("./rembo_restart.nc",time_init)


contains

    subroutine load_topo_rtopo2(z_srf,H_ice,reg_mask,path,domain,grid_name)
        ! Load the data into the rembo_class object

        implicit none 

        real(8),           intent(INOUT) :: z_srf(:,:) 
        real(8),           intent(INOUT) :: H_ice(:,:) 
        real(8),           intent(INOUT) :: reg_mask(:,:) 
        character(len=*),  intent(IN)    :: path 
        character(len=*),  intent(IN)    :: domain
        character(len=*),  intent(IN)    :: grid_name 

        ! Local variables
        character(len=512) :: filename
        real(8) :: region_number 

        filename = trim(path)//"/"//trim(grid_name)//"_TOPO-RTOPO-2.0.1.nc"

        ! Static fields

        ! ## Surface elevation ##
        call nc_read(filename,"z_srf",z_srf)
        
        ! ## Ice thickness ##
        call nc_read(filename,"H_ice",H_ice)
        
        ! Load regions to delete regions out of interest 
        filename = trim(path)//"/"//trim(grid_name)//"_REGIONS.nc"
        call nc_read(filename,"mask",reg_mask)
        
        if (trim(domain) .eq. "Greenland") then 
            region_number = 1.3 
        else if (trim(domain) .eq. "Antarctica") then  
            region_number = 2.11
        else 
            write(*,*) "Domain not recognized: "//trim(domain)
            stop 
        end if 
         
        where (reg_mask .ne. region_number) 
            reg_mask = 0.0 
        elsewhere 
            reg_mask = 1.0 
        end where
        
        return 

    end subroutine load_topo_rtopo2

end program emb_driver
      

