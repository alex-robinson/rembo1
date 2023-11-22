program emb_driver
    
    use nml 
    use rembo_sclimate 
    
    implicit none

    integer :: n_step, n_step_max 
    double precision :: timer_start, timer_tot

    real(8) :: year0
    real(8) :: yearf 
    real(8) :: T_summer 
    real(8) :: T_mon(12)
    real(8) :: co2 

    write(*,*)
    write(*,*) "  *** 2D Diffusion: Temperature and precipitation ***"
    write(*,*) "                   over Greenland" 
    write(*,*)
    
    ! Get start time in seconds
    call cpu_time(timer_start) 
  
    ! Initialize the climate model REMBO, including loading parameters from options_rembo 
    call rembo_init(year0)
    call timing(0,timer_start,timer_tot)
    
    ! Read in parameters
    call nml_read("rembo_Greenland.nml","ctrl","year0",      year0)
    call nml_read("rembo_Greenland.nml","ctrl","yearf",      yearf)
    call nml_read("rembo_Greenland.nml","ctrl","T_summer",   T_summer)
    
    !year0    = 0.0 
    !yearf    = 10.0 
    !T_summer = 5.0 
    T_mon    = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,5.0,4.0,3.0,2.0,1.0]
    !T_summer = 0.0
    !T_mon    = 0.0
    co2      = 350.0 

    ! Update REMBO, with ice sheet topography    
    !call rembo_update(year0,T_summer)
    call rembo_update(year0,T_summer,dT_mon=T_mon,co2=co2)
    
    n_step_max = yearf - year0
  
    ! ### Run iterations ###
    do n_step = 1, n_step_max    ! in years

        ! call REMBO1
        call rembo_update(year0+n_step,T_summer) 
        !call rembo_update(year0+n_step,T_summer,dT_mon=T_mon,co2=co2)
            
        ! Update the timers for each timestep and output
        call timing(n_step,timer_start,timer_tot)

    end do 

    ! Write a restart file
    call rembo_write_restart("./rembo_restart.nc",year0)

    
end program emb_driver
      

