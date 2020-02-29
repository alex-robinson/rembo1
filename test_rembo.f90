program emb_driver
  
    use rembo_sclimate 
    
    implicit none

    integer :: n_step, n_step_max 
    double precision :: timer_start, timer_tot
    real(8) :: year0
    real(8) :: yearf 
    real(8) :: T_summer 

    write(*,*)
    write(*,*) "  *** 2D Diffusion: Temperature and precipitation ***"
    write(*,*) "                   over Greenland" 
    write(*,*)
    
    T_summer = 0.0 

    year0 = 0.0 
    yearf = 10.0 

    ! Get start time in seconds
    call cpu_time(timer_start) 
  
    ! Initialize the climate model REMBO, including loading parameters from options_rembo 
    call rembo_init()
    call timing(0,timer_start,timer_tot)

    ! Update REMBO, with ice sheet topography    
    call rembo_update(0,T_summer)
     
    n_step_max = yearf - year0
  
    ! ### Run iterations ###
    do n_step = 1, n_step_max    ! in years

        ! call REMBO1     
        call rembo_update(n_step,T_summer)
            
        ! Update the timers for each timestep and output
        call timing(n_step,timer_start,timer_tot)

    end do 

  
end program emb_driver
      

