program emb_driver
  
  use emb_global
  
  implicit none

  integer :: n_step
  double precision :: timer_start, timer_tot

  write(*,*)
  write(*,*) "  *** 2D Diffusion: Temperature and precipitation ***"
  write(*,*) "                   over Greenland" 
  write(*,*)
      
  ! Get start time in seconds
  call cpu_time(timer_start) 
  
  ! Initialize the model
  call sclimate(0)
  call timing(0,timer_start,timer_tot)
  
  nstep_max = yearf - year0
  
  ! ### Run iterations ###
  do n_step = 1, nstep_max    ! in years

    call sclimate(n_step)

    ! Update the timers for each timestep and output
    call timing(n_step,timer_start,timer_tot)
  
  end do 

  
end program emb_driver
      

