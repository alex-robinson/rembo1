program remboyelmo_driver

    use emb_global
    use yelmo 

    implicit none

    integer :: n_step
    double precision :: timer_start, timer_tot

    ! Yelmo variables
    character(len=512) :: path_out, path_par, path_const 
    character(len=512) :: path_file1D, path_file2D 

    type(yelmo_class)  :: yelmo1 
    
    write(*,*)
    write(*,*) "                         ===== Greenland simulation ====="
    write(*,*) "Atmosphere (temperature and precipitation): 2D energy-moisture balance model (REMBO)"
    write(*,*) "Snowpack: One-layer snowpack model with ITM or PDD melt scheme (inside REMBO)"
    write(*,*) "Ice sheet: 3D thermomechanical hybrid ice sheet model (Yelmo)" 
    write(*,*)

    ! Get start time in seconds
    call cpu_time(timer_start) 

    ! Initialize the climate model REMBO
    call sclimate(0)
    call timing(0,timer_start,timer_tot)


    ! Initialize the ice sheet model Yelmo 
    path_out    = "."
    path_par    = trim(path_out)//"/yelmo_Greenland.nml"
    path_const  = trim(path_out)//"/yelmo_const_Earth.nml"
    path_file1D = trim(path_out)//"yelmo1D.nc"
    path_file2D = trim(path_out)//"yelmo2D.nc"
  
    ! General initialization of yelmo constants (used globally)
    call yelmo_global_init(path_const)

    ! Initialize data objects
    call yelmo_init(yelmo1,filename=path_par)
    


    nstep_max = yearf - year0
  
    ! ### Run iterations ###
    do n_step = 1, nstep_max    ! in years

        call sclimate(n_step)

        ! Update the timers for each timestep and output
        call timing(n_step,timer_start,timer_tot)
  
    end do 

  
end program remboyelmo_driver
      

