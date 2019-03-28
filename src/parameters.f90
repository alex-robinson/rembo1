
!! ### FORTRAN90 library
!! ### PARAMETER LOADING FUNCTIONS

module parameters

  implicit none
  
  !! ### PARAMETER variables
  type param_type
    character (len=30) :: name
    real (8) :: val
  end type
  
  type cparam_type
    character (len=30) :: name
    character (len=256) :: val
  end type
  
  type aparam_type
    character (len=30) :: name
    real (8), allocatable :: val(:)
  end type
  
  ! All necessary program params
  type(param_type),  target, allocatable :: params(:), glob_p(:)
  type(cparam_type), target, allocatable :: cparams(:), glob_cp(:)
  
contains
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  e m b _ a r g s
  ! Author     :  Alex Robinson
  ! Purpose    :  Get command line arguments from program call
  !               *only obtains character arguments
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  subroutine args(folder)
    
    implicit none
    
    integer :: i, narg, check
    
    character (len=256)  :: string
    character (len=*)    :: folder

    ! Get the number of arguments total
    narg = command_argument_count()
    
    ! Set default values (in case no arguments provided)
    !folder = "output/test/"
    folder = "./"
    
    if (narg .gt. 0) then
      
      ! Get the last argument, the rest are ignored
      ! note: this should be arg #1 (after the executable)
      call get_command_argument(narg,string)
      folder = trim(adjustl(string))
    end if
    
    i = len(trim(folder))
    if ( scan(trim(folder),"/",back=.TRUE.) .ne. i ) folder = trim(folder)//"/"
    
    write(*,*) "len(folder): ",len(trim(folder)),scan(trim(folder),"/",back=.TRUE.)
    write(*,"(a1,5x,a,a)") "e","folder: ", trim(folder)
    
    return
    
  end subroutine args
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  g e t _ p a r a m s
  ! Author     :  Alex Robinson
  ! Purpose    :  Get all parameter values from specified options file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine get_params(fnm,io_in,glob)
    
    implicit none
    
    integer :: i, k, stat
    integer, parameter :: nmax = 300
    double precision :: param, tmp, val(nmax)
    character (len=*) :: fnm
    character (len=30) :: tmpc, nm(nmax), cnm(nmax)
    character (len=256) :: tmpc1, cval(nmax)
    character (len=1) :: com, split
    integer, optional :: io_in, glob
    integer :: io, global
    
    io = 6
    if (present(io_in)) io = io_in
    
    global = 0
    if (present(glob)) global = glob
    
    write(io,*) "Parameter file: ",trim(fnm)
    
    ! ################################################## 
    ! First obtain all character parameters
    ! ################################################## 
    
    open(unit=29,file=trim(fnm),status="old")
    
    k = 0
    stat = 0
    ! Loop through file, only take parameters from lines
    ! with an "=" as the 31st character
    do i = 1, nmax*2
      
      read(29,"(a1,9x,a30,a1,a)",iostat=stat) com, tmpc, split, tmpc1
      if (stat .eq. -4) exit
      
      if (com .ne. "#" .and. split .eq. ":") then
        k = k + 1
        nm(k)   = trim(adjustl(tmpc))
        cval(k) = trim(adjustl(tmpc1))
      end if
    
    end do
    
    close(29)
    
    write(*,*) k, "character parameters found."
    
    ! Allocate the global variable params to correct length
    if (allocated(cparams)) deallocate(cparams)
    allocate(cparams(k))
    
    ! Fill in cparams with the values obtained from file
    do i = 1, k
      cparams(i)%name = nm(i)
      cparams(i)%val  = cval(i)
    end do
    
    ! ################################################## 
    ! Now do the same to obtain all numerical parameters
    ! ################################################## 
    
    open(unit=29,file=trim(fnm),status="old")

    k = 0
    stat = 0
    ! Loop through file, only take parameters from lines
    ! with an "=" as the 31st character
    do i = 1, nmax*2
      
      read(29,"(a1,9x,a30,a1,f15.0)",iostat=stat) com, tmpc, split, tmp
      if (stat .eq. -1) exit
      
      if (com .ne. "#" .and. split .eq. "=") then
        k = k + 1
        nm(k)  = trim(adjustl(tmpc))
        val(k) = tmp      
      end if
    
    end do
    
    close(29)
    
    write(*,*) k, "numerical parameters found."
    
    ! Allocate the global variable params to correct length
    if (allocated(params)) deallocate(params)
    allocate(params(k))
    
    ! Fill in params with the values obtained from file
    do i = 1, k
      params(i)%name = nm(i)
      params(i)%val  = val(i)
    end do
    
    ! Save as global values if desired
    if ( global .eq. 1 ) then
      
      if ( allocated(glob_p) ) then
        write(*,*) "glob_p parameter list already allocated."
        write(*,*) "Something is happening out of order."
        stop
      end if
    
      ! Allocate the global parameter arrays to store the global variables
      !  (that way they can replace any module variables with the same names)
      k = size(params)
      allocate(glob_p(k))
      glob_p = params
      
      k = size(cparams)
      allocate(glob_cp(k))
      glob_cp = cparams
      
      write(io,*) "(saved as global parameters)"
      
    end if
    
    return
    
  end subroutine get_params
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  o u t _ p a r a m s
  ! Author     :  Alex Robinson
  ! Purpose    :  Write parameter values to a file or output number
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine out_params(fnm)
    
    implicit none
    
    integer :: i, io
    character (len=*), optional :: fnm
    
    io = 6
    if (present(fnm)) then
      io = 70
      open(unit=io,file=trim(fnm),status="unknown")
    end if
    
    write(*,*) ""
    write(io,"(a1,2x,i5,2x,a)") "#",size(cparams), "character parameters."
    do i = 1, size(cparams)
      write(io,"(a30,a1,1x,a)") cparams(i)%name,"=",'"'//trim(cparams(i)%val)//'"'
    end do
    write(*,*) ""
    write(io,"(a1,2x,i5,2x,a)") "#",size(params), "numerical parameters."
    do i = 1, size(params)
      write(io,"(a30,a1,1x,g15.4)") params(i)%name,"=",params(i)%val
    end do
    
    if (present(fnm)) close(io)
    
    return
    
  end subroutine out_params
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  s u b _ p a r a m s
  ! Author     :  Alex Robinson
  ! Purpose    :  Substitute current parameter values with the global
  !               ones (if any are duplicates)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sub_params()
    
    implicit none
    
    integer :: i, j
    
    do i = 1, size(cparams)
      
      do j = 1, size(glob_cp)
        
        if ( trim(adjustl(cparams(i)%name)) .eq. trim(adjustl(glob_cp(j)%name)) ) then
          
          cparams(i)%val = glob_cp(j)%val
          write(*,"(2x,a30,a)") trim(adjustl(cparams(i)%name)), " using global value."
          exit
          
        end if
      
      end do
    
    end do
    
    do i = 1, size(params)
      
      do j = 1, size(glob_p)
        
        if ( trim(adjustl(params(i)%name)) .eq. trim(adjustl(glob_p(j)%name)) ) then
          
          params(i)%val = glob_p(j)%val
          write(*,"(2x,a30,a)") trim(adjustl(params(i)%name)), " using global value."
          exit
          
        end if
      
      end do
    
    end do
    
    return
    
  end subroutine sub_params
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  p a r a m
  ! Author     :  Alex Robinson
  ! Purpose    :  Output a parameter value based on the string name
  !               of the parameter
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  function param(nm,val,set)
    
    implicit none
    
    integer :: i
    double precision :: param
    double precision, optional :: val
    character (len=*) :: nm
    
    type(param_type), target, optional :: set(:)
    type(param_type), pointer :: parameters(:)
    
    parameters => params
    if (present(set)) parameters => set
    
    param = -9999.d0   ! Error value, in case param name not found
    
    ! Cycle through all parameters and find the index of the current
    ! one specified by nm
    do i = 1, size(parameters)
      if (trim(nm) .eq. trim(adjustl(parameters(i)%name))) then
        ! If val is an argument, set the value of the parameters
        if (present(val)) parameters(i)%val = val   
        ! Set the output value equal to the value of the parameter
        param = parameters(i)%val
        exit
      end if
    end do

    if (param .eq. -9999.d0) then
      write(*,*) "parameter error: param doesn't exist, ",trim(nm)
      stop
    end if
    
    return
    
  end function param
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  p a r a m
  ! Author     :  Alex Robinson
  ! Purpose    :  Output a parameter value based on the string name
  !               of the parameter
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  function cparam(nm,val,set)
    
    implicit none
    
    integer :: i
    character (len=256) :: cparam
    character (len=*), optional :: val
    character (len=*) :: nm
    
    type(cparam_type), target, optional :: set(:)
    type(cparam_type), pointer :: parameters(:)
    
    parameters => cparams
    if (present(set)) parameters => set
    
    cparam = "NA"  ! Error value, in case param name not found
    
    ! Cycle through all parameters and find the index of the current
    ! one specified by nm
    do i = 1, size(parameters)
      if (trim(nm) .eq. trim(adjustl(parameters(i)%name))) then
        ! If val is an argument, set the value of the parameters
        if (present(val)) parameters(i)%val = val   
        ! Set the output value equal to the value of the parameter
        cparam = trim(adjustl(parameters(i)%val))
        exit
      end if
    end do

    if (trim(cparam) .eq. "NA") then
      write(*,*) "parameter error: param doesn't exist, ",trim(nm)
      stop
    end if
    
    return
    
  end function cparam



  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  g e t _ o t i m e
  ! Author     :  Alex Robinson
  ! Purpose    :  Load a vector of "output times" from a file
  !               in order to standardize output across modules
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  subroutine get_output_times(nm,otimes,ntot,y0,yf,dt)
    
    implicit none
    
    character (len=*) :: nm
    character (len=256), parameter :: fnm = "out.times.txt"
    character (len=256) :: tmpc1, tmpc2
    integer :: i, stat, n, ntot
    double precision :: y0, yf, t0, tf, dt1
    double precision, optional :: dt
    double precision :: tmp(500)
    double precision, pointer :: otimes(:)
    !double precision, pointer, allocatable :: get_otime(:)
    
    ! Check for presence of dt, adjust dt1 as necessary
    dt1 = -1.d0
    if ( present(dt) ) dt1 = dt*1d-3
    
    ! Adjust years to ka for input times
    t0 = y0*1d-3; tf = yf*1d-3
    
    ! If dt1 is greater than zero, generate evenly spaced
    ! output times
    if ( dt1 .gt. 0.d0 ) then
      
      write(*,*) trim(nm)//": generating evenly-spaced output times..."
      
      !! Determine at which times 2d output will be generated
      !! (useful to know for netcdf file)
      n = int( (tf - t0) / dt1 )+1
      allocate(otimes(n))
      
      otimes(1) = t0; do i = 2, n; otimes(i) = otimes(i-1) + dt1; end do   ! ka
!       do i = 1, n
!         otimes(i) = t0 + (i-1)*dt1  ! ka
!       end do
      
    else if (dt1 .lt. 0.d0 ) then !! If dt < 0, then load values from file
      
      write(*,*) trim(nm)//": loading output times from file "//trim(fnm)
      
      ! Open the file containing the output times
      open(unit=30,file=trim(fnm),status="old")
      
      stat = 0
      ! Loop through file, only take parameters from lines
      ! with an "=" as the 31st character
      do i = 1, 200
        
        ! Read the lines that of eg, "nm = 10"
        read(30,*,iostat=stat) tmpc1, tmpc2, n
        
        ! If end of file is reached, exit
        if (stat .eq. -4) exit
        
        ! If the line of the desired array is found, exit this loop
        if ( trim(tmpc1) .eq. trim(nm) ) exit
        
      end do
      
      if (stat .eq. -4) then
        write(*,*) "parameters::get_otime:: End of file reached in: "//trim(fnm)
        write(*,*) "   Output array not found: "//trim(nm)
        write(*,*) "   Stopping."
        stop
      else
        
        ! Allocate the output array
        if (associated(otimes)) deallocate(otimes)
        allocate(otimes(n))
        
        ! Next line in the file should contain the n parameter values
        read(30,*) otimes
      
      end if
      
      close(30)
      
      ! Figure out which indices are in range of simulation
      ! Remove unused values
      do i = 1, n; if ( otimes(i) .ge. t0 ) exit; end do
      otimes = cshift(otimes,(i-1)); n = n - (i-1)
      i = 1
      do while ( otimes(i) .lt. tf .and. i .le. n )
        i = i+1
      end do
      n = i-1
      
      ! Add one index for final year if n==0 are found
      if (n .eq. 0) then
        n = 1
        otimes(1) = tf 
      end if

      write(*,*) "dt =",dt1
      write(*,*) "n  =",n
      
    else
      write(*,*) "parameters::get_output_times:: incorrect dt! dt =", dt1
      stop
    end if  !! end loading section (from file) !!
    
    ! Store otimes temporarily to large array to add new values
    tmp(1:n) = otimes(1:n)
    
    ! Make sure last year of simulation will be part of array
    if ( tmp(n) .lt. tf ) then
      n = n+1; tmp(n) = tf
    end if
    
    ! Send the final array in the function output
    if (associated(otimes)) deallocate(otimes)
    allocate(otimes(n)); otimes = tmp(1:n)
    
    ntot = n
    
    ! Write the output times to standard out
    write(*,"(a,i4,a,1000f10.3)") trim(nm)//" =",ntot,":",otimes
    
    ! Check to make sure output is not excessive!
    if ( ntot .gt. 500 ) then
      write(*,*) "parameters::get_output_times:: Too much output desired!! Modify program choices and try again."
      stop
    end if
    
    return
    
  end subroutine get_output_times
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  match
  ! Author     :  Alex Robinson
  ! Purpose    :  find the index of a double within a vector of doubles
  !               Return the index, or 0 if not found
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  function match(t0,t)
    
    integer :: i, match
    double precision :: t0, t(:)    

    ! Determine the current output index
    match = 0
    do i=1, size(t)
      if ( dabs(t0-t(i)) .lt. 1d-3 ) match = i
    end do
    
  end function match



end module parameters
