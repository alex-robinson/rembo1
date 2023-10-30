module emb_functions
  
  use emb_global
  use ncio
  use interp1D

  implicit none
  
 contains
  
  subroutine rembo_update_boundary_forcing_monthly(T_anomaly,T_warming,T_monthly,day_month,day_year)
  
    implicit none

    real(8), intent(OUT) :: T_anomaly(:)    ! Daily temperature anomaly values
    real(8), intent(OUT) :: T_warming       ! Summer average temperature anomaly
    
    real(8), intent(IN)  :: T_monthly(12)   ! Monthly temperature anomalies to interpolate
    real(8), intent(IN) :: day_month 
    real(8), intent(IN) :: day_year 

    ! Calculate temperature anomaly for each day of the year
    call calc_daily_values(T_anomaly,T_monthly,day_month,day_year)

    ! Store summer mean value too
    T_warming = sum(T_monthly(6:8)) / 3d0

    return

  end subroutine rembo_update_boundary_forcing_monthly

  subroutine rembo_update_boundary_forcing_sin(T_anomaly,T_warming,T_summer,T_winter,day_year)
  
    implicit none

    real(8), intent(OUT) :: T_anomaly(:)    ! Daily temperature anomaly values
    real(8), intent(OUT) :: T_warming       ! Summer average temperature anomaly
    
    real(8), intent(IN)  :: T_summer
    real(8), intent(IN)  :: T_winter
    real(8), intent(IN) :: day_year 

    ! Local variables

    ! Calculate temperature anomaly for each day of the year
    call calc_deltaT_sin(T_anomaly,T_summer,T_winter,day_year)

    ! Store summer mean value too
    T_warming = T_summer 

    return

  end subroutine rembo_update_boundary_forcing_sin

  subroutine calc_deltaT_sin(dT,T_summer,T_winter,day_year)
    ! Typically T_winter = T_summer * T_wintfac, with T_wintfac=2 (winter 2x warmer than summer)
    
    real(8), intent(OUT) :: dT(:)
    real(8), intent(IN)  :: T_summer
    real(8), intent(IN)  :: T_winter
    real(8), intent(IN)  :: day_year    ! total days per year (day_year=360 typically)
    ! Local variables
    real (8) :: Tamp, Tmean
    integer  :: k 

    Tamp    = (T_winter - T_summer) / 2.d0
    Tmean   = (T_winter + T_summer) / 2.d0
    
    ! Get the temperature anomaly for each day of the year
    do k = 1, int(day_year)
      dT(k) = Tmean + Tamp * dcos(2.d0*pi*(k-15)/day_year)
    end do 

    return

  end subroutine calc_deltaT_sin

  subroutine calc_daily_values(var_daily,var_mon,day_month,day_year)

    implicit none

    real(8), intent(OUT) :: var_daily(:)    ! size=nk
    real(8), intent(IN)  :: var_mon(:)      ! size=nm
    real(8), intent(IN)  :: day_month 
    real(8), intent(IN)  :: day_year 

    ! Local variables
    integer :: dnow, q, k, nm
    integer :: q0, q1, q2 
    integer :: mid, day 
    real(8) :: wt0, wt1, wt2, wttot

    nm = size(var_mon,1)

    ! ####################################################################
    ! Interpolate data in time: monthly => daily
    ! ####################################################################
    dnow = 0

    do q = 1, nm

      do k = 1, int(day_month)

        dnow = dnow+1

        q1 = q                     ! q1 is current month

        q0 = q1-1                  ! q0 is previous month
        if (q0 .eq. 0) q0 = 12     ! Loop to december

        if (q1 .eq. 13) q1 = 1     ! Loop to january for current month

        q2 = q1+1                  ! q2 is next month
        if (q2 .eq. 13) q2 = 1     ! Loop to january for next month

        mid = day_month / 2   ! Halfway point of current month (+/- 1 day)
        day = k - mid

        ! When day < mid, weight shared between wt0 and wt1
        ! When day = mid, weight goes to current month
        ! When day > mid, weight shared between wt1 and wt2
        wt0 = dble (- min(0,day))
        wt1 = day_month - abs(day)
        wt2 = dble (  max(0,day))

        ! Normalize weights
        wttot = dble(wt0 + wt1 + wt2)
        wt0 = wt0 / wttot
        wt1 = wt1 / wttot
        wt2 = wt2 / wttot
        
        var_daily(dnow) = var_mon(q0)*wt0 + &
                          var_mon(q1)*wt1 + &
                          var_mon(q2)*wt2
     
      end do
    end do
    
    return

  end subroutine calc_daily_values


! ==========


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  l o a d _ f o r c i n g
  ! Author     :  Alex Robinson
  ! Purpose    :  Load relevant forcing data:
  !               1. Monthly temperature anomalies
  !               2. Atmospheric co2
  !               3. sea level
  !               (for eg paleo runs or future scenarios)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine load_forcing(fnm,t_start,t_end)
    
    implicit none
    
    integer, parameter :: nt0 = 50000
    double precision, parameter :: time_err = 9.9e12
    integer :: i, j, q, k, ios, ntime, ntf
    character(len=*) :: fnm
    character(len=10) :: tmp
    character(len=2)  :: ext
    double precision, dimension(nt0) :: time, seal, co2e, co2, dTann
    double precision, dimension(nt0,12) :: tgrl(nt0,12), dTseas
    double precision, dimension(nt0) :: dataset
    double precision :: t_start, t_end
    double precision :: T0(12)

    ! Enter some default values to be safe
    time    = time_err
    seal    = 0.d0
    co2e    = 0.d0
    co2     = 280.d0
    dataset = 1 
    dTann   = 0.d0 
    dTseas  = 0.d0 

    ext = fnm((len(fnm)-1):len(fnm))

    ! Now read in these values from the forcing input file if they exist...
    if ( ext .eq. "nc" .and. index(fnm,"hybrid") .gt. 0) then
      call nc_read(fnm,"time",time(1:6001))
      time(1:6001) = time(1:6001)*1d3
      call nc_read(fnm,"rsl",seal)
      call nc_read(fnm,"co2",co2)
      co2e = 0.d0 
      call nc_read(fnm,"dTann",dTann)
      call nc_read(fnm,"dTseas",dTseas)
      call nc_read(fnm,"dataset",dataset)
      write(*,*) "FORCING: READ FROM HYBRID NC FILE."

      ! Scale temperature anomalies based on parameter scaling factors
      where(dataset .eq. 1.d0 .and. dTann .gt. 0.d0) dTann = dTann * f_hol 
      where(dataset .eq. 2.d0 .and. dTann .lt. 0.d0) dTann = dTann * f_glac

      ! Make sure the Eemian factor is applied to all positive temperatures during the Eemian!
!       where(dataset .eq. 3.d0 .and. dTann .gt. 0.d0) dTann = dTann * f_eem 
      where(time .ge. -130d3 .and. time .le. -115d3 .and. dTann .gt. 0.d0) dTann = dTann * f_eem 

      if (f_eem .eq. -1.d0) then 
        call schema_interglacial_mis5(time,dT=dTann,dT_max=dT_factor,dT_width=dT_width,a1=0.6d0)
      end if 

      do i = 1,12
        tgrl(:,i) = dTann + (dTseas(:,i) * f_seas) 
      end do

!       ! Now modify interglacial summer to match schematic forcing
!       ! (use dTann as object to hold summer average)
!       dTann = sum(tgrl(:,6:8),dim=2) / 3.d0
!       do i = 1,12
!         tgrl(:,i) = tgrl(:,i) - dTann 
!       end do

!       do i = 1,12
!         tgrl(:,i) = tgrl(:,i) + dTann 
!       end do


    else if ( ext .eq. "nc") then 
      call nc_read(fnm,"TIME",time)
      call nc_read(fnm,"gmsl",seal)
      call nc_read(fnm,"co2emi",co2e)
      call nc_read(fnm,"co2", co2)
      call nc_read(fnm,"tgrl",tgrl)
      write(*,*) "FORCING: READ FROM NC FILE."
    else
      call read_forcing_ascii(fnm,time,seal,co2e,co2,tgrl)
      write(*,*) "FORCING: READ FROM ASCII."
    end if
    
    ! Determine end of the actual data contained in arrays
    do k = 1, nt0
      if ( time(k) .eq. time_err ) exit
    end do
    ntf = k - 1
    write(*,*) "ntf = ",ntf

    ! #### Diagnostic output...
!     do k = 1, 100
!      write(*,"(17f10.2)") time(k),seal(k),co2e(k),co2(k),dataset(k),tgrl(k,:)
!     end do
!     write(*,*) "min/max seal =",minval(seal(1:ntf)),maxval(seal(1:ntf))
!     write(*,*) "min/max co2  =",minval(co2(1:ntf)),maxval(co2(1:ntf))
!     stop
    ! ####

    ! Allocate the forcing array for the
    ! length of current simulation time, in years
    ntime = int(t_end - t_start) + 1
    if (allocated(forcing)) deallocate(forcing)
    allocate(forcing(ntime))
    ! Fill time array with yearly values from start to finish
    do k = 1, ntime
      forcing(k)%time = t_start + dble(k-1)
    end do
    
    ! Interpolate the input data, so that a value exists
    ! for every simulation year 

    ! 1d fields
    forcing%rsl  = interp_spline(x=time(1:ntf),y=seal(1:ntf),xout=forcing%time)
    forcing%co2e = interp_spline(x=time(1:ntf),y=co2e(1:ntf),xout=forcing%time)
    forcing%co2  = interp_spline(x=time(1:ntf),y=co2(1:ntf), xout=forcing%time)
     
    ! Now interpolate each month of temperatures
    do q = 1, 12
      forcing%dT(q)  = interp_spline(x=time(1:ntf),y=tgrl(1:ntf,q), xout=forcing%time)
    end do
    
    ! Offset the temperatures to get anomalies
    T0 = forcing_base(forcing_yr0,forcing_yrf)

    do k = 1, size(forcing)
      forcing(k)%dT = forcing(k)%dT - T0 
    end do
    
    if (dT_min .ne. 100.d0) then 
      do i = 1,12

        ! Call separately for MIS-11 and MIS-5 (hard coded time ranges for now...)
        call scale_interglacial_cos(forcing%time,forcing%dT(i),scale=dT_factor, &
                                    dT_min=dT_min,time_min=-425.d3,time_max=-395.d3)

        ! For MIS-5, dT_min is hard-coded to 0deg C, since scaling only applies to positive temp anomalies
        ! Note that the scaling factor dT_factor is applied as (1.0+dT_factor)*dT in this case
        call scale_interglacial(forcing%time,forcing%dT(i),scale=dT_factor, &
                                dT_min=0.d0,time_min=-135.d3,time_max=-115.d3)

      end do

    else if (dT_width .ne. 0.d0 .and. f_eem .ne. -1.d0) then 

      do i = 1,12
        call schema_interglacial(forcing%time,forcing%dT(i),scale=dT_factor,width=dT_width, &
                                 t0=-443.d3,t1=-403.d3,month=i)
!         call schema_interglacial(forcing%time,forcing%dT(i),scale=dT_factor,width=dT_width, &
!                                  t0=-135.d3,t1=-110.d3,month=i)
      end do

    else if (f_eem .ne. -1.d0) then

      ! Adjust the positive paleo boundary temperatures by the dT_factor
      do i = 1,12
        where( forcing%dT(i) .gt. 0.d0 .and. forcing%time .lt. 0.d0 ) forcing%dT(i) = forcing%dT(i) * dT_factor
      end do
 
    end if  

    ! Ouput to ascii file in output folder for checking (every ka)
    open(2,file=trim(outfldr)//"forcing.dat",status="unknown")
    
    write(2,"(a1,16a12)") "#","kyr","m","Gt","ppm", &
                                    "degC","degC","degC","degC","degC","degC",&
                                    "degC","degC","degC","degC","degC","degC"
                                    
    write(2,"(1x,16a12)") "time","rsl","co2e","co2",&
                                  "jan","feb","mar","apr","may","jun",&
                                  "jul","aug","sep","oct","nov","dec"
    do k = 1, ntime
      if ( k .eq. 1 .or. k .eq. ntime .or. &
           mod(forcing(k)%time,1e2) .eq. 0.d0 ) &
             write(2,"(1x,16f12.3)") &
                          forcing(k)%time, forcing(k)%rsl, &
                          forcing(k)%co2e, forcing(k)%co2, &
                          forcing(k)%dT
    end do
    close(2)   

    return
  
  end subroutine load_forcing
  
  subroutine scale_interglacial(time,dT,scale,dT_min,time_min,time_max)
    ! This subroutine will replace the input temperature anomaly time series
    ! within a window of times t0 and t1 with a scaled 'interglacial' forcing
    ! curve. 

    implicit none 

    double precision :: time(:), dT(:)
    double precision :: scale, dT_min, time_min, time_max

!     double precision, parameter :: dT_min   = -5.d0 
!     double precision, parameter :: time_min = -135.d3 
!     double precision, parameter :: time_max = -115.d3 
    
    integer :: k0, k1, k 

    ! Only apply this if interglacial exists in time window 
    if (minval(time) .le. time_min .and. maxval(time) .ge. time_max) then 

      ! Get maximum width of window
      k0 = minloc(dabs(time-time_min),1)
      k1 = minloc(dabs(time-time_max),1)

      ! Get indices above minimum temperature
      do k = k0, k1 
        k0 = k 
        if (dT(k) .ge. dT_min) exit 
      end do 

      do k = k1, k0, -1  
        k1 = k 
        if (dT(k) .ge. dT_min) exit 
      end do 

      if (k1 .lt. k0) then 
        write(*,*) "scale_interglacial:: Error: incorrect indices!"
        write(*,*) k0, time(k0), dT(k0)
        write(*,*) k1, time(k1), dT(k1)
        stop
      end if 
      
      do k = k0, k1

        if (dT(k) .ge. 0.d0) then 
          dT(k) = dT(k) * (1.d0 + scale) 
        else if (dT(k) .ge. dT_min) then 
          dT(k) = dT(k) * (1.d0 + scale*(dT(k)-dT_min)**2/(0.d0-dT_min)**2)
        end if 

      end do 

    end if 

    return 

  end subroutine scale_interglacial 

  subroutine scale_interglacial_cos(time,dT,scale,dT_min,time_min,time_max)
    ! This subroutine will replace the input temperature anomaly time series
    ! within a window of times t0 and t1 with a scaled 'interglacial' forcing
    ! curve. 

    implicit none 

    double precision :: time(:), dT(:)
    double precision :: scale, dT_min, time_min, time_max

!     double precision, parameter :: dT_min   = -2.d0 
!     double precision, parameter :: time_min = -425.d3 
!     double precision, parameter :: time_max = -395.d3 
    
    integer :: k0, k1, k 
    double precision, allocatable :: x1(:), y1(:) 
    integer :: n1

    ! Only apply this if interglacial exists in time window 
    if (minval(time) .le. time_min .and. maxval(time) .ge. time_max) then 

      ! Get maximum width of window
      k0 = minloc(dabs(time-time_min),1)
      k1 = minloc(dabs(time-time_max),1)

      ! Get indices above minimum temperature
      do k = k0, k1 
        k0 = k 
        if (dT(k) .ge. dT_min) exit 
      end do 

      do k = k1, k0, -1  
        k1 = k 
        if (dT(k) .ge. dT_min) exit 
      end do 

      if (k1 .lt. k0) then 
        write(*,*) "scale_interglacial:: Error: incorrect indices!"
        write(*,*) k0, time(k0), dT(k0)
        write(*,*) k1, time(k1), dT(k1)
        stop
      end if 
      
      n1 = k1-k0+1
      allocate(x1(n1),y1(n1))

      do k = 1, n1 
        x1(k) = -pi/2.d0 + (k-1)*pi/dble(n1)
      end do 
      y1 = cos(x1)

      dT(k0:k1) = dT(k0:k1) + maxval(dT(k0:k1))*y1*scale 

    end if 

    return 

  end subroutine scale_interglacial_cos 

  subroutine schema_interglacial(time,dT,scale,width,t0,t1,month)
    ! This subroutine will replace the input temperature anomaly time series
    ! within a window of times t0 and t1 with a schematic and symmetrical 
    ! 'interglacial' forcing curve centered on the original peak of the time series

    implicit none 

    double precision :: time(:), dT(:)
    double precision :: scale, width, t0, t1 
    integer :: month 

    double precision :: dt_inter 

    integer :: k0, k1, kp, k 
    double precision :: x(7), y(7)
    double precision :: tp, dx, dx_max, fp_fill 
    double precision, allocatable :: x1(:), y1(:) 
    integer :: n1, nadd

    character(len=512) :: fname 

    dt_inter = t1 - t0 

    ! Only apply this if interglacial exists in time window 
    if (minval(time) .le. t0 .and. maxval(time) .ge. t0+dt_inter) then 

      ! Specify the initial time and temp anomaly
      x(2) = t0 
      y(2) = interp_linear(time,dT,xout=x(2))

      ! Find indices of temporary time window of interest 
      ! to determine peak time and temp anomaly
      k0 = minloc(abs(time-x(2)),dim=1)
      k1 = minloc(abs(time-(x(2)+dt_inter)),dim=1)
      kp = k0-1 + maxloc(dT(k0:k1),dim=1)
      x(4) = time(kp) !+2.d3
      y(4) = dT(kp)
      if (y(4) .gt. 0.d0) y(4) = y(4)*scale 

      ! Determine actual mirrored final time point (symmetrical to initial point)
      x(6) = x(4) + (x(4)-x(2))
      y(6) = y(2) 

      ! Determine wider buffer points 
      x(1) = x(2) - 3.d3 
      x(7) = x(6) + 5.d3 
      y(1) = interp_linear(time,dT,xout=x(1))
      y(7) = interp_linear(time,dT,xout=x(7))

      ! Determine zero points
      ! Max half-width of interglacial with 1 ka reduced for smoothness
      ! Make sure applied width is <= the max width 
      dx_max = (x(6)-x(2))/2.d0 - 1.d3 
      dx     = min(width*1d3/2.d0,dx_max)    

      x(3) = x(4) - dx 
      y(3) = 0.0
      x(5) = x(4) + dx 
      y(5) = 0.0 

      ! ###### With 'designer' points #######

      ! Get fractional adjustment of distance of fill points
      ! next to peak (proportional to width of interglacial)
      fp_fill = (dx*2.d0)/160.d3
      fp_fill = min(fp_fill,0.8d0)
      fp_fill = max(fp_fill,0.05d0) 

      ! Add fill points to new vector
      nadd = 4
      n1 = size(x) + nadd
      allocate(x1(n1),y1(n1))

      ! Left pinning points
      x1(1:2) = x(1:2)
      y1(1:2) = y(1:2) 

      ! New fill point
      x1(3) = x(2) + (x(3)-x(2))*0.99d0
      y1(3) = y(2) + (y(3)-y(2))*0.99d0  

      ! Left zero point
      x1(4) = x(3)
      y1(4) = y(3)

      ! New fill point
      x1(5) = x(3) + (x(4)-x(3))*(1.d0 - fp_fill)
      y1(5) = y(3) + (y(4)-y(3))*0.985d0  

      ! Peak point 
      x1(6) = x(4)
      y1(6) = y(4)

      ! New fill point
      x1(7) = x(4) + (x(4)-x1(5))
      y1(7) = y1(5)

      ! Right zero point
      x1(8) = x(5)
      y1(8) = y(5)

      ! New fill point
      x1(9) = x(5) + (x(3)-x1(3))
      y1(9) = y1(3)

      ! Right pinning points
      x1(10:11) = x(6:7)
      y1(10:11) = y(6:7)

      ! Reduce the magnitude of the right pinning point and fill point 
      y1(9:10)  = y1(9:10)*0.5d0 

      if (.FALSE.) then 
        ! Ouput to ascii file in output folder for checking (every ka)
        write(fname,"(a,i1,a4)") trim(outfldr)//"forcing_ties_0", month, ".txt"
        if (month .ge. 10) write(fname,"(a,i2,a4)") trim(outfldr)//"forcing_ties_", month, ".txt"
        open(3,file=fname,status="unknown")
        write(3,"(2a12)") "time", "dT" 
        do k = 1, size(x1)
          write(3,"(2f12.2)") x1(k)*1d-3, y1(k) 
        end do 
        close(3)   
      end if 

      ! Find index of initial time point
      do k = 1, size(time)
        k0 = k 
        if (time(k0) .ge. x1(1)) exit 
      end do 

      ! Find index of actual end point to define time window of interest
      k1 = k0
      do k = k0, size(time)
        k1 = k 
        if (time(k1) .ge. x1(size(x1))) exit 
      end do 

      dT(k0:k1) = interp_spline(x1,y1,xout=time(k0:k1))

    end if 

    return 

  end subroutine schema_interglacial

    subroutine schema_interglacial_mis5(time,dT,dT_max,dT_width,a1)
    ! This subroutine will replace the input temperature anomaly time series
    ! (annual mean, precip-weighted annual mean or jja, etc) within a window 
    ! of times t0 and t1 with a schematic 'interglacial' forcing curve

    implicit none 

    double precision :: time(:), dT(:)
    double precision :: dT_max, dT_width, a1 
    integer :: month 

    integer :: k0, k1, k 
    double precision :: x(10), y(10)
    double precision :: alpha

    character(len=512) :: fname 

    x(1) = -131.7d3
    y(1) =   -8.8d0

    x(2) = -130.0d3
    y(2) =   -7.0d0

    x(3) = -128.0d3
    y(3) =   -4.0d0 

    x(4) = -127.2d3 
    y(4) =    0.0d0 

    y(5) = dT_max 
    alpha = a1* (y(4)-y(3)) / (x(4)-x(3))
    x(5) = x(4) + (y(5)-y(4)) * 1/alpha

    x(6) = x(5)   + 0.20d3 
    y(6) = dT_max - 0.02d0

    x(7) = x(4) + min(dT_width*1d3,9.d3)
    x(7) = max(x(7),x(6)+0.8d3)
    y(7) = 0.0 

    x(9) = -116.0d3
    y(9) =   -5.0d0

    x(10) = -113.0d3 
    y(10) =   -8.0d0 

!     x(8) = 0.5d0*x(7) + 0.5d0*x(9)
!     y(8) = 0.5d0*y(7) + 0.5d0*y(9)

    x(8) = x(7) + 0.010d3 
    y(8) = y(7) - 0.02d0 

    ! Only apply this if interglacial exists in time window 
    if (minval(time) .le. minval(x) .and. maxval(time) .ge. maxval(x)) then 

      if (.TRUE.) then 
        ! Ouput to ascii file in output folder for checking
        write(fname,"(a,a)") trim(outfldr)//"forcing_ties.txt"
        open(3,file=fname,status="unknown")
        write(3,"(2a12)") "time", "dT" 
        do k = 1, size(x)
          write(3,"(2f12.2)") x(k)*1d-3, y(k) 
        end do 
        close(3)   
      end if 

      ! Find index of initial time point
      do k = 1, size(time)
        k0 = k 
        if (time(k0) .ge. x(1)) exit 
      end do 

      ! Find index of actual end point to define time window of interest
      k1 = k0
      do k = k0, size(time)
        k1 = k 
        if (time(k1) .ge. maxval(x)) exit 
      end do 

      dT(k0:k1) = interp_spline(x,y,xout=time(k0:k1))

    end if 

    return 

  end subroutine schema_interglacial_mis5

  subroutine stretch_interglacial(time,dT,stretch,t0,t1)
    ! This subroutine will adjust the temperature anomaly time series
    ! within a window of times t0 and t1, such that the points at which 
    ! the anomalies cross to positive anomalies and back down to negative
    ! anomalies are shifted left or right by the number of years 'stretch'

    implicit none 

    double precision :: time(:), dT(:)
    double precision :: stretch, t0, t1 
    
    integer :: k0, k1, kp, k 
    double precision :: x(6), y(6)
    double precision :: tmp(2) 
    double precision, allocatable :: offset(:) 
    double precision, allocatable :: time1(:), dT1(:) 
    integer :: n1

    ! Define end time points (where scaling dissipates to zero)
    x(1) = t0 
    x(6) = t1

    ! Find indices of end time points
    do k = 1, size(time)
      k0 = k 
      if (time(k0) .ge. x(1)) exit 
    end do 

    k1 = k0
    do k = k0, size(time)
      k1 = k 
      if (time(k1) .ge. x(6)) exit 
    end do 

    ! ###### With pre-interpolation #######

    ! Interpolate to higher resolution 
    n1 = (k1-k0+1)*10
    allocate(time1(n1),dT1(n1))
    do k = 1, size(time1)
      time1(k) = time(k0) + (k-1)*(time(k1)-time(k0))/size(time1)
    end do 
    dT1 = interp_linear(time,dT,xout=time1)

    ! Find interglacial maximum within end points 
    ! (within a little spread to avoid a spike at peak)
    kp = maxloc(dT1,dim=1)
    x(3) = time1(kp)-1.0d3 
    x(4) = time1(kp)+1.0d3 

    ! Find first inflection point 
    do k = kp, 1, -1 
      if (dT1(k) .le. 0.d0) exit 
    end do 
    x(2) = time1(k)-stretch*1d3 

    ! Find second inflection point 
    do k = kp, size(time1) 
      if (dT1(k) .le. 0.d0) exit 
    end do 
    x(5) = time1(k)+stretch*1d3

    ! Now define offset values at key points 
    y = 0.d0 

    ! Inflection points 
    tmp = interp_linear(time1,dT1,xout=[x(2),x(5)])
    y(2) = -tmp(1)
    y(5) = -tmp(2) 
    
!     ###### No pre-interpolation #######
!     ! Find interglacial maximum within end points 
!     ! (within a little spread to avoid a spike at peak)
!     k = maxloc(dT(k0:k1),dim=1)
!     kp = (k0-1) + k 
!     x(3) = time(kp)-0.5d3 
!     x(4) = time(kp)+1.0d3 

!     ! Find first inflection point 
!     do k = kp, k0, -1 
!       if (dT(k) .le. 0.d0) exit 
!     end do 
!     x(2) = time(k)-stretch*1d3 

!     ! Find second inflection point 
!     do k = kp, k1 
!       if (dT(k) .le. 0.d0) exit 
!     end do 
!     x(5) = time(k)+stretch*1d3

!     ! Now define offset values at key points 
!     y = 0.d0 

!     ! Inflection points 
!     tmp = interp_linear(time(k0:k1),dT(k0:k1),xout=[x(2),x(5)])
!     y(2) = -tmp(1)
!     y(5) = -tmp(2) 

    write(*,"(a10,6f10.1)") "times: ", x 
    write(*,"(a10,6f10.1)") "  dTs: ", y 
    write(*,*) 

    ! Get offset vector using linear fitting 
    allocate(offset(k1-k0+1))
    offset = interp_linear(x,y,xout=time(k0:k1))

    ! Apply offset 
    dT(k0:k1) = dT(k0:k1) + offset 

    return 

  end subroutine stretch_interglacial

  subroutine read_forcing_ascii(fnm,time,seal,co2e,co2,tgrl)
    
    implicit none

    character(len=* ) :: fnm
    character(len=10) :: tmp
    double precision, dimension(:) :: time,seal,co2e,co2
    double precision, dimension(:,:) :: tgrl 
    integer :: n, k, ios

    double precision :: time0,seal0,co2e0,co20,ecc0,obl0,pre0,tgrl0(12)

    !! READ IN boundary data
    open(13,file=trim(fnm),status="old")
    
    ! Skip through header
    read(13,*) tmp
    read(13,*) tmp
    
    ! Read through file until end of array or end-of-file is reached
    do k = 1, size(time)
      read(13,iostat=ios,fmt=*) &
                              time0, seal0, co2e0, co20, ecc0, obl0, pre0, tgrl0 
      
      if ( ios < 0 ) then
        exit
      else
        time(k)   = time0*1d3
        seal(k)   = seal0
        co2(k)    = co20
        tgrl(k,:) = tgrl0
      end if
      
    end do
    
    ! Adjust co2e 
    co2e = 0.d0 

    n = k-1
    close(13)
    
    return

  end subroutine read_forcing_ascii
  
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  g e t _ f o r c i n g _ a n o m a l i e s
  ! Author     :  Alex Robinson
  ! Purpose    :  Determine the average value of the forcing to
  !               generate temp. anomalies instead of abs. values
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  function forcing_base(time0,timef)
    
    implicit none
    
    double precision :: time0, timef
    double precision :: forcing_base(12)
    double precision :: T0(12)
    integer :: k, nt 
    
    T0 = 0.d0
    nt = 0
    do k = 1, size(forcing)
      if (forcing(k)%time .ge. time0 .and. &
          forcing(k)%time .le. timef) then
        T0 = T0 + forcing(k)%dT
        nt = nt + 1
      end if
    end do
    
    if (nt .gt. 0) T0 = T0 / dble(nt)
    
    write(*,"(a20,12f8.2)") "forcing: base temps =",T0

    forcing_base = T0 
    
    return

  end function forcing_base

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  g e t _ f o r c i n g _ n o w
  ! Author     :  Alex Robinson
  ! Purpose    :  Determine the current values of boundary forcing
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine get_forcing_now(time,mask)

    implicit none
    
    integer :: k, q, dnow, q0, q1, q2, mid, day
    double precision :: time
    double precision :: wt0, wt1, wt2, wttot
    
    double precision, dimension(:,:), optional :: mask
    double precision :: f0, f1, frac
    
    ! Determine index of current year in forcing
    do k = 1, size(forcing)
      if ( forcing(k)%time .ge. time ) exit
    end do
    
    ! Assign current year's data to forcing_now
    forcing_now = forcing(k)
    
    ! Overwrite forcing fields with climber data if coupled
    !if (boundary_forcing .eq. 3 .and. nstep .gt. 0) then
    !  call read_from_climber(trim(outfldr),forcing_now%dT,forcing_now%rsl, &
    !                         forcing_now%co2e,forcing_now%co2)
    !end if

    ! Calculate the fraction of Greenland that's ice covered
    ! and scale the temperature anomaly if desired (if not, set paleo_frac_dT=0)
    if ( present(mask) ) then
      f0 = count(mask .eq. 0.d0)
      f1 = count(mask .le. 1.d0)
      frac = f0/f1
      forcing_now%dTbnd = forcing_now%dT 
      forcing_now%dTamp = (1.d0-frac)*paleo_frac_dT
      forcing_now%dT    = forcing_now%dTbnd + forcing_now%dTamp
    end if
    
    ! Assign temperature anomaly to Tanomaly vector
    ! ####################################################################
    ! Interpolate data in time: monthly => daily
    ! ####################################################################
    dnow = 0

    do q = 1, nm

      do k = 1, int(day_month)

        dnow = dnow+1

        q1 = q                     ! q1 is current month

        q0 = q1-1                  ! q0 is previous month
        if (q0 .eq. 0) q0 = 12     ! Loop to december

        if (q1 .eq. 13) q1 = 1     ! Loop to january for current month

        q2 = q1+1                  ! q2 is next month
        if (q2 .eq. 13) q2 = 1     ! Loop to january for next month

        mid = day_month / 2   ! Halfway point of current month (+/- 1 day)
        day = k - mid

        ! When day < mid, weight shared between wt0 and wt1
        ! When day = mid, weight goes to current month
        ! When day > mid, weight shared between wt1 and wt2
        wt0 = dble (- min(0,day))
        wt1 = day_month - abs(day)
        wt2 = dble (  max(0,day))

        ! Normalize weights
        wttot = dble(wt0 + wt1 + wt2)
        wt0 = wt0 / wttot
        wt1 = wt1 / wttot
        wt2 = wt2 / wttot
        
        Tanomaly(dnow) = forcing_now%dT(q0)*wt0 + &
                         forcing_now%dT(q1)*wt1 + &
                         forcing_now%dT(q2)*wt2
     
      end do
    end do
    
    write(*,"(a,f12.1,2f10.2)") "forcing (time,co2,dTann): ", &
                time, forcing_now%co2, sum(forcing_now%dT)/size(forcing_now%dT)
    return
  
  end subroutine get_forcing_now


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  d e l t a T
  ! Author     :  Alex Robinson
  ! Purpose    :  Determine seasonal value of global warming
  !               around Greenland
  !               summer warming = global warming
  !               winter warming = higher warming
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
  function deltaT_orig(day) result(deltaT)
    
    integer :: day
    real (8) :: deltaT, Twrmng, Tw, Tamp, Tmean
    
    real (8) :: step_frac
    
    ! Fraction of warming to apply until maximum step as a transition from 0
    ! (so that T_warming level applied is not instantaneous)
    step_frac = 1.d0
    if ( T_trans_max .gt. 0 ) &
      step_frac = max( min( dble(nstep-T_warming_delay)/T_trans_max, 1.d0 ), 0.d0 )

    Twrmng = 0.d0  ! Initialize Twrmng as zero
    
    select case(slow_hyst)
      case(1,-1)
        Twrmng = dTtrans(nstep)
      case(2,-2)
        Twrmng = T_warming + dTtrans0(nstep)
      case(3)
        Twrmng = T_warming + Tnoise(nstep)
      case(4,-4)
        Twrmng = dTtrans1(nstep)
      case DEFAULT
        if ( nstep .ge. T_warming_delay ) Twrmng = T_warming *step_frac

    end select

    Tw      = Twrmng * T_wintfac     ! Default T_wintfac is 2 (winter 2x warmer than summer)
    Tamp    = (Tw - Twrmng) / 2.d0
    Tmean   = (Tw + Twrmng) / 2.d0
    
    if ( day .eq. 0 ) then  ! provide the global (Greenland summer) value
      deltaT = Twrmng  
      
      ! Add summer value from paleo if needed
      if ( boundary_forcing .gt. 0 ) deltaT = deltaT + sum(forcing_now%dT(6:8)) / 3d0
        
    else                    ! provide the local Greenland day's warming
      deltaT = Tmean + Tamp * dcos(2.d0*pi*(day-15)/day_year)
       
      ! Add daily value from paleo if needed
      if ( day .gt. nk ) day = day - nk
      if ( boundary_forcing .gt. 0 ) deltaT = deltaT + Tanomaly(day)
    end if
    
    return

  end function deltaT_orig
  
  function deltaT(day) result(dT_out)
    
    integer  :: day
    real (8) :: dT_out      ! Warming for the given day, or summer value

    ! Local variables
    real (8) :: T_summer  
    real (8) :: Tw, Tamp, Tmean
    real (8) :: step_frac
    
    ! Update T_summer with global forcing dT value 
    T_summer = T_warming_in 
    
    Tw      = T_summer * T_wintfac     ! Default T_wintfac is 2 (winter 2x warmer than summer)
    Tamp    = (Tw - T_summer) / 2.d0
    Tmean   = (Tw + T_summer) / 2.d0
    
    if ( day .eq. 0 ) then  ! provide the global (Greenland summer) value
      dT_out = T_summer  
      
      ! Add summer value from paleo if needed
      if ( boundary_forcing .gt. 0 ) dT_out = dT_out + sum(forcing_now%dT(6:8)) / 3d0
        
    else                    ! provide the local Greenland day's warming
      dT_out = Tmean + Tamp * dcos(2.d0*pi*(day-15)/day_year)
       
      ! Add daily value from paleo if needed
      if ( day .gt. nk ) day = day - nk
      if ( boundary_forcing .gt. 0 ) dT_out = dT_out + Tanomaly(day)
    end if
    
    return

  end function deltaT
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  d t T r a n s 0
  ! Author     :  Alex Robinson
  ! Purpose    :  Generate correct T_warming for gradual changes over
  !               time (continuous stability diagram!)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
  function dTtrans0(n_step)
    
    integer :: n_step
    real (8) :: Twrmng, year
    double precision :: rate
    double precision :: dTtrans0
    
    year = dble(n_step)
    
    !rate = T_diff / 1d6    ! 2.0 deg C / 1 million years
    rate = dT_rate ! K / a
       
    if (slow_hyst .eq. -2) rate = -rate  ! cooling for ice-free init.
    
    ! Make sure warming should be applied...
    dTtrans0 = 0.d0
    if ( n_step .ge. T_warming_delay ) dTtrans0 = rate*(year - T_warming_delay)
    
    if ( trim(domain) .eq. "calov-island" ) then
    ! (reinhard)
    ! Faster increasing temperature for the first 100 years)
      rate = 5.d0 / 100.d0   ! 5 deg C / 100 years, misused for transient global warming runs

      if(n_step.le.100) then
        dTtrans0 = rate*dble(n_step)
      else
        dTtrans0 = rate*100.d0
      end if
    end if
    
    return

  end function dTtrans0
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  t r a n s T w a r m i n g 2
  ! Author     :  Alex Robinson
  ! Purpose    :  Generate correct T_warming for gradual changes over
  !               time (continuous stability diagram!)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
  function dTtrans(n_step)
    
    integer :: n_step, n_tot, n_yr    
    double precision :: y, dx
    double precision :: dT1, dT2, dT_mid, dT_min, dT_max, dT_diff
    double precision :: dTtrans
    
    ! Get the current n and the total n
    n_yr  = n_step
    n_tot = yearf-year0
    
    ! Set the temperature range to be applied (hard-coded for now),
    ! and set the threshold to the input T_warming parameter
    dT_min  = -1.5d0; dT_max = 3.5d0
    dT_mid  = T_warming
    dT_diff = T_diff
    
    ! Set up current generic cubic y value that fits inside a range from -1 to 1
    dx = (1.d0-(-1.d0))/dble(n_tot)
    y = ( -1.d0 + max(n_yr-1,0)*dx )**3

    if (slow_hyst .eq. -1) then
      y = -y                   ! cooling for ice-free init.
      dT_mid = dT_mid - dT_diff   ! Threshold is 1.45 degree lower (for ice-free init.)
      dT_min = -0.5d0
    end if
    
    if ( dT_mid .lt. dT_min .or. dT_mid .gt. dT_max ) then
      write(stdout,*) "dTtrans:: threshold dT_mid is not inside bounds of temperatures... try again!"
      call embout(1)
    end if
    
    ! Get scaling ranges for y >= 0 and y < 0
    ! (based on where the threshold is, dT_mid)
    dT2 = dT_max - dT_mid
    dT1 = dT_mid - dT_min

    ! Get the temperature anomaly for the current year
    ! (scaled by whether we are before or after threshold)
    if ( y .ge. 0.d0 ) then
      dTtrans = dT2*y + dT_mid
    else
      dTtrans = dT1*y + dT_mid
    end if
    
    return

  end function dTtrans
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  d t T r a n s 1
  ! Author     :  Alex Robinson
  ! Purpose    :  Generate correct T_warming for gradual changes over
  !               time (continuous stability diagram!)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
  function dTtrans1(n_step)
    
    integer :: n_step
    double precision :: dTtrans1
    
!     double precision, parameter :: dTdt_min = 1.d0/1d6      ! 1deg/million years
!     double precision, parameter :: dTdt_max = 4.d0/1d6      ! 2deg/million years
!     double precision, parameter :: dVdt_max = 30.d0  ! Gt/a
    double precision :: dTdt_min, dTdt_max
    double precision :: dVdt_now, dVdt_max, dTdt, Tlim
    double precision :: dVdtm_now
    
    double precision :: fac, expfac
!     double precision, parameter :: fac = 5.d0
!     double precision, parameter :: expfac = dexp(fac)
    
    fac = h_fac
    expfac = dexp(fac)
    
    dVdt_max = h_dVdt_max               ! 30 Gt/a, more or less...
    dTdt_min = h_dTdt_min/1d6           ! 1deg/million years
    dTdt_max = h_dTdt_max/1d6           ! 4deg/million years
    
    if (n_step .eq. 0) then
      
      transT%T_warming = T_warming
      dTtrans1         = transT%T_warming
      transT%nstep = 0
      transT%dVdt = 0.d0
      
      transT%dVdtm = 0.d0
      transT%mstep = 0
      
    else if ( n_step .eq. transT%nstep ) then
    
      dTtrans1         = transT%T_warming
      
    else
      
      transT%mstep = transT%mstep+1
      if (transT%mstep .gt. size(transT%dVdtm)) transT%mstep = 1
      transT%dVdtm(transT%mstep) = transT%dVdt
      
      ! Get the absolute value of the current change in vol (Gt/a)
      ! (not greater than the max allowed)
      !       dVdt_now = min(dabs(transT%dVdt),dVdt_max)
      dVdtm_now = sum(transT%dVdtm,transT%dVdtm .ne. 0.d0)/max(1,count(transT%dVdtm .ne. 0.d0))
      dVdt_now  = min(dabs(dVdtm_now),dVdt_max)
 
      ! Calculate the current dTdt based on allowed values and rate of smb 
      ! BASED ON COS (smoother transition)
      !dTdt = (dTdt_max-dTdt_min)*0.5d0*(dcos(pi*dVdt_now/dVdt_max)+1.d0) + dTdt_min
      
      ! BASED ON EXPONENTIAL (sharper transition, tuneable)
      dTdt = (dTdt_max-dTdt_min)*(1.d0-(dexp(dVdt_now/dVdt_max*fac)-1.d0)/(expfac-1.d0)) + dTdt_min
      
      ! Cooling for ice-free init.
      if (slow_hyst .lt. 0) dTdt = -dTdt
      
      
      ! Make sure warming should be applied...
      if ( nstep .ge. T_warming_delay ) then
        ! Assign new T_warming values
        transT%T_warming = transT%T_warming + dTdt*(n_step-transT%nstep)
      end if
      
      ! Output current T_warming value
      dTtrans1 = transT%T_warming
      
      ! Activate kill switch if max/min temperature has been reached
      if (slow_hyst .gt. 0) then  ! T_warming is increasing
        Tlim = T_warming + T_diff
        if ( dTtrans1 .gt. Tlim ) kill = 1
      else                        ! T_warming is decreasing
        Tlim = T_warming - T_diff
        if ( dTtrans1 .lt. Tlim ) kill = 1
      end if
      
      write(stdout,"(a5,3f10.3,3g15.3)") "dTdt ", n_step*1d-3, dTtrans1, dTdt*1e6, &
                              transT%dVdt, dVdt_now, dVdtm_now
    
      transT%nstep = n_step
      transT%dTdt  = dTdt*1e6
      
    end if
    
    return

  end function dTtrans1
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  T n o i s e
  ! Author     :  Alex Robinson
  ! Purpose    :  Generate sinusoidal temperature noise to 
  !               mimic natural variability
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
  function Tnoise(n_step)
    
    integer :: n_step
    real (8) :: Tnoise, year
    
    year = dble(n_step)
    
    Tnoise = T_noise * dsin(2.d0*pi*year/T_noise_period)
      
    return

  end function Tnoise
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c o 2
  ! Author     :  Alex Robinson
  ! Purpose    :  Determine atmos. concentration of co2 based
  !               on temperature change relative to preindustrial
  !               and chosen climate sensitivity
  !               If boundary_forcing is active (1), then dT is ignored
  !               and co2 conc. is obtained from boundary forcing data
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
  function co2(dT)
    
    real (8) :: co2, dT
    
    real (8), parameter :: co2_0    = 280.d0
    real (8), parameter :: ln2      = dlog(2.d0)

    if ( boundary_forcing .eq. 1 .or. boundary_forcing .eq. 3 ) then ! (co2 from file or climber)
      
      ! Get amount of CO2 from anomaly data
      co2 = forcing_now%co2
      
    else if ( boundary_forcing .eq. 0 .or. boundary_forcing .eq. 2) then  ! (prescribed T_anomaly but not co2)
      
      ! First get amount of CO2 based on dT
      co2 = co2_0 * dexp ( dT * ln2 / clim_sens )
      
    end if
    
    return

  end function co2
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  R c o 2
  ! Author     :  Alex Robinson
  ! Purpose    :  Determine the radiative forcing for certain 
  !               co2 concentration
  !               (Based on paper by ? - taken from Andrey)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
  function Rco2(co2)
    
    real (8) :: Rco2, co2
    
    real (8), parameter :: co2_0    = 280.d0
    real (8), parameter :: Rco2_fac = 5.35d0 
    
    !if ( co2_0 .ne. 280.d0 ) write(*,*) "co2_0: ",co2_0
    !write(*,*) "co2: ",co2,", co2_0: ",co2_0
    !write(*,*) "dlog( co2/co2_0 ): ",dlog( co2 / co2_0 )
    Rco2 = Rco2_fac * dlog( co2 / co2_0 )
  
    return

  end function Rco2
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  a i r d e n s i t y
  ! Author     :  Alex Robinson
  ! Purpose    :  Determine the air density depending on elevation
  ! zs [m]
  ! rho [kg/m3]
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
  function airdens(zs)
    
    real (8) :: zs, airdens
    
    airdens = 1.3d0 * exp(-zs/8600.d0)
  
    return
  end function airdens
        

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  e m b _ F Q S A T
  ! Purpose    :  Same as FQSAT in climber, except real (8),
  !               returns sat. specific humidity, as a function 
  !               of °C
  ! Author     :  Alex Robinson (24. June 2008)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function emb_FQSAT(T)  
    
    implicit none
        
    real (8) :: T, emb_FQSAT
    
    emb_FQSAT=3.8d-3*EXP(17.67d0*T/(T+243.5d0))

    return
    
  end function emb_FQSAT
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  s i n s o l 2 d
  ! Purpose    :  Same as climber module sinsol, except it is
  !               based on a 2d grid 
  ! Author     :  Alex Robinson (24. June 2008)
  ! =Input======  
  !  S0        - solar constant, normally 1365.d0 W / m2
  !  BTIME     - astronomical time (=NYRA), present day 1950 = 0.0
  !  LATS      - 2d array of latitude values, in radians!!
  ! =Output======
  !  SOLARM2D  - daily solar radiation at the top of the atmosphere
  !	 COSZM2D   - daily averaged solar zenith angle
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sinsol2d(solarm2d,lats,btime0)

    implicit none
    
    integer, parameter :: nx = nxe, ny= nye, nl = 50
    real (8) :: lats(ny,nx)
    real (8) :: lats0(nl), solarm(nk,nl), coszm(nk,nl)
    real (8) :: solarm2d(nk,ny,nx), coszm2d(nk,ny,nx)
    real (8) :: btime0
    
    integer :: i, j, k, h, jj1, jj0
    real (8) :: s, cosn, cosp, fi
    
    real (4) :: BTIME, ECC, XOBCH, TPERI, ZAVEXPE
    real (4) :: PCLOCK, PYTIME, PDISSE, PZEN1, PZEN2, PZEN3, PRAE

    ! If incoming time is in calendar years, shift it for input to 
    ! sinsol, otherwise leave it as it is
    BTIME = btime0
    if (timetype .eq. 1) BTIME = btime0 - 1950.d0
    
    if ( init_sinsol .eq. 0 ) then
      call INI_SINSOL("./input/")
      init_sinsol = 1
    end if
    
    ! Populate the lats lookup table    
    lats0(1)  = floor(minval(lats))
    lats0(nl) = ceiling(maxval(lats))
    
    fi = ( lats0(nl) - lats0(1) ) / nl
    
    do j = 2, nl-1
      lats0(j) = lats0(j-1) + fi
    end do
    
!...1) Berger program calculates orbital parameters
!   ===========================================   
    call BERGOR(BTIME,ECC,XOBCH,TPERI,ZAVEXPE)
!   =========================================== 

!...2) Daily insolation is calculated by hourly integration for each day
     
    ! Loop over interpolation matrix (lats0, solarm, coszm)
    do j = 1, nl

        fi=lats0(j)
        do k = 1, nk
          PYTIME=k*2.d0*pi/day_year
          solarm(k,j)=0.d0
          coszm(k,j) =0.d0        
          do h = 1, 24
            PCLOCK=h*2.d0*pi/24.d0   
!     =================================================================          
            call ORBIT(ECC,XOBCH,TPERI,ZAVEXPE, &
                 PCLOCK,PYTIME,PDISSE,PZEN1,PZEN2,PZEN3,PRAE)
!     =================================================================                    
            cosp=PZEN1*dsin(fi)+PZEN2*dcos(fi)
            cosn=max(cosp,0.0)
      
            s=s0*cosn*PDISSE
            solarm(k,j)=solarm(k,j)+s
            coszm(k,j)=coszm(k,j)+s*cosn
          enddo
        enddo

      
!...  Daily insolation and zenite angle
                                 

        do k = 1, nk
          solarm(k,j)=solarm(k,j)/24.d0
          if (solarm(k,j) .gt. 0.) then
            coszm(k,j)=coszm(k,j)/(solarm(k,j)*24.d0)
          else
            coszm(k,j)=0.d0
          endif
        enddo
      
    enddo ! end j-loop
    
    ! Now interpolate to fill in the gaps
    do i = 1, nx
      do j = 1, ny
      
        fi = lats(j,i)   ! Store current latitude
        
        do k = 1, nl
          if ( lats0(k) .gt. fi ) then
            jj0 = k-1; exit
          end if
        end do
        
        do k = jj0, nl
          if ( lats0(k) .ge. fi ) then
            jj1 = k; exit
          end if
        end do
        
        do k = 1, nk
          call interp1(fi, solarm2d(k,j,i), lats0(jj0), solarm(k,jj0), lats0(jj1), solarm(k,jj1))
          if ( solarm2d(k,j,i) .lt. 1000.d0 .and. solarm2d(k,j,i) .ge. 0.d0 ) then
          
          else
            write(*,*) "solarm2d not in range",k,j,i
            write(*,"(4f12.2)") lats0(jj0), solarm(k,jj0), lats0(jj1),solarm(k,jj1)
            write(*,"(2f12.2)") fi, solarm2d(k,j,i)
            stop
          end if
          
        end do
        
      end do
    end do

    return
  
  end subroutine sinsol2d 

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  t o h i r e s
  ! Author     :  Alex Robinson
  ! Purpose    :  interpolate a cartesian grid to a higher resolution
  !               Note: assumes lower res is a multiple of higher res!
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  subroutine tohires(v_in,v_out,ratio,m)
    
    implicit none
    
    integer,  parameter :: nx = nxs, ny = nys
    
    integer,  parameter :: nxl = nxe, nyl = nye
    
    real (8) :: v_in(:,:), v_out(:,:)
!     real (8) :: v_in(nyl,nxl), v_out(ny,nx)
    
    integer :: nin, i, j, inow, jnow, ih, jh, jh2, q, m
    real (8) :: ratio, wt, frac
    
    nin = ceiling(1/ratio) - 1 
    frac = 1.d0 / (nin+1)
    
    v_out = 0.d0
    
    ! Loop through lo-res grid
    do inow = 1, nxl-1
      ih = 1 + ceiling((inow-1)/ratio)
      
      do jnow = 1, nyl
        jh = 1 + ceiling((jnow-1)/ratio)
        
        ! Save matching points horizontally
        v_out(jh,ih)       = v_in(jnow,inow)
        v_out(jh,ih+nin+1) = v_in(jnow,inow+1)
        
        ! Fill in gaps horizontally
        do q = 1, nin
          wt = q*frac
          v_out(jh,ih+q) = (1.d0-wt)*v_in(jnow,inow) + wt*v_in(jnow,inow+1)
        end do
        
      end do
    end do
    
    ! Loop through hi-res grid
    do ih= 1, nx
      do jnow= 1, nyl-1
        
        jh  = 1 + ceiling((jnow-1)/ratio)
        jh2 = 1 + ceiling(jnow/ratio)
        
        ! Fill in gaps vertically
        do q = 1, nin
          wt = q*frac
          v_out(jh+q,ih) = (1.d0-wt)*v_out(jh,ih) + wt*v_out(jh2,ih)
        end do
        
      end do
    end do      
    
    if (m .eq. 1) then
      write(stdout,*) "Code to fix mask!" ! Not yet needed
    end if  
    
    return
  end subroutine tohires
     
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  t o l o r e s
  ! Author     :  Alex Robinson
  ! Purpose    :  Aggregate a cartesian grid to a lower resolution
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine tolores(v_in,v_out,wts,nr,ratio,m)
      
    implicit none
    
    integer :: nx, ny, nxl, nyl
    
    real (8) :: v_in(:,:), v_out(:,:)
    
    integer :: nd, nr, n0
    integer :: i, j, inow, jnow, x, y, cx, cy
    character (len=*), optional :: m
    real (8) :: wts(:,:), ratio
    real (8), allocatable :: mask(:,:), pts(:,:)
    
    integer :: n_wtr, n_land, n_ice
    
    nx  = size(v_in,2)
    ny  = size(v_in,1)
    nxl = size(v_out,2)
    nyl = size(v_out,1)
    
    ! Diameter of neighborhood
    nd = nr*2 + 1
    
    ! Center index of neighborhood
    n0 = nr+1
    
    if (allocated(mask)) deallocate(mask)
    if (allocated(pts)) deallocate(pts)
    allocate(mask(nd,nd), pts(nd,nd))
          
    ! Average points in neighborhood
    
    ! Cycle through low res grid
    do inow = 1, nxl
      do jnow = 1, nyl
          
        cx = 1 + int((inow-1)/ratio)
        cy = 1 + int((jnow-1)/ratio)

        ! Fill in neighbors of current point from high res grid
        do i=1,nd
          
          x = cx - n0 + i
          x = max(1,x)
          x = min(nx,x)
          
          do j=1,nd

            y = cy - n0 + j
            y = max(1,y)
            y = min(ny,y)

            pts(j,i) = v_in(y,x)
          
          end do
        end do
        
        ! Get the value of the current point
        v_out(jnow,inow) = sum(pts*wts)

      end do
    end do
    
    if (present(m)) then
    ! If manipulating a mask, ensure that integer
    ! quality is retained.
    ! Simply set current point on lores grid to corresponding
    ! point on hires grid       
      do inow = 1, nxl
        do jnow = 1, nyl
        
          cx = 1 + int((inow-1)/ratio)
          cy = 1 + int((jnow-1)/ratio)
          
          v_out(jnow,inow) = v_in(cy,cx)
          
          !if (v_out(jnow,inow) .lt. 0.5d0) then
          !  v_out(jnow,inow) = 0.d0
          !else if (v_out(jnow,inow) .lt. 1.5d0) then
          !  v_out(jnow,inow) = 1.d0
          !else
          !  v_out(jnow,inow) = 2.d0
          !end if
          
        end do
      end do
    end if
    
    return
      
  end subroutine tolores

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine :  t o l o w t s
    ! Author     :  Alex Robinson
    ! Purpose    :  Get neighbor weights by inverse distance
    !               on a cartesian grid
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine tolowts(wts,nr,method)
      
      implicit none
      
      integer :: i, j, c, nr, nd    ! Radius in points
      real (8) :: wts(:,:)

      real (8) :: dx, dist, distmax, factor
      character (len=*), optional :: method
      character (len=10) :: meth

      meth = "exp"
      if (present(method)) meth = method
      write(stdout,*) "weighting method: ", trim(meth)
      
      c = nr+1
      nd = nr*2+1
      
      dx = 1.d0
      distmax = dsqrt( (nr*dx)**2 + (nr*dx)**2 )
      factor = 4.d0 / distmax
      
      if (trim(meth) .eq. "idw") then
        do i = 1, nd
          do j = 1, nd
          
            dist = dsqrt( ((i-c)*dx)**2 + ((j-c)*dx)**2 )
            
            if (dist .lt. dx) dist = dx/5.d0    ! Make sure dist is small, but > 0
            wts(j,i) = 1 / dist**2           ! Assign a weight

          end do
        end do
        
        !wts(c,c) = 0.d0
        !wts(c,c) = sum(wts)*factor
      
      else if (trim(meth) .eq. "exp") then
        do i = 1, nd
          do j = 1, nd
          
            dist = dsqrt( ((i-c)*dx)**2 + ((j-c)*dx)**2 )
            
            wts(j,i) = exp(-dist*factor) ! Assign a weight
            
          end do
        end do
      
      
      else
        write(stdout,*) "No appropriate weighting method given... try again!"
        call embout(1)
      end if
      
      wts = wts / sum(wts)

      return
      
    end subroutine tolowts
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  m o n t h l y 2 d a i l y
  ! Author     :  Alex Robinson
  ! Purpose    :  Return indices of months and weights
  !               to be used for linear daily interpolation
  !               from monthly values
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  subroutine monthly2daily(m0,wt0,m1,wt1)
  
    implicit none
    
    integer :: k, m, n, mid
    integer, dimension(:) :: m0, m1
    double precision, dimension(:) :: wt0, wt1
    double precision :: wttot
    
    ! Get length of arrays to fill, midpoint of month (in days)
    ! and the total weight (total number of days in a month)
    n = size(m0)
    mid = ndm / 2
    
    ! ####################################################################
    ! Interpolate data in time: monthly => daily
    ! ####################################################################   
    do k = 1, n
      
      do m = 1, nm+1
        if ( m*ndm-mid .gt. k ) exit
      end do
      m1(k) = m; m0(k) = m1(k)-1
      
      if ( m1(k) .gt. 12 ) m1(k) =  1
      if ( m0(k) .lt.  1 ) m0(k) = 12
      
      wt1(k) = dble ( abs( mod(k-mid,ndm) )  )
      wt0(k) = day_month - wt1(k)
      
      wttot = wt0(k) + wt1(k)
      wt0(k) = wt0(k)/wttot; wt1(k) = wt1(k)/wttot
      
    end do
            
    return
  
  end subroutine monthly2daily

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  h g r a d
  ! Author     :  Alex Robinson
  ! Purpose    :  Horizontal gradient of u0
  !               Returns gradient in km / km or m / m, depends on input
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
  subroutine hgrad(u0,dx,uu)  
  
    implicit none
    
    integer, parameter :: nx = nxe, ny = nye
    
    integer :: i, j
    real (8), dimension(ny,nx), intent(IN) :: u0
    real (8), dimension(ny,nx) :: uu
    real (8) :: dx, inv_2dx, dudx, dudy
    
    inv_2dx = 1.d0 / (2.d0 * dx)
    
    ! Assign boundary values
    uu(1:ny,1)  = 0.d0
    uu(1:ny,nx) = 0.d0
    uu(1,1:nx)  = 0.d0
    uu(ny,1:nx) = 0.d0
    
    do i = 2, nx-1
      do j = 2, ny-1
      
        !uu(j,i) = dsqrt( inv_4dx2 * (u0(j,i+1) - u0(j,i-1))**2 + &
        !                 (u0(j+1,i) - u0(j-1,i))**2 )
        
        dudx = (u0(j,i+1) - u0(j,i-1)) * inv_2dx
        dudy = (u0(j+1,i) - u0(j-1,i)) * inv_2dx
        
        uu(j,i) = dsqrt( dudx**2 + dudy**2 )
        
      end do
    end do
    
    return
    
  end subroutine hgrad
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  d 2 u _ d x 2
  ! Author     :  Alex Robinson
  ! Purpose    :  Quick calculation of 2nd derivative in xy
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
  subroutine d2u_dx2(u0,dx,uu)  
  
    implicit none
    
    integer, parameter :: nx = nxe, ny = nye
    
    integer :: i, j
    real (8), dimension(ny,nx), intent(IN) :: u0
    real (8), dimension(ny,nx) :: uu
    real (8) :: dx, inv_dx2
    
    inv_dx2 = 1.d0 / (dx**2)
    
    ! Assign boundary values
    uu(1:ny,1)  = 0.d0
    uu(1:ny,nx) = 0.d0
    uu(1,1:nx)  = 0.d0
    uu(ny,1:nx) = 0.d0
    
    do i = 2, nx-1
      do j = 2, ny-1
      
        uu(j,i) = inv_dx2 * ( u0(j,i+1) + u0(j,i-1) + u0(j+1,i) &
                  + u0(j-1,i) - 4.d0 * u0(j,i) ) 
      
      end do
    end do
  
    return
  end subroutine d2u_dx2

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i n t e r p 1
  !   Author     :  Alex Robinson
  !   Purpose    :  Interpolates for the y value at the desired x value, 
  !                 given x and y values around the desired point.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine interp1(x, y, x0, y0, x1, y1)
  
    implicit none
    
    integer :: i
    real (8) :: x, y, x0, y0, x1, y1, alph
    
    if ( x .lt. x0 .or. x .gt. x1 ) then
      write(*,*) "interp1: x < x0 or x > x1 !"
      write(*,*) "x = ",x
      write(*,*) "x0 = ",x0
      write(*,*) "x1 = ",x1
      stop
    end if

    alph = (x - x0) / (x1 - x0)
    
    y = y0 + alph*(y1 - y0)
    
    return

  end subroutine interp1 
  
  subroutine interp_series(x0,y0,x1,y1)
  
    implicit none
    
    integer :: k, n0, n1, q
    double precision, dimension(:) :: x0, y0, x1, y1
    
    n0 = size(x0)
    n1 = size(x1)

    ! Interpolate x0,y0 data to fill in the data on new values, x1,y1
    do k = 1, n1

      if ( x1(k) .le. x0(1) ) then        ! Use first value
      
        y1(k) = y0(1)
      
      else if ( x1(k) .ge. x0(n0) ) then   ! Use last value

        y1(k) = y0(n0)

      else                                ! Interpolate to get value
        
        ! Find the index equal to or greater than current value
        do q = 1, n0
          if ( x0(q) .ge. x1(k) ) exit
        end do
        
        ! Interpolate to get y1 values at x1
        call interp1(x1(k),y1(k), x0(q-1),y0(q-1),x0(q),y0(q) )

      end if     
      
    end do
    
    ! Output some checks
    write(*,*) "Finished interpolating time series:"
    write(*,"(a20,2g12.3,1x,a1,1x,2g12.3)") "min/max x0 ; y0 :",minval(x0),maxval(x0),";",minval(y0),maxval(y0)
    write(*,"(a20,2g12.3,1x,a1,1x,2g12.3)") "min/max x1 ; y1 :",minval(x1),maxval(x1),";",minval(y1),maxval(y1)

    return

  end subroutine interp_series

!   ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   !   Subroutine :  p o l y f i t
!   !   Author     :  http://rosettacode.org/wiki/Polynomial_regression
!   !   Purpose    :  Returns polynomial coefficients to a polynomial 
!   !                 regression
!   ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   function polyfit(vx, vy, d)
!     
!     implicit none
!     
!     integer, intent(in)                   :: d
!     real (8), dimension(d+1)              :: polyfit
!     real (8), dimension(:), intent(in)    :: vx, vy
!  
!     real (8), dimension(:,:), allocatable :: X
!     real (8), dimension(:,:), allocatable :: XT
!     real (8), dimension(:,:), allocatable :: XTX
!  
!     integer :: i, j
!  
!     integer :: n, lda, lwork
!     integer :: info
!     integer,  dimension(:), allocatable :: ipiv
!     real (8), dimension(:), allocatable :: work
!  
!     n = d+1
!     lda = n
!     lwork = n
!  
!     allocate(ipiv(n))
!     allocate(work(lwork))
!     allocate(XT(n, size(vx)))
!     allocate(X(size(vx), n))
!     allocate(XTX(n, n))
!  
!     ! prepare the matrix
!     do i = 0, d
!        do j = 1, size(vx)
!           X(j, i+1) = vx(j)**i
!        end do
!     end do
!  
!     XT  = transpose(X)
!     XTX = matmul(XT, X)
!  
!     ! calls to LAPACK subs DGETRF and DGETRI
!     call DGETRF(n, n, XTX, lda, ipiv, info)
!     if ( info /= 0 ) then
!        print *, "problem"
!        return
!     end if
!     call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
!     if ( info /= 0 ) then
!        print *, "problem"
!        return
!     end if
!  
!     polyfit = matmul( matmul(XTX, XT), vy)
!  
!     deallocate(ipiv)
!     deallocate(work)
!     deallocate(X)
!     deallocate(XT)
!     deallocate(XTX)
!  
!   end function
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  x y 2 i n d e x
  ! Author     :  Alex Robinson
  ! Purpose    :  Get the i,j indices of a specific x,y value
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine xy2index(grd,x,y,inow,jnow)
    
    implicit none
  
    type(location) :: grd(nys,nxs)
    integer :: nx, ny
    integer :: inow, jnow, k
    real (8) :: x, y, dx, dy
    
    nx = size(grd,2)
    ny = size(grd,1)
    
    do inow = 1, nx
      if (grd(1,inow)%x .ge. x) exit
    end do
    
    do jnow = 1, ny
      if (grd(jnow,1)%y .ge. y) exit
    end do   
    
    if ( inow .eq. 1 ) write(*,*) grd(1,inow)%x, ":",x
    if ( jnow .eq. 1 ) write(*,*) grd(jnow,1)%y, ":",y
    
    if (inow .gt. 1) then 
      dx = (grd(1,inow)%x-x) / (grd(1,inow)%x-grd(1,inow-1)%x)
      if (dx .ge. 0.5) inow = inow-1
    end if 

    if (jnow .gt. 1) then 
      dy = (grd(jnow,1)%y-y) / (grd(jnow,1)%y-grd(jnow-1,1)%y)
      if (dy .ge. 0.5) jnow = jnow-1
    end if 

    return
  
  end subroutine xy2index
        
! #### WRITING OUTPUT FUNCTIONS ####

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Program    :  w r i t e c h e c k
  ! Purpose    :  Write an array to stdout for diagnostics
  ! Author     :  Alex Robinson (03 Jul 2008)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
  subroutine writecheck(unt,var,fmt1,title)

    real (8) :: var(:,:)
    integer :: nx, ny
    character (len=*) :: fmt1, title
    character (len=15) :: fmt2, c_nx
    
    integer :: i, j, unt
    
    nx = size(var,2)
    ny = size(var,1)
    
    write(unt,*) nx
    write(c_nx,"(i3)") nx
    c_nx=adjustl(c_nx)
    fmt2 = "("//trim(c_nx)//trim(fmt1)//")"
    
    write(unt,*)
    write(unt,*) title
    do j = ny, 1, -1
      write(unt,fmt2) (var(j,i), i=1, nx)
    end do

    return
  
  end subroutine writecheck
  
  
end module emb_functions

! #### READING INPUT FUNCTIONS ####
      
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l o a d _ v a r
  !   Author     :  Alex Robinson
  !   Purpose    :  load a generic variable from standard ascii file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  subroutine load_var(ny,nx,var,dattype,fnm)
  
    implicit none
    
    integer :: ny,nx
    real (8) :: var(ny,nx)
    integer :: i, j, check, intvar(ny,nx)
    character(len=50) :: fmt1
    character(len=*) :: dattype, fnm

    open(24,iostat=check,file=trim(fnm),status="old")
    if (check .ne. 0) then
      write(*,*) "Error opening file: ",trim(fnm)
    end if

    do j=1,6
      read(24,*)
    end do
    
    if (dattype .eq. "dble") then
      write(fmt1,"(a1,i3,a6)") "(",nx,"f12.3)"
      do j=ny, 1, -1
        read(24,fmt1) (var(j,i), i=1,nx)
      end do
    else if (dattype .eq. "int") then
      write(fmt1,"(a1,i3,a3)") "(",nx,"i1)"
      do j=ny, 1, -1
        read(24,fmt1) (intvar(j,i), i=1,nx)
      end do
      var = dble(intvar)
    else
      write(fmt1,"(a1,i3,a6,a1)") "(",nx,dattype,")"
      do j=ny, 1, -1
        read(24,fmt1) (var(j,i), i=1,nx)
      end do
    end if 
    
    close(24)
  
    return
  
  end subroutine load_var

