module rembo_functions
  
  use emb_global
  use emb_functions
  use ncio 

  implicit none

  type rembo_type
    double precision, dimension(nys,nxs) :: &
             tt, tte, pp, snow, snowh, ice, runoff, runoff_snow, runoff_rain, S, dh, &
             rf, as, ap, refrozen, melt, melted_ice, melt_insol, melt_S0, pdd_corr
    double precision, dimension(nys,nxs) :: &
             ap0, S0
  end type

  type emb_out
    double precision, dimension(nye,nxe) ::   &
        tt, pp, snow, qq, rhum, &
        S, co2, ap, apS, ABT, d2T, Lr, Ls, Ldh
  end type
  
  double precision, dimension (nys,nxs) :: veg_pdds
  
contains 
  
  ! ##### REMBO OUTPUT FUNCTIONS ######
  
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  e m b _ n c
  ! Author   :  Alex Robinson
  ! Purpose  :  Initialize netcdf output for rembo
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine emb_nc(fnm,units,t0,lat,lon,xx,yy,mask,zs,zsp,dzsp,pkap,ptau)
    
    implicit none
    
    double precision, dimension(:,:) :: lat, lon, xx, yy
    double precision, dimension(:,:),optional :: mask, zs, zsp,dzsp,pkap,ptau
    character (len=*) :: fnm, units
    double precision :: t0
    
    call nc_create(fnm)

    ! Create the netcdf file and the dimension variables
    call nc_write_dim(fnm,"x",x=x0*1e-3,dx=dxe*1e-3,nx=nxe)
    call nc_write_dim(fnm,"y",x=y0*1e-3,dx=dxe*1e-3,nx=nye)
    call nc_write_dim(fnm,"time",x=[t0],units=units,unlimited=.TRUE.)
    
    ! Add grid variables 
    call nc_write_t(fnm,"lat",lat,dim1="x",dim2="y")
    call nc_write_t(fnm,"lon",lon,dim1="x",dim2="y")
    call nc_write_t(fnm,"xx",  xx,dim1="x",dim2="y")
    call nc_write_t(fnm,"yy",  yy,dim1="x",dim2="y")
    
    if (present(mask)) call nc_write_t(fnm,"mask",   mask,dim1="x",dim2="y")
    if (present(zs))   call nc_write_t(fnm,"zs",       zs,dim1="x",dim2="y")
    if (present(zsp))  call nc_write_t(fnm,"zsp",     zsp,dim1="x",dim2="y")
    if (present(dzsp)) call nc_write_t(fnm,"dzsp",   dzsp,dim1="x",dim2="y")
    if (present(pkap)) call nc_write_t(fnm,"pkappa", pkap,dim1="x",dim2="y")
    if (present(ptau)) call nc_write_t(fnm,"ptau",   ptau,dim1="x",dim2="y")
    
    return
  
  end subroutine emb_nc
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  e m b _ n c _ w r i t e
  ! Author   :  Alex Robinson
  ! Purpose  :  Output timestep of netcdf output for rembo
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine emb_nc_write(fnm,vars,ndat)
    
    implicit none
    
    type(emb_out) :: vars
    
    character (len=*) :: fnm
    integer :: ndat
    
    call nc_write(fnm,"time",ndat,dim1="time",start=[ndat],count=[1])   ! ka => years
    
    call nc_write(fnm,"tt",       vars%tt,  dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])
    call nc_write(fnm,"pp",       vars%pp,  dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])
    call nc_write(fnm,"snow",     vars%snow,dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])
    call nc_write(fnm,"qq",       vars%qq,  dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])
    call nc_write(fnm,"rhum",     vars%rhum,dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])
    call nc_write(fnm,"S",        vars%S,   dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])
    call nc_write(fnm,"co2",      vars%co2, dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])
    call nc_write(fnm,"ap",       vars%ap,  dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])
    call nc_write(fnm,"apS",      vars%apS, dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])
    call nc_write(fnm,"ABT",      vars%ABT, dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])
    call nc_write(fnm,"d2T",      vars%d2T, dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])
    call nc_write(fnm,"Lr",       vars%Lr,  dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])
    call nc_write(fnm,"Ls",       vars%Ls,  dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])
    call nc_write(fnm,"Ldh",      vars%Ldh, dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])
    call nc_write(fnm,"snow",     vars%snow,dim1="x",dim2="y",dim3="time",start=[1,1,ndat],count=[nxe,nye,1])

    return
    
  end subroutine emb_nc_write
   
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  r e m b o _ n c
  ! Author   :  Alex Robinson
  ! Purpose  :  Initialize netcdf output for rembo
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine rembo_nc(fnm,units,t0,lat,lon,xx,yy,mask,mask_hydro,zs)
    
    implicit none
    
    double precision, dimension(:,:) :: lat, lon, xx, yy
    double precision, dimension(:,:),optional :: mask, mask_hydro, zs
    character (len=*) :: fnm, units
    double precision :: t0
    
    call nc_create(fnm)

    ! Create the netcdf file and the dimension variables
    call nc_write_dim(fnm,"x",    x=x0*1e-3,dx=dxs*1e-3,nx=nxs)
    call nc_write_dim(fnm,"y",    x=y0*1e-3,dx=dxs*1e-3,nx=nys)
    call nc_write_dim(fnm,"month",x=1.d0,dx=1.d0,nx=12)
    call nc_write_dim(fnm,"time", x=[t0],units=units,unlimited=.TRUE.)
    
    ! Add grid variables 
    call nc_write_t(fnm,"lat",lat,dim1="x",dim2="y")
    call nc_write_t(fnm,"lon",lon,dim1="x",dim2="y")
    call nc_write_t(fnm,"xx", xx, dim1="x",dim2="y")
    call nc_write_t(fnm,"yy", yy, dim1="x",dim2="y")
    
    if (present(mask))       call nc_write_t(fnm,"mask",mask,dim1="x",dim2="y")
    if (present(mask_hydro)) call nc_write_t(fnm,"mask_hydro",mask_hydro,dim1="x",dim2="y")
    if (present(zs))         call nc_write_t(fnm,"zs",  zs,dim1="x",dim2="y")
    
    return
  
  end subroutine rembo_nc
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  r e m b o _ n c _ w r i t e
  ! Author   :  Alex Robinson
  ! Purpose  :  Output timestep of netcdf output for rembo
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine rembo_nc_write(fnm,vars,ndat,time,mask,mask_hydro,zs,tjja,tdjf,ttp,rf, &
                            refrozen,pkappa,qqfac,pdds,write_mon)
    
    implicit none
    
    type(rembo_type) :: vars
    double precision, dimension(:,:), optional :: mask, mask_hydro, zs, tjja, tdjf, ttp
    double precision, dimension(:,:), optional :: rf, refrozen,pkappa,qqfac, pdds
    logical, optional :: write_mon

    character (len=*) :: fnm
    integer :: ndat
    double precision :: time
    character (len=256), allocatable :: dims(:)

    call nc_write(fnm,"time",time,dim1="time",start=[ndat],count=[1])   ! years
    
    ! Add variables
    if (present(mask))       call nc_write_t(fnm,"mask",      mask,      dim1="x",dim2="y",dim3="time", &
                                            start=[1,1,1],count=[nxs,nys,1])
    if (present(mask_hydro)) call nc_write_t(fnm,"mask_hydro",mask_hydro,dim1="x",dim2="y")
    if (present(zs))         call nc_write_t(fnm,"zs",        zs  ,      dim1="x",dim2="y",dim3="time", &
                                            start=[1,1,1],count=[nxs,nys,1])
    
    if (allocated(dims)) deallocate(dims); allocate(dims(3))
    dims(1) = "x"; dims(2) = "y"; dims(3) = "time"

    call nc_write_t(fnm,"tt",       vars%tt,dims=dims,start=[1,1,ndat],count=[nxs,nys,1])
    call nc_write_t(fnm,"tte",      vars%tte,dims=dims,start=[1,1,ndat],count=[nxs,nys,1]) 
    call nc_write_t(fnm,"pp",       vars%pp,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm/a")
    call nc_write_t(fnm,"snow",     vars%snow,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm/a")
    call nc_write_t(fnm,"snowh",    vars%snowh,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm")
    call nc_write_t(fnm,"ice",      vars%ice,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm/a")
    call nc_write_t(fnm,"runoff",   vars%runoff,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm/a")
    call nc_write_t(fnm,"melt",     vars%melt,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm/a")
    call nc_write_t(fnm,"smb",      vars%pp-vars%runoff,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm/a")
    call nc_write_t(fnm,"smbi",     vars%ice-vars%melted_ice,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm/a")
    call nc_write_t(fnm,"S",        vars%S,  dims=dims,start=[1,1,ndat],count=[nxs,nys,1])
    if (melt_choice .le. 0) call nc_write_t(fnm,"S0", vars%S0,  dims=dims,start=[1,1,ndat],count=[nxs,nys,1])
    call nc_write_t(fnm,"as",       vars%as, dims=dims,start=[1,1,ndat],count=[nxs,nys,1])
    call nc_write_t(fnm,"ap",       vars%ap, dims=dims,start=[1,1,ndat],count=[nxs,nys,1])
    call nc_write_t(fnm,"ap0",      vars%ap0,dims=dims,start=[1,1,ndat],count=[nxs,nys,1])
    call nc_write_t(fnm,"rf",       vars%rf,dims=dims,start=[1,1,ndat],count=[nxs,nys,1])
    call nc_write_t(fnm,"refrozen", vars%refrozen,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm/a")
    call nc_write_t(fnm,"melt_insol", vars%melt_insol,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm/a")
    call nc_write_t(fnm,"melt_S0",    vars%melt_S0,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm/a")

    ! If desired, write mean fields too
    if (present(tjja))     call nc_write_t(fnm,"tjja",tjja,dims=dims,start=[1,1,ndat],count=[nxs,nys,1]) 
    if (present(tdjf))     call nc_write_t(fnm,"tdjf",tdjf,dims=dims,start=[1,1,ndat],count=[nxs,nys,1]) 
    if (present(ttp))      call nc_write_t(fnm,"ttp",ttp,dims=dims,start=[1,1,ndat],count=[nxs,nys,1]) 
    if (present(rf))       call nc_write_t(fnm,"jh_rf",rf,dims=dims,start=[1,1,ndat],count=[nxs,nys,1]) 
    if (present(refrozen)) call nc_write_t(fnm,"jh_refrozen",refrozen,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm/a") 
    if (present(pkappa))   call nc_write_t(fnm,"pkappa",pkappa,dims=dims,start=[1,1,ndat],count=[nxs,nys,1])
    if (present(qqfac))    call nc_write_t(fnm,"qqfac",qqfac,dims=dims,start=[1,1,ndat],count=[nxs,nys,1])
    if (present(pdds))     call nc_write_t(fnm,"pdds",pdds,dims=dims,start=[1,1,ndat],count=[nxs,nys,1])

    if (present(write_mon)) then 
      if (write_mon) then
        ! Set up dims for 3d vars, x,y,time
        if (allocated(dims)) deallocate(dims); allocate(dims(4))
        dims(1) = "x"; dims(2) = "y"; dims(3) = "month"; dims(4) = "time"

        call nc_write_t(fnm,"m_S",       savedm%S,dims=dims,start=[1,1,1,ndat],count=[nxs,nys,nm,1],units="W m-2")
        call nc_write_t(fnm,"m_S0",      savedm%S0,dims=dims,start=[1,1,1,ndat],count=[nxs,nys,nm,1],units="W m-2")
        call nc_write_t(fnm,"m_tt",      savedm%tt,dims=dims,start=[1,1,1,ndat],count=[nxs,nys,nm,1],units="degrees celcius")
        call nc_write_t(fnm,"m_melt",    savedm%melt/30,dims=dims,start=[1,1,1,ndat],count=[nxs,nys,nm,1],units="mm d-1")
        call nc_write_t(fnm,"m_as",      savedm%as,dims=dims,start=[1,1,1,ndat],count=[nxs,nys,nm,1],units="1")
        call nc_write_t(fnm,"m_pdd_corr",savedm%pdd_corr,dims=dims,start=[1,1,1,ndat],count=[nxs,nys,nm,1],units="mm d-1 / W m-2")
      end if
    end if 

    return
    
  end subroutine rembo_nc_write
  
  subroutine rembo_nc_write_small(fnm,ann,mon,ndat,time,mask,zs,tjja,tdjf,ttp,pdds)
    
    implicit none
    
    type(rembo_type) :: ann, mon(:)
    double precision, dimension(:,:) :: mask, zs, tjja, tdjf, ttp, pdds

    character (len=*) :: fnm
    integer :: ndat
    double precision :: time
    character (len=256), allocatable :: dims(:)

    if (allocated(dims)) deallocate(dims); allocate(dims(3))
    dims(1) = "x"; dims(2) = "y"; dims(3) = "time"

    call nc_write(fnm,"time",time,dim1="time",start=[ndat],count=[1])   ! years
    
    ! Add variables
    call nc_write_t(fnm,"mask",     mask,dims=dims,start=[1,1,ndat],count=[nxs,nys,1])
    call nc_write_t(fnm,"zs",       zs  ,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="m")
    
    call nc_write_t(fnm,"tann",     ann%tt,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="degrees C")
    call nc_write_t(fnm,"tjan",     mon(1)%tt,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="degrees C")
    call nc_write_t(fnm,"tjul",     mon(7)%tt,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="degrees C")
    call nc_write_t(fnm,"pp",       ann%pp,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm/a")
    call nc_write_t(fnm,"snow",     ann%snow,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm/a")
    call nc_write_t(fnm,"smb",      ann%pp-ann%runoff,dims=dims,start=[1,1,ndat],count=[nxs,nys,1],units="mm/a")

    call nc_write_t(fnm,"pdds",pdds,dims=dims,start=[1,1,ndat],count=[nxs,nys,1])
    call nc_write_t(fnm,"tjja",tjja,dims=dims,start=[1,1,ndat],count=[nxs,nys,1]) 
    call nc_write_t(fnm,"tdjf",tdjf,dims=dims,start=[1,1,ndat],count=[nxs,nys,1]) 
    call nc_write_t(fnm,"ttp", ttp,dims=dims,start=[1,1,ndat],count=[nxs,nys,1]) 

    return
    
  end subroutine rembo_nc_write_small
    
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  r e m b o _ r e s t a r t
  ! Author   :  Alex Robinson
  ! Purpose  :  Make a restart file for rembo
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine rembo_restart(yearnow0,mask,zs,in,pdds,file)
    
    implicit none
    
    integer :: n
    double precision :: yearnow, yearnow0
    double precision, dimension(:,:) :: mask, zs,pdds
    type(rembo_type) :: in
    character (len=256) :: fnm, c_yr
    character (len=*), optional :: file
    
    yearnow = dabs(yearnow0*1d-3) ! Get abs value of current year in ka
    
    if ( yearnow .lt. 1.d0 .and. yearnow .ne. 0.d0 ) then
      write(c_yr,"(a4,i1)") "000.",int(yearnow*10.d0)
    else if ( yearnow .lt. 10.d0 ) then
      write(c_yr,"(a2,i1)") "00",int(yearnow)
    else if ( yearnow .lt. 100.d0 ) then
      write(c_yr,"(a1,i2)") "0",int(yearnow)
    else if ( yearnow .lt. 1000.d0 ) then
      write(c_yr,"(i3)") int(yearnow)
    else
      write(c_yr,"(i4)") int(yearnow)
    end if
    
    ! Add on bp designation if needed
    if ( yearnow0 .lt. 0.d0 ) c_yr = "bp"//trim(c_yr)
    
    fnm = trim(outfldr)//"restart."//trim(c_yr)//".nc"
    
    ! Overwrite local naming convention if filename is given as argument
    if ( present(file) ) fnm = trim(file)
    
    ! Initialize the file (if running in stand-alone mode)
    if ( clim_coupled .eq. -1 ) then
      call rembo_nc(fnm,"years",yearnow0, &
                    sico_grid%lat*todegs,sico_grid%lon*todegs,sico_grid%x,sico_grid%y, &
                    mask=mask,zs=zs)
    end if
    
    ! Write the fields
    call nc_write_t(fnm,"tt",   in%tt,   dim1="x",dim2="y",long_name="surface temperature (last day of year)")
    call nc_write_t(fnm,"snowh",in%snowh,dim1="x",dim2="y",long_name="snow height (last day of year)")
    call nc_write_t(fnm,"dh",   in%dh,   dim1="x",dim2="y",long_name="net melt (last day of year)",units="mm")
    call nc_write_t(fnm,"pdds", pdds,    dim1="x",dim2="y",long_name="annual PDDs (for vegetation type)",units="PDD")
        
    if ( ap_fixed .eq. 1 ) then
      call nc_write_t(fnm,"ap",   in%ap0,  dim1="x",dim2="y",long_name="planetary albedo (last day of year)")
    else
      call nc_write_t(fnm,"ap",   in%ap,   dim1="x",dim2="y",long_name="planetary albedo (last day of year)")
    end if

    return
  
  end subroutine rembo_restart

  
  ! ##### REMBO MATH FUNCTIONS ######
  ! Note: overloading = makes it much slower!! Why??? Instead, made
  !       function to perform individual calculations
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  r e m b o _ a v e
  ! Author   :  Alex Robinson
  ! Purpose  :  Average (or sum) the rembo fields as needed
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine rembo_ave(d,ave,clim,smb)
    
    implicit none
    
    type(rembo_type) :: d(:), ave
    
    integer :: k, n
    double precision :: div
    
    logical, optional :: clim, smb
    logical :: calc_clim, calc_smb
    
    n = size(d)
    div = dble(n)
    
    calc_clim = .FALSE.
    if ( present(clim) ) calc_clim = clim
    
    calc_smb  = .TRUE.
    if ( present(smb) ) calc_smb = smb
    
    if ( calc_clim ) then  ! Calculate averages for climate variables
      
      ! Set all rembo values to zero
      ave%tt      = 0.d0
      ave%tte     = 0.d0
      ave%pp      = 0.d0
      ave%snow    = 0.d0
      
      ! Loop over the time indices to sum up
      do k = 1, n
        ave%tt      = ave%tt         + d(k)%tt
        ave%tte     = ave%tte        + d(k)%tte
        ave%pp      = ave%pp         + d(k)%pp
        ave%snow    = ave%snow       + d(k)%snow
      end do
      
      ! Get the average where necessary
      ave%tt      = ave%tt      / div
      ave%tte     = ave%tte     / div
      
    end if
    
    
    if ( calc_smb ) then  ! Now calculate averages for smb variables
    
      ! Set all rembo values to zero
      ave%snowh   = 0.d0;    ave%ice       = 0.d0
      ave%runoff  = 0.d0;    ave%melt      = 0.d0
      ave%S       = 0.d0;    ave%dh        = 0.d0
      ave%rf      = 0.d0;    ave%as        = 0.d0
      ave%ap      = 0.d0;    ave%refrozen  = 0.d0
      ave%runoff_snow  = 0.d0
      ave%runoff_rain  = 0.d0
      ave%melted_ice   = 0.d0
      ave%melt_insol   = 0.d0
      ave%melt_S0      = 0.d0
      ave%pdd_corr     = 0.d0
      ave%S0           = 0.d0 

      ! Loop over the time indices to sum up
      do k = 1, n

        ave%snowh   = ave%snowh      + d(k)%snowh
        ave%ice     = ave%ice        + d(k)%ice
        ave%runoff  = ave%runoff     + d(k)%runoff
        ave%melt    = ave%melt       + d(k)%melt
        ave%S       = ave%S          + d(k)%S
        ave%S0      = ave%S0         + d(k)%S0
        ave%dh      = ave%dh         + d(k)%dh
        ave%rf      = ave%rf         + d(k)%rf
        ave%as      = ave%as         + d(k)%as
        ave%ap      = ave%ap         + d(k)%ap
        ave%refrozen  = ave%refrozen + d(k)%refrozen
        
        ave%runoff_snow  = ave%runoff_snow + d(k)%runoff_snow 
        ave%runoff_rain  = ave%runoff_rain + d(k)%runoff_rain
        ave%melted_ice   = ave%melted_ice  + d(k)%melted_ice
        ave%melt_insol   = ave%melt_insol  + d(k)%melt_insol
        ave%melt_S0      = ave%melt_S0     + d(k)%melt_S0
        ave%pdd_corr     = ave%pdd_corr    + d(k)%pdd_corr
      end do

      ! Get the average where necessary
      ave%snowh   = ave%snowh   / div
      ave%S       = ave%S       / div
      ave%S0      = ave%S0      / div
      ave%rf      = ave%rf      / div
      ave%as      = ave%as      / div
      ave%ap      = ave%ap      / div
      
      ave%pdd_corr = ave%pdd_corr / div

      ! Other fields are sums: pp, snow, ice, melt, melt, dh, frozen
      
      ! Note: including S0 in averaging procedure here is a repetetive, redundant
      ! calculation, however it was the cleanest incorporation of the calculation.
      ! The cost is minimal, but should be improved in future versions. [ajr, 2013-03-08]

    end if
    
      
    return
    
  end subroutine rembo_ave
  
  subroutine rembo_array(var,d,varname)

    type(rembo_type) :: d(:)
    character (len=*) :: varname

    double precision, dimension(:,:,:) :: var

    integer:: im, nm

    nm = size(d)

!     if (allocated(rembo_array)) deallocate(rembo_array)
!     allocate(rembo_array(nm,ny,nx))

    if (trim(varname)=="pdd_corr") then
      do im = 1, nm 
        var(im,:,:) = d(im)%pdd_corr 
      end do 
    else if (trim(varname)=="S") then 
      do im = 1, nm 
        var(im,:,:) = d(im)%S
      end do
    else if (trim(varname)=="S0") then 
      do im = 1, nm 
        var(im,:,:) = d(im)%S0
      end do
    else if (trim(varname)=="tt") then 
      do im = 1, nm 
        var(im,:,:) = d(im)%tt
      end do
    else if (trim(varname)=="melt") then 
      do im = 1, nm 
        var(im,:,:) = d(im)%melt
      end do
    else if (trim(varname)=="as") then 
      do im = 1, nm 
        var(im,:,:) = d(im)%as
      end do

    else
      write(*,*) "rembo:: rembo_array error: variable with this name not allowed: "//trim(varname) 
      stop
    end if 

    return

  end subroutine rembo_array

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : m a v e 
  ! Author     : Alex Robinson, 27. Oct 2008
  ! Purpose    : make an area sum or average filtered by a mask
  !              return value as an average of the input units
  !              or in [1e12 Gt/a]
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function mave(ma,dx,var_in,ave,area)
    
    implicit none
    
    integer:: ny, nx
    real (8), dimension(:,:) :: ma, var_in
    real (8), dimension(:,:), allocatable :: var
    real (8) :: dx, mave
    
    character (len=*), optional :: ave
    real (8), dimension(:,:), optional :: area
    
    ny = size(var_in,1)
    nx = size(var_in,2)
    
    allocate(var(ny,nx))
    
    ! Reset var to zero
    var = 0.d0
    
    ! Use variable values only where mask is equal to one
    where (ma .eq. 1.d0)  var = var_in
    
    if (present(ave)) then
      if (ave .eq. "ave") then
        mave = sum(var) / sum(ma)
      else
        write(stdout,*) "Incorrect choice for mave."
        write(stdout,*) "Should be 'ave' or nothing. Try again!"
        call embout(1)
      end if
    else
      if (present(area)) then
        mave = sum( var*area )  / 1.d12
      else
        mave = sum(var) * dx**2 / 1.d12
      end if
    end if
    
    deallocate(var)
    
    return
    
  end function mave  
  
  subroutine tune_qqfac(pp,pp0,mask,lons,qqfac)
  
    implicit none
    
    integer:: ny, nx
    double precision, parameter :: pp_max = 3000.d0
    double precision, parameter :: dp_max =  300d0  ! mm
    double precision, parameter :: qqfac_min = 1.0d0
    double precision, parameter :: qqfac_max = 3.0d0
    double precision, parameter :: dfac = 0.1d0 * (qqfac_max-qqfac_min)
    
    real (8), dimension(:,:) :: pp, pp0, mask, lons, qqfac
    real (8), dimension(:,:), allocatable :: qqfac_new, pp1, pp00, resid
    real (8) :: dx, mave
    
    real (8) :: klon, llmin, llmax

    llmin = minval(lons); llmax = maxval(lons)
    klon = -1.d0
    
    ny = size(pp,1)
    nx = size(pp,2)
    
    allocate(qqfac_new(ny,nx),pp00(ny,nx),resid(ny,nx),pp1(ny,nx))
    
    ! If it's the first timestep, initialize kappa field
    if ( nstep .ge. 0 .or. equili .ne. 0 .or. tuning .ne. 1 ) then
    
      qqfac = 1.d0
      
!       qqfac = ( 1.d0 + klon*(2.d0*(lons-llmin)/(llmax-llmin) - 1.d0) )
!       where ( qqfac .lt. 1.d0 ) qqfac = 1.d0
      
      
    else  ! Adjust kappa locally
    
      ! Store old kappa as starting point for new kappa
      qqfac_new = qqfac
      
      ! Get new array of precip data (limit to 1000mm)
      pp00 = pp0
      where( pp00 .gt. pp_max ) pp00 = pp_max
      
      ! Store precip and modify it
      pp1 = pp
      
      ! Determine the residuals with precip data and normalize
      resid = (pp1 - pp00) / dp_max
      where ( resid .gt.  1.d0 ) resid =  1.d0
      where ( resid .lt. -1.d0 ) resid = -1.d0
      
      ! Artificially reduce positive residual
      ! (since I can't improve this much...)
      where ( resid .gt. 0.d0 ) resid = resid*0.2d0
      
      ! Make sure only ocean points recieve more moisture
      where ( mask .ne. 2.d0 ) resid = 0.d0
      
      ! Adjust the kappa field by desired amount
      ! (scaled by residual at that point)
      qqfac_new = qqfac - resid*dfac

      ! Scale the new kappa field to make sure it's within desired bounds
      where ( qqfac_new .gt. qqfac_max ) qqfac_new = qqfac_max
      where ( qqfac_new .lt. qqfac_min ) qqfac_new = qqfac_min
      
      !! Ignore previous, just calculate adjustment based on 
      !! magnitude of pp0 precip
      resid = pp0
      where ( pp0 .lt. 1500.d0 ) resid = 1.d0
      where ( pp0 .ge. 1500.d0 ) resid = resid/1500.d0*2.d0
      qqfac_new = resid
      
      ! Save new kappa to original
      qqfac = qqfac_new
      
    end if
    
    return
    
  end subroutine tune_qqfac
  
  subroutine tune_pkappa(pp,pp0,mask,kappa)
  
    implicit none
    
    integer:: ny, nx
    double precision, parameter :: pp_max = 800.d0
    double precision, parameter :: dp_max = 300d0  ! mm
    double precision, parameter :: kap_min =  10d4
    double precision, parameter :: kap_max =  60d4
    double precision, parameter :: dkap = 0.1d0 * (kap_max-kap_min)
    
    real (8), dimension(:,:) :: pp, pp0, mask, kappa
    real (8), dimension(:,:), allocatable :: kappa_new, pp00, resid
    real (8) :: dx, mave
    
    ny = size(pp,1)
    nx = size(pp,2)
    
    allocate(kappa_new(ny,nx),pp00(ny,nx),resid(ny,nx))
    
    ! If it's the first timestep, initialize kappa field
    if ( nstep .eq. 0 .or. equili .ne. 0 .or. tuning .ne. 1 ) then
    
      kappa = pkappa
      
    else  ! Adjust kappa locally
    
      ! Store old kappa as starting point for new kappa
      kappa_new = kappa
      
      ! Get new array of precip data (limit to 1000mm)
      pp00 = pp0
      where( pp00 .gt. pp_max ) pp00 = pp_max
      
      ! Determine the residuals with precip data and normalize
      resid = (pp - pp00) / dp_max
      where ( resid .gt.  1.d0 ) resid = 1.d0
      where ( resid .lt. -1.d0 ) resid = -1.d0
      
      ! Artificially reduce negative residual
      ! (since I can't improve this much...)
      !where ( resid .lt. 0.d0 ) resid = resid*0.2d0
      
      ! Make sure only ocean points are less weighted
      !where ( mask .eq. 2.d0 ) resid = resid*0.1d0
      
      ! Adjust the kappa field by desired amount
      ! (scaled by residual at that point)
      kappa_new = kappa - resid*dkap
      
      ! Scale the new kappa field to make sure it's within desired bounds
      where ( kappa_new .gt. kap_max ) kappa_new = kap_max
      where ( kappa_new .lt. kap_min ) kappa_new = kap_min
      
      ! Set ocean points equal to original kappa
      !where ( mask .eq. 2.d0 ) kappa_new = pkappa
      
      ! Save new kappa to original
      kappa = kappa_new
      
    end if
    
    return
    
  end subroutine tune_pkappa
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c l i m c h e c k e r
  ! Author     :  Alex Robinson
  ! Purpose    :  Output means/sums of variables for a given region
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine climchecker_new(fnm,d,m,a,mask,zs,lats,lons,T_warming,T_anomaly,yearnow)
    
    implicit none
    
    integer,  parameter :: nx = nxs, ny = nys
    real (8), parameter :: dx = dxs
    
    character (len=*) :: fnm
    type(rembo_type)  :: d(:), m(:), a
    real (8), dimension(ny,nx) :: mask, zs, lats, lons
    real(8), intent(IN) :: T_warming        ! Summer anomaly value
    real(8), intent(IN) :: T_anomaly(:)     ! Daily anomaly values

    real (8) :: yearnow
    
    real (8), parameter :: ela_ebs = 100.d0   ! mm/a
    real (8), parameter :: area_conv = 1e-6*1e-6  ! m2=>km2=>1e6km2
    real (8), parameter :: vol_conv  = 1e-9*1e-6  ! m3=>km3=>1e6km3
    real (8), parameter :: Gt_conv   = 1e-12      ! m3 w.e. => Gt
    real (8), parameter, dimension(3) :: melt_min = (/ 7.75d0,8.5d0,9.25d0 /)   ! mm/d
    
    real (8), dimension(ny,nx) :: ice, land, ima, inmb, melt
    
    integer :: nm13 = 13 
    real (8) :: dT, aco2, dT_all(13), aco2_all(13), dT_bnd(13)
    real (8) :: cma775(13), cma850(13), cma925(13)
    real (8), dimension(nk) :: dma775,dma850,dma925

    integer  :: i, j, kd, km, k1, k2, k, n_tot, nrec
    real (8) :: area, time, tmp(1)
    integer  :: j65_0,j65_1
    real (8) :: frac65, S65, lat65
    
    character (len=256), allocatable :: dims(:)
    
    type output1
      double precision, dimension(13) :: &
              tt, pp, tte, as, ap, snow, snowh, ice, runoff,  &
              melt, smb, smbi, pmb, nmb, ma, refrozen, &
              S, S65, melt_insol, melt_S0, pdd_corr, &
              ela, aar  
    end type
  
    type(output1) :: gis

    if ( sum(mask) .eq. 0.d0 ) then
      write(stdout,*) "**Currently no ice available, using fixed ice sheet size"
      write(stdout,*) "  equal to present-day in climcheck to compare accum/abl patterns!"
      
      ! Load present-day ice mask
      call nc_read_t(trim(topo_file),"H_ice",mask)
      
      where(mask .gt. 0.d0) 
        mask = 1.d0 
      elsewhere 
        mask = 0.d0 
      end where 

    end if
    
    ! Determine index of 65 deg North Latitude
    lat65 = 65.d0*torads
    i = 1
    if ( lat65 .le. maxval(lats(:,i)) .and. lat65 .gt. minval(lats(:,i)) ) then

      do j = 2, ny
        if (lats(j,i) .ge. lat65) then
          j65_0 = j-1
          j65_1 = j
          exit
        end if
      end do
    
      ! Determine interpolation fraction to get exactly 65 deg North
      frac65 = (1.d0 - (lats(j65_1,i) - lat65 ) / &
                     (lats(j65_1,i) - lats(j65_0,i) )  )
    
    else
      j65_0 = 1
      j65_1 = 2
      frac65 = 1
    end if

    ! Determine areas of current mask of interest
    ! (assuming ice is only land-based)
    area  = sum(mask) *dx**2 *area_conv   
    n_tot = int(sum(mask))
    
    km = 13

    ! Determine daily melt extent
    dma775 = 0.d0; dma850 = 0.d0; dma925 = 0.d0
    do k = 1, nk
      
      melt = 0.d0
      where ( mask .eq. 1.d0 .and. d(k)%melt .gt. melt_min(1) ) melt = 1.d0
      dma775(k) = sum(melt)*dx**2 * area_conv   ! scale from m2=>km2 and then from km2 => 1e6km2
      
      melt = 0.d0
      where ( mask .eq. 1.d0 .and. d(k)%melt .gt. melt_min(2) ) melt = 1.d0
      dma850(k) = sum(melt)*dx**2 * area_conv   ! scale from m2=>km2 and then from km2 => 1e6km2
      
      melt = 0.d0
      where ( mask .eq. 1.d0 .and. d(k)%melt .gt. melt_min(3) ) melt = 1.d0
      dma925(k) = sum(melt)*dx**2 *area_conv   ! scale from m2=>km2 and then from km2 => 1e6km2
    
    end do
    
    cma775(km) = sum(dma775)   
    cma850(km) = sum(dma850) 
    cma925(km) = sum(dma925) 
    
    ! Get insolation at 65 deg North at summer solstice (max value)
    S65 = frac65*maxval(d%S(j65_1,1)) + (1.d0-frac65)*maxval(d%S(j65_0,1))
    
    ! ma (annual melt area - not cumulative, in the sense that melt is only counted once per place)
    ! nmb (negative annual surface mass balance)
    ! pmb (positive annual surface mass balance)
    ima = 0.d0; inmb = 0.d0
    where (mask .eq. 1.d0 .and. a%melt .gt. melt_min(2) )   ima  = 1.d0
    where (mask .eq. 1.d0 .and. a%pp - a%runoff .lt. 0.d0 ) inmb = 1.d0
    
    gis%ma(km)  = sum(ima) *dx**2 *area_conv
    gis%nmb(km) = sum(inmb)*dx**2 *area_conv
    gis%pmb(km) = area - gis%nmb(km)
    
    if (gis%pmb(km) + gis%nmb(km) .ne. 0.0) then 
      gis%aar(km) = gis%pmb(km) / (gis%pmb(km) + gis%nmb(km))
    else
      gis%aar(km) = 0.0 
    end if 

    gis%tt(km)       = mave(mask,dx,a%tt,"ave")
    gis%pp(km)       = mave(mask,dx,a%pp)
    gis%tte(km)      = mave(mask,dx,a%tte,"ave")
    gis%S(km)        = mave(mask,dx,a%S,"ave")
    gis%S65(km)      = S65
    gis%as(km)       = mave(mask,dx,a%as,"ave")
    gis%ap(km)       = mave(mask,dx,a%ap,"ave")
    gis%snow(km)     = mave(mask,dx,a%snow)
    gis%snowh(km)    = mave(mask,dx,a%snow,"ave")
    gis%ice(km)      = mave(mask,dx,a%ice)
    gis%runoff(km)   = mave(mask,dx,a%runoff)
    gis%melt(km)     = mave(mask,dx,a%melt)
    gis%smb(km)      = mave(mask,dx,a%pp-a%runoff)
    gis%smbi(km)     = mave(mask,dx,a%ice-a%melted_ice)
    gis%refrozen(km) = mave(mask,dx,a%refrozen)
    
    gis%melt_insol(km)  = mave(mask,dx,a%melt_insol)
    gis%melt_S0(km)     = mave(mask,dx,a%melt_S0)
    gis%pdd_corr(km)    = mave(mask,dx,a%pdd_corr)

    ! Additionally get Equilibrium Line Altitude (ELA) -----
    call find_ela(gis%ela(km),zs,a%pp-a%runoff,mask=mask.eq.1.d0)

    do km = 1, nm
      
      ! Get indices of days in current month
      k1 = int( (km-1)*day_month) + 1 
      k2 = k1 + int(day_month) - 1
      
      ! Determine monthly melt extent from daily contribution
      cma775(km) = sum(dma775(k1:k2))
      cma850(km) = sum(dma850(k1:k2))
      cma925(km) = sum(dma925(k1:k2))
      
      ! Get insolation at 65 deg North for this time
      S65 = frac65*m(km)%S(j65_1,1) + (1.d0-frac65)*m(km)%S(j65_0,1)
   
      ima = 0.d0; inmb = 0.d0
      where (mask .eq. 1.d0 .and. m(km)%melt .gt. melt_min(2) )       ima  = 1.d0
      where (mask .eq. 1.d0 .and. m(km)%pp - m(km)%runoff .lt. 0.d0 ) inmb = 1.d0
      
      gis%ma(km)  = sum(ima) *dx**2 *area_conv
      gis%nmb(km) = sum(inmb)*dx**2 *area_conv
      gis%pmb(km) = area - gis%nmb(km) 
      
      if (gis%pmb(km) + gis%nmb(km) .ne. 0.0) then 
        gis%aar(km) = gis%pmb(km) / (gis%pmb(km) + gis%nmb(km))
      else
        gis%aar(km) = 0.0 
      end if 
    
      ! Store monthly means in global output arrays
      gis%tt(km)       = mave(mask,dx,m(km)%tt,"ave")
      gis%pp(km)       = mave(mask,dx,m(km)%pp)          ! [1d12 Gt/month]
      gis%tte(km)      = mave(mask,dx,m(km)%tte,"ave")
      gis%S(km)        = mave(mask,dx,m(km)%S,"ave")
      gis%S65(km)      = S65
      gis%as(km)       = mave(mask,dx,m(km)%as,"ave")
      gis%ap(km)       = mave(mask,dx,m(km)%ap,"ave")
      gis%snow(km)     = mave(mask,dx,m(km)%snow)
      gis%snowh(km)    = mave(mask,dx,m(km)%snow,"ave")
      gis%ice(km)      = mave(mask,dx,m(km)%ice)
      gis%runoff(km)   = mave(mask,dx,m(km)%runoff)
      gis%melt(km)     = mave(mask,dx,m(km)%melt)
      gis%smb(km)      = mave(mask,dx,m(km)%pp-m(km)%runoff)
      gis%smbi(km)     = mave(mask,dx,m(km)%ice-m(km)%melted_ice)
      gis%refrozen(km) = mave(mask,dx,m(km)%refrozen)
      
      gis%melt_insol(km)  = mave(mask,dx,m(km)%melt_insol)
      gis%melt_S0(km)     = mave(mask,dx,m(km)%melt_S0)
      gis%pdd_corr(km)    = mave(mask,dx,m(km)%pdd_corr)

      ! Additionally get Equilibrium Line Altitude (ELA) -----
      call find_ela(gis%ela(km),zs,m(km)%pp-m(km)%runoff,mask=mask.eq.1.d0) 

    end do
    
    ! Initialize output file (netcdf)
    if (init_summary2 .eq. 0) then
      
      call nc_create(fnm)
      call nc_write_dim(trim(fnm),"x",    x=x0*1e-3,nx=1)
      call nc_write_dim(trim(fnm),"y",    x=y0*1e-3,nx=1)
      call nc_write_dim(trim(fnm),"time", x=[year0], units="years",unlimited=.TRUE. )
      call nc_write_dim(trim(fnm),"month",x=1.d0,dx=1.d0,nx=13,units="month",axis="Z")
      call nc_write_dim(trim(fnm),"day",  x=1.d0,dx=1.d0,nx=nk,units="day")
      call nc_write_dim(trim(fnm),"point",x=1.d0,dx=1.d0,nx=1,units="n/a")
      
      call nc_write_dim(trim(fnm),"nrec", x=0.d0,dx=1.d0,nx=1,units="n/a")

      init_summary2 = 1
      
    end if
    
    ! Determine the number of records in this file, add 1
    call nc_read(trim(fnm),"nrec",nrec)
    nrec = nrec + 1
    
    do km = 1, 13
      
      if (km .eq. 13) dT = T_warming
      if (km .lt. 13) dT = T_anomaly(km*int(day_month)-15)
      aco2 = co2(dT)
      
      dT_all(km)   = dT    
      aco2_all(km) = aco2
    end do
    
    ! Use the forcing_now object to get dT without albedo amplification
    dT_bnd(1:12) = forcing_now%dTbnd 
    dT_bnd(13)   = sum(forcing_now%dTbnd(6:8)) / 3d0
    
    ! Update time variables
    call nc_write(trim(fnm),"time",yearnow,dim1="time",start=[nrec],count=[1])

    call nc_write(trim(fnm),"nrec",nrec,   dim1="nrec",start=[1],count=[1])
    
    write(*,"(a,a2,a6,i5)") trim(fnm),":","nrec=",nrec
    
    ! Set up dims for 1d vars, dependent on time
    if (allocated(dims)) deallocate(dims); allocate(dims(2))
    dims(1) = "month"; dims(2) = "time"
    
    ! GIS output (just ice sheet)
    call nc_write(trim(fnm),"tt",        gis%tt,dims=dims,start=[1,nrec],count=[nm13,1])
    call nc_write(trim(fnm),"pp",        gis%pp,dims=dims,start=[1,nrec],count=[nm13,1],units="Gt/a")
    call nc_write(trim(fnm),"as",        gis%as,dims=dims,start=[1,nrec],count=[nm13,1])
    call nc_write(trim(fnm),"ap",        gis%ap,dims=dims,start=[1,nrec],count=[nm13,1])
    call nc_write(trim(fnm),"tte",       gis%tte,dims=dims,start=[1,nrec],count=[nm13,1])
    call nc_write(trim(fnm),"S65",       gis%S65,dims=dims,start=[1,nrec],count=[nm13,1])
    call nc_write(trim(fnm),"S",         gis%S,dims=dims,start=[1,nrec],count=[nm13,1])
    call nc_write(trim(fnm),"dT",        dT_all,dims=dims,start=[1,nrec],count=[nm13,1],units="degrees celcius")
    call nc_write(trim(fnm),"dTb",       dT_bnd,dims=dims,start=[1,nrec],count=[nm13,1],units="degrees celcius")
    call nc_write(trim(fnm),"aco2",      aco2_all,dims=dims,start=[1,nrec],count=[nm13,1],units="ppm")
    call nc_write(trim(fnm),"snowh",     gis%snowh,dims=dims,start=[1,nrec],count=[nm13,1],units="mm")
    call nc_write(trim(fnm),"pmb",       gis%pmb,dims=dims,start=[1,nrec],count=[nm13,1],long_name="pos. smb area")
    call nc_write(trim(fnm),"nmb",       gis%nmb,dims=dims,start=[1,nrec],count=[nm13,1],long_name="neg. smb area")
!     call nc_write(trim(fnm),"ma",        gis%ma,dims=dims,start=[1,nrec],count=[nm13,1],long_name="max. melt area")
!     call nc_write(trim(fnm),"cma775",    cma775,dims=dims,start=[1,nrec],count=[nm13,1],long_name="cum. melt area, min 7.75mm/d")
    call nc_write(trim(fnm),"cma850",    cma850,dims=dims,start=[1,nrec],count=[nm13,1],long_name="cum. melt area, min 8.50mm/d")
!     call nc_write(trim(fnm),"cma925",    cma925,dims=dims,start=[1,nrec],count=[nm13,1],long_name="cum. melt area, min 9.25mm/d")
    call nc_write(trim(fnm),"rf",        gis%refrozen,dims=dims,start=[1,nrec],count=[nm13,1],long_name="refrozen",units="Gt/a")
    call nc_write(trim(fnm),"ice",       gis%ice,dims=dims,start=[1,nrec],count=[nm13,1],units="Gt/a")
    call nc_write(trim(fnm),"runoff",    gis%runoff,dims=dims,start=[1,nrec],count=[nm13,1],units="Gt/a")
    call nc_write(trim(fnm),"snow",      gis%snow,dims=dims,start=[1,nrec],count=[nm13,1],units="Gt/a")
    call nc_write(trim(fnm),"melt",      gis%melt,dims=dims,start=[1,nrec],count=[nm13,1],units="Gt/a")
    call nc_write(trim(fnm),"smb",       gis%smb,dims=dims,start=[1,nrec],count=[nm13,1],units="Gt/a")
    call nc_write(trim(fnm),"smbi",      gis%smbi,dims=dims,start=[1,nrec],count=[nm13,1],units="Gt/a")

    call nc_write(trim(fnm),"melt_insol",      gis%melt_insol,dims=dims,start=[1,nrec],count=[nm13,1],units="Gt/a")
    call nc_write(trim(fnm),"melt_S0",         gis%melt_S0,   dims=dims,start=[1,nrec],count=[nm13,1],units="Gt/a")
    call nc_write(trim(fnm),"pdd_corr",        gis%pdd_corr,  dims=dims,start=[1,nrec],count=[nm13,1],units="Gt/a")

    call nc_write(trim(fnm),"dTdt",      transT%dTdt,dim1="time",start=[nrec],count=[1],units="deg/1e6yr")
    call nc_write(trim(fnm),"mask_area", area,       dim1="time",start=[nrec],count=[1],units="1e6km3")
    call nc_write(trim(fnm),"nmask",     dble(n_tot),dim1="time",start=[nrec],count=[1],units="na")
    
    call nc_write(trim(fnm),"ela",       gis%ela,dims=dims,start=[1,nrec],count=[nm13,1], &
                                            long_name="Equil. line altitude (ELA)",units="m")
    call nc_write(trim(fnm),"aar",       gis%aar,dims=dims,start=[1,nrec],count=[nm13,1], &
                                            long_name="Accum. area ratio (AAR)",units="1")
    
    call nc_write(trim(fnm),"aar_ann",   gis%aar(13),dim1="time",start=[nrec],count=[1], &
                                            long_name="Accum. area ratio (AAR)",units="1")
    
    ! Boundary forcing
    ! call nc_write(trim(fnm),"dT_jja",sum(forcing_now%dTbnd(6:8))/3.d0,dim1="time", &
    !               start=[nrec],count=[1],units="degrees celcius")
    ! if (paleo_frac_dT .ne. 0.d0) then 
    !   call nc_write(trim(fnm),"dT_amp",forcing_now%dTamp,dims=dims,start=[1,nrec],count=[nm13,1],units="degrees celcius")
    ! end if 

    ! Now write daily output (for gis) - only if running rembo solo!!
    if ( clim_coupled .lt. 0 ) then
      dims(1) = "day"; dims(2) = "time"
      call nc_write(trim(fnm),"dma775",  dma775,dims=dims,start=[1,nrec],count=[nm13,1],long_name="daily melt area, min. 7.75mm/d")
      call nc_write(trim(fnm),"dma850",  dma850,dims=dims,start=[1,nrec],count=[nm13,1],long_name="daily melt area, min. 8.50mm/d")
      call nc_write(trim(fnm),"dma925",  dma925,dims=dims,start=[1,nrec],count=[nm13,1],long_name="daily melt area, min. 9.25mm/d")
    end if
    
    return

  end subroutine climchecker_new
  
  

  subroutine find_ela(ela,z_srf,smb,mask)
    ! Determine the ELA from smb values.

    implicit none 

    double precision, intent(OUT) :: ela 
    double precision, intent(IN)  :: z_srf(:,:) 
    double precision, intent(IN)  :: smb(:,:) 
    logical,          intent(IN)  :: mask(:,:) 

    ! Local variables 
    integer :: nx, ny 
    integer :: npts 
    double precision :: npts_dble 
    double precision, allocatable :: wts(:,:) 
    double precision :: xmin, xmax 
    double precision :: ymin, ymax
    double precision :: xbar, ybar 
    double precision :: b, a 

    nx = size(mask,1)
    ny = size(mask,2) 
    allocate(wts(nx,ny)) 

    ! Get mask of current points and count them
    npts      = count(mask)
    npts_dble = real(npts,8) 

    ! Now, find the ELA 

    if (npts .eq. 0) then 
      ! No points available, set ela to zero for now 
    
      ela = 0.0d0 

    else 
      ! Calculate ELA

      ! Determine range of input data points 
      xmin = minval(z_srf,mask=mask)
      xmax = maxval(z_srf,mask=mask)
      ymin = minval(smb,  mask=mask)
      ymax = maxval(smb,  mask=mask)
      
      if (ymin .gt. 0.d0) then 
        ela = maxval(z_srf,mask=mask) 
      else if (ymax .lt. 0.d0) then 
        ela = minval(z_srf,mask=mask) 
      else 
        ! Calculate ELA

        wts = 0.0
        where (mask .and. smb .ne. 0.0) wts = 1.0/smb 
        where (mask .and. smb .eq. 0.0) wts = maxval(wts)
        wts = wts / sum(wts) 

        ela = sum(wts*z_srf)

      end if 

    end if 

    return

  end subroutine find_ela


end module rembo_functions


