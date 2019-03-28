! ####################################################################
! Module     : exchange
! Author     : Alex Robinson
! Purpose    : To facilitate transfer of information between models
! ####################################################################
module exchange
  
  use parameters
  
  implicit none
  
  !! Structure containing driver information needed by main module
  !! (but obtained from individual module)
  type driver_info
    integer :: init, to_ice, to_clim
    integer :: ny, nx
    double precision :: t0, t1, dt_ice, dt_surf, dt_clim
  end type
  
  type(driver_info) :: info
  
  !! Structures containing information that each module will want
  type vars_to_ice
    double precision :: tt, tdjf, tjja, ttp, tts, tjan, tjul
    double precision :: smb, melted_ice, evap, refrozen
    double precision :: precip, snow, runoff_snow, runoff_rain
    double precision :: h_snow
  end type
  
  type vars_to_clim
    double precision,  allocatable, dimension(:,:) :: zs
    integer,  allocatable, dimension(:,:)          :: m2
    double precision :: dVdt
  end type

  ! Define the structures
  type(vars_to_ice),  allocatable, dimension(:,:) :: to_ice
  type(vars_to_clim) :: to_clim

contains
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : init_exchange
  ! Author     : Alex Robinson
  ! Purpose    : Initialize structures that will be used to transfer
  !              information between modules
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine init_exchange(coupled,ny,nx,t0,t1,dt_ice,dt_clim,dt_surf)
                           
   
    implicit none
    
    integer, optional :: coupled, ny, nx
    double precision, optional :: t0, t1, dt_ice, dt_clim, dt_surf
    
    if ( present(t0) )       info%t0 = t0
    if ( present(t1) )       info%t1 = t1
    if ( present(dt_ice) )   info%dt_ice = dt_ice
    if ( present(dt_clim) )  info%dt_clim = dt_clim
    if ( present(dt_surf) )  info%dt_surf = dt_surf
    
    if ( present(ny) .and. present(nx) .and. .not. allocated(to_ice) ) then
      info%ny = ny
      info%nx = nx
      allocate( to_ice(ny,nx) )
      allocate( to_clim%zs(ny,nx), to_clim%m2(ny,nx) )
      
      info%init    = 1
      info%to_clim = 0
      info%to_ice  = 0
      write(*,*) "Allocated to_ice and to_clim:", ny, " x ",nx
    end if
    
    ! Ensure that if models are not coupled, then no information 
    ! exchange should take place
    if ( present(coupled) ) then
      if ( coupled .eq. 0 ) info%init = -1
    end if
    
    return
    
  end subroutine init_exchange
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : vars2clim
  ! Author     : Alex Robinson
  ! Purpose    : Send variables to climate module. Arguments are climate
  !              module variables. Information the climate module would
  !              like to retreive is contained in the to_clim structure
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine vars2clim(zs, m2, dVdt, todo)
    
    implicit none
    
    character(len=*), optional :: todo
    character(len=4) :: c
    double precision, dimension(:,:) :: zs, m2
    double precision :: dVdt
    
    c = ""
    if ( present(todo) ) c = todo    
    
    if ( c .eq. "save" ) then
      
      if ( info%init .eq. 1 ) then
      
        to_clim%zs = zs
        to_clim%m2 = int(m2)
        to_clim%dVdt = dVdt
        
        info%to_clim = 1
        
        !write(*,*) "vars2clim: variables updated."
        
      end if
    
    else
    
      if ( info%to_clim .eq. 1 ) then
      
        zs   = to_clim%zs
        m2   = dble(to_clim%m2)
        dVdt = to_clim%dVdt
        info%to_clim = 0
        
        !write(*,*) "vars2clim: variables retrieved."
        
      end if
    
    end if
    
    return

  end subroutine vars2clim
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : vars2ice
  ! Author     : Alex Robinson
  ! Purpose    : Send variables to ice module. Arguments are ice
  !              module variables. Information the ice module would
  !              like to retreive is contained in the to_ice structure
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine vars2ice(tt,tdjf,tjja,ttp,tts,tjan,tjul, &
                      precip, snow, runoff_snow, runoff_rain, melted_ice, &
                      refrozen, evap, smb, h_snow, todo)
                       
    implicit none
    
    character(len=*), optional :: todo
    character(len=4) :: c
    double precision, dimension(:,:) :: precip, snow, runoff_snow, runoff_rain, melted_ice, &
                                        refrozen, evap, smb, tt, tdjf, tjja, ttp, tts, tjan, tjul    
    double precision, dimension(:,:) :: h_snow 

    c = ""
    if ( present(todo) ) c = todo
    
    if ( c .eq. "save" ) then
      
      if ( info%init .eq. 1 ) then
        
        to_ice%tt          = tt
        to_ice%tjan        = tjan
        to_ice%tjul        = tjul
        to_ice%tdjf        = tdjf
        to_ice%tjja        = tjja
        to_ice%ttp         = ttp
        to_ice%tts         = tts
        to_ice%precip      = precip
        to_ice%snow        = snow
        to_ice%runoff_snow = runoff_snow
        to_ice%runoff_rain = runoff_rain
        to_ice%melted_ice  = melted_ice
        to_ice%refrozen    = refrozen
        to_ice%evap        = evap
        to_ice%smb         = smb
        to_ice%h_snow      = h_snow
        info%to_ice        = 1
        
        write(*,*) "vars2ice: variables updated."
        
      end if
    
    else
    
      if ( info%to_ice .eq. 1 ) then
        
        tt          = to_ice%tt
        tjan        = to_ice%tjan
        tjul        = to_ice%tjul
        tdjf        = to_ice%tdjf
        tjja        = to_ice%tjja
        ttp         = to_ice%ttp
        tts         = to_ice%tts
        precip      = to_ice%precip
        snow        = to_ice%snow
        runoff_snow = to_ice%runoff_snow
        runoff_rain = to_ice%runoff_rain
        melted_ice  = to_ice%melted_ice
        refrozen    = to_ice%refrozen
        evap        = to_ice%evap
        smb         = to_ice%smb
        h_snow      = to_ice%h_snow
        info%to_ice = 0
        
        write(*,*) "vars2ice: variables retrieved."
        
      end if
      
    end if
    
    return

  end subroutine vars2ice



end module exchange
