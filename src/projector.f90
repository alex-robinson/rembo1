module projector

  use emb_global
  
  implicit none
    
contains

  !  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !  Subroutine :  c r e a t e g r i d
  !  Purpose    :  Create a cartesian grid, find the lat/lon
  !                values that correspond to the grid (via 
  !                inverse stereographic projection
  !  Author     :  Alex Robinson (15. Apr 2008)
  !  updated    :  
  !  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine creategrid(grid,x0,y0,dx,dy,nx,ny)
      
    implicit none
    
    type(location) :: grid(:,:)
    
    integer :: i, j, nx, ny
    real (8) :: x0, y0, dx, dy, tmpx, tmpy
    
    write(*,*) "Creating cartesian grid..."
        
    do i = 1, nx
      grid(i,:)%x = x0 + (i-1)*dx*10**3
    end do
    do j = 1, ny
      grid(:,j)%y = y0 + (j-1)*dx*10**3
    end do
    
    ! Perform transformation to sphere/ellipsoid from plane
    ! for all points on the grid
    do i = 1, nx
      do j = 1, ny
        call plane_ellipsoid(grid(i,j)%lon, grid(i,j)%lat, &
                             grid(i,j)%x, grid(i,j)%y, -1)
        grid(i,j)%lon = grid(i,j)%lon * todegs
        grid(i,j)%lat = grid(i,j)%lat * todegs
      end do
    end do
    
    call check_grid(grid)
    
    return
    
  end subroutine creategrid

  !  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !  Subroutine : c h e c k _ g r i d
  !  Purpose    : Output some information about the current grid
  !  Author     : Alex Robinson (16. Apr 2008)
  !  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine check_grid(grid)
  
    implicit none
  
    type(location) :: grid(:,:)
  
    integer :: imax, jmax
    real (8) :: tmpx, tmpy
    
    imax = size(grid,1)
    jmax = size(grid,2)
    
    ! Calculate the size of the grid in meters
    !tmpx = calcdistw(grid(1,1)%lon,grid(1,1)%lat,grid(nx,1)%lon,grid(nx,1)%lat)
    !tmpy = calcdistw(grid(1,1)%lon,grid(1,1)%lat,grid(1,ny)%lon,grid(1,ny)%lat)
    tmpx =0.d0
    tmpy =0.d0
    
    write(*,"(2x,a)") "========================================================"
    write(*,"(5x,a,i6,a2,i6)") "Grid points imax x jmax:        ", imax,"x",jmax
    !write(*,"(5x,a,f12.2,a2,f12.2)") "Grid size (length x height):", tmpx,"x",tmpy
    write(*,"(5x,a)") "    Corners of grid grid lat,lon:"
    write(*,*)
    write(*,"(5x,f12.1,a1,f12.1,10x,f12.1,a1,f12.1)") &
                                grid(1,jmax)%lon, ",", grid(1,jmax)%lat, &
                                grid(imax,jmax)%lon, ",", grid(imax,jmax)%lat
    write(*,*)
    write(*,"(5x,f12.1,a1,f12.1,10x,f12.1,a1,f12.1)") &
                                grid(1,1)%lon, ",", grid(1,1)%lat, &
                                grid(imax,1)%lon, ",", grid(imax,1)%lat            
    write(*,*)
    write(*,*) "    Corners of grid grid x,y:"
    write(*,*)  
    write(*,"(5x,f12.1,a1,f12.1,10x,f12.1,a1,f12.1)") &
                               grid(1,jmax)%x, ",", grid(1,jmax)%y, &
                               grid(imax,jmax)%x, ",", grid(imax,jmax)%y                                              
    write(*,*)
    write(*,"(5x,f12.1,a1,f12.1,10x,f12.1,a1,f12.1)") &
                               grid(1,1)%x,",", grid(1,1)%y, &
                               grid(imax,1)%x,",", grid(imax,1)%y                                                           
    write(*,*)
    return
    
  end subroutine
  
  !  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !  Subroutine :  p l a n e _ e l l i p s o i d
  !  Purpose    :  Computation of position (x,y) for latitude phi 
  !                and longitude lambda in a polar stereographics
  !                projection with reference to the WGS84 ellipsoid.
  !  Author     :  Reinhard Calov
  !  updated    :  Alex Robinson (17. Mar 2008)
  !  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine plane_ellipsoid(lambda, phi, x, y, proj)

    implicit none
    
    integer :: proj, l
    real (8) :: e, mc, t, tc, kp, rho, phi_p, eps
    real (8) :: lambda, phi, x, y
    real (8) :: sinphi0, sinlambda0, cosphi0, coslambda0
    
    e=sqrt((a_earth**2-b_earth**2)/(a_earth**2))
    
    sinphi0    = dsin(phi0);  sinlambda0 = dsin(lambda0)
    cosphi0    = dcos(phi0);  coslambda0 = dcos(lambda0)
  
    if (proj .eq. 1) then
      ! Perform stereographic projection: from ellipsoid to plane

      mc=cosphi0/dsqrt(1.d0-e*e*sinphi0*sinphi0)
      t=dsqrt(((1.d0-dsin(phi))/(1.d0+dsin(phi)))* &
             ((1.d0+e*dsin(phi))/(1.d0-e*dsin(phi)))**e)
      tc=dsqrt(((1.d0-sinphi0)/(1.d0+sinphi0))* &
              ((1.d0+e*sinphi0)/(1.d0-e*sinphi0))**e)
      rho=a_earth*mc*t/tc

      x= rho*dsin(lambda-lambda0)
      y=-rho*dcos(lambda-lambda0)

      !kp=0.5d0*mc*dsqrt((1.0d0+e)**(1.0d0+e)*(1.0d0-e)**(1.0d0-e))/ &
      !                        (a*tc)
      !write(6,'(a,1pe11.4,a)') 'kp= ', kp
    
    else if (proj .eq. -1) then
      ! Perform inverse projection: from plane to ellipsoid

      tc=dsqrt(((1.d0-sinphi0)/(1.d0+sinphi0))* &
             ((1.d0+e*sinphi0)/(1.d0-e*sinphi0))**e)
      mc=cosphi0/dsqrt(1.d0-e*e*sinphi0*sinphi0)
      rho=dsqrt(x*x+y*y)
      t=rho*tc/(a_earth*mc)

      lambda=lambda0+datan(x/(-y))

      !  fix point iteration

      phi_p=0.5d0*pi-2.d0*datan(t)
      l=0
      eps=3600.d0
      do while(eps.ge.1d-9)
        l=l+1
        phi=0.5d0*pi-2.d0*datan(t*((1.d0-e*dsin(phi_p))/ &
                 (1.d0+e*dsin(phi_p)))**(0.5d0*e))
        eps=dabs(phi-phi_p)
        phi_p=phi
      end do
      
    else
    
      stop "Error in arguments for plane_ellipsoid subroutine: proj must be -1 or 1"
      
    end if
    
    return
    
  end subroutine plane_ellipsoid 

  !  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !  Subroutine :  t e s t t r a n s
  !  Purpose    :  used to test projection algorithms
  !  Author     :  Alex Robinson (17. Mar 2008)
  !  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine testtrans()
  
    implicit none
    
    integer :: i, proj
    real (8) :: lambda, phi, x, y
    real (8) :: xy(4,2), lamphi(4,2)
    
    xy(1,:) = (/ -802500.d0, -3402500.d0 /)
    xy(2,:) = (/  702500.d0, -3402500.d0 /)
    xy(3,:) = (/ -802500.d0,  -597500.d0 /)
    xy(4,:) = (/  702500.d0,  -597500.d0 /)
    lamphi(1,:) = (/ -52.2710201d0, 58.6035491d0 /)
    lamphi(2,:) = (/ -27.3342981d0, 58.7883271d0 /)
    lamphi(3,:) = (/ -92.3305365d0, 80.8106358d0 /)
    lamphi(4,:) = (/  10.6177100d0, 81.5269340d0 /)
    
    write(*,*) "\n == Output from testtrans =="
    
    proj = 1       ! Perform transformation to plane from sphere/ellipsoid
    
    write(*,*) "\nplane_ellipsoid"
    do i=1,4
      write(6,*) i
      call plane_ellipsoid(lamphi(i,1)*torads, lamphi(i,2)*torads, x, y, proj)
      write(6,*) "x: ",xy(i,1)*0.001,"=",x*0.001, "(",xy(i,1)-x,")"
      write(6,*) "y: ",xy(i,2)*0.001,"=",y*0.001, "(",xy(i,2)-y,")"
      
    end do
    
    proj = -1      ! Perform transformation to sphere/ellipsoid from plane
    
    write(*,*) "\n\nplane_ellipsoid"
    do i=1,4
      write(6,*) i
      call plane_ellipsoid(lambda, phi, xy(i,1), xy(i,2) ,proj)    
      write(6,*) "lambda: ", lamphi(i,1),"=",lambda*todegs, &
               "(",lambda*todegs - lamphi(i,1),")" 
      write(6,*) "phi: ", lamphi(i,2),"=",phi*todegs, &
               "(",phi*todegs - lamphi(i,2),")"  
    end do
        
    return 
    
  end subroutine testtrans

end module projector
