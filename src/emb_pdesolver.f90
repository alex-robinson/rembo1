module emb_pdesolver
   
  use emb_global
  
  implicit none
  
contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  e x p l i c i t
  ! Purpose    :  Perform 1 timestep of explicit routine on 2d grid
  ! Author     :  Alex Robinson (24. June 2008)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine explicit(u0,uu,F,alpha,gama,beta,relax)

    implicit none
    
    integer, parameter :: nx = nxe, ny = nye
    
    real (8), Intent(IN), dimension(ny,nx) :: u0, F, alpha, relax
    real (8), Intent(IN) :: gama, beta
    real (8), dimension(ny,nx) :: uu, u, bet, alph, rel
    real (8) :: gam
    integer :: i, j, k    
    
    ! Constants
    alph = alpha
    gam  = gama
    bet  = 1.d0-4.d0*alph - gam*beta + gam*relax
    
    rel = gam*relax*u0
    
    ! Update boundary conditions
    uu(1:ny,1)  = u0(1:ny,1)
    uu(1:ny,nx) = u0(1:ny,nx)
    uu(1,1:nx)  = u0(1,1:nx)
    uu(ny,1:nx) = u0(ny,1:nx)
    
    ! Define matrix (to fill in boundary regions)
    u = uu
    
    ! explicit scheme
    do j= 2, ny-1
      do i= 2, nx-1
        u(j,i)= bet(j,i) * uu(j,i)   &
                + alph(j,i) * (uu(j,i+1)+uu(j,i-1)+uu(j+1,i)+uu(j-1,i)) &
                + gam*F(j,i) - rel(j,i)
      end do
    end do
      
    ! Send to output matrix
    uu = u
    
    return
    
  end subroutine explicit

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  a d i
  ! Purpose    :  Perform 1 timestep of adi routine on 2d grid
  !               kappa varying in space (prescribed)
  ! Author     :  Alex Robinson (01. Aug 2008)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine adi(u0,uu,F,alpha,gama,beta,relax)

    implicit none
    
    integer, parameter :: nx = nxe, ny = nye
    
    real (8), intent(IN), dimension(ny,nx) :: u0, F, alpha, relax
    real (8), intent(IN) :: gama, beta
    real (8), dimension(ny,nx) :: uu, alph, coeffs, rel
    real (8), dimension(nx) :: ax, bx, cx, uprimex, cprimex, mx, rx
    real (8), dimension(ny) :: ay, by, cy, uprimey, cprimey, my, ry
    real (8) :: gam, bet0, bet
    integer :: i, j, k
    
    ! Constants, divide by two because it is 1/2 timestep
    alph = alpha / 2.d0
    gam = gama / 2.d0
    bet = 0.d0             ! Explicit coefficient, 0 if using implicit
    bet0 = beta            ! Implicit coefficient, 0 if using explicit
    
    rel = gam*relax*u0     ! Assign 'boundary' values to relaxation matrix
    coeffs = 1.d0 - 2.d0*alph - gam*bet + gam*relax
    
    ! Define the three diagonal rows of the matrix: x
    ax = 0.d0
    bx = 1.d0
    cx = 0.d0

    ! Define the three diagonal rows of the matrix: y
    ay = 0.d0
    by = 1.d0
    cy = 0.d0
    
    ! Update boundary conditions
    uu(1,1:nx)  = u0(1,1:nx)
    uu(ny,1:nx) = u0(ny,1:nx)
    uu(1:ny,1)  = u0(1:ny,1)
    uu(1:ny,nx) = u0(1:ny,nx)
    
    ! First implicit x, explicit y - column by column
    do j = 2, ny-1
      
      cprimex(1) = cx(1) / bx(1)     
      
      ! RHS
      uprimex(1) = ( coeffs(j,1) * uu(j,1)     &
                    + alph(j,1) * (uu(j+1,1) + uu(j-1,1))    &
                    + gam*F(j,1) - rel(j,1) ) / bx(1)     
      
      ! Define the three diagonal rows of the matrix: x
      ax(2:nx-1) = -alph(j,2:nx-1)
      bx(2:nx-1) = 1 + 2.d0 * alph(j,2:nx-1) + gam*bet0
      cx(2:nx-1) = -alph(j,2:nx-1) 
                
      do i = 2, nx

        cprimex(i) = cx(i) / (bx(i) - ax(i) * cprimex(i-1))
        
        ! RHS
        uprimex(i) = ( ( coeffs(j,i) * uu(j,i)      &
                    + alph(j,i) * (uu(j+1,i) + uu(j-1,i))      & 
                    + gam*F(j,i) - rel(j,i) ) - ax(i)    &
                    * uprimex(i-1)) / (bx(i) - ax(i) * cprimex(i-1))        
        
      end do 
      
      do i = nx-1, 2, -1
        uu(j,i) = uprimex(i) - cprimex(i) * uu(j,i+1)
      end do
      
    end do ! End of j loop
    
    ! Now explicit x, implicit y, row by row
    do i = 2, nx-1
      
      cprimey(1) = cy(1) / by(1)
      
      ! RHS
      uprimey(1) = ( coeffs(1,i) * uu(1,i)      &
                    + alph(1,i) * (uu(1,i+1) + uu(1,i-1))   &  
                    + gam*F(1,i) - rel(1,i) ) / by(1)
      
      ! Define the three diagonal rows of the matrix: y
      ay(2:ny-1) = -alph(2:ny-1,i)
      by(2:ny-1) = 1 + 2.d0 * alph(2:ny-1,i) + gam*bet0
      cy(2:ny-1) = -alph(2:ny-1,i)
    
      do j = 2, ny

        cprimey(j) = cy(j) / ( by(j) - ay(j) * cprimey(j-1) )
        
        ! RHS               
        uprimey(j) = ( ( coeffs(j,i) * uu(j,i)      &
                    + alph(j,i) * (uu(j,i+1) + uu(j,i-1))     &  
                    + gam*F(j,i) - rel(j,i) ) - ay(j)    &
                    * uprimey(j-1)) / (by(j) - ay(j) * cprimey(j-1))                  

      end do
      
      do j = ny-1, 2, -1
        uu(j,i) = uprimey(j) - cprimey(j) * uu(j+1,i)
      end do
      
    end do  ! End of i-loop
    
    return
    
  end subroutine adi
      
end module emb_pdesolver
