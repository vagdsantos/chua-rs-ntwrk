!! Written by Vagner dos Santos
!! Professor at State Univertsity of Ponta Grossa, Paraná State,  Brazil.
!!
!! Compile and execute using the following commands:
!!
!! INTEL Compiler:
!! ifort -fast chua_chimera.f90 -o x.x ; time ./x.x
!!
!! GNU Compiler:
!! gfortran -O3 -cpp chua_chimera.f90 -o x.x ; time ./x.x
!!
module precisionparams
  integer, parameter :: bit=8!    floating-point precision (4 = single precision; 8 = double precision)
  real(kind=bit), parameter :: pi = 4.e0_bit*atan(1.e0_bit)!      Pi constant
end module precisionparams
!
!_____________________________________________________________________________________
!
module netwrkparams
  USE precisionparams
!------------------------------- Network Parameters ----------------------------------
  integer, parameter :: n = 300,&!                        n = num. of nodes
                        neq = 3,&!                        neq = num. of equation per node
                        nt = n*neq!                       nt = total number of equations
  !real(kind=bit), parameter :: r = 1.e0_bit/3.e0_bit!    r = coupling range
  real(kind=bit) :: sigma = 0.75e0_bit!                   sigma = ring coupling strength
  real(kind=bit) :: kparam = 0.1e0_bit!                    kparam = central node coupling param
  integer, parameter :: p = 100!                          p = number of coupled neighbours
end module netwrkparams
!
!_____________________________________________________________________________________
!
module rk4params
  USE precisionparams
  implicit none
!------------------------ Numerical integration Parameters ---------------------------
  real(kind=bit), parameter :: tf = 5.0e3_bit,&!        tf = total integration time
                               trans = 3.5e3_bit,&!     trans = transient time
                               tstep = 1.e-2_bit!       tstep = time step
  integer, parameter :: nsteps = nint(tf/tstep)!        nsteps = number of integration steps
end module rk4params
!
!______________________________________________________________________________
!
module chua_params
  USE precisionparams
  implicit none
!--------------------------------- Chua Parameters -----------------------------------
  real(kind=bit), parameter :: alpha = 9.4e0_bit,&
                               beta = 14.28e0_bit,&
                               aparam = -1.143e0_bit,&
                               bparam = -0.714e0_bit
end module chua_params
!
!______________________________________________________________________________
!
!=========================== MAIN PROGRAM =====================================
program rede_chua
  use netwrkparams
  implicit none
  real(kind=bit) :: x(1:nt+neq), x0(1:nt+neq)!  x=dynamical variables; x0= initial conditions
  integer, allocatable :: seed(:)!        seed for the random number generator
  integer :: i, nrand
!---------------------- Random Number Generator -------------------------------
  call random_seed(size=nrand)!   defines the size of the seed variable
  allocate(seed(nrand))!          allocates memory for the 'seed' variable
  do i = 1, nrand
    seed(i) = 4243*(-1)**i + i
  enddo
  call random_seed(put=seed)!     initializes the random number generator using the specified seed
!------------------------------------------------------------------------------
  x0 = 0.e0_bit
  call random_number(x0)!         Generates a uniform random initial condition in the interval (0,1)

  x = x0
  call rk4sys(x,1)!               Calls the forth order Runge-Kutta numerical integrator
!------------------------------------------------------------------------------
  !print*, 'x1=',x0(1),'y1=',x0(2),'z1=',x0(3)
  !print*, ''
  print*,'FINISHED SUCCESSFULLY'
  print*, 'network size=',n
  print*, ''
  print*, 'p=',p,'sigma=',sigma,'k',kparam
  call system("gnuplot graph_png.plt")!     Plots the spacetime using GNUPLOT according the the script file "graph_png.plt"
end program rede_chua
!==============================================================================
!
!-------------------- NUMERICAL INTEGRATION SUBROUTINE ------------------------
subroutine rk4sys(x,l)
  USE netwrkparams
  USE rk4params
  real(kind=bit) :: x(1:nt+neq),y(1:nt+neq),f(1:nt+neq,4), lapl(1:nt)
  integer :: i,k,l
  character*60 :: filename(2)
!------------------------------------------------------------------------------
    write(unit=filename(1),fmt=666) l !open distinct files
    write(unit=filename(2),fmt=888) l !open distinct files
    open(Unit=3, File=filename(1), Status='unknown')
    open(Unit=5, File=filename(2), Status='unknown')
666 format('lapl',I0.2,'.dat')
888 format('at',I0.2,'.dat')
!------------------------------- Time Loop ------------------------------------
time:do k = 1, nsteps
!----------------------------- RK4 Integrator ---------------------------------
      call xpsys(x,f(:,1))
      y = x + 0.5e0_bit*tstep*f(:,1)
      call xpsys(y,f(:,2))
      y = x + 0.5e0_bit*tstep*f(:,2)
      call xpsys(y,f(:,3))
      y = x + tstep*f(:,3)
      call xpsys(y,f(:,4))
      x = x + (tstep/6.0e0_bit) * (f(:,1) + 2.0e0_bit*(f(:,2) + f(:,3)) + f(:,4))
!------------------------------AVOID_DIVERGENCE--------------------------------
      if(any(abs(x) .gt. (5.e2_bit))) then

        print*, "divergence"
        do i=1,nt+neq
          if(abs(x(i)) .gt. (5.e2_bit)) print*, 'i=',i
        enddo
        print*, 'div_time=',k*tstep
        exit
      endif
!------------------------------- DATA OUTPUT ----------------------------------
      if(k*tstep .gt. trans) then
        if(modulo(k,nint(1/tstep)) .eq. 0) write(5,333) x! Writes dynamical state of all nodes to file unit=5

        lapl = abs(-2e0_bit*x(1:nt) + cshift(x(1:nt),3,1) + cshift(x(1:nt),-3,1))!  Calculates the Laplacian of the ring nodes
        if(modulo(k,nint(1/tstep)) .eq. 0) write(3,333) lapl! Writes the laplacian of the ring nodes to file unit=3
      endif
    end do time
!
333 format(2000(f12.4,x))
end subroutine rk4sys
!
!==============================================================================
!
!--------------------------- EQUATIONS SUBROUTINE -----------------------------
!
subroutine xpsys(x,f)
  USE netwrkparams
  USE rk4params
  real(kind=bit) :: x(1:nt+neq), f(1:nt+neq), coup(2,n), coupsum(2)
  real(kind=bit) :: chua_x,chua_y,chua_z
  integer i,j,k,l
!----------------------- RING NETWORK EQUATIONS -----------------------------
  coup = 0.e0_bit
  coupsum = 0.e0_bit
net:do i = 1, nt, neq
      l = (i + 2)/3
      if(l .eq. 1) then!                Ring Coupling Function
acop:   do k=i-neq*p, i+neq*p,neq
          j = k
          if(k .le. 0) then
            j = nt + k
          elseif(j .gt. nt) then
            j = k - nt
          endif
          coupsum(1) = coupsum(1) + x(j)!       Summation for the x-variable from 1-p to 1+p (index=1 icluded)
          coupsum(2) = coupsum(2) + x(j+1)!     Summation for the y-variable from 1-p to 1+p (index=1 icluded)
        enddo acop
        coup(1,l) = coupsum(1) - (2*p + 1)*x(i)!       Coupling term of the x-variable from 1-p to 1+p 
        coup(2,l) = coupsum(2) - (2*p + 1)*x(i+1)!     Coupling term of the y-variable from 1-p to 1+p 
      else
        coupsum(1) = coupsum(1) + x(modulo(i+neq*p,nt)) - x(modulo(nt+(i-neq*(p+1)),nt))!       Summation for the x-variable from i-p to i+p (index=i icluded)
        coupsum(2) = coupsum(2) + x(modulo(i+neq*p,nt)+1) - x(modulo(nt+(i-neq*(p+1)),nt)+1)!     Summation for the y-variable from i-p to i+p (index=i icluded)

        coup(1,l) = coupsum(1) - (2*p + 1)*x(i)!           Coupling term of the x-variable from i-p to i+p
        coup(2,l) = coupsum(2) - (2*p + 1)*x(i+1)!         Coupling term of the y-variable from i-p to i+p
      endif!            End of Ring Coupling Function of node i
!
!---------------------------- RING NODES EQUATIONS -----------------------------
      f(i) = chua_x(x(i),x(i+1)) &!                   chua x-equation
      &    + sigma*coup(1,l)/(2.e0_bit*real(p)) &!    ring coupling function
      &    + (kparam/2.e0_bit)*(x(nt+1) - x(i))!      star coupling function

      f(i+1) = chua_y(x(i),x(i+1),x(i+2)) &!          chua y-equation
      &      + sigma*coup(2,l)/(2.e0_bit*real(p))!    ring coupling function

      f(i+2) = chua_z(x(i+1))!                        chua z-equation
    enddo net
!-------------------------- CENTRAL NODE EQUATION -----------------------------
    f(nt+1) = chua_x(x(nt+1),x(nt+2)) &!        x-variable of the central node
    &     + (kparam/2.e0_bit)*( sum(x(1:nt:neq)) - n*x(nt+1) )!  central node coupling function
  
    f(nt+2) = chua_y(x(nt+1),x(nt+2),x(nt+3))!       y-variable of the central node
  
    f(nt+3) = chua_z(x(nt+2))!                   z-variable of the central node
end subroutine xpsys
!==============================================================================
!
!------------------------------ CHUA FUNCTIONS --------------------------------
!
function chua_x(x,y) result(fx)
  use chua_params
  real(kind=bit), intent(in) :: x, y!  input
  real(kind=bit) :: fx!  output

  fx = alpha*( y - x - ( bparam*x + 0.5e0_bit*( aparam - bparam )*( abs(x + 1.e0_bit) - abs(x - 1.e0_bit))))
end function
!
function chua_y(x,y,z) result(fy)
  use chua_params
  real(kind=bit), intent(in) :: x, y, z!  input
  real(kind=bit) :: fy!  output

  fy = x - y + z
end function
!
function chua_z(y) result(fz)
  use chua_params
  real(kind=bit), intent(in) :: y!  input
  real(kind=bit) :: fz!  output

  fz = -beta*y
end function
!==============================================================================
