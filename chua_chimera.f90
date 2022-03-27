!! Written by Vagner dos Santos
!! Professor at State Univertsity of Ponta Grossa, Paraná State,  Brazil.
!!
!! Compile and execute using the following commands:
!!
!! INTEL Compiler:
!! ifort -qopenmp -qmkl chua_chimera.f90 -o x.x ; time ./x.x
!!
!!
!______________________________________________________________________________
!
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
    integer, parameter :: n = 150,&!                        n = num. of nodes
                          neq = 3,&!                        neq = num. of equation per node
                          nt = n*neq!                       nt = total number of equations
    integer, parameter :: p = 50!                           p = number of coupled neighbours
    !real(kind=bit), parameter :: r = 1.e0_bit/3.e0_bit!    r = coupling range
    real(kind=bit), parameter :: sigma = 0.68e0_bit!        sigma = ring coupling strength
    real(kind=bit) :: kparam = 0.0e0_bit!                   kparam = central node coupling param
    real(kind=bit) :: x0ini = -1.e0_bit,&!		              x0ini = beginning of the of the x0 range
                      x0range = 2.e0_bit,&!		            x0range = range of x0
                      y0ini = -1.e0_bit,&!		              y0ini = beginning of the of the y0 range
                      y0range = 2.e0_bit!                  y0range = range of y0  
    integer, parameter :: grid = 512,&!                      grid = grid for the basin of attraction
                          num_threads = 8!			     num_threads = number of threads to parallelize the program.
  end module netwrkparams
  !
  !_____________________________________________________________________________________
  !
  module rk4params
    USE precisionparams
    implicit none
  !------------------------ Numerical integration Parameters ---------------------------
    real(kind=bit), parameter :: tf = 4.6e3_bit,&!        tf = total integration time
                                 trans = tf-0.1e3_bit,&!     trans = transient time
                                 tstep = 1.e-2_bit,&!       tstep = time step
                                 deltamin = 0.02e0_bit,&!      deltamin = threshold for chimera detection
                                 deltamax = 2.e0_bit!      deltaax = threshold SW-DW detection
    integer, parameter :: nsteps = nint(tf/tstep),&!        nsteps = number of integration steps
                          samplesize = nint((tf-trans)/tstep)!    samplesize = size of the timeseries used for analysis
    integer :: state = 1!                state = 1 after the time numerical integrations means no problem. state = -1 means the trajectory goes to infinity.
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
                                 aparam = -1.143_bit,&
                                 bparam = -0.714_bit
  end module chua_params
  !
  !=========================== MAIN PROGRAM =====================================
  program chua_network
    use netwrkparams
    use rk4params
    use omp_lib
    implicit none
    real(kind=bit) :: x(1:nt+neq),&!    x=dynamical variables
                      x0(1:nt+neq),&!   x0= initial conditions
                      xdata(1:nt+neq,samplesize)!   xdata = data output
  !  integer, allocatable :: seed(:)!        seed for the random number generator
    integer :: k,l,id, &!      general purpose integers
  !             nrand, &!    size of the seed() variable (to be defined by the code)
  !             usrseed, &!   user defined seed
               diagnostic(0:grid,0:grid)!	diagnostic = -1(divergence), 1(sync/cluster), 2(SW-chimera), 3(DW-chimera), 4 (incoherence)
    character*100 :: filename(1)
  !---------------------- Random Number Generator -------------------------------
    !print*, "Enter a seed (integer):"
    !read(*,*) usrseed
    !read(*,*) kparam
    !kparam = kparam/10.e0_bit
  
  
  !  usrseed = 84031470
  !  usrseed = 999016234
  !  call random_seed(size=nrand)!   defines the size of the seed variable
  !  allocate(seed(nrand))!          allocates memory for the 'seed' variable
  !  do i = 1, nrand
  !    seed(i) = usrseed*(-1)**i + i
  !  enddo
  !  call random_seed(put=seed)!     initializes the random number generator using the specified seed
  !  call random_seed()
  !  
  !  call random_number(x0)!         Generates a uniform random initial condition in the interval (0,1)
  !
  !
    call OMP_SET_NUM_THREADS(num_threads)
  
  !$OMP PARALLEL PRIVATE(k,l,x,x0,id,state,xdata,filename)
    x0 = 0.e0_bit! 				All nodes at an unstable fixed point at the origin
    id = omp_get_thread_num()
    
  
  !------------------- Output files ---------------------------------------------- 
    write(unit=filename(1),fmt=777) sigma, kparam, x0(3),deltamin,tf!       open distinct output files
    open(Unit=id+10, File=filename(1), Status='unknown')
  !
  777 format('./basin-sigma=',f4.2,'-kparam=',f4.2,'-z0=',f4.2,'-delta=',f4.2,'-tf=',f5.0,'dat')!         Format specifiers fot the name of the dynamical variables timeseries output file
    
  !------------------- IC_Loop ---------------------------------------------- 
  !$OMP DO
  p1:do k = 1, grid
      x0(1) = x0ini + k*x0range/grid!    x_1 initial condition
      !if(id .eq. 0) print*, 'thread id:',id, 'loop index:', k,'of',grid/num_threads
  p2: do l = 1, grid
        x0(2) = y0ini + l*y0range/grid!  y_1 Initial condition
  !      print*, 'k=',k,'l=',l, 'id=',id
  
        x = x0
        state = 1
        call rk4sys(x,xdata)!               Calls the forth order Runge-Kutta numerical integrator
  !
        call final_dyn_state(xdata,diagnostic(k,l))
  !    
  !------------------- TIMESERIES Output ---------------------------------------------- 
  !      write(unit=filename(1),fmt=234) k,l,deltamin!       open distinct output files
  !      open(Unit=id+10, File=filename(1), Status='unknown')
  !
  !234   format('./at-k=',I0,'-l=',I0,'-deltamin=',f4.2'.dat')!         Format specifiers fot the name of the dynamical variables timeseries output file
  !
        !do i = 1, samplesize
        !  write(id+10,333) xdata(:,i)
        !enddo 
        !close(id+10)
      enddo p2
    enddo p1
  !$OMP ENDDO
  
  
o1:do k = 1, grid
    x0(1) = x0ini + k*x0range/grid!    x_1 initial condition
o2: do l = 1, grid
      x0(2) = y0ini + l*y0range/grid!  y_1 Initial condition
      write(id+10,222) x0(1), x0(2), diagnostic(k,l)!         Write output to file
    enddo o2
  enddo o1
  
  
  
  close(id+10)
  !$OMP END PARALLEL
  !
  
  !333 format(2000(f12.4,x))!  Format of the output data
  222 format(2(f10.6,x), I2)!  Format of the output data
  !------------------------------------------------------------------------------
  !  print*,'FINISHED SUCCESSFULLY'
    !print*, 'network size=',n, 'user seed=', usrseed
    !print*, ''
    !print*, 'p=',p,'sigma=',sigma,'k',kparam
  ! call system("gnuplot graph_svg.plt")!     Plots the spacetime using GNUPLOT according the the script file "graph_png.plt"
  ! call system("cat basin-sigma=* >> basin.dat")
  end program chua_network
  !==============================================================================
  !
  !-------------------- FINAL STATE DIAGNOSTIC SUBROUTINE ------------------------
  subroutine final_dyn_state(xdata,diagnostic)
    use netwrkparams
    use rk4params
    implicit none
    real(kind=bit) :: xdata(1:nt+neq,samplesize),&!   xdata = data output
                      distmat(1:n+1,1:n+1),&!   distmat = Euclidean distance matrix
                      recmat(1:n+1,1:n+1),& !         recmat = recurrence matrix
                      evals(n),&!
                      evalsum(n),&!
                      evecs(n,n),&!
                      signvec(n)!                signvec = signed matrix
    real(kind=bit) :: delta,&!                delta = threshold to identify coherent regions
                      csum!       csum = sum of the signed vector components
    integer :: i,j,m, &!      general purpose integers
               diagnostic!    diagnostic = -1(divergence), 1(sync/SW-cluster), 2(DW-cluster), 3(SW-chimera), 4(DW-chimera), 5(incoherence)
  !
    if(state .eq. 1) then
        
  
  th: do m = 1,2
  
        delta = deltamin
        if(m .eq. 2) delta = deltamax
  
  

        csum = 0
        evalsum = 0
ts:     do i=1, samplesize, nint(1/tstep)
    
          do j=1, n+1
            distmat(:,j) = xdata(1:nt+neq:neq,i)
          enddo
          distmat(:,:) = abs(distmat - transpose(distmat))
    
          recmat = 0.0
          where(distmat(:,:) < delta)
            recmat = 1.0
          endwhere
    
          call eigens(recmat(1:n,1:n),evals,evecs)
    
          j = 1
          signvec = 0.e0
          do while(evals(j) .gt. 1.0)
              signvec = signvec + abs(evecs(:,j))
              j = j + 1
          enddo
    
          where(signvec > 1.e-8)
              signvec = 1.0
          endwhere
    
          csum = csum + sum(signvec)
          evalsum = evalsum + evals
        enddo ts
        csum = csum/real(samplesize*tstep)
        
  
        if(m == 1) then
          select case(nint(csum))
            case(1:149) 
              diagnostic = 3!'SW-chimera'
            case(0)
              diagnostic = 5!'Incoherence'
            case(150) 
              diagnostic = 1!'Sync/Sw-cluster'
          end select
        else if (m == 2) then
          if((evals(2) > 1.e0_bit) .and. (diagnostic .eq. 3)) diagnostic = 4!'DW-chimera'
          if((evals(2) > 1.e0_bit) .and. (diagnostic .eq. 1)) diagnostic = 2!'DW-cluster'
        endif
  
      enddo th
    elseif(state == -1) then
      diagnostic = -1!'Divergence'
  endif
  
  end subroutine final_dyn_state
  !==============================================================================
  !
  !-------------------- NUMERICAL INTEGRATION SUBROUTINE ------------------------
  subroutine rk4sys(x,xout)
    USE netwrkparams
    USE rk4params
    real(kind=bit) :: x(1:nt+neq),y(1:nt+neq),f(1:nt+neq,4)
    real(kind=bit),intent(out) :: xout(1:nt+neq,samplesize)
    integer :: k
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
        if(any(abs(x) .gt. (1.e5_bit))) then
  
          !print*, "ERROR!!", "DIVERGENT NODE!!"
          !do i=1,nt+neq
          !  if(abs(x(i)) .gt. (5.e2_bit)) print*, 'node =',i,'diverged at time T=',k*tstep
          !enddo
          state = -1
  
          exit
        endif
  !------------------------------- DATA OUTPUT ----------------------------------
  !      if((k*tstep .gt. trans) .and. (modulo(k,nint(1/tstep)) .eq. 0)) xout(:,nint(k*tstep-trans)) = x! Writes dynamical state of all nodes to output variable!
        if(k*tstep .gt. trans) xout(:,k-nint(trans/tstep)) = x! Writes dynamical state of all nodes to output variable!
      end do time
  !
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
          coup(1,l) = coupsum(1) - (2.e0_bit*p + 1)*x(i)!       Coupling term of the x-variable from 1-p to 1+p 
          coup(2,l) = coupsum(2) - (2.e0_bit*p + 1)*x(i+1)!     Coupling term of the y-variable from 1-p to 1+p 
        else
          coupsum(1) = coupsum(1) + x(modulo(i+neq*p,nt)) - x(modulo(nt+(i-neq*(p+1)),nt))!       Summation for the x-variable from i-p to i+p (index=i icluded)
          coupsum(2) = coupsum(2) + x(modulo(i+neq*p,nt)+1) - x(modulo(nt+(i-neq*(p+1)),nt)+1)!     Summation for the y-variable from i-p to i+p (index=i icluded)
  
          coup(1,l) = coupsum(1) - (2.e0_bit*p + 1)*x(i)!           Coupling term of the x-variable from i-p to i+p
          coup(2,l) = coupsum(2) - (2.e0_bit*p + 1)*x(i+1)!         Coupling term of the y-variable from i-p to i+p
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
  !--------------------------- EQUATIONS SUBROUTINE -----------------------------
  !
  subroutine eigens(A,evals,evecs)
    USE precisionparams
    use netwrkparams
    implicit none
    integer, parameter :: LDA = n, &!         LDA: The leading dimension of the array A. LDA >= max(1,N).
                          LWORK = 4*n,&!      LWORK = The length of the array WORK.  LWORK >= max(1,3*N-1).
                          LDVL = n, &!           LDVL: The leading dimension of the array VL
                          LDVR = n!           LDVR: The leading dimension of the array VR
    INTEGER :: i,j, &
               INFO!                          INFO = 0: Successful exit; INFO /= 0: Error.
    REAL(KIND=bit) :: A(n,n), &!            A: The symmetric matrix
                      WR(n), &!             WR: If INFO = 0, the RIGHT eigenvalues in ascending order.
                      WI(n), &!             WI: Imaginaty part of the eigenvalues
                      VL(LDVL,n), &!        VL: If JOBVL = 'V', the LEFT eigenvectors u(j) are stored one after another in the columns of VL, in the same order as their eigenvalues.
                      VR(LDVR,n), &!        VR: If JOBVR = 'V', the RIGHT eigenvectors u(j) are stored one after another in the columns of VL, in the same order as their eigenvalues.
                      WORK(LWORK)!      WORK: On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
    
    REAL(KIND=bit) :: evals(n),evecs(n,n)
  
    CHARACTER :: JOBVL, &!   JOBVL = 'n': Compute LEFT eigenvalues; JOBVL = 'V': Compute LEFT eigenvalues and LEFT eigenvectors.
                 JOBVR, &!   JOBVR = 'n': Compute RIGHT eigenvalues; JOBVR = 'V': Compute RIGHT eigenvalues and RIGHT eigenvectors.
                 UPLO!      UPLO = 'U': Upper triangle of A is stored; JOBVL = 'L': Lower triangle of A is stored.
  !  CHARACTER*60 :: FILENAME(3)
    !______________________________________________________________________________
    JOBVL = 'V'; JOBVR = 'V'
    UPLO = 'U'
  
  
    call dgeev(JOBVL,JOBVR, N, A, LDA,WR,WI,VL, LDVL,VR, LDVR,WORK, LWORK, INFO)! calculates the eigenvalues and eigenvectors
  
  !  write(unit=filename(2),fmt=777) delta!open distinct output files
  !  write(unit=filename(3),fmt=888) delta!open distinct output files
  !  open(Unit=11, File=filename(2), Status='unknown')
  !  open(Unit=13, File=filename(3), Status='unknown')
  !  777 format('./data/eigenvals-delta=',f5.3,'.dat')!       Format specifiers fot the name of the eigenvalues output file
  !  888 format('./data/eigenvecs-delta=',f5.3,'.dat')!       Format specifiers fot the name of the eigenvectors output file
  
    j=1
    do while(maxval(WR,1) .ge. -100.e0)
  
      i = maxloc(WR,1)
      !print*, maxloc(WR),i, WR(i)
    
      evals(j) = WR(i)
      evecs(:,j) = VR(:,i)
      j = j + 1
  
  !    write(11,333) real(i), WR(i), WI(i)
  !    write(13,333) (abs(VR(k, i)), k=1, N)
  
      WR(i) = -1000.0
    enddo
  
  !  close(11)
  !  close(13)
  
  !  333 format(2000(f12.4,x))!  Format of the output data
  end subroutine eigens
  
