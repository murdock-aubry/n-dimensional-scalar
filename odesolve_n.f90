!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for solving nonlinear systems of ordinary differential 
!  equations via a highly robust, but not terrible efficient, adaptive Chebyshev
!  spectral method. odesolve_nonlinear1
!
!  The following subroutines are publicly callable:
!
!    odesolve_init - initialize the solver by populating a structure which stores
!      data used by the other subroutines
!
!    odesolve_nonlinear - solve a "specified value problem" for a system of nonlinear
!      ordinary differential equations; that is, solve a problem of the form
!
!        {  y'(t) = f(t,y(t)),            a < t < b,                                      (1)
!        {  y(c)  = yc
!
!      where a <= c <= b, f is a smooth function R x C^n --> C^n and y is a
!      mapping R --> C^n. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module odesolve_n

    use utils
    use chebyshev
    use linalg0
    use chebpw
    use iso_c_binding
    use ieee_arithmetic
    
    type     odesolve_vars_t
    integer                                 :: kcheb
    integer                                 :: ntest
    integer                                 :: maxstack
    integer                                 :: maxints
    integer                                 :: maxiters
    integer                                 :: maxsteps
    
    double precision, allocatable           :: xscheb(:)
    double precision, allocatable           :: whtscheb(:)
    double precision, allocatable           :: acoefs(:,:)
    double precision, allocatable           :: aintl(:,:)
    double precision, allocatable           :: aintr(:,:)
    end type odesolve_vars_t
    
    
    interface          odesolve_nonlinear
    module procedure   odesolve_nonlinear1
    module procedure   odesolve_nonlinear2
    end interface      odesolve_nonlinear
    
    interface
    
    subroutine odesolve_fun1(n,t,y,f,df,userptr, sum_set, sum_coefs)
    use iso_c_binding
    implicit double precision (a-h,p-z)
    integer           :: n
    type(c_ptr)       :: userptr
    double complex    :: y(:), f(:), df(:,:)
    double precision  :: sum_coefs(:)
    integer           :: sum_set(:, :)
    end subroutine
    
    subroutine odesolve_fun2(n,t,y,f,df,pars)
    use iso_c_binding
    implicit double precision (a-h,p-z)
    integer                    :: n
    double complex, pointer    :: pars(:)
    double complex             :: y(:), f(:), df(:,:)
    end subroutine
    
    
    end interface
    
    contains
    
    subroutine odesolve_init(vars,kcheb0,ntest0,maxstack0,maxints0,maxsteps0,maxiters0)
    implicit double precision (a-h,o-z)
    type(odesolve_vars_t)                       :: vars
    integer, optional                           :: kcheb0, ntest0, maxstack0, maxsteps0
    integer, optional                           :: maxiters0, maxints0
    
    !
    !  Initialize the structure containing any data needed by the other procedures 
    !  in this module.
    !
    !  If vars is the only parameter passed to this routine, then reasonable defaults
    !  are chosen.
    !
    !  Input parameters:
    !    kcheb0 - the number of terms in the local Chebyshev expansions used to represent
    !      solutions
    !    ntest0 - the number of trailing terms in the local Chebyshev expansions to 
    !      check when testing the fitness of a purported solution
    !    maxstack0 - the length of the stack used by the adaptive procedure
    !    maxints0 - the maximum number of permissible intervals for the discretization
    !      scheme used to represent solutions
    !    maxsteps0 - the maximum number of iterations for the trapezoid rule used in the
    !                ODE solver
    !    maxiters0 - the maximum number of iterations for Newton's method
    !
    !  Output parameters:
    !    vars - the structure containing all of the data needed by the other procedures in
    !      this module
    !
    
    if (.not. present(kcheb0) ) then
    kcheb    = 16
    ntest    = 2
    maxstack = 100
    maxints  = 1000
    maxsteps = 8
    maxiters = 8
    else
    kcheb    = kcheb0
    ntest    = ntest0
    maxstack = maxstack0
    maxints  = maxints0
    maxsteps = maxsteps0
    maxiters = maxiters0
    endif
    
    vars%kcheb    = kcheb
    vars%ntest    = ntest
    vars%maxstack = maxstack
    vars%maxints  = maxints
    vars%maxsteps = maxsteps
    vars%maxiters = maxiters
    
    call chebyshev_quad(kcheb,vars%xscheb,vars%whtscheb)
    call chebyshev_coefsmatrix(kcheb,vars%acoefs)
    call chebyshev_intlmatrix(kcheb,vars%aintl)
    call chebyshev_intrmatrix(kcheb,vars%aintr)
    
    end subroutine
    
    
    
    subroutine odesolve_fit(vars,ifit,vals,dcoefs)
    implicit double precision (a-h,o-z)
    type(odesolve_vars_t)         :: vars
    double complex                :: vals(:,:), coefs(vars%kcheb)
    double precision              :: dd1(size(vals,2)),dd2(size(vals,2))
    
    double complex, allocatable   :: coefs0(:,:)
    
    n      = size(vals,2)
    kcheb  = vars%kcheb
    ntest  = vars%ntest
    dcoefs = 0
    
    
    if (ifit .eq. 1) then
    
    allocate(coefs0(kcheb,n))
    coefs0     = matmul(vars%acoefs,vals)
    
    dcoefs  = 0 
    ddd1 = maxval(abs(coefs0))
    if (ddd1 .eq. 0) ddd1 = 1
    
    do i=1,n
    ddd2    = maxval(abs(coefs0(kcheb-ntest+1:kcheb,i)))
    dcoefs  = max(dcoefs, ddd2/ddd1)
    end do
    
    else if (ifit .eq. 2) then
    coefs     = matmul(vars%acoefs,vals(:,1))
    ddd1      = maxval(abs(coefs)) + 1
    ddd2      = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
    
    dcoefs    = ddd2/ddd1
    
    else if(ifit .eq. 3) then
    
    allocate(coefs0(kcheb,n))
    coefs0     = matmul(vars%acoefs,vals)
    dcoefs     = 0 
    
    do i=1,n
    ddd2    = maxval(abs(coefs0(kcheb-ntest+1:kcheb,i)))
    dcoefs  = max(dcoefs, ddd2)
    end do
    
    else if (ifit .eq. 4) then
    
    coefs     = matmul(vars%acoefs,vals(:,1))
    dcoefs    = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
    
    endif
    
    end subroutine
    
    
    subroutine odesolve_nonlinear1(vars,ier,ifit,eps,a,b,c,yc,n,odefun,userptr,chebsol,ys,yders,sum_set, sum_coefs)
    implicit double precision (a-h,o-z)
    type(odesolve_vars_t)                       :: vars
    type(c_ptr)                                 :: userptr
    type(chebpw_scheme)                         :: chebsol
    procedure(odesolve_fun1)                    :: odefun
    double complex, allocatable, intent(out)    :: ys(:,:), yders(:,:)
    double complex                              :: yc(:)
    !
    !  Adaptively solve a problem of the form
    !
    !    { y'(t) = f( t, y(t) ),    a < t < b,
    !    {  y(c) = yc
    !
    !  where a <= c <= b and f: R x C^n \to C^n is a smooth function supplied
    !  via an external subroutine.   This routine uses a straightforward and rather slow spectral 
    !  method.  Each component of the solution is represented via a piecewise Chebyshev expansion.
    !
    !  Input parameters:
    !    vars - the structure populated by odesolve_init
    !
    !    ifit - an integer parameter which determines the criterion used for adaptive
    !      subdivision;
    !
    !      ifit = 1   means measure the relative accuracy of every component of the solution 
    !                 for fit
    !      ifit = 2   means only measure the relative accuracy of the first component for
    !                 fit
    !      ifit = 3   means measure the absolute accuracy of every component of the solution 
    !                 for fit
    !      ifit = 4   means only measure the absolute accuracy of the first component for
    !                 fit
    !
    !    eps - the desired precision for the solution
    !    (a,b) - the interval over which the problem is given
    !    c - the point at which the conditions for y(t) are specified
    !    n - the number of equations in the system
    !    odesolve - an external subroutine conforming to the odesolve_fun1
    !     interface which supplies the 
    !    userptr - a user-supplied "void *" pointer which is passed to odesolve
    !    yc - a vector giving the value of y at the point c
    !    
    !  Output parameters:
    !    ier - an error return code;
    !       ier  = 0   indicates successful execution
    !       ier  = 4   means the maximum number of intervals was exceeded
    !       ier  = 8   means the stack overflowed
    !
    !    chebsol - a structure specifying the piecewise discretization scheme
    !      used to represent the solutions
    !    ys - an array whose jth column gives the values of the jth component
    !      y_j of the solution at the discretization nodes
    !    yders - an array whose jth column gives the values of the derivative y_j'
    !     of the jth component of the solution at the discretization nodes     
    !
    
    double precision, allocatable       :: stack(:,:),  ts0(:), ab2(:,:), ab1(:,:), ab(:,:)
    double complex, allocatable         :: ys1(:,:,:), yders1(:,:,:)
    double complex, allocatable         :: ys2(:,:,:), yders2(:,:,:)
    double complex, allocatable         :: f0(:,:), df0(:,:,:), hs0(:,:), hders0(:,:)
    
    double precision                    :: sum_coefs(:)
    integer                             :: sum_set(:, :)
    
    
    ier        = 0
    
    k          = vars%kcheb
    maxsteps   = vars%maxsteps
    maxiters   = vars%maxiters
    maxstack   = vars%maxstack
    maxints    = vars%maxints
    
    
    allocate(stack(2,maxstack), ts0(k))
    allocate(ys1(k,n,maxints), ab1(2,maxints), yders1(k,n,maxints) )
    allocate(ys2(k,n,maxints), ab2(2,maxints), yders2(k,n,maxints) )
    allocate(f0(k,n), df0(k,n,n), hs0(k,n), hders0(k,n)  )
    
    !
    !  Solve backward from b to a
    !
    nints1     = 0
    
    if (a < c) then
    
        ind        = 0
        nstack     = 1
        stack(1,1) = a
        stack(2,1) = c
    
    
        do while (nstack .gt. 0)
    
            a0       = stack(1,nstack)
            b0       = stack(2,nstack)
            nstack   = nstack-1
            ts0      = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
    
            ifaccept = 0
            jer      = -1
            dy       = -1
            dyp      = -1
    
            if (nints1 .eq. 0 ) then
                ys1(k,:,nints1+1) = yc
            else
                ys1(k,:,nints1+1) = ys1(1,:,nints1)
            endif
    
            !
            !  Use the trapezoidal rule to build an initial approximation
            !
    
            call odesolve_traptvp1(jer,eps,maxsteps,a0,b0,k,ts0,n,odefun,userptr,ys1(:,:,nints1+1),&
                    yders1(:,:,nints1+1),sum_set,sum_coefs)
    
    
            ! Use spectral integration to compute the y from y'
    
            ! do i=1,k
            ! call odefun(n,ts0(i),ys1(i,:,nints1+1),f0(i,:),df0(i,:,:),userptr)
            ! yders1(i,:,nints1+1) = f0(i,:)
            ! end do
    
            do j=1,n
                ys1(:,j,nints1+1) = ys1(k,j,nints1+1)+(b0-a0)/2 * matmul(vars%aintr,yders1(:,j,nints1+1))
            end do
    
    
            if (jer .eq. 0) then
    
                do iter=1,maxiters
    
                    do i=1,k
                        call odefun(n,ts0(i),ys1(i,:,nints1+1),f0(i,:),df0(i,:,:),userptr, sum_set, sum_coefs)
                    end do
    
    
                    f0 = f0 - yders1(:,:,nints1+1)
                    hs0      = 0
                    call odesolve_lineartvp(vars,a0,b0,k,ts0,n,df0,f0,hs0,hders0)
    
                    dd1 = norm2(abs(ys1(:,:,nints1+1)))
                    if (dd1 .eq. 0) dd1 = 1
                    dd2 = norm2(abs(hs0))
    
                    ys1(:,:,nints1+1)    = ys1(:,:,nints1+1) + hs0
                    yders1(:,:,nints1+1) = yders1(:,:,nints1+1) + hders0
    
                    if (dd2 .lt. eps*dd1) exit
    
                end do
    
                call odesolve_fit(vars,ifit,ys1(:,:,nints1+1),dy)
                ! call odesolve_fit0(vars,ys1(:,:,nints1+1),dy)
    
                ! if (dy .lt. eps .AND. dyp .lt. eps) ifaccept = 1
    
                if (dy .lt. eps) ifaccept = 1
    
    
            endif
    
            ! ind = ind+1
            ! write(*, "(4(I6,1X),4(D15.8,1X),I4)")  ind,nints1,nstack,jer,a0,b0,dy,dyp,ifaccept
            ! write(13,"(4(I6,1X),4(D15.8,1X),I4)")  ind,nints1,nstack,jer,a0,b0,dy,dyp,ifaccept
    
            if (ifaccept .eq. 0) then
    
                if (nstack+2 .ge. maxstack) then
                    ier = 8
                    return
                endif
    
    
                nstack          = nstack+1
                stack(1,nstack) = a0
                stack(2,nstack) = (a0+b0)/2
    
                nstack          = nstack+1
                stack(1,nstack) = (a0+b0)/2
                stack(2,nstack) = b0
    
                else
    
                if (nints1+1 .ge. maxints) then
                    ier = 4
                    return
                endif
    
    
                nints1        = nints1+1
                ab1(1,nints1) = a0
                ab1(2,nints1) = b0
    
            endif
    
        end do
    
    endif
    
    
    
    !
    !  Solve forward from c to b.
    !
    
    nints2     = 0
    
    
    if (c < b) then
        ind        = 0
        nstack     = 1
        stack(1,1) = c
        stack(2,1) = b
    
    
        do while (nstack .gt. 0)
    
            a0       = stack(1,nstack)
            b0       = stack(2,nstack)
            nstack   = nstack-1
            ts0      = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
    
            ifaccept = 0
            jer      = -1
            dy       = -1
            dyp      = -1
    
            ys2(:,:,nints2+1) = 0
    
            if (nints2 .eq. 0 ) then
                ys2(1,:,nints2+1) = yc
            else
                ys2(1,:,nints2+1) = ys2(k,:,nints2)
            endif
    
            !
            !  Use the trapezoidal rule to build an initial approximation
            !
    
    
            call odesolve_trapivp1(jer,eps,maxsteps,a0,b0,k,ts0,n,odefun,userptr,ys2(:,:,nints2+1),&
                    yders2(:,:,nints2+1),sum_set,sum_coefs)
    
            
    
            ! Use spectral integration to compute the y from y'
            ! do i=1,k
            ! call odefun(n,ts0(i),ys2(i,:,nints2+1),f0(i,:),df0(i,:,:),userptr)
            ! yders2(i,:,nints2+1) = f0(i,:)
            ! end do
    
            do j=1,n
                ys2(:,j,nints2+1) = ys2(1,j,nints2+1)+(b0-a0)/2 * matmul(vars%aintl,yders2(:,j,nints2+1))
            end do
    
            if (jer .eq. 0) then
    
                do iter=1,maxiters
    
                    do i=1,k
                        call odefun(n,ts0(i),ys2(i,:,nints2+1),f0(i,:),df0(i,:,:),userptr, sum_set, sum_coefs)
                    end do
    
                    f0       = f0 - yders2(:,:,nints2+1)
                    hs0      = 0
                    call odesolve_linearivp(vars,a0,b0,k,ts0,n,df0,f0,hs0,hders0)
    
                    dd1 = norm2(abs(ys2(:,:,nints2+1)))
                    if (dd1 .eq. 0) dd1 = 1
                    dd2 = norm2(abs(hs0))
    
                    ys2(:,:,nints2+1)    = ys2(:,:,nints2+1) + hs0
                    yders2(:,:,nints2+1) = yders2(:,:,nints2+1) + hders0
    
                    if (dd2 .lt. eps*dd1) exit
    
                end do
    
                call odesolve_fit(vars,ifit,ys2(:,:,nints2+1),dy)
                !call odesolve_fit0(vars,ys2(:,:,nints2+1),dy)
                ! if (dy .lt. eps .AND. dyp .lt. eps) ifaccept = 1
                if (dy .lt. eps) ifaccept = 1
    
            endif
    
            
    
            ! ind = ind+1
            !write(*, "(4(I6,1X),4(D15.8,1X),I4)")  ind,nints2,nstack,jer,a0,b0,dy,dyp,ifaccept
            !write(13,"(4(I6,1X),4(D15.8,1X),I4)")  ind,nints2,nstack,jer,a0,b0,dy,dyp,ifaccept
    
    
            if (ifaccept .eq. 0) then
    
                if (nstack+2 .ge. maxstack) then
                    ier = 8
                    return
                endif
    
                nstack          = nstack+1
                stack(1,nstack) = (a0+b0)/2
                stack(2,nstack) = b0
    
                nstack          = nstack+1
                stack(1,nstack) = a0
                stack(2,nstack) = (a0+b0)/2
    
            else
    
                if (nints2+1 .ge. maxints) then
                    ier = 4
                    return
                endif
    
                nints2        = nints2+1
                ab2(1,nints2) = a0
                ab2(2,nints2) = b0
            endif
    
        end do
    
    
    endif
    
    !
    !  Copy out the solution
    !
    
    nints = nints1+nints2
    
    
    allocate(ab(2,1:nints))
    ab(:,1:nints1)       = ab1(:,nints1:1:-1)
    ab(:,nints1+1:nints) = ab2(:,1:nints2)
    call chebpw_specified(chebsol,k,nints,ab)
    allocate(ys(k*nints,n), yders(k*nints,n) )
    
    ys = 0 
    
    n1 = k*nints1
    n2 = k*nints
    
    do j=1,n
        ys(1:n1,j)       = reshape(ys1(1:k,j,nints1:1:-1), [k*nints1] )
        ys(n1+1:n2,j)    = reshape(ys2(1:k,j,1:nints2), [k*nints2] )
        yders(1:n1,j)    = reshape(yders1(1:k,j,nints1:1:-1), [k*nints1] )
        yders(n1+1:n2,j) = reshape(yders2(1:k,j,1:nints2), [k*nints2] )
    end do
    
    end subroutine
    
    
    subroutine odesolve_linearivp(vars,a,b,k,ts,n,as,fs,ys,yders)
    implicit double precision (a-h,o-z)
    type(odesolve_vars_t)        :: vars
    type(c_ptr)                  :: userptr
    double precision             :: ts(:)
    double complex               :: ys(:,:), yders(:,:), fs(:,:)
    double complex               :: as(:,:,:)
    !
    !  Solve an initial value problem for the linear problem
    !
    !    { y'(t) = A(t) y(t) + f(t),     a < t < b,
    !    { y(a)  = ya 
    !
    !  Input parameters:
    !    vars - the structure populated by the initialization routine
    !    (a,b) - the the interval over which the problem is given
    !    k - the number of Chebyshev nodes on the interval (a,b)
    !    ts - an array specifying the Chebyshev discretization nodes t_1, ..., t_k on the interval
    !       (a,b)
    !    n - the number of differential equations in the system
    !    as - a (k,n,n) matrix such that as(:,:,j) gives the values of the matrix
    !      A(t_j)
    !    fs - an (k,n) matrix such that fs(:,j) gives the values of the jth component
    !      of the vector-valued function f at the discretization nodes 
    !
    !    ys(1,:) - specifies the  vector ya
    !
    !  Output parameters:
    !    ys - an (k,n) array such that ys(i,j) gives the value y_j(t_i)
    !    yders - an (k,n) array such that yders(i,j) gives the value y_j'(t_i)
    
    !
    
    double complex, allocatable :: amatr0(:,:), rhs0(:), x0(:), y0(:)
    double complex, allocatable :: amatr(:,:), rhs(:), bmatr(:,:)
    
    k         = vars%kcheb
    
    allocate(amatr(n*k,n*k), rhs(n*k), bmatr(k,k) )
    
    i1 = 1
    i2 = k
    
    do l=1,n
    rhs(i1:i2) = ys(1,l) + (b-a)/2*matmul(vars%aintl,fs(:,l))
    i1         = i1+k
    i2         = i2+k
    end do
    
    i1 = 1
    i2 = k
    
    amatr = 0
    
    do i=1,n
    j1 = 1
    j2 = k
    do j=1,n
    
    bmatr              = -(b-a)/2 * vars%aintl
    do l=1,k
    bmatr(:,l)         = bmatr(:,l) * as(l,i,j)
    end do
    amatr(i1:i2,j1:j2) = bmatr
    
    j1 = j1+k
    j2 = j2+k
    
    end do
    
    i1 = i1+k
    i2 = i2+k
    end do
    
    do i=1,n*k
    amatr(i,i) = 1.0d0+amatr(i,i)
    end do
    
    call linalg0_solve(n*k,amatr,rhs)
    
    i1 = 1
    i2 = k
    
    do i=1,n
    ys(:,i) = rhs(i1:i2)
    i1      = i1+k
    i2      = i2+k
    end do
    
    
    do i=1,k
    yders(i,:) = matmul(as(i,:,:),ys(i,:)) + fs(i,:)
    end do
    
    end subroutine
    
    
    subroutine odesolve_lineartvp(vars,a,b,k,ts,n,as,fs,ys,yders)
    implicit double precision (a-h,o-z)
    type(odesolve_vars_t)        :: vars
    type(c_ptr)                  :: userptr
    double precision             :: ts(:)
    double complex               :: ys(:,:), yders(:,:), fs(:,:)
    double complex               :: as(:,:,:)
    !
    !  Solve the problem
    !
    !    { y'(t) = A(t) y(t) + f(t),     a < t < b,
    !    { y(b)  = yb 
    !
    !  Input parameters:
    !    vars - the structure populated by the initialization routine
    !    (a,b) - the the interval over which the problem is given
    !    k - the number of Chebyshev nodes on the interval (a,b)
    !    ts - an array specifying the Chebyshev discretization nodes t_1, ..., t_k on the interval
    !       (a,b)
    !    n - the number of differential equations in the system
    !    as - a (k,n,n) matrix such that as(:,:,j) gives the values of the matrix
    !      A(t_j)
    !    fs - an (k,n) matrix such that fs(:,j) gives the values of the jth component
    !      of the vector-valued function f at the discretization nodes 
    !
    !    ys(k,:) - specifies the  vector yb
    !
    !  Output parameters:
    !    ys - an (k,n) array such that ys(i,j) gives the value y_j(t_i)
    !    yders - an (k,n) array such that yders(i,j) gives the value y_j'(t_i)
    
    !
    
    double complex, allocatable :: amatr0(:,:), rhs0(:), x0(:), y0(:)
    double complex, allocatable :: amatr(:,:), rhs(:), bmatr(:,:)
    
    k         = vars%kcheb
    
    allocate(amatr(n*k,n*k), rhs(n*k), bmatr(k,k) )
    
    i1 = 1
    i2 = k
    
    do l=1,n
    rhs(i1:i2) = ys(k,l) + (b-a)/2*matmul(vars%aintr,fs(:,l))
    i1         = i1+k
    i2         = i2+k
    end do
    
    i1 = 1
    i2 = k
    
    amatr = 0
    
    do i=1,n
    j1 = 1
    j2 = k
    do j=1,n
    bmatr              = -(b-a)/2 * vars%aintr
    do l=1,k
    bmatr(:,l)         = bmatr(:,l) * as(l,i,j)
    end do
    amatr(i1:i2,j1:j2) = bmatr
    
    j1 = j1+k
    j2 = j2+k
    end do
    
    i1 = i1+k
    i2 = i2+k
    end do
    
    do i=1,n*k
    amatr(i,i) = 1.0d0+amatr(i,i)
    end do
    
    call linalg0_solve(n*k,amatr,rhs)
    
    i1 = 1
    i2 = k
    
    do i=1,n
    ys(:,i) = rhs(i1:i2)
    i1      = i1+k
    i2      = i2+k
    end do
    
    
    do i=1,k
    yders(i,:) = matmul(as(i,:,:),ys(i,:)) + fs(i,:)
    end do
    
    end subroutine
    
    
    subroutine odesolve_trapivp1(jer,eps,maxsteps,a,b,k,ts,n,odefun,userptr,ys,yders,sum_set, sum_coefs)
    implicit double precision (a-h,o-z)
    procedure(odesolve_fun1)     :: odefun
    type(c_ptr)                  :: userptr
    double complex               :: ys(:,:), yders(:,:)
    double precision             :: ts(:)
    !
    !  Use the trapezoidal rule to approximate the solution of the problem
    !
    !   {  y'(t) = f(t, y(t)),    a < t < b
    !   {  y(a)  = ya
    !
    !  Input parameters:
    !    eps - the desired precision for the computations
    !    (a,b) - the interval overwhich the problem is given
    !    k - the number of discretization nodes of the interval
    !    ts - the discretization points
    !    n - the number of differential equations in the system
    !    odefun - a user-supplied subroutine for evaluating the function
    !      f and its derivative with respect to y
    !    userptr - a "void *" pointer which is passed to the user-supplied subroutine
    !
    !    ys(1,:) - upon input, the vector ya giving the initial values of the
    !      components of the solution y
    !
    !  Output parameters:
    !    jer - an error return code;
    !      jer = 0   indicates successful execution
    !      jer = 4   indicates the maximum number of Newton steps was exceeded
    !      jer = 16  indicates that NaN or Inf was encountered
    !
    !    ys - an (k,n) array such that ys(i,j) gives the value of y_j(t_i)
    !    yders - an (k,n) array such that yders(i,j) gives the value of
    !      y_j'(t_i)
    !  
    
    double complex, allocatable   :: amatr(:,:), delta(:), rhs(:)
    double complex, allocatable   :: df0(:,:), df1(:,:)
    
    double precision              :: sum_coefs(:)
    integer                       :: sum_set(:, :)
    
    eps0 = epsilon(0.0d0)
    
    allocate(amatr(n,n), df0(n,n), df1(n,n))
    allocate(delta(n), rhs(n))
    
    jer     = 0
    
    do i=2,k
        t0  = ts(i-1)
        t1  = ts(i)
        h   = t1-t0
    
        call odefun(n,t0,ys(i-1,:),yders(i-1,:),df0,userptr,sum_set, sum_coefs)
        ys(i,:) = ys(i-1,:)
        ! + h * yders(i-1,:)
    
    
    
        do iter=1,maxsteps
    
            call odefun(n,t1,ys(i,:),yders(i,:),df1,userptr,sum_set, sum_coefs)
    
            amatr = -h/2*df1
            
            do j=1,n
                amatr(j,j) = amatr(j,j)+1
            end do
            delta = ys(i-1,:)-ys(i,:)+h/2 * (yders(i-1,:)+yders(i,:))
    
            
    
            call linalg0_solve(n,amatr,delta)
    
            ! rhs = delta
            ! call linalg0_qrsolve(1.0d-12,n,n,amatr,delta,rhs)
    
            if (eps0 .gt. 1.0d-16) then
                do j=1,k
                    ddd = real(delta(j))
                    if (ieee_is_nan(ddd) .or. .not. ieee_is_finite(ddd))  then
                        jer = 16
                        return
                    endif
                end do
            endif
    
            dd1 = norm2(abs(delta))
            dd2 = norm2(abs(ys(i,:)))
            if (dd2 .eq. 0) dd2 = 1
    
            ys(i,:) = ys(i,:) + delta
    
            if (dd1 .lt. eps * dd2) exit
    
    
        end do
    
    
        if (iter .gt. maxsteps) then
            jer = 4
            return
        endif
    
    end do
    
    
    end subroutine
    
    
    subroutine odesolve_traptvp1(jer,eps,maxsteps,a,b,k,ts,n,odefun,userptr,ys,yders,sum_set,sum_coefs)
    implicit double precision (a-h,o-z)
    procedure(odesolve_fun1)     :: odefun
    type(c_ptr)                  :: userptr
    double complex               :: ys(:,:), yders(:,:)
    double precision             :: ts(:)
    !
    !  Use the trapezoidal rule to approximate the solution of the problem
    !
    !   {  y'(t) = f(t, y(t)),    a < t < b
    !   {  y(b)  = yb
    !
    !  Input parameters:
    !    eps - the desired precision for the computations
    !    (a,b) - the interval overwhich the problem is given
    !    k - the number of discretization nodes of the interval
    !    ts - the discretization points
    !    n - the number of differential equations in the system
    !    odefun - a user-supplied subroutine for evaluating the function
    !      f and its derivative with respect to y
    !    userptr - a "void *" pointer which is passed to the user-supplied subroutine
    !
    !    ys(k,:) - upon input, the vector ya giving the terminal values of the
    !      components of the solution y
    !
    !  Output parameters:
    !    jer - an error return code;
    !      jer = 0   indicates successful execution
    !      jer = 4   indicates the maximum number of Newton steps was exceeded
    !      jer = 16  indicates that NaN or Inf was encountered
    !
    !    ys - an (k,n) array such that ys(i,j) gives the value of y_j(t_i)
    !    yders - an (k,n) array such that yders(i,j) gives the value of
    !      y_j'(t_i)
    !  
    
    double complex, allocatable   :: amatr(:,:), delta(:), rhs(:)
    double complex, allocatable   :: df0(:,:), df1(:,:)
    
    integer                       :: sum_set(:, :)
    double precision              :: sum_coefs(:)
    
    allocate(amatr(n,n), df0(n,n), df1(n,n))
    allocate(delta(n), rhs(n))
    eps0    = epsilon(0.0d0)
    jer     = 0
    
    do i=k-1,1,-1
    t0  = ts(i)
    t1  = ts(i+1)
    h   = t1-t0
    
    call odefun(n,t1,ys(i+1,:),yders(i+1,:),df1,userptr,sum_set, sum_coefs)
    ys(i,:) = ys(i+1,:)
    ! - h * yders(i+1,:)
    
    do iter=1,maxsteps
    
    call odefun(n,t0,ys(i,:),yders(i,:),df0,userptr,sum_set,sum_coefs)
    
    amatr = h/2*df0
    do j=1,n
    amatr(j,j) = amatr(j,j)+1
    end do
    delta = ys(i+1,:)-ys(i,:)-h/2 * (yders(i,:)+yders(i+1,:))
    
    call linalg0_solve_c(n,amatr,delta)
    
    if (eps0 .gt. 1.0d-16) then
    do j=1,k
    ddd = real(delta(j))
    if (ieee_is_nan(ddd) .or. .not. ieee_is_finite(ddd))  then
    jer = 16
    return
    endif
    end do
    endif
    
    dd1 = norm2(abs(delta))
    dd2 = norm2(abs(ys(i,:)))
    if (dd2 .eq. 0) dd2 = 1
    if (dd1 .lt. eps * dd2) exit
    
    ys(i,:) = ys(i,:) + delta
    
    end do
    
    if (iter .gt. maxsteps) then
    jer = 4
    return
    endif
    
    end do
    
    
    end subroutine
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  Second interface
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine odesolve_nonlinear2(vars,ier,ifit,eps,a,b,c,yc,n,odefun,pars,chebsol,ys,yders)
    implicit double precision (a-h,o-z)
    type(odesolve_vars_t)                       :: vars
    double complex, pointer                     :: pars(:)
    type(chebpw_scheme)                         :: chebsol
    procedure(odesolve_fun2)                    :: odefun
    double complex, allocatable, intent(out)    :: ys(:,:), yders(:,:)
    double complex                              :: yc(:)
    !
    !  Adaptively solve a problem of the form
    !
    !    { y'(t) = f( t, y(t) ),    a < t < b,
    !    {  y(c) = yc
    !
    !  where a <= c <= b and f: R x C^n \to C^n is a smooth function supplied
    !  via an external subroutine.   This routine uses a straightforward and rather slow spectral 
    !  method.  Each component of the solution is represented via a piecewise Chebyshev expansion.
    !
    !  Input parameters:
    !    vars - the structure populated by odesolve_init
    !
    !    ifit - an integer parameter which determines the criterion used for adaptive
    !      subdivision;
    !
    !      ifit = 1   means test every component of the solution for fit
    !      ifit = 2   means test the first component of the solution for fit (this is meant to 
    !                 be used when solving nth order scalar equations)
    !
    !    eps - the desired precision for the solution
    !    (a,b) - the interval over which the problem is given
    !    c - the point at which the conditions for y(t) are specified
    !    n - the number of equations in the system
    !    odesolve - an external subroutine conforming to the odesolve_fun2
    !     interface which supplies the 
    !    pars - array of doubles passed to odesolve
    !    yc - a vector giving the value of y at the point c
    !    
    !  Output parameters:
    !    chebsol - a structure specifying the piecewise discretization scheme
    !      used to represent the solutions
    !    ys - an array whose jth column gives the values of the jth component
    !      y_j of the solution at the discretization nodes
    !    yders - an array whose jth column gives the values of the derivative y_j'
    !     of the jth component of the solution at the discretization nodes     
    !
    
    double precision, allocatable       :: stack(:,:),  ts0(:), ab2(:,:), ab1(:,:), ab(:,:)
    double complex, allocatable         :: ys1(:,:,:), yders1(:,:,:)
    double complex, allocatable         :: ys2(:,:,:), yders2(:,:,:)
    double complex, allocatable         :: f0(:,:), df0(:,:,:), hs0(:,:), hders0(:,:)
    
    ier        = 0
    
    k          = vars%kcheb
    maxsteps   = vars%maxsteps
    maxiters   = vars%maxiters
    maxstack   = vars%maxstack
    
    allocate(stack(2,maxstack), ts0(k))
    allocate(ys1(k,n,maxstack), ab1(2,maxstack), yders1(k,n,maxstack) )
    allocate(ys2(k,n,maxstack), ab2(2,maxstack), yders2(k,n,maxstack) )
    allocate(f0(k,n), df0(k,n,n), hs0(k,n), hders0(k,n)  )
    
    !
    !  Solve backward from b to a
    !
    nints1     = 0
    
    if (a < c) then
    
    ind        = 0
    nstack     = 1
    stack(1,1) = a
    stack(2,1) = c
    
    do while (nstack .gt. 0)
    
    a0       = stack(1,nstack)
    b0       = stack(2,nstack)
    nstack   = nstack-1
    ts0      = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
    
    ifaccept = 0
    jer      = -1
    dy       = -1
    dyp      = -1
    
    if (nints1 .eq. 0 ) then
    ys1(k,:,nints1+1) = yc
    else
    ys1(k,:,nints1+1) = ys1(1,:,nints1)
    endif
    
    !
    !  Use the trapezoidal rule to build an initial approximation
    !
    
    call odesolve_traptvp2(jer,eps,maxsteps,a0,b0,k,ts0,n,odefun,pars,ys1(:,:,nints1+1),yders1(:,:,nints1+1))
    
    
    ! Use spectral integration to compute the y from y'
    
    ! do i=1,k
    ! call odefun(n,ts0(i),ys1(i,:,nints1+1),f0(i,:),df0(i,:,:),userptr)
    ! yders1(i,:,nints1+1) = f0(i,:)
    ! end do
    
    do j=1,n
    ys1(:,j,nints1+1) = ys1(k,j,nints1+1)+(b0-a0)/2 * matmul(vars%aintr,yders1(:,j,nints1+1))
    end do
    
    
    if (jer .eq. 0) then
    
    do iter=1,maxiters
    
    do i=1,k
    call odefun(n,ts0(i),ys1(i,:,nints1+1),f0(i,:),df0(i,:,:),pars)
    end do
    
    
    f0 = f0 - yders1(:,:,nints1+1)
    hs0      = 0
    call odesolve_lineartvp(vars,a0,b0,k,ts0,n,df0,f0,hs0,hders0)
    
    dd1 = norm2(abs(ys1(:,:,nints1+1)))
    if (dd1 .eq. 0) dd1 = 1
    dd2 = norm2(abs(hs0))
    
    ys1(:,:,nints1+1)    = ys1(:,:,nints1+1) + hs0
    yders1(:,:,nints1+1) = yders1(:,:,nints1+1) + hders0
    
    if (dd2 .lt. eps*dd1) exit
    
    end do
    
    call odesolve_fit(vars,ifit,ys1(:,:,nints1+1),dy)
    ! call odesolve_fit0(vars,ys1(:,:,nints1+1),dy)
    
    ! if (dy .lt. eps .AND. dyp .lt. eps) ifaccept = 1
    
    if (dy .lt. eps) ifaccept = 1
    
    
    endif
    
    ! ind = ind+1
    ! write(*, "(2(I6,1X),4(D15.8,1X),I4)")  ind,jer,a0,b0,dy,dyp,ifaccept
    ! write(13,"(2(I6,1X),4(D15.8,1X),I4)")  ind,jer,a0,b0,dy,dyp,ifaccept
    
    if (ifaccept .eq. 0) then
    
    if (nstack+2 .ge. maxstack) then
    ier = 8
    return
    endif
    
    
    nstack          = nstack+1
    stack(1,nstack) = a0
    stack(2,nstack) = (a0+b0)/2
    
    nstack          = nstack+1
    stack(1,nstack) = (a0+b0)/2
    stack(2,nstack) = b0
    
    else
    
    if (nints1+1 .ge. maxstack) then
    ier = 4
    return
    endif
    
    
    nints1        = nints1+1
    ab1(1,nints1) = a0
    ab1(2,nints1) = b0
    
    endif
    
    end do
    
    endif
    
    
    !
    !  Solve forward from c to b.
    !
    
    nints2     = 0
    
    if (c < b) then
    ind        = 0
    nstack     = 1
    stack(1,1) = c
    stack(2,1) = b
    
    
    do while (nstack .gt. 0)
    
    a0       = stack(1,nstack)
    b0       = stack(2,nstack)
    nstack   = nstack-1
    ts0      = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
    
    ifaccept = 0
    jer      = -1
    dy       = -1
    dyp      = -1
    
    ys2(:,:,nints2+1) = 0
    if (nints2 .eq. 0 ) then
    ys2(1,:,nints2+1) = yc
    else
    ys2(1,:,nints2+1) = ys2(k,:,nints2)
    endif
    
    !
    !  Use the trapezoidal rule to build an initial approximation
    !
    call odesolve_trapivp2(jer,eps,maxsteps,a0,b0,k,ts0,n,odefun,pars,ys2(:,:,nints2+1),yders2(:,:,nints2+1))
    
    
    ! Use spectral integration to compute the y from y'
    ! do i=1,k
    ! call odefun(n,ts0(i),ys2(i,:,nints2+1),f0(i,:),df0(i,:,:),userptr)
    ! yders2(i,:,nints2+1) = f0(i,:)
    ! end do
    
    do j=1,n
    ys2(:,j,nints2+1) = ys2(1,j,nints2+1)+(b0-a0)/2 * matmul(vars%aintl,yders2(:,j,nints2+1))
    end do
    
    
    if (jer .eq. 0) then
    
    do iter=1,maxiters
    
    do i=1,k
    call odefun(n,ts0(i),ys2(i,:,nints2+1),f0(i,:),df0(i,:,:),pars)
    end do
    
    f0       = f0 - yders2(:,:,nints2+1)
    hs0      = 0
    call odesolve_linearivp(vars,a0,b0,k,ts0,n,df0,f0,hs0,hders0)
    
    dd1 = norm2(abs(ys2(:,:,nints2+1)))
    if (dd1 .eq. 0) dd1 = 1
    dd2 = norm2(abs(hs0))
    
    ys2(:,:,nints2+1)    = ys2(:,:,nints2+1) + hs0
    yders2(:,:,nints2+1) = yders2(:,:,nints2+1) + hders0
    
    if (dd2 .lt. eps*dd1) exit
    
    end do
    
    call odesolve_fit(vars,ifit,ys2(:,:,nints2+1),dy)
    !call odesolve_fit0(vars,ys2(:,:,nints2+1),dy)
    ! if (dy .lt. eps .AND. dyp .lt. eps) ifaccept = 1
    
    
    if (dy .lt. eps) ifaccept = 1
    
    endif
    
    ! ind = ind+1
    ! write(*, "(2(I6,1X),4(D15.8,1X),I4)")  ind,jer,a0,b0,dy,dyp,ifaccept
    ! write(13,"(2(I6,1X),4(D15.8,1X),I4)")  ind,jer,a0,b0,dy,dyp,ifaccept
    
    
    if (ifaccept .eq. 0) then
    
    if (nstack+2 .ge. maxstack) then
    ier = 8
    return
    endif
    
    nstack          = nstack+1
    stack(1,nstack) = (a0+b0)/2
    stack(2,nstack) = b0
    
    nstack          = nstack+1
    stack(1,nstack) = a0
    stack(2,nstack) = (a0+b0)/2
    
    else
    
    if (nints2+1 .ge. maxstack) then
    ier = 4
    return
    endif
    
    
    nints2        = nints2+1
    ab2(1,nints2) = a0
    ab2(2,nints2) = b0
    
    endif
    
    end do
    
    
    endif
    
    
    
    !
    !  Copy out the solution
    !
    
    nints = nints1+nints2
    
    
    allocate(ab(2,1:nints))
    ab(:,1:nints1)       = ab1(:,nints1:1:-1)
    ab(:,nints1+1:nints) = ab2(:,1:nints2)
    call chebpw_specified(chebsol,k,nints,ab)
    allocate(ys(k*nints,n), yders(k*nints,n) )
    
    ys = 0 
    
    n1 = k*nints1
    n2 = k*nints
    
    do j=1,n
    
    ys(1:n1,j)       = reshape(ys1(1:k,j,nints1:1:-1), [k*nints1] )
    ys(n1+1:n2,j)    = reshape(ys2(1:k,j,1:nints2), [k*nints2] )
    yders(1:n1,j)    = reshape(yders1(1:k,j,nints1:1:-1), [k*nints1] )
    yders(n1+1:n2,j) = reshape(yders2(1:k,j,1:nints2), [k*nints2] )
    end do
    
    
    end subroutine
    
    
    
    subroutine odesolve_trapivp2(jer,eps,maxsteps,a,b,k,ts,n,odefun,pars,ys,yders)
    implicit double precision (a-h,o-z)
    procedure(odesolve_fun2)     :: odefun
    double complex, pointer      :: pars(:)
    double complex               :: ys(:,:), yders(:,:)
    double precision             :: ts(:)
    !
    !  Use the trapezoidal rule to approximate the solution of the problem
    !
    !   {  y'(t) = f(t, y(t)),    a < t < b
    !   {  y(a)  = ya
    !
    !  Input parameters:
    !    eps - the desired precision for the computations
    !    (a,b) - the interval overwhich the problem is given
    !    k - the number of discretization nodes of the interval
    !    ts - the discretization points
    !    n - the number of differential equations in the system
    !    odefun - a user-supplied subroutine for evaluating the function
    !      f and its derivative with respect to y
    !    userptr - a "void *" pointer which is passed to the user-supplied subroutine
    !
    !    ys(1,:) - upon input, the vector ya giving the initial values of the
    !      components of the solution y
    !
    !  Output parameters:
    !    jer - an error return code;
    !      jer = 0   indicates successful execution
    !      jer = 4   indicates the maximum number of Newton steps was exceeded
    !      jer = 16  indicates that NaN or Inf was encountered
    !
    !    ys - an (k,n) array such that ys(i,j) gives the value of y_j(t_i)
    !    yders - an (k,n) array such that yders(i,j) gives the value of
    !      y_j'(t_i)
    !  
    
    double complex, allocatable   :: amatr(:,:), delta(:), rhs(:)
    double complex, allocatable   :: df0(:,:), df1(:,:)
    
    allocate(amatr(n,n), df0(n,n), df1(n,n))
    allocate(delta(n), rhs(n))
    eps0    = epsilon(0.0d0)
    jer     = 0
    
    do i=2,k
    t0  = ts(i-1)
    t1  = ts(i)
    h   = t1-t0
    
    call odefun(n,t0,ys(i-1,:),yders(i-1,:),df0,pars)
    ys(i,:) = ys(i-1,:)
    ! + h * yders(i-1,:)
    
    do iter=1,maxsteps
    
    call odefun(n,t1,ys(i,:),yders(i,:),df1,pars)
    
    amatr = -h/2*df1
    do j=1,n
    amatr(j,j) = amatr(j,j)+1
    end do
    delta = ys(i-1,:)-ys(i,:)+h/2 * (yders(i-1,:)+yders(i,:))
    
    call linalg0_solve(n,amatr,delta)
    
    if (eps0 .gt. 1.0d-16) then
    do j=1,k
    ddd = real(delta(j))
    if (ieee_is_nan(ddd) .or. .not. ieee_is_finite(ddd))  then
    jer = 16
    return
    endif
    end do
    endif
    
    dd1 = norm2(abs(delta))
    dd2 = norm2(abs(ys(i,:)))
    if (dd2 .eq. 0) dd2 = 1
    if (dd1 .lt. eps * dd2) exit
    !if (dd1 .lt. eps) exit
    
    
    ys(i,:) = ys(i,:) + delta
    
    end do
    
    if (iter .gt. maxsteps) then
    jer = 4
    return
    endif
    
    end do
    
    
    end subroutine
    
    
    subroutine odesolve_traptvp2(jer,eps,maxsteps,a,b,k,ts,n,odefun,pars,ys,yders)
    implicit double precision (a-h,o-z)
    procedure(odesolve_fun2)     :: odefun
    double complex, pointer      :: pars(:)
    double complex               :: ys(:,:), yders(:,:)
    double precision             :: ts(:)
    !
    !  Use the trapezoidal rule to approximate the solution of the problem
    !
    !   {  y'(t) = f(t, y(t)),    a < t < b
    !   {  y(b)  = yb
    !
    !  Input parameters:
    !    eps - the desired precision for the computations
    !    (a,b) - the interval overwhich the problem is given
    !    k - the number of discretization nodes of the interval
    !    ts - the discretization points
    !    n - the number of differential equations in the system
    !    odefun - a user-supplied subroutine for evaluating the function
    !      f and its derivative with respect to y
    !    userptr - a "void *" pointer which is passed to the user-supplied subroutine
    !
    !    ys(k,:) - upon input, the vector ya giving the terminal values of the
    !      components of the solution y
    !
    !  Output parameters:
    !    jer - an error return code;
    !      jer = 0   indicates successful execution
    !      jer = 4   indicates the maximum number of Newton steps was exceeded
    !      jer = 16  indicates that NaN or Inf was encountered
    !
    !    ys - an (k,n) array such that ys(i,j) gives the value of y_j(t_i)
    !    yders - an (k,n) array such that yders(i,j) gives the value of
    !      y_j'(t_i)
    !  
    
    double complex, allocatable   :: amatr(:,:), delta(:), rhs(:)
    double complex, allocatable   :: df0(:,:), df1(:,:)
    
    allocate(amatr(n,n), df0(n,n), df1(n,n))
    allocate(delta(n), rhs(n))
    
    eps0    = epsilon(0.0d0)
    jer     = 0
    
    do i=k-1,1,-1
    t0  = ts(i)
    t1  = ts(i+1)
    h   = t1-t0
    
    call odefun(n,t1,ys(i+1,:),yders(i+1,:),df1,pars)
    ys(i,:) = ys(i+1,:)
    ! - h * yders(i+1,:)
    
    do iter=1,maxsteps
    
    call odefun(n,t0,ys(i,:),yders(i,:),df0,pars)
    
    amatr = h/2*df0
    do j=1,n
    amatr(j,j) = amatr(j,j)+1
    end do
    delta = ys(i+1,:)-ys(i,:)-h/2 * (yders(i,:)+yders(i+1,:))
    
    call linalg0_solve_c(n,amatr,delta)
    
    if (eps0 .gt. 1.0d-16) then
    do j=1,k
    ddd = real(delta(j))
    if (ieee_is_nan(ddd) .or. .not. ieee_is_finite(ddd))  then
    jer = 16
    return
    endif
    end do
    endif
    
    
    dd1 = norm2(abs(delta))
    dd2 = norm2(abs(ys(i,:)))
    if (dd2 .eq. 0) dd2 = 1
    if (dd1 .lt. eps * dd2) exit
    
    ys(i,:) = ys(i,:) + delta
    
    
    end do
    
    if (iter .gt. maxsteps) then
    jer = 4
    return
    endif
    
    end do
    
    end subroutine
    
    
    end module
    