!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains a highly-robust but fairly slow solver for second order ordinary
!  differential equations.  It operates via an adaptive integration method and solutions
!  are represented as piecewise Chebyshev expansions. 
!
!  The following routines should be regarded as public:
! 
!    odetwo_init - initialize the solver by populating a structure which is passed to
!      the other subroutines
!
!    odetwo_nonlinear - adaptively solve a problem of the form
!
!        {  y'(t) = f(t,y(t))         a <= t <= b
!        {  y(c)  = yc
!        {  y'(c) = ypc
! 
!      where  a <= c <= b and f:R x C^2 --> C^2 is a smooth function supplied by 
!      an external subroutine
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module odetwo

use utils
use linalg0
use chebyshev
use chebpw
use iso_c_binding
use ieee_arithmetic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The structure which stores any data needed by the procedures in this modules and
!  which is populated by the initialization routine.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type       odetwo_vars_t

integer                                    :: kcheb
integer                                    :: ntest
integer                                    :: maxstack
integer                                    :: maxints
integer                                    :: maxiters
integer                                    :: maxsteps

double precision, allocatable              :: xscheb(:)
double precision, allocatable              :: whtscheb(:)
double precision, allocatable              :: acoefs(:,:)
double precision, allocatable              :: aintl(:,:)
double precision, allocatable              :: aintr(:,:)
double precision, allocatable              :: aintl2(:,:)
double precision, allocatable              :: aintr2(:,:)

end type   odetwo_vars_t


interface


subroutine odetwo_fun1(t,y,yp,f,dfdy,dfdyp,userptr)
use iso_c_binding
implicit double precision (a-h,o-z)
double complex              :: y, yp, f, dfdy, df, dfdyp
type(c_ptr)                 :: userptr
end subroutine

subroutine odetwo_fun2(t,y,yp,f,dfdy,dfdyp,pars)
use iso_c_binding
implicit double precision (a-h,o-z)
double complex              :: y, yp, f, dfdy, df, dfdyp
double complex, pointer     :: pars(:)
end subroutine

end interface

interface          odetwo_nonlinear
module procedure   odetwo_nonlinear1
module procedure   odetwo_nonlinear2
end interface      odetwo_nonlinear

contains


subroutine odetwo_init(vars,kcheb0,ntest0,maxstack0,maxints0,maxsteps0,maxiters0)
implicit double precision (a-h,o-z)
type(odetwo_vars_t)                         :: vars
integer, optional                           :: kcheb0, ntest0, maxstack0, maxsteps0
integer, optional                           :: maxints0, maxiters0
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
!    maxints0 - the maximum number of intervals for the discretization scheme used
!      to represent solutions
!    maxsteps0 - the maximum number of iterations for the trapezoid rule
!    maxiters0 - the maximum number of iterations for Newton's method
!
!  Output parameters:
!    vars - the structure containing all of the data needed by the other procedures in
!      this module
!

if (.not. present(kcheb0) ) then
kcheb    = 16
ntest    = 2
maxstack = 10000
maxints  = 10000
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

allocate(vars%aintl2(kcheb,kcheb))
allocate(vars%aintr2(kcheb,kcheb))

vars%aintl2 = matmul(vars%aintl,vars%aintl)
vars%aintr2 = matmul(vars%aintr,vars%aintr)

end subroutine


subroutine odetwo_fit(vars,vals,dcoefs)
implicit double precision (a-h,o-z)
type(odetwo_vars_t)         :: vars
double complex              :: vals(:), coefs(vars%kcheb)

kcheb  = vars%kcheb
ntest  = vars%ntest
coefs  = matmul(vars%acoefs,vals)
dd1    = norm2(abs(coefs))
dd2    = norm2(abs(coefs(kcheb-ntest+1:kcheb)))
if (dd1 .eq. 0) dd1=1
dcoefs = dd2/dd1

end subroutine



subroutine odetwo_nonlinear1(vars,ier,eps,a,b,c,yc,ypc,chebsol,ys,yders,yder2s,odefun,userptr)
implicit double precision (a-h,o-z)
type(odetwo_vars_t)                                 :: vars
type(chebpw_scheme), intent(out)                    :: chebsol
double complex, allocatable, intent(out)            :: ys(:), yders(:), yder2s(:)
double complex                                      :: yc, ypc
procedure(odetwo_fun1)                              :: odefun
type(c_ptr)                                         :: userptr
!  
!  Adaptively solve the problem
!
!      {  y'(t) = f(t,y(t)),         a <= t <= b,
!      {  y(c)  = yc
!      {  y'(c) = ypc
! 
!  where a <= c <= b and f:R x C^2 --> C^2  is a smooth function suppplied
!  via an external subroutine.  The solution and its first two derivatives
!  are represented via piecewise Chebyshev expansions
!
!  Input parameters:
!    vars - the structure populated by odetwo_init
!    eps - the desired precision for the calculations
!    (a,b) - the interval over which the problem is given
!    c - the point c at which the values of y and y' are given
!    yc - the desired value of y(c)
!    ypc - the desired value of y'(c)
!
!  Output parameters:
!    chebsol - a structure specifying the piecewise discretization scheme used to
!      represent the solutions
!    ys - the values of the obtained solution at the discretization nodes
!    yders - the values of the derivative of the obtained solution at the discretization 
!       nodes
!    yder2s - the values of the second derivative of the obtained solution at the 
!      discretization nodes
!

double precision, allocatable                       :: stack(:,:), ab1(:,:), ab2(:,:), ts0(:), ab(:,:)
double complex, allocatable                         :: ps0(:), qs0(:), fs0(:)
double complex, allocatable                         :: hs(:), hders(:), hder2s(:)
double complex, allocatable                         :: ys1(:,:), yders1(:,:), yder2s1(:,:)
double complex, allocatable                         :: ys2(:,:), yders2(:,:), yder2s2(:,:)

double complex                                      :: f, dfdy, dfdyp, ima

ier     = 0
epstrap = eps
epsnewt = eps
ima     = (0.0d0,1.0d0)

!
!  Read the parameters from the vars structure
! 

k             = vars%kcheb
ntest         = vars%ntest
maxstack      = vars%maxstack
maxints       = vars%maxints
maxsteps      = vars%maxsteps
maxiters      = vars%maxiters

!
!  Allocate temporary arrays
!


allocate(stack(2,maxstack))
 
allocate(ts0(k),ps0(k),qs0(k),fs0(k))
allocate(hs(k),hders(k),hder2s(k))

allocate(ab1(2,maxints), ys1(k,maxints), yders1(k,maxints), yder2s1(k,maxints) )
allocate(ab2(2,maxints), ys2(k,maxints), yders2(k,maxints), yder2s2(k,maxints) )

!
!  Solve going backward from the point c
!


nints1 = 0
ind    = 0

if (a < c) then

nstack     = 1
stack(1,1) = a
stack(2,1) = c

do while (nstack > 0 )

a0         = stack(1,nstack)
b0         = stack(2,nstack)
nstack     = nstack-1

ts0        = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
ifaccept   = 0
dcoefs     = -1
dy         = -1
dy0        = -1
dyp        = -1

! Evaluate the coefficients 

!
!  Use the implicit trapezoidal method to construct an approximation
!  to use as an initial guess for Newton's method.
!

if (nints1 .eq. 0) then
ys1(k,nints1+1)    = yc
yders1(k,nints1+1) = ypc
else
ys1(k,nints1+1)    = ys1(1,nints1)
yders1(k,nints1+1) = yders1(1,nints1)
endif


! Build an initial solution using the trapezoid rule
call odetwo_traptvp1(jer,epstrap,maxsteps,k,ts0,odefun,ys1(:,nints1+1),yders1(:,nints1+1),yder2s1(:,nints1+1),userptr)

if (jer .eq. 0) then

!  The values of y, y' and y'' calculated using the implicit trapezoidal method satisfy
!  the differential equation at the specified nodes, but they are not consistent with eachother.  
!  That is, y' is not the derivative of y in any reasonable sense.
!
!  We integrate the second derivative twice in order to produce consistent
!  approximations.
!
!  Skipping this step is fatal.
!

yders1(:,nints1+1)  = matmul(vars%aintr*(b0-a0)/2,yder2s1(:,nints1+1)) + yders1(k,nints1+1)
ys1(:,nints1+1)     = matmul(vars%aintr*(b0-a0)/2,yders1(:,nints1+1))  + ys1(k,nints1+1)

!
!  Perform Newton iterations.
!

do inewt = 1, maxiters

!
!  Form the coefficients for the linearized problem and compute
!  the error in the current solution.
!
do i=1,k
call odefun(ts0(i),ys1(i,nints1+1),yders1(i,nints1+1),f,dfdy,dfdyp,userptr)
ps0(i)    = -dfdyp
qs0(i)    = -dfdy
fs0(i)    = f-yder2s1(i,nints1+1)
end do

!
!  Call the code_linear_ivp routine to solve the linearized system.
!

hs(k)    = 0
hders(k) = 0
call odetwo_lineartvp1(vars,a0,b0,k,ps0,qs0,fs0,hs,hders,hder2s)


if ( norm2(abs(hs)) .lt. eps*norm2(abs(ys1(:,nints1+1)))) exit

!
!  Update the solution
!

ys1(:,nints1+1)     = ys1(:,nints1+1) + hs
yders1(:,nints1+1)  = yders1(:,nints1+1) + hders
yder2s1(:,nints1+1) = yder2s1(:,nints1+1) + hder2s

end do


 call odetwo_fit(vars,ys1(:,nints1+1),dy)
!call odetwo_fit(vars,yders1(:,nints1+1),dyp)
if (dy .lt. eps ) ifaccept = 1
endif

! ind = ind+1
! write (*,"(2(I6),3(D20.10,2X),1X,I1)")    ind,jer,a0,b0,dyp,ifaccept
! write (13,"(2(I6),3(D20.10,2X),1X,I1)")   ind,jer,a0,b0,dyp,ifaccept

if (ifaccept .eq. 1) then

if (nints1+1 .gt. maxints) then
ier    = 4
return
endif

nints1        = nints1+1
ab1(1,nints1) = a0
ab1(2,nints1) = b0

else

if (nstack+2 .gt. maxstack) then
ier = 8
return
endif
 
nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

endif


end do

endif

!
!  Solve going forward from the point c
!

nints2 = 0
ind    = 0

if (b > c) then

nstack     = 1
stack(1,1) = c
stack(2,1) = b

do while (nstack > 0 )

a0         = stack(1,nstack)
b0         = stack(2,nstack)
nstack     = nstack-1

ts0        = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
ifaccept   = 0
dcoefs     = -1
dy         = -1
dy0        = -1
dyp        = -1

! Evaluate the coefficients 

!
!  Use the implicit trapezoidal method to construct an approximation
!  to use as an initial guess for Newton's method.
!

if (nints2 .eq. 0) then
ys2(1,nints2+1)    = yc
yders2(1,nints2+1) = ypc
else
ys2(1,nints2+1)    = ys2(k,nints2)
yders2(1,nints2+1) = yders2(k,nints2)
endif



! Build an initial solution using the trapezoid rule
call odetwo_trapivp1(jer,epstrap,maxsteps,k,ts0,odefun,ys2(:,nints2+1),yders2(:,nints2+1),yder2s2(:,nints2+1),userptr)

if (jer .eq. 0) then

!  The values of y, y' and y'' calculated using the implicit trapezoidal method satisfy
!  the differential equation at the specified nodes, but they are not consistent with eachother.  
!  That is, y' is not the derivative of y in any reasonable sense.
!
!  We integrate the second derivative twice in order to produce consistent
!  approximations.
!
!  Skipping this step is fatal.
!

yders2(:,nints2+1)  = matmul(vars%aintl*(b0-a0)/2,yder2s2(:,nints2+1)) + yders2(1,nints2+1)
ys2(:,nints2+1)     = matmul(vars%aintl*(b0-a0)/2,yders2(:,nints2+1))  + ys2(1,nints2+1)


!
!  Perform Newton iterations.
!

do inewt = 1, maxiters

!
!  Form the coefficients for the linearized problem and compute
!  the error in the current solution.
!
do i=1,k
call odefun(ts0(i),ys2(i,nints2+1),yders2(i,nints2+1),f,dfdy,dfdyp,userptr)
ps0(i)    = -dfdyp
qs0(i)    = -dfdy
fs0(i)    = f-yder2s2(i,nints2+1)
end do

!
!  Call the code_linear_ivp routine to solve the linearized system.
!


hs(1)    = 0
hders(1) = 0
call odetwo_linearivp1(vars,a0,b0,k,ps0,qs0,fs0,hs,hders,hder2s)

if ( norm2(abs(hs)) .lt. eps*norm2(abs(ys2(:,nints2+1)))) exit

!
!  Update the solution
!

ys2(:,nints2+1)     = ys2(:,nints2+1) + hs
yders2(:,nints2+1)  = yders2(:,nints2+1) + hders
yder2s2(:,nints2+1) = yder2s2(:,nints2+1) + hder2s

end do

call odetwo_fit(vars,ys2(:,nints2+1),dy)
!call odetwo_fit(vars,yders2(:,nints2+1),dyp)

if (dy .lt. eps) ifaccept = 1

endif


! ind = ind+1
! write (*,"(2(I6),3(D20.10,2X),1X,I1)")    ind,jer,a0,b0,dyp,ifaccept
! write (13,"(2(I6),3(D20.10,2X),1X,I1)")   ind,jer,a0,b0,dyp,ifaccept

if (ifaccept .eq. 1) then

if (nints2+1 .gt. maxints) then
ier    = 4
return
endif

nints2        = nints2+1
ab2(1,nints2) = a0
ab2(2,nints2) = b0

else

if (nstack+2 .gt. maxstack) then
ier = 8
return
endif
 
nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2
endif


end do

endif


!
!  Form the piecewise Chebyshev discretization scheme and copy out the solution
!

nints = nints1+nints2

allocate(ab(2,1:nints))
ab(:,1:nints1)       = ab1(:,nints1:1:-1)
ab(:,nints1+1:nints) = ab2(:,1:nints2)

!call prin2("ab = ",ab)
call chebpw_specified(chebsol,k,nints,ab)

allocate(ys(k*nints), yders(k*nints), yder2s(k*nints) )

n1 = k*nints1
n2 = k*nints

ys(1:n1)         = reshape(ys1(1:k,nints1:1:-1), [k*nints1] )
ys(n1+1:n2)      = reshape(ys2(1:k,1:nints2), [k*nints2] )

yders(1:n1)      = reshape(yders1(1:k,nints1:1:-1), [k*nints1] )
yders(n1+1:n2)   = reshape(yders2(1:k,1:nints2), [k*nints2] )

yder2s(1:n1)     = reshape(yder2s1(1:k,nints1:1:-1), [k*nints1] )
yder2s(n1+1:n2)  = reshape(yder2s2(1:k,1:nints2), [k*nints2] )

end subroutine




subroutine odetwo_nonlinear2(vars,ier,eps,a,b,c,yc,ypc,chebsol,ys,yders,yder2s,odefun,pars)
implicit double precision (a-h,o-z)
type(odetwo_vars_t)                                 :: vars
type(chebpw_scheme), intent(out)                    :: chebsol
double complex, allocatable, intent(out)            :: ys(:), yders(:), yder2s(:)
double complex                                      :: yc, ypc
procedure(odetwo_fun2)                              :: odefun
double complex, pointer                             :: pars(:)
!  
!  Adaptively solve the problem
!
!      {  y'(t) = f(t,y(t)),         a <= t <= b,
!      {  y(c)  = yc
!      {  y'(c) = ypc
! 
!  where a <= c <= b and f:R x C^2 --> C^2  is a smooth function suppplied
!  via an external subroutine.  The solution and its first two derivatives
!  are represented via piecewise Chebyshev expansions
!
!  Input parameters:
!    vars - the structure populated by odetwo_init
!    eps - the desired precision for the calculations
!    (a,b) - the interval over which the problem is given
!    c - the point c at which the values of y and y' are given
!    yc - the desired value of y(c)
!    ypc - the desired value of y'(c)
!
!  Output parameters:
!    chebsol - a structure specifying the piecewise discretization scheme used to
!      represent the solutions
!    ys - the values of the obtained solution at the discretization nodes
!    yders - the values of the derivative of the obtained solution at the discretization 
!       nodes
!    yder2s - the values of the second derivative of the obtained solution at the 
!      discretization nodes
!

double precision, allocatable                       :: stack(:,:), ab1(:,:), ab2(:,:), ts0(:), ab(:,:)
double complex, allocatable                         :: ps0(:), qs0(:), fs0(:)
double complex, allocatable                         :: hs(:), hders(:), hder2s(:)
double complex, allocatable                         :: ys1(:,:), yders1(:,:), yder2s1(:,:)
double complex, allocatable                         :: ys2(:,:), yders2(:,:), yder2s2(:,:)

double complex                                      :: f, dfdy, dfdyp, ima

ier     = 0
epstrap = eps
epsnewt = eps
ima     = (0.0d0,1.0d0)

!
!  Read the parameters from the vars structure
! 

k             = vars%kcheb
ntest         = vars%ntest
maxstack      = vars%maxstack
maxints       = vars%maxints
maxsteps      = vars%maxsteps
maxiters      = vars%maxiters

!
!  Allocate temporary arrays
!


allocate(stack(2,maxstack))
 
allocate(ts0(k),ps0(k),qs0(k),fs0(k))
allocate(hs(k),hders(k),hder2s(k))

allocate(ab1(2,maxints), ys1(k,maxints), yders1(k,maxints), yder2s1(k,maxints) )
allocate(ab2(2,maxints), ys2(k,maxints), yders2(k,maxints), yder2s2(k,maxints) )

!
!  Solve going backward from the point c
!


nints1 = 0
ind    = 0

if (a < c) then

nstack     = 1
stack(1,1) = a
stack(2,1) = c

do while (nstack > 0 )

a0         = stack(1,nstack)
b0         = stack(2,nstack)
nstack     = nstack-1

ts0        = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
ifaccept   = 0
dcoefs     = -1
dy         = -1
dy0        = -1
dyp        = -1

! Evaluate the coefficients 

!
!  Use the implicit trapezoidal method to construct an approximation
!  to use as an initial guess for Newton's method.
!

if (nints1 .eq. 0) then
ys1(k,nints1+1)    = yc
yders1(k,nints1+1) = ypc
else
ys1(k,nints1+1)    = ys1(1,nints1)
yders1(k,nints1+1) = yders1(1,nints1)
endif


! Build an initial solution using the trapezoid rule
call odetwo_traptvp2(jer,epstrap,maxsteps,k,ts0,odefun,ys1(:,nints1+1),yders1(:,nints1+1),yder2s1(:,nints1+1),pars)

if (jer .eq. 0) then

!  The values of y, y' and y'' calculated using the implicit trapezoidal method satisfy
!  the differential equation at the specified nodes, but they are not consistent with eachother.  
!  That is, y' is not the derivative of y in any reasonable sense.
!
!  We integrate the second derivative twice in order to produce consistent
!  approximations.
!
!  Skipping this step is fatal.
!

yders1(:,nints1+1)  = matmul(vars%aintr*(b0-a0)/2,yder2s1(:,nints1+1)) + yders1(k,nints1+1)
ys1(:,nints1+1)     = matmul(vars%aintr*(b0-a0)/2,yders1(:,nints1+1))  + ys1(k,nints1+1)

!
!  Perform Newton iterations.
!

do inewt = 1, maxiters

!
!  Form the coefficients for the linearized problem and compute
!  the error in the current solution.
!
do i=1,k
call odefun(ts0(i),ys1(i,nints1+1),yders1(i,nints1+1),f,dfdy,dfdyp,pars)
ps0(i)    = -dfdyp
qs0(i)    = -dfdy
fs0(i)    = f-yder2s1(i,nints1+1)
end do

!
!  Call the code_linear_ivp routine to solve the linearized system.
!

hs(k)    = 0
hders(k) = 0
call odetwo_lineartvp1(vars,a0,b0,k,ps0,qs0,fs0,hs,hders,hder2s)


if ( norm2(abs(hs)) .lt. eps*norm2(abs(ys1(:,nints1+1)))) exit

!
!  Update the solution
!

ys1(:,nints1+1)     = ys1(:,nints1+1) + hs
yders1(:,nints1+1)  = yders1(:,nints1+1) + hders
yder2s1(:,nints1+1) = yder2s1(:,nints1+1) + hder2s

end do


 call odetwo_fit(vars,ys1(:,nints1+1),dy)
!call odetwo_fit(vars,yders1(:,nints1+1),dyp)
if (dy .lt. eps ) ifaccept = 1
endif

! ind = ind+1
! write (*,"(2(I6),3(D20.10,2X),1X,I1)")    ind,jer,a0,b0,dyp,ifaccept
! write (13,"(2(I6),3(D20.10,2X),1X,I1)")   ind,jer,a0,b0,dyp,ifaccept

if (ifaccept .eq. 1) then

if (nints1+1 .gt. maxints) then
ier    = 4
return
endif

nints1        = nints1+1
ab1(1,nints1) = a0
ab1(2,nints1) = b0

else

if (nstack+2 .gt. maxstack) then
ier = 8
return
endif
 
nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

endif


end do

endif

!
!  Solve going forward from the point c
!

nints2 = 0
ind    = 0

if (b > c) then

nstack     = 1
stack(1,1) = c
stack(2,1) = b

do while (nstack > 0 )

a0         = stack(1,nstack)
b0         = stack(2,nstack)
nstack     = nstack-1

ts0        = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
ifaccept   = 0
dcoefs     = -1
dy         = -1
dy0        = -1
dyp        = -1

! Evaluate the coefficients 

!
!  Use the implicit trapezoidal method to construct an approximation
!  to use as an initial guess for Newton's method.
!

if (nints2 .eq. 0) then
ys2(1,nints2+1)    = yc
yders2(1,nints2+1) = ypc
else
ys2(1,nints2+1)    = ys2(k,nints2)
yders2(1,nints2+1) = yders2(k,nints2)
endif



! Build an initial solution using the trapezoid rule
call odetwo_trapivp2(jer,epstrap,maxsteps,k,ts0,odefun,ys2(:,nints2+1),yders2(:,nints2+1),yder2s2(:,nints2+1),pars)

if (jer .eq. 0) then

!  The values of y, y' and y'' calculated using the implicit trapezoidal method satisfy
!  the differential equation at the specified nodes, but they are not consistent with eachother.  
!  That is, y' is not the derivative of y in any reasonable sense.
!
!  We integrate the second derivative twice in order to produce consistent
!  approximations.
!
!  Skipping this step is fatal.
!

yders2(:,nints2+1)  = matmul(vars%aintl*(b0-a0)/2,yder2s2(:,nints2+1)) + yders2(1,nints2+1)
ys2(:,nints2+1)     = matmul(vars%aintl*(b0-a0)/2,yders2(:,nints2+1))  + ys2(1,nints2+1)


!
!  Perform Newton iterations.
!

do inewt = 1, maxiters

!
!  Form the coefficients for the linearized problem and compute
!  the error in the current solution.
!
do i=1,k
call odefun(ts0(i),ys2(i,nints2+1),yders2(i,nints2+1),f,dfdy,dfdyp,pars)
ps0(i)    = -dfdyp
qs0(i)    = -dfdy
fs0(i)    = f-yder2s2(i,nints2+1)
end do

!
!  Call the code_linear_ivp routine to solve the linearized system.
!


hs(1)    = 0
hders(1) = 0
call odetwo_linearivp1(vars,a0,b0,k,ps0,qs0,fs0,hs,hders,hder2s)

if ( norm2(abs(hs)) .lt. eps*norm2(abs(ys2(:,nints2+1)))) exit

!
!  Update the solution
!

ys2(:,nints2+1)     = ys2(:,nints2+1) + hs
yders2(:,nints2+1)  = yders2(:,nints2+1) + hders
yder2s2(:,nints2+1) = yder2s2(:,nints2+1) + hder2s

end do

call odetwo_fit(vars,ys2(:,nints2+1),dy)
!call odetwo_fit(vars,yders2(:,nints2+1),dyp)

if (dy .lt. eps) ifaccept = 1

endif


! ind = ind+1
! write (*,"(2(I6),3(D20.10,2X),1X,I1)")    ind,jer,a0,b0,dyp,ifaccept
! write (13,"(2(I6),3(D20.10,2X),1X,I1)")   ind,jer,a0,b0,dyp,ifaccept

if (ifaccept .eq. 1) then

if (nints2+1 .gt. maxints) then
ier    = 4
return
endif

nints2        = nints2+1
ab2(1,nints2) = a0
ab2(2,nints2) = b0

else

if (nstack+2 .gt. maxstack) then
ier = 8
return
endif
 
nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2
endif


end do

endif


!
!  Form the piecewise Chebyshev discretization scheme and copy out the solution
!

nints = nints1+nints2

allocate(ab(2,1:nints))
ab(:,1:nints1)       = ab1(:,nints1:1:-1)
ab(:,nints1+1:nints) = ab2(:,1:nints2)

!call prin2("ab = ",ab)
call chebpw_specified(chebsol,k,nints,ab)

allocate(ys(k*nints), yders(k*nints), yder2s(k*nints) )

n1 = k*nints1
n2 = k*nints

ys(1:n1)         = reshape(ys1(1:k,nints1:1:-1), [k*nints1] )
ys(n1+1:n2)      = reshape(ys2(1:k,1:nints2), [k*nints2] )

yders(1:n1)      = reshape(yders1(1:k,nints1:1:-1), [k*nints1] )
yders(n1+1:n2)   = reshape(yders2(1:k,1:nints2), [k*nints2] )

yder2s(1:n1)     = reshape(yder2s1(1:k,nints1:1:-1), [k*nints1] )
yder2s(n1+1:n2)  = reshape(yder2s2(1:k,1:nints2), [k*nints2] )

end subroutine




subroutine odetwo_linearivp1(vars,a,b,k,ps,qs,fs,ys,yders,yder2s)
implicit double precision (a-h,o-z)
type(odetwo_vars_t)           :: vars
integer, intent(in)           :: k
double complex, intent(in)    :: ps(k),qs(k),fs(k)
double complex, intent(out)   :: ys(k),yders(k),yder2s(k)
!
!  Solve an initial value for the ordinary differential equation
!
!     y''(t) + p(t) y'(t) + q(t) y(t) = f(t)                                              (3)
!
!  on the interval [a,b] using a standard Chebyshev spectral method.
!
!  Input parameters:
!    (a,b) - the interval on which the ODE (3) is given
!    k - the number of Chebyshev points on the interval [a,b]
!    ps - an array specifying the values of the function p(t) appearing in (3)
!      at the k Chebyshev nodes on [a,b]
!    qs - an array specifying the values of the function q(t) appearing in (3)
!      at the k Chebyshev node on [a,b]
!    fs - an array speciying the values of the function f(t) appearing in (3)
!      at the k Chebyshev nodes on [a,b]
!
!    ys(1) - the value of y(a)
!    yders(1) - the value of y'(a)
!
!  Output parameters:
!
!    ys - the values of the solution y of (3) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yders - the values of the solution y' of (3) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yder2s - the values of the solution y'' of (3) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!

double complex, allocatable :: amatr(:,:),xs(:),sigma(:),rhs(:)
double complex              :: alpha,beta

!
!  Allocate memory for the procedure and setup some parameters.
!

allocate(amatr(k,k),xs(k),sigma(k),rhs(k))


xs       = max(a,min(b,(b-a)/2 *vars%xscheb + (b+a)/2))
alpha    = ys(1)
beta     = yders(1)

!
!  We represent the solution in the form
!
!      y(t) = alpha + beta (t-a) + \int_a^t (t-s) sigma(s) ds,
!
!  insert this representation into (1), and solve the resulting system of
!  linear equations in order to obtain the values of sigma.  
!    

amatr  = 0
do i=1,k
amatr(i,i) = 1.0d0
sigma(i)   = 0.0d0
end do

!
! Handle the p(t) * y'(t) term.
!
do i=1,k
amatr(i,:) = amatr(i,:) + ps(i) * vars%aintl(i,:)*(b-a)/2
sigma(i)   = sigma(i) - ps(i)*beta
end do

!
!  Handle the q(t) y(t) term
!

do i=1,k
amatr(i,:) = amatr(i,:) + qs(i) * vars%aintl2(i,:)*((b-a)/2)**2
sigma(i)   = sigma(i) - qs(i) * (alpha + beta*(xs(i)-a))
end do


!
!  Form the right-hand side.
!
do i=1,k
sigma(i) = sigma(i) + fs(i)
end do

!
!  Use a QR decomposition to invert the linear system
!

call linalg0_solve_c(k,amatr,sigma)


!
!  Calculate y(t) and y'(t) from sigma.
!

yder2s = sigma
yders  = (b-a)/2*matmul(vars%aintl,sigma)
ys     = ((b-a)/2)**2*matmul(vars%aintl2,sigma)

do i=1,k
ys(i)     = ys(i) + alpha + beta*(xs(i)-a)
yders(i)  = yders(i) + beta
end do

end subroutine



subroutine odetwo_lineartvp1(vars,a,b,k,ps,qs,fs,ys,yders,yder2s)
implicit double precision (a-h,o-z)
type(odetwo_vars_t)           :: vars
integer, intent(in)           :: k
double complex, intent(in)    :: ps(k),qs(k),fs(k)
double complex, intent(out)   :: ys(k),yders(k),yder2s(k)

!
!  Solve a terminal value for the ordinary differential equation
!
!     y''(t) + p(t) y'(t) + q(t) y(t) = f(t)                                              (4)
!
!  on the interval [a,b] using a standard Chebyshev spectral method.
!
!  Input parameters:
!
!    (a,b) - the interval on which the ODE (4) is given
!    k - the number of Chebyshev points on the interval [a,b]
!    ps - an array specifying the values of the function p(t) appearing in (4)
!      at the k Chebyshev nodes on [a,b]
!    qs - an array specifying the values of the function q(t) appearing in (4)
!      at the k Chebyshev node on [a,b]
!    fs - an array speciying the values of the function f(t) appearing in (4)
!      at the k Chebyshev nodes on [a,b]
!
!    ys(k) - the value of y(b)
!    yders(k) - the value of y'(b)
!
!  Output parameters:
!
!    ys - the values of the solution y of (4) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yders - the values of the solution y' of (4) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yder2s - the values of the solution y'' of (4) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!

double complex, allocatable :: amatr(:,:),xs(:),sigma(:),rhs(:)
double complex              :: alpha,beta
integer, allocatable        :: ipiv(:)

!
!  Allocate memory for the procedure and setup some parameters.
!
allocate(amatr(k,k),xs(k),sigma(k))

xs       = max(a,min(b,(b-a)/2 * vars%xscheb + (b+a)/2))
alpha    = ys(k)
beta     = yders(k)

!
!  We represent the solution in the form
!
!      y(t) = alpha + beta (t-b) + \int_b^t (t-s) sigma(s) ds,
!
!  insert this representation into (2), and solve the resulting system of
!  linear equations in order to obtain the values of sigma.  
!    

amatr  = 0
do i=1,k
amatr(i,i) = 1.0d0
sigma(i)   = 0.0d0
end do

!
! Handle the p(t) * y'(t) term.
!
do i=1,k
amatr(i,:) = amatr(i,:) + ps(i) * vars%aintr(i,:)*(b-a)/2
sigma(i)   = sigma(i) - ps(i)*beta
end do

!
!  Handle the q(t) y(t) term
!

do i=1,k
amatr(i,:) = amatr(i,:) + qs(i) * vars%aintr2(i,:)*((b-a)/2)**2
sigma(i)   = sigma(i) - qs(i) * (alpha + beta*(xs(i)-b))
end do


!
!  Form the right-hand side.
!
do i=1,k
sigma(i) = sigma(i) + fs(i)
end do

!
!  Use a QR decomposition to invert the linear system
!

call linalg0_solve_c(k,amatr,sigma)


!
!  Calculate y(t) and y'(t) from sigma.
!

yder2s = sigma
yders  = (b-a)/2*matmul(vars%aintr,sigma)
ys     = ((b-a)/2)**2*matmul(vars%aintr2,sigma)

do i=1,k
ys(i)     = ys(i) + alpha + beta*(xs(i)-b)
yders(i)  = yders(i) + beta
end do

end subroutine



subroutine odetwo_trapivp1(ier,eps,maxsteps,k,ts,odefun,ys,yders,yder2s,userptr)
implicit double precision (a-h,o-z)
integer, intent(out)          :: ier
integer, intent(in)           :: k
double precision, intent(in)  :: ts(k)
double complex, intent(out)   :: ys(k),yders(k),yder2s(k)
procedure (odetwo_fun1)       :: odefun
type(c_ptr)                   :: userptr

!
!  Use the implicit trapezodial method to crudely approximate the solution of an
!  initial value problem for the nonlinear ordinary differential equation
!
!    y''(t) = f(t,y(t),y'(t))                                                             (5)
! 
!  on the interval [a,b].  The user specifies the nodes on which the solution
!  of is to be computed.
!
!  Input parameters:
!
!    eps - precision for the calculation
!    maxsteps - the maximum number of Newton steps
!    k - the number of nodes at which the solution of (5) is to be computed
!    ts - an array of length k which supplies the a sorted list of nodes at which 
!     (4) is to be solved
!    odefun - the user-supplied external routine conforming to the interface odetwo_fun1
!
!    ys(1) - the first entry of this array is the initial value for the solution y
!    yders(1) - the first entry of this array is the initial value for the solution y'
!    userptr - a "void *" pointer which is passed as an argument to the user-specified
!     external function odefun
!
!  Output parameters:
!
!    ier - an error return code;
!      ier = 0    indicates successful execution
!      ier = 16   means that the Newton procedure did not converge
!      ier = 64   means that NaN or Inf was encountered
!
!    ys - the values of the obtained approximation of the solution of (5) at
!     the specified nodes
!    yders - the values of the obtained approximation of the derivative of the solution of 
!     (5) at the specified nodes
!    yder2s - the values of the obtained approximation of the second derivative of the solution 
!     of  (5) at the specified nodes
!  

double complex :: y0,yp0,ypp0,ypp1,yp1,y1,dd,dd0,dfdy,dfdyp,val,der,delta

ier      = 0
! eps0     = epsilon(0.0d0)
! eps      = eps0*10

!
!  Evaluate the second derivative at the left-hand endpoint of the interval.
!

call odefun(ts(1),ys(1),yders(1),yder2s(1),dfdy,dfdyp,userptr)


do i=2,k
t0 = ts(i-1)
t  = ts(i)
h  = t-t0

y0   = ys(i-1)
yp0  = yders(i-1)
ypp0 = yder2s(i-1)

!
!  Set the initial guess.
!

!yp1 = yp0 + h *ypp0
yp1 = yp0 
y1  = y0 + h/2*(yp0+yp1)
call odefun(t,y1,yp1,ypp1,dfdy,dfdyp,userptr)
dd    = yp1 - yp0 - h/2*(ypp0+ypp1)

!
!  Conduct Newton iterations in an attempt to improve the siutation.
!

do iter=1,maxsteps

!
!  Record the current approximation.
!
! dd0       = dd
! ys(i)     = y1
! yders(i)  = yp1
! yder2s(i) = ypp1

!
!  Take a Newton step.
!

val   = dd
der   = 1.0d0-h/2*(dfdy*h/2+dfdyp)
delta = val/der

yp1   = yp1-delta
y1    = y0 + h/2*(yp0+yp1)
call odefun(t,y1,yp1,ypp1,dfdy,dfdyp,userptr)
dd    = yp1 - yp0 - h/2*(ypp0+ypp1)


ddd = abs(yp1)
if (ddd .eq. 0) ddd=1
if (abs(delta) .lt. eps * ddd)  exit

end do


if (ieee_is_nan(real(y1)) .OR. ieee_is_nan(imag(y1)) ) then
ier = 64
return
endif

if (.not. ieee_is_finite(real(y1)) .OR. .not. ieee_is_finite(imag(y1))) then
ier = 64
return
endif

ys(i)     = y1
yders(i)  = yp1
yder2s(i) = ypp1


if (iter .gt. maxsteps) then
ier = 16
return
endif

end do




end subroutine


subroutine odetwo_traptvp1(ier,eps,maxsteps,k,ts,odefun,ys,yders,yder2s,userptr)
implicit double precision (a-h,o-z)
integer, intent(out)           :: ier
integer, intent(in)            :: k
double precision, intent(in)   :: ts(k)
double complex, intent(out)    :: ys(k),yders(k),yder2s(k)
procedure (odetwo_fun1)        :: odefun
type(c_ptr)                    :: userptr
!
!  Use the implicit trapezoidal method to crudely approximate the solution of a
!  terminal value problem for the nonlinear ordinary differential equation
!
!    y''(t) = f(t,y(t),y'(t))                                                             (5)
! 
!  on the interval [a,b].  The user specifies the nodes on which the solution
!  of is to be computed.
!
!  Input parameters:
!
!    k - the number of nodes at which the solution of (5) is to be computed
!    ts - an array of length k which supplies the a sorted list of nodes at which 
!     (5) is to be solved
!    odefun - a user-supplied external procedure of type "codefunction" which
!     supplied the values of the function f as well as the derivatives of
!     f w.r.t. y' and y  (see the definition of codefunction above)
!
!    ys(k) - the last entry of this array is the terminal value for the solution y
!    yders(k) - the first entry of this array is the terminal value for the solution y'
!    userptr - a "void *" pointer which is passed as an argument to the user-specified
!     external function odefun
!
!  Output parameters:
!
!    ier - an error return code;
!       ier = 0    indicates successful execution
!       ier = 1024 means NaN or Infinity was encountered; this usually means
!                  that a spectacular overflow took place
!
!    ys - the values of the obtained approximation of the solution of (5) at
!     the specified nodes
!    yders - the values of the obtained approximation of the derivative of the solution of 
!     (5) at the specified nodes
!    yder2s - the values of the obtained approximation of the second derivative of the solution 
!     of  (5) at the specified nodes
!  

double complex :: y0,yp0,ypp0,ypp1,yp1,y1,dd,dd0,dfdy,dfdyp,val,der,delta

ier        = 0

call odefun(ts(k),ys(k),yders(k),yder2s(k),dfdy,dfdyp,userptr)

do i=k-1,1,-1
t  = ts(i)
h  = ts(i+1)-ts(i)

y1   = ys(i+1)
yp1  = yders(i+1)
ypp1 = yder2s(i+1)

!
!  Set the initial guess.
!

yp0 = yp1
y0  = y1 - h/2*(yp0+yp1)

call odefun(t,y0,yp0,ypp0,dfdy,dfdyp,userptr)
dd  = yp1 - yp0 - h/2*(ypp0+ypp1)

!
!  Try to improve it via Newton iterations
!

do iter=1,maxsteps

!
!  Record the current guess
!

! dd0       = dd
! ys(i)     = y0
! yders(i)  = yp0
! yder2s(i) = ypp0

!
!  Take a Newton step
!
val   = dd
der   = -1.0d0-h/2*(-dfdy*h/2+dfdyp)
delta = val/der

yp0   = yp0 -delta
y0    = y1 - h/2*(yp0+yp1)
call odefun(t,y0,yp0,ypp0,dfdy,dfdyp,userptr)
dd  = yp1 - yp0 - h/2*(ypp0+ypp1)


ddd = abs(yp0)
if (ddd .eq. 0) ddd=1
if (abs(delta) .lt. eps * ddd)  exit

end do



if (iter .gt. maxsteps) then
ier = 16
return
endif

if (ieee_is_nan(real(y0)) .OR. ieee_is_nan(imag(y0)) ) then
ier = 64
return
endif

if (.not. ieee_is_finite(real(y0)) .OR. .not. ieee_is_finite(imag(y0))) then
ier = 64
return
endif

ys(i)     = y0
yders(i)  = yp0
yder2s(i) = ypp0

end do

end subroutine




subroutine odetwo_trapivp2(ier,eps,maxsteps,k,ts,odefun,ys,yders,yder2s,pars)
implicit double precision (a-h,o-z)
integer, intent(out)          :: ier
integer, intent(in)           :: k
double precision, intent(in)  :: ts(k)
double complex, intent(out)   :: ys(k),yders(k),yder2s(k)
procedure (odetwo_fun2)       :: odefun
double complex, pointer       :: pars(:)

!
!  Use the implicit trapezodial method to crudely approximate the solution of an
!  initial value problem for the nonlinear ordinary differential equation
!
!    y''(t) = f(t,y(t),y'(t))                                                             (5)
! 
!  on the interval [a,b].  The user specifies the nodes on which the solution
!  of is to be computed.
!
!  Input parameters:
!
!    eps - precision for the calculation
!    maxsteps - the maximum number of Newton steps
!    k - the number of nodes at which the solution of (5) is to be computed
!    ts - an array of length k which supplies the a sorted list of nodes at which 
!     (4) is to be solved
!    odefun - the user-supplied external routine conforming to the interface odetwo_fun1
!
!    ys(1) - the first entry of this array is the initial value for the solution y
!    yders(1) - the first entry of this array is the initial value for the solution y'
!    userptr - a "void *" pointer which is passed as an argument to the user-specified
!     external function odefun
!
!  Output parameters:
!
!    ier - an error return code;
!      ier = 0    indicates successful execution
!      ier = 16   means that the Newton procedure did not converge
!      ier = 64   means that NaN or Inf was encountered
!
!    ys - the values of the obtained approximation of the solution of (5) at
!     the specified nodes
!    yders - the values of the obtained approximation of the derivative of the solution of 
!     (5) at the specified nodes
!    yder2s - the values of the obtained approximation of the second derivative of the solution 
!     of  (5) at the specified nodes
!  

double complex :: y0,yp0,ypp0,ypp1,yp1,y1,dd,dd0,dfdy,dfdyp,val,der,delta

ier      = 0
! eps0     = epsilon(0.0d0)
! eps      = eps0*10

!
!  Evaluate the second derivative at the left-hand endpoint of the interval.
!

call odefun(ts(1),ys(1),yders(1),yder2s(1),dfdy,dfdyp,pars)


do i=2,k
t0 = ts(i-1)
t  = ts(i)
h  = t-t0

y0   = ys(i-1)
yp0  = yders(i-1)
ypp0 = yder2s(i-1)

!
!  Set the initial guess.
!

!yp1 = yp0 + h *ypp0
yp1 = yp0 
y1  = y0 + h/2*(yp0+yp1)
call odefun(t,y1,yp1,ypp1,dfdy,dfdyp,pars)
dd    = yp1 - yp0 - h/2*(ypp0+ypp1)

!
!  Conduct Newton iterations in an attempt to improve the siutation.
!

do iter=1,maxsteps

!
!  Record the current approximation.
!
! dd0       = dd
! ys(i)     = y1
! yders(i)  = yp1
! yder2s(i) = ypp1

!
!  Take a Newton step.
!

val   = dd
der   = 1.0d0-h/2*(dfdy*h/2+dfdyp)
delta = val/der

yp1   = yp1-delta
y1    = y0 + h/2*(yp0+yp1)
call odefun(t,y1,yp1,ypp1,dfdy,dfdyp,pars)
dd    = yp1 - yp0 - h/2*(ypp0+ypp1)


ddd = abs(yp1)
if (ddd .eq. 0) ddd=1
if (abs(delta) .lt. eps * ddd)  exit

end do


if (ieee_is_nan(real(y1)) .OR. ieee_is_nan(imag(y1)) ) then
ier = 64
return
endif

if (.not. ieee_is_finite(real(y1)) .OR. .not. ieee_is_finite(imag(y1))) then
ier = 64
return
endif

ys(i)     = y1
yders(i)  = yp1
yder2s(i) = ypp1


if (iter .gt. maxsteps) then
ier = 16
return
endif

end do




end subroutine


subroutine odetwo_traptvp2(ier,eps,maxsteps,k,ts,odefun,ys,yders,yder2s,pars)
implicit double precision (a-h,o-z)
integer, intent(out)           :: ier
integer, intent(in)            :: k
double precision, intent(in)   :: ts(k)
double complex, intent(out)    :: ys(k),yders(k),yder2s(k)
procedure (odetwo_fun2)        :: odefun
double complex, pointer        :: pars(:)
!
!  Use the implicit trapezoidal method to crudely approximate the solution of a
!  terminal value problem for the nonlinear ordinary differential equation
!
!    y''(t) = f(t,y(t),y'(t))                                                             (5)
! 
!  on the interval [a,b].  The user specifies the nodes on which the solution
!  of is to be computed.
!
!  Input parameters:
!
!    k - the number of nodes at which the solution of (5) is to be computed
!    ts - an array of length k which supplies the a sorted list of nodes at which 
!     (5) is to be solved
!    odefun - a user-supplied external procedure of type "codefunction" which
!     supplied the values of the function f as well as the derivatives of
!     f w.r.t. y' and y  (see the definition of codefunction above)
!
!    ys(k) - the last entry of this array is the terminal value for the solution y
!    yders(k) - the first entry of this array is the terminal value for the solution y'
!    userptr - a "void *" pointer which is passed as an argument to the user-specified
!     external function odefun
!
!  Output parameters:
!
!    ier - an error return code;
!       ier = 0    indicates successful execution
!       ier = 1024 means NaN or Infinity was encountered; this usually means
!                  that a spectacular overflow took place
!
!    ys - the values of the obtained approximation of the solution of (5) at
!     the specified nodes
!    yders - the values of the obtained approximation of the derivative of the solution of 
!     (5) at the specified nodes
!    yder2s - the values of the obtained approximation of the second derivative of the solution 
!     of  (5) at the specified nodes
!  

double complex :: y0,yp0,ypp0,ypp1,yp1,y1,dd,dd0,dfdy,dfdyp,val,der,delta

ier        = 0

call odefun(ts(k),ys(k),yders(k),yder2s(k),dfdy,dfdyp,pars)

do i=k-1,1,-1
t  = ts(i)
h  = ts(i+1)-ts(i)

y1   = ys(i+1)
yp1  = yders(i+1)
ypp1 = yder2s(i+1)

!
!  Set the initial guess.
!

yp0 = yp1
y0  = y1 - h/2*(yp0+yp1)

call odefun(t,y0,yp0,ypp0,dfdy,dfdyp,pars)
dd  = yp1 - yp0 - h/2*(ypp0+ypp1)

!
!  Try to improve it via Newton iterations
!

do iter=1,maxsteps

!
!  Record the current guess
!

! dd0       = dd
! ys(i)     = y0
! yders(i)  = yp0
! yder2s(i) = ypp0

!
!  Take a Newton step
!
val   = dd
der   = -1.0d0-h/2*(-dfdy*h/2+dfdyp)
delta = val/der

yp0   = yp0 -delta
y0    = y1 - h/2*(yp0+yp1)
call odefun(t,y0,yp0,ypp0,dfdy,dfdyp,pars)
dd  = yp1 - yp0 - h/2*(ypp0+ypp1)


ddd = abs(yp0)
if (ddd .eq. 0) ddd=1
if (abs(delta) .lt. eps * ddd)  exit

end do



if (iter .gt. maxsteps) then
ier = 16
return
endif

if (ieee_is_nan(real(y0)) .OR. ieee_is_nan(imag(y0)) ) then
ier = 64
return
endif

if (.not. ieee_is_finite(real(y0)) .OR. .not. ieee_is_finite(imag(y0))) then
ier = 64
return
endif

ys(i)     = y0
yders(i)  = yp0
yder2s(i) = ypp0

end do

end subroutine


end module
