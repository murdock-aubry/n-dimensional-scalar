!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for solving inital and boundary value problems for
!  inhomogeneous nth order scalar ordinary differential equations of the form
!
!    y^(n)(t) + q_{n-1} y^(n-1)(t) + ... + q_1 y'(t)  + q_0(t) y(t) = 0                      (1)
!
!  over an interval I.  It operates by constructing a basis in the space of solutions 
!  of (1) consisting of functions of the form
!
!    y_j(t) = exp(r_j(t))                                                                    (2)
!
!  with the r_j(t) slowly-varying.  Solutions of this type exist assuming the
!  coefficients q_j(t) are smooth and slowly-varying, and the eigenvalues of the 
!  coefficient matrix corresponding to (1) are distinct for each t in I.
!
!  Once the functions r_j(t) have been constructed, any reasonable boundary value 
!  problem for (1) can then be easily solved. 
!
!  For the time-being, the code only handles equations of orders 2, 3 and 4.  Code
!  for general n will eventually be written.
!
!  The following subroutines should be regarded as publicly callable:
!
!    scalar_init - initialize the solver by populating a data structure which is then
!      passed to the other subroutines in this module
!
!    scalar_homotope - compute the derivatives of the r_j(t) at a point t using the 
!      homotopy method
!
!    scalar_levin - compute the derivative of the r_j(t) at a point t using the
!      Levin method
!
!    scalar_riccati - given the values of one of the r_j(t) at a point, solve the 
!      Riccati equation associated with (1) in order to calculate it over the
!      whole interval
!
!    scalar_levinode - use the Levin method to compute a collection of nonoscilatory
!      solutions of the Riccati equation which represent a basis in the space of
!      solutions of (1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module scalar

use iso_c_binding
use utils
use linalg0
use chebyshev
use chebpw
use odesolve
use odetwo


type       scalar_vars_t
type(odesolve_vars_t)                              :: odesolve_vars
type(odetwo_vars_t)                                :: odetwo_vars

integer                                            :: kcheb
integer                                            :: ntest
integer                                            :: maxleviters
integer                                            :: maxstack
integer                                            :: maxints

double precision, allocatable                      :: xscheb(:)
double precision, allocatable                      :: whtscheb(:)
double precision, allocatable                      :: acoefs(:,:)
double precision, allocatable                      :: aintl(:,:)

double precision, allocatable                      :: adiff(:,:)
double precision, allocatable                      :: adiff2(:,:)
double precision, allocatable                      :: adiff3(:,:)
double precision, allocatable                      :: adiff4(:,:)
double precision, allocatable                      :: adiff5(:,:)

double precision, allocatable                      :: aintr(:,:)

end type   scalar_vars_t



interface

subroutine      scalar_fun(n,t,qs,pars)
implicit double precision (a-h,o-z)
double precision        :: t
double complex          :: qs(:)
double complex, pointer :: pars(:)
end subroutine

end interface


type       scalar_odedata
integer                                            :: n
procedure(scalar_fun), pointer, nopass             :: fun
double complex, pointer                            :: pars(:)
double complex, allocatable                        :: dks(:)
integer                                            :: ifleft
double precision                                   :: wa, wb
end type   scalar_odedata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  There are for the version of 3rd order code which uses the scalar ODE solver
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

interface 

subroutine      scalar_fun33(t,q0,q1,q2,par1,par2,par3)
implicit double precision (a-h,o-z)
double precision        :: t
double complex          :: q0, q1, q2
double complex          :: par1, par2, par3
!double complex, pointer :: pars(:)
!  y'''(t) + q2(t) y''(t) + q1(t) y'(t) + q0(t) y(t) = 0
end subroutine
end interface

type       scalar_odedata33
procedure(scalar_fun33), pointer, nopass           :: fun
double complex                                     :: par1, par2, par3
double complex                                     :: dk0, dk1, dk2
integer                                            :: ifleft
double precision                                   :: wa, wb
end type   scalar_odedata33

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

subroutine scalar_init(vars,kcheb0,ntest0,maxstack0,maxints0,maxsteps0,maxiters0,&
  maxleviters0)
implicit double precision (a-h,o-z)
type(scalar_vars_t)                         :: vars
integer, optional                           :: kcheb0, ntest0, maxstack0, maxsteps0
integer, optional                           :: maxints0, maxiters0, maxleviters0
!
!  Initialize the structure containing any data needed by the other procedures 
!  in this module.
!
!  If vars is the only parameter passed to this routine, then reasonable defaults
!  are chosen.
!
!  Input parameters:
!    kcheb0 - the number of terms in the local Chebyshev expansions used to represent
!      solutions of ODEs
!    ntest0 - the number of trailing terms in the local Chebyshev expansions to 
!      check when testing the fitness of a purported solution
!    maxstack0 - the length of the stack used by the adaptive ODE solver
!    maxsteps0 - the maximum number of iterations for the trapezoid rule used by the ODE 
!      solver
!    maxiters0 - the maximum number of iterations for Newton procedure used by the ODE 
!      solver
!    maxlevinter0 - the maximum number of iterations for Levin's method
!
!  Output parameters:
!    vars - the structure containing all of the data needed by the other procedures in
!      this module
!

if (.not. present(kcheb0) ) then
kcheb       = 16
ntest       = 2
maxstack    = 1000
maxints     = 10000
maxsteps    = 8
maxiters    = 8
maxleviters = 30
else
kcheb       = kcheb0
ntest       = ntest0
maxstack    = maxstack0
maxints     = maxints0
maxsteps    = maxsteps0
maxiters    = maxiters0
maxleviters = maxleviters0
endif

call odesolve_init(vars%odesolve_vars,kcheb,ntest,maxstack,maxints,maxsteps,maxiters)
call odetwo_init(vars%odetwo_vars,kcheb,ntest,maxstack,maxints,maxsteps,maxiters)

vars%kcheb       = kcheb
vars%ntest       = ntest
vars%maxleviters = maxleviters
vars%maxstack    = maxstack
vars%maxints     = maxints

call chebyshev_quad(kcheb,vars%xscheb,vars%whtscheb)
call chebyshev_coefsmatrix(kcheb,vars%acoefs)
call chebyshev_intlmatrix(kcheb,vars%aintl)
call chebyshev_intrmatrix(kcheb,vars%aintr)
call chebyshev_diffmatrix(kcheb,vars%adiff)

allocate(vars%adiff2(kcheb,kcheb))
allocate(vars%adiff3(kcheb,kcheb))
allocate(vars%adiff4(kcheb,kcheb))
allocate(vars%adiff5(kcheb,kcheb))

vars%adiff2 = matmul(vars%adiff,vars%adiff)
vars%adiff3 = matmul(vars%adiff2,vars%adiff)
vars%adiff4 = matmul(vars%adiff3,vars%adiff)
vars%adiff5 = matmul(vars%adiff4,vars%adiff)

end subroutine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Order 2 using the standard ODE solver
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine scalar_homotope2(vars,ier,eps,a,b,ifleft,fun,pars,rvals)
implicit double precision (a-h,o-z)
type(scalar_vars_t)              :: vars
procedure(scalar_fun)            :: fun
double complex, pointer          :: pars(:)
double complex                   :: rvals(:,:)
!
!  Use the homotopy method to compute the derivatives of the r_j(t) at a specified 
!  point.  More explicitly, this routine homotopes the coefficient matrix for the user's 
!  nth order scalar system to the constant system
!
!          (  0     1   )
!  A0    = (  0     0   )
!          ( -dk0 -dk1  )
!
!  in order to calculate appropriate initial conditions for the user's system.
!  Reasonable values of dk? are determined by calling the user-supplied external
!  subroutine.
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    eps - precision for the calculations
!    (a,b) - the interval on which to perform the procedure
!    ifleft - an integer parameter indicating whether the values of the r_j(t) should
!      be computed at the left-hand endpoint of the interval a or the right-hand
!      endpoint b;
!
!        ifleft = 1   means calculate the values at a
!        ifleft = 0   means calculate the values at b
!
!    fun - an external subroutine supplying the coefficients of the user's scalar
!      equation
!    par? - user-supplied parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0    indicates successful execution
!       ier = 4    indicates the ODE solver failed
!       ier = 256  indicates that the eigensolver failed
!
!    rvals - the (i,j) entry contains the value of r_j^(i-1)(c), where c=a if
!      ifleft =1  and c=b if ifleft=0
!    

type(chebpw_scheme)               :: chebsol
type(c_ptr)                       :: userptr
type(scalar_odedata), pointer     :: odedata
double complex, allocatable       :: ys(:,:), yders(:,:)
double complex                    :: ima, yc(1), y0, y0p, w(2), dks(2)
double precision                  :: rcoefs(2), icoefs(2), reigs(2), ieigs(2)
integer                           :: its(2), flag

ier  = 0

k    = vars%kcheb
eps0 = epsilon(0.0d0)


wa   = 11d0/(b-a)
wb   = (a+b)/2.0d0
ima  = (0.0d0,1.0d0)

n    = 2

call fun(n,(a+b)/2,dks,pars)

allocate(odedata)
allocate(odedata%dks(n))

odedata%n       = n-1
odedata%pars   => pars
odedata%dks     = dks
odedata%fun    => fun
odedata%wa      = wa
odedata%wb      = wb
odedata%ifleft  = ifleft
userptr         = c_loc(odedata)


if (ifleft .eq. 0) then
c            = a
else
c            = b
endif

!
!  Find the roots of the polynomial
!

its        = 200

rcoefs(1) = real(dks(2))
icoefs(1) = imag(dks(2))

rcoefs(2) = real(dks(1))
icoefs(2) = imag(dks(1))

call init_random_seed()
call zamvw(n,rcoefs,icoefs,reigs,ieigs,its,flag)

if (flag .ne. 0) then
ier = 256
return
endif

w(1) = reigs(1) + ima*ieigs(1)
w(2) = reigs(2) + ima*ieigs(2)

if (norm2(abs(w)) .eq. 0) then
ier = 256
return
endif

!call scalar_sortw(2,w)

do i=1,n
yc(1)      = w(i)
ifit       = 2

call odesolve_nonlinear1(vars%odesolve_vars,jer,ifit,eps,a,b,c,yc,n-1, &
 scalar_ode2,userptr,chebsol,ys,yders)

if (jer .ne. 0) then
ier = 4
return
endif

nn = chebsol%nints*k
if (ifleft .eq. 1) then
rvals(1,i)  = ys(1,1)
else
rvals(1,i)  = ys(nn,1)
endif


end do

end subroutine


subroutine scalar_levin2(vars,ier,eps,a,b,ifleft,fun,pars,rvals)
   implicit double precision (a-h,o-z)
   type(scalar_vars_t)      :: vars
   procedure(scalar_fun)    :: fun
   double complex, pointer  :: pars(:)
   double complex           :: rvals(:,:)
   !
   !  Use the Levin method to compute the derivatives of the r_j(t) at a specified 
   !  point.
   !
   !  Input parameters:
   !    vars - the structure prepared by scalar_init
   !    eps - precision for the Newton iteration
   !    (a,b) - the interval on which to perform the procedure
   !    ifleft - an integer parameter indicating whether the values of the r_j(t) should
   !      be computed at the left-hand endpoint of the interval a or the right-hand
   !      endpoint b;
   !
   !        ifleft = 1   means calculate the values at a
   !        ifleft = 0   means calculate the values at b
   !
   !    fun - an external subroutine supplying the coefficients of the 
   !      user's scalar system 
   !    par? - user-supplied parameters which are passed to fun
   !
   !  Output parameters:
   !    ier - an error return code;
   !       ier = 0     indicates successful execution
   !       ier = 4     indicates the ODE solver failed
   !       ier = 256   means that the eigensolver failed
   !       ier = 1024  indicates the coefficients were not represented resolved
   !                   by Chebyshev expansions of order vars%kcheb over [a,b]
   !
   !    rvals - the (i,j) entry contains the value of r_j^(i)(c), where c=a if 
   !      ifleft =1 and c=b if ifleft=0
   !

   double precision, allocatable            :: ts(:)
   double precision                         :: rcoefs(2), icoefs(2), reigs(2), ieigs(2)
   integer                                  :: its(2), flag

   double complex, allocatable              :: q0s(:), q1s(:), q2s(:), qs(:), coefs(:), ws(:,:)

   double complex, allocatable              :: r(:), rp(:), rpp(:), rhs(:), amatr(:,:), delta(:)
   double complex                           :: ima, rval, rder, w(2)
   double complex                           :: rp0, rpp0

   double complex, allocatable             :: roots(:)


   

   kcheb    = vars%kcheb
   ntest    = vars%ntest
   maxiters = vars%maxleviters

   ier      = 0
   ima      = (0.0d0,1.0d0)
   eps0     = epsilon(0.0d0)
   epssol   = eps0*10

   

   n        = 2
   allocate( ts(kcheb), q0s(kcheb), q1s(kcheb), qs(n), coefs(kcheb)  )
   allocate( r(kcheb), rp(kcheb), rpp(kcheb), rhs(kcheb), amatr(kcheb,kcheb), delta(kcheb) )

   ts = (b-a)/2*vars%xscheb + (b+a)/2

   do i=1,kcheb
      call fun(n,ts(i),qs,pars)
      q0s(i) = qs(1)
      q1s(i) = qs(2)
   end do


   

   ! 
   !  Test the coefficinets for fit
   !

   coefs = matmul(vars%acoefs,q0s)
   dd1   = maxval(abs(coefs))
   dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
   if (dd1 .eq. 0) dd1 = 1
   dd    = dd2/dd1
   if (dd .gt. eps) then
      ier = 1024
      return
   endif

   coefs = matmul(vars%acoefs,q1s)
   dd1   = maxval(abs(coefs))
   dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
   if (dd1 .eq. 0) dd1 = 1
   dd    = dd2/dd1
   if (dd .gt. eps) then
      ier = 1024
      return
   endif



   !
   !  Use the eigensolver
   !

   call init_random_seed()
   allocate(ws(kcheb,2))

   n          = 2
   its        = 200

   ! do i=1,kcheb

   !    rcoefs(1) = real(q1s(i))
   !    icoefs(1) = imag(q1s(i))

   !    rcoefs(2) = real(q0s(i))
   !    icoefs(2) = imag(q0s(i))

   !    call zamvw(n,rcoefs,icoefs,reigs,ieigs,its,flag)

   !    if (flag .ne. 0) then
   !       ier = 256
   !       return
   !    endif

   !    w(1) = reigs(1) + ima*ieigs(1)
   !    w(2) = reigs(2) + ima*ieigs(2)


   !    if (i .eq. 1) then
   !          call scalar_sortw(2,w)
   !          idx1 = 1
   !          idx2 = 2
   !       else
   !          idx1 = minloc(abs(w(1)-ws(i-1,:)),1)
   !          idx2 = minloc(abs(w(2)-ws(i-1,:)),1)
   !    endif


   !    ws(i,idx1) = w(1)
   !    ws(i,idx2) = w(2)

   ! end do

   dk = real(pars(1))
   call root_estimates(2, dk, roots)



   do isol=1,n

      dd   = 2/(b-a)
      r    = roots(isol)
      rp   = dd * matmul(vars%adiff,r)


      do iter=1,maxiters
         dd      = 2/(b-a)
         ifdone  = 0
         rhs     = -(rp+r**2+q1s*r+q0s)


         amatr   = dd * vars%adiff


         do i=1,kcheb
            amatr(i,i) = amatr(i,i) + (2*r(i) + q1s(i))
         end do

      
         
         !call linalg0_utvsolve(epssol,kcheb,kcheb,amatr,delta,rhs)
         call linalg0_qrtsolve(epssol,kcheb,kcheb,amatr,delta,rhs)

      

         ! perform at least one iteration
         if (iter .gt. 1) then
         dd1 = norm2(abs(delta))
         dd0 = norm2(abs(r))+1.0d0
         dd2  = dd1/dd0
         !print *, dd2
         if (dd2 .lt. eps) ifdone = 1
         endif

         r    = r + delta

         rp   = dd * matmul(vars%adiff,r)

         if (ifdone .eq. 1) exit

      end do


      
      if (iter .gt. maxiters) then
         ier = 4
         return
      endif

      

      if (ifleft .eq. 1) then
            rvals(1,isol) = r(1)
         else
            rvals(1,isol) = r(kcheb)
      endif

   end do

    
end subroutine


subroutine scalar_riccati2(vars,ier,eps,a,b,c,rvals,fun,pars,chebsol,rs)
implicit double precision (a-h,o-z)
type(scalar_vars_t)                         :: vars
procedure(scalar_fun)                       :: fun
double complex, pointer                     :: pars(:)
double complex                              :: rvals(:)
type(chebpw_scheme)                         :: chebsol
double complex, allocatable, intent(out)    :: rs(:,:)
!
!  Solve the Riccati equation corresponding to (1) over the interval [a,b] 
!  given the values of r_j(t) and its first 2 derivatives at a point t in [a,b].
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    eps - precision for the Newton iteration
!    (a,b) - the interval on which to perform the procedure
!    c - the point at which the values of r_j(t) and its derivatives are given
!    rvals - an array of length 3 which gives the values of r_j(c), r_j'(c), and r_j''(c)
!    fun - an external subroutine supplying the coefficients of the 
!      user's scalar system 
!    par? - user-supplied parameters which are passed to fun
!
!  Output parameters:
!    ier  - a error return code;
!      ier = 0  indicates successful execution
!      ier = 4  the ODE solver failed
!
!    chebsol - a structure describing the discretization scheme used to represent
!       the solution
!    rs - the jth column of this array gives the values of r^(j-1) at the discretization 
!         nodes
!
!
type(c_ptr)                                :: userptr
type(scalar_odedata), pointer              :: odedata
double complex, allocatable                :: ys(:,:), yders(:,:), rr(:)
double complex                             :: ima


ier   = 0
eps0  = epsilon(0.0d0)
ima   = (0.0d0,1.0d0)

allocate(odedata)
userptr = c_loc(odedata)

n              = 2
odedata%n      = n-1
odedata%fun   => fun
odedata%pars  => pars
odedata%wa     = 0
odedata%wb     = 0
odedata%ifleft = -1
ifit           = 2

call odesolve_nonlinear1(vars%odesolve_vars,jer,ifit,eps,a,b,c,rvals(2:2),n-1,&
  scalar_ode2,userptr,chebsol,ys,yders)
if (jer .ne. 0) then
ier = 4
return
endif




call chebpw_info(chebsol,k,nints)
allocate(rs(k*nints,3))

rs(:,3) = yders(:,1)
rs(:,2) = ys(:,1)
call chebpw_int(chebsol,rs(:,2),c,rvals(1),rr)
rs(:,1) = rr

end subroutine


subroutine scalar_levinode2(vars,ier,eps,a,b,c,rvals,fun,pars,chebsol,rs)
implicit double precision (a-h,p-z)
type(scalar_vars_t)                         :: vars
type(chebpw_scheme), intent(out)            :: chebsol
double complex, allocatable, intent(out)    :: rs(:,:,:)
double complex                              :: rvals(:)
procedure(scalar_fun)                       :: fun
double complex, pointer                     :: pars(:)
!
!  Use the global Levin method to compute a set of phase functions.
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    eps - precision for the Newton iteration
!    (a,b) - the interval on which to perform the procedure
!    c - the point at which the values of r_j(t) and its derivatives are given
!    rvals - the desired values of the r_j
!    fun - an external subroutine supplying the coefficients of the 
!      user's scalar system 
!    par? - user-supplied parameters which are passed to fun
!  
!
!  Output parameters:
!    ier  - a error return code;
!      ier = 0   indicates successful execution
!      ier = 4   the Newton iterations did not converge
!      ier = 256 the eigensolver failed
!
!    chebsol - a structure describing the discretization scheme used to represent
!      the solution
!    rs - a (*,n,n) such that r(:,i,j) gives the values of r_j^(i-1) at the 
!      discretization nodes
!
!
double precision, allocatable               :: stack(:,:), ts(:)
double complex, allocatable                 :: qs(:,:), r(:), rp(:), rpp(:), coefs(:), ws(:,:)
double complex, allocatable                 :: delta(:), rhs(:), amatr(:,:)

double complex, allocatable                 :: rs0(:,:,:,:), rvals0(:)

double precision, allocatable               :: ab(:,:)

double complex                              :: ima, qs0(2), w(2)
double precision                            :: rcoefs(2), icoefs(2), ieigs(2), reigs(2), drs(2)
integer                                     :: its(2), flag

n        = 2
ima      = (0.0d0,1.0d0)
eps0     = epsilon(0.0d0)
epssol   = eps0*10

ier      = 0
ind      = 0
nints    = 0
kcheb    = vars%kcheb
ntest    = vars%ntest
maxiters = vars%maxleviters
maxints  = vars%maxints
maxstack = vars%maxstack

allocate(ts(kcheb), qs(kcheb,2), r(kcheb), rp(kcheb), rpp(kcheb), coefs(kcheb) )
allocate(delta(kcheb), rhs(kcheb), amatr(kcheb,kcheb))
allocate(rs0(kcheb,2,2,maxints), stack(2,maxstack), ab(2,maxints))

call init_random_seed()
allocate(ws(kcheb,2))

nstack     = 1
stack(1,1) = a
stack(2,1) = b

do while(nstack > 0)

   !
   !  Pop an interval off of the stack
   !
    
   a0     = stack(1,nstack)
   b0     = stack(2,nstack)
   nstack = nstack-1
   ind    = ind+1

   ifsplit = 1
   dcoefs  = 1d300
   drs     = 1d300
   jer     = 0

   !
   !  Fetch the values of the coefficients
   !
    
   ts = vars%xscheb*(b0-a0)/2 + (b0+a0)/2
   do i=1,kcheb
   call fun(n,ts(i),qs0,pars)
   qs(i,1) = qs0(1)
   qs(i,2) = qs0(2)
   end do    

   !
   !  Check that the coefficients are properly discretized
   !

   dcoefs = 0
   do i=1,2
      coefs  = matmul(vars%acoefs,qs(:,i))
      dd1    = maxval(abs(coefs))+1
      dd2    = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
      dd     = dd2/dd1
      dcoefs = max(dd,dcoefs)
   end do

   !
   !  Compute the eigenvalues to use as inital guesses if this is the first interval
   !

   if (dcoefs .lt. eps .AND. nints .eq. 0) then

      n          = 2
      its        = 200
       
      do i=1,kcheb
       
         rcoefs(1) = real(qs(i,2))
         icoefs(1) = imag(qs(i,2))
         
         rcoefs(2) = real(qs(i,1))
         icoefs(2) = imag(qs(i,1))
         
         call zamvw(n,rcoefs,icoefs,reigs,ieigs,its,flag)
         if (flag .ne. 0) then
            ier = 256
            return
         endif
         
         w(1) = reigs(1) + ima*ieigs(1)
         w(2) = reigs(2) + ima*ieigs(2)
         
         
         if (i .eq. 1) then
         call scalar_sortw(2,w)
         idx1 = 1
         idx2 = 2
         else
         
         idx1 = minloc(abs(w(1)-ws(i-1,:)),1)
         idx2 = minloc(abs(w(2)-ws(i-1,:)),1)
         
         endif
         
         ws(i,idx1) = w(1)
         ws(i,idx2) = w(2)
       
      end do


   endif
    
    
   !
   !  If the coefficients are properly discretized, proceed
   !
    
   if (jer .eq. 0 .AND. dcoefs .lt. eps) then
         
      iffail  = 0

      do isol=1,2
         
         if (nints .eq. 0 ) then
               dd     = 2/(b0-a0)
               r      = ws(:,isol)
               rp     = dd*matmul(vars%adiff,r)
            else
               r      = rs0(kcheb,1,isol,nints)
               rp     = rs0(kcheb,2,isol,nints)
         endif

         ifdone = 0

         do iter=1,maxiters
            dd      = 2/(b0-a0)
            rhs     = -(rp+r**2+qs(:,2)*r+qs(:,1))

            amatr   = dd*vars%adiff

            
            do i=1,kcheb
            amatr(i,i) = amatr(i,i) + 2*r(i) + qs(i,2)
         end do

         call linalg0_qrtsolve(epssol,kcheb,kcheb,amatr,delta,rhs)

         
         !perform at least one iteration
         if (iter .gt. 1) then
         dd1 = norm2(abs(delta))
         dd0 = norm2(abs(r))+1.0d0
         dd2  = dd1/dd0
         if (dd2 .lt. eps) ifdone=1
         endif
         
         r    = r + delta
         rp   = dd * matmul(vars%adiff,r)


         ! write(*,"(3X,3(I3,1X),D15.7,2(I3,1X),1(D15.7,1X),(I3,1X))") ind,nints,jer,dcoefs,isol,iter,dd2,ifdone
         ! write(13,"(3X,3(I3,1X),D15.7,2(I3,1X),1(D15.7,1X),(I3,1X))") ind,nints,jer,dcoefs,isol,iter,dd2,ifdone

         if (ifdone .eq. 1) exit
         end do
         
       
         if (iter .gt. maxiters) then
         iffail = 1
         exit
         endif

         rs0(:,1,isol,nints+1) = r
         rs0(:,2,isol,nints+1) = rp

         coefs     = matmul(vars%acoefs,r)
         dd1       = maxval(abs(coefs))+1
         dd2       = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
         dd        = dd2/dd1
         drs(isol) = max(dd,dcoefs)
              
      end do

      if (iffail .eq. 0 .AND. maxval(drs) .lt. eps) ifsplit = 0
    
   endif

   ! write(*,"(3(I5,1X),5(D15.7,1X),I2)")  ind,nints,jer,a0,b0,dcoefs,drs,ifsplit
   ! write(13,"(3(I5,1X),5(D15.7,1X),I2)") ind,nints,jer,a0,b0,dcoefs,drs,ifsplit

   !
   !  Either split the interval or add it to the list of "accepted" intervals
   !
   if (ifsplit .eq. 1) then
         if (nstack +2 .ge. maxstack) then
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
         if (nints+1 .ge. maxints) then
            ier = 4
            return
         endif

         nints       = nints+1
         ab(1,nints) = a0
         ab(2,nints) = b0

   endif

end do

!
!  Copy out the solution
!

call chebpw_specified(chebsol,kcheb,nints,ab)

nn = nints*kcheb
allocate(rs(nn,3,2))

i1 = 1
i2 = kcheb

do int=1,nints
   a0  = ab(1,int)
   b0  = ab(2,int)

   rs(i1:i2,2,1) = rs0(1:kcheb,1,1,int)
   rs(i1:i2,3,1) = rs0(1:kcheb,2,1,int)

   rs(i1:i2,2,2) = rs0(1:kcheb,1,2,int)
   rs(i1:i2,3,2) = rs0(1:kcheb,2,2,int)

   i1 = i1+kcheb
   i2 = i2+kcheb

end do

call chebpw_int2(chebsol,rs(:,2,1),c,rvals(1),rvals0)
rs(:,1,1) = rvals0
call chebpw_int2(chebsol,rs(:,2,2),c,rvals(2),rvals0)
rs(:,1,2) = rvals0

end subroutine



subroutine scalar_residual2(vars,fun,pars,chebsol,rs,dres)
implicit double precision (a-h,o-z)
type(scalar_vars_t)                         :: vars
procedure(scalar_fun)                       :: fun
double complex, pointer                     :: pars(:)
type(chebpw_scheme)                         :: chebsol
double complex                              :: rs(:,:)
!
!  Return a measure of the error obtained when a solution is inserted into 
!  the Riccati equation.
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    (a,b) - the interval on which to perform the procedure
!    fun - an external subroutine supplying the coefficients of the 
!      user's scalar system 
!    par? - user-supplied parameters which are passed to fun
!    chebsol - a structure describing the discretization scheme used to represent
!       the solution
!    rs - the values of the solution at the discretization nodes
!
!  Output parameters:
!    dres - a measure of the residual
!

double precision, allocatable :: ts(:),errs(:)
double complex                :: r, rp, qs(2), q0, q1


eps0  = epsilon(0.0d0)
ima   = (0.0d0,1.0d0)

n     = 2
nn    = 10000
a     = chebsol%a
b     = chebsol%b
allocate(ts(nn),errs(nn))


do i=1,nn
ts(i) = a+(b-a)*(i-1.0d0)/(nn-1.0d0)
end do

errs  = 1d300

do i=1,nn
t  = ts(i)
call fun(n,t,qs,pars)
q0 = qs(1)
q1 = qs(2)

call chebpw_interp(chebsol,rs(:,2),t,r)
call chebpw_interp(chebsol,rs(:,3),t,rp)

dd      = max(abs(q0),abs(q1))
errs(i) = abs(rp + r**2 + q1*r + q0) / dd

end do


!call prin2("errs = ",errs) 
dres = maxval(errs)

end subroutine


subroutine scalar_ode2(n,t,y,f,df,userptr)
implicit double precision (a-h,p-z)
type(c_ptr)                    :: userptr
double complex                 :: y(:), f(:), df(:,:)
type(scalar_odedata), pointer  :: odedata
double complex                 :: us(2)
double complex                 :: q0, q1, q2, r, rp, rpp
double complex                 :: dk0, dk1, dk2

call c_f_pointer(userptr,odedata)

ifleft = odedata%ifleft
wa     = odedata%wa
wb     = odedata%wb
phi    = (erf(wa*(t-wb))+1)/2




call odedata%fun(n+1,t,us,odedata%pars)



r      = y(1)



if (ifleft .eq. 1) then
dk0 = odedata%dks(1)
dk1 = odedata%dks(2)

q0  = phi*dk0 + (1-phi)*us(1)
q1  = phi*dk1 + (1-phi)*us(2)
elseif (ifleft .eq. 0) then
dk0 = odedata%dks(1)
dk1 = odedata%dks(2)

q0  = (1-phi)*dk0 + phi*us(1)
q1  = (1-phi)*dk1 + phi*us(2)
else
q0  = us(1)
q1  = us(2)
endif

f(1)    = -(r**2+q1*r+q0)


df(1,1) = -(2*r+q1)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Order 3 using the standard ODE solver
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine scalar_homotope3(vars,ier,eps,a,b,ifleft,fun,pars,rvals)
implicit double precision (a-h,o-z)
type(scalar_vars_t)              :: vars
procedure(scalar_fun)            :: fun
double complex, pointer          :: pars(:)
double complex                   :: rvals(:,:)
!
!  Use the homotopy method to compute the values of the derivatives of the r_j(t)
!  at a specified point.  More explicitly, this routine homotopes the coefficient 
!  matrix for the user's nth order scalar system to the constant system
!
!          (  0     1    0  )
!  A0    = (  0     0    1  )
!          ( -dk0 -dk1 -dk2 )
!
!  in order to calculate appropriate initial conditions for the user's system.
!  Reasonable values of dk? are determined by calling the user-supplied external
!  subroutine.
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    eps - precision for the calculations
!    (a,b) - the interval on which to perform the procedure
!    ifleft - an integer parameter indicating whether the values of the r_j(t) should
!      be computed at the left-hand endpoint of the interval a or the right-hand
!      endpoint b;
!
!        ifleft = 1   means calculate the values at a
!        ifleft = 0   means calculate the values at b
!
!    fun - an external subroutine supplying the coefficients of the user's scalar
!      equation
!    par? - user-supplied parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0    indicates successful execution
!       ier = 4    indicates the ODE solver failed
!       ier = 256  indicates that the eigensolver failed
!
!    rvals - the (i,j) entry contains the value of r_j^(i)(c), where c=a if ifleft =1
!      and c=b if ifleft=0
!    

type(chebpw_scheme)               :: chebsol
type(c_ptr)                       :: userptr
type(scalar_odedata), pointer     :: odedata
double complex, allocatable       :: ys(:,:), yders(:,:)
double complex                    :: ima, yc(2), y0, y0p, w(3), dks(3)
double precision                  :: rcoefs(3), icoefs(3), reigs(3), ieigs(3)
integer                           :: its(3), flag

ier  = 0

k    = vars%kcheb
eps0 = epsilon(0.0d0)


wa   = 11d0/(b-a)
wb   = (a+b)/2.0d0
ima  = (0.0d0,1.0d0)

n    = 3

call fun(n,(a+b)/2,dks,pars)

allocate(odedata)
allocate(odedata%dks(n))
odedata%n       = n-1
odedata%pars   => pars
odedata%dks     = dks
odedata%fun    => fun
odedata%wa      = wa
odedata%wb      = wb
odedata%ifleft  = ifleft
userptr         = c_loc(odedata)


if (ifleft .eq. 0) then
c            = a
else
c            = b
endif

!
!  Find the roots of the polynomial
!

n          = 3
its        = 200

rcoefs(1) = real(dks(3))
icoefs(1) = imag(dks(3))

rcoefs(2) = real(dks(2))
icoefs(2) = imag(dks(2))

rcoefs(3) = real(dks(1))
icoefs(3) = imag(dks(1))

! call prin2("rcoefs = ",rcoefs)
! call prin2("icoefs = ",icoefs)

call init_random_seed()
call zamvw(n,rcoefs,icoefs,reigs,ieigs,its,flag)

! call prini("flag = ",flag)

if (flag .ne. 0) then
ier = 256
return
endif

w(1) = reigs(1) + ima*ieigs(1)
w(2) = reigs(2) + ima*ieigs(2)
w(3) = reigs(3) + ima*ieigs(3)

!call scalar_sortw(3,w)

! call prinz("in scalar_homotopy3, ws = ",w)

if (norm2(abs(w)) .eq. 0) then
ier = 256
return
endif

do i=1,n
yc(1)      = w(i)
yc(2)      = 0
ifit       = 2

call odesolve_nonlinear1(vars%odesolve_vars,jer,ifit,eps,a,b,c,yc,n-1, &
  scalar_ode3,userptr,chebsol,ys,yders)

if (jer .ne. 0) then
ier = 4
return
endif


nn = chebsol%nints*k
if (ifleft .eq. 1) then
rvals(1,i)  = ys(1,1)
rvals(2,i)  = ys(1,2)
else
rvals(1,i)  = ys(nn,1)
rvals(2,i)  = ys(nn,2)
endif


end do

end subroutine


subroutine scalar_levin3(vars,ier,eps,a,b,ifleft,fun,pars,rvals)
implicit double precision (a-h,o-z)
type(scalar_vars_t)      :: vars
procedure(scalar_fun)    :: fun
double complex, pointer  :: pars(:)
double complex           :: rvals(:,:)
!
!  Use the Levin method to compute the values of the derivatives of the r_j(t)
!  at a specified point.
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    eps - precision for the Newton iteration
!    (a,b) - the interval on which to perform the procedure
!    ifleft - an integer parameter indicating whether the values of the r_j(t) should
!      be computed at the left-hand endpoint of the interval a or the right-hand
!      endpoint b;
!
!        ifleft = 1   means calculate the values at a
!        ifleft = 0   means calculate the values at b
!
!    fun - an external subroutine supplying the coefficients of the 
!      user's scalar system 
!    par? - user-supplied parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0     indicates successful execution
!       ier = 4     indicates the ODE solver failed
!       ier = 256   means that the eigensolver failed
!       ier = 1024  indicates the coefficients were not represented resolved
!                   by Chebyshev expansions of order vars%kcheb over [a,b]
!
!    rvals - the (i,j) entry contains the value of r_j^(i)(c), where c=a if ifleft =1
!      and c=b if ifleft=0
!

double precision, allocatable            :: ts(:)
double precision                         :: rcoefs(3), icoefs(3), reigs(3), ieigs(3)
integer                                  :: its(3), flag

double complex, allocatable              :: q0s(:), q1s(:), q2s(:), qs(:), coefs(:), ws(:,:)

double complex, allocatable              :: r(:), rp(:), rpp(:), rhs(:), amatr(:,:), delta(:)
double complex                           :: ima, rval, rder, w(3)
double complex                           :: rp0, rpp0

double complex, allocatable             :: roots(:)

kcheb    = vars%kcheb
ntest    = vars%ntest
maxiters = vars%maxleviters


ier      = 0
ima      = (0.0d0,1.0d0)
eps0     = epsilon(0.0d0)
epssol   = eps0*10

n        = 3
allocate( ts(kcheb), q0s(kcheb), q1s(kcheb), q2s(kcheb),  qs(n), coefs(kcheb)  )
allocate( r(kcheb), rp(kcheb), rpp(kcheb), rhs(kcheb), amatr(kcheb,kcheb), delta(kcheb) )

!allocate( rs(kcheb), rders(kcheb), delta(kcheb), rhs(kcheb), amatr(kcheb,kcheb) )

ts = (b-a)/2*vars%xscheb + (b+a)/2

do i=1,kcheb
call fun(n,ts(i),qs,pars)
q0s(i) = qs(1)
q1s(i) = qs(2)
q2s(i) = qs(3)
end do

! 
!  Test the coefficinets for fit
!

coefs = matmul(vars%acoefs,q0s)
dd1   = maxval(abs(coefs))
dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
if (dd1 .eq. 0) dd1 = 1
dd    = dd2/dd1
if (dd .gt. eps) then
ier = 1024
return
endif



coefs = matmul(vars%acoefs,q1s)
dd1   = maxval(abs(coefs))
dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
if (dd1 .eq. 0) dd1 = 1
dd    = dd2/dd1
if (dd .gt. eps) then
ier = 1024
return
endif

coefs = matmul(vars%acoefs,q2s)
dd1   = maxval(abs(coefs))
dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
if (dd1 .eq. 0) dd1 = 1
dd    = dd2/dd1
if (dd .gt. eps) then
ier = 1024
return
endif



!
!  Use the eigensolver
!

call init_random_seed()
allocate(ws(kcheb,3))

n          = 3
its        = 200

! do i=1,kcheb

! rcoefs(1) = real(q2s(i))
! icoefs(1) = imag(q2s(i))

! rcoefs(2) = real(q1s(i))
! icoefs(2) = imag(q1s(i))

! rcoefs(3) = real(q0s(i))
! icoefs(3) = imag(q0s(i))


! call zamvw(n,rcoefs,icoefs,reigs,ieigs,its,flag)

! if (flag .ne. 0) then
! ier = 256
! return
! endif

! w(1) = reigs(1) + ima*ieigs(1)
! w(2) = reigs(2) + ima*ieigs(2)
! w(3) = reigs(3) + ima*ieigs(3)



! if (i .eq. 1) then
! call scalar_sortw(3,w)

! idx1 = 1
! idx2 = 2
! idx3 = 3
! else

! idx1 = minloc(abs(w(1)-ws(i-1,:)),1)
! idx2 = minloc(abs(w(2)-ws(i-1,:)),1)
! idx3 = minloc(abs(w(3)-ws(i-1,:)),1)

! endif



! ws(i,idx1) = w(1)
! ws(i,idx2) = w(2)
! ws(i,idx3) = w(3)


! end do

dk = real(pars(1))
call root_estimates(3, dk, roots)


do isol=1,n

dd   = 2/(b-a)
r    = roots(isol)
rp   = dd * matmul(vars%adiff,r)
rpp  = dd**2 * matmul(vars%adiff2,r)


do iter=1,maxiters
dd      = 2/(b-a)
ifdone  = 0
rhs     = -(rpp+3*r*rp+r**3+q2s*rp+q2s*r**2+q1s*r+q0s)
amatr   = dd**2 * vars%adiff2 


do i=1,kcheb
   amatr(i,i) = amatr(i,i) + 3*rp(i) + 3*r(i)**2 + 2*q2s(i)*r(i)+q1s(i)
end do


do i=1,kcheb
amatr(i,:) = amatr(i,:) + (3*r(i)+q2s(i))*dd*vars%adiff(i,:)
end do


!call linalg0_utvsolve(epssol,kcheb,kcheb,amatr,delta,rhs)
call linalg0_qrtsolve(epssol,kcheb,kcheb,amatr,delta,rhs)


! perform at least one iteration
if (iter .gt. 1) then
dd1 = norm2(abs(delta))
dd0 = norm2(abs(r))+1.0d0
dd2  = dd1/dd0
if (dd2 .lt. eps) ifdone = 1
endif

r    = r + delta
rp   = dd * matmul(vars%adiff,r)
rpp  = dd * matmul(vars%adiff,rp)


!print *,isol,iter,dd2,eps,ifdone

if (ifdone .eq. 1) exit

end do


if (iter .gt. maxiters) then
ier = 4
return
endif


if (ifleft .eq. 1) then
rvals(1,isol) = r(1)
rvals(2,isol) = rp(1)
else
rvals(1,isol) = r(kcheb)
rvals(2,isol) = rp(kcheb)
endif



end do
    
end subroutine


subroutine scalar_riccati3(vars,ier,eps,a,b,c,rvals,fun,pars,chebsol,rs)
implicit double precision (a-h,o-z)
type(scalar_vars_t)                         :: vars
procedure(scalar_fun)                       :: fun
double complex, pointer                     :: pars(:)
double complex                              :: rvals(:)
type(chebpw_scheme)                         :: chebsol
double complex, allocatable, intent(out)    :: rs(:,:)
!
!  Solve the Riccati equation corresponding to (1) over the interval [a,b] 
!  given the values of r_j(t) and its first 2 derivatives at a point t in [a,b].
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    eps - precision for the Newton iteration
!    (a,b) - the interval on which to perform the procedure
!    c - the point at which the values of r_j(t) and its derivatives are given
!    rvals - an array of length 3 which gives the values of r_j(c), r_j'(c), and r_j''(c)
!    fun - an external subroutine supplying the coefficients of the 
!      user's scalar system 
!    par? - user-supplied parameters which are passed to fun
!
!  Output parameters:
!    ier  - a error return code;
!      ier = 0  indicates successful execution
!      ier = 4  the ODE solver failed
!
!    chebsol - a structure describing the discretization scheme used to represent
!       the solution
!    rs - the values of the solution at the discretization nodes
!
!
type(c_ptr)                                :: userptr
type(scalar_odedata), pointer              :: odedata
double complex, allocatable                :: ys(:,:), yders(:,:), rr(:)
double complex                             :: ima


ier   = 0
eps0  = epsilon(0.0d0)
ima   = (0.0d0,1.0d0)

allocate(odedata)
userptr = c_loc(odedata)

n              = 3
odedata%n      = n-1
odedata%fun   => fun
odedata%pars  => pars
odedata%wa     = 0
odedata%wb     = 0
odedata%ifleft = -1
ifit           = 2



call odesolve_nonlinear1(vars%odesolve_vars,jer,ifit,eps,a,b,c,rvals(2:3),n-1,&
  scalar_ode3,userptr,chebsol,ys,yders)


! call prinz("", yders)
! stop

if (jer .ne. 0) then
ier = 4
return
endif



call chebpw_info(chebsol,k,nints)


allocate(rs(k*nints,4))

rs(:,4) = yders(:,2)
rs(:,3) = ys(:,2)
rs(:,2) = ys(:,1)

call chebpw_int(chebsol,rs(:,2),c,rvals(1),rr)
rs(:,1) = rr

end subroutine


subroutine scalar_levinode3(vars,ier,eps,a,b,c,rvals,fun,pars,chebsol,rs)
implicit double precision (a-h,p-z)
type(scalar_vars_t)                         :: vars
type(chebpw_scheme), intent(out)            :: chebsol
double complex, allocatable, intent(out)    :: rs(:,:,:)
double complex                              :: rvals(:)
procedure(scalar_fun)                       :: fun
double complex, pointer                     :: pars(:)

!
!  Use the global Levin method to compute a set of phase functions.
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    eps - precision for the Newton iteration
!    (a,b) - the interval on which to perform the procedure
!    c - the point at which the values of r_j(t) and its derivatives are given
!    rvals - the desired values of the r_j
!    fun - an external subroutine supplying the coefficients of the 
!      user's scalar system 
!    par? - user-supplied parameters which are passed to fun
!  
!
!  Output parameters:
!    ier  - a error return code;
!      ier = 0  indicates successful execution
!      ier = 4  the Newton iterations did not converge
!
!    chebsol - a structure describing the discretization scheme used to represent
!      the solution
!    rs - a (*,n,n) such that r(:,i,j) gives the values of r_j^(i-1) at the 
!      discretization nodes
!
!

double precision, allocatable               :: stack(:,:), ts(:)
double complex, allocatable                 :: qs(:,:), r(:), rp(:), rpp(:), coefs(:), ws(:,:)
double complex, allocatable                 :: delta(:), rhs(:), amatr(:,:), rvals0(:)

double complex, allocatable                 :: rs0(:,:,:,:)

double precision, allocatable               :: ab(:,:)

double complex                              :: ima, qs0(3), w(3)
double precision                            :: rcoefs(3), icoefs(3), ieigs(3), reigs(3), drs(3)
integer                                     :: its(3), flag

n        = 3
ima      = (0.0d0,1.0d0)
eps0     = epsilon(0.0d0)
epssol   = eps0*10

ier      = 0
ind      = 0
nints    = 0
kcheb    = vars%kcheb
ntest    = vars%ntest
maxiters = vars%maxleviters
maxints  = vars%maxints
maxstack = vars%maxstack

allocate(ts(kcheb), qs(kcheb,3), r(kcheb), rp(kcheb), rpp(kcheb), coefs(kcheb) )
allocate(delta(kcheb), rhs(kcheb), amatr(kcheb,kcheb))
allocate(rs0(kcheb,3,3,maxints), stack(2,maxstack), ab(2,maxints))

call init_random_seed()
allocate(ws(kcheb,3))

nstack     = 1
stack(1,1) = a
stack(2,1) = b

do while(nstack > 0)

   !
   !  Pop an interval off of the stack
   !
    
   a0     = stack(1,nstack)
   b0     = stack(2,nstack)
   nstack = nstack-1
   ind    = ind+1

   ifsplit = 1
   dcoefs  = 1d300
   drs     = 1d300
   jer     = 0

   !
   !  Fetch the values of the coefficients
   !
    
   ts = vars%xscheb*(b0-a0)/2 + (b0+a0)/2
   do i=1,kcheb
      call fun(n,ts(i),qs0,pars)
      qs(i,1) = qs0(1)
      qs(i,2) = qs0(2)
      qs(i,3) = qs0(3)
   end do    

   !
   !  Check that the coefficients are properly discretized
   !

   dcoefs = 0
   do i=1,3
   coefs  = matmul(vars%acoefs,qs(:,i))
   dd1    = maxval(abs(coefs))+1
   dd2    = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
   dd     = dd2/dd1

   dcoefs = max(dd,dcoefs)
   end do
   !
   !  Compute the eigenvalues to use as inital guesses if this is the first interval
   !

   if (dcoefs .lt. eps .AND. nints .eq. 0) then
      its        = 200
       
      do i=1,kcheb
       
      rcoefs(1) = real(qs(i,3))
      icoefs(1) = imag(qs(i,3))
       
      rcoefs(2) = real(qs(i,2))
      icoefs(2) = imag(qs(i,2))
       
      rcoefs(3) = real(qs(i,1))
      icoefs(3) = imag(qs(i,1))
       
      call zamvw(n,rcoefs,icoefs,reigs,ieigs,its,flag)
      if (flag .ne. 0) then
      ier = 256
      return
      endif
       
      w(1) = reigs(1) + ima*ieigs(1)
      w(2) = reigs(2) + ima*ieigs(2)
      w(3) = reigs(3) + ima*ieigs(3)
       
      if (i .eq. 1) then
      call scalar_sortw(3,w)
       
      idx1 = 1
      idx2 = 2
      idx3 = 3
      else
       
      idx1 = minloc(abs(w(1)-ws(i-1,:)),1)
      idx2 = minloc(abs(w(2)-ws(i-1,:)),1)
      idx3 = minloc(abs(w(3)-ws(i-1,:)),1)
       
      endif
       
      ws(i,idx1) = w(1)
      ws(i,idx2) = w(2)
      ws(i,idx3) = w(3)
       
      end do
       

   endif
   
   !
   !  If the coefficients are properly discretized, proceed
   !
    
   if (jer .eq. 0 .AND. dcoefs .lt. eps) then

         
      iffail  = 0

      do isol=1,3
         
         if (nints .eq. 0 ) then
         dd   = 2/(b0-a0)
         r    = ws(:,isol)
         rp   = dd    * matmul(vars%adiff,r)
         rpp  = dd**2 * matmul(vars%adiff2,r)
         else
         r      = rs0(kcheb,1,isol,nints)
         rp     = rs0(kcheb,2,isol,nints)
         rpp    = rs0(kcheb,3,isol,nints)
         endif

         ifdone = 0

         do iter=1,maxiters
         dd      = 2/(b0-a0)
         rhs     = -(rpp+3*r*rp+r**3+qs(:,3)*rp+qs(:,3)*r**2+qs(:,2)*r+qs(:,1))

         amatr   = dd**2 * vars%adiff2 

         do i=1,kcheb
         amatr(i,:) = amatr(i,:) + (3*r(i)+qs(i,3))*dd*vars%adiff(i,:)
         end do
         
         do i=1,kcheb
         amatr(i,i) = amatr(i,i) + 3*rp(i) + 3*r(i)**2 + 2*qs(i,3)*r(i)+qs(i,2)
         end do

         call linalg0_qrtsolve(epssol,kcheb,kcheb,amatr,delta,rhs)

         
         !perform at least one iteration
         if (iter .gt. 1) then
         dd1 = norm2(abs(delta))
         dd0 = norm2(abs(r))+1.0d0
         dd2  = dd1/dd0
         if (dd2 .lt. eps) ifdone=1
         endif
         
         r    = r + delta
         rp   = dd * matmul(vars%adiff,r)
         rpp  = dd * matmul(vars%adiff,rp)

         ! write(*,"(3X,3(I3,1X),D15.7,2(I3,1X),1(D15.7,1X),(I3,1X))") ind,nints,jer,dcoefs,isol,iter,dd2,ifdone
         ! write(13,"(3X,3(I3,1X),D15.7,2(I3,1X),1(D15.7,1X),(I3,1X))") ind,nints,jer,dcoefs,isol,iter,dd2,ifdone

         if (ifdone .eq. 1) exit

         end do
         
         if (iter .gt. maxiters) then
         iffail = 1
         exit
         endif

         rs0(:,1,isol,nints+1) = r
         rs0(:,2,isol,nints+1) = rp
         rs0(:,3,isol,nints+1) = rpp

         coefs     = matmul(vars%acoefs,r)
         dd1       = maxval(abs(coefs))+1
         dd2       = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
         dd        = dd2/dd1
         drs(isol) = max(dd,dcoefs)
              
      end do

      if (iffail .eq. 0 .AND. maxval(drs) .lt. eps) ifsplit = 0
    
   endif

   ! write(*,"(3(I5,1X),6(D15.7,1X),I2)")  ind,nints,jer,a0,b0,dcoefs,drs,ifsplit
   ! write(13,"(3(I5,1X),6(D15.7,1X),I2)") ind,nints,jer,a0,b0,dcoefs,drs,ifsplit

   !
   !  Either split the interval or add it to the list of "accepted" intervals
   !
   if (ifsplit .eq. 1) then
   if (nstack +2 .ge. maxstack) then
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
   if (nints+1 .ge. maxints) then
   ier = 4
   return
   endif

   nints       = nints+1
   ab(1,nints) = a0
   ab(2,nints) = b0

   endif

end do

!
!  Copy out the solution
!

call chebpw_specified(chebsol,kcheb,nints,ab)

allocate(rs(nints*kcheb,4,3))


i1 = 1
i2 = kcheb

do int=1,nints
a0  = ab(1,int)
b0  = ab(2,int)

rs(i1:i2,2,1) = rs0(1:kcheb,1,1,int)
rs(i1:i2,3,1) = rs0(1:kcheb,2,1,int)
rs(i1:i2,4,1) = rs0(1:kcheb,3,1,int)


rs(i1:i2,2,2) = rs0(1:kcheb,1,2,int)
rs(i1:i2,3,2) = rs0(1:kcheb,2,2,int)
rs(i1:i2,4,2) = rs0(1:kcheb,3,2,int)


rs(i1:i2,2,3) = rs0(1:kcheb,1,3,int)
rs(i1:i2,3,3) = rs0(1:kcheb,2,3,int)
rs(i1:i2,4,3) = rs0(1:kcheb,3,3,int)

i1 = i1+kcheb
i2 = i2+kcheb

end do


call chebpw_int2(chebsol,rs(:,2,1),c,rvals(1),rvals0)
rs(:,1,1) = rvals0

call chebpw_int2(chebsol,rs(:,2,2),c,rvals(2),rvals0)
rs(:,1,2) = rvals0

call chebpw_int2(chebsol,rs(:,2,3),c,rvals(3),rvals0)
rs(:,1,3) = rvals0


end subroutine


subroutine scalar_residual3(vars,fun,pars,chebsol,rs,dres)
implicit double precision (a-h,o-z)
type(scalar_vars_t)                         :: vars
procedure(scalar_fun)                       :: fun
double complex, pointer                     :: pars(:)
type(chebpw_scheme)                         :: chebsol
double complex                              :: rs(:,:)
!
!  Return a measure of the error obtained when a solution is inserted into 
!  the Riccati equation.
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    (a,b) - the interval on which to perform the procedure
!    fun - an external subroutine supplying the coefficients of the 
!      user's scalar system 
!    par? - user-supplied parameters which are passed to fun
!    chebsol - a structure describing the discretization scheme used to represent
!       the solution
!    rs - the values of the solution at the discretization nodes
!
!  Output parameters:
!    dres - a measure of the residual
!
!

double precision, allocatable :: ts(:),errs(:)
double complex                :: r, rp, rpp, qs(3), q0, q1, q2, ima


eps0  = epsilon(0.0d0)
ima   = (0.0d0,1.0d0)

n     = 3
nn    = 10000
a     = chebsol%a
b     = chebsol%b
allocate(ts(nn),errs(nn))

do i=1,nn
ts(i) = a+(b-a)*(i-1.0d0)/(nn-1.0d0)
end do

errs  = 1d300

do i=1,nn
t  = ts(i)
call fun(n,t,qs,pars)
q0 = qs(1)
q1 = qs(2)
q2 = qs(3)

dd = max(abs(q0),abs(q1))
dd = max(abs(q2),dd)

call chebpw_interp(chebsol,rs(:,2),t,r)
call chebpw_interp(chebsol,rs(:,3),t,rp)
call chebpw_interp(chebsol,rs(:,4),t,rpp)

errs(i)    = abs(rpp+3*r*rp+r**3+q2*rp+q2*r**2+q1*r+q0) / dd

end do

dres = maxval(errs)
end subroutine



subroutine scalar_ode3(n,t,y,f,df,userptr)
implicit double precision (a-h,p-z)
type(c_ptr)                    :: userptr
double complex                 :: y(:), f(:), df(:,:)
type(scalar_odedata), pointer  :: odedata
double complex                 :: us(3)
double complex                 :: q0, q1, q2, r, rp, rpp
double complex                 :: dk0, dk1, dk2

call c_f_pointer(userptr,odedata)

ifleft = odedata%ifleft
wa     = odedata%wa
wb     = odedata%wb
phi    = (erf(wa*(t-wb))+1)/2


call odedata%fun(n+1,t,us,odedata%pars)

r      = y(1)
rp     = y(2)


if (ifleft .eq. 1) then
dk0 = odedata%dks(1)
dk1 = odedata%dks(2)
dk2 = odedata%dks(3)
q0  = phi*dk0 + (1-phi)*us(1)
q1  = phi*dk1 + (1-phi)*us(2)
q2  = phi*dk2 + (1-phi)*us(3)
elseif (ifleft .eq. 0) then
dk0 = odedata%dks(1)
dk1 = odedata%dks(2)
dk2 = odedata%dks(3)
q0  = (1-phi)*dk0 + phi*us(1)
q1  = (1-phi)*dk1 + phi*us(2)
q2  = (1-phi)*dk2 + phi*us(3)
else
q0  = us(1)
q1  = us(2)
q2  = us(3)
endif



f(1)    = rp
f(2)    = - ( r**3 + 3*r*rp + q2*rp + q2*r**2 + q1*r + q0)


df(1,1) = 0
df(1,2) = 1

df(2,1) = -( 3*r**2 + 3*rp + 2*q2*r + q1 )
df(2,2) = -( 3*r + q2 )


end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Scalar equations of order 4
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine scalar_homotope4(vars,ier,eps,a,b,ifleft,fun,pars,rvals)
implicit double precision (a-h,o-z)
type(scalar_vars_t)              :: vars
procedure(scalar_fun)            :: fun
double complex, pointer          :: pars(:)
double complex                   :: rvals(:,:)
!
!  Use the homotopy method to compute the values of the derivatives of the r_j(t)
!  at a specified point.  More explicitly, this routine homotopes the coefficient 
!  matrix for the user's nth order scalar system to the constant system of the form
!
!          (  0     1    0    0  )
!  A0    = (  0     0    1    0  )
!          (  0     0    0    1  ) 
!          ( -dk0 -dk1 -dk2 -dk4 ) 
!
!  in order to calculate appropriate initial conditions for the user's system.  
!  Reasonable values of dk? are determined by calling the user-supplied external
!  subroutine.
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    eps - precision for the calculations
!    (a,b) - the interval on which to perform the procedure
!    ifleft - an integer parameter indicating whether the values of the r_j(t) should
!      be computed at the left-hand endpoint of the interval a or the right-hand
!      endpoint b;
!
!        ifleft = 1   means calculate the values at a
!        ifleft = 0   means calculate the values at b
!
!    fun - an external subroutine supplying the coefficients of the user's scalar
!      equation
!    par? - user-supplied parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   indicates the ODE solver failed
!       ier = 256 means the eigensolver failed
!
!    rvals - the (i,j) entry contains the value of r_j^(i)(c), where c=a if ifleft =1
!      and c=b if ifleft=0
!    

type(chebpw_scheme)               :: chebsol
type(c_ptr)                       :: userptr
type(scalar_odedata), pointer     :: odedata
double complex, allocatable       :: ys(:,:), yders(:,:)
double complex                    :: ima, yc(3)
double complex                    :: y0, y0p, amatr(4,4),u(4,4), w(4), dks(4)
! double complex, allocatable       :: work(:)
! double precision, allocatable     :: rwork(:)

double precision                  :: rcoefs(4), icoefs(4), reigs(4), ieigs(4)
integer                           :: its(4), flag

ier  = 0

k    = vars%kcheb
eps0 = epsilon(0.0d0)


n    = 4
wa   = 11d0/(b-a)
wb   = (a+b)/2.0d0
ima  = (0.0d0,1.0d0)

call fun(n,(a+b)/2.0d0,dks,pars)


n    = 4
allocate(odedata)
allocate(odedata%dks(n))

odedata%n       = n-1
odedata%pars   => pars
odedata%dks     = dks
odedata%fun    => fun
odedata%wa      = wa
odedata%wb      = wb
odedata%ifleft  = ifleft
userptr         = c_loc(odedata)


if (ifleft .eq. 0) then
c            = a
else
c            = b
endif

its       = 200

rcoefs(1) = real(dks(4))
icoefs(1) = imag(dks(4))

rcoefs(2) = real(dks(3))
icoefs(2) = imag(dks(3))

rcoefs(3) = real(dks(2))
icoefs(3) = imag(dks(2))

rcoefs(4) = real(dks(1))
icoefs(4) = imag(dks(1))

call init_random_seed()
call zamvw(n,rcoefs,icoefs,reigs,ieigs,its,flag)
if (flag .ne. 0) then
ier = 256
return
endif

w(1) = reigs(1) + ima*ieigs(1)
w(2) = reigs(2) + ima*ieigs(2)
w(3) = reigs(3) + ima*ieigs(3)
w(4) = reigs(4) + ima*ieigs(4)

!call scalar_sortw(4,w)
!call prinz("in scalar_homotopy4, ws = ",w)


do i=1,n
yc(1)      = w(i)
yc(2)      = 0
yc(3)      = 0
ifit       = 2

call odesolve_nonlinear1(vars%odesolve_vars,jer,ifit,eps,a,b,c,yc,n-1,scalar_ode4,userptr,chebsol,ys,yders)
if (jer .ne. 0) then
ier = 4
return
endif

nn = chebsol%nints*k

if (ifleft .eq. 1) then
rvals(1,i)  = ys(1,1)
rvals(2,i)  = ys(1,2)
rvals(3,i)  = ys(1,3)
else
rvals(1,i)  = ys(nn,1)
rvals(2,i)  = ys(nn,2)
rvals(3,i)  = ys(nn,3)
endif

end do

end subroutine


subroutine scalar_levin4(vars,ier,eps,a,b,ifleft,fun,pars,rvals)
implicit double precision (a-h,o-z)
type(scalar_vars_t)      :: vars
procedure(scalar_fun)    :: fun
double complex, pointer  :: pars(:)
double complex           :: rvals(:,:)
!
!  Use the Levin method to compute the values of the derivatives of the r_j(t)
!  at a specified point.
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    eps - precision for the calculations
!    (a,b) - the interval on which to perform the procedure
!    ifleft - an integer parameter indicating whether the values of the r_j(t) should
!      be computed at the left-hand endpoint of the interval a or the right-hand
!      endpoint b;
!
!        ifleft = 1   means calculate the values at a
!        ifleft = 0   means calculate the values at b
!
!    fun - an external subroutine supplying the coefficients of the 
!      user's scalar system 
!    par? - user-supplied parameters which are passed to fun
!
!  Output parameters:
!    ierc - an error return code;
!       ier = 0    indicates successful execution
!       ier = 4    indicates the ODE solver failed
!       ier = 256  indicates that the eigensolver failed
!       ier = 1024 means that the coefficients are not resolved by the Chebyshev
!                  discretization
!

double precision, allocatable            :: ts(:)
double complex, allocatable              :: q0s(:), q1s(:), q2s(:), q3s(:), qs(:), coefs(:)
double precision                         :: rcoefs(4), icoefs(4), reigs(4), ieigs(4) 
integer                                  :: its(4), flag

double complex, allocatable              :: r(:), rp(:), rpp(:), rppp(:), rhs(:), amatr(:,:), delta(:)
double complex                           :: ima, rval, rder, w(4)
double complex                           :: rp0, rpp0

double complex, allocatable              :: ws(:,:), roots(:)

kcheb    = vars%kcheb
ntest    = vars%ntest
maxiters = vars%maxleviters


ier      = 0
ima      = (0.0d0,1.0d0)
eps0     = epsilon(0.0d0)
epssol   = eps0*10
epsnewt  = eps0*1000
pi       = acos(-1.0d0)

n        = 4
allocate( ts(kcheb), q0s(kcheb), q1s(kcheb), q2s(kcheb),  q3s(kcheb), qs(n), coefs(kcheb)  )
allocate( r(kcheb), rp(kcheb), rpp(kcheb), rppp(kcheb), rhs(kcheb), amatr(kcheb,kcheb), delta(kcheb) )

!allocate( rs(kcheb), rders(kcheb), delta(kcheb), rhs(kcheb), amatr(kcheb,kcheb) )

ts = (b-a)/2*vars%xscheb + (b+a)/2

do i=1,kcheb
call fun(n,ts(i),qs,pars)
q0s(i) = qs(1)
q1s(i) = qs(2)
q2s(i) = qs(3)
q3s(i) = qs(4)
end do



! 
!  Test the coefficinets for fit
!

coefs = matmul(vars%acoefs,q0s)
dd1   = maxval(abs(coefs))
dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
if (dd1 .eq. 0) dd1 = 1
dd    = dd2/dd1
if (dd .gt. eps) then
ier = 1024
return
endif

coefs = matmul(vars%acoefs,q1s)
dd1   = maxval(abs(coefs))
dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
if (dd1 .eq. 0) dd1 = 1
dd    = dd2/dd1
if (dd .gt. eps) then
ier = 1024
return
endif

coefs = matmul(vars%acoefs,q2s)
dd1   = maxval(abs(coefs))
dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
if (dd1 .eq. 0) dd1 = 1
dd    = dd2/dd1
if (dd .gt. eps) then
ier = 1024
return
endif

coefs = matmul(vars%acoefs,q3s)
dd1   = maxval(abs(coefs))
dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
if (dd1 .eq. 0) dd1 = 1
dd    = dd2/dd1
if (dd .gt. eps) then
ier = 1024
return
endif

!
!  Use the eigensolver
!

call init_random_seed()
allocate(ws(kcheb,4))

n          = 4
! its        = 20000

! do i=1,kcheb
!    rcoefs(1) = real(q3s(i))
!    icoefs(1) = imag(q3s(i))

!    rcoefs(2) = real(q2s(i))
!    icoefs(2) = imag(q2s(i))

!    rcoefs(3) = real(q1s(i))
!    icoefs(3) = imag(q1s(i))

!    rcoefs(4) = real(q0s(i))
!    icoefs(4) = imag(q0s(i))



!    call zamvw(n,rcoefs,icoefs,reigs,ieigs,its,flag)
!    if (flag .ne. 0) then
!    ier = 256
!    return
!    endif

!    w(1) = reigs(1) + ima*ieigs(1)
!    w(2) = reigs(2) + ima*ieigs(2)
!    w(3) = reigs(3) + ima*ieigs(3)
!    w(4) = reigs(4) + ima*ieigs(4)


!    if (i .eq. 1) then
!    call scalar_sortw(4,w)
!    idx1 = 1
!    idx2 = 2
!    idx3 = 3
!    idx4 = 4
!    else

!    idx1 = minloc(abs(w(1)-ws(i-1,:)),1)
!    idx2 = minloc(abs(w(2)-ws(i-1,:)),1)
!    idx3 = minloc(abs(w(3)-ws(i-1,:)),1)
!    idx4 = minloc(abs(w(4)-ws(i-1,:)),1)
!    endif



!    ws(i,idx1) = w(1)
!    ws(i,idx2) = w(2)
!    ws(i,idx3) = w(3)
!    ws(i,idx4) = w(4)
! end do


dk = real(pars(1))
call root_estimates(4, dk, roots)


do isol=1,n

dd   = 2/(b-a)
!r    = ws(:,isol)
r = roots(isol)




rp   = dd * matmul(vars%adiff,r)
rpp  = dd**2 * matmul(vars%adiff2,r)
rppp = dd**3 * matmul(vars%adiff3,r)


do iter=1,maxiters
dd     = 2/(b-a)
ifdone = 0

rhs     = -(q0s+r**4+6*r**2*rp+3*rp**2+4*r*rpp+rppp+q2s*r**2+q2s*rp+q1s*r+ q3s*rpp+q3s*r**3+3*q3s*r*rp)
amatr   = dd**3 * vars%adiff3

do i=1,kcheb
   amatr(i,i) = amatr(i,i) + (4*rpp(i) + 12*r(i)*rp(i) + 4*r(i)**3 + 2*q2s(i)*r(i) + q1s(i) + 3*q3s(i)*r(i)**2 + 3*q3s(i)*rp(i))
end do


do i=1,kcheb
   amatr(i,:) = amatr(i,:) + (6*rp(i)+6*r(i)**2+q2s(i)+3*q3s(i)*r(i))*dd*vars%adiff(i,:)
end do


do i=1,kcheb
amatr(i,:) = amatr(i,:) + (4*r(i)+q3s(i))*dd**2*vars%adiff2(i,:)
end do

! if (isol == 3) then
!    call prinz("", amatr)
!    stop 
! end if

call linalg0_qrtsolve(epssol,kcheb,kcheb,amatr,delta,rhs)




! perform at least one iteration
if (iter .gt. 1) then
dd1 = norm2(abs(delta))
dd0 = norm2(abs(r))+1.0d0
dd2  = dd1/dd0
if (dd2 .lt. epsnewt) ifdone = 1
endif

r    = r + delta
rp   = dd * matmul(vars%adiff,r)
rpp  = dd * matmul(vars%adiff,rp)
rppp = dd * matmul(vars%adiff,rpp)


!print *,iter,norm2(abs(delta))/norm2(abs(r))

if (ifdone .eq. 1) exit

end do





if (iter .gt. maxiters) then
ier = 4
return
endif



if (ifleft .eq. 1) then
rvals(1,isol) = r(1)
rvals(2,isol) = rp(1)
rvals(3,isol) = rpp(1)
else
rvals(1,isol) = r(kcheb)
rvals(2,isol) = rp(kcheb)
rvals(3,isol) = rpp(kcheb)
endif


end do

end subroutine


subroutine scalar_riccati4(vars,ier,eps,a,b,c,rvals,fun,pars,chebsol,rs)
implicit double precision (a-h,o-z)
type(scalar_vars_t)                         :: vars
procedure(scalar_fun)                       :: fun
double complex, pointer                     :: pars(:)
double complex                              :: rvals(:)
type(chebpw_scheme)                         :: chebsol
double complex, allocatable, intent(out)    :: rs(:,:)

!
!  Solve the Riccati equation corresponding to (1) over the interval [a,b] 
!  given the values of r_j(t) and its first 2 derivatives at a point t in [a,b].
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    eps - precision for the Newton iteration
!    (a,b) - the interval on which to perform the procedure
!    c - the point at which the values of r_j(t) and its derivatives are given
!    rvals - an array of length 3 which gives the values of r_j(c), r_j'(c), and r_j''(c)
!    fun - an external subroutine supplying the coefficients of the 
!      user's scalar system 
!    par? - user-supplied parameters which are passed to fun
!
!  Output parameters:
!    ier  - a error return code;
!      ier = 0  indicates successful execution
!      ier = 4  the ODE solver failed
!
!    chebsol - a structure describing the discretization scheme used to represent
!       the solution
!    rs - the values of the solution at the discretization nodes
!
!
type(c_ptr)                                :: userptr
type(scalar_odedata), pointer              :: odedata
double complex, allocatable                :: ys(:,:), yders(:,:), rr(:)
double complex                             :: ima


ier   = 0
eps0  = epsilon(0.0d0)
ima   = (0.0d0,1.0d0)

allocate(odedata)
userptr = c_loc(odedata)

n              = 4
odedata%n      = n-1
odedata%fun   => fun
odedata%pars  => pars
odedata%wa     = 0
odedata%wb     = 0
odedata%ifleft = -1
ifit           = 2


call odesolve_nonlinear1(vars%odesolve_vars,jer,ifit,eps,a,b,c,rvals(2:n),n-1,scalar_ode4,userptr,chebsol,ys,yders)




if (jer .ne. 0) then
ier = 4
return
endif



call chebpw_info(chebsol,k,nints)
allocate(rs(k*nints,5))

rs(:,2) = ys(:,1)
rs(:,3) = ys(:,2)
rs(:,4) = ys(:,3)
rs(:,5) = yders(:,3)


call chebpw_int(chebsol,rs(:,2),c,rvals(1),rr)
rs(:,1) = rr

end subroutine



subroutine scalar_levinode4(vars,ier,eps,a,b,c,rvals,fun,pars,chebsol,rs)
implicit double precision (a-h,p-z)
type(scalar_vars_t)                         :: vars
type(chebpw_scheme), intent(out)            :: chebsol
double complex                              :: rvals(:)
double complex, allocatable, intent(out)    :: rs(:,:,:)
procedure(scalar_fun)                       :: fun
double complex, pointer                     :: pars(:)
!
!  Use the global Levin method to compute a set of phase functions.
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    eps - precision for the Newton iteration
!    (a,b) - the interval on which to perform the procedure
!    c - the point at which the values of r_j(t) and its derivatives are given
!    rvals - the desired values of the r_j
!    fun - an external subroutine supplying the coefficients of the 
!      user's scalar system 
!    par? - user-supplied parameters which are passed to fun
!  
!
!  Output parameters:
!    ier  - a error return code;
!      ier = 0  indicates successful execution
!      ier = 4  the Newton iterations did not converge
!
!    chebsol - a structure describing the discretization scheme used to represent
!      the solution
!    rs - a (*,n,n) such that r(:,i,j) gives the values of r_j^(i-1) at the 
!      discretization nodes
!
!

double precision, allocatable               :: stack(:,:), ts(:)
double complex, allocatable                 :: qs(:,:), r(:), rp(:), rpp(:), coefs(:), ws(:,:)
double complex, allocatable                 :: delta(:), rhs(:), amatr(:,:), rppp(:)

double complex, allocatable                 :: rs0(:,:,:,:), rvals0(:)

double precision, allocatable               :: ab(:,:)

double complex                              :: ima, qs0(4), w(4)
double precision                            :: rcoefs(4), icoefs(4), ieigs(4), reigs(4), drs(4)
integer                                     :: its(4), flag

n        = 4
ima      = (0.0d0,1.0d0)
eps0     = epsilon(0.0d0)
epssol   = eps0*10

ier      = 0
ind      = 0
nints    = 0
kcheb    = vars%kcheb
ntest    = vars%ntest
maxiters = vars%maxleviters
maxstack = vars%maxstack
maxints  = vars%maxints


allocate(ts(kcheb), qs(kcheb,4), r(kcheb), rp(kcheb), rpp(kcheb), coefs(kcheb) )
allocate(delta(kcheb), rhs(kcheb), amatr(kcheb,kcheb), rppp(kcheb))
allocate(rs0(kcheb,4,4,maxints), stack(2,maxstack), ab(2,maxints))

call init_random_seed()
allocate(ws(kcheb,4))


nstack     = 1
stack(1,1) = a
stack(2,1) = b

do while(nstack > 0)

   !
   !  Pop an interval off of the stack
   !
    
   a0     = stack(1,nstack)
   b0     = stack(2,nstack)
   nstack = nstack-1
   ind    = ind+1

   ifsplit = 1
   dcoefs  = 1d300
   drs     = 1d300
   jer     = 0

   !
   !  Fetch the values of the coefficients
   !
    
   ts = vars%xscheb*(b0-a0)/2 + (b0+a0)/2
   do i=1,kcheb
   call fun(n,ts(i),qs0,pars)
   qs(i,1) = qs0(1)
   qs(i,2) = qs0(2)
   qs(i,3) = qs0(3)
   qs(i,4) = qs0(4)
   end do    
    
   !
   !  Check that the coefficients are properly discretized
   !

   dcoefs = 0
   do i=1,4
   coefs  = matmul(vars%acoefs,qs(:,i))
   dd1    = maxval(abs(coefs))+1
   dd2    = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
   dd     = dd2/dd1
   dcoefs = max(dd,dcoefs)
   end do

   !
   !  Compute the eigenvalues to use as inital guesses if this is the first interval
   !

   if (dcoefs .lt. eps .AND. nints .eq. 0) then
    
   n          = 4
   its        = 200
    
   do i=1,kcheb
    
   rcoefs(1) = real(qs(i,4))
   icoefs(1) = imag(qs(i,4))
    
   rcoefs(2) = real(qs(i,3))
   icoefs(2) = imag(qs(i,3))
    
   rcoefs(3) = real(qs(i,2))
   icoefs(3) = imag(qs(i,2))
    
   rcoefs(4) = real(qs(i,1))
   icoefs(4) = imag(qs(i,1))    
    
   call zamvw(n,rcoefs,icoefs,reigs,ieigs,its,flag)
   if (flag .ne. 0) then
   ier = 256
   return
   endif
    
   w(1) = reigs(1) + ima*ieigs(1)
   w(2) = reigs(2) + ima*ieigs(2)
   w(3) = reigs(3) + ima*ieigs(3)
   w(4) = reigs(4) + ima*ieigs(4)
    
    
    
   if (i .eq. 1) then
   call scalar_sortw(4,w)
    
   idx1 = 1
   idx2 = 2
   idx3 = 3
   idx4 = 4
   else
    
   idx1 = minloc(abs(w(1)-ws(i-1,:)),1)
   idx2 = minloc(abs(w(2)-ws(i-1,:)),1)
   idx3 = minloc(abs(w(3)-ws(i-1,:)),1)
   idx4 = minloc(abs(w(4)-ws(i-1,:)),1)
    
   endif
    
   ws(i,idx1) = w(1)
   ws(i,idx2) = w(2)
   ws(i,idx3) = w(3)
   ws(i,idx4) = w(4)
    
   end do
    
    
   endif
   

   !
   !  If the coefficients are properly discretized, proceed
   !
    
   if (jer .eq. 0 .AND. dcoefs .lt. eps) then
         
      iffail  = 0

      do isol=1,4
         
         if (nints .eq. 0 ) then    
         dd   = 2/(b0-a0)
         r    = ws(:,isol)
         rp   = dd * matmul(vars%adiff,r)
         rpp  = dd**2 * matmul(vars%adiff2,r)
         rppp = dd**3 * matmul(vars%adiff3,r)
         else
         r      = rs0(kcheb,1,isol,nints)
         rp     = rs0(kcheb,2,isol,nints)
         rpp    = rs0(kcheb,3,isol,nints)
         rppp   = rs0(kcheb,4,isol,nints)
         endif

         ifdone = 0

         do iter=1,maxiters

         dd      = 2/(b0-a0)
          
         rhs     = -(qs(:,1)+r**4+6*r**2*rp+3*rp**2+4*r*rpp+rppp+qs(:,3)*r**2+&
                    qs(:,3)*rp+qs(:,2)*r+ qs(:,4)*rpp+qs(:,4)*r**3+3*qs(:,4)*r*rp)

         amatr   = dd**3 * vars%adiff3
          
         do i=1,kcheb
         amatr(i,:) = amatr(i,:) + (4*r(i)+qs(i,4))*dd**2*vars%adiff2(i,:)
         end do
          
         do i=1,kcheb
         amatr(i,:) = amatr(i,:) + (6*rp(i)+6*r(i)**2+qs(i,3)+&
                                   3*qs(i,4)*r(i))*dd*vars%adiff(i,:)
         end do
          
         do i=1,kcheb
         amatr(i,i) = amatr(i,i) + (4*rpp(i) + 12*r(i)*rp(i) + 4*r(i)**3 + &
                      2*qs(i,3)*r(i) + qs(i,2) + 3*qs(i,4)*r(i)**2 + 3*qs(i,4)*rp(i))
         end do
          
         call linalg0_qrtsolve(epssol,kcheb,kcheb,amatr,delta,rhs)

         
         !perform at least one iteration
         if (iter .gt. 1) then
         dd1 = norm2(abs(delta))
         dd0 = norm2(abs(r))+1.0d0
         dd2  = dd1/dd0
         if (dd2 .lt. eps) ifdone=1
         endif
         
         r     = r + delta
         rp    = dd * matmul(vars%adiff,r)
         rpp   = dd * matmul(vars%adiff,rp)
         rppp  = dd * matmul(vars%adiff,rpp)

        ! write(*,"(3X,3(I3,1X),D15.7,2(I3,1X),1(D15.7,1X),(I3,1X))") ind,nints,jer,dcoefs,isol,iter,dd2,ifdone
        ! write(13,"(3X,3(I3,1X),D15.7,2(I3,1X),1(D15.7,1X),(I3,1X))") ind,nints,jer,dcoefs,isol,iter,dd2,ifdone

         if (ifdone .eq. 1) exit

         end do
         
       
         if (iter .gt. maxiters) then
         iffail = 1
         exit
         endif

         rs0(:,1,isol,nints+1) = r
         rs0(:,2,isol,nints+1) = rp
         rs0(:,3,isol,nints+1) = rpp
         rs0(:,4,isol,nints+1) = rppp

         coefs     = matmul(vars%acoefs,r)
         dd1       = maxval(abs(coefs))+1
         dd2       = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
         dd        = dd2/dd1
         drs(isol) = max(dd,dcoefs)
              
      end do

      if (iffail .eq. 0 .AND. maxval(drs) .lt. eps) ifsplit = 0
    
   endif

  ! write(*,"(3(I3,1X),D15.7,2(I3,1X),4(D15.7,1X),2(I3,1X))")  ind,nints,jer,dcoefs,isol,iter,drs,iffail,ifsplit
  ! write(13,"(3(I3,1X),D15.7,2(I3,1X),4(D15.7,1X),2(I3,1X))") ind,nints,jer,dcoefs,isol,iter,drs,iffail,ifsplit

   !
   !  Either split the interval or add it to the list of "accepted" intervals
   !
   if (ifsplit .eq. 1) then
   if (nstack +2 .ge. maxstack) then
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
   if (nints+1 .ge. maxints) then
   ier = 4
   return
   endif

   nints       = nints+1
   ab(1,nints) = a0
   ab(2,nints) = b0

   endif

end do

!
!  Copy out the solution
!

call chebpw_specified(chebsol,kcheb,nints,ab)

allocate(rs(nints*kcheb,5,4))

i1 = 1
i2 = kcheb

do int=1,nints
a0  = ab(1,int)
b0  = ab(2,int)

rs(i1:i2,2,1) = rs0(1:kcheb,1,1,int)
rs(i1:i2,3,1) = rs0(1:kcheb,2,1,int)
rs(i1:i2,4,1) = rs0(1:kcheb,3,1,int)
rs(i1:i2,5,1) = rs0(1:kcheb,4,1,int)

rs(i1:i2,2,2) = rs0(1:kcheb,1,2,int)
rs(i1:i2,3,2) = rs0(1:kcheb,2,2,int)
rs(i1:i2,4,2) = rs0(1:kcheb,3,2,int)
rs(i1:i2,5,2) = rs0(1:kcheb,4,2,int)

rs(i1:i2,2,3) = rs0(1:kcheb,1,3,int)
rs(i1:i2,3,3) = rs0(1:kcheb,2,3,int)
rs(i1:i2,4,3) = rs0(1:kcheb,3,3,int)
rs(i1:i2,5,3) = rs0(1:kcheb,4,3,int)

rs(i1:i2,2,4) = rs0(1:kcheb,1,4,int)
rs(i1:i2,3,4) = rs0(1:kcheb,2,4,int)
rs(i1:i2,4,4) = rs0(1:kcheb,3,4,int)
rs(i1:i2,5,4) = rs0(1:kcheb,4,4,int)

i1 = i1+kcheb
i2 = i2+kcheb

end do


call chebpw_int2(chebsol,rs(:,2,1),c,rvals(1),rvals0)
rs(:,1,1) = rvals0

call chebpw_int2(chebsol,rs(:,2,2),c,rvals(2),rvals0)
rs(:,1,2) = rvals0

call chebpw_int2(chebsol,rs(:,2,3),c,rvals(3),rvals0)
rs(:,1,3) = rvals0

call chebpw_int2(chebsol,rs(:,2,4),c,rvals(4),rvals0)
rs(:,1,4) = rvals0

end subroutine



subroutine scalar_ode4(n,t,y,f,df,userptr)
implicit double precision (a-h,p-z)
type(c_ptr)                    :: userptr
double complex                 :: y(:), f(:), df(:,:)
type(scalar_odedata), pointer  :: odedata
double complex                 :: us(4)
double complex                 :: r, rp, rpp
double complex                 :: q0, q1, q2, q3
double complex                 :: dk0, dk1, dk2, dk3

call c_f_pointer(userptr,odedata)

ifleft = odedata%ifleft
wa     = odedata%wa
wb     = odedata%wb
phi    = (erf(wa*(t-wb))+1)/2

call odedata%fun(n+1,t,us,odedata%pars)

r      = y(1)
rp     = y(2)
rpp    = y(3)



if (ifleft .eq. 1) then
dk0 = odedata%dks(1)
dk1 = odedata%dks(2)
dk2 = odedata%dks(3)
dk3 = odedata%dks(4)

q0  = phi*dk0 + (1-phi)*us(1)
q1  = phi*dk1 + (1-phi)*us(2)
q2  = phi*dk2 + (1-phi)*us(3)
q3  = phi*dk3 + (1-phi)*us(4)

elseif (ifleft .eq. 0) then
dk0 = odedata%dks(1)
dk1 = odedata%dks(2)
dk2 = odedata%dks(3)
dk3 = odedata%dks(4)

q0  = (1-phi)*dk0 + phi*us(1)
q1  = (1-phi)*dk1 + phi*us(2)
q2  = (1-phi)*dk2 + phi*us(3)
q3  = (1-phi)*dk3 + phi*us(4)

else
q0  = us(1)
q1  = us(2)
q2  = us(3)
q3  = us(4)
endif

f(1)    = rp
f(2)    = rpp
f(3)    = -q0 -r**4-6*r**2*rp -3*rp**2 -4*r*rpp -q1*r - q2*r**2 - q2*rp -q3*r**3-3*q3*r*rp-q3*rpp



df(1,1) = 0
df(1,2) = 1
df(1,3) = 0

df(2,1) = 0
df(2,2) = 0
df(2,3) = 1

df(3,1) = -4*r**3 - 12*r*rp - 4*rpp - q1 -2*q2*r - 3*q3*r**2 -3*q3*rp
df(3,2) = -6*r**2 - 6*rp -q2 -3*q3*r
df(3,3) = -4*r-q3

end subroutine



subroutine scalar_sortw(n,w)
implicit double precision (a-h,o-z)
double complex               :: w(:), w2(n), u1, u2
double precision             :: v(n)
integer                      :: idxs(n)

eps = epsilon(0.0d0)*10

do i=1,n
v(i)    = abs(w(i))
idxs(i) = i
end do

call insort2(n,v,idxs)

w2 = w(idxs)
w  = w2


do i=1,n-1
dd1 = abs(w(i))
dd2 = abs(w(i+1))

if( (dd1-dd2)/dd1 .lt. eps) then

if (real(w(i)) .lt. real(w(i+1))) then
u1     = w(i)
u2     = w(i+1)
w(i)   = u2
w(i+1) = u1
endif

endif

end do

end subroutine



subroutine scalar_residual4(vars,fun,pars,chebsol,rs,dres)
implicit double precision (a-h,o-z)
type(scalar_vars_t)                         :: vars
procedure(scalar_fun)                       :: fun
double complex, pointer                     :: pars(:)
type(chebpw_scheme)                         :: chebsol
double complex                              :: rs(:,:)
!
!  Return a measure of the error obtained when a solution is inserted into 
!  the Riccati equation.
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    (a,b) - the interval on which to perform the procedure
!    fun - an external subroutine supplying the coefficients of the 
!      user's scalar system 
!    par? - user-supplied parameters which are passed to fun
!    chebsol - a structure describing the discretization scheme used to represent
!       the solution
!    rs - the values of the solution at the discretization nodes
!
!  Output parameters:
!    dres - a measure of the residual
!
!

double precision, allocatable :: ts(:),errs(:)
double complex                :: r, rp, rpp, rppp, qs(4), q0, q1, q2, q3, ima


eps0  = epsilon(0.0d0)
ima   = (0.0d0,1.0d0)

n     = 4
nn    = 10000
a     = chebsol%a
b     = chebsol%b
allocate(ts(nn),errs(nn))

do i=1,nn
ts(i) = a+(b-a)*(i-1.0d0)/(nn-1.0d0)
end do

errs  = 1d300

do i=1,nn
t  = ts(i)
call fun(n,t,qs,pars)
q0 = qs(1)
q1 = qs(2)
q2 = qs(3)
q3 = qs(4)

call chebpw_interp(chebsol,rs(:,2),t,r)
call chebpw_interp(chebsol,rs(:,3),t,rp)
call chebpw_interp(chebsol,rs(:,4),t,rpp)
call chebpw_interp(chebsol,rs(:,5),t,rppp)

errs(i)     = abs(q0+r**4+6*r**2*rp+3*rp**2+4*r*rpp+rppp+q2*r**2+q2*rp+q1*r+ q3*rpp+q3*r**3+3*q3*r*rp)

end do

dres = maxval(errs)

end subroutine





! subroutine scalar_levin5(vars,ier,eps,a,b,ifleft,fun,pars,rvals)
!    implicit double precision (a-h,o-z)
!    type(scalar_vars_t)      :: vars
!    procedure(scalar_fun)    :: fun
!    double complex, pointer  :: pars(:)
!    double complex           :: rvals(:,:)
!    !
!    !  Use the Levin method to compute the derivatives of the r_j(t) at a specified 
!    !  point.
!    !
!    !  Input parameters:
!    !    vars - the structure prepared by scalar_init
!    !    eps - precision for the Newton iteration
!    !    (a,b) - the interval on which to perform the procedure
!    !    ifleft - an integer parameter indicating whether the values of the r_j(t) should
!    !      be computed at the left-hand endpoint of the interval a or the right-hand
!    !      endpoint b;
!    !
!    !        ifleft = 1   means calculate the values at a
!    !        ifleft = 0   means calculate the values at b
!    !
!    !    fun - an external subroutine supplying the coefficients of the 
!    !      user's scalar system 
!    !    par? - user-supplied parameters which are passed to fun
!    !
!    !  Output parameters:
!    !    ier - an error return code;
!    !       ier = 0     indicates successful execution
!    !       ier = 4     indicates the ODE solver failed
!    !       ier = 256   means that the eigensolver failed
!    !       ier = 1024  indicates the coefficients were not represented resolved
!    !                   by Chebyshev expansions of order vars%kcheb over [a,b]
!    !
!    !    rvals - the (i,j) entry contains the value of r_j^(i)(c), where c=a if 
!    !      ifleft =1 and c=b if ifleft=0
!    !

!    double precision, allocatable            :: ts(:)
!    double precision                         :: rcoefs(5), icoefs(5), reigs(5), ieigs(5)
!    integer                                  :: its(5), flag

!    double complex, allocatable              :: q0s(:), q1s(:), q2s(:), qs(:), coefs(:), ws(:,:)

!    double complex, allocatable              :: r(:), rp(:), rpp(:), rhs(:), amatr(:,:), delta(:)
!    double complex                           :: ima, rval, rder, w(5)
!    double complex                           :: rp0, rpp0

!    kcheb    = vars%kcheb
!    ntest    = vars%ntest
!    maxiters = vars%maxleviters

!    ier      = 0
!    ima      = (0.0d0,1.0d0)
!    eps0     = epsilon(0.0d0)
!    epssol   = eps0*10

!    n        = 2
!    allocate( ts(kcheb), q0s(kcheb), q1s(kcheb), q2s(kcheb), q3s(kcheb), q4s(kcheb), qs(n), coefs(kcheb)  )
!    allocate(r(kcheb), rp(kcheb), rpp(kcheb), rppp(kcheb), rpppp(kcheb) rhs(kcheb), amatr(kcheb,kcheb), delta(kcheb))

!    ts = (b-a)/2*vars%xscheb + (b+a)/2

!    do i=1,kcheb
!       call fun(n,ts(i),qs,pars)
!       q0s(i) = qs(1)
!       q1s(i) = qs(2)
!       q2s(i) = qs(3)
!       q3s(i) = qs(4)
!       q4s(i) = qs(5)
!    end do

!    ! 
!    !  Test the coefficinets for fit
!    !

!    coefs = matmul(vars%acoefs,q0s)
!    dd1   = maxval(abs(coefs))
!    dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
!    if (dd1 .eq. 0) dd1 = 1
!    dd    = dd2/dd1
!    if (dd .gt. eps) then
!       ier = 1024
!       return
!    endif

!    coefs = matmul(vars%acoefs,q1s)
!    dd1   = maxval(abs(coefs))
!    dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
!    if (dd1 .eq. 0) dd1 = 1
!    dd    = dd2/dd1
!    if (dd .gt. eps) then
!       ier = 1024
!       return
!    endif

!    coefs = matmul(vars%acoefs,q2s)
!    dd1   = maxval(abs(coefs))
!    dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
!    if (dd1 .eq. 0) dd1 = 1
!    dd    = dd2/dd1
!    if (dd .gt. eps) then
!       ier = 1024
!       return
!    endif

!    coefs = matmul(vars%acoefs,q3s)
!    dd1   = maxval(abs(coefs))
!    dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
!    if (dd1 .eq. 0) dd1 = 1
!    dd    = dd2/dd1
!    if (dd .gt. eps) then
!       ier = 1024
!       return
!    endif

!    coefs = matmul(vars%acoefs,q4s)
!    dd1   = maxval(abs(coefs))
!    dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))
!    if (dd1 .eq. 0) dd1 = 1
!    dd    = dd2/dd1
!    if (dd .gt. eps) then
!       ier = 1024
!       return
!    endif



!    !
!    !  Use the eigensolver
!    !

!    call init_random_seed()
!    allocate(ws(kcheb,2))

!    n          = 5
!    its        = 2000

!    do i=1,kcheb

!       rcoefs(1) = real(q4s(i))
!       icoefs(1) = imag(q4s(i))

!       rcoefs(2) = real(q3s(i))
!       icoefs(2) = imag(q3s(i))

!       rcoefs(3) = real(q2s(i))
!       icoefs(3) = imag(q2s(i))

!       rcoefs(4) = real(q1s(i))
!       icoefs(4) = imag(q1s(i))

!       rcoefs(5) = real(q0s(i))
!       icoefs(5) = imag(q0s(i))

!       call zamvw(n,rcoefs,icoefs,reigs,ieigs,its,flag)
!       if (flag .ne. 0) then
!          ier = 256
!          return
!       endif

!       w(1) = reigs(1) + ima*ieigs(1)
!       w(2) = reigs(2) + ima*ieigs(2)
!       w(3) = reigs(3) + ima*ieigs(3)
!       w(4) = reigs(4) + ima*ieigs(4)
!       w(5) = reigs(5) + ima*ieigs(5)

!       if (i .eq. 1) then
!             call scalar_sortw(2,w)
!             idx1 = 1
!             idx2 = 2
!             idx3 = 3 
!             idx4 = 4
!             idx5 = 5
!          else
!             idx1 = minloc(abs(w(1)-ws(i-1,:)),1)
!             idx2 = minloc(abs(w(2)-ws(i-1,:)),1)
!             idx3 = minloc(abs(w(3)-ws(i-1,:)),1)
!             idx4 = minloc(abs(w(4)-ws(i-1,:)),1)
!             idx5 = minloc(abs(w(5)-ws(i-1,:)),1)
!       endif

!       ws(i,idx1) = w(1)
!       ws(i,idx2) = w(2)
!       ws(i,idx3) = w(3)
!       ws(i,idx4) = w(4)
!       WS(i,idx5) = w(5)
!    end do

!    do isol=1,n

!       dd    = 2/(b-a)
!       r     = ws(:,isol)
!       rp    = dd * matmul(vars%adiff,r)
!       rpp   = dd**2 * matmul(vars%adiff2,r)
!       rppp  = dd**3 * matmul(vars%adiff3,r)
!       rpppp = dd**4 * matmul(vars%adiff4,r)


!       do iter=1,maxiters
!          dd      = 2/(b-a)
!          ifdone  = 0
!          rhs     = -(rp+r**2+q1s*r+q0s)
!          amatr   = dd * vars%adiff

!          do i=1,kcheb
!             amatr(i,i) = amatr(i,i) + (2*r(i) + q1s(i))
!          end do

!          !call linalg0_utvsolve(epssol,kcheb,kcheb,amatr,delta,rhs)
!          call linalg0_qrtsolve(epssol,kcheb,kcheb,amatr,delta,rhs)

!          ! perform at least one iteration
!          if (iter .gt. 1) then
!          dd1 = norm2(abs(delta))
!          dd0 = norm2(abs(r))+1.0d0
!          dd2  = dd1/dd0
!          if (dd2 .lt. eps) ifdone = 1
!          endif

!          r    = r + delta
!          rp   = dd * matmul(vars%adiff,r)

!          ! print *,isol,iter,dd2,eps

!          if (ifdone .eq. 1) exit

!       end do


!       if (iter .gt. maxiters) then
!          ier = 4
!          return
!       endif


!       if (ifleft .eq. 1) then
!             rvals(1,isol) = r(1)
!          else
!             rvals(1,isol) = r(kcheb)
!       endif

!    end do

    
! end subroutine

















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Scalar equations of order 3 using specialized second order solver
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine scalar_homotope33(vars,ier,eps,a,b,ifleft,fun,par1,par2,par3, &
  r1p,r1pp,r2p,r2pp,r3p,r3pp)
implicit double precision (a-h,o-z)
type(scalar_vars_t)      :: vars
procedure(scalar_fun33)  :: fun
double complex           :: par1,par2,par3
double complex           :: r1p,r1pp
double complex           :: r2p,r2pp
double complex           :: r3p,r3pp
!
!  Use the homotopy method to compute the values of the derivatives of the r_j(t)
!  at a specified point.  More explicitly, this routine homotopes the coefficient 
!  matrix for the user's nth order scalar system to the constant system
!
!          (  0     1    0  )
!  A0    = (  0     0    1  )
!          ( -dk0 -dk1 -dk2 )
!
!  in order to calculate appropriate initial conditions for the user's system.
!
!  Be warned that if the constants dk0, dk1 and dk2 are poorly chosen, or the resulting
!  homotopy has a tunring point, this will likely fail.
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    eps - precision for the calculations
!    (a,b) - the interval on which to perform the procedure
!    ifleft - an integer parameter indicating whether the values of the r_j(t) should
!      be computed at the left-hand endpoint of the interval a or the right-hand
!      endpoint b;
!
!        ifleft = 1   means calculate the values at a
!        ifleft = 0   means calculate the values at b
!
!    fun - an external subroutine supplying the coefficients of the user's scalar
!      equation
!    par? - user-supplied parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   indicates the ODE solver failed
!       ier = 256 the eigensolver failed
!
!   r?p, r?pp - the values of r_?'(t) and r_?''(t) at the point a (ifleft = 1)
!    or b (ifleft = 0)
!
!    
type(chebpw_scheme)               :: chebsol
type(c_ptr)                       :: userptr
type(scalar_odedata33), pointer   :: odedata
double complex, allocatable       :: ys(:,:), yders(:,:)
double complex                    :: ima, yc(2)

double complex, allocatable       :: ys0(:), yders0(:), yder2s0(:)
double complex                    :: w(3), y0, y0p
double complex                    :: dk0, dk1, dk2

! double complex                    :: y0, y0p, amatr(3,3),u(3,3),w(3)
! double complex, allocatable       :: work(:)
! double precision, allocatable     :: rwork(:)

double precision                  :: icoefs(3), rcoefs(3), ieigs(3), reigs(3)
integer                           :: its(3), flag


ier  = 0

k    = vars%kcheb
eps0 = epsilon(0.0d0)

wa   = 11.0d0/(b-a)
wb   = (a+b)/2
ima  = (0.0d0,1.0d0)

call fun((a+b)/2,dk0,dk1,dk2,par1,par2,par3)

allocate(odedata)

odedata%par1    = par1
odedata%par2    = par2
odedata%par3    = par3
odedata%fun    => fun
odedata%dk0     = dk0
odedata%dk1     = dk1
odedata%dk2     = dk2
odedata%wa      = wa
odedata%wb      = wb
odedata%ifleft  = ifleft
userptr         = c_loc(odedata)

if (ifleft .eq. 1) then
c            = b
else
c            = a
endif

n          = 3
its        = 200

rcoefs(1) = real(dk2)
icoefs(1) = imag(dk2)

rcoefs(2) = real(dk1)
icoefs(2) = imag(dk1)

rcoefs(3) = real(dk0)
icoefs(3) = imag(dk0)

call prin2("rcoefs = ",rcoefs)
call prin2("icoefs = ",icoefs)

call init_random_seed()
call zamvw(n,rcoefs,icoefs,reigs,ieigs,its,flag)

call prin2("reigs = ",reigs)
call prin2("ieigs = ",ieigs)

if (flag .ne. 0) then
ier = 256
return
endif

w(1) = reigs(1) + ima*ieigs(1)
w(2) = reigs(2) + ima*ieigs(2)
w(3) = reigs(3) + ima*ieigs(3)

!call scalar_sortw(n,w)

call prinz("in scalar_homotopy3, ws = ",w)

y0         = w(1)
y0p        = 0
call odetwo_nonlinear1(vars%odetwo_vars,jer,eps,a,b,c,y0,y0p,chebsol,ys0,yders0,yder2s0, &
  scalar_ode33,userptr)
if (jer .ne. 0) then
ier = 4
return
endif

nn = chebsol%nints*k
if (ifleft .eq. 1) then
r1p  = ys0(1)
r1pp = yders0(1)
else
r1p  = ys0(nn)
r1pp = yders0(nn)
endif

call prini("in scalar_homotopy3, nints1 = ",chebsol%nints)


!!!!!!!!!!!!!!!!

y0         = w(2)
y0p        = 0
call odetwo_nonlinear1(vars%odetwo_vars,jer,eps,a,b,c,y0,y0p,chebsol,ys0,yders0,yder2s0, &
  scalar_ode33,userptr)
if (jer .ne. 0) then
ier = 4
return
endif

nn = chebsol%nints*k
if (ifleft .eq. 1) then
r2p  = ys0(1)
r2pp = yders0(1)
else
r2p  = ys0(nn)
r2pp = yders0(nn)
endif

call prini("in scalar_homotopy3, nints2 = ",chebsol%nints)

!!!!!!!!!!!!!!!!

y0         = w(3)
y0p        = 0
call odetwo_nonlinear1(vars%odetwo_vars,jer,eps,a,b,c,y0,y0p,chebsol,ys0,yders0,yder2s0, &
  scalar_ode33,userptr)
if (jer .ne. 0) then
ier = 4
return
endif

call prini("in scalar_homotopy3, nints3 = ",chebsol%nints)

nn = chebsol%nints*k
if (ifleft .eq. 1) then
r3p  = ys0(1)
r3pp = yders0(1)
else
r3p  = ys0(nn)
r3pp = yders0(nn)
endif

end subroutine


subroutine scalar_levin33(vars,ier,eps,a,b,ifleft,fun,par1,par2,par3, &
  r1p,r1pp,r2p,r2pp,r3p,r3pp)
implicit double precision (a-h,o-z)
type(scalar_vars_t)      :: vars
procedure(scalar_fun33)  :: fun
double complex           :: par1,par2,par3
double complex           :: r1p,r1pp
double complex           :: r2p,r2pp
double complex           :: r3p,r3pp
!
!  Use the Levin method to compute the values of the derivatives of the r_j(t)
!  at a specified point.
!
!  Input parameters:
!    vars - the structure prepared by scalar_init
!    eps - precision for the calculations
!    (a,b) - the interval on which to perform the procedure
!    ifleft - an integer parameter indicating whether the values of the r_j(t) should
!      be computed at the left-hand endpoint of the interval a or the right-hand
!      endpoint b;
!
!        ifleft = 1   means calculate the values at a
!        ifleft = 0   means calculate the values at b
!
!    fun - an external subroutine supplying the coefficients of the 
!      user's scalar system 
!    par? - user-supplied parameters which are passed to fun
!
!  Output parameters:
!    ier - an error return code;
!       ier = 0   indicates successful execution
!       ier = 4   indicates the ODE solver failed
!
!   r?p, r?pp - the values of r_?'(t) and r_?''(t) at the point a (ifleft = 1)
!    or b (ifleft = 0)
!


double precision, allocatable            :: ts(:)
double complex, allocatable              :: q0s(:), q1s(:), q2s(:)


double complex, allocatable              :: r(:), rp(:), rpp(:), rhs(:), amatr(:,:), delta(:)
double complex                           :: ima, rval, rder, bmatr(3,3), work(1000), w(3)
double complex                           :: rp0, rpp0

double precision                  :: icoefs(3), rcoefs(3), ieigs(3), reigs(3)
integer                           :: its(3), flag

kcheb    = vars%kcheb
maxiters = vars%maxleviters

ier      = 0
ima      = (0.0d0,1.0d0)
eps0     = epsilon(0.0d0)
epssol   = eps0*10
epsnewt  = eps0*100

n        = 3
allocate( ts(kcheb), q0s(kcheb), q1s(kcheb), q2s(kcheb) )
allocate( r(kcheb), rp(kcheb), rpp(kcheb), rhs(kcheb), amatr(kcheb,kcheb), delta(kcheb) )

! allocate(qs(kcheb))
! allocate(rs(kcheb), rders(kcheb), delta(kcheb), rhs(kcheb), amatr(kcheb,kcheb))

ts = (b-a)/2*vars%xscheb + (b+a)/2
do i=1,kcheb
call fun(ts(i),q0s(i),q1s(i),q2s(i),par1,par2,par3)
end do




rcoefs(1) = real(q2s(1))
icoefs(1) = imag(q2s(1))

rcoefs(2) = real(q1s(1))
icoefs(2) = imag(q1s(1))

rcoefs(3) = real(q0s(1))
icoefs(3) = imag(q0s(1))

call prin2("rcoefs = ",rcoefs)
call prin2("icoefs = ",icoefs)

call init_random_seed()
call zamvw(n,rcoefs,icoefs,reigs,ieigs,its,flag)
if (flag .ne. 0) then
ier = 256
return
endif

w(1) = reigs(1) + ima*ieigs(1)
w(2) = reigs(2) + ima*ieigs(2)
w(3) = reigs(3) + ima*ieigs(3)

!call scalar_sortw(n,w)

call prinz("in scalar_homotopy3, ws = ",w)

! n          = 3
! bmatr      = 0
! bmatr(2,1) = 1
! bmatr(3,2) = 1
! bmatr(1,3) = -q0s(1)
! bmatr(2,3) = -q1s(1)
! bmatr(3,3) = -q2s(1)

! call nonsym_eigvals(jer,bmatr,n,w,nroots,work)



do isol=1,3

r    = w(isol)
rp   = 0
rpp  = 0


do iter=1,maxiters
dd     = 2/(b-a)
ifdone = 0
rhs     = -(rpp+3*r*rp+r**3+q2s*rp+q2s*r**2+q1s*r+q0s)


amatr = dd**2 * vars%adiff2 

do i=1,kcheb
amatr(i,:) = amatr(i,:) + (3*r(i)+q2s(i))*dd*vars%adiff(i,:)
end do

do i=1,kcheb
amatr(i,i) = amatr(i,i) + 3*rp(i) + 3*r(i)**2 + 2*q2s(i)*r(i)+q1s(i)
end do


! !call linalg0_utvsolve(epssol,kcheb,kcheb,amatr,delta,rhs)
call linalg0_qrtsolve(epssol,kcheb,kcheb,amatr,delta,rhs)


! perform at least one iteration
if (iter .gt. 1) then
dd1 = norm2(abs(delta))
dd0 = norm2(abs(r))+1.0d0
dd2  = dd1/dd0
if (dd2 .lt. epsnewt) ifdone = 1
endif

r    = r + delta
rp   = dd * matmul(vars%adiff,r)
rpp  = dd * matmul(vars%adiff,rp)

!print *,"   ",norm2(abs(rhs)),norm2(abs(delta)), norm2(abs(r))

if (ifdone .eq. 1) exit

end do


if (iter .gt. maxiters) then
ier = 4
return
endif

if (ifleft .eq. 1) then
rp0  = r(1)
rpp0 = rp(1)
else
rp0  = r(kcheb)
rpp0 = rp(kcheb)
endif


if (isol .eq. 1) then
r1p  = rp0
r1pp = rpp0
elseif (isol .eq. 2) then
r2p  = rp0
r2pp = rpp0
else
r3p  = rp0
r3pp = rpp0
endif


end do

    
end subroutine



subroutine scalar_riccati33(vars,ier,eps,a,b,c,r,rp,rpp,fun,par1,par2,par3,chebsol,rs,rders,rder2s)
implicit double precision (a-h,o-z)
type(scalar_vars_t)                         :: vars
procedure(scalar_fun33)                     :: fun
double complex                              :: par1,par2,par3
double complex                              :: r, rp, rpp
type(chebpw_scheme)                         :: chebsol
double complex, allocatable, intent(out)    :: rs(:), rders(:), rder2s(:)
!
!  Solve the Riccati equation corresponding to (1) given the 
!
!  Input parameters:
!
!  Output parameters:
!
!
type(c_ptr)                                :: userptr
type(scalar_odedata33), pointer            :: odedata
double complex, allocatable                :: rder3s(:)
! double complex                             :: yc(2), ima
! double complex, allocatable                :: ys0(:), yders0(:), yder2s0(:)

ier   = 0
eps0  = epsilon(0.0d0)
ima   = (0.0d0,1.0d0)

allocate(odedata)
userptr = c_loc(odedata)

n              = 2
odedata%fun   => fun
odedata%par1   = par1
odedata%par2   = par2
odedata%par3   = par3
odedata%wa     = 1
odedata%wb     = 1
odedata%ifleft = -1



call elapsed(t1)
call odetwo_nonlinear1(vars%odetwo_vars,jer,eps,a,b,c,rp,rpp,chebsol,rders,rder2s,rder3s, &
  scalar_ode33,userptr)
call elapsed(t2)


if (jer .ne. 0) then
ier = 4
return
endif

call chebpw_int(chebsol,rders,c,r,rs)

end subroutine


subroutine scalar_ode33(t,y,yp,f,dfdy,dfdyp,userptr)
implicit double precision (a-h,p-z)
type(c_ptr)                     :: userptr
double complex                  :: y, yp, f, dfdy, df, dfdyp
type(scalar_odedata33), pointer :: odedata

double complex                 :: q2, q1, q0
double complex                 :: u2, u1, u0


call c_f_pointer(userptr,odedata)

wa     = odedata%wa
wb     = odedata%wb

dk0    = odedata%dk0
dk1    = odedata%dk1
dk2    = odedata%dk2

ifleft = odedata%ifleft
phi    = (erf(wa*(t-wb))+1)/2
call odedata%fun(t,u0,u1,u2,odedata%par1,odedata%par2,odedata%par3)

if (ifleft .eq. 1) then
q0 = phi*dk0 + (1-phi)*u0
q1 = phi*dk1 + (1-phi)*u1
q2 = phi*dk2 + (1-phi)*u2
elseif (ifleft .eq. 0) then
q0 = (1-phi)*dk0 + phi*u0
q1 = (1-phi)*dk1 + phi*u1
q2 = (1-phi)*dk2 + phi*u2
else
q0 = u0
q1 = u1
q2 = u2
endif

f     = -q0-q1*y-q2*y**2-q2*yp-y**3-3*y*yp
dfdy  = -q1-2*y*q2-3*y**2-3*yp
dfdyp = -q2-3*y


end subroutine



subroutine root_estimates(n, dk, roots)
   implicit double precision (a-h, o-z)
   integer                                     :: n
   double precision                            :: dk
   double complex, allocatable                 :: roots(:)
   double complex                              :: ima

   ima = (0.0d0, 1.0d0)
   pi      = acos(-1.0d0)


   allocate(roots(n))
   
   do i = 0, n-1
       roots(i+1) = dk * exp(2 * pi * ima * i / n)
   end do
end subroutine



end module
