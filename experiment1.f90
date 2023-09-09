!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Experiment1:  Solve an initial value problem for a 2x2 system
!
!  lambda1 = I*k
!  lambda2 = -I*k
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module experiment1_functions
use utils
use linalg0
use chebyshev
use chebpw
use odesolve
use scalar
use legendre
use odetwo

contains

subroutine test0(eps,dk,ntest,vals,dtime)
implicit double precision (a-h,o-z)

type(odesolve_vars_t)            :: vars
type(odetwo_vars_t)              :: odetwo_vars

type(chebpw_scheme)              :: chebsol
double complex, pointer          :: pars(:)
double complex, allocatable      :: ys(:,:),yders(:,:)
double complex                   :: ima, yc(2), val


double complex                   :: vals(:)

double complex, allocatable      :: y0(:),yp0(:),ypp0(:)
double complex                   :: yc0,ypc0

call odesolve_init(vars,30,4,1000000,1000000,12,12)
call odetwo_init(odetwo_vars,30,4,1000000,1000000,12,12)

call elapsed(t1)
n     = 2
a     = -1.0d0
b     =  1.0d0
c     = a
!eps   = 1.0d-12
ima   = (0.0d0,1.0d0)
pi    = acos(-1.0d0)

yc(1) = 1.0d0
yc(2) = 1.0d0
ifit  = 2

allocate(pars(1))
pars(1) = dk

yc0   = 1
ypc0  = 1
call  odetwo_nonlinear2(odetwo_vars,ier,eps,a,b,c,yc0,ypc0,chebsol,y0,yp0,ypp0,testfun0,pars)

if (ier .ne. 0) then
print *,"odetwo failed, ier = ",ier
stop
endif

nn = chebsol%nints*chebsol%k
allocate(ys(nn,3))
ys(:,1) = y0
ys(:,2) = yp0
ys(:,3) = ypp0

!call odesolve_nonlinear2(vars,ier,ifit,eps,a,b,c,yc,n,testfun1,pars,chebsol,ys,yders)
if (ier .ne. 0) then
print *,"odesolve failed"
stop
endif

do i=1,ntest
t = a + (b-a)*(i-1.0d0)/(ntest-1.0d0)
call chebpw_interp(chebsol,ys(:,1),t,val)
vals(i)  = val
end do
call elapsed(t2)

dtime=t2-t1

end subroutine


subroutine test1(dk,ntest,vals,dtime,ncoefs)
implicit double precision (a-h,o-z)

type(scalar_vars_t)              :: vars
double complex                   :: ima
type(chebpw_scheme)              :: chebsol1, chebsol2
double complex, allocatable      :: r1(:,:), r2(:,:)
double complex, pointer          :: pars(:), rvals(:,:), rs(:,:,:)
double complex                   :: a11,a12,a22,a21,c1,c2,det
double complex                   :: val1,val2,der1,der2,val0,der0
double complex                   :: vals(:)

ima     = (0.0d0,1.0d0)
pi      = acos(-1.0d0)
eps0    = epsilon(0.0d0)
eps     = 1.0d-12
errmax  = 0.0d0

call scalar_init(vars)

n = 2
allocate(pars(1),rvals(2,2))
call elapsed(t1)

rvals    = 0

a        = -1.0d0
b        =  1.0d0

ifleft   =  1
pars(1)  = dk

a0       = -1.00d0
b0       =  -0.90d0
call scalar_levin2(vars,ier,eps,a0,b0,ifleft,testfun2,pars,rvals(2:2,:))

if (ier .ne. 0) then
call prini("levin2 ier = ",ier)
stop
endif


call prinz("", rvals(2:2, :))

stop


if (ifleft .eq. 1) c = a0
if (ifleft .eq. 0) c = b0

call scalar_riccati2(vars,ier1,eps,a,b,c,rvals(:,1),testfun2,pars,chebsol1,r1)
call scalar_riccati2(vars,ier2,eps,a,b,c,rvals(:,2),testfun2,pars,chebsol2,r2)

call prinz("", r1)
call prinz("", r2)

ncoefs = chebsol1%nints*chebsol1%k + chebsol2%nints*chebsol2%k
!
!  Find the coefficients
!

t    = a
val0 = 1
der0 = 1

call chebpw_interp(chebsol1,r1(:,1),t,val1)
call chebpw_interp(chebsol1,r1(:,2),t,der1)

call chebpw_interp(chebsol2,r2(:,1),t,val2)
call chebpw_interp(chebsol2,r2(:,2),t,der2)

a11 = exp(val1)
a12 = exp(val2)
a21 = exp(val1)*der1
a22 = exp(val2)*der2

det = a11*a22-a12*a21
c1  = ( a22*val0-a12*der0)/det
c2  = (-a21*val0+a11*der0)/det

call elapsed(t2)
dtime=t2-t1

do i=1,ntest
t    = a + (b-a)*(i-1.0d0)/(ntest-1.0d0)

call chebpw_interp(chebsol1,r1(:,1),t,val1)
call chebpw_interp(chebsol2,r2(:,1),t,val2)

vals(i) = c1*exp(val1)+c2*exp(val2)
end do

end subroutine




subroutine test2(dk,ntest,vals,dtime, ncoefs)
implicit double precision (a-h,o-z)

type(scalar_vars_t)              :: vars
double complex                   :: ima
type(chebpw_scheme)              :: chebsol
double complex, allocatable      :: rs(:,:,:), rvals(:)
double complex, pointer          :: pars(:)
double complex                   :: a11,a12,a22,a21,c1,c2,det
double complex                   :: val1,val2,der1,der2,val0,der0
double complex                   :: vals(:)

ima     = (0.0d0,1.0d0)
pi      = acos(-1.0d0)
eps0    = epsilon(0.0d0)
eps     = 1.0d-12
errmax  = 0.0d0

call scalar_init(vars)

n = 2
allocate(pars(1),rvals(2))
call elapsed(t1)

rvals    = 0
a        = -1.0d0
b        =  1.0d0
c        = a
ifleft   = 1
pars(1)  = dk

call scalar_levinode2(vars,ier,eps,a,b,c,rvals,testfun2,pars,chebsol,rs)
if (ier .ne. 0) then
call prini("levin2 ier = ",ier)
stop
endif

ncoefs = 2 * chebsol%nints * chebsol%k

!
!  Find the coefficients
!

t    = a
val0 = 1
der0 = 1

call chebpw_interp(chebsol,rs(:,1,1),t,val1)
call chebpw_interp(chebsol,rs(:,2,1),t,der1)

call chebpw_interp(chebsol,rs(:,1,2),t,val2)
call chebpw_interp(chebsol,rs(:,2,2),t,der2)

a11 = exp(val1)
a12 = exp(val2)
a21 = exp(val1)*der1
a22 = exp(val2)*der2

det = a11*a22-a12*a21
c1  = ( a22*val0-a12*der0)/det
c2  = (-a21*val0+a11*der0)/det

call elapsed(t2)
dtime=t2-t1


do i=1,ntest
t    = a + (b-a)*(i-1.0d0)/(ntest-1.0d0)

call chebpw_interp(chebsol,rs(:,1,1),t,val1)
call chebpw_interp(chebsol,rs(:,1,2),t,val2)

vals(i) = c1*exp(val1)+c2*exp(val2)
end do

end subroutine


subroutine plot2(dk)
    implicit double precision (a-h,o-z)
    
    type(scalar_vars_t)              :: vars
    double complex                   :: ima, val
    type(chebpw_scheme)              :: chebsol
    double complex, allocatable      :: rs(:,:,:), rvals(:)
    double complex, pointer          :: pars(:)
    
    
    double precision, allocatable     :: xs(:), ys1(:), ys2(:)
    
    ima     = (0.0d0,1.0d0)
    pi      = acos(-1.0d0)
    eps0    = epsilon(0.0d0)
    eps     = 1.0d-12
    errmax  = 0.0d0
    
    call scalar_init(vars)
    
    n = 2
    allocate(pars(1),rvals(2))
    call elapsed(t1)
    
    rvals    = 0
    a        = -1.0d0
    b        =  1.0d0
    c        = a
    pars(1)  = dk
    
    call scalar_levinode2(vars,ier,eps,a,b,c,rvals,testfun2,pars,chebsol,rs)
    if (ier .ne. 0) then
    call prini("levin2 ier = ",ier)
    stop
    endif

   
    
    
    nn = 10000
    allocate(xs(nn),ys1(nn),ys2(nn))
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,nn
    t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
    call chebpw_interp(chebsol,rs(:,2,1),t,val)
    xs(i)  = t
    ys1(i) = real(val)
    ys2(i) = imag(val)
    end do
    
    c =  0 
    d  = -1

    call plot_functions("experiment1-plot11.py","experiment1-plot11.pdf",1,"t","",0,0,a,b,0.0d0,0.3d0, &
       nn,xs,ys1,"","","")
    
    
    call plot_functions("experiment1-plot21.py","experiment1-plot21.pdf",1,"t","",0,0,a,b,-1.0d6,-1.7d5, &
       nn,xs,ys2,"","","")
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,nn
    t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
    call chebpw_interp(chebsol,rs(:,2,2),t,val)
    xs(i)  = t
    ys1(i) = real(val)
    ys2(i) = imag(val)
    end do
    
    call plot_functions("experiment1-plot12.py","experiment1-plot12.pdf",1,"t","",0,0,a,b,c,d, &
       nn,xs,ys1,"","","")
    
    
    call plot_functions("experiment1-plot22.py","experiment1-plot22.pdf",1,"t","",0,0,a,b,c,d, &
       nn,xs,ys2,"","","")
    
end subroutine



subroutine testfun0(t,y,yp,f,dfdy,dfdyp,pars)
use iso_c_binding
implicit double precision (a-h,o-z)
double complex              :: y, yp, f, dfdy, df, dfdyp
double complex, pointer     :: pars(:)
double complex              :: qs(2), ima
data ima /(0.0d0,1.0d0)/

dk     = pars(1)
call testfun2(n,t,qs,pars)


f      = -qs(2)*yp -qs(1)*y
dfdy   = -qs(1)
dfdyp  = -qs(2)

end subroutine


subroutine testfun1(n,t,y,f,df,pars)
implicit double precision (a-h,p-z)
integer                    :: n
double complex, pointer    :: pars(:)
double complex             :: y(:), f(:), df(:,:)
double complex             :: qs(2), ima
data ima / (0.0d0,1.0d0) /

dk     = pars(1)
call testfun2(n,t,qs,pars)


f(1)    = y(2)
f(2)    = -qs(1)*y(1)-qs(2)*y(2)

df(1,1) = 0
df(1,2) = 1

df(2,1) = -qs(1)
df(2,2) = -qs(2)

end subroutine


subroutine testfun2(n,t,qs,pars)
implicit double precision (a-h,o-z)
double complex                  :: qs(:)
double complex, pointer         :: pars(:)
double complex                  :: dk0, dk1, dk2, ima
data ima / (0.0d0,1.0d0) /

dk     = pars(1)

!qs(1)  = dk**2/(1+t**2)*(1+cos(t)**2)
! (1.5d0 + cos(log(dk)*t))

qs(1)  = dk ** n * (1 + t ** 2)
qs(2)  = 0!-ima*dk*(1+t**2)

end subroutine


end module


program experiment1

use utils
use linalg0
use chebyshev
use chebpw
use odesolve
use scalar
use experiment1_functions

implicit double precision (a-h,o-z)

double precision, allocatable :: ts(:)
double complex, allocatable   :: vals0(:)
double complex, allocatable   :: vals1(:)
double complex, allocatable   :: vals2(:)
double complex                :: ima

double precision, allocatable :: times(:, :), errs(:, :), dks(:), dcoefs(:, :)

double complex              :: test(2), w2(2)
integer                       :: index(4)

ima = (0.0d0, 1.0d0)
eps0  = epsilon(0.0d0)
ntest = 10000
eps   = 1.0d-14
ii1 = 8
ii2 = 14

nn = ii2 - ii1 + 1

allocate(ts(ntest),vals0(ntest),vals1(ntest),vals2(ntest))

allocate(times(2, nn), errs(2, nn), dks(nn), dcoefs(2, nn))

dk = 200

call test1(dk,ntest,vals1,dtime1,ncoefs1) 

stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Just do it in double precision
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ii=ii1,ii2
   mm = ii - ii1 + 1
   dk  = 2.0d0**ii
   dtime0 = -1
   dtime1 = -1
   dtime2 = -1
   dtime3 = -1
   call test0(eps,dk,ntest,vals0,dtime0)
   call test1(dk,ntest,vals1,dtime1,ncoefs1)                 ! local Levin
   call test2(dk,ntest,vals2,dtime2,ncoefs2)                 ! global Levin
   errmax1 = maxval(abs(vals0-vals1))
   errmax2 = maxval(abs(vals0-vals2))
   write (*,"(4(D15.7,1X),4(D15.4,1X))") dk,errmax1,errmax2,dtime0,dtime1,dtime2
   
   dks(mm)       = dk
   times(2, mm)  = dtime1 * 1000.0d0
   times(1, mm)  = dtime2 * 1000.0d0
   errs(2, mm)   = errmax1
   errs(1, mm)   = errmax2
   dcoefs(2, mm) = ncoefs1
   dcoefs(1, mm) = ncoefs2

end do
a = 2 ** ii1
b = 2 ** ii2

tmin = 0.0d0
tmax = 10.0d0

errmin = 1.0d-16
errmax = 1.0d-6


call plot_functions("experiment1-times.py","experiment1-times.pdf",2,"$k$","Time (millisecond)",2,0, &
        a, b, tmin, tmax, nn, dks, times(:, :),"best","global Levin method*local Levin method*","b*k--*g:*")

call plot_functions("experiment1-errs.py","experiment1-errs.pdf",2,"$k$","Maximum Absolute Error",2,1, &
        a, b, errmin, errmax, nn, dks, errs(:, :),"best","global Levin method*local Levin method*","b*k--*g:*")


call plot_functions("experiment1-coefs.py","experiment1-coefs.pdf",2,"$k$","Chebyshev Coefficients",2,0, &
        a,b,0.0d0,200.0d0,nn,dks,dcoefs,"best","global Levin method*local Levin method*","b*k--*g:*")


dk = 2.0d0 ** 20

call plot2(dk)
stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Comparison with extended precision version of the experiment
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! iw = 101

! if (eps0 .lt. 1.0d-16) then
!    eps = 1.0d-16
    
!    open(iw,FILE="experiment1.save")
    
!    do ii=ii1,ii2
!    dk = 2.0d0**ii
!    call test0(eps,dk,ntest,vals0,dtime)
!    do i=1,ntest
!    write(iw,"(3(D44.26,1X))") ts(i),vals0(i)
!    end do

!    print *,dk,dtime
!    end do

!    close(iw)

! else
!    eps = 1.0d-12
    
!    open(iw,FILE="experiment1.save",STATUS='old')

!    do ii=ii1,ii2
!    dk = 2.0d0**ii
    
!    do i=1,ntest
!    read(iw,"(3(D44.26,1X))") ts(i),vals0(i)
!    end do
    
!    call test1(dk,ntest,vals1,dtime1)
!    call test2(dk,ntest,vals2,dtime2)
    
!    errmax1 = maxval(abs(vals0-vals1))
!    errmax2 = maxval(abs(vals0-vals2))
    
!    write (*,"(4(D15.7,1X),3(D15.4,1X))") dk,errmax1,errmax2,dtime0,dtime1,dtime2
   
!    end do

!    close(iw)


! endif


end program
