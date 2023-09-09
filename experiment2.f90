!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Experiment2:  Solve an initial value problem for a 3x3 system
!
!
!  lambda1 = -I*k
!  lambda2 = I*k
!  lambda3 = -1
!
!  **PERTURBED AWAY ENOUGH THAT ALL EIGENVALUES ARE LARGE ***
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module experiment2_functions
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

double complex, pointer          :: pars(:)

type(chebpw_scheme)              :: chebsol
double complex, allocatable      :: ys(:,:),yders(:,:)



double complex                   :: ima, yc(3), val
double complex                   :: vals(:)

double complex, allocatable      :: y0(:),yp0(:),ypp0(:)
double complex                   :: yc0,ypc0

call odesolve_init(vars,30,4,1000000,1000000,12,12)
call odetwo_init(odetwo_vars,30,4,1000000,1000000,12,12)

call elapsed(t1)
n     = 3
a     = -1.0d0
b     =  1.0d0
ima   = (0.0d0,1.0d0)
pi    = acos(-1.0d0)
eps0  = epsilon(0.0d0)

c     = 0
yc(1) = 1.0d0
yc(2) = -ima*dk
yc(3) = -dk**2

ifit  = 2

allocate(pars(1))
pars(1) = dk


call odesolve_nonlinear2(vars,ier,ifit,eps,a,b,c,yc,n,testfun1,pars,chebsol,ys,yders)
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

! call chebpw_plot(chebsol,ys(:,1),a,b,"y.py","y.pdf")
! stop

end subroutine


subroutine test1(dk,ntest,vals,dtime)
implicit double precision (a-h,o-z)

type(scalar_vars_t)              :: vars
double complex                   :: ima
type(chebpw_scheme)              :: chebsol1, chebsol2, chebsol3
double complex, allocatable      :: r1(:,:), r2(:,:), r3(:,:)
double complex, pointer          :: pars(:), rvals(:,:), rs(:,:,:)

double complex                   :: amatr(3,3), cs(3), yc(3)
double complex                   :: val, der, der2
double complex                   :: val1, val2, val3, vals(:)

ima     = (0.0d0,1.0d0)
pi      = acos(-1.0d0)
eps0    = epsilon(0.0d0)
eps1    = 1.0d-13
eps2    = 1.0d-13
errmax  = 0.0d0

call scalar_init(vars)

n = 3
allocate(pars(1),rvals(3,3))
call elapsed(t1)

rvals    = 0

a        = -1.0d0
b        =  1.0d0

ifleft   =  1
pars(1)  = dk

a0       = -1.00d0
b0       = -0.90d0

call elapsed(t11)

call scalar_levin3(vars,ier,eps1,a0,b0,ifleft,testfun2,pars,rvals(2:3,:))


if (ier .ne. 0) then
call prini("levin ier = ",ier)
stop
endif

if (ifleft .eq. 1) c = a0
if (ifleft .eq. 0) c = b0


call scalar_riccati3(vars,ier1,eps2,a,b,c,rvals(:,1),testfun2,pars,chebsol1,r1)
call scalar_riccati3(vars,ier2,eps2,a,b,c,rvals(:,2),testfun2,pars,chebsol2,r2)
call scalar_riccati3(vars,ier2,eps2,a,b,c,rvals(:,3),testfun2,pars,chebsol3,r3)


call elapsed(t22)

dtime = t22 - t11

print *, "time = ", dtime

!call prinz("", r3)
stop


!
!  Find the coefficients
!

t     = 0
yc(1) = 1.0d0
yc(2) = -ima*dk
yc(3) = -dk**2

call chebpw_interp(chebsol1,r1(:,1),t,val)
call chebpw_interp(chebsol1,r1(:,2),t,der)
call chebpw_interp(chebsol1,r1(:,3),t,der2)

amatr(1,1) = exp(val)
amatr(2,1) = exp(val)*der
amatr(3,1) = exp(val)*(der**2 + der2)

call chebpw_interp(chebsol2,r2(:,1),t,val)
call chebpw_interp(chebsol2,r2(:,2),t,der)
call chebpw_interp(chebsol2,r2(:,3),t,der2)

amatr(1,2) = exp(val)
amatr(2,2) = exp(val)*der
amatr(3,2) = exp(val)*der**2 + exp(val)*der2

call chebpw_interp(chebsol3,r3(:,1),t,val)
call chebpw_interp(chebsol3,r3(:,2),t,der)
call chebpw_interp(chebsol3,r3(:,3),t,der2)


amatr(1,3) = exp(val)
amatr(2,3) = exp(val)*der
amatr(3,3) = exp(val)*der**2 + exp(val)*der2

call linalg0_qrsolve(eps0*10,n,n,amatr,cs,yc)

call elapsed(t2)
dtime=t2-t1

do i=1,ntest
t    = a + (b-a)*(i-1.0d0)/(ntest-1.0d0)

call chebpw_interp(chebsol1,r1(:,1),t,val1)
call chebpw_interp(chebsol2,r2(:,1),t,val2)
call chebpw_interp(chebsol3,r3(:,1),t,val3)

val = cs(1)*exp(val1)+cs(2)*exp(val2)+cs(3)*exp(val3)
vals(i) = val
end do

end subroutine



subroutine test2(dk,ntest,vals,dtime)
implicit double precision (a-h,o-z)

type(scalar_vars_t)              :: vars
double complex                   :: ima
type(chebpw_scheme)              :: chebsol
double complex, allocatable      :: rs(:,:,:), rvals(:)
double complex, pointer          :: pars(:)

double complex                   :: amatr(3,3), cs(3), yc(3)
double complex                   :: val, der, der2
double complex                   :: val1, val2, val3, vals(:)


ima     = (0.0d0,1.0d0)
pi      = acos(-1.0d0)
eps0    = epsilon(0.0d0)
eps     = 1.0d-12
errmax  = 0.0d0

call scalar_init(vars)

n = 3
allocate(pars(1),rvals(3))
call elapsed(t1)

rvals    = 0
a        = -1.0d0
b        =  1.0d0
c        = a
pars(1)  = dk

call scalar_levinode3(vars,ier,eps,a,b,c,rvals,testfun2,pars,chebsol,rs)
if (ier .ne. 0) then
call prini("levin2 ier = ",ier)
stop
endif

! call chebpw_plot(chebsol,rs(:,2,1),a,b,"r1.py","r1.pdf")
! call chebpw_plot(chebsol,rs(:,2,2),a,b,"r2.py","r2.pdf")
! call chebpw_plot(chebsol,rs(:,2,3),a,b,"r3.py","r3.pdf")

!
!  Find the coefficients
!

t     = 0
yc(1) = 1.0d0
yc(2) = -ima*dk
yc(3) = (-ima*dk)**2

call chebpw_interp(chebsol,rs(:,1,1),t,val)
call chebpw_interp(chebsol,rs(:,2,1),t,der)
call chebpw_interp(chebsol,rs(:,3,1),t,der2)

amatr(1,1) = exp(val)
amatr(2,1) = exp(val)*der
amatr(3,1) = exp(val)*(der**2 + der2)

call chebpw_interp(chebsol,rs(:,1,2),t,val)
call chebpw_interp(chebsol,rs(:,2,2),t,der)
call chebpw_interp(chebsol,rs(:,3,2),t,der2)

amatr(1,2) = exp(val)
amatr(2,2) = exp(val)*der
amatr(3,2) = exp(val)*der**2 + exp(val)*der2

call chebpw_interp(chebsol,rs(:,1,3),t,val)
call chebpw_interp(chebsol,rs(:,2,3),t,der)
call chebpw_interp(chebsol,rs(:,3,3),t,der2)


amatr(1,3) = exp(val)
amatr(2,3) = exp(val)*der
amatr(3,3) = exp(val)*der**2 + exp(val)*der2

cs = yc
!call linalg0_solve_c(n,amatr,cs)
call linalg0_qrsolve(eps0*10,n,n,amatr,cs,yc)


call elapsed(t2)
dtime=t2-t1


do i=1,ntest

t    = a + (b-a)*(i-1.0d0)/(ntest-1.0d0)

call chebpw_interp(chebsol,rs(:,1,1),t,val1)
call chebpw_interp(chebsol,rs(:,1,2),t,val2)
call chebpw_interp(chebsol,rs(:,1,3),t,val3)

val     = cs(1)*exp(val1)+cs(2)*exp(val2)+cs(3)*exp(val3)
vals(i) = val
end do

end subroutine




subroutine testfun1(n,t,y,f,df,pars)
implicit double precision (a-h,p-z)
integer                    :: n
double complex, pointer    :: pars(:)
double complex             :: y(:), f(:), df(:,:)
double complex             :: qs(2), ima
data ima / (0.0d0,1.0d0) /

call testfun2(n,t,qs,pars)

f(1)    = y(2)
f(2)    = y(3)
f(3)    = -qs(1)*y(1)-qs(2)*y(2)-qs(3)*y(3)

df(1,1) = 0
df(1,2) = 1
df(1,3) = 0

df(2,1) = 0
df(2,2) = 0
df(2,3) = 1

df(3,1) = -qs(1)
df(3,2) = -qs(2)
df(3,3) = -qs(3)

end subroutine

subroutine testfun2(n,t,qs,pars)
implicit double precision (a-h,o-z)
double complex                  :: qs(:)
double complex, pointer         :: pars(:)
double complex                  :: dk0, dk1, dk2, ima
data ima / (0.0d0,1.0d0) /

dk     = pars(1)

! qs(1)  = -ima*dk**2 * (2+exp(ima*t))/(1+t**2)
! qs(2)  = dk**2
! qs(3)  = -ima

qs(1)  = -dk ** 3 * (1 + t ** 2)
qs(2)  = 0!dk**2  * (1+sin(3*t)**2)
qs(3)  = 0!-1  * (cos(6*t)**2+2)
end subroutine

end module


program experiment2

use utils
use linalg0
use chebyshev
use chebpw
use odesolve
use scalar
use experiment2_functions

implicit double precision (a-h,o-z)

double precision, allocatable :: ts(:)
double complex, allocatable   :: vals0(:)
double complex, allocatable   :: vals1(:)
double complex, allocatable   :: vals2(:)
double complex, allocatable   :: vals3(:)


eps0  = epsilon(0.0d0)
ntest = 10000
eps   = 1.0d-14
ii1   = 8
ii2   = 20

allocate(ts(ntest),vals0(ntest),vals1(ntest),vals2(ntest),vals3(ntest))

dk = 2000
call test1(dk,ntest,vals1,dtime1)

stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Just do it in double precision
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ii=ii1,ii2
dk  = 2.0d0**ii

dtime0  = -1
dtime1  = -1
dtime2  = -1
dtime3  = -1
errmax1 = -1
errmax2 = -1
errmax3 = -1
vals0   =  0
vals1   =  0
vals2   =  0

call test0(eps,dk,ntest,vals0,dtime0)             ! brute force solver
call test1(dk,ntest,vals1,dtime1)                 ! local Levin
call test2(dk,ntest,vals2,dtime2)                 ! global Levin

errmax1 = maxval(abs(vals0-vals1))
errmax2 = maxval(abs(vals0-vals2))

write (*,"(4(D15.7,1X),4(D15.4,1X))") dk,errmax1,errmax2,dtime0,dtime1,dtime2
end do

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
