!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Experiment5:  an example with large real eigenvalues
!
!  lambda1: Ik
!  lambda2: -Ik
!  lambda3: k
!  lambda4: -k
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module experiment5_functions
use utils
use linalg0
use chebyshev
use chebpw
use odesolve
use scalar
use legendre
use odetwo

contains



subroutine plot1(dk)
implicit double precision (a-h,o-z)

type(scalar_vars_t)              :: vars
double complex                   :: ima, val
type(chebpw_scheme)              :: chebsol1, chebsol2, chebsol3, chebsol4
double complex, allocatable      :: r1(:,:), r2(:,:), r3(:,:), r4(:,:)
double complex, pointer          :: pars(:), rvals(:,:), rs(:,:,:)

double precision, allocatable    :: xs(:), ys1(:), ys2(:)


ima     = (0.0d0,1.0d0)
pi      = acos(-1.0d0)
eps0    = epsilon(0.0d0)
eps1    = 1.0d-14
eps2    = 1.0d-12
errmax  = 0.0d0

call scalar_init(vars)

n = 4
allocate(pars(1),rvals(4,4))
call elapsed(t1)

rvals    = 0

a        = -1.0d0
b        =  1.0d0

ifleft   =  1
pars(1)  = dk


a0       = -1.00d0
b0       = -0.90d0

call elapsed(t11)
call scalar_levin4(vars,ier,eps1,a0,b0,ifleft,testfun2,pars,rvals(2:4,:))


if (ier .ne. 0) then
call prini("levin ier = ",ier)
stop
endif

!call prinz("rvals = ",rvals)

if (ifleft .eq. 1) c = a0
if (ifleft .eq. 0) c = b0

call scalar_riccati4(vars,ier1,eps2,a,b,c,rvals(:,1),testfun2,pars,chebsol1,r1)
call scalar_riccati4(vars,ier2,eps2,a,b,c,rvals(:,2),testfun2,pars,chebsol2,r2)
call scalar_riccati4(vars,ier3,eps2,a,b,c,rvals(:,3),testfun2,pars,chebsol3,r3)
call scalar_riccati4(vars,ier4,eps2,a,b,c,rvals(:,4),testfun2,pars,chebsol4,r4)
call elapsed(t22)

print *, t22 - t11
stop
! call prinz("", r1)
! stop





!print *,chebsol1%nints, chebsol2%nints, chebsol3%nints,chebsol4%nints
if (ier1+ier2+ier3+ier4 .gt. 0)then
print *,ier1,ier2,ier3,ier4
stop
endif




stop


nn = 10000
allocate(xs(nn),ys1(nn),ys2(nn))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nn
t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
call chebpw_interp(chebsol1,r1(:,2),t,val)
xs(i)  = t
ys1(i) = real(val)
ys2(i) = imag(val)

if (abs(ys1(i)) .lt. 1.0d-16) ys1(i)=0
if (abs(ys2(i)) .lt. 1.0d-16) ys2(i)=0

end do

c =   0
d  =  -1
call plot_functions("experiment5-plota11.py","experiment5-plota11.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys1,"","","")


c =  0
d  = -1

call plot_functions("experiment5-plota21.py","experiment5-plota21.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys2,"","","")



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nn
t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
call chebpw_interp(chebsol2,r2(:,2),t,val)
xs(i)  = t
ys1(i) = real(val)
ys2(i) = imag(val)

if (abs(ys1(i)) .lt. 1.0d-16) ys1(i)=0
if (abs(ys2(i)) .lt. 1.0d-16) ys2(i)=0

end do

c =  0 
d  =  -1
call plot_functions("experiment5-plota12.py","experiment5-plota12.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys1,"","","")


c =   0
d  =  -1

call plot_functions("experiment5-plota22.py","experiment5-plota22.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys2,"","","")



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nn
t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
call chebpw_interp(chebsol3,r3(:,2),t,val)
xs(i)  = t
ys1(i) = real(val)
ys2(i) = imag(val)
if (abs(ys1(i)) .lt. 1.0d-16) ys1(i)=0
if (abs(ys2(i)) .lt. 1.0d-16) ys2(i)=0

end do

c =  0
d  =  -1
call plot_functions("experiment5-plota13.py","experiment5-plota13.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys1,"","","")


c =   0
d  =  -1

call plot_functions("experiment5-plota23.py","experiment5-plota23.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys2,"","","")



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nn
t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
call chebpw_interp(chebsol4,r4(:,2),t,val)
xs(i)  = t
ys1(i) = real(val)
ys2(i) = imag(val)

if (abs(ys1(i)) .lt. 1.0d-16) ys1(i)=0
if (abs(ys2(i)) .lt. 1.0d-16) ys2(i)=0

end do

c =   0
d  =  -1
call plot_functions("experiment5-plota14.py","experiment5-plota14.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys1,"","","")


c =   0
d  =  -1

call plot_functions("experiment5-plota24.py","experiment5-plota24.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys2,"","","")





end subroutine



subroutine plot2(dk)
implicit double precision (a-h,o-z)

type(scalar_vars_t)              :: vars
double complex                   :: ima, val
type(chebpw_scheme)              :: chebsol
double complex, allocatable      :: rs(:,:,:), rvals(:)
double complex, pointer          :: pars(:)


double precision, allocatable    :: xs(:), ys1(:), ys2(:)

ima     = (0.0d0,1.0d0)
pi      = acos(-1.0d0)
eps0    = epsilon(0.0d0)
eps     = 1.0d-12
errmax  = 0.0d0

call scalar_init(vars)

n = 4
allocate(pars(1),rvals(4))
call elapsed(t1)

rvals    = 0
a        = -1.0d0
b        =  1.0d0
c        = a
pars(1)  = dk


call scalar_levinode4(vars,ier,eps,a,b,c,rvals,testfun2,pars,chebsol,rs)
if (ier .ne. 0) then
call prini("levin2 ier = ",ier)
stop
endif

call elapsed(t2)
dtime = t2-t1


nn = 10000
allocate(xs(nn),ys1(nn),ys2(nn))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nn
t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
call chebpw_interp(chebsol,rs(:,2,1),t,val)
xs(i)  = t
ys1(i) = real(val)
ys2(i) = imag(val)

if (abs(ys1(i)) .lt. 1.0d-16) ys1(i)=0
if (abs(ys2(i)) .lt. 1.0d-16) ys2(i)=0

end do

c =   0
d  =  -1
call plot_functions("experiment5-plotb11.py","experiment5-plotb11.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys1,"","","")


c =  0
d  = -1

call plot_functions("experiment5-plotb21.py","experiment5-plotb21.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys2,"","","")



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nn
t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
call chebpw_interp(chebsol,rs(:,2,2),t,val)
xs(i)  = t
ys1(i) = real(val)
ys2(i) = imag(val)

if (abs(ys1(i)) .lt. 1.0d-16) ys1(i)=0
if (abs(ys2(i)) .lt. 1.0d-16) ys2(i)=0

end do

c =  0 
d  =  -1
call plot_functions("experiment5-plotb12.py","experiment5-plotb12.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys1,"","","")


c =   0
d  =  -1

call plot_functions("experiment5-plotb22.py","experiment5-plotb22.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys2,"","","")



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nn
t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
call chebpw_interp(chebsol,rs(:,2,3),t,val)
xs(i)  = t
ys1(i) = real(val)
ys2(i) = imag(val)
if (abs(ys1(i)) .lt. 1.0d-16) ys1(i)=0
if (abs(ys2(i)) .lt. 1.0d-16) ys2(i)=0

end do

c =  0
d  =  -1
call plot_functions("experiment5-plotb13.py","experiment5-plotb13.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys1,"","","")


c =   0
d  =  -1

call plot_functions("experiment5-plotb23.py","experiment5-plotb23.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys2,"","","")



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nn
t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
call chebpw_interp(chebsol,rs(:,2,4),t,val)
xs(i)  = t
ys1(i) = real(val)
ys2(i) = imag(val)

if (abs(ys1(i)) .lt. 1.0d-16) ys1(i)=0
if (abs(ys2(i)) .lt. 1.0d-16) ys2(i)=0

end do

c =   0
d  =  -1
call plot_functions("experiment5-plotb14.py","experiment5-plotb14.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys1,"","","")


c =   0
d  =  -1

call plot_functions("experiment5-plotb24.py","experiment5-plotb24.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys2,"","","")



end subroutine



subroutine plot3(dk)
implicit double precision (a-h,o-z)

type(scalar_vars_t)              :: vars
double complex                   :: ima, val
type(chebpw_scheme)              :: chebsol
double complex, allocatable      :: rs(:,:,:), rvals(:)
double complex, pointer          :: pars(:)


double precision, allocatable    :: xs(:), ys1(:), ys2(:)

ima     = (0.0d0,1.0d0)
pi      = acos(-1.0d0)
eps0    = epsilon(0.0d0)
eps     = 1.0d-12
errmax  = 0.0d0

call scalar_init(vars)

n = 4
allocate(pars(1),rvals(4))
call elapsed(t1)

rvals    = 0
a        = -1.0d0
b        =  1.0d0
c        = a
pars(1)  = dk


call scalar_levinode4(vars,ier,eps,a,b,c,rvals,testfun2,pars,chebsol,rs)
if (ier .ne. 0) then
call prini("levin2 ier = ",ier)
stop
endif

call elapsed(t2)
dtime = t2-t1


nn = 10000
allocate(xs(nn),ys1(nn),ys2(nn))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nn
t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
call chebpw_interp(chebsol,rs(:,2,1),t,val)
xs(i)  = t
ys1(i) = real(val)
ys2(i) = imag(val)

if (abs(ys1(i)) .lt. 1.0d-16) ys1(i)=0
if (abs(ys2(i)) .lt. 1.0d-16) ys2(i)=0

end do

c =   0
d  =  -1
call plot_functions("experiment5-plotc11.py","experiment5-plotc11.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys1,"","","")


c =  0
d  = -1

call plot_functions("experiment5-plotc21.py","experiment5-plotc21.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys2,"","","")



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nn
t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
call chebpw_interp(chebsol,rs(:,2,2),t,val)
xs(i)  = t
ys1(i) = real(val)
ys2(i) = imag(val)

if (abs(ys1(i)) .lt. 1.0d-16) ys1(i)=0
if (abs(ys2(i)) .lt. 1.0d-16) ys2(i)=0

end do

c =  0 
d  =  -1
call plot_functions("experiment5-plotc12.py","experiment5-plotc12.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys1,"","","")


c =   0
d  =  -1

call plot_functions("experiment5-plotc22.py","experiment5-plotc22.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys2,"","","")



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nn
t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
call chebpw_interp(chebsol,rs(:,2,3),t,val)
xs(i)  = t
ys1(i) = real(val)
ys2(i) = imag(val)
if (abs(ys1(i)) .lt. 1.0d-16) ys1(i)=0
if (abs(ys2(i)) .lt. 1.0d-16) ys2(i)=0

end do

c =  0
d  =  -1
call plot_functions("experiment5-plotc13.py","experiment5-plotc13.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys1,"","","")


c =   0
d  =  -1

call plot_functions("experiment5-plotc23.py","experiment5-plotc23.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys2,"","","")



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nn
t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
call chebpw_interp(chebsol,rs(:,2,4),t,val)
xs(i)  = t
ys1(i) = real(val)
ys2(i) = imag(val)

if (abs(ys1(i)) .lt. 1.0d-16) ys1(i)=0
if (abs(ys2(i)) .lt. 1.0d-16) ys2(i)=0

end do

c =   0
d  =  -1
call plot_functions("experiment5-plotc14.py","experiment5-plotc14.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys1,"","","")


c =   0
d  =  -1

call plot_functions("experiment5-plotc24.py","experiment5-plotc24.pdf",1,"","",0,0,a,b,c,d, &
   nn,xs,ys2,"","","")



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
f(3)    = y(4)
f(4)    = -qs(1)*y(1)-qs(2)*y(2)-qs(3)*y(3)-qs(4)*y(4)

df(1,1) = 0
df(1,2) = 1
df(1,3) = 0
df(1,4) = 0

df(2,1) = 0
df(2,2) = 0
df(2,3) = 1
df(2,4) = 0

df(3,1) = 0
df(3,2) = 0
df(3,3) = 0
df(3,4) = 1

df(4,1) = -qs(1)
df(4,2) = -qs(2)
df(4,3) = -qs(3)
df(5,4) = -qs(4)

end subroutine


subroutine testfun2(n,t,qs,pars)
implicit double precision (a-h,o-z)
double complex                  :: qs(:)
double complex, pointer         :: pars(:)
double complex                  :: dk0, dk1, dk2, ima
data ima / (0.0d0,1.0d0) /

dk     = pars(1)

qs(1)  = - dk**4  * (1+t**2)
qs(2)  = 0
qs(3)  = 0
qs(4)  = 0

end subroutine

end module



program experiment5

use utils
use linalg0
use chebyshev
use chebpw
use odesolve
use scalar
use experiment5_functions

implicit double precision (a-h,o-z)

double precision, allocatable :: ts(:)
double complex, allocatable   :: vals0(:)
double complex, allocatable   :: vals1(:)
double complex, allocatable   :: vals2(:)
double complex, allocatable   :: vals3(:)

double precision, allocatable :: dks(:), residues(:,:), times(:,:)

eps0  = epsilon(0.0d0)
ntest = 10000



dk = 2000
call plot1(dk)
stop

dk = 2**18
call plot2(dk)

dk = 2**20
call plot3(dk)

stop

end program
