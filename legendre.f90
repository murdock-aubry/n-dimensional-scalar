!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for constructing Guass-Legendre quadrature rules, and for
!  forming and manipulating univariate Legendre expansions.  A univariate Legendre 
!  expansion is a sum of the form 
!
!              n-1
!      f(x) = \sum   a_i P_i(x),                                                            (1)
!              i=0
!
!  where P_i denotes the Legendre polynomial of degree i.  This code provides 
!  routines for both real-valued and complex-valued expansions.
!
!  Expansions of the form (1) are represented either via the vector of coefficients
!   
!    ( a_0     )
!    ( a_1     )
!    ( ...     )                                                                            (2)
!    ( a_{n-1} )
!
!  or via the vector 
!
!    ( f(x_1)  )
!    ( f(x_2)  )
!    (  ...    )                                                                           (3)
!    ( f(x_n)  )
!
!  of their values at the nodes of the n-point Gauss-Legendre quadrature rule.  
!
!  The following routines should be regarded as public:
!      
!    legendre_quad - return the nodes and weights of the n-point Gauss-Legendre 
!      quadrature rule; note that this routine uses an O(n^2) algorithm and is
!      intended for constructing quadrature rules of relatively small sizes (n < 200)
!
!    legendre_interp - use barycentric interpolation to evaluate an expansion 
!     of the form (1) given its values at nodes of the n-point Gauss-Legendre
!     quadrature rule
!
!    legendre_eval - evaluate one or more expansions of the form (1) at a specified
!      point given their coefficient vectors (2)
!
!    legendre_evalder - evaluate one or more expansion of the form (1) and their
!      derivatives at a specified point given the vector (2) of coefficients
!
!    legendre_coefsmatrix - return the (n,n) matrix which takes the vector (3) of 
!      the values of an expansion of the form (1) to the vector (2) of its
!      coefficients
!
!    legendre_interpmatrix - return an (m,n) matrix which takes the vector (3) of the
!      values of an expansion at the nodes of the n-point Gauss-Legendre rule
!      to its values at a user-specified collection of points
!
!    legendre_diffmatrix - return the (n,n) ``spectral differentiation matrix'' which
!      takes the vector (3) of values of an expansion of the form (1) to
!      the vector of values of the derivative of the expansion
!
!      **WARNING** SPECTRAL DIFFERENTIATION BECOME INCREASINGLY ILL-CONDITIONED AS THE
!                  ORDER OF THE EXPANSIONS GROW
!
!    legendre_intlmatrix - return the (n,n) spectral integration matrix which
!      takes the vector (3) of values of an expansion of the form (1) to
!      the vector of  values of the function
!
!                        x
!            g(x)  = \int f(x) dt
!                       -1
!
!      at the Legendre nodes
!
!    legendre_intrmatrix - return the (n,n) spectral integration matrix which
!      takes the vector (3) of values of an expansion of the form (1) to
!      the vector of values of the function
!
!                        x
!            g(x)  = \int f(x) dt
!                        1
!
!      at the Legendre nodes 
!
!  The following are low-level routines should be used with caution by the user:
!
!    leges - evaluate the first n Legendre polynomials at a specified point 
!      via the three-term recurrence relation
!
!    legesder - evaluate the first n Legendre polynomials and their 
!      derivatives at a specified point via the three-term recurrence relation
!
!    legesq - evaluate the first n Legendre functions of the second kind
!
!    legendre_pdnu - evaluate the Legendre function of the first kind of a noninteger
!      degree at a specified point using a combination of Taylor expansions and the 
!      three-term recurrence relation satsified by the Legendre functions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module legendre

interface     legendre_interp
module procedure legendre_interp1
module procedure legendre_interp2
module procedure legendre_interp3
module procedure legendre_interp4
end interface legendre_interp

interface     legendre_eval
module procedure legendre_eval1
module procedure legendre_eval2
module procedure legendre_eval3
module procedure legendre_eval4
end interface legendre_eval

interface     legendre_evalder
module procedure legendre_evalder1
module procedure legendre_evalder2
module procedure legendre_evalder3
module procedure legendre_evalder4
end interface legendre_evalder

contains


subroutine legendre_quad(n,xslege,whtslege)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)            :: xslege(:),whtslege(:)
!
!  Form the n-point Gauss-Legendre quadrature rule.  Note that this routine uses a
!  O(n^2) algorithm and is inteded for use in the case in which n is relatively
!  small (n < 100 or so).
!
!  Input parameters:
!    n - the length of the desired quadrature rule
!
!  Output parameters:
!    xslege - the nodes of the desired quadrature rule
!    whtlege - the weights of the desired quadrature rule
!

double precision :: pols(n+1)

data pi / 3.14159265358979323846264338327950288d0 /
allocate( xslege(n), whtslege(n) )

maxiters = 12
eps0     = epsilon(0.0d0)
xslege=0

!
!   Use Newton's method and the recurrence relation to find half of the
!   roots --- the others are obtained via symmetry.
!
!   Note that we also store the value of the derivative at each of the obtained
!   roots for use in computing the weights below.
!
ifodd = 0
nn = (n+1)/2
if (nn /= n/2) then
ifodd=1
nn=nn-1
endif

!
!  Use Tricomi's formula to estimate the roots of P_n
!

do i =nn+1,n
   dk = i-nn-1
   dn = n
   theta = (4*(ceiling(dn/2)-dk)-1) / (4*dn+2) * pi
   x0 = 1.0d0 - (dn-1)/(8*dn**3)-1.0d0/(384.0d0*dn**4)* &
        (39.0d0-28.0d0/sin(theta)**2)
   xslege(i)=x0*cos(theta)
enddo

!
!  Use Chebyshev nodes to estimate the roots of P_n
!
! do i=nn+1,n
!    j=n-i+1
!    xslege(i) = cos(pi*(j+0.0d0)/(n+1.0d0))
! end do


!
!  Perform Newton iterations in order to refine the estimates.
!
do iter = 1,maxiters

!
!  Evaluate the Legendre polynomial of degree n at each point; save
!  the values of the derivative in the whts array.
!
do i=nn+1,n
call legeder(n,xslege(i),pols(i),whtslege(i))
end do

!
!  Perform one Newton iteration
!
pols(nn+1:n)     = pols(nn+1:n)/whtslege(nn+1:n)
xslege(nn+1:n)   = xslege(nn+1:n) - pols(nn+1:n)

if(norm2(pols(nn+1:n)) < eps0) then
exit
endif

end do


if (iter == maxiters)  then
print *,"legemdre_quad: newton iterations failed!"
stop
end if

!
! Compute the weights using the derivatives we stored above.
!
do j=nn+1,n
x           = xslege(j)
dd          = 2.0d0/(1.0d0-x**2)
whtslege(j) = dd/(whtslege(j)**2) 
end do

! do j=nn+1,n
! x = xslege(j)
! call legeder(n-1,x,val,der)
! whtslege(j) = 2*(1-x**2)/(n**2*val**2)
! end do

!
! Reflect the quadrature nodes.
!
do j=1,nn
xslege(j)   = -xslege(n-j+1)
whtslege(j) = whtslege(n-j+1)
end do

!
! Handle the root at 0 if n is odd.
!

if (ifodd .eq. 1) then
x0          = 0
call legeder(n,x0,pol,der)
xslege(nn+1)   = x0
whtslege(nn+1) = 2.0d0/(der**2) 
endif

end subroutine


subroutine legendre_interp1(n,xslege,whtslege,vals,x,valout)
implicit double precision (a-h,o-z)
double precision           :: xslege(n), whtslege(n)
double complex             :: vals(:), valout
!
!  Use barycentric Lagrange interpolation to evaluate a complex-valued expansion of
!  the form (1) at a specified points.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    (xslege,whtslege) - the nodes and weights of the n-point Guass-Legendre
!      quadrature rule
!    vals - the *scaled* values of (1) at the n Gauss-Legendre nodes
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    valout - the value of (1) at the point x
!
double precision :: whts(n)
double complex   :: dsum1

eps0  = epsilon(0.0d0)

dsign = 1.0d0
do i=1,n
dd      = x - xslege(i)
if ( abs(dd) .lt. eps0*10) then
valout = vals(i)
!/sqrt(whtslege(i))
return
endif
whts(i) = dsign/dd * sqrt((1.0d0-xslege(i)**2) * whtslege(i) )
dsign   = -dsign
end do


valout = dot_product(whts,vals) / sum(whts)

end subroutine


subroutine legendre_interp2(n,xslege,whtslege,vals,x,valout)
implicit double precision (a-h,o-z)
double precision           :: xslege(n), whtslege(n)
double precision           :: vals(:), valout
!
!  Use barycentric Lagrange interpolation to evaluate a real-valued expansion of
!  the form (1) at a specified points.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    (xslege,whtslege) - the nodes and weights of the n-point Guass-Legendre
!      quadrature rule
!    vals - the *scaled* values of (1) at the n Gauss-Legendre nodes
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    valout - the value of (1) at the point x
!
double precision :: whts(n)
double precision :: dsum1

eps0  = epsilon(0.0d0)
dsign = 1.0d0
do i=1,n
dd      = x - xslege(i)
if ( abs(dd) .lt. eps0) then
valout = vals(i)
return
endif
whts(i) = dsign/dd * sqrt((1.0d0-xslege(i)**2)*whtslege(i))
dsign   = -dsign
end do

valout = dot_product(whts,vals) / sum(whts)

end subroutine

subroutine legendre_interp3(n,xslege,whtslege,vals,x,valsout)
implicit double precision (a-h,o-z)
double precision           :: xslege(n), whtslege(n)
double precision           :: vals(:,:), valsout(:)
!
!  Use barycentric Lagrange interpolation to evaluate a collection of real-valued
!  expansions of the form (1) at a specified point.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    (xslege,whtslege) - the nodes and weights of the n-point Guass-Legendre
!      quadrature rule
!    vals - a matrix whose jth column gives the *scaled* values (3) of the jth
!      input expanion
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    valsout - an array whose jth entry gives the value of the jth input expansion
!      at the point x
!
double precision              :: whts(n)
!double precision, allocatable :: sums(:)


l = size(vals,2)
!allocate(sums(l))

! compute the barcyentric weights and sum them
eps0  = epsilon(0.0d0)
dsign = 1.0d0
do i=1,n
dd      = x - xslege(i)
if ( abs(dd) .lt. eps0) then
valsout(1:l) = vals(i,1:l)
!/sqrt(whtslege(i))
return
endif
whts(i) = dsign/dd * sqrt((1.0d0-xslege(i)**2)*whtslege(i))
dsign   = -dsign
end do

valsout = matmul(whts,vals) / sum(whts)

end subroutine


subroutine legendre_interp4(n,xslege,whtslege,vals,x,valsout)
implicit double precision (a-h,o-z)
double precision           :: xslege(n), whtslege(n)
double complex             :: vals(:,:), valsout(:)
!
!  Use barycentric Lagrange interpolation to evaluate a collection of complex-valued
!  expansions of the form (1) at a specified point.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    (xslege,whtslege) - the nodes and weights of the n-point Guass-Legendre
!      quadrature rule
!    vals - a matrix whose jth column gives the *scaled* values (3) of the jth
!      input expanion
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    valsout - an array whose jth entry gives the value of the jth input expansion
!      at the point x
double precision              :: whts(n)
!double complex, allocatable   :: sums(:)

l = size(vals,2)
!allocate(sums(l))

! compute the barcyentric weights and sum them
eps0  = epsilon(0.0d0)
dsign = 1.0d0
do i=1,n
dd      = x - xslege(i)
if ( abs(dd) .lt. eps0) then
valsout = vals(i,:)
!/sqrt(whtslege(i))
return
endif
whts(i) = dsign/dd * sqrt((1.0d0-xslege(i)**2)*whtslege(i))
dsign   = -dsign
end do

! compute the sums
valsout = matmul(whts,vals) / sum(whts)


end subroutine


subroutine legendre_eval1(n,coefs,x,val)
implicit double precision (a-h,o-z)
double complex                 :: coefs(:), val
!
!  Evaluate a single complex-valued expansion of the form (1) given the vector
!  (2) of its expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - an array specifying the coefficients
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the expansion at the point x
! 

double precision :: pols(n)
double complex   :: p0, p1, p2

! call leges(n,x,pols)
! val = dot_product(pols,coefs)


!
!  Clenshaw algorithm 
!


p2 = coefs(n)
p1 = coefs(n-1) + 17.0d0/9.0d0 * coefs(n)

do k=n-3,1,-1
p0 = coefs(k+1) + (2*k+1.0d0)/(k+1.0d0)*x*p1 - (k+1.0d0)/(k+2.0d0)*p2
p2 = p1
p1 = p0
end do

val = coefs(1) + x*p1 - p2/2

end subroutine


subroutine legendre_eval2(n,coefs,x,val)
implicit double precision (a-h,o-z)
double precision                 :: coefs(:), val
!
!  Evaluate a single real-valued expansion of the form (1) given the vector
!  (2) of its expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - an array specifying the coefficients
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the expansion at the point x
! 

double precision :: pols(n)

! call leges(n,x,pols)
! val = dot_product(pols,coefs)

!
!  Clenshaw algorithm 
!


p2 = coefs(n)
p1 = coefs(n-1) + 17.0d0/9.0d0 * coefs(n)

do k=n-3,1,-1
p0 = coefs(k+1) + (2*k+1.0d0)/(k+1.0d0)*x*p1 - (k+1.0d0)/(k+2.0d0)*p2
p2 = p1
p1 = p0
end do

val = coefs(1) + x*p1 - p2/2

end subroutine


subroutine legendre_eval3(n,nfuns,coefs,x,vals)
implicit double precision (a-h,o-z)
double precision                 :: coefs(:,:), vals(:)
!
!  Evaluate a collection of real-valued expansions of the form (1) given
!  their expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    nfuns - the number of expansions to evaluate
!    coefs - a matrix whose jth column gives the coefficient vector (2)
!      for the jth input expansion
!    x - the point at which to evaluate the expansions
!
!  Output parameters:
!    vals - the jth entry of this array will give the value of the jth
!       input expansion
! 


!double precision :: p0(nfuns), p1(nfuns), p2(nfuns)

! double precision :: pols(n+1)
! call leges(n,x,pols)
! vals = matmul(pols,coefs)


!
!  Vectorized Clenshaw algorithm 
!

! p2 = coefs(n,:)
! p1 = coefs(n-1,:) + 17.0d0/9.0d0 * coefs(n,:)

! do k=n-3,1,-1
! p0 = coefs(k+1,:) + (2*k+1.0d0)/(k+1.0d0)*x*p1 - (k+1.0d0)/(k+2.0d0)*p2
! p2 = p1
! p1 = p0
! end do

! vals = coefs(1,:) + x*p1 - 0.5d0*p2

!
!
!

do j=1,nfuns
p2 = coefs(n,j)
p1 = coefs(n-1,j) + 17.0d0/9.0d0 * coefs(n,j)

do k=n-3,1,-1
p0 = coefs(k+1,j) + (2*k+1.0d0)/(k+1.0d0)*x*p1 - (k+1.0d0)/(k+2.0d0)*p2
p2 = p1
p1 = p0
end do

vals(j) = coefs(1,j) + x*p1 - 0.5d0*p2
end do

end subroutine


subroutine legendre_eval4(n,nfuns,coefs,x,vals)
implicit double precision (a-h,o-z)
double complex                  :: coefs(:,:), vals(:)
!
!  Evaluate a collection of complex-valued expansions of the form (1) given
!  their expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - a matrix whose jth column gives the coefficient vector (2)
!      for the jth input expansion
!    x - the point at which to evaluate the expansions
!
!  Output parameters:
!    vals - the jth entry of this array will give the value of the jth
!       input expansion
! 

!double complex :: p0(nfuns), p1(nfuns), p2(nfuns)

! double precision :: pols(n)

! call leges(n,x,pols)
! vals = matmul(pols,coefs)


!
!  Vectorized Clenshaw algorithm 
!


! p2 = coefs(n,:)
! p1 = coefs(n-1,:) + 17.0d0/9.0d0 * coefs(n,:)

! do k=n-3,1,-1
! p0 = coefs(k+1,:) + (2*k+1.0d0)/(k+1.0d0)*x*p1 - (k+1.0d0)/(k+2.0d0)*p2
! p2 = p1
! p1 = p0
! end do

! vals = coefs(1,:) + x*p1 - 0.5d0*p2

!
!  Clenshaw algorithm
!

do j=1,nfuns

p2 = coefs(n,j)
p1 = coefs(n-1,j) + 17.0d0/9.0d0 * coefs(n,j)

do k=n-3,1,-1
p0 = coefs(k+1,j) + (2*k+1.0d0)/(k+1.0d0)*x*p1 - (k+1.0d0)/(k+2.0d0)*p2
p2 = p1
p1 = p0
end do

vals(j) = coefs(1,j) + x*p1 - 0.5d0*p2

end do

end subroutine


subroutine legendre_evalder1(n,coefs,x,val,der)
implicit double precision (a-h,o-z)
double complex                 :: coefs(:), val, der
!
!  Evaluate an expansion of the form (1) and its derivative given the expansion
!  coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - an array specifying the coefficients
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the expansion
! 

double precision :: pols(n), ders(n)

call legeders(n,x,pols,ders)

val = dot_product(pols,coefs)
der = dot_product(ders,coefs)

end subroutine


subroutine legendre_evalder2(n,coefs,x,val,der)
implicit double precision (a-h,o-z)
double precision                 :: coefs(:)
!
!  Evaluate an expansion of the form (1) and its derivative given the expansion
!  coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - an array specifying the coefficients
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the expansion
! 

double precision :: pols(n), ders(n)

call legeders(n,x,pols,ders)

val = 0
der = 0

val = dot_product(pols,coefs)
der = dot_product(ders,coefs)

end subroutine


subroutine legendre_evalder3(n,nfuns,coefs,x,vals,ders)
implicit double precision (a-h,o-z)
double precision                 :: coefs(:,:), vals(:), ders(:)
!
!  Evaluate a collection of expansions of the form (1) and their derivatives at
!  a specified point x given the vector (3) of coefficients for each expansion.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - a matrix whose jth column gives the expansion coefficients for the
!      jth input expansion
!    x - the point at which to evaluate the expansions and their derivatives
!
!  Output parameters:
!    vals - an array whose jth entry will contain the value of the jth input
!      expansion at x
!    ders - an array whose jth entry will contain the value of the derivative
!      of the jth expansion at x
! 

double precision :: pols(n), ders0(n)

call legeders(n,x,pols,ders0)

vals = matmul(pols,coefs)
ders = matmul(pols,coefs)

end subroutine


subroutine legendre_evalder4(n,nfuns,coefs,x,vals,ders)
implicit double precision (a-h,o-z)
double complex                 :: coefs(:,:), vals(:), ders(:)
!
!  Evaluate a collection of expansions of the form (1) and their derivatives at
!  a specified point x given the vector (3) of coefficients for each expansion.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - a matrix whose jth column gives the expansion coefficients for the
!      jth input expansion
!    x - the point at which to evaluate the expansions and their derivatives
!
!  Output parameters:
!    vals - an array whose jth entry will contain the value of the jth input
!      expansion at x
!    ders - an array whose jth entry will contain the value of the derivative
!      of the jth expansion at x
! 

double precision :: pols(n), ders0(n)

call legeders(n,x,pols,ders0)

vals = matmul(pols,coefs)
ders = matmul(ders0,coefs)

end subroutine


subroutine legendre_coefsmatrix(n,xslege,whtslege,umatr)
implicit double precision (a-h,o-z)
double precision                           :: xslege(n), whtslege(n)
double precision, allocatable, intent(out) :: umatr(:,:)
!
!  Return the (n,n) matrix which takes the values of an expansion of the
!  form (1) at the nodes of the n-point Gauss-Legendre quadrature to the
!  expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    xslege - the nodes of the n-point Gauss-Legendre quadrature rule
!    whtslege - the weights of the n-point Gauss-Legendre rule
!
!  Output parameters:
!    umatr - the (n,n) matrix which takes values to coefficients
! 

double precision, allocatable :: dnorms(:)


allocate(umatr(n,n),dnorms(n))

do i=1,n
call leges(n,xslege(i),umatr(:,i))
end do

do i=1,n
umatr(i,:) = umatr(i,:) * (i-0.5d0) * whtslege
end do

end subroutine


subroutine legendre_interpmatrix(n,xslege,whtslege,m,xsout,ainterp)
implicit double precision (a-h,o-z)
dimension xslege(n), whtslege(n), xsout(n)
double precision, allocatable, intent(out)  :: ainterp(:,:)
!
!  Construct the matrix which takes the values of (1) at the Gauss-Legendre
!  nodes to its values at the nodes of a user-specified quadrature rule.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    (xslege,whtslege) - the nodes and weights of the n-point Guass-Legendre
!      quadrature rule
!    m - the number of points in the output target quadrature rule
!    (xsout,whtsout) - the nodes and weights of the output quadrature rule
!
!  Output parameters:
!    ainterp - the (m,n) interpolation matrix
!

dimension whts(n)

eps0  = epsilon(0.0d0)

allocate(ainterp(m,n))
ainterp = 0

do j=1,m
x     = xsout(j)
idx   = 0

dsign = 1.0d0
do i=1,n
dd      = x - xslege(i)
if ( abs(dd) .lt. eps0) then
idx = i 
exit
endif

whts(i) = dsign/dd * sqrt( ( 1.0d0-xslege(i)**2 ) * whtslege(i) )
dsign   = -dsign
end do

if (idx .ne. 0) then
ainterp(j,idx) = 1
cycle
endif

dsum         = sum(whts)
ainterp(j,:) = whts/dsum

end do

end subroutine

function eye(n)
implicit double precision (a-h,o-z)
double precision, allocatable :: eye(:,:)

allocate(eye(n,n))
eye =0 
do i=1,n
eye(i,i) = 1
end do

end function


subroutine legendre_diffmatrix(n,xslege,whtslege,adiff)
implicit double precision (a-h,o-z)
double precision                           :: xslege(:), whtslege(:)
double precision, allocatable, intent(out) :: adiff(:,:)
!
!  Return the spectral differentiation matrix which takes the scaled
!  values of an expansion of the form (1) at the Gauss-Legendre nodes
!  to the values of its derivatves at the same nodes.
!
!  Input parameters:
!    n - the number of terms in the Legendre expansion
!    (xs,whts) - the n-point Gauss-Legendre quadrature rule
!
!  Output parameters:
!    adiff - the (n,n) spectral differentiation matrix
!

double precision, allocatable :: u(:,:), vals(:,:), pols(:), ders(:)
allocate(adiff(n,n), u(n,n), vals(n,n), pols(0:n-1), ders(0:n-1) )

do i=1,n
x   = xslege(i)
wht = whtslege(i)
call legeders(n,x,u(:,i),vals(i,:))
end do

do i=1,n
u(i,:) = u(i,:) * (i-0.5d0) * whtslege
end do

adiff = matmul(vals,u)

end subroutine



subroutine legendre_intlmatrix(n,xslege,whtslege,aint)
implicit double precision (a-h,o-z)
double precision, intent(in)               :: xslege(:), whtslege(:)
double precision, allocatable, intent(out) :: aint(:,:)
!
!  Construct the (n,n) spectral integration matrix which takes the
!  values of the an expansion f(x) of the form (1) at the Gauss-Legendre
!  nodes to the values of the function
!
!                x
!    g(x) = \int   f(t) dt
!               -1
!
!   at the Gauss-Legendre nodes 
!
!  Input parameters:
!    n - the number of terms in the Legendre expansions
!    (xs,whts) - the nodes and weights of the n-point Legendre quadrature
!
!  Output parameters:
!    aint - the (n,n) spectral integration matrix
!

double precision, allocatable :: u(:,:), vals(:,:), pols(:)

allocate(aint(n,n))
allocate(pols(0:n), u(n,n), vals(n,n) )

!call legendre_coefsmatrix(n,xslege,whtslege,u)

do i=1,n

x   = xslege(i)
wht = whtslege(i)

call leges(n+1,x,pols)

u(:,i) = pols(0:n-1)

vals(i,1) = x+1

do j=1,n-1
vals(i,j+1) = (pols(j+1)-pols(j-1)) / (2*j+1.0d0)
end do

end do

do i=1,n
u(i,:) = u(i,:) * (i-0.5d0) * whtslege
end do

aint = matmul(vals,u)

end subroutine



subroutine legendre_intrmatrix(n,xslege,whtslege,aint)
implicit double precision (a-h,o-z)
double precision, intent(in)               :: xslege(:), whtslege(:)
double precision, allocatable, intent(out) :: aint(:,:)
!
!  Construct the (n,n) spectral integration matrix which takes the
!  scaled values of the an expansion f(x) of the form (1) at the Gauss-Legendre
!  nodes to the values of the function
!
!                x
!    g(x) = \int   f(t) dt
!                1
!
!   at the Gauss-Legendre nodes.
!
!  Input parameters:
!    n - the number of terms in the Legendre expansions
!    (xs,whts) - the nodes and weights of the n-point Legendre quadrature
!
!  Output parameters:
!    aint - the (n,n) spectral integration matrix
!

double precision, allocatable :: u(:,:), vals(:,:), pols(:)

allocate(aint(n,n))
allocate(pols(0:n), u(n,n), vals(n,n) )

!call legendre_coefsmatrix(n,xslege,whtslege,u)

do i=1,n

x   = xslege(i)
wht = whtslege(i)

call leges(n+1,x,pols)

u(:,i) = pols(0:n-1)

vals(i,1) = x-1

do j=1,n-1
vals(i,j+1) = (pols(j+1)-pols(j-1)) / (2*j+1.0d0)
end do

end do

do i=1,n
u(i,:) = u(i,:) * (i-0.5d0) * whtslege
end do

aint = matmul(vals,u)


end subroutine

subroutine legeder(n,x,pol,der)
implicit double precision (a-h,o-z)
!
!  Evaluate the Legendre polynomial of degree n and its derivative
!  at the point x.
!
!  Input parameters:
!    n - the degree of the polynomial to evaluate
!    x - the point at which the polynomial is to be evaluated
!
!  Output parameters:
!    pol - the value of P_n(x)
!    der - the value of P_n'(x)
!


if ( x == 1.0d0) then
pol = 1
der = (n)*(n+1.0d0)/2
return
endif

if ( x == -1.0d0) then
pol = (-1.0d0)**n
der = -pol*(n)*(n+1.0d0)/2
return
endif


if (n == 0) then
pol = 1
der = 0
else if (n == 1) then
pol = x
der = 1.0d0
else
p1 = 1
p2 = x

do j=2,n
   p  = ((2*j-1)*x*p2-(j-1)*p1)/j
   p1 = p2
   p2 = p
end do
!
pol = p

!
! Compute the derivative using another well-known formula.
!

der=n*(x*p2-p1)/(x**2-1.0d0)
endif

end subroutine


subroutine legeqder(n,x,val,der)
implicit double precision (a-h,o-z)
!
!  Evaluate the Legendre function of the second kind  of degree n and its derivative
!  at the point x.
!
!  Input parameters:
!    n - the degree of the polynomial to evaluate
!    x - the point at which the polynomial is to be evaluated
!
!  Output parameters:
!    pol - the value of Q_n(x)
!    der - the value of Q_n'(x)
!


q1 = 0.5d0*(log(1+x)-log(1-x))
q2 = -1.0d0 + x * 0.5d0*(log(1+x)-log(1-x))

d1 = 0.5d0*(1/(1+x) + 1/(1-x))
d2 = 0.5d0*(2*x/(1-x**2)-log(1-x)+log(1+x))

if (n .eq. 0) then
val = q1
der = d1
return
endif

if (n .eq. 1) then
val = q2
der = d2
return
endif


do j=2,n
q = ((2*j-1)*x*q2-(j-1)*q1)/j
d = ((2*j-1)*q2+(2*j-1)*x*d2-(j-1)*d1)/j

q1 = q2
q2 = q

d1 = d2
d2 = d

end do

val = q
der = d


end subroutine


subroutine leges(n,x,pols)
implicit double precision (a-h,o-z)
double precision           :: pols(:)
!
!  Return the values of the Legendre polynomials of the
!  first kind of degrees 0 through n-1 at a specified point.
!
!  Input parameters:
!     n - the number of polynomials to evaluate
!     x - the point at which to evaluate them
!
!  Output parameters:
!     pols - the array of length n containing the values of the first n 
!      Legendre polynomials
!

if (x == 1.0d0) then
do i=1,n
pols(i) = 1.0d0
end do
return
endif

if (x == -1.0d0) then
dsign = 1.0d0
do i=1,n
pols(i) = dsign
dsign   = -dsign
end do
return
endif

pols(1) = 1.0d0
if (n == 1) return
pols(2) = x
if (n == 2) return

do j=2,n-1
pols(j+1) = ((2*j-1)*x*pols(j)-(j-1)*pols(j-1))/j
end do

end subroutine


subroutine legeqs(n,x,vals)
implicit double precision (a-h,o-z)
double precision           :: vals(:)
!
!  Return the values of the first n functions of the second kind at a specified
!  point x. 
!
!  Input parameters:
!     n - the number of polynomials to evaluate
!     x - the point at which to evaluate them
!
!  Output parameters:
!     vals - the array of length n containing the values of the first n 
!

vals(1) = 0.5d0*(log(1+x)-log(1-x))
if (n == 1) return
vals(2) = -1.0d0 + x*0.5d0*(log(1+x)-log(1-x))
if (n == 2) return

do j=2,n-1
vals(j+1) = ((2*j-1)*x*vals(j)-(j-1)*vals(j-1))/j
end do

end subroutine





subroutine legeders(n,x,pols,ders)
implicit double precision (a-h,o-z)
integer                       :: n
double precision              :: x
double precision, intent(out) :: pols(:),ders(:)
!
!  Evaluate the n Legendre polynomials of degree 0 through n-1 at the
!  point x using the standard 3-term recurrence relation.  Return the values
!  of their derivative at the point x as well.
!
!  Input parameters:
!    n - an integer specifying the number of polynomials to be evaluated
!    x - the point at which to evaluate the polynomials and their derivatives 
!
!  Output parameters:
!    pols - the ith entry of this user-supplied array of length n
!      will contain the value of normalized Legendre polynomial of degree
!      i-1 at x
!    ders - the ith entry of this user-supplied array will contain the
!      value of the derivative of the normalized Legendre polynomial at x
!
!
double precision :: pols2(n)


if ( x == 1.0d0) then
do i=1,n
pols(i) = 1.0d0
ders(i) = (i-1)*i/2
end do
return
endif

if ( x == -1.0d0) then
dsign = 1.0d0
do i=1,n
pols(i) = dsign
ders(i) = -dsign*(i-1)*i/2
dsign   = -dsign
end do
return
endif

pols(1)  = 1
ders(1)  = 0

if (n .gt. 1) then
pols(2)  = x
ders(2)  = 1
endif

!
!  Calculate the values of the polynomials
!
do j=2,n-1
pols(j+1)=((2.0d0*j-1.0d0)*x*pols(j)-(j-1.0d0)*pols(j-1))/j
end do

!
!  Compute the derivatives of the polynomials
!
d=x**2.0d0-1.0d0
do j=3,n
ders(j)=(j-1.0d0)*(x*pols(j)-pols(j-1))/d
end do


end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Code for evaluating Legendre functions of noninteger degrees
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine legendre_pdnu(dnu,t,val)
implicit double precision (a-h,o-z)
!
!  Use the recurrence relation satisfied by solutions of Legendre's equation
!  to evaluate P_dnu(t) with t in the interval [-1,1].
!
!  Input parameters:
!    dnu - the degree of the Legendre function to evaluation
!    t - the point at which to evaluate it in the interval [-1,1]
!
!  Output parameters:
!    val - the value of P_dnu(t)
!

if (dnu .le. 2) then
call legendre_ptaylor(dnu,t,val)
return
endif

n     = floor(dnu)
delta = dnu-n

call legendre_ptaylor(delta,t,val0)
call legendre_ptaylor(delta+1.0d0,t,val1)

dd = delta+1

do i=2,n
val  = (2*dd+1)*t*val1 - dd * val0
val  = val / (dd+1)
dd   = dd+1
val0 = val1
val1 = val
end do

end subroutine




subroutine legendre_ptaylor(dnu,t,val)
implicit double precision (a-h,o-z)
!
!  Use a Taylor series in order to evaluate the the Legendre function of the first kind of
!  order dnu at a point t in the interval [-1,1].  The degree dnu must be in the
!  interval [0,2].
!

maxterms = 100

if (dnu < 0 .OR. dnu >  2) then
print *,"legendre_taylor called with dnu not in the range 0 < dnu < 2"
stop
endif

a   = 1
val = 0
dd  = 1

do i=1,maxterms

delta = a*dd
if (abs(delta) .lt. 1.0d-30) exit

val  = val + delta
dd  = dd  * (1-t)/2
a   = a * (-dnu+i-1)*(dnu+i)/(i**2)

end do

end subroutine


end module
