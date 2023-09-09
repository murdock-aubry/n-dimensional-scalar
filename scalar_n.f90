module scalar_n

    use iso_c_binding
    use utils
    use linalg0
    use chebyshev
    use chebpw
    use odesolve_n
    use odetwo

    type       scalar_vars_n
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

        double precision, allocatable                      :: adiffs(:,:,:)
        double precision, allocatable                      :: adiff(:, :)


        double precision, allocatable                      :: aintr(:,:)

    end type   scalar_vars_n

    type       scalar_odedata_n
        integer                                            :: n
        procedure(scalar_fun), pointer, nopass             :: fun
        double precision, pointer                          :: pars(:)
        double precision                                   :: dk
        integer                                            :: ifleft
        double precision                                   :: wa, wb
    end type   scalar_odedata_n

    interface 

    subroutine      scalar_fun(n,t,qs,par)
        implicit double precision (a-h,o-z)
        double precision        :: t
        double complex          :: qs
        double precision        :: dk
    end subroutine

    end interface

    contains

    subroutine scalar_init_n(n,vars,kcheb0,ntest0,maxstack0,maxints0,maxsteps0,maxiters0,&
        maxleviters0)
        implicit double precision (a-h,o-z)
        type(scalar_vars_n)                         :: vars
        integer, optional                           :: kcheb0, ntest0, maxstack0, maxsteps0
        integer, optional                           :: maxints0, maxiters0, maxleviters0
        double precision, allocatable               :: adiff(:, :)
        !
        !  Initialize the structure containing any data needed by the other procedures 
        !  in this module.
        !
        !  If vars is the only parameter passed to this routine, then reasonable defaults
        !  are chosen.
        !
        !  Input parameters:
        !    n - the order of the differential equation
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
        maxints     = 100000
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
        call chebyshev_diffmatrix(kcheb,adiff)
        
        allocate(vars%adiffs(kcheb,kcheb,n))

        vars%adiffs(:, :, 1) = adiff

        do i = 2, n
            vars%adiffs(:, :, i) = matmul(vars%adiffs(:, :, i-1),adiff)
        end do
    end subroutine

    subroutine scalar_levin_n(n, sum_set, sum_coefs, vars, ier, eps, a, b, ifleft, fun, dk, rvals)
        implicit double precision (a-h, o-z)
        type(scalar_vars_n)                     :: vars
        procedure(scalar_fun)                   :: fun
        double complex                          :: rvals(:,:), ima
        integer                                 :: sum_set(:, :)
        double precision                        :: sum_coefs(:)
        integer, allocatable                    :: its(:)

        double complex, allocatable             :: qs0(:), coefs(:)
        double complex, allocatable             :: ws(:, :), w(:)
        double precision, allocatable           :: ts(:)
        double precision, allocatable           :: rcoefs(:), icoefs(:), reigs(:), ieigs(:)
        integer, allocatable                    :: idxs(:)

        double complex, allocatable             :: r(:), rps(:, :), rhs(:), rhs_sum(:)

        double precision, allocatable           :: adiffs(:, :, :)
    
        double precision, allocatable           :: terms(:)
        double complex, allocatable             :: amatr(:, :), delta(:)
        double complex, allocatable             :: prod(:), sum1(:), sum2(:)

        double complex, allocatable             :: roots(:)

    

        kcheb    = vars%kcheb
        ntest    = vars%ntest
        maxiters = vars%maxleviters
        adiffs   = vars%adiffs


        ier      = 0
        ima      = (0.0d0,1.0d0)
        eps0     = epsilon(0.0d0)
        epssol   = eps0*10
        epsnewt  = eps0*1000
        pi       = acos(-1.0d0)


        m = size(sum_set(1, :))


        allocate(prod(kcheb), sum1(kcheb), sum2(kcheb))
        allocate(ts(kcheb), qs0(kcheb), coefs(kcheb))
        allocate(r(kcheb), rps(n, kcheb), rhs(kcheb), rhs_sum(kcheb), delta(kcheb))
        allocate(amatr(kcheb, kcheb))
        allocate(terms(n), its(n))

        ts = (b-a)/2*vars%xscheb + (b+a)/2
        do i = 1, kcheb
            call fun(n, ts(i), qs0(i), dk)
        end do


        coefs = matmul(vars%acoefs,qs0)
        dd1   = maxval(abs(coefs))
        dd2   = maxval(abs(coefs(kcheb-ntest+1:kcheb)))

        if (dd1 .eq. 0) dd1 = 1
        dd    = dd2/dd1
        if (dd .gt. eps) then
            ier = 1024
            return
        endif




        !
        ! Use the eigensolver
        !

        ! call init_random_seed()
        ! allocate(ws(kcheb,n), w(n))
        ! allocate(reigs(n), ieigs(n), idxs(n))
        ! allocate(rcoefs(n), icoefs(n))
        
        ! its = 20000
        ! do i = 1, kcheb

        !     do j = 1, n-1
        !         rcoefs(j) = 0
        !         icoefs(j) = 0
        !     end do

        !     rcoefs(n) = real(qs0(i))
        !     icoefs(n) = imag(qs0(i))

            


        !     call zamvw(n,rcoefs,icoefs,reigs,ieigs,its,flag)



        !     if (flag .ne. 0) then
        !         ier = 256
        !         return
        !     endif

        !     do j = 1, n
        !         w(j) = reigs(j) + ima * ieigs(j)
        !     end do




        !     if (i .eq. 1) then
        !         call scalar_sortw(n,w)
        !         do k = 1, n
        !             idxs(k) = k
        !         end do
        !     else
        !         do k = 1, n
        !             idxs(k) = minloc(abs(w(k) - ws(i-1, :)), 1)
        !         end do
        !     endif


        !     do j = 1, n
        !         ws(i, idxs(j)) = w(j)
        !     end do

        ! end do

        call root_estimates(n, dk, roots)

        !
        ! Now compute all n phase functions
        !
        

        do isol = 1, n
            dd = 2 / (b - a)

            ! Use eigenvalues as initial estimates
            ! r = ws(:, isol)
            r = roots(isol)


            ! Estimate derivatives
            do j = 1, n-1
                rps(j, :) = (dd ** j) * matmul(adiffs(:, :, j), r)
            end do

            ! Run Levin method
            do iter = 1, maxiters
                dd = 2 / (b-a)
                ifdone = 0
                prod = 1
                sum1 = 0
                sum2 = 0

                ! First construct the rhs
                rhs_sum = 0 ! Reset summation value every iteration
                ! exterior sum
                do i = 1, m
                    terms = sum_set(:, i)
                    prod = 1 ! Need to reset product value
                    ! Compute the product nested in sum on RHS
                    do j = 1, n
                        if (j == 1) then
                            prod = r ** sum_set(j, i)
                        else
                            prod = prod * (rps(j-1, :)) ** sum_set(j, i)
                        end if
                    end do
                    rhs_sum = rhs_sum + sum_coefs(i) * prod
                end do

                rhs = -rhs_sum - qs0 

                


                ! Now construct left hand side

                amatr = 0
                sum1 = 0
                
                ! Contribute the first term with D^{(n-1)}
                amatr = dd ** (n - 1) * adiffs(:, :, n-1)

                ! Contribute the final term where there is no differential matrix
                do j = 1, m 
                    prod = 1 ! Need to reset nested product
                    do k = 2, n
                        prod = prod * (rps(k-1, :) ** sum_set(k, j))
                    end do 
                    sum1 = sum1 + sum_set(1, j) * sum_coefs(j) * prod * (r ** (sum_set(1, j) - 1))
                end do

                do l = 1, kcheb
                    amatr(l, l) = amatr(l, l) + sum1(l)
                end do


                ! Do all other terms
                
                do i = 2, n - 1 !outer sum 
                    sum1 = 0 
                    sum2 = 0 ! Need to reset inner summation 
                    do j = 1, m !second sum
                        prod = 1 ! Need to reset nested product
                        do k = 1, n !nested product
                            if (k == i) then
                                prod = prod ! This is the term we don't consider, map via identity
                            else 
                                if (k == 1) then
                                    prod = prod * r ** sum_set(k, j)
                                else 
                                    prod = prod * (rps(k-1, :)) ** sum_set(k, j) 
                                end if
                            end if
                        end do
                        sum1 = sum1 + sum_set(i, j) * sum_coefs(j) * (rps(i-1,:)) ** (sum_set(i, j)-1) * prod
                    end do
                    do l = 1, kcheb
                        amatr(l, :) = amatr(l, :) + sum1(l) * dd**(i-1) * adiffs(l, :, i-1) 
                    end do
                end do 

                ! Solve the system A * delta = rhs for delta
                call linalg0_qrtsolve(epssol,kcheb,kcheb,amatr,delta,rhs)

                ! perform at least one iteration
                if (iter .gt. 1) then
                    dd1 = norm2(abs(delta))
                    dd0 = norm2(abs(r))+1.0d0
                    dd2  = dd1/dd0
                    if (dd2 .lt. eps) ifdone = 1
                endif

                ! Update the current guess and associated derivatives

                r = r + delta

                do ider = 1, n-1
                    if (ider == 1) then
                        rps(ider, :) = dd * matmul(adiffs(:, :, 1), r)
                    else 
                        rps(ider, :) = dd * matmul(adiffs(:, :, 1), rps(ider - 1, :))
                    end if 
                end do

                if (ifdone .eq. 1) exit
            end do
            

            if (iter .gt. maxiters) then
                ier = 4
                return
            endif


            if (ifleft .eq. 1) then
                do ider = 1, n-1
                    if (ider == 1) then
                        rvals(ider,isol) = r(1)
                    else
                        rvals(ider,isol) = rps(ider-1, 1)
                    end if
                end do
            else
                do ider = 1, n-1
                    if (ider == 1) then
                        rvals(ider,isol) = r(kcheb)
                    else
                        rvals(ider,isol) = rps(ider-1,kcheb)
                    end if
                end do
            endif

        end do

    end subroutine 

    subroutine scalar_riccati_n(n, vars, ier, eps, a, b, c, rvals, fun, pars, chebsol, rs, sum_set, sum_coefs)
        ! Final index indicated basis function
        implicit double precision (a-h, o-z)
        type(scalar_vars_n)                                 :: vars
        procedure(scalar_fun)                               :: fun
        double precision                                    :: dk
        double precision, pointer                           :: pars(:)
        double complex                                      :: rvals(:)
        type(chebpw_scheme)                                 :: chebsol
        double complex, allocatable, intent(out)            :: rs(:,:)

        type(c_ptr)                                         :: userptr
        type(scalar_odedata_n), pointer                     :: odedata
        double complex, allocatable                         :: ys(:,:), yders(:,:), rr(:)
        double complex                                      :: ima

        integer                                :: sum_set(:, :)
        double precision                       :: sum_coefs(:)

        ier = 0
        eps0 = epsilon(0.0d0)
        ima = (0.0d0, 1.0d0)

        allocate(odedata)
        userptr = c_loc(odedata)
        odedata%n      = n-1
        odedata%fun   => fun
        odedata%pars  => pars
        odedata%wa     = 0
        odedata%wb     = 0
        odedata%ifleft = -1
        ifit           = 2


        
        
        call odesolve_nonlinear1(vars%odesolve_vars,jer,ifit,eps,a,b,c,rvals(2:n),n-1,&
          scalar_ode_n,userptr,chebsol,ys,yders,sum_set,sum_coefs)

        ! call prinz("", yders)

          
        if (jer .ne. 0) then
            ier = 4
            return
        endif

        call chebpw_info(chebsol,k,nints)

        
        allocate(rs(k*nints,n+1))

        do j = 2, n
            rs(:, j) = ys(:, j-1)
        end do


        rs(:,n+1) = yders(:,n-1)


        call chebpw_int(chebsol,rs(:,2),c,rvals(1),rr)
        rs(:,1) = rr


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

    subroutine scalar_ode_n(n,t,y,f,df,userptr,sum_set,sum_coefs)
        implicit double precision (a-h,p-z)
        type(c_ptr)                             :: userptr
        double complex                          :: y(:), f(:), df(:,:)
        type(scalar_odedata_n), pointer         :: odedata
        double complex                          :: prod1, prod2, sum1, sum2
        double complex, allocatable             :: rps(:)
        double complex                          :: us, q0
        double precision                        :: sum_coefs(:)
        integer                                 :: sum_set(:, :) 
 
        call c_f_pointer(userptr,odedata)

        allocate(rps(n))
        
        ifleft = odedata%ifleft
        wa     = odedata%wa
        wb     = odedata%wb
        phi    = (erf(wa*(t-wb))+1)/2
        
        
        call odedata%fun(n+1,t,us,odedata%pars(1))

        do i = 1, n
            rps(i) = y(i)
        end do

        

   
        if (ifleft .eq. 1) then
                dk = odedata%dk
            
                q0  = phi*dk + (1-phi)*us

        elseif (ifleft .eq. 0) then
                dk = odedata%dk

                q0  = (1-phi)*dk + phi*us
        else
                q0  = us
        endif


        n1 = n+1
        ! call partitions(n1, sum_set, sum_coefs)
        m = size(sum_set(1, :))


        !!!!!! NOTE COEFFS AND N-TUPLES MUST BE IN CORRECT ORDER
        !!!!!! MUST HAVE (0,0,0,....,0,1) FIRST

        ! Construct vector [f] and matrix Df
        !             [         r'          ]               [0   1     0     0   ...    0   ]
        !             |         r''         |               |0   0     1     0   ...    0   |
        !             |         ...         |               |.   .     .     .   .      .   |
        !       [f] = |         ...         |       Df =    |.   .     .     .    .     .   |
        !             |         ...         |               |.   .     .     .     .    .   |
        !             |       r^(n-1)       |               |0   0     0     0   ...    1   |
        !             [RHS of Riccati = f(r)]               [f'  f'' f^(3)  ...  ... f^(n-1)]
        

        if (n1 .gt. 2) then
            do i = 1, n1-2
                f(i) = rps(i+1)
            end do
        end if 


        sum2 = 0

        ! do the first and last terms first since they have largest magnitude
        sum2 = -q0
        sum2 = sum2 - rps(1) ** sum_set(1, m)

        do i = 2, m-1 !Summation loop
            prod1 = 1
            do j = 1, n
                 prod1 = prod1 * rps(j) ** sum_set(j, i)
            end do
            sum2 = sum2 - sum_coefs(i) * prod1
        end do

        f(n1-1) = sum2




        ! Construct matrix Df

        df = 0

        ! Contruct upper off-diagonal
        do i = 1, n1-2
            df(i, i+1) = 1
        end do


        ! Construct final row
        do k = 1, n1-1
            sum1 = 0
            prod1 = 1
            prod2 = 1
            ! Costruct sum for each k value
            do i = 2, m
                ! First contruct nested n-tuple product

                prod1 = 1
                do j = 0, k-1
                    prod1 = prod1 * (sum_set(1, i) - j)
                end do
                prod1 = prod1 * rps(1) ** (sum_set(1, i) - k)

                ! Now contruct other nested product

                prod2 = 1
                do j = 2, n1
                    prod2 = prod2 * rps(j) ** sum_set(j, i)
                end do

                sum1 = sum1 + sum_coefs(i) * prod1 * prod2
            end do

            call factorial(k, k2)

            df(n1-1, k) = -sum1 / k2


        end do

    end subroutine

    subroutine factorial(n, fact)
        integer                     :: n, fact, prod
        prod = 1
        do i = 0, n-1
            prod = prod * (n - i)
        end do
        fact = prod 
    end subroutine

    subroutine partitions(n, sum_set, sum_coefs)
        implicit double precision (a-h, o-z)
        integer, allocatable                        :: sum_set(:, :), sum_set1(:, :), maxcolumns(:)
        double precision, allocatable               :: sum_coefs(:)
        integer                                     :: fact1, fact2, fact3

        maxn = 20
        
        allocate(maxcolumns(20))

        nn1 = maxn * (maxn - 1) / 2 - 1

        open(unit=20, file="partitions.txt")
        open(unit=210, file = "partitions_num.txt")

        do i = 1,maxn - 1
            read(210, *) maxcolumns(i) 
        end do

        m1 = maxcolumns(maxn - 1)
        m2 = maxcolumns(n-1)
        
        allocate(sum_set1(nn1, m1), sum_set(n, m2), sum_coefs(m2))

        do irow = 1, nn1
            read(20, *) (sum_set1(irow, icol), icol= 1, m1)
        end do

        close(20)
        close(210)

        nn2 = n * (n-1) / 2 - 1

        sum_set = sum_set1(nn2 + 1: nn2 + n, 1:m2)
        

        ! Compute coefficients


        do i = 1, m2
            prod = 1
            do j = 1, n
                call factorial(sum_set(j, i), fact1)
                call factorial(j, fact2)
                prod = prod * fact1 * (fact2 ** sum_set(j, i))
            end do
            call factorial(n, fact3)
            sum_coefs(i) = fact3 / prod
        end do

    

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