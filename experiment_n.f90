module experiment_n_functions
    use utils
    use linalg0
    use chebyshev
    use chebpw
    use legendre
    use odetwo

    use odesolve_n
    use scalar_n

    use odesolve
    use scalar

    contains

    subroutine test_local(eps, dk, n, dtime)
        implicit double precision (a-h,o-z)
        type(scalar_vars_n)                         :: vars
        double complex                              :: ima
        type(chebpw_scheme)                         :: chebsols
        double complex, allocatable                 :: rn(:,:)
        double precision                            :: dk
        double complex, pointer                     :: rvals(:,:)
        double precision, pointer                   :: pars(:)
        double complex, allocatable                 :: phase(:, :, :)
        
        integer, allocatable                        :: sum_set(:, :)
        double precision, allocatable               :: sum_coefs(:)

        double precision                            :: eps
        
        eps1     = 1.0d-12
        eps2     = eps
        ima = (0.0d0, 1.0d0)


        call partitions(n, sum_set, sum_coefs)
        m = size(sum_set(1, :))


        allocate(phase(ntest, n-1, n))
        allocate(rvals(n, n))
        allocate(pars(1))

        ! sum_set is the set of n-tuples which one of the summations is taken over.
        ! It takes on the following structure:
        !               <--- different tuples --->
        !               ^
        !        S =    |   elements of each tuple
        !               _


        call scalar_init_n(n, vars)

        a        = -1.0d0
        b        = 1.0d0
        ifleft   = 1
        a0       = -1.00d0
        b0       = -0.90d0

        pars(1) = dk

        call elapsed(t1)

        call scalar_levin_n(n, sum_set, sum_coefs, vars, ier, eps1, a0, b0, ifleft, testfun2, dk, rvals(2:n, :))


        if (ier .ne. 0) then
            call prini("levin_n ier = ",ier)
            stop
        endif

        if (ifleft .eq. 1) c = a0
        if (ifleft .eq. 0) c = b0
        

        do i = 1, n
            call scalar_riccati_n(n,vars,ier1,eps2,a,b,c,rvals(:, i),testfun2,pars,chebsols,rn, sum_set, sum_coefs)
        end do
        call elapsed(t2)

        dtime = t2 - t1

        ! print *, "time = ", dtime

    end subroutine

    subroutine test_precision(eps, dk, n, ntest, phase, dtime)
        implicit double precision (a-h,o-z)
        type(scalar_vars_n)                         :: vars
        double complex                              :: ima
        type(chebpw_scheme)                         :: chebsols
        double complex, allocatable                 :: rn(:,:)
        double precision                            :: dk
        double complex, pointer                     :: rvals(:,:)
        double precision, pointer                   :: pars(:)
        double complex, allocatable                 :: phase(:, :)
        double complex                              :: val
        
        integer, allocatable                        :: sum_set(:, :)
        double precision, allocatable               :: sum_coefs(:)
        double precision                            :: eps
        
        eps1     = 1.0d-10
        eps2     = eps
        ima = (0.0d0, 1.0d0)
        ntest = 1000
        call partitions(n, sum_set, sum_coefs)
        m = size(sum_set(1, :))


        if (allocated(phase)) deallocate(phase)

        allocate(phase(ntest, n))
        allocate(rvals(n, n))
        allocate(pars(1))
        ! sum_set is the set of n-tuples which one of the summations is taken over.
        ! It takes on the following structure:
        !               <--- different tuples --->
        !               ^
        !        S =    |   elements of each tuple
        !               _
        call scalar_init_n(n, vars)
        a        = -1.0d0
        b        = 1.0d0
        ifleft   = 1
        a0       = -1.00d0
        b0       = -0.90d0
        pars(1) = dk
        call elapsed(t1)

        call scalar_levin_n(n, sum_set, sum_coefs, vars, ier, eps1, a0, b0, ifleft, testfun2, dk, rvals(2:n, :))

        if (ier .ne. 0) then
            call prini("levin_n ier = ",ier)
            stop
        endif

        if (ifleft .eq. 1) c = a0
        if (ifleft .eq. 0) c = b0
        
        do i = 1, n
            call scalar_riccati_n(n,vars,ier1,eps2,a,b,c,rvals(:, i),testfun2,pars,chebsols,rn, sum_set, sum_coefs)

            do j=1,ntest
                t    = a + (b-a)*(j-1.0d0)/(ntest-1.0d0)
                print *, t
                
                call chebpw_interp(chebsols,rn(:,1),t,val)
                
                phase(j, i) = val
            end do

        end do

        call elapsed(t2)
        
        dtime = t2 - t1
        print *, "time = ", dtime

    end subroutine

    subroutine test_lowdim(eps, dk, ntest, err2, err3, err4, dtime)
        implicit double precision (a-h,o-z)
        type(scalar_vars_n)                         :: vars22, vars23, vars24
        type(scalar_vars_t)                         :: vars12, vars13, vars14
        double complex                              :: ima
        type(chebpw_scheme)                         :: chebsol12, chebsol13, chebsol14
        type(chebpw_scheme)                         :: chebsol22, chebsol23, chebsol24
        
        !double precision                            :: dk

        double precision, pointer                   :: pars(:)
        double complex, pointer                     :: pars2(:)

        double complex, allocatable                 :: phase(:, :)
        double complex                              :: val1, val2
        
        integer, allocatable                        :: sum_set2(:, :), sum_set3(:, :), sum_set4(:, :)
        double precision, allocatable               :: sum_coefs2(:), sum_coefs3(:), sum_coefs4(:)
        double precision                            :: eps

        double complex, allocatable                 :: r12(:,:), r13(:,:), r14(:,:)
        double complex, allocatable                 :: r22(:,:), r23(:,:), r24(:,:)

        double precision, allocatable               :: err2(:, :), err3(:, :), err4(:, :)


        double complex, pointer                     :: rvals11(:,:), rvals12(:,:)
        double complex, pointer                     :: rvals21(:,:), rvals22(:,:)
        double complex, pointer                     :: rvals31(:,:), rvals32(:,:)

        
        eps1     = 1.0d-10
        eps2     = eps
        ima = (0.0d0, 1.0d0)
        ntest = 10000
        errmax = 0.0d0

        call partitions(2, sum_set2, sum_coefs2)
        m2 = size(sum_set2(1, :))

        call partitions(3, sum_set3, sum_coefs3)
        m3 = size(sum_set3(1, :))

        call partitions(4, sum_set4, sum_coefs4)
        m4 = size(sum_set4(1, :))

        

        allocate(rvals11(2, 2))
        allocate(rvals12(2, 2))

        allocate(rvals21(3, 3))
        allocate(rvals22(3, 3))

        allocate(rvals31(4, 4))
        allocate(rvals32(4, 4))



        rvals11 = 0
        rvals12 = 0
        rvals21 = 0
        rvals22 = 0
        rvals31 = 0
        rvals32 = 0




        if (allocated(err2)) deallocate(err2)
        if (allocated(err3)) deallocate(err3)
        if (allocated(err4)) deallocate(err4)

        allocate(err2(ntest, 2), err3(ntest, 3), err4(ntest, 4))
        
        allocate(pars(1), pars2(1))

        ! sum_set is the set of n-tuples which one of the summations is taken over.
        ! It takes on the following structure:
        !               <--- different tuples --->
        !               ^
        !        S =    |   elements of each tuple
        !               _
        
        call scalar_init(vars12)
        call scalar_init(vars13)
        call scalar_init(vars14)

        call scalar_init_n(2, vars22)
        call scalar_init_n(3, vars23)
        call scalar_init_n(4, vars24)


        a        = -1.0d0
        b        = 1.0d0
        ifleft   = 1
        a0       = -1.00d0
        b0       = -0.90d0
        pars(1)  = dk
        pars2(1) = dk

        

        call scalar_levin_n(2, sum_set2, sum_coefs2, vars22, ier, eps1, a0, b0, ifleft, testfun2, dk, rvals12(2:2, :))
        call scalar_levin_n(3, sum_set3, sum_coefs3, vars23, ier, eps1, a0, b0, ifleft, testfun2, dk, rvals22(2:3, :))
        call scalar_levin_n(4, sum_set4, sum_coefs4, vars24, ier, eps1, a0, b0, ifleft, testfun2, dk, rvals32(2:4, :))
        
        
        call scalar_levin2(vars12,ier,eps1,a0,b0,ifleft,fun2,pars2,rvals11(2:2,:))
        call scalar_levin3(vars13,ier,eps1,a0,b0,ifleft,fun3,pars2,rvals21(2:3,:))
        call scalar_levin4(vars14,ier,eps1,a0,b0,ifleft,fun4,pars2,rvals31(2:4,:))


        

        if (ier .ne. 0) then
            call prini("levin_n ier = ",ier)
            stop
        endif


        if (ifleft .eq. 1) c = a0
        if (ifleft .eq. 0) c = b0

        
        
        do i = 1, 1
            call scalar_riccati2(vars12,ier1,eps2,a,b,c,rvals11(:,i),fun2,pars2,chebsol12,r12)
            call scalar_riccati_n(2,vars22,ier1,eps2,a,b,c,rvals12(:, i),testfun2,pars,chebsol22,r22, sum_set2, sum_coefs2)

            do j=1,ntest
                t    = a + (b-a)*(j-1.0d0)/(ntest-1.0d0)
                
                call chebpw_interp(chebsol12,r12(:,i),t,val1)
                call chebpw_interp(chebsol22,r22(:,i),t,val2)

                
                err2(j, i) = abs(val1 - val2)
            end do

        end do

        !print *, "n = 2 max err:", maxval(err2(:, 1))

    


        do i = 1, 1
            call scalar_riccati3(vars13,ier1,eps2,a,b,c,rvals21(:,i),fun3,pars2,chebsol13,r13)


            call scalar_riccati_n(3,vars23,ier1,eps2,a,b,c,rvals22(:, i),testfun2,pars,chebsol23,r23, sum_set3, sum_coefs3)


            do j=1,ntest
                t    = a + (b-a)*(j-1.0d0)/(ntest-1.0d0)
                
                call chebpw_interp(chebsol13,r13(:,i),t,val1)
                call chebpw_interp(chebsol23,r23(:,i),t,val2)


                ! if (j == 900) then 
                !     print *, abs(val1 - val2)
                !     stop
                ! end if
                
                err3(j, i) = abs(val1 - val2)
            end do

        end do

        

        !print *, "n = 3 max err:", maxval(err3(:, 1))




        do i = 1, 1
            call scalar_riccati4(vars14,ier1,eps2,a,b,c,rvals31(:,i),fun4,pars2,chebsol14,r14)

            call scalar_riccati_n(4,vars24,ier1,eps2,a,b,c,rvals32(:, i),testfun2,pars,chebsol24,r24, sum_set4, sum_coefs4)


            do j=1,ntest
                t    = a + (b-a)*(j-1.0d0)/(ntest-1.0d0)
                
                call chebpw_interp(chebsol14,r14(:,i),t,val1)
                call chebpw_interp(chebsol24,r24(:,i),t,val2)
                
                err4(j, i) = abs(val1 - val2)
            end do

        end do

        !print *, "n = 4 max err:", maxval(err4(:, 1))

        

    end subroutine

    subroutine plotting(eps, dk, n, ntest, mm, xs, ys1, ys2)
        implicit double precision (a-h,o-z)
        type(scalar_vars_n)                         :: vars
        double complex                              :: ima
        type(chebpw_scheme)                         :: chebsols
        double complex, allocatable                 :: rn(:,:)
        double precision                            :: dk
        double complex, pointer                     :: rvals(:,:)
        double precision, pointer                   :: pars(:)
        double complex, allocatable                 :: phase(:, :, :)
        
        integer, allocatable                        :: sum_set(:, :)
        double precision, allocatable               :: sum_coefs(:)
    
        double precision                            :: eps
        double complex                              :: val


        double precision, allocatable               :: xs(:), ys1(:), ys2(:)
        
        eps1     = 1.0d-12
        eps2     = eps
        ima = (0.0d0, 1.0d0)
    
    
        call partitions(n, sum_set, sum_coefs)
        m = size(sum_set(1, :))
    
    
        allocate(phase(ntest, n-1, n))
        allocate(rvals(n, n))
        allocate(pars(1))
    
    
        call scalar_init_n(n, vars)
    
        a        = -1.0d0
        b        = 1.0d0
        ifleft   = 1
        a0       = -1.00d0
        b0       = -0.90d0
    
        pars(1) = dk
    
        call elapsed(t1)
    
        call scalar_levin_n(n, sum_set, sum_coefs, vars, ier, eps1, a0, b0, ifleft, testfun2, dk, rvals(2:n, :))
    
    
        if (ier .ne. 0) then
            call prini("levin_n ier = ",ier)
            stop
        endif
    
        if (ifleft .eq. 1) c = a0
        if (ifleft .eq. 0) c = b0
        
    
        call scalar_riccati_n(n,vars,ier1,eps2,a,b,c,rvals(:, 1),testfun2,pars,chebsols,rn, sum_set, sum_coefs)

        nn = ntest

        if (allocated(xs)) deallocate(xs)
        if (allocated(ys1)) deallocate(ys1)
        if (allocated(ys2)) deallocate(ys2)
        allocate(xs(nn),ys1(nn),ys2(nn))

        do i=1,nn
            t = a + (b-a)*(i-1.0d0)/(nn-1.0d0)
            call chebpw_interp(chebsols,rn(:,mm),t,val)
            xs(i)  = t
            ys1(i) = real(val)
            ys2(i) = imag(val)
            
            if (abs(ys1(i)) .lt. 1.0d-16) ys1(i)=0
            if (abs(ys2(i)) .lt. 1.0d-16) ys2(i)=0
        end do

    


    end subroutine

    subroutine testfun1(n,t,y,f,df,dk)
        implicit double precision (a-h,p-z)
        integer                    :: n
        double precision           :: par
        double complex             :: y(:), f(:), df(:,:)
        double complex             :: qs, ima
        data ima / (0.0d0,1.0d0) /
        
        call testfun2(n,t,qs,dk)

        do i = 1, n-1
            f(i) = y(i)
            df(i, i+1) = 1
            df(n, i+1) = 0
        end do

        f(n) = 0

        do i = 1, n
            f(n) = f(n) - qs * y(1)
        end do 
        df(n,1) = -qs

    end subroutine

    subroutine testfun2(n,t,qs,dk)
        implicit double precision (a-h,o-z)
        double complex                  :: qs
        double precision                :: dk
        data ima / (0.0d0,1.0d0) /
        
        qs = - dk ** n * (1 + t ** 2)
        
    end subroutine

    ! Functions for comparison test

    subroutine fun4(n,t,qs,pars)
    
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

    subroutine fun3(n,t,qs,pars)
        implicit double precision (a-h,o-z)
        double complex                  :: qs(:)
        double complex, pointer         :: pars(:)
        double complex                  :: dk0, dk1, dk2, ima
        data ima / (0.0d0,1.0d0) /
        
        dk     = pars(1)
        
        qs(1)  = -dk ** 3 * (1 + t ** 2)
        qs(2)  = 0
        qs(3)  = 0
    end subroutine

    subroutine fun2(n,t,qs,pars)
        implicit double precision (a-h,o-z)
        double complex                  :: qs(:)
        double complex, pointer         :: pars(:)
        double complex                  :: dk0, dk1, dk2, ima
        data ima / (0.0d0,1.0d0) /
        
        dk     = pars(1)
        
        qs(1)  = -dk ** n * (1 + t ** 2)
        qs(2)  = 0
        
    end subroutine

end module

program experimentn
    use utils
    use linalg0
    use chebyshev
    use chebpw
    use odesolve

    use scalar
    use scalar_n
    use experiment_n_functions

    implicit double precision (a-h,o-z)
    integer                       :: ntest
    double precision              :: dtime
    double precision, allocatable :: eps(:)
    double complex, allocatable     :: phase1(:, :), phase2(:, :), vals(:)
    double precision, allocatable   :: ts(:), errmax(:)

    double precision, allocatable   :: err2(:, :), err3(:, :), err4(:, :)
    double precision, allocatable   :: errs(:, :)

    double precision, allocatable   :: dks(:)

    double precision, allocatable   :: times(:), orders(:)

    double precision, allocatable   :: xs(:), ys1(:), ys2(:)

    allocate(eps(20))
    eps = 1.0d0



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                    !
    !                       Quick Test Board                             !
    !                                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! ntest = 10000
        ! eps0 = 1.0d-8
        ! n = 8
        ! m = 7
        ! dk = 300

        ! a = -1.0d0
        ! b = 1.0d0

        ! c = 0
        ! d = -1

        ! call plotting(eps0, dk, n, ntest, m, xs, ys1, ys2)

        

        ! call plot_functions("experiment_n-phase51r.py","experiment_n-phase51r.pdf",1,"t","",0,0,a,b,c,d, &
        !    ntest,xs,ys1,"","","")

        ! call plot_functions("experiment_n-phase51i.py","experiment_n-phase51i.pdf",1,"t","",0,0,a,b,c,d, &
        !    ntest,xs,ys2,"","","")

        

    !                                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                    !
    !                       Test 1: k = 300                              !
    !                                                                    !
    !                        n = 1,..., 9                                !
    !                                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

        ! eps(1) = 1.0d-12 ! precision for n = 2
        ! eps(2) = 1.0d-12 ! precision for n = 3
        ! eps(3) = 1.0d-10
        ! eps(4) = 1.0d-11
        ! eps(5) = 1.0d-10
        ! eps(6) = 1.0d-9
        ! eps(7) = 1.0d-8
        ! eps(8) = 5.0d-8

        ! dk = 300
        ! n = 9
        ! eps0  = epsilon(0.0d0)
        ! allocate(times(n-1), orders(n-1))

        ! do i = 2, n
        !     call test_local(eps(i-1), dk, i, times(i-1))
        !     orders(i-1) = i*1.0d0
        ! end do

        ! times = times

        ! a = 1.0d0
        ! b = (n+1)*1.0d0
        ! c = 0.0d0
        ! d = maxval(times) * 5/4

        ! call plot_functions("experiment_n-times1.py","experiment_n-times1.pdf",1,"n","Time (in milliseconds)",0,0,a,b,c,d, &
        !     n-1,orders,times,"best","k = 300*","")


    !                                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                    !
    !                       Test 2: k = 1000                             !
    !                                                                    !
    !                        n = 1,..., 8                                !
    !                                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! eps(1) = 1.0d-12
        ! eps(2) = 1.0d-12
        ! eps(3) = 1.0d-11
        ! eps(4) = 1.0d-11
        ! eps(5) = 1.0d-9
        ! eps(6) = 1.0d-8
        ! eps(7) = 1.0d-6

        ! dk = 1000
        ! n = 8
        ! eps0  = epsilon(0.0d0)

        ! if (allocated(times)) deallocate(times)
        ! if (allocated(orders)) deallocate(orders)
        ! allocate(times(n-1), orders(n-1))

        ! do i = 2, n
        !     call test_local(eps(i), dk, i, times(i-1))
        !     orders(i-1) = i*1.0d0
        ! end do

        ! times = times * 1000

        ! a = 1.0d0
        ! b = (n+1)*1.0d0
        ! c = 0.0d0
        ! d = maxval(times) * 5/4

        ! call plot_functions("experiment_n-times2.py","experiment_n-times2.pdf",1,"n","Time (in milliseconds)",0,0,a,b,c,d, &
        !     n-1,orders,times,"best","k = 1000*","")



    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                    !
    !                       Test 3: k = 10000                            !
    !                                                                    !
    !                        n = 1,..., 6                                !
    !                                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! eps(1) = 1.0d-12
        ! eps(2) = 1.0d-12
        ! eps(3) = 2.0d-10
        ! eps(4) = 4.0d-11
        ! eps(5) = 5.0d-9


        ! dk = 10000
        ! n = 6
        ! eps0  = epsilon(0.0d0)

        ! if (allocated(times)) deallocate(times)
        ! if (allocated(orders)) deallocate(orders)
        ! allocate(times(n-1), orders(n-1))
        

        ! do i = 2, n
        !     call test_local(eps(i-1), dk, i, times(i-1))
        !     orders(i-1) = i*1.0d0
        ! end do

        ! times = times * 1000

        ! a = 1.0d0
        ! b = (n+1)*1.0d0
        ! c = 0.0d0
        ! d = maxval(times) * 5/4

        ! call plot_functions("experiment_n-times3.py","experiment_n-times3.pdf",1,"n","Time (in milliseconds)",0,0,a,b,c,d, &
        !     n-1,orders,times,"best","k = 10000*","")


    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                    !
    !                       Test 4: k = 100000                           !
    !                                                                    !
    !                        n = 1,..., 5                                !
    !                                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! eps(1) = 1.0d-12
        ! eps(2) = 1.0d-12
        ! eps(3) = 2.0d-11
        ! eps(4) = 4.0d-10


        ! dk = 10000
        ! n = 5
        ! eps0  = epsilon(0.0d0)

        ! if (allocated(times)) deallocate(times)
        ! if (allocated(orders)) deallocate(orders)
        ! allocate(times(n-1), orders(n-1))
        

        ! do i = 2, n
        !     call test_local(eps(i-1), dk, i, times(i-1))
        !     orders(i-1) = i*1.0d0
        ! end do

        ! times = times * 1000

        ! a = 1.0d0
        ! b = (n+1)*1.0d0
        ! c = 0.0d0
        ! d = maxval(times) * 5/4

        ! call plot_functions("experiment_n-times4.py","experiment_n-times4.pdf",1,"n","Time (in milliseconds)",0,0,a,b,c,d, &
        !     n-1,orders,times,"best","k = 100000*","")


    !                                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                    !
    !                       Double precision test                        !
    !                                                                    !
    !                        n = 5, k = 2^8,..., 2^15                    !
    !                                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! iw = 101

        ! a = -1.0d0
        ! b = 1.0d0
        ! ntest = 10000
        ! ii1 = 8
        ! ii2 = 11
        ! n = 5

        ! ! determine whether to use double or extended precision arithemetic
        ! eps0 = 1.0d-15

        ! allocate(vals(ntest))

        
        ! ! phase1 = 0.0d0
        ! ! phase2 = 0.0d0

        ! vals = 0

        

        ! allocate(ts(ntest), errmax(ii2 - ii1))

        ! do i = 1, ntest
        !     ts(i) = a + (b-a) * (i - 1.0d0) / (ntest - 1.0d0)
        ! end do
    
        ! if (eps0 .lt. 1.0d-16) then
        !     eps1 = 1.0d-15

        !     open(iw,FILE="experiment_n.save")

        !     do ii=ii1,ii2
        !         dk = 2.0d0**ii + 1
        
        !         call test_precision(eps1, dk, n, ntest, phase1, dtime)

        !         do i=1,ntest
        !             write(iw,"(3(D44.26,1X))") ts(i), phase1(i, 1)
        !         end do

        !     end do

        !     close(iw)

        ! else
        !     eps1 = 1.0d-12


        !     open(iw,FILE="experiment_n.save",STATUS='old')

        !     do ii=ii1,ii2
        !         dk = 2.0d0**ii + 1

        !         do i=1,ntest
        !             read(iw,"(3(D44.26,1X))") ts(i), vals(i)
        !         end do
                

        !         call test_precision(eps1, dk, n, ntest, phase2, dtime)

        !         errmax(ii - ii1 + 1) = maxval(abs(vals-phase2(:, 1)))

        !         !write (*,"(4(D15.7,1X),3(D15.4,1X))") dk,errmax1,errmax2,dtime0,dtime1,dtime2
            
        !     end do

        !     close(iw)

        !     print *, errmax


        ! endif
    
    !                                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                    !
    !                       Low Dim Comparison                           !
    !                                                                    !
    !                  n = 3 and 4,     k = 2^8,..., 2^15                !
    !                                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        ! eps0 = 1.0d-12

        ! ii1 = 8
        ! ii2 = 15
        ! mm = ii2 - ii1 + 1

        ! allocate(dks(mm))




        ! allocate(errs(ii2 - ii1+1, 2))
        ! errs = 0

        ! do ii = ii1, ii2
        !     nn = ii - ii1 + 1

        !     dk = 2.0d0 ** ii+2
            

        !     if (dk .lt. 2.0d0 ** 12) then
        !         eps0 = 1.0d-12
        !     elseif (dk .lt. 2.0d0 ** 16) then 
        !         eps0 = 1.0d-12
        !     else 
        !         eps0 = 1.0d-12
        !     end if


        !     call test_lowdim(eps0, dk, ntest, err2, err3, err4, dtime)

            
        !     dks(nn) = dk
            
        !     errs(nn, 1) = maxval(abs(err3(:, 1)))
        !     errs(nn, 2) = maxval(abs(err4(:, 1)))
        ! end do

        ! c = 0
        ! d = maxval(errs) * 5.0d0 / 4

        ! a = 2.0d0 ** (ii1 - 1.0d0 / 10)
        ! b = 2.0d0 ** (ii2 + 1.0d0 / 10)
        

        ! call plot_functions("experiment_n-errs.py","experiment_n-errs.pdf",2,"k","Maximum Absolute Error",1,0,a,b,c,d, &
        !     mm,dks,errs,"best","n = 3*n = 4*","")



    !                                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                    !
    !                       Plotting Functions                           !
    !                                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ntest = 10000
        eps0 = 1.0d-10
        n = 8
        m = 6
        dk = 10000

        a = -1.0d0
        b = 1.0d0

        c = 0
        d = -1

        ! n = 5

        call plotting(eps0, dk, 5, ntest, 3, xs, ys1, ys2)
        call plot_functions("experiment_n-phase53.py","experiment_n-phase53.pdf",1,"t","",0,0,a,b,c,d, &
            ntest,xs,ys1,"","","")


        call plotting(eps0, dk, 5, ntest, 5, xs, ys1, ys2)
        call plot_functions("experiment_n-phase55.py","experiment_n-phase55.pdf",1,"t","",0,0,a,b,c,d, &
            ntest,xs,ys1,"","","")


        
        ! n = 9
        
        eps0 = 1.0d-7
        dk = 300

        call plotting(eps0, dk, 9, ntest, 7, xs, ys1, ys2)
        call plot_functions("experiment_n-phase97.py","experiment_n-phase97.pdf",1,"t","",0,0,a,b,c,d, &
            ntest,xs,ys1,"","","")


        call plotting(eps0, dk, 9, ntest, 8, xs, ys1, ys2)
        call plot_functions("experiment_n-phase99.py","experiment_n-phase99.pdf",1,"t","",0,0,a,b,c,d, &
            ntest,xs,ys1,"","","")

        
end program







