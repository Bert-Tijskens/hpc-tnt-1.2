! ============================================================================
! Name        : ppmd01.f90
! Author      : E. Tijskens
! Description : Driver program for mymd module
! ============================================================================

program ppmd01
    implicit none
    logical,parameter :: test   =.true.
    logical,parameter :: verbose=.false.

    integer,parameter :: n=29 ! 29
    integer :: m=2**n,k=1,i

    if (test) then
        call test_lj_pot     (verbose)
        call test_util_random(verbose)
    else
        print *, "Skipping tests"
    endif

    open(10,file='ppmd01.txt',status='REPLACE')
    close(10)

    do i=1,1!n-8
        call experiment1(m,k)
        m=m/2
        k=k*2
    enddo

    print *, "*** program ppmd01 finished ***"
end program

subroutine dummy(p)
    use wprec
    real(wp),intent(inout) :: p(:)
    real(wp) :: w
    w    = p(1)
    p(1) = p(2)
    p(3) = w
end subroutine dummy

  ! lennard-jones potential as a function of r
    elemental function lj_pot(r)
        use wprec
        implicit none
        real(wp)             :: lj_pot ! type declaration of return variable
        real(wp), intent(in) :: r      ! type declaration of dummy argument

        real(wp)             :: rr,rr6 ! ! type declaration of dummy argument

        rr = 1.0d0/(r*r)
        rr6 = rr*rr*rr;
        lj_pot = 4.0d0*rr6*(rr6-1.0d0);
    end function

  ! lennard-jones potential as a function of r**2
  ! It is a waste of resources to compute the square root in the interatomic
  ! distance if we need to square it again.
    elemental function lj_pot2(r2)
        use wprec
        implicit none
        real(wp), intent(in) :: r2
        real(wp)             :: lj_pot2

        real(wp)             :: rr,rr6

        rr = 1.0d0/r2
        rr6 = rr*rr*rr;
        lj_pot2 = 4.0d0*rr6*(rr6-1.0d0);
    end function

subroutine experiment1(m,k)
    use wprec
    use util_random
    implicit none
    integer, intent(in) :: m ! length of the position arrays and number of potentials computed
    integer, intent(in) :: k ! number of times that the experiment is repeated

    integer :: im,ik,j(m)

    real(wp) :: x0,y0,z0, xyz0(3),p(3*m+3),r,lj_pot2

    real :: start_time, stop_time, cpu_ordered, cpu_random, cpu_strided

    print *,"m=",m," k=",k," m*k=",m*k
    call cpu_time(start_time)
    call fill_uniform_real(xyz0)
    x0 = xyz0(1)
    y0 = xyz0(2)
    z0 = xyz0(3)
    call fill_uniform_real(p,0._wp,3._wp)
    call fill_uniform_integer(j,1,m)
    call cpu_time(stop_time)
    print *, "setup  : cpu time:", stop_time - start_time, "seconds"

    ! ordered access, SoA stride = 1
    call cpu_time(start_time)
    do ik=1,k
        do im=1,m
            r = (p(im)-x0)*(p(im)-x0) + (p(m+im)-y0)*(p(m+im)-y0) + (p(2*m+im)-z0)*(p(2*m+im)-z0)
            r = lj_pot2(r)
        enddo
!        print *,ik
        if ( r < -100._wp ) then
            ! Prevents smarty-pants compilers from doing “clever” stuff, e.g remove the outer loop
            call dummy(p)
        endif
    enddo
    call cpu_time(stop_time)
    cpu_ordered = stop_time - start_time
    print *, "ordered: cpu time:", cpu_ordered, "seconds"

    ! ordered access, AoS stride = 3
    call cpu_time(start_time)
    do ik=1,k
        do im=1,m
            r = (p(im)-x0)**2 +(p(1+im)-y0)**2 +(p(2+im)-z0)**2
            r = lj_pot2(r)
        enddo
!        print *,ik
        if ( r < -100._wp ) then
            ! Prevents smarty-pants compilers from doing “clever” stuff, e.g remove the outer loop
            call dummy(p)
        endif
    enddo
    call cpu_time(stop_time)
    cpu_strided = stop_time - start_time
    print *, "strided : cpu time:", cpu_strided, "seconds"

    ! ramdom access
    call cpu_time(start_time)
    do ik=1,k
        do im=1,m
            r = (p(j(im))-x0)**2 +(p(j(im)+m)-y0)**2 +(p(j(im)+2*m)-z0)**2
            r = lj_pot2(r)
        enddo
    !        print *,ik
        if ( r < -100._wp ) then
            ! Prevents smarty-pants compilers from doing “clever” stuff, e.g remove the outer loop
            call dummy(p)
        endif
    enddo
    call cpu_time(stop_time)
    cpu_random = stop_time - start_time
    print *, "random : cpu time:", cpu_random, "seconds"
    print *, "random/ordered", cpu_random/cpu_ordered

    open(10,file='ppmd01.txt',status='OLD',position='APPEND')
        write(10,*) m,k,cpu_ordered,cpu_strided,cpu_random
    close(10)

end subroutine experiment1

!
! test routines
!
subroutine test_lj_pot(verbose)
    use wprec
    implicit none
    logical, intent(in) :: verbose
    ! local variables
    real(wp) :: r, pot, pot2, lj_pot,lj_pot2
    r = 1.2_wp
    pot = lj_pot(r)
    pot2 = lj_pot2(r*r)
    if (verbose) then
        print *, pot
        print *, pot2
    endif
    if (pot2==pot) then
        print *, "test_lj_pot OK"
    else
        print *, "test_lj_pot FAILED"
    endif
end subroutine test_lj_pot

subroutine test_util_random(verbose)
    use util_random
    implicit none
    logical, intent(in) :: verbose
    ! local variables
    logical :: ok = .true.
    real(wp) :: x(1),x10(10)
    integer  :: i,j(100)

    ! test fill_uniform_real
    call fill_uniform_real(x)
    if ( x(1)<0_wp .or. x(1)>=1_wp ) then
        ok = .false.
    endif
    call fill_uniform_real(x10)
    do i=1,10
        if (verbose) then
            print *, x10(i)
        endif
        if ( x10(i)<0._wp .or. x10(i)>=1.0_wp ) then
            ok = .false.
        endif
    enddo
    call fill_uniform_real(x10,1.0_wp,1.5_wp)
    do i=1,10
        if (verbose) then
            print *, x10(i)
        endif
        if ( x10(i)<1_wp .or. x10(i)>=1.5_wp ) then
            ok = .false.
        endif
    enddo
    if (.not.ok) then
        print *,"failed at fill_uniform_real "
    endif

    ! test fill_uniform_integer
    call fill_uniform_integer(j,1,10)
    do i=1,100
        if (verbose) then
            print *, j(i)
        endif
        if ( J(i)<1 .or. J(i)>10 ) then
            ok = .false.
        endif
    enddo

    if (ok) then
        print *, "test_util_random OK"
    else
        print *, "test_util_random FAILED"
    end if
end subroutine test_util_random
