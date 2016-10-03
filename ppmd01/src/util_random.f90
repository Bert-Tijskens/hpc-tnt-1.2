! utilities for random numbers

module util_random
    use wprec
    implicit none
contains
    subroutine fill_uniform_real(x,opt_x1,opt_x2)
        ! fill x with uniform real(wp) numbers in the range [opt_x1,opt_x2[
        real(wp), intent(out) :: x(:)
        real(wp), intent(in), optional :: opt_x1, opt_x2

        real(wp) :: x1, x2
        logical  :: transform = .false.

        call random_number(x)

        if (present(opt_x1)) then
            x1 = opt_x1
            if (x1/=real(0,kind(x))) then
                transform = .true.
            endif
        else
            x1 = real(0,kind(x))
        end if

        if (present(opt_x2)) then
            x2 = opt_x2
            if (x2/=real(1,kind(x))) then
                transform = .true.
            endif
        else
            x2 = real(1,kind(x))
        end if

        if (transform) then
            x2 = x2 - x1
            x = x1 + x*x2
        endif
    end subroutine

    subroutine fill_uniform_integer(i,i1,i2)
        ! fill i with uniform integer numbers in the range [opt_i1,opt_i2]
        integer, intent(out) :: i(:)
        integer, intent(in)  :: i1, i2

        integer :: n
        real(wp) :: x1, x2, xi(size(i))

        n = size(i)
        x1 = real(i1,wp)
        x2 = real(i2+1,wp)

        call fill_uniform_real(xi,x1,x2)
        i = floor(xi)

    end subroutine fill_uniform_integer

end module util_random
