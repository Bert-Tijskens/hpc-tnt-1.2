module prec
    integer, parameter :: sp = kind(1.e0) ! single precision identifier
      ! declare single precision variables as
      !   real(sp) :: a
      ! express single precision numbers as
      !   1.23e8_sp, 3.14_sp, 4_sp
    integer, parameter :: dp = kind(1.d0) ! double precision identifier
      ! declare double precision variables as
      !   real(sp) :: a
      ! express double precision numbers as
      !   1.23e8_dp, 3.14_dp, 4_dp
end module prec
