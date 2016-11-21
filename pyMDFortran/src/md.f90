! A module for molecular dynamics using the Lennard-Jones_potential
! in REDUCED UNITS
! source : http://www.phys.ubbcluj.ro/~tbeu/MD/C2_for.pdf

! if the simulation would be for argon:

! quantity      unit                    value(Ar)
!-------------------------------------------------------
! length        sigma                   3.40e-10  m
! energy        epsilon                 1.65e-21 J
! mass          m                       6.69e-26 kg
! time          sigma*sqrt(m/epsilon)   2.17e-12 s      roughly timescale of atomic vibration
! velocity      sqrt(epsilon/m)         1.57e2   m/s
! force         epsilon/sigma           4.85e-12 N
! pressure      epsilon/sigma**3        4.20e7   N/m**2
! temperature   epsilon/kB              120      K

module md
  ! use wprec ! modules too hard for f2py
    implicit none
    integer, parameter :: wp = kind(1.d0) ! double precision identifier
contains

    elemental function lj_pot2(r2)
      ! lennard-jones potential as a function of r**2
      ! It is a waste of resources to compute the square root in the interatomic
      ! distance if we need to square it again.

        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        real(wp), intent(in) :: r2
        real(wp)             :: lj_pot2

        real(wp)             :: rr,rr6

        rr = 1.0d0/r2
        rr6 = rr*rr*rr;
        lj_pot2 = 4.0d0*rr6*(rr6-1.0d0);
    end function

  !----------------------------------------
  ! subprograms related to particle motion
  !----------------------------------------
  ! in a larger code you should move0 this part to different file and module.

  ! The lennard-jones interaction force between two particles separated by a
  ! distance vector (dx,dy,dz) is computed as
  !     (dx,dy,dz)*lj_force_factor2(r2)
  ! with
  !     r2 = dx*dx+dy*dy+dz*dz
  ! The force factor for ij-pair is obtained by differentiating the potential energy of
  ! atom i w.r.t. the position vector of atom j.

  ! elemental : http://fortranwiki.org/fortran/show/elemental
  ! The main benefit of elemental procedures is that knowledge
  ! that a function is elemental simplifies parallel execution.

    elemental function lj_force_factor2(r2)
      ! multiply inter-particle distance vector with this factor to get the force

        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        real(wp), intent(in) :: r2 !squared distance
        real(wp)             :: lj_force_factor2

        real(wp)             :: ri2,ri6

        ri2 = 1.0d0/r2
        ri6 = ri2*ri2*ri2;
        lj_force_factor2 = 48.0d0*ri2*ri6*(ri6-0.5d0);
    end function

    function ppmd01ff(n,rx,ry,rz,ax,ay,az,n_used,n_iter)
      ! Compute all interactions between atom 1 and atoms 2:n and update the acceleations.
      ! (This is the force equivalent of the ppmd01 benchmark)
      ! Repeat n_iter times (for accurate timing results)
      ! Return the cputime
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        integer                ,intent(in)    :: n        ! # atoms
        real(wp),dimension(n)  ,intent(in)    :: rx,ry,rz ! positions
        real(wp),dimension(n)  ,intent(inout) :: ax,ay,az ! accelerations
        integer                ,intent(in)    :: n_iter
        integer                ,intent(in)    :: n_used

        real :: ppmd01ff  ! the cpu time consumed by the subprogram

        real(wp)  :: aij,dx,dy,dz,xi,yi,zi
        real :: cpustart,cpustop
        integer :: ja,iter

        call cpu_time(cpustart)
        do iter=1,n_iter
            xi = rx(1)
            yi = ry(1)
            zi = rz(1)
            !DIR$ SIMD
            do ja=2,n_used
                dx = rx(ja)-xi
                dy = ry(ja)-yi
                dz = rz(ja)-zi
                aij = lj_force_factor2( dx**2 + dy**2 + dz**2 )
              ! update particle i acceleration
                ax(1) = ax(1) + aij*dx
                ay(1) = ay(1) + aij*dy
                az(1) = az(1) + aij*dz
              ! update particle j acceleration
                ax(ja) = ax(ja) - aij*dx
                ay(ja) = ay(ja) - aij*dy
                az(ja) = az(ja) - aij*dz
              ! print *,ia,ja,ax(ia),ay(ia),az(ia),ax(ja),ay(ja),az(ja)
            enddo
        enddo
        call cpu_time(cpustop)
        ppmd01ff = cpustop-cpustart

    end function
!-------------------------------------------------------------------------------
!   Verlet list traversal subprograms

    subroutine print_verlet_list(n,verlet_list,m)
      ! print the verlet list entries one per row
      ! (for testing purposes)
        integer                ,intent(in)    :: n,m
        integer ,dimension(m,n),intent(in)    :: verlet_list

        integer :: ia,j,ni

        print *,'md.f90: print_verlet_list( verlet_linear(m=',m,',n_atoms=',n,') )'
        do ia=1,n
            ni = verlet_list(1,ia)
            do j = 2,ni+1
                print *,ia,verlet_list(j,ia)+1 ! +1 since fortran starts counting from 1 !
            enddo
         enddo
    end subroutine

    subroutine test_interactions_verlet(n,rx,ry,rz,ax,ay,az,verlet_list,m)
      ! Computes the forces between all atom pairs in the verlet list using a dummy force factor
      ! of 1 and updates the accelerations correspondingly. Because of the dummy force factor
      ! the accelerations are just counting the number of interactions.
      !     ax(i) counts the number of entries in the verlet list of atom (i),
      !         I.e. for every interacting atom pair
      !             ax(i) = ax(i) + aij ! aij=1
      !             ax(j) = ax(j) ! unchanged
      !     ay(i) counts the total number of interactions, i.e. the number ofentries in the
      !         verlet list of atom (i) + the number of times that i occurs in the verlet list
      !         of the other atoms.
      !         I.e. for every interacting atom pair
      !             ay(i) = ay(i) + aij ! aij=1
      !             ay(j) = ay(j) + aij
      !     ax(i) counts the number of entries in the verlet list of atom (i) - the number of
      !         times that i occurs in the verlet list of the other atoms.
      !         I.e. for every interacting atom pair
      !             az(i) = az(i) + aij ! aij=1
      !             az(j) = az(j) - aij
      !
      ! It is assumed that the accelerations are appropriately initialized (typically by zeros).
      !
      ! The following post-conditions hold:
      !     ax(i) == len_verlet_list_i
      !     ay(i) == len_verlet_list_i + n_occurences_i
      !     az(i) == len_verlet_list_i - n_occurences_i
      ! where
      !     len_verlet_list_i is the length of the verlet list of atom i,
      !     n_occurences_i is the number of times that atom i occurs all other verlet lists of atoms j!=i.

        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        integer                ,intent(in)    :: n,m
        real(wp),dimension(n)  ,intent(in)    :: rx,ry,rz ! positions
        real(wp),dimension(n)  ,intent(inout) :: ax,ay,az ! accelerations
        integer ,dimension(m,n),intent(in)    :: verlet_list

        real(wp) :: aij
        integer  :: ia,j,ja,ni

        print *,'md.f90: test_interactions_verlet(rx(n),ry(n),rz(n),ax(n),ay(n),az(n),verlet_list(n,m)): n=',n,', m=',m
        do ia=1,n
            ni = verlet_list(1,ia)
            do j = 2,ni+1
                ja = verlet_list(j,ia) + 1 ! +1 since fortran starts counting from 1 !
                aij = 1
              ! update particle i acceleration
                ax(ia) = ax(ia) + aij
                ay(ia) = ay(ia) + aij
                az(ia) = az(ia) + aij
              ! update particle j acceleration
              ! ax(ja) = ax(ja)
                ay(ja) = ay(ja) + aij
                az(ja) = az(ja) - aij
              ! print *,ia,ja,ax(ia),ay(ia),az(ia),ax(ja),ay(ja),az(ja)
            enddo
         enddo
    end subroutine

    function compute_interactions_verlet(n,rx,ry,rz,ax,ay,az,verlet_list,m,n_iter)
      ! Computes the lennard-jones forces between all atom pairs in the verlet list and update
      ! the accelerations.
      ! It is assumed that the accelerations are appropriately initialized (typically by zeros).
      ! Repeat niter times (for accurate timing results)
      ! Return the cputime
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        integer                ,intent(in)    :: n,m
        real(wp),dimension(n)  ,intent(in)    :: rx,ry,rz ! positions
        real(wp),dimension(n)  ,intent(inout) :: ax,ay,az ! accelerations
        integer ,dimension(m,n),intent(in)    :: verlet_list
        integer                ,intent(in)    :: n_iter

        real :: compute_interactions_verlet  ! the cpu time consumed by the subprogram

        real(wp)  :: aij,dx,dy,dz
        real :: cpustart,cpustop
        integer :: ia,j,ja,ni,iter

        call cpu_time(cpustart)
        do iter=1,n_iter
        do ia=1,n
!            print *,ia
            ni = verlet_list(1,ia)
            do j = 2,ni+1
                ja = verlet_list(j,ia) + 1 ! +1 since fortran starts counting from 1 !
                dx = rx(ja)-rx(ia)
                dy = ry(ja)-ry(ia)
                dz = rz(ja)-rz(ia)
                aij = lj_force_factor2( dx**2 + dy**2 + dz**2 )
              ! update particle i acceleration
                ax(ia) = ax(ia) + aij*dx
                ay(ia) = ay(ia) + aij*dy
                az(ia) = az(ia) + aij*dz
              ! update particle j acceleration
                ax(ja) = ax(ja) - aij*dx
                ay(ja) = ay(ja) - aij*dy
                az(ja) = az(ja) - aij*dz
!                print *,ia,ja,ax(ia),ay(ia),az(ia),ax(ja),ay(ja),az(ja)
            enddo
        enddo
        enddo
        call cpu_time(cpustop)
        compute_interactions_verlet = cpustop-cpustart
!        print *,cpustart,cpustop,cpustop-cpustart
    end function

!------------------------------------------------------------------------------
!   The same subprograms for linear verlet list
!   The linear verlet list is a data structure for storing verlet lists in a 1D
!   Array.
!       [n_pairs_i1 j1_i1 j2_i1 .. n_pairs_i2 j1_i2 j2_i2 .. i3 ..]
!   This data structure allows contiguous data access during traversal.

    subroutine print_verlet_linear(n_atoms,verlet_linear,m)
      ! print the verlet list entries one per row
      ! (for testing purposes)
        integer              ,intent(in) :: n_atoms,m
        integer ,dimension(m),intent(in) :: verlet_linear

        integer :: ia,j,ia_pairs,k

        print *,'md.f90: print_verlet_linear( n_atoms=',n_atoms,', verlet_linear(m=',m,') )'
        j=1
        do ia=1,n_atoms
            ia_pairs = verlet_linear(j)
            do k = j+1,j+ia_pairs
                print *,ia,verlet_linear(k)+1 ! +1 since fortran starts counting from 1 !
            enddo
            j = j + 1 + ia_pairs
        enddo
    end subroutine

    subroutine test_interactions_verlet_linear(n_atoms,rx,ry,rz,ax,ay,az,verlet_linear,m)
      ! A test_interactions_verlet but for a linear verlet list data structure

        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        integer                    ,intent(in)    :: n_atoms,m
        real(wp),dimension(n_atoms),intent(in)    :: rx,ry,rz ! positions
        real(wp),dimension(n_atoms),intent(inout) :: ax,ay,az ! accelerations
        integer ,dimension(m)      ,intent(in)    :: verlet_linear

        real(wp)  :: aij
        integer :: ia,j,ja,k,ia_pairs

        print *,'md.f90: test_interactions_verlet_linear(rx(n),ry(n),rz(n),ax(n),ay(n),az(n),verlet_linear(m)): n=',n_atoms,', m=',m
        j=1
        do ia=1,n_atoms
            ia_pairs = verlet_linear(j)
            do k = j+1,j+ia_pairs
                ja = verlet_linear(k)+1 ! +1 since fortran starts counting from 1 !
                aij = 1
              ! update particle i acceleration
                ax(ia) = ax(ia) + aij
                ay(ia) = ay(ia) + aij
                az(ia) = az(ia) + aij
              ! update particle j acceleration
              ! ax(ja) = ax(ja)
                ay(ja) = ay(ja) + aij
                az(ja) = az(ja) - aij
              ! print *,ia,ja,ax(ia),ay(ia),az(ia),ax(ja),ay(ja),az(ja)
            enddo
            j = j + 1 + ia_pairs
        enddo
    end subroutine

    function compute_interactions_verlet_linear(n_atoms,rx,ry,rz,ax,ay,az,verlet_linear,m,n_iter)
      ! Computes the lennard-jones forces between all atom pairs in the linear verlet list and update
      ! the accelerations.
      ! It is assumed that the accelerations are appropriately initialized (typically by zeros).
      ! Repeat n_iter times (for accurate timing results)
      ! Return the cputime
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        integer                    ,intent(in)    :: n_atoms,m
        real(wp),dimension(n_atoms),intent(in)    :: rx,ry,rz ! positions
        real(wp),dimension(n_atoms),intent(inout) :: ax,ay,az ! accelerations
        integer ,dimension(m)      ,intent(in)    :: verlet_linear
        integer                    ,intent(in)    :: n_iter   ! number of iterations, usually 1, >1 to get meaningfull results for cputime measurement

        real :: compute_interactions_verlet_linear  ! the cpu time consumed by this subprogram

        real(wp)  :: aij,dx,dy,dz
        real :: cpustart,cpustop
        integer :: iter,ia,ja,ia_pairs,j,k

        call cpu_time(cpustart)
        do iter=1,n_iter
        j=1 ! 0+1 since fortran starts counting from 1 !
        do ia=1,n_atoms
            ia_pairs = verlet_linear(j)
            !DIR$ SIMD
            do k = j+1,j+ia_pairs
                ja = verlet_linear(k)+1 ! +1 since fortran starts counting from 1 !
                dx = rx(ja)-rx(ia)
                dy = ry(ja)-ry(ia)
                dz = rz(ja)-rz(ia)
                aij = lj_force_factor2( dx**2 + dy**2 + dz**2 )
              ! update particle ia acceleration
                ax(ia) = ax(ia) + aij*dx
                ay(ia) = ay(ia) + aij*dy
                az(ia) = az(ia) + aij*dz
              ! update particle ja acceleration
                ax(ja) = ax(ja) - aij*dx
                ay(ja) = ay(ja) - aij*dy
                az(ja) = az(ja) - aij*dz
!                print *,ia,ja,ax(ia),ay(ia),az(ia),ax(ja),ay(ja),az(ja)
            enddo
            j = j + 1 + ia_pairs
        enddo
        enddo
        call cpu_time(cpustop)
        compute_interactions_verlet_linear = cpustop-cpustart
!        print *,cpustart,cpustop,cpustop-cpustart
    end function

!   untested so far...
    function compute_interactions_verlet_linear2(n_atoms,rx,ry,rz,ax,ay,az,verlet_n,verlet_pairs,m,n_iter)
      ! Computes the lennard-jones forces between all atom pairs in the linear verlet list and update
      ! the accelerations.
      ! In this version we have two arrays for the linear verlet list:
      !   . verlet_n(n_atoms) containing the number of verlet pairs for each atom,
      !   . verlet_pairs(m) containing all the pairs, starting with verlet_n(1) pairs with atom 1,
      !     followed by, verlet_n(2) pairs with atom 2 and so on.
      !
      ! It is assumed that the accelerations are appropriately initialized (typically by zeros).
      ! Repeat n_iter times (for accurate timing results)
      ! Return the cputime
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        integer                    ,intent(in)    :: n_atoms,m
        real(wp),dimension(n_atoms),intent(in)    :: rx,ry,rz ! positions
        real(wp),dimension(n_atoms),intent(inout) :: ax,ay,az ! accelerations
        integer ,dimension(n_atoms),intent(in)    :: verlet_n
        integer ,dimension(m)      ,intent(in)    :: verlet_pairs
        integer                    ,intent(in)    :: n_iter   ! number of iterations, usually 1, >1 to get meaningfull results for cputime measurement

        real :: compute_interactions_verlet_linear2  ! the cpu time consumed by this subprogram

        real(wp)  :: aij,dx,dy,dz
        real :: cpustart,cpustop
        integer :: iter,ia,ja,ia_pairs,j,k

        call cpu_time(cpustart)
        do iter=1,n_iter
        j=1 ! 0+1 since fortran starts counting from 1 !
        do ia=1,n_atoms
            ia_pairs = verlet_n(ia)
            !DIR$ SIMD
            do k = j,j+ia_pairs
                ja = verlet_pairs(k)+1 ! +1 since fortran starts counting from 1 !
                dx = rx(ja)-rx(ia)
                dy = ry(ja)-ry(ia)
                dz = rz(ja)-rz(ia)
                aij = lj_force_factor2( dx**2 + dy**2 + dz**2 )
              ! update particle ia acceleration
                ax(ia) = ax(ia) + aij*dx
                ay(ia) = ay(ia) + aij*dy
                az(ia) = az(ia) + aij*dz
              ! update particle ja acceleration
                ax(ja) = ax(ja) - aij*dx
                ay(ja) = ay(ja) - aij*dy
                az(ja) = az(ja) - aij*dz
!                print *,ia,ja,ax(ia),ay(ia),az(ia),ax(ja),ay(ja),az(ja)
            enddo
            j = j + 1 + ia_pairs
        enddo
        enddo
        call cpu_time(cpustop)
        compute_interactions_verlet_linear2 = cpustop-cpustart
!        print *,cpustart,cpustop,cpustop-cpustart
    end function

end module
