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

    subroutine print_verlet_list(n,verlet_list,m)
      ! print the verlet list entries one per row
        integer                ,intent(in)    :: n,m
        integer ,dimension(m,n),intent(in)    :: verlet_list

        integer :: ia,j,ja,ni

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

        real(wp)  :: aij,dx,dy,dz
        integer :: ia,j,ja,ni

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


end module
