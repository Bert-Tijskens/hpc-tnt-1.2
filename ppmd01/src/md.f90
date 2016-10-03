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

module mymd
  ! use wprec ! modules too hard for f2py
    implicit none
    integer, parameter :: wp = kind(1.d0) ! double precision identifier
contains

  ! lennard-jones potential as a function of r
    elemental function lj_pot(r)
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
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
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        real(wp), intent(in) :: r2
        real(wp)             :: lj_pot2

        real(wp)             :: rr,rr6

        rr = 1.0d0/r2
        rr6 = rr*rr*rr;
        lj_pot2 = 4.0d0*rr6*(rr6-1.0d0);
    end function

    ! compute the interaction energy of an atom at (rx0ry0,rz0) with atoms at
    ! (rx(i),ry(i),rz(i)), i = 1..n
    function interaction_energy(n,rx,ry,rz,rx0,ry0,rz0)
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        integer, intent(in) ::n
        real(wp), intent(in), dimension(n) :: rx,ry,rz
        real(wp), intent(in) :: rx0,ry0,rz0
        real(wp)             :: interaction_energy

        integer :: i
        real(wp) :: r
        interaction_energy = 0.0
        do i=1,n
            r = (rx(i)-rx0)**2 + (ry(i)-ry0)**2 + (ry(i)-ry0)**2
            interaction_energy = interaction_energy + lj_pot2(r)
        enddo
    end function

  ! the potential energy of a system with particle coordinates (rx,ry,rz)
  ! (this a O(n2) loop)
    function Epot(n,rx,ry,rz)
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        real(wp),intent(in),dimension(n) :: rx,ry,rz
        real(wp)                         :: Epot,r2
        integer :: i,j,n

        Epot = 0.d0
        do i=1,n
            do j=1,i-1
                r2 = (rx(j)-rx(i))**2 + (ry(j)-ry(i))**2 + (rz(j)-rz(i))**2
                Epot = Epot+lj_pot2(r2)
            end do
        end do
    end function

  ! the kinetic energy of a system with particle velocities (vx,vy,vz)
    function Ekin(n,vx,vy,vz)
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        real(wp),intent(in),dimension(n) :: vx,vy,vz
        real(wp)                         :: Ekin
        integer :: i,n

        Ekin = 0.d0
        do i=1,n
            Ekin = Ekin + vx(i)**2 + vy(i)**2 + vz(i)**2
        end do
        Ekin = Ekin*0.5d0
    end function

  !----------------------------------------
  ! subprograms related to particle motion
  !----------------------------------------
  ! in a larger code you should move0 this part to different file and module.

  ! The lennard-jones interaction force between two particles separated by a
  ! distance vector (dx,dy,dz) is computed as
  !     (dx,dy,dz)*lj_force_factor(r2)
  ! with
  !     r2 = dx*dx+dy*dy+dz*dz
  ! The force factor for ij-pair is obtained by differentiating the potential energy of
  ! atom i w.r.t. the position vector of atom j.

    elemental function lj_force_factor(r2)
        ! multiply inter-particle distance vector with this factor to get the force

        ! elemental : http://fortranwiki.org/fortran/show/elemental
        ! The main benefit of elemental procedures is that advance knowledge
        ! that a function is elemental simplifies parallel execution.

        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        real(wp), intent(in) :: r2 !squared distance
        real(wp)             :: lj_force_factor
        real(wp)             :: ri2,ri6

        ri2 = 1.0d0/r2
        ri6 = ri2*ri2*ri2;
        lj_force_factor = 48.0d0*ri2*ri6*(ri6-0.5d0);
    end function
    ! we choose to implement lj_force_factor as a function because the result
    ! is a single value.

    ! the subprogram compute_accelerations below computes n*3 result values and
    ! it is more convenient to implement it as a subroutine
    subroutine compute_accelerations(n,rx,ry,rz,ax,ay,az)
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        integer :: n
        real(wp),dimension(n),intent(in)    :: rx,ry,rz ! positions
        real(wp),dimension(n),intent(inout) :: ax,ay,az ! accelerations

        real(wp)  :: aij,dx,dy,dz
        integer :: i,j

        do i=1,n
          ! update particle i acceleration
            ax(i) = 0d0
            ay(i) = 0d0
            az(i) = 0d0
        enddo
        do i=1,n
            do j=1,i-1
                dx = rx(i)-rx(j)
                dy = ry(i)-ry(j)
                dz = rz(i)-rz(j)
                aij = lj_force_factor( dx**2 + dy**2 + dz**2 )
              ! update particle i acceleration
                ax(i) = ax(i) + aij*dx
                ay(i) = ay(i) + aij*dy
                az(i) = az(i) + aij*dz
              ! update particle j acceleration
                ax(j) = ax(j) - aij*dx
                ay(j) = ay(j) - aij*dy
                az(j) = az(j) - aij*dz
            enddo
        enddo
    end subroutine

    ! velocity verlet after http://en.wikipedia.org/wiki/Verlet_integration
    ! time integration consists of
    !
    !   vv_update_positions  ( [n,] rx,ry,rz,vx,vy,vz,ax,ay,az,dt)
    !   compute_accelerations( [n,] rx,ry,rz,         ax,ay,az)
    !   vv_update_velocities ( [n,]          vx,vy,vz,ax,ay,az,dt)
    !
    ! Note:
    ! . the accelerations are obviously depending on the interaction model, i.c. LJ
    ! . the argument n, representing the length of the arrays, is only
    !   needed in Fortran. In Python a numpy array is a data structure that
    !   knows its length and the f2py automatically extracts it and passes it
    !   to the fortran function

    subroutine vv_update_positions(n,rx,ry,rz,vx,vy,vz,ax,ay,az,dt)
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        integer            ,intent(in)    :: n
        real(wp),dimension(n),intent(inout) :: rx,ry,rz,vx,vy,vz
        real(wp),dimension(n),intent(in)    :: ax,ay,az
        real(wp)             ,intent(in)    :: dt

        ! treat every space dimension in the same way, accessing the arrays contiguously (=efficiently)
        ! don't repeat code
        call vv_update_position(rx,vx,ax,dt)
        call vv_update_position(ry,vy,ay,dt)
        call vv_update_position(rz,vz,az,dt)
        ! here we make efficient use of the fact that elemental functions
        ! can accept array arguments

    contains
        ! a nested subprogram definition :
        ! this is good practice if we need it only inside the parent function
        ! do not expose what the client (our python code) does not need
        elemental subroutine vv_update_position(r,v,a,dt)
            real(wp),intent(inout) :: r,v
            real(wp),intent(in)    :: a,dt

            !halfstep update of the velocity
            v = v + 0.5*a*dt
            !update position
            r = r + dt*v

        end subroutine
    end subroutine

    ! more of the same
    subroutine vv_update_velocities(n,vx,vy,vz,ax,ay,az,dt)
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        integer            ,intent(in)    :: n
        real(wp),dimension(n),intent(inout) :: vx,vy,vz
        real(wp),dimension(n),intent(in)    :: ax,ay,az
        real(wp)             ,intent(in)    :: dt

        call vv_update_velocity(vx,ax,dt)
        call vv_update_velocity(vy,ay,dt)
        call vv_update_velocity(vz,az,dt)
    contains
        elemental subroutine vv_update_velocity(v,a,dt)
            real(wp),intent(inout) :: v
            real(wp),intent(in)    :: a,dt

            !halfstep update of the velocity
            v = v + 0.5*a*dt

        end subroutine
    end subroutine

    ! A routine that wraps time integration over <ntimesteps> time steps using
    ! the standard velocity verlet formulation
    ! It returns the cputime comsumed. (It could also return nothing, in
    ! which case it would become a subroutine).
    function move0(n,rx,ry,rz,vx,vy,vz,ax,ay,az,dt,ntimesteps)
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        integer            ,intent(in)    :: n,ntimesteps
        real(wp),dimension(n),intent(inout) :: rx,ry,rz,vx,vy,vz,ax,ay,az
        real(wp)             ,intent(in)    :: dt
        real(wp)                            :: move0

        real(wp)  :: cput0,cput1
        integer :: it

        call cpu_time(cput0)

        do it = 1,ntimesteps
            call vv_update_positions  (n,rx,ry,rz,vx,vy,vz,ax,ay,az,dt)
            call compute_accelerations(n,rx,ry,rz,         ax,ay,az)
            call vv_update_velocities (n,         vx,vy,vz,ax,ay,az,dt)
        enddo
!       !note that this could be rewritten as:
!        call vv_update_positions  (n,rx,ry,rz,vx,vy,vz,ax,ay,az,dt)
!        do it = 1,ntimesteps-1
!            call compute_accelerations(n,rx,ry,rz,         ax,ay,az)
!            call vv_update_velocities (n,         vx,vy,vz,ax,ay,az,dt)
!            call vv_update_positions  (n,rx,ry,rz,vx,vy,vz,ax,ay,az,dt)
!            ! if the two calls above are merged the two halfstep velocity
!            ! updates
!            !   v = v + 0.5*a*dt
!            ! can be be merged too:
!            !   v = v + a*dt
!            ! thus avoiding an extra loop over all particles
!            ! (of course the critical efficiency concern is the
!            ! compute_accelerations call which is O(N**2) )
!        enddo
!        call compute_accelerations(n,rx,ry,rz,         ax,ay,az)
!        call vv_update_velocities (n,         vx,vy,vz,ax,ay,az,dt)
!       !see subprogram move
        call cpu_time(cput1)
        move0 = cput1-cput0 !return value
    end function

    ! a subroutine merging vv_update_velocities and vv_update_positions
    subroutine vv_update_velocities_positions(n,rx,ry,rz,vx,vy,vz,ax,ay,az,dt)
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        integer             ,intent(in)    :: n
        real(wp),dimension(n),intent(inout) :: rx,ry,rz,vx,vy,vz
        real(wp),dimension(n),intent(in)    :: ax,ay,az
        real(wp)             ,intent(in)    :: dt

        !elemental functions can accept array arguments
        call vv_update_velocity_position_1d(rx,vx,ax,dt)
        call vv_update_velocity_position_1d(ry,vy,ay,dt)
        call vv_update_velocity_position_1d(rz,vz,az,dt)
    contains
        elemental subroutine vv_update_velocity_position_1d(r,v,a,dt)
            real(wp),intent(inout) :: r,v
            real(wp),intent(in)    :: a,dt

           !update velocity
            v = v + a*dt
           !update position
            r = r + v*dt

        end subroutine
    end subroutine

    ! a more efficient move
    function move(n,rx,ry,rz,vx,vy,vz,ax,ay,az,dt,ntimesteps)
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        integer             ,intent(in)    :: n,ntimesteps
        real(wp),dimension(n),intent(inout) :: rx,ry,rz,vx,vy,vz,ax,ay,az
        real(wp)             ,intent(in)    :: dt
        real(wp)                            :: move

        real(wp)  :: cput0,cput1
        integer :: it

        call cpu_time(cput0)

        call vv_update_positions               (n,rx,ry,rz,vx,vy,vz,ax,ay,az,dt)
        do it = 1,ntimesteps-1
            call compute_accelerations         (n,rx,ry,rz,         ax,ay,az)
            call vv_update_velocities_positions(n,rx,ry,rz,vx,vy,vz,ax,ay,az,dt)
        enddo
        call compute_accelerations             (n,rx,ry,rz,         ax,ay,az)
        call vv_update_velocities              (n,         vx,vy,vz,ax,ay,az,dt)

        call cpu_time(cput1)
        move = cput1-cput0 !return value
    end function

end module
