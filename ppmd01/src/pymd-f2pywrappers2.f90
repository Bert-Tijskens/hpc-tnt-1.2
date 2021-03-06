!     -*- f90 -*-
!     This file is autogenerated with f2py (version:2)
!     It contains Fortran 90 wrappers to fortran functions.

      subroutine f2pywrap_mymd_lj_pot (lj_potf2pywrap, r)
      use mymd, only : lj_pot
      real(kind=wp) r
      real(kind=wp) lj_potf2pywrap
      lj_potf2pywrap = lj_pot(r)
      end subroutine f2pywrap_mymd_lj_pot
      subroutine f2pywrap_mymd_lj_pot2 (lj_pot2f2pywrap, r2)
      use mymd, only : lj_pot2
      real(kind=wp) r2
      real(kind=wp) lj_pot2f2pywrap
      lj_pot2f2pywrap = lj_pot2(r2)
      end subroutine f2pywrap_mymd_lj_pot2
      subroutine f2pywrap_mymd_Epot (Epotf2pywrap, n, rx, ry, rz)
      use mymd, only : Epot
      integer n
      real(kind=wp) rx(n)
      real(kind=wp) ry(n)
      real(kind=wp) rz(n)
      real(kind=wp) Epotf2pywrap
      Epotf2pywrap = Epot(n, rx, ry, rz)
      end subroutine f2pywrap_mymd_Epot
      subroutine f2pywrap_mymd_Ekin (Ekinf2pywrap, n, vx, vy, vz)
      use mymd, only : Ekin
      integer n
      real(kind=wp) vx(n)
      real(kind=wp) vy(n)
      real(kind=wp) vz(n)
      real(kind=wp) Ekinf2pywrap
      Ekinf2pywrap = Ekin(n, vx, vy, vz)
      end subroutine f2pywrap_mymd_Ekin
      subroutine f2pywrap_mymd_lj_force_factor (lj_force_factorf2pywrap,&
     & r2)
      use mymd, only : lj_force_factor
      real(kind=wp) r2
      real(kind=wp) lj_force_factorf2pywrap
      lj_force_factorf2pywrap = lj_force_factor(r2)
      end subroutine f2pywrap_mymd_lj_force_factor
      subroutine f2pywrap_mymd_move0 (move0f2pywrap, n, rx, ry, rz, vx, &
     &vy, vz, ax, ay, az, dt, ntimesteps)
      use mymd, only : move0
      integer n
      real(kind=wp) dt
      integer ntimesteps
      real(kind=wp) rx(n)
      real(kind=wp) ry(n)
      real(kind=wp) rz(n)
      real(kind=wp) vx(n)
      real(kind=wp) vy(n)
      real(kind=wp) vz(n)
      real(kind=wp) ax(n)
      real(kind=wp) ay(n)
      real(kind=wp) az(n)
      real(kind=wp) move0f2pywrap
      move0f2pywrap = move0(n, rx, ry, rz, vx, vy, vz, ax, ay, az, dt, n&
     &timesteps)
      end subroutine f2pywrap_mymd_move0
      subroutine f2pywrap_mymd_move (movef2pywrap, n, rx, ry, rz, vx, vy&
     &, vz, ax, ay, az, dt, ntimesteps)
      use mymd, only : move
      integer n
      real(kind=wp) dt
      integer ntimesteps
      real(kind=wp) rx(n)
      real(kind=wp) ry(n)
      real(kind=wp) rz(n)
      real(kind=wp) vx(n)
      real(kind=wp) vy(n)
      real(kind=wp) vz(n)
      real(kind=wp) ax(n)
      real(kind=wp) ay(n)
      real(kind=wp) az(n)
      real(kind=wp) movef2pywrap
      movef2pywrap = move(n, rx, ry, rz, vx, vy, vz, ax, ay, az, dt, nti&
     &mesteps)
      end subroutine f2pywrap_mymd_move
      
      subroutine f2pyinitmymd(f2pysetupfunc)
      use mymd, only : wp
      use mymd, only : compute_accelerations
      use mymd, only : vv_update_positions
      use mymd, only : vv_update_velocities
      use mymd, only : vv_update_velocities_positions
      interface 
      subroutine f2pywrap_mymd_lj_pot (lj_potf2pywrap, lj_pot, r)
      real(kind=wp) lj_pot
      real(kind=wp) r
      real(kind=wp) lj_potf2pywrap
      end subroutine f2pywrap_mymd_lj_pot 
      subroutine f2pywrap_mymd_lj_pot2 (lj_pot2f2pywrap, lj_pot2, r2)
      real(kind=wp) lj_pot2
      real(kind=wp) r2
      real(kind=wp) lj_pot2f2pywrap
      end subroutine f2pywrap_mymd_lj_pot2 
      subroutine f2pywrap_mymd_Epot (Epotf2pywrap, Epot, n, rx, ry, rz)
      real(kind=wp) Epot
      integer n
      real(kind=wp) rx(n)
      real(kind=wp) ry(n)
      real(kind=wp) rz(n)
      real(kind=wp) Epotf2pywrap
      end subroutine f2pywrap_mymd_Epot 
      subroutine f2pywrap_mymd_Ekin (Ekinf2pywrap, Ekin, n, vx, vy, vz)
      real(kind=wp) Ekin
      integer n
      real(kind=wp) vx(n)
      real(kind=wp) vy(n)
      real(kind=wp) vz(n)
      real(kind=wp) Ekinf2pywrap
      end subroutine f2pywrap_mymd_Ekin 
      subroutine f2pywrap_mymd_lj_force_factor (lj_force_factorf2pywrap,&
     & lj_force_factor, r2)
      real(kind=wp) lj_force_factor
      real(kind=wp) r2
      real(kind=wp) lj_force_factorf2pywrap
      end subroutine f2pywrap_mymd_lj_force_factor 
      subroutine f2pywrap_mymd_move0 (move0f2pywrap, move0, n, rx, ry, r&
     &z, vx, vy, vz, ax, ay, az, dt, ntimesteps)
      real(kind=wp) move0
      integer n
      real(kind=wp) dt
      integer ntimesteps
      real(kind=wp) rx(n)
      real(kind=wp) ry(n)
      real(kind=wp) rz(n)
      real(kind=wp) vx(n)
      real(kind=wp) vy(n)
      real(kind=wp) vz(n)
      real(kind=wp) ax(n)
      real(kind=wp) ay(n)
      real(kind=wp) az(n)
      real(kind=wp) move0f2pywrap
      end subroutine f2pywrap_mymd_move0 
      subroutine f2pywrap_mymd_move (movef2pywrap, move, n, rx, ry, rz, &
     &vx, vy, vz, ax, ay, az, dt, ntimesteps)
      real(kind=wp) move
      integer n
      real(kind=wp) dt
      integer ntimesteps
      real(kind=wp) rx(n)
      real(kind=wp) ry(n)
      real(kind=wp) rz(n)
      real(kind=wp) vx(n)
      real(kind=wp) vy(n)
      real(kind=wp) vz(n)
      real(kind=wp) ax(n)
      real(kind=wp) ay(n)
      real(kind=wp) az(n)
      real(kind=wp) movef2pywrap
      end subroutine f2pywrap_mymd_move
      end interface
      external f2pysetupfunc
      call f2pysetupfunc(wp,f2pywrap_mymd_lj_pot,f2pywrap_mymd_lj_pot2,f&
     &2pywrap_mymd_Epot,f2pywrap_mymd_Ekin,f2pywrap_mymd_lj_force_factor&
     &,compute_accelerations,vv_update_positions,vv_update_velocities,f2&
     &pywrap_mymd_move0,vv_update_velocities_positions,f2pywrap_mymd_mov&
     &e)
      end subroutine f2pyinitmymd


