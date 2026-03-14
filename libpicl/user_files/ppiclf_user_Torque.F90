!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for computing the torque terms on the
!    RHS of the angular velocity equations
!
!
! if collisional_flag = 1  F = Fn
!                     = 2  F = Fn + Ft + Tt
!                     = 3  F = Fn + Ft + Tt + Th + Tr
! where Tt = collisional torque
!       Th = hydrodynamic torque
!       Tr = rolling torque
!
! Note that the collisional and rolling torques are due to
!     particle-particle interactions and thus are evaulated
!     in ppiclf_user_EvalNearestNeighbor. Therefore, only the
!     hydrodynamic torque is left to be calculated.
!
!-----------------------------------------------------------------------
!
#include "PPICLF_STD.h"
submodule (ppiclf_m_user_ForceModels) ppiclf_m_user_ForceModels_Torque
    ! particle data
    use ppiclf_data, only:
    use ppiclf_m_particledata, only:
    ! grid data
    use ppiclf_data, only:
    use ppiclf_data, only:
    use ppiclf_data, only:
    ! particle options variables
    use ppiclf_data, only:
    use ppiclf_data, only: 
    use ppiclf_data, only: ppiclf_dt
    ! use ppiclf_data, only:
    ! comm variables
    use ppiclf_data, only: 
    ! binning variables
    use ppiclf_data, only: 
    ! ghost particle variables
    use ppiclf_data, only: 
    ! wall support variables
    use ppiclf_data, only:
    ! AngularPeriodic variables (?)(SEE NOTE IN ppiclf_data)
    use ppiclf_data, only:


    use ppiclf_m_user_data
    use ppiclf_m_user_RFLUdata

    use ppiclf_op, only: ppiclf_exittr
    implicit none


    contains
    module procedure ppiclf_user_Torque_driver
        !
        ! Input:
        !
        ! type(PPICLF_U_t_particle), intent(inout) :: particle
        ! type(ppiclf_u_t_interp), intent(in) :: interp    
        ! integer*4, intent(in) :: iStage
        ! integer*4, intent(out) :: ierr
        ! type(PPICLF_t_realNVec), intent(inout) :: tau
        ! type(PPICLF_t_realNVec), intent(inout) :: tau_hydro
        ! real*8, intent(in) :: dp, rhop, rhof, rmass, rmu, rpi
        !
        ! Internal:
        !
      
        ! real*8 taux_undist, tauy_undist, tauz_undist
        type(PPICLF_t_realNVec) tau_undist
        real*8 rmass_local

        !
        ! Code:
        !
        ! taux_hydro = 0.0d0
        ! tauy_hydro = 0.0d0
        ! tauz_hydro = 0.0d0
        tau_hydro%vec = 0.0d0
        ! taux_undist = 0.0d0
        ! tauy_undist = 0.0d0
        ! tauz_undist = 0.0d0
        tau_undist%vec = 0.0d0

        if (collisional_flag >= 3) then
            call Torque_Hydro(particle, interp, tau_hydro, dp, rhop, rhof, rmass, rmu, rpi, ierr)
            if (ierr .ne. 0) return
            call Torque_Undisturbed(particle, interp, tau_undist, dp, rhof)
        endif

        ! taux = taux + taux_hydro + taux_undist
        ! tauy = tauy + tauy_hydro + tauy_undist
        ! tauz = tauz + tauz_hydro + tauz_undist
        tau = tau + tau_hydro + tau_undist

        return
    end procedure ppiclf_user_Torque_driver
    !
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !
    ! Created June 17, 2024
    !
    ! Subroutine for hydrodynamic torque
    !
    !
    !-----------------------------------------------------------------------
    !
    pure subroutine Torque_Hydro(particle, interp, tau_hydro, dp, rhop, rhof, rmass, rmu, rpi, ierr)
        type(PPICLF_U_t_particle), intent(inout) :: particle
            type(ppiclf_u_t_interp), intent(in) :: interp    
        type(PPICLF_t_realNVec), intent(inout) :: tau_hydro
        integer*4, intent(inout) :: ierr
        real*8, intent(in) :: dp, rhop, rhof, rmass, rmu, rpi
        !
        ! Internal:
        !
        ! real*8 omgrx, omgry, omgrz
        type(PPICLF_t_realNVec) omgr
        real*8 omgr_mag
        real*8 Ct1, Ct2, Ct3, Ct
        real*8 reyr, beta, rIp, factor

        !
        ! Code:
        !
        ! Compute relative angular velocity components
        !    and magnitude
        ! omgrx = 0.5d0*(@{USEPARTICLE(ppiclf_parts(i)%rprop%VOR%X)}@) - (@{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%X)}@)
        ! omgry = 0.5d0*(@{USEPARTICLE(ppiclf_parts(i)%rprop%VOR%Y)}@) - (@{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%Y)}@)
        ! omgrz = 0.5d0*(@{USEPARTICLE(ppiclf_parts(i)%rprop%VOR%Z)}@) - (@{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%Z)}@)
        omgr = (interp%VOR - particle%y%ang_vel) * 0.5d0
        omgr_mag = nvecMagnitude(omgr) ! sqrt(omgrx*omgrx + omgry*omgry + omgrz*omgrz)

        ! Particle rotational Reynolds number
        reyr = rhof*dp*dp*omgr_mag/(4.0d0*rmu)
        reyr = max(reyr,0.001d0)

        ! Compute the hydrodynamic torque parameter Ct=Ct(Re_r)
        if (reyr < 1) then
            Ct1 = 0.0d0
            Ct2 = 16.0d0*rpi
            Ct3 = 0.0d0
        elseif (reyr < 10) then
            Ct1 = 0.0d0
            Ct2 = 16.0d0*rpi
            Ct3 = 0.0418d0
        elseif (reyr < 20) then
            Ct1 = 5.32d0
            Ct2 = 37.2d0
            Ct3 = 0.0d0
        elseif (reyr < 50) then
            Ct1 = 6.44d0
            Ct2 = 32.2d0
            Ct3 = 0.0d0
        elseif (reyr < 100) then
            Ct1 = 6.45d0
            Ct2 = 32.1d0
            Ct3 = 0.0d0
        else
            ! call ppiclf_exittr('Re rotational too large$', reyr, 0)
            ierr = 1
            return
        endif

        Ct = Ct1/sqrt(reyr) + Ct2/Reyr + Ct3*reyr

        ! Now compute hydrodynamic torque components
        beta = rhop/rhof
        rIp  = rmass*dp*dp/10.0d0
        factor = rIp*60.0d0*Ct*omgr_mag/(64.0d0*rpi*beta)

        ! taux_hydro = factor*omgrx
        ! tauy_hydro = factor*omgry
        ! tauz_hydro = factor*omgrz
        tau_hydro = omgr * factor


        return
    end subroutine Torque_Hydro
    !
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !
    ! Created April 01, 2025
    !
    ! Subroutine for undisturbed torque
    !
    !
    !-----------------------------------------------------------------------
    !
    pure subroutine Torque_Undisturbed(particle, interp, tau_undist, dp, rhof)
        !
        ! Input:
        !
        type(PPICLF_U_t_particle), intent(in) :: particle
        type(ppiclf_u_t_interp), intent(in) :: interp
        type(PPICLF_t_realNVec), intent(inout) :: tau_undist
        real*8, intent(in) :: dp, rhof
        !
        ! Internal:
        !
        real*8 rIf

        !
        ! Code:
        !

        ! Moment of interia with respect to gas
        rIf = rhof*dp*dp*(particle%rprop%VOLp)/10.0d0

        ! Undisturbed torque component
        ! Written using angular velocity = 0.5*vorticity
        ! taux_undist = 0.5d0*rIf*(@{USEPARTICLE(ppiclf_parts(i)%rprop%SDO%X)}@)
        ! tauy_undist = 0.5d0*rIf*(@{USEPARTICLE(ppiclf_parts(i)%rprop%SDO%Y)}@)
        ! tauz_undist = 0.5d0*rIf*(@{USEPARTICLE(ppiclf_parts(i)%rprop%SDO%Z)}@)
        tau_undist = interp%SDO * rIF * 0.5d0


        return
    end subroutine Torque_Undisturbed
end submodule ppiclf_m_user_ForceModels_Torque
