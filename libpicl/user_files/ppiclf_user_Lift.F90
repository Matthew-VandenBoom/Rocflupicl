!-----------------------------------------------------------------------
!
! Created April 01, 2025
!
! Subroutine for computing the lift terms
!    Lift components first requires computing gas-phase vorticity
!    and particle angular velocity
!
!
! if collisional_flag = 4  Add Saffman and Magnus lift
!
!
!-----------------------------------------------------------------------
!
#include "PPICLF_STD.h"
submodule (ppiclf_m_user_ForceModels) ppiclf_m_user_ForceModels_Lift
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
    module procedure ppiclf_user_Lift_driver
        !
        ! Input:
        !
        ! type(PPICLF_U_t_particle), intent(inout) :: particle
        ! type(ppiclf_u_t_interp), intent(in) :: interp    
        ! integer*4, intent(in) :: iStage
        ! type(PPICLF_t_realNVec), intent(inout) :: lift
        ! integer*4, intent(out) :: ierr
        ! type(PPICLF_t_realNVec), intent(in) :: v
        ! real*8, intent(in) :: vmag, dp, rep, rhof, rmu, rnu
        !
        ! Code:
        !
        lift%vec = 0.0d0

        if (collisional_flag >= 4) then
            call Lift_Saffman(particle, interp, lift, v, vmag, dp, rep, rhof, rmu, rnu)
            call Lift_Magnus (particle, interp, lift, v, vmag, dp, rep, rhof, rmu)
        endif

        ierr = 0
        return
    end procedure ppiclf_user_Lift_driver
    !
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !
    ! Created April 01, 2025
    !
    ! Subroutine for Saffman lift - shear-induced lift
    !
    ! Requires gas-phase vorticity to be computed
    ! Valid for Rep < 50 and omg* < 0.8 (see Loft, "Lift of a spherical
    !    particle subject to vorticity and/or spin", AIAA J., 
    !    Vol. 46,  pp. 801-809, 2008)
    !      
    ! References:
    ! 1) Fundamentals of Dispersed Multiphase Flows (S.Balachandar), Chap.5
    !
    !-----------------------------------------------------------------------
    !
    pure subroutine Lift_Saffman(particle, interp, lift, v, vmag, dp, rep, rhof, rmu, rnu)
        !
        ! Parameters
        !
        type(PPICLF_U_t_particle), intent(inout) :: particle
        type(ppiclf_u_t_interp), intent(in) :: interp    
        type(PPICLF_t_realNVec), intent(inout) :: lift
        type(PPICLF_t_realNVec), intent(in) :: v
        real*8, intent(in) :: vmag, dp, rep, rhof, rmu, rnu

        !
        ! Internal:
        !
        ! real*8 omgx, omgy, omgz
        type(PPICLF_t_realNVec) omg
        real*8 omg_mag, omg_star
        real*8 epi, Jepi
        real*8 d1, d2, d3
        real*8 factor
        ! real*8 elx, ely, elz
        type(PPICLF_t_realNVec) el
        real*8 elm, ielm

        !
        ! Code:
        !
        if (vmag .lt. 1.d-8) return

        ! Compute gas-phase vorticity components and magnitude
        ! omgx = @{USEPARTICLE(particle%rprop%VOR%X)}@
        ! omgy = @{USEPARTICLE(particle%rprop%VOR%Y)}@
        ! omgz = @{USEPARTICLE(particle%rprop%VOR%Z)}@
        omg = interp%VOR
        omg_mag = nvecMagnitude(omg) ! sqrt(omgx*omgx + omgy*omgy + omgz*omgz)

        ! Compute Mei correction
        omg_star = omg_mag*dp/vmag
        epi = sqrt(omg_star/rep)

        d1 = 1.0d0 + tanh(2.5d0*(log10(epi)+0.191d0))
        d2 = 0.667d0 + tanh(6.0d0*epi-1.92d0)
        Jepi = 0.3d0*d1*d2

        factor = 1.615d0*rmu*(dp*dp)*vmag*sqrt(omg_mag/rnu)

        ! Compute lift components
        ! elx = vy*omgz - vz*omgy
        ! ely = vz*omgx - vx*omgz
        ! elz = vx*omgy - vy*omgx
        el%vec(1) = v%vec(2)*omg%vec(3) - v%vec(3)*omg%vec(2)
        el%vec(2) = v%vec(3)*omg%vec(1) - v%vec(1)*omg%vec(3)
        el%vec(3) = v%vec(1)*omg%vec(2) - v%vec(2)*omg%vec(1)

        elm = nvecMagnitude(el) ! sqrt(elx*elx + ely*ely +elz*elz)
        elm = max(1.0d-20,elm)
        ielm = 1.0d0/elm

        ! liftx = liftx + factor*Jepi*elx*ielm
        ! lifty = lifty + factor*Jepi*ely*ielm
        ! liftz = liftz + factor*Jepi*elz*ielm
        lift = lift + el*(factor*Jepi*ielm)


        return
    end subroutine Lift_Saffman
    !
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !
    ! Created April 01, 2025
    !
    ! Subroutine for Magnus lift - lift induced by particle rotation
    !
    ! Requires particle angular velocity to be calculated
    !      
    ! References:
    ! 1) Fundamentals of Dispersed Multiphase Flows (S.Balachandar), Chap.5
    ! 2) Loth and Drogan, "An equation of motion for particles of finite
    ! Reynolds number and size", (2009). 
    !      
    !-----------------------------------------------------------------------
    !
    pure subroutine Lift_Magnus(particle, interp, lift, v, vmag, dp, rep, rhof, rmu)
        !
        ! Parameters
        !
        type(PPICLF_U_t_particle), intent(inout) :: particle
        type(ppiclf_u_t_interp), intent(in) :: interp    
        type(PPICLF_t_realNVec), intent(inout) :: lift
        real*8, intent(in) :: vmag, dp, rep, rhof, rmu
        type(PPICLF_t_realNVec), intent(in) :: v
        !
        ! Internal:
        !
        ! real*8 omgx, omgy, omgz
        type(PPICLF_t_realNVec) omg
        real*8 omg_mag, omg_star
        real*8 epi, CL
        real*8 d1
        real*8 factor
        ! real*8 elx, ely, elz
        type(PPICLF_t_realNVec) el

        !
        ! Code:
        !
        if (vmag .lt. 1.d-8) return

        ! Compute particle angular velocity
        ! omgx = @{USEPARTICLE(particle%y%ang_vel%X)}@
        ! omgy = @{USEPARTICLE(particle%y%ang_vel%Y)}@
        ! omgz = @{USEPARTICLE(particle%y%ang_vel%Z)}@
        omg = particle%y%ang_vel
        omg_mag = nvecMagnitude(omg) ! sqrt(omgx*omgx + omgy*omgy + omgz*omgz)

        ! Correction to lift
        omg_star = omg_mag*dp/vmag
        epi = omg_star
        d1 = 0.675d0 + 0.15d0*(1.0d0 + tanh(0.28d0*(epi-2.0d0)))
        CL = 1.0d0 - d1*tanh(0.18*sqrt(rep))

        factor = 0.125d0*dp*dp*dp*rhof

        ! Compute lift components
        ! elx = vy*omgz - vz*omgy
        ! ely = vz*omgx - vx*omgz
        ! elz = vx*omgy - vy*omgx
        el%vec(1) = v%vec(2)*omg%vec(3) - v%vec(3)*omg%vec(2)
        el%vec(2) = v%vec(3)*omg%vec(1) - v%vec(1)*omg%vec(3)
        el%vec(3) = v%vec(1)*omg%vec(2) - v%vec(2)*omg%vec(1)

        ! liftx = liftx + factor*CL*elx
        ! lifty = lifty + factor*CL*ely
        ! liftz = liftz + factor*CL*elz
        lift = lift + (el*(factor*CL))

        return    
    end subroutine Lift_Magnus
end submodule ppiclf_m_user_ForceModels_Lift
