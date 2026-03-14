module ppiclf_m_user_ForceModels
    use ppiclf_m_types, only : PPICLF_t_realNVec
    use ppiclf_user_particle
    implicit none
    interface
        pure module subroutine ppiclf_user_BR_driver(particle, interp, iStage,burnrate_model,qq,mdot_me,mdot_ox, rpi, vmag, ierr)
            type(PPICLF_U_t_particle), intent(inout) :: particle
            type(ppiclf_u_t_interp), intent(in) :: interp
            integer*4, intent(in) :: iStage, burnrate_model
            integer*4, intent(out) :: ierr
            real*8, intent(inout) :: qq
            real*8, intent(inout) :: mdot_me, mdot_ox
            real*8, intent(in) :: rpi, vmag
        end subroutine ppiclf_user_BR_driver
        pure module subroutine ppiclf_user_Lift_driver(particle, interp, iStage,lift, v, vmag, dp, rep, rhof, rmu, rnu, ierr)
            type(PPICLF_U_t_particle), intent(inout) :: particle
            type(ppiclf_u_t_interp), intent(in) :: interp    
            integer*4, intent(in) :: iStage
            integer*4, intent(out) :: ierr
            type(PPICLF_t_realNVec), intent(inout) :: lift
            type(PPICLF_t_realNVec), intent(in) :: v
            real*8, intent(in) :: vmag, dp, rep, rhof, rmu, rnu
        end subroutine ppiclf_user_Lift_driver
        pure module subroutine ppiclf_user_Torque_driver(particle, interp, iStage,tau,tau_hydro, dp, rhop, rhof, rmass, rmu, rpi, ierr)
            type(PPICLF_U_t_particle), intent(inout) :: particle
            type(ppiclf_u_t_interp), intent(in) :: interp    
            integer*4, intent(in) :: iStage
            integer*4, intent(out) :: ierr
            type(PPICLF_t_realNVec), intent(inout) :: tau
            type(PPICLF_t_realNVec), intent(inout) :: tau_hydro
            real*8, intent(in) :: dp, rhop, rhof, rmass, rmu, rpi
        end subroutine ppiclf_user_Torque_driver
        pure module subroutine ppiclf_user_HT_driver(particle, interp,qq, rkappa, dp, rep, rpr, rmachp, rphif, ierr)
            type(PPICLF_U_t_particle), intent(inout) :: particle
            type(ppiclf_u_t_interp), intent(in) :: interp  
            integer*4, intent(out) :: ierr
            real*8, intent(inout) :: qq
            real*8, intent(in) :: rkappa, dp, rep, rpr, rmachp, rphif
        end subroutine ppiclf_user_HT_driver

        pure module subroutine ppiclf_user_QS_Parmar(mp, phi, re, cd, beta, rmachp, rphip, rep, dp, rmu)
            real*8, intent(inout)   :: mp, phi, re, cd, beta
            real*8, intent(in)      :: rmachp, rphip, rep, dp, rmu
        end subroutine ppiclf_user_QS_Parmar
        pure module subroutine ppiclf_user_QS_ModifiedParmar(mp, phi, re, cd, beta, rmachp, rphip, rep, dp, rmu)
            real*8, intent(inout)   :: mp, phi, re, cd, beta
            real*8, intent(in)      :: rmachp, rphip, rep, dp, rmu
        end subroutine ppiclf_user_QS_ModifiedParmar
        pure module subroutine ppiclf_user_QS_Osnes(mp, phi, re, cd, beta, rmachp, rphip, rep, dp, rmu)
            real*8, intent(inout)   :: mp, phi, re, cd, beta
            real*8, intent(in)      :: rmachp, rphip, rep, dp, rmu
        end subroutine ppiclf_user_QS_Osnes
        pure module subroutine ppiclf_user_QS_Gidaspow(phi, re, cd, beta, rphip, rep, dp, rmu, rphif, rhof, vmag)
            real*8, intent(inout)   :: phi, re, cd, beta
            real*8, intent(in)      :: rphip, rep, dp, rmu, rphif, rhof, vmag
        end subroutine ppiclf_user_QS_Gidaspow

        pure module subroutine ppiclf_user_QS_fluct_Lattanzi(particle, interp, fqs_fluct, UnifRnd, upmean, u2pmean, icpmean, fac, rphip, vmag, rmachp, rmu, dp, rep, phipmean)
            type(PPICLF_U_t_particle), intent(inout) :: particle
            type(ppiclf_u_t_interp), intent(in) :: interp 
            real*8, intent(inout) :: fqs_fluct(3), UnifRnd(6)
            type(PPICLF_t_realNVec), intent(inout) :: upmean, u2pmean
            integer*4, intent(in) :: icpmean
            real*8, intent(in) :: fac, rphip, vmag, rmachp, rmu, dp, rep, phipmean
        end subroutine ppiclf_user_QS_fluct_Lattanzi

        pure module subroutine ppiclf_user_QS_fluct_Osnes(particle, interp, fqs_fluct, UnifRnd, rsg, t_par, mp, phi, re, rem, upmean, u2pmean, &
                icpmean, rphip, vmag, rmachp, rmu, dp, rep, phipmean, fac, fqs, v, xi_par, xi_perp, xi_T)
            ! 
            ! Input:
            type(PPICLF_U_t_particle), intent(inout) :: particle
            type(ppiclf_u_t_interp), intent(in) :: interp 
            real*8, intent(inout) :: fqs_fluct(3), UnifRnd(6), Rsg(3,3), T_par(3), mp, phi, re, rem
            type(PPICLF_t_realNVec), intent(inout) :: upmean, u2pmean
            integer*4, intent(in) :: icpmean
            real*8, intent(in) :: rphip, vmag, rmachp, rmu, dp, rep, phipmean, fac
            type(PPICLF_t_realNVec), intent(in) :: fqs, v
            !
            ! Output:
            real*8, intent(out) :: xi_par, xi_perp, xi_T
        end subroutine ppiclf_user_QS_fluct_Osnes

        pure module subroutine ppiclf_user_VU_Rocflu(particle, interp,fvu, reyL, rnu, dp, rphip, vmag)
            type(PPICLF_U_t_particle), intent(inout) :: particle
            type(ppiclf_u_t_interp), intent(in) :: interp 
            type(PPICLF_t_realNVec), intent(inout) :: fvu
            real*8, intent(in) :: reyL, rnu, dp, rphip, vmag
        end subroutine ppiclf_user_VU_Rocflu
        module subroutine ppiclf_user_ShiftUnsteadyData
        end subroutine ppiclf_user_ShiftUnsteadyData

        pure module subroutine ppiclf_user_UpdatePlag(particle, interp, v, rphif, rhof)
            type(PPICLF_U_t_particle), intent(inout) :: particle
            type(ppiclf_u_t_interp), intent(in) :: interp 
            type(PPICLF_t_realNVec), intent(in) :: v
            real*8, intent(in) :: rphif, rhof
        end subroutine ppiclf_user_UpdatePlag
        ! module subroutine ppiclf_user_prop2plag
        ! end subroutine ppiclf_user_prop2plag

        ! module subroutine ppiclf_user_plag2prop
        ! end subroutine ppiclf_user_plag2prop

        ! Added mass functions
        pure module subroutine ppiclf_user_AM_Parmar(particle, interp, fam, rmass_add, rmachp, rphip, rhof, rphif, v)
            type(PPICLF_U_t_particle), intent(inout) :: particle
            type(ppiclf_u_t_interp), intent(in) :: interp 
            type(PPICLF_t_realNVec), intent(inout) :: fam
            real*8, intent(inout) :: rmass_add
            real*8, intent(in) :: rmachp, rphip, rhof, rphif
            type(PPICLF_t_realNVec), intent(in) :: v
        end subroutine ppiclf_user_AM_Parmar
        pure module subroutine ppiclf_user_AM_Briney_Unary(particle, interp, fam, FamUnary, rmass_add, rmachp, rphif, rhof, v)
            type(PPICLF_U_t_particle), intent(inout) :: particle
            type(ppiclf_u_t_interp), intent(in) :: interp 
            type(PPICLF_t_realNVec), intent(inout) :: fam, FamUnary
            real*8, intent(inout) :: rmass_add
            real*8, intent(in) :: rmachp, rphif, rhof
            type(PPICLF_t_realNVec), intent(in) :: v
        end subroutine ppiclf_user_AM_Briney_Unary

        pure module subroutine ppiclf_user_AM_Briney_Binary(particle, interp, fam, rphip, Wdot_neighbor_mean, nneighbors)
            type(PPICLF_U_t_particle), intent(inout) :: particle
            type(ppiclf_u_t_interp), intent(in) :: interp 
            type(PPICLF_t_realNVec), intent(inout) :: fam
            real*8, intent(in) :: rphip
            type(PPICLF_t_realNVec), intent(in) :: Wdot_neighbor_mean
            integer*4, intent(in) :: nneighbors
        end subroutine ppiclf_user_AM_Briney_Binary
    end interface
end module ppiclf_m_user_ForceModels