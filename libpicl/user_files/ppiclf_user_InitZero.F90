!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine to set user-defined values at time t=0
!
!-----------------------------------------------------------------------
!

submodule (ppiclf_user) ppiclf_user_InitZero_imp
    use ppiclf_m_user_RFLUdata, only : PPICLF_TIMEBH !, PPICLF_DRUDTPLAG, PPICLF_DRUDTMIXT
    implicit none
    contains

    module procedure ppiclf_user_InitZero
        !
        ! Internal:
        !
        !
        ! Code:
        !
        ppiclf_TimeBH = 0.0d0

        ! ppiclf_drudtMixt = 0.0d0
        ! ppiclf_drudtPlag = 0.0d0
        return
    end procedure ppiclf_user_InitZero
    
    module procedure ppiclf_user_ZeroParticle
        call ZeroSolutionProps(particle%y)
        call ZeroSolutionProps(particle%ydot)
        call ZeroSolutionProps(particle%ydotc)
        call ZeroSolutionProps(particle%y1(1))
        call Zerorprop(particle%rprop)
#ifdef PPICLF_VU
        particle%DRUDTPLAG(:) = 0
        particle%DRUDTMIXT(:) = 0
#endif
    end procedure ppiclf_user_ZeroParticle

    module procedure ppiclf_user_ZeroInterp
        interp%RHOf         = 0.0d0
        interp%PHIp         = 0.0d0
        interp%U%vec        = 0.0d0
        interp%CS           = 0.0d0
        interp%DPDX%vec     = 0.0d0
        interp%SDR%vec      = 0.0d0
        interp%RHSR         = 0.0d0
        interp%PGC%vec      = 0.0d0
        interp%Ft           = 0.0d0
        interp%VOR%vec      = 0.0d0
        interp%P            = 0.0d0
        interp%RHOG%vec     = 0.0d0
        interp%DPVDX%vec    = 0.0d0
        interp%SDO%vec      = 0.0d0
    end procedure ppiclf_user_ZeroInterp

    pure subroutine ZeroSolutionProps(solProps)
        type(PPICLF_t_solutionProps), intent(inout) :: solProps
        solProps%pos%vec        = 0.0d0
        solProps%vel%vec        = 0.0d0
        solProps%T              = 0.0d0
        solProps%ang_vel%vec    = 0.0d0
        solProps%metal          = 0.0d0
        solProps%oxide          = 0.0d0
    end subroutine ZeroSolutionProps

    pure subroutine Zerorprop(rprop)
        type(PPICLF_t_rprop), intent(inout) :: rprop
        rprop%RHOp          =0.0d0
        rprop%Dp            =0.0d0
        rprop%VOLp          =0.0d0
        rprop%JDPi          =0.0d0
        rprop%JDPe          =0.0d0
        rprop%JSPL          =0.0d0
        rprop%JSPT          =0.0d0
        rprop%FLUCTF%vec    =0.0d0
        rprop%WDOT%vec      =0.0d0
        rprop%IDp           =0.0d0
        rprop%BRNT          =0.0d0
        rprop%XIPAR         =0.0d0
        rprop%XIPERP        =0.0d0
        rprop%XIT           =0.0d0
    end subroutine Zerorprop
end submodule ppiclf_user_InitZero_imp

