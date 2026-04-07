#include <PPICLF_STD.h>
#define PPICLF_N_PREV_SOLS 1

module PPICLF_user_particle
    use ppiclf_m_types
    use ppiclf_op, only: ppiclf_exittr
    use mpi
    implicit none
    ! type PPICLF_t_realNVec

    ! end type PPICLF_t_realNVec

    

    type PPICLF_t_solutionProps
        type(PPICLF_t_realNVec)                 :: pos
        type(PPICLF_t_realNVec)                 :: vel
        real*8                                  :: T        ! temperature
        type(PPICLF_t_realNVec)                 :: ang_vel
        real*8                                  :: metal
        real*8                                  :: oxide
    end type PPICLF_t_solutionProps

    type PPICLF_t_rprop
        real*8                                  :: RHOp     ! particle density
        real*8                                  :: Dp       ! particle diameter
        real*8                                  :: VOLp     ! particle volume
        real*8                                  :: JDPi     ! no longer used
        real*8                                  :: JDPe     ! no longer used
        real*8                                  :: JSPL     ! = constant 1
        real*8                                  :: JSPT     ! no longer used. ! Might have been repurposed to something else and is actually used?
        type(PPICLF_t_realNVec)                 :: FLUCTF   ! Fluctuating QS force
        type(PPICLF_t_realNVec)                 :: WDOT     ! Relative acceleration
        real*8                                  :: IDp      ! Initial particle diameter
        real*8                                  :: BRNT     ! Total burn time
        real*8                                  :: XIPAR        
        real*8                                  :: XIPERP
        real*8                                  :: XIT

    end type PPICLF_t_rprop

    type PPICLF_U_t_interp
        real*8                                  :: RHOf     ! fluid density
        real*8                                  :: PHIp     ! particle volume fraction
        type(PPICLF_t_realNVec)                 :: U        ! gas velocity vector
        real*8                                  :: CS       ! gas velocity sound speed
        type(PPICLF_t_realNVec)                 :: DPDX     ! gas pressure derivative
        type(PPICLF_t_realNVec)                 :: SDR      ! Gas substantial derivative
        real*8                                  :: RHSR     ! gas density time derivative
        type(PPICLF_t_realNVec)                 :: PGC      ! Gas density velocity
        real*8                                  :: Ft       ! gas temperature
        type(PPICLF_t_realNVec)                 :: VOR      ! gas vorticity
        real*8                                  :: P        ! Fluid pressure
        type(PPICLF_t_realNVec)                 :: RHOG     ! gradient of rho^g (fluid)
        type(PPICLF_t_realNVec)                 :: DPVDX    ! ????
        type(PPICLF_t_realNVec)                 :: SDO      ! DOmega/Dt (???)
    end type PPICLF_U_t_interp

    type PPICLF_U_t_feedback
        !--- Particle Volume Fraction Feedback
        real*8                                  :: PHIP
        !--- x,y,z Forces Feedback
        ! x feedback force
        real*8                                  :: FX
        ! y feedback force
        real*8                                  :: FY
        ! z feedback force
        real*8                                  :: FZ
        !---Energy Feedback
        real*8                                  :: E
        !--- More VF quanities. ***NEED TO CONFIRM THEY ARE USED ***
        real*8                                  :: PHIPD
        real*8                                  :: PHIPU
        real*8                                  :: PHIPV
        real*8                                  :: PHIPW
        real*8                                  :: PHIPT
        !--- Reynolds Subgrid Stress Tensor (RSG)
        ! real*8                                  :: JRSG(3,3)
        real*8                                  :: RSG11
        real*8                                  :: RSG12
        real*8                                  :: RSG13
        real*8                                  :: RSG21
        real*8                                  :: RSG22
        real*8                                  :: RSG23
        real*8                                  :: RSG31
        real*8                                  :: RSG32
        real*8                                  :: RSG33
        !--- Pseudo Turbulent Kinetic Energy
        ! real*8                                  :: JTSG(3)
        real*8                                  :: TSG1
        real*8                                  :: TSG2
        real*8                                  :: TSG3
    end type PPICLF_U_t_feedback

    type PPICLF_t_iprop
        integer*4                               :: ParticleRank
        integer*4                               :: xBin
        integer*4                               :: yBin
        integer*4                               :: zBin
        integer*4                               :: binNum
        integer*4                               :: RemoveParticle
    end type PPICLF_t_iprop

#ifdef PPICLF_VU
    type PPICLF_t_DRUDTPLAG
        ! TODO: figure out what data is actually stored in here
        type(PPICLF_t_realNVec)                 :: unknown
    end type PPICLF_t_DRUDTPLAG
    type PPICLF_t_DRUDTMIXT
        ! TODO: figure out what data is actually stored in here
        type(PPICLF_t_realNVec)                 :: unknown
    end type PPICLF_t_DRUDTMIXT
#endif

    type PPICLF_U_t_particle
        type(PPICLF_t_tag)                      :: tag
        type(PPICLF_t_solutionProps)            :: Y
        type(PPICLF_t_solutionProps)            :: YDot
        type(PPICLF_t_solutionProps)            :: YDotc
        type(PPICLF_t_solutionProps)            :: Y1(PPICLF_N_PREV_SOLS)

        type(PPICLF_t_rprop)                    :: rprop
        type(PPICLF_t_iprop)                    :: iprop

! #ifdef PPICLF_VU
!         type(PPICLF_t_DRUDTPLAG)                :: DRUDTPLAG(PPICLF_VU)
!         type(PPICLF_t_DRUDTMIXT)                :: DRUDTMIXT(PPICLF_VU)
! #endif
        type(PPICLF_t_realNVec)                 :: drudtPlag(PPICLF_VU)
        type(PPICLF_t_realNVec)                 :: drudtMixt(PPICLF_VU)

        ! type(PPICLF_t_interp), pointer          :: interp
        ! type(PPICLF_t_feedback), pointer        :: feedback
        integer*8                               :: rngState
    end type PPICLF_U_t_particle

    type PPICLF_U_t_ghostParticle
        type(PPICLF_t_tag)                      :: tag
        type(PPICLF_t_solutionProps)            :: Y
        type(PPICLF_t_rprop)                    :: rprop
        type(PPICLF_t_iprop)                    :: iprop
    end type PPICLF_U_t_ghostParticle

    ! Operator Definitions for the solutions

    interface operator (+)
        module procedure solADDVec
        module procedure solADDScalar
        module procedure interpADDInterp
        module procedure interpADDScalar
        module procedure feedbackADDfeedback
    end interface
    interface operator(-)
        module procedure solSUBVec
        module procedure solSUBScalar
        module procedure interpSUBInterp
        module procedure interpSUBScalar
    end interface
    interface operator(*)
        module procedure solMULTScalar
        module procedure interpMULTScalar
        module procedure feedbackMULTscalar
    end interface
    interface operator(/)
        module procedure solDIVScalar
        module procedure interpDIVScalar
    end interface
    interface operator(.eq.)
        module procedure partEQpart
    end interface

    ! Handles to the MPI derived types
    ! The elements are the handles for the following types, in order:
    ! t_solutionProps, t_rprop, t_interp, t_feedback, t_iprop, t_drudtplag, t_drudtmixt, t_particle, t_ghostParticle
    integer, save :: user_MPI_TYPE_Handles(9)
    contains

    ! solution addition overloads
    function solADDVec(s1, s2) result(resSol)
        type(PPICLF_t_solutionProps), intent(in)    :: s1, s2
        type(PPICLF_t_solutionProps)                :: resSol
        resSol%pos      = s1%pos        + s2%pos
        resSol%vel      = s1%vel        + s2%vel
        resSol%T        = s1%T          + s2%T
        resSol%ang_vel  = s1%ang_vel    + s2%ang_vel
        resSol%metal    = s1%metal      + s2%metal
        resSol%oxide    = s1%oxide      + s2%oxide
    end function solADDVec
    function solADDScalar(s1, s2) result(resSol)
        type(PPICLF_t_solutionProps), intent(in)    :: s1
        real*8, intent(in)                          :: s2
        type(PPICLF_t_solutionProps)                :: resSol
        resSol%pos      = s1%pos        + s2
        resSol%vel      = s1%vel        + s2
        resSol%T        = s1%T          + s2
        resSol%ang_vel  = s1%ang_vel    + s2
        resSol%metal    = s1%metal      + s2
        resSol%oxide    = s1%oxide      + s2
    end function solADDScalar
    
    ! solution subtraction overloads
    function solSUBVec(s1, s2) result(resSol)
        type(PPICLF_t_solutionProps), intent(in)    :: s1, s2
        type(PPICLF_t_solutionProps)                :: resSol
        resSol%pos      = s1%pos        - s2%pos
        resSol%vel      = s1%vel        - s2%vel
        resSol%T        = s1%T          - s2%T
        resSol%ang_vel  = s1%ang_vel    - s2%ang_vel
        resSol%metal    = s1%metal      - s2%metal
        resSol%oxide    = s1%oxide      - s2%oxide
    end function solSUBVec
    function solSUBScalar(s1, s2) result(resSol)
        type(PPICLF_t_solutionProps), intent(in)    :: s1
        real*8, intent(in)                          :: s2
        type(PPICLF_t_solutionProps)                :: resSol
        resSol%pos      = s1%pos        - s2
        resSol%vel      = s1%vel        - s2
        resSol%T        = s1%T          - s2
        resSol%ang_vel  = s1%ang_vel    - s2
        resSol%metal    = s1%metal      - s2
        resSol%oxide    = s1%oxide      - s2
    end function solSUBScalar
    
    ! solution multiplication overloads
    function solMULTScalar(s1, s2) result(resSol)
        type(PPICLF_t_solutionProps), intent(in)    :: s1
        real*8, intent(in)                          :: s2
        type(PPICLF_t_solutionProps)                :: resSol
        resSol%pos      = s1%pos        * s2
        resSol%vel      = s1%vel        * s2
        resSol%T        = s1%T          * s2
        resSol%ang_vel  = s1%ang_vel    * s2
        resSol%metal    = s1%metal      * s2
        resSol%oxide    = s1%oxide      * s2
    end function solMULTScalar

    ! solution division overloads
    function solDIVScalar(s1, s2) result(resSol)
        type(PPICLF_t_solutionProps), intent(in)    :: s1
        real*8, intent(in)                          :: s2
        type(PPICLF_t_solutionProps)                :: resSol
        resSol%pos      = s1%pos        / s2
        resSol%vel      = s1%vel        / s2
        resSol%T        = s1%T          / s2
        resSol%ang_vel  = s1%ang_vel    / s2
        resSol%metal    = s1%metal      / s2
        resSol%oxide    = s1%oxide      / s2
    end function solDIVScalar

    function interpADDInterp(i1, i2)    result(resInterp)
        type(PPICLF_U_t_interp), intent(in) :: i1, i2
        type(PPICLF_U_t_interp)             :: resInterp
        resInterp%RHOf  = i1%RHOf   + i2%RHOf
        resInterp%PHIp  = i1%PHIp   + i2%PHIp
        resInterp%U     = i1%U      + i2%U
        resInterp%CS    = i1%CS     + i2%CS
        resInterp%DPDX  = i1%DPDX   + i2%DPDX
        resInterp%SDR   = i1%SDR    + i2%SDR
        resInterp%RHSR  = i1%RHSR   + i2%RHSR
        resInterp%PGC   = i1%PGC    + i2%PGC
        resInterp%Ft    = i1%Ft     + i2%Ft
        resInterp%VOR   = i1%VOR    + i2%VOR
        resInterp%P     = i1%P      + i2%P
        resInterp%RHOG  = i1%RHOG   + i2%RHOG
        resInterp%DPVDX = i1%DPVDX  + i2%DPVDX
        resInterp%SDO   = i1%SDO    + i2%SDO
    end function interpADDInterp
    function interpADDScalar(i1, i2)    result(resInterp)
        type(PPICLF_U_t_interp), intent(in) :: i1
        real*8, intent(in)                  :: i2
        type(PPICLF_U_t_interp)             :: resInterp
        resInterp%RHOf  = i1%RHOf   + i2
        resInterp%PHIp  = i1%PHIp   + i2
        resInterp%U     = i1%U      + i2
        resInterp%CS    = i1%CS     + i2
        resInterp%DPDX  = i1%DPDX   + i2
        resInterp%SDR   = i1%SDR    + i2
        resInterp%RHSR  = i1%RHSR   + i2
        resInterp%PGC   = i1%PGC    + i2
        resInterp%Ft    = i1%Ft     + i2
        resInterp%VOR   = i1%VOR    + i2
        resInterp%P     = i1%P      + i2
        resInterp%RHOG  = i1%RHOG   + i2
        resInterp%DPVDX = i1%DPVDX  + i2
        resInterp%SDO   = i1%SDO    + i2
    end function interpADDScalar
    function interpSUBInterp(i1, i2)    result(resInterp)
        type(PPICLF_U_t_interp), intent(in) :: i1, i2
        type(PPICLF_U_t_interp)             :: resInterp
        resInterp%RHOf  = i1%RHOf   - i2%RHOf
        resInterp%PHIp  = i1%PHIp   - i2%PHIp
        resInterp%U     = i1%U      - i2%U
        resInterp%CS    = i1%CS     - i2%CS
        resInterp%DPDX  = i1%DPDX   - i2%DPDX
        resInterp%SDR   = i1%SDR    - i2%SDR
        resInterp%RHSR  = i1%RHSR   - i2%RHSR
        resInterp%PGC   = i1%PGC    - i2%PGC
        resInterp%Ft    = i1%Ft     - i2%Ft
        resInterp%VOR   = i1%VOR    - i2%VOR
        resInterp%P     = i1%P      - i2%P
        resInterp%RHOG  = i1%RHOG   - i2%RHOG
        resInterp%DPVDX = i1%DPVDX  - i2%DPVDX
        resInterp%SDO   = i1%SDO    - i2%SDO
    end function interpSUBInterp
    function interpSUBScalar(i1, i2)    result(resInterp)
        type(PPICLF_U_t_interp), intent(in) :: i1
        real*8, intent(in)                  :: i2
        type(PPICLF_U_t_interp)             :: resInterp
        resInterp%RHOf  = i1%RHOf   - i2
        resInterp%PHIp  = i1%PHIp   - i2
        resInterp%U     = i1%U      - i2
        resInterp%CS    = i1%CS     - i2
        resInterp%DPDX  = i1%DPDX   - i2
        resInterp%SDR   = i1%SDR    - i2
        resInterp%RHSR  = i1%RHSR   - i2
        resInterp%PGC   = i1%PGC    - i2
        resInterp%Ft    = i1%Ft     - i2
        resInterp%VOR   = i1%VOR    - i2
        resInterp%P     = i1%P      - i2
        resInterp%RHOG  = i1%RHOG   - i2
        resInterp%DPVDX = i1%DPVDX  - i2
        resInterp%SDO   = i1%SDO    - i2
    end function interpSUBScalar
    function interpMULTScalar(i1, i2)   result(resInterp)
        type(PPICLF_U_t_interp), intent(in) :: i1
        real*8, intent(in)                  :: i2
        type(PPICLF_U_t_interp)             :: resInterp
        resInterp%RHOf  = i1%RHOf   * i2
        resInterp%PHIp  = i1%PHIp   * i2
        resInterp%U     = i1%U      * i2
        resInterp%CS    = i1%CS     * i2
        resInterp%DPDX  = i1%DPDX   * i2
        resInterp%SDR   = i1%SDR    * i2
        resInterp%RHSR  = i1%RHSR   * i2
        resInterp%PGC   = i1%PGC    * i2
        resInterp%Ft    = i1%Ft     * i2
        resInterp%VOR   = i1%VOR    * i2
        resInterp%P     = i1%P      * i2
        resInterp%RHOG  = i1%RHOG   * i2
        resInterp%DPVDX = i1%DPVDX  * i2
        resInterp%SDO   = i1%SDO    * i2
    end function interpMULTScalar
    function interpDIVScalar(i1, i2)    result(resInterp)
        type(PPICLF_U_t_interp), intent(in) :: i1
        real*8, intent(in)                  :: i2
        type(PPICLF_U_t_interp)             :: resInterp
        resInterp%RHOf  = i1%RHOf   / i2
        resInterp%PHIp  = i1%PHIp   / i2
        resInterp%U     = i1%U      / i2
        resInterp%CS    = i1%CS     / i2
        resInterp%DPDX  = i1%DPDX   / i2
        resInterp%SDR   = i1%SDR    / i2
        resInterp%RHSR  = i1%RHSR   / i2
        resInterp%PGC   = i1%PGC    / i2
        resInterp%Ft    = i1%Ft     / i2
        resInterp%VOR   = i1%VOR    / i2
        resInterp%P     = i1%P      / i2
        resInterp%RHOG  = i1%RHOG   / i2
        resInterp%DPVDX = i1%DPVDX  / i2
        resInterp%SDO   = i1%SDO    / i2
    end function interpDIVScalar

    ! feedback addition overloads
    pure function feedbackADDfeedback(f1, f2) result(resFDBK)
        type(PPICLF_U_t_feedback), intent(in) :: f1, f2
        type(PPICLF_U_t_feedback) :: resFDBK
        resFDBK%PHIP    = f1%PHIP   + f2%PHIP
        resFDBK%FX      = f1%FX     + f2%FX
        resFDBK%FY      = f1%FY     + f2%FY
        resFDBK%FZ      = f1%FZ     + f2%FZ
        resFDBK%E       = f1%E      + f2%E
        resFDBK%PHIPD   = f1%PHIPD  + f2%PHIPD
        resFDBK%PHIPU   = f1%PHIPU  + f2%PHIPU
        resFDBK%PHIPV   = f1%PHIPV  + f2%PHIPV
        resFDBK%PHIPW   = f1%PHIPW  + f2%PHIPW
        resFDBK%PHIPT   = f1%PHIPT  + f2%PHIPT
        resFDBK%RSG11   = f1%RSG11  + f2%RSG11
        resFDBK%RSG12   = f1%RSG12  + f2%RSG12
        resFDBK%RSG13   = f1%RSG13  + f2%RSG13
        resFDBK%RSG21   = f1%RSG21  + f2%RSG21
        resFDBK%RSG22   = f1%RSG22  + f2%RSG22
        resFDBK%RSG23   = f1%RSG23  + f2%RSG23
        resFDBK%RSG31   = f1%RSG31  + f2%RSG31
        resFDBK%RSG32   = f1%RSG32  + f2%RSG32
        resFDBK%RSG33   = f1%RSG33  + f2%RSG33
        resFDBK%TSG1    = f1%TSG1   + f2%TSG1
        resFDBK%TSG2    = f1%TSG2   + f2%TSG2
        resFDBK%TSG3    = f1%TSG3   + f2%TSG3
    end function feedbackADDfeedback

    ! feedback multiplication overloads
    pure function feedbackMULTscalar(f1, s) result(resFDBK)
        type(PPICLF_U_t_feedback), intent(in) :: f1
        real*8, intent(in) :: s
        type(PPICLF_U_t_feedback) :: resFDBK
        resFDBK%PHIP    = f1%PHIP   * s
        resFDBK%FX      = f1%FX     * s
        resFDBK%FY      = f1%FY     * s
        resFDBK%FZ      = f1%FZ     * s
        resFDBK%E       = f1%E      * s
        resFDBK%PHIPD   = f1%PHIPD  * s
        resFDBK%PHIPU   = f1%PHIPU  * s
        resFDBK%PHIPV   = f1%PHIPV  * s
        resFDBK%PHIPW   = f1%PHIPW  * s
        resFDBK%PHIPT   = f1%PHIPT  * s
        resFDBK%RSG11   = f1%RSG11  * s
        resFDBK%RSG12   = f1%RSG12  * s
        resFDBK%RSG13   = f1%RSG13  * s
        resFDBK%RSG21   = f1%RSG21  * s
        resFDBK%RSG22   = f1%RSG22  * s
        resFDBK%RSG23   = f1%RSG23  * s
        resFDBK%RSG31   = f1%RSG31  * s
        resFDBK%RSG32   = f1%RSG32  * s
        resFDBK%RSG33   = f1%RSG33  * s
        resFDBK%TSG1    = f1%TSG1   * s
        resFDBK%TSG2    = f1%TSG2   * s
        resFDBK%TSG3    = f1%TSG3   * s
    end function feedbackMULTscalar

    pure function partEQpart(p1, p2) result(res)
        type(PPICLF_U_t_particle), intent(in) :: p1, p2
        logical res
        res = .true.

        if(any(p1%y%pos%vec .ne. p2%y%pos%vec)) res = .false.
        if(any(p1%y%vel%vec .ne. p2%y%vel%vec)) res = .false.
        if(p1%y%oxide .ne. p2%y%oxide ) res = .false.
        if(p1%tag%partNum .ne. p2%tag%partNum ) res = .false.
        if(p1%tag%rankNum .ne. p2%tag%rankNum ) res = .false.
        if(p1%tag%cycleNum .ne. p2%tag%cycleNum ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
        ! if(p1% .ne. p2% ) res = .false.
    end function partEQpart
    ! checks for any nan values in an interp object
    pure function interpHasNAN(interp)  result(b)
        type(PPICLF_U_t_interp), intent(in) :: interp
        logical b
        b = &
            isnan(interp%RHOf)              .or. &
            isnan(interp%PHIp)              .or. &
            any(isnan(interp%U%vec))        .or. &
            isnan(interp%CS)                .or. &
            any(isnan(interp%DPDX%vec))     .or. &
            any(isnan(interp%SDR%vec))      .or. &
            isnan(interp%RHSR)              .or. &
            any(isnan(interp%PGC%vec))      .or. &
            isnan(interp%Ft)                .or. &
            any(isnan(interp%VOR%vec))      .or. &
            isnan(interp%P)                 .or. &
            any(isnan(interp%RHOG%vec))     .or. &
            any(isnan(interp%DPVDX%vec))    .or. &
            any(isnan(interp%SDO%vec))
    end function


    ! sets up the derived types in mpi, and returns the handles to the ones ppiclf needs (Particle, GhostParticle, Interp, feedback)
    subroutine ppiclf_user_Create_MPI_Derivedtypes(ParticleStuctHandle, GhostParticleStructHandle, InterpStructHandle, FeedbackStructHandle, realNVecHandle, tagHandle)
        !
        ! Input:
        !
        integer, intent(in) :: realNVecHandle, tagHandle
        integer, intent(out) :: ParticleStuctHandle, GhostParticleStructHandle, InterpStructHandle, FeedbackStructHandle
        !
        ! Internal:
        !
        integer :: extent_tag, extent_realNVec, extent_Real8, extent_Int4, extent_solutionProps, extent_rprop, extent_iprop, prev_extent
        integer :: i, ierr
        integer :: oldtypes(20), blockcounts(20), offsets(20)
        ! ! variables for faking a pointer type
        ! integer :: fakePointerType, test_extent, extent_pointer
        ! type(PPICLF_U_t_interp), pointer :: testPointer
        user_MPI_TYPE_Handles = -1
        ! t_solutionProps, t_rprop, t_interp, t_feedback, t_iprop, t_drudtplag, t_drudtmixt, t_particle, t_ghostParticle
        ! get the extent of realNVec and real8 mpi types, for use in calculating offsets for solutionprops
        call MPI_TYPE_EXTENT(tagHandle, extent_tag, ierr)
        call MPI_TYPE_EXTENT(realNVecHandle, extent_realNVec, ierr)
        call MPI_TYPE_EXTENT(MPI_REAL8, extent_Real8, ierr)
        call MPI_TYPE_EXTENT(MPI_INTEGER4, extent_Int4, ierr)

        i = 1
        ! create and commit MPI derived type for PPICLF_t_solutionProps
        ! 2 nVecs (pos, vel)
        oldtypes(i) = realNVecHandle
        blockcounts(i) = 2
        offsets(i) = 0
        i = i + 1

        ! 1 real8 (temp)
        oldtypes(i) = MPI_REAL8
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_realNVec
        i = i + 1

        ! 1 nvec (ang_vel)
        oldtypes(i) = realNVecHandle
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_Real8
        i = i + 1

        ! 2 real8 (metal, oxide)
        oldtypes(i) = MPI_REAL8
        blockcounts(i) = 2
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_realNVec

        call MPI_TYPE_STRUCT(i, blockcounts(1:i), offsets(1:i), oldtypes(1:i), user_MPI_TYPE_Handles(1), ierr)
        call MPI_TYPE_COMMIT(user_MPI_TYPE_Handles(1), ierr)

        ! get the extent of solutionProps, for use later
        call MPI_TYPE_EXTENT(user_MPI_TYPE_Handles(1), extent_solutionProps, ierr)

        i = 1
        ! create and commit MPI derived type for t_rprop
        ! 7 real8s (RHOp, Dp, VOLp, JDPi, JDPe, JSPL, JSPT)
        oldtypes(i) = MPI_REAL8
        blockcounts(i) = 7
        offsets(i) = 0
        i = i + 1

        ! 2 nvec (FLUCTF, WDOT)
        oldtypes(i) = realNVecHandle
        blockcounts(i) = 2
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_Real8
        i = i + 1

        ! 5 real8 (IDp, BRNT, XIPAR, XIPERP, XIT)
        oldtypes(i) = MPI_REAL8
        blockcounts(i) = 5
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_realNVec

        call MPI_TYPE_STRUCT(i, blockcounts(1:i), offsets(1:i), oldtypes(1:i), user_MPI_TYPE_Handles(2), ierr)
        call MPI_TYPE_COMMIT(user_MPI_TYPE_Handles(2), ierr)

        ! get the extent of rprop, for use later
        call MPI_TYPE_EXTENT(user_MPI_TYPE_Handles(2), extent_rprop, ierr)
        
        ! create and commit MPI derived type for t_interp
        i = 1
        ! 2 real8s (RHOf, PHIp)
        oldtypes(i) = MPI_REAL8
        blockcounts(i) = 2
        offsets(i) = 0
        i = i + 1

        ! 1 nvec (U)
        oldtypes(i) = realNVecHandle
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_Real8
        i = i + 1

        ! 1 real (CS)
        oldtypes(i) = MPI_REAL8
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_realNVec
        i = i + 1

        ! 2 nvecs (DPDX, SDR)
        oldtypes(i) = realNVecHandle
        blockcounts(i) = 2
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_Real8
        i = i + 1

        ! 1 real (RHSR)
        oldtypes(i) = MPI_REAL8
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_realNVec
        i = i + 1

        ! 1 nvec (PGC)
        oldtypes(i) = realNVecHandle
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_Real8
        i = i + 1

        ! 1 real (Ft)
        oldtypes(i) = MPI_REAL8
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_realNVec
        i = i + 1

        ! 1 nvec (VOR)
        oldtypes(i) = realNVecHandle
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_Real8
        i = i + 1

        ! 1 real (P)
        oldtypes(i) = MPI_REAL8
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_realNVec
        i = i + 1

        ! 3 nvec (RHOG, DPVDX, SDO)
        oldtypes(i) = realNVecHandle
        blockcounts(i) = 3
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_Real8

        call MPI_TYPE_STRUCT(i, blockcounts(1:i), offsets(1:i), oldtypes(1:i), user_MPI_TYPE_Handles(3), ierr)
        call MPI_TYPE_COMMIT(user_MPI_TYPE_Handles(3), ierr)
        InterpStructHandle = user_MPI_TYPE_Handles(3)

        ! create and commit MPI derived type for t_feedback
        ! need to verify this actually works
        call MPI_TYPE_CONTIGUOUS(22, MPI_REAL8, user_MPI_TYPE_Handles(4), ierr)
        call MPI_TYPE_COMMIT(user_MPI_TYPE_Handles(4), ierr)
        FeedbackStructHandle = user_MPI_TYPE_Handles(4)

        ! create and commit MPI derived type for t_iprop
        call MPI_TYPE_CONTIGUOUS(6, MPI_INTEGER4, user_MPI_TYPE_Handles(5), ierr)
        call MPI_TYPE_COMMIT(user_MPI_TYPE_Handles(5), ierr)

        ! get the extent of iprop
        call MPI_TYPE_EXTENT(user_MPI_TYPE_Handles(5), extent_iprop, ierr)

#ifdef PPICLF_VU
        ! create and commit MPI derived type for t_drudtplag
        call MPI_TYPE_CONTIGUOUS(1, realNVecHandle, user_MPI_TYPE_Handles(6), ierr)
        call MPI_TYPE_COMMIT(user_MPI_TYPE_Handles(6), ierr)
        ! create and commit MPI derived type for t_drudtmixt
        call MPI_TYPE_CONTIGUOUS(1, realNVecHandle, user_MPI_TYPE_Handles(7), ierr)
        call MPI_TYPE_COMMIT(user_MPI_TYPE_Handles(7), ierr)
#endif
        i = 1
        ! create and commit MPI derived type for t_particle
        oldtypes(i) = tagHandle
        blockcounts(i) = 1
        offsets(i) = 0
        i = i + 1

        ! 3 + PPICLF_N_PREV_SOLS solutionProps (Y, Ydot, Ydotc, Y1)
        oldtypes(i) = user_MPI_TYPE_Handles(1) ! handle for solution prop
        blockcounts(i) = 3 + PPICLF_N_PREV_SOLS
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_tag
        i = i + 1

       ! 1 rprop 
        oldtypes(i) = user_MPI_TYPE_Handles(2) ! handle for rprop
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_solutionProps
        i = i + 1

        ! 1 iprop
        oldtypes(i) = user_MPI_TYPE_Handles(5) ! handle for iprop
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_rprop
        i = i + 1
        prev_extent = extent_iprop
#ifdef PPICLF_VU
        ! PPICLF_VU t_drudtplag 
        oldtypes(i) = user_MPI_TYPE_Handles(6) ! handle for drudtplag
        blockcounts(i) = PPICLF_VU
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_iprop
        i = i + 1
        call MPI_TYPE_EXTENT(user_MPI_TYPE_Handles(6), prev_extent, ierr)

        ! PPICLF_VU t_drudtmixt
        oldtypes(i) = MPI_REAL8
        blockcounts(i) = PPICLF_VU
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * prev_extent
        i = i + 1
        call MPI_TYPE_EXTENT(user_MPI_TYPE_Handles(7), prev_extent, ierr)

        
#endif
        ! 1 int8 (rngState)
        oldtypes(i) = MPI_INTEGER8
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * prev_extent
        ! ! find a type with the same extent as a pointer
        ! fakePointerType = -1
        ! extent_pointer = STORAGE_SIZE(testPointer) / 8

        ! ! test integer8
        ! if (fakePointerType .eq. -1) then
        !     call MPI_TYPE_EXTENT(MPI_INTEGER8, test_extent, ierr)
        !     if (extent_pointer .eq. test_extent) then
        !         fakePointerType = MPI_INTEGER8
        !     end if
        ! endif


        ! if (fakePointerType .eq. -1) then
        !     call ppiclf_exittr("Failed to find a compatable type for a pointer", 0.0d0, 0)
        ! endif
        ! ! 2 pointers
        ! oldtypes(i) = fakePointerType
        ! blockcounts(i) = 2
        ! offsets(i) = offsets(i - 1) + blockcounts(i - 1) * prev_extent

        call MPI_TYPE_STRUCT(i, blockcounts(1:i), offsets(1:i), oldtypes(1:i), user_MPI_TYPE_Handles(8), ierr)
        call MPI_TYPE_COMMIT(user_MPI_TYPE_Handles(8), ierr)
        ParticleStuctHandle = user_MPI_TYPE_Handles(8)

        i = 1
        ! create and commit MPI derived type for t_ghostParticle
        oldtypes(i) = tagHandle
        blockcounts(i) = 1
        offsets(i) = 0
        i = i + 1

        ! 1 solutionProps (Y)
        oldtypes(i) = user_MPI_TYPE_Handles(1) ! handle for solution prop
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_tag
        i = i + 1

       ! 1 rprop 
        oldtypes(i) = user_MPI_TYPE_Handles(2) ! handle for rprop
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_solutionProps
        i = i + 1

        ! 1 iprop
        oldtypes(i) = user_MPI_TYPE_Handles(5) ! handle for iprop
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_rprop
        

        call MPI_TYPE_STRUCT(i, blockcounts(1:i), offsets(1:i), oldtypes(1:i), user_MPI_TYPE_Handles(9), ierr)
        call MPI_TYPE_COMMIT(user_MPI_TYPE_Handles(9), ierr)
        GhostParticleStructHandle = user_MPI_TYPE_Handles(9)
    end subroutine ppiclf_user_Create_MPI_Derivedtypes

    subroutine ppiclf_user_Destroy_MPI_Derivedtypes
        integer ierr, i
        do i = 1, 9
            if (user_MPI_TYPE_Handles(i) .ne. -1) then
                call MPI_TYPE_FREE(user_MPI_TYPE_Handles(i), ierr)
            end if
        end do
    end subroutine ppiclf_user_Destroy_MPI_Derivedtypes
end module PPICLF_user_particle



! Custom attributes
! these tell the preprocessor how the user file is going to use each of the properties. They are valid between substructures

! GHOST - 
! NOGHOST - This property does not need to be included in ghost particles. 

! SINGLETS - this property is only used within a single cycle. It doesn't need to be saved between timesteps, or when the particle is transfered to another rank
! MULTITS - This property is used across multiple timestep. 

! INTERNAL - This property is only used within the user files. We won't generate macros to access it from outside the user files, ie in the ROCPICL driver code.
! EXTERNAL - This property will be accessed outsdie the user macros. We will generate macros to allow for this access. Mutually exclusive with SINGLETS.

! INTERP - This property is an interpolated property. This implies NOGHOST, SINGLETS, and INTERNAL, but these can be overridden.
