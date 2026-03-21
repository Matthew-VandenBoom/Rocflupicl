#include "PPICLF_STD.h"
!-----------------------------------------------------------------------
!
! Sam - Code for solving ydot = F(y,t)
!
! Set internal flags
!   These flags allows the user to turn on or
!      off various force terms, or to select
!      between different versions of each force model
!   To turn off all particle motion, set stationary = 0
!   To turn off a force set the flag to zero
!   For two-way coupling set collisional_flag = 0
!   For four-way coupling set collisional_flag = 1
!   To turn off feedback force set feedback_flag = 0
!   To use fluctuations turn corresponding flag = 1
!
!
!   stationary = 0 if 1, do not move particles but do
!         calculate drag forces; feedback force can also
!         be turned on
!
!   qs_flag = 2  ! none = 0; Parmar = 1; Osnes = 2
!   am_flag = 1  ! Parmar = 1; Briney = 2
!   pg_flag = 1
!   collisional_flag = 1 ! two way coupled = 0 ; four-way = 1
!
!   heattransfer_flag = 1
!
!   ViscousUnsteady_flag = 0 no viscous unsteady drag
!                        = 1 history kernal for visc. unsteady drag
!
!   feedback_flag = 1
!
!   qs_fluct_flag = 1 ! None = 0 ; Lattanzi = 1 ; Osnes = 2 
!
!
!
!-----------------------------------------------------------------------
!
submodule (ppiclf_user) ppiclf_user_SetYdot_imp
    ! particle data
    use ppiclf_data, only: ppiclf_npart
    use ppiclf_m_particledata, only: ppiclf_parts, ppiclf_gparts
    ! grid data
    use ppiclf_data, only:
    use ppiclf_data, only:
    use ppiclf_data, only:
    ! particle options variables
    use ppiclf_data, only:
    use ppiclf_data, only: ppiclf_ndim
    use ppiclf_data, only: ppiclf_nndist, ppiclf_dt, ppiclf_time, ppiclf_rk3ark, ppiclf_filter
    ! use ppiclf_data, only:
    ! comm variables
    use ppiclf_data, only: ppiclf_nid
    ! binning variables
    use ppiclf_data, only: ppiclf_n_bins, ppiclf_bins_dx
    ! ghost particle variables
    use ppiclf_data, only: ppiclf_npart_gp, PPICLF_PPInteractions
    ! wall support variables
    use ppiclf_data, only:
    ! AngularPeriodic variables (?)(SEE NOTE IN ppiclf_data)
    use ppiclf_data, only:
#ifdef TEST
    use ppiclf_data, only: particle_nn, PPICLF_TOTNNDIST
#endif 

    use ppiclf_m_user_data
    use ppiclf_m_user_RFLUdata
    use ppiclf_m_types
    use ppiclf_m_wrapped_types, only: ppiclf_t_neighborInfo

    use mpi, only: MPI_WTIME
    use ppiclf_op, only: ppiclf_exittr
    use ppiclf_solve, only: ppiclf_solve_FindNearestNeighborSB_Start, ppiclf_solve_FindNearestNeighborSB

    use ppiclf_m_user_ForceModels, only: ppiclf_user_BR_driver, ppiclf_user_Lift_driver, ppiclf_user_Torque_driver, ppiclf_user_HT_driver
    use ppiclf_m_user_ForceModels, only: ppiclf_user_QS_Parmar, ppiclf_user_QS_ModifiedParmar, ppiclf_user_QS_Osnes, ppiclf_user_QS_Gidaspow
    use ppiclf_m_user_ForceModels, only: ppiclf_user_QS_fluct_Lattanzi, ppiclf_user_QS_fluct_Osnes
    use ppiclf_m_user_ForceModels, only: ppiclf_user_VU_Rocflu, ppiclf_user_UpdatePlag, ppiclf_user_ShiftUnsteadyData
    use ppiclf_m_user_ForceModels, only: ppiclf_user_AM_Briney_Unary, ppiclf_user_AM_Briney_Binary, ppiclf_user_AM_Parmar
    use ppiclf_m_user_SubbinMap, only: ppiclf_user_subbinMap
    use ppiclf_m_user_EvalNearestNeighbor, only: ppiclf_user_EvalNearestNeighbor

    use ppiclf_user_random
    implicit none
    contains

module procedure ppiclf_user_SetYdotInit
    integer*4 i
#ifdef TEST
    type(PPICLF_t_NNSB_Search_Data) searchInfo
    type(ppiclf_t_neighborInfo) searchResult
#endif
    if (allocated(SBin_map)) deallocate(SBin_map)
    if (allocated(SBin_counter)) deallocate(SBin_counter)

    allocate(SBin_map( 0 : ( (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_nndist)/ppiclf_nndist)         + 1) * (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_nndist)/ppiclf_nndist)        + 1) * (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_nndist)/ppiclf_nndist)        + 1) - 1), (ppiclf_npart+ppiclf_npart_gp)))
    allocate(SBin_counter( 0 : ( (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_nndist)/ppiclf_nndist)         + 1) * (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_nndist)/ppiclf_nndist)        + 1) * (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_nndist)/ppiclf_nndist)        + 1) - 1)))
    
    !-----------------------------------------------------------------------
    !
    ! Avery added 10/10/2024 - Map particles to subbins if collisional force, 
    ! Briney Added Mass force, QS fluctation force, or pseudo turbulence is flagged
    !
    if(PPICLF_PPInteractions) then
#ifdef PERF
        tstart = MPI_WTIME()
#endif
        call ppiclf_user_subbinMap(i_Bin, n_SBin, tot_SBin,SBin_counter ,SBin_map)
#ifdef PERF
        tfinal = MPI_WTIME()
        PPICLF_TParticleParticleModels = tfinal - tstart
#endif
    endif ! Collisions, QS Fluct, Briney AM, or pseudoTurb flags on

    icallb = icallb + 1
    nstage = 3
    istage = mod(icallb,nstage)
    if (istage .eq. 0) istage = 3  

    ! Count every iStage=1 for debug output
    if (iStage .eq. 1) idebug = idebug + 1

#ifndef TEST
    ! Print dt and time every time step
    if (ppiclf_nid==0) then
        if (istage .eq. 1) then
            write(6,'(a,2x,2(1pe14.6),2x,i3)') '*** PPICLF dt, time = ', ppiclf_dt,ppiclf_time
        endif
    endif
#endif

    burnrate_model = 0
    if (burnrate_flag .gt. 0) then
        if ( TRIM(ppiclf_matname)=='AL'.or. TRIM(ppiclf_matname)=='Al' ) then
            burnrate_model = 1
        elseif ( TRIM(ppiclf_matname)=='Mg' ) then
            burnrate_model = 2
        elseif ( TRIM(ppiclf_matname)=='C' ) then
            burnrate_model = 3
        else
            print*,'Error: no burn rate model'
            stop
        endif
    endif

    ! we generate a max of 12 random numbers per particle
    call ppiclf_user_random_setup(ppiclf_npart, 12)
    do i = 1, ppiclf_npart
        ppiclf_parts(i)%rngState = (i - 1) * 12
    end do

#ifdef TEST
    ! unit_test only tests the SB nearest neighbor search in this
    ! subroutine.  The full subroutine is called to ensure that
    ! the array initialization is correct.
    DO i = 1,ppiclf_npart
        PARTICLE_NN(i) = 0 
        PPICLF_TOTNNDIST(i) = 0.0D0
        call ppiclf_solve_FindNearestNeighborSB_Start(i,tot_SBin,SBin_counter,SBin_map,n_SBin,i_Bin, searchInfo)
        do
            searchResult = ppiclf_solve_FindNearestNeighborSB(searchInfo)
            if (.not. searchResult%exists) exit ! exit when the nearest neighbor search doesn't find another neighbor
            if (searchResult%j .eq. 0) cycle ! skip boundary
            ! if the search did find a neigbor save the resulting info for the unit tests to retrieve later
            PARTICLE_NN(i) = PARTICLE_NN(i) + 1
            ! PPICLF_TOTNNDIST(i) = PPICLF_TOTNNDIST(i) + dist_total
            PPICLF_TOTNNDIST(i) = PPICLF_TOTNNDIST(i) + nvecMagnitudeSQ((ppiclf_parts(i)%y%pos - searchResult%neighbor%y%pos))
        end do
    END DO
    RETURN
#endif

    ! Reapply axi-sym collision correction
    ! Right now hard coding smallest radius  
    do i=1,ppiclf_npart
       ppiclf_parts(i)%rprop%JDPe = (0.00005/(ppiclf_parts(i)%rprop%JSPT)) * (ppiclf_parts(i)%rprop%DP)
    end do 

    ! Set initial max values - must be done npart loop
    if (ppiclf_debug >= 1) then
        phimax    = 0.d0
        fqsx_max  = 0.d0; fqsy_max  = 0.d0; fqsz_max  = 0.d0
        famx_max  = 0.d0; famy_max  = 0.d0; famz_max  = 0.d0
        fdpdx_max = 0.d0; fdpdy_max = 0.d0; fdpdz_max = 0.d0
        fcx_max   = 0.d0; fcy_max   = 0.d0; fcz_max   = 0.d0
        fvux_max  = 0.d0; fvuy_max  = 0.d0; fvuz_max  = 0.d0
        qq_max    = 0.d0;
        fqsx_fluct_max = 0.d0; fqsy_fluct_max = 0.d0
        fqsz_fluct_max = 0.d0
        fqsx_total_max = 0.d0; fqsy_total_max = 0.d0
        fqsz_total_max = 0.d0
        fqs_mag = 0.0d0; fam_mag = 0.0d0; fdp_mag = 0.0d0
        fc_mag  = 0.0d0
        umean_max = 0.d0; vmean_max = 0.d0; wmean_max = 0.d0
    endif
end procedure ppiclf_user_SetYdotInit


module procedure ppiclf_user_YdotParticle
    !
    ! Input:
    !
    ! integer*4, intent(in) :: ip ! particle index into the internal ppiclf array, shouldn't be directly used by user. Given for use when calling other ppiclf functions
    ! type(PPICLF_U_t_particle), intent(inout) :: particle ! particle structure, as defined by the user in ppiclf_user_particle.F90
    ! type(PPICLF_U_t_interp), intent(in) :: interp ! Interpolated data from fluid solver
    ! type(PPICLF_U_t_feedback), intent(out) :: feedback ! feedback data, to be filled in ydotParticle
    ! integer*4, intent(out) :: ierr ! use to communicate errors back up the call stack from the force models to ppiclf_solve_setydot, which can actually call ppiclf_exittr

    !
    ! Internal:
    !
    ! real*8 fqsx, fqsy, fqsz
    type(PPICLF_t_realNVec) fqs
    real*8 fqsforce
    real*8 fqs_fluct(3)
    real*8 xi_par, xi_perp, xi_T
    ! real*8 famx, famy, famz ! see below 
    
    ! real*8 fdpdx, fdpdy, fdpdz
    type(PPICLF_t_realNVec) fdpd
    ! real*8 fdpvdx, fdpvdy, fdpvdz
    type(PPICLF_t_realNVec) fdpvd
    ! real*8 fcx, fcy, fcz
    type(PPICLF_t_realNVec) fc ! collision force (?)
    ! real*8 fbx, fby, fbz 
    type(PPICLF_t_realNVec) fb
    ! real*8 fvux, fvuy, fvuz
    type(PPICLF_t_realNVec) fvu

    ! real*8 ug, vg, wg
    type(PPICLF_t_realNVec) ug

    real*8 beta,cd

    real*8 factor, rcp_fluid, rmass_add

    real*8 gkern

    ! Needed for Added mass calculation
    integer*4 j, l
    real*8 SDrho, maxFilter

    real*8 vgradrhog
    integer*4 i, n, ic, k
    integer*4 store_forces

    ! Needed for heat transfer
    real*8 qq, rmass_therm, temp

    ! Needed for reactive particles
    real*8    mdot_me, mdot_ox
    real*8    Pres

    ! Needed for angular velocity
    ! real*8 taux, tauy, tauz
    type(PPICLF_t_realNVec) tau
    real*8 tau_mag

    real*8 rmass_omega
    
    ! real*8 taux_hydro, tauy_hydro, tauz_hydro 
    type(PPICLF_t_realNVec) tau_hydro

    ! real*8 liftx, lifty, liftz
    type(PPICLF_t_realNVec) lift
    real*8 lift_mag

    !
    ! For ppiclf_user_Fluctuations.f
    !
    ! COMMON Name: user_flcut01, user_fluct02, user_fluct03
    integer*4 icpmean
    real*8 phipmean
    ! real*8 upmean, vpmean, wpmean
    ! real*8 u2pmean, v2pmean, w2pmean
    type(PPICLF_t_realNVec) upmean, u2pmean

    real*8 UnifRnd(6), Rsg(3,3), T_par(3)

    !
    ! For ppiclf_user_AddedMass.f
    !
    ! COMMON Name: user_AddedMass01, user_AddedMass02
    integer*4 nneighbors
    ! real*8 Fam(3) ! replaced the redundant(?) Fam(3) and famx, famy, famz, with one realNvec
    type(PPICLF_t_realNVec) fam, Famunary, FamBinary
    ! real*8 FamUnary(3), FamBinary(3)

    ! real*8 Wdot_neighbor_mean(3)
    type(PPICLF_t_realNVec) Wdot_neighbor_mean
    real*8 R_pair(6,6)


    real*8 rmu,rkappa,rmass,rhof,dp,rep,rphip,                  &
        rphif,asndf,rmachp,rhop,rhoMixt,reyL,rnu,fac,               &
        rcp_part,rpr,                                               &
        phi, mp, re, rem



    real*8 vmag !, vx, vy, vz
    type(PPICLF_t_realNVec) v

    ! internal data for the nearest neighbor SB search function
    ! stored here so that it can be reused between calls to the search function
    ! without causing issues when we try to parallelize it
    type(PPICLF_t_NNSB_Search_Data) searchInfo
    type(ppiclf_t_neighborInfo) searchResult ! the result from one call to ppiclf_solve_findnearestneighbor

#ifdef PERF
    REAL*8 tstart, tfinal
#endif

    !
    ! Code:
    !


    ! rpi        = acos(-1.0d0)
    rcp_part   = ppiclf_rcp_part
    rpr        = 0.70d0
    rcp_fluid  = 1004.64d0

    fac = ppiclf_rk3ark(iStage)*ppiclf_dt
    if (1==2) then
        if (ppiclf_nid==0) print*,'dt,fac=',istage,ppiclf_dt,fac,stationary, qs_flag, am_flag, pg_flag,collisional_flag, heattransfer_flag, feedback_flag,qs_fluct_flag, ppiclf_debug, ppiclf_nTimeBH,ppiclf_nUnsteadyData
    endif

    !
    ! Evaluate ydot, the rhs of the equations of motion
    ! for the particles
    !

    ! do i=1,ppiclf_npart

    ! Choose viscosity law
    if (rmu_flag==rmu_fixed_param) then
        ! Constant viscosity law
        rmu = rmu_ref
    elseif (rmu_flag==rmu_suth_param) then
        ! Sutherland law
        temp    = interp%Ft
        rmu     = rmu_ref*sqrt(temp/tref)*(1.0d0+suth/tref)/(1.0d0+suth/temp)
    else
        call ppiclf_exittr('Unknown viscosity law$', 0.0d0, 0)
        return ! dummy return, but lets the compiler know that rmu is always initialized pass this point
    endif
    rkappa = rcp_fluid*rmu/rpr


    ! Useful values
    rmass  = particle%rprop%VOLP * particle%rprop%RHOP
    ! vx     = interp%U%X - particle%y%vel%x
    ! vy     = interp%U%Y - particle%y%vel%y
    ! vz     = interp%U%Z - particle%y%vel%z
    v = interp%U - particle%y%vel

    vmag   = nvecMagnitude(v)
    rhof   = interp%RHOf  
    dp     = particle%rprop%Dp
    rep    = vmag*dp*rhof/rmu
    rphip  = interp%PHIp
    rphif  = 1.0d0-(interp%PHIp)
    asndf  = interp%CS
    rmachp = vmag/asndf
    rhop   = particle%rprop%RHOP

    ! TLJ - 04/03/2025; Do not calculate forces if vmag = 0
    !       Otherwise the particles might move before the 
    !       shock arrives
    if (vmag <= 1.d-3) return ! cycle

    ! Thierry - at initial times when rmachp and vmag are very
    ! small, we would get very large CD (>100)
    ! This leads to NaN projected force when PseudoTurbulence is 
    ! enabled. One fix I'm trying is to cycle when rmachp and vmag
    ! are small

    ! Thierry - seeing if that fixes the issue of very large CD
    ! initially 
    if(vmag .lt. 1.0 .or. rmachp .lt. 1.d-3)  return !cycle

    ! TLJ - redefined rprop(PPICLF_R_JSPT,i) to be the particle
    !   velocity magnitude for plotting purposes - 01/03/2025

    ! ***This is bad practice and leads to difficult code to debug -Avery***
    particle%rprop%JSPT = nvecMagnitude(particle%y%vel)

    rep = max(0.1d0,rep)

    ! Redefine volume fractions
    ! Need to make sure phi_p + phi_f = 1
    rphip = interp%PHIp
    rphif = 1.0d0-rphip

    ! TLJ: Needed for viscous unsteady force
    !      Using same nomenclature as rocinteract subroutines
    reyL = dp*vmag*rhof/rmu
    rnu = rmu/rhof

    ! Zero out for each particle i
    ! famx = 0.0d0; famy = 0.0d0; famz = 0.0d0; 
    rmass_add = 0.0d0;
    Fam%vec = 0.0d0
    FamUnary%vec =  0.0d0;
    FamBinary%vec =0.0d0;
    Wdot_neighbor_mean%vec = 0.0d0
    nneighbors = 0.0d0
    ! fqsx = 0.0d0; fqsy = 0.0d0; fqsz = 0.0d0; beta = 0.0d0;
    fqs%vec = 0
    fqs_fluct = 0.0d0;
    ! fdpdx = 0.0d0; fdpdy = 0.0d0; fdpdz = 0.0d0;
    fdpd%vec = 0
    ! fcx = 0.0d0; fcy = 0.0d0; fcz = 0.0d0;
    fc%vec = 0
    ! taux = 0.0d0; tauy = 0.0d0; tauz = 0.0d0;
    tau%vec = 0
    ! liftx = 0.0d0; lifty = 0.0d0; liftz = 0.0d0;
    lift%vec = 0
    ! fvux = 0.0d0; fvuy = 0.0d0; fvuz = 0.0d0;
    fvu%vec = 0
    qq=0.0d0
    mdot_me = 0.0d0; mdot_ox = 0.0d0;
    ! upmean = 0.0; vpmean = 0.0; wpmean = 0.0;
    ! u2pmean = 0.0; v2pmean = 0.0; w2pmean = 0.0;
    upmean%vec = 0
    u2pmean%vec = 0
    ! fdpvdx = 0.0d0; fdpvdy = 0.0d0; fdpvdz = 0.0d0;
    fdpd%vec = 0
    !--- Added for PseudoTurbulence
    Rsg = 0.0d0; T_par = 0.0d0



    !
    ! Step 1a: New Added-Mass model of Briney
    !
    ! 07/15/2024 - If am_flag = 2, then we need to call
    !   the Unary term before we make any calls to nearest
    !   neighbor
    ! 06/05/2024 - Thierry - For each particle i, initialize
    ! variables to be used in nearest neighbors to zero
    ! before looping over particle j (j neq i)
    ! Briney Added Mass flag
#ifdef PERF
    tstart = MPI_WTIME()
#endif
    if (am_flag == 2) then 
        ! 07/14/24 - Thierry - If Briney Algorithm flag and fluct_flag
        !   are ON -> evaluate added-mass unary term before evaluating
        !   neighbor-induced acceleration in EvalNearestNeighbor
        call ppiclf_user_AM_Briney_Unary(particle, interp, fam, FamUnary, rmass_add, rmachp, rphif, rhof, v)
    endif ! end am_flag = 2


    !
    ! Step 1b: Call NearestNeighbor if particles i and j interact
    !
    if(PPICLF_PPInteractions) THEN

        !AVERY - we should fix vmag ~0 bug and remove conditional check
        if ((qs_fluct_flag>=1) .and. (vmag .gt. 1.d-8)) then
            ! Compute mean for particle i
            !    add neighbor particle j afterward
            ! Box filter is used if qs_fluct_filter_flag=0
            !   The box filter used here is a simple cube centered
            !     at particle i with half-width dist2 (see
            !     ppiclf_user_EvalNearestNeighbor.f for definition)
            !   We use a simple arithmetic mean
            !   phipmean is not used
            ! Gaussian filter is used if qs_fluct_filter_flag=1
            !   We use the value of the Gaussian times the volume
            !     of the particle to get the filtered particle volume

            if (qs_fluct_filter_flag==0) then
                ! box filter
                !***AVERY - I think phipmean is using wrong index ***
                ! Also, none of the below are means...
                phipmean = particle%rprop%VOLP
                ! upmean   = particle%y%vel%vec(1)
                ! vpmean   = particle%y%vel%vec(2)
                ! wpmean   = particle%y%vel%vec(3)
                upmean = particle%y%vel
                ! u2pmean  = upmean**2
                ! v2pmean  = vpmean**2
                ! w2pmean  = wpmean**2
                u2pmean = upmean * upmean
                icpmean  = 1
            else if (qs_fluct_filter_flag==1) then
                ! gaussian kernel
                ! r = 0
                !***AVERY - verify max filter dimension should be used***
                maxFilter = MAX(ppiclf_filter(1),ppiclf_filter(2), ppiclf_filter(3))
                gkern = sqrt(rpi*maxFilter**2/ (4.0d0*log(2.0d0)))**(-ppiclf_ndim)
                !***AVERY - I think phipmean is using wrong index ***
                ! Also, none of the below are means...
                phipmean = gkern * particle%rprop%VOLp
                ! upmean   = gkern * (particle%y%vel%vec(1)) * particle%rprop%VOLP
                ! vpmean   = gkern * (particle%y%vel%vec(2)) * particle%rprop%VOLP
                ! wpmean   = gkern * (particle%y%vel%vec(3)) * particle%rprop%VOLP
                upmean = particle%y%vel * gkern * particle%rprop%volp
                ! u2pmean  = gkern * ((particle%y%vel%vec(1))**2)* particle%rprop%VOLP
                ! v2pmean  = gkern * ((particle%y%vel%vec(2))**2)* particle%rprop%VOLP
                ! w2pmean  = gkern * ((particle%y%vel%vec(3))**2)* particle%rprop%VOLP
                u2pmean = (particle%y%vel * particle%y%vel) * gkern * particle%rprop%volp
                icpmean = 1
            end if
        end if

        ! add neighbors
        ! AVERY - always use subbins
        call ppiclf_solve_FindNearestNeighborSB_Start(i,tot_SBin,SBin_counter,SBin_map,n_SBin,i_Bin, searchInfo)
        do
            searchResult = ppiclf_solve_FindNearestNeighborSB(searchInfo)
            if (.not. searchResult%exists) exit ! exit when the nearest neighbor search doesn't find another neighbor
            
            ! if the search did find a neigbor, run EvalNearestNeighbor on it
            call ppiclf_user_EvalNearestNeighbor(ip, particle, interp, searchResult%j, searchResult%neighbor, &
                fam, Wdot_neighbor_mean, R_pair, upmean, u2pmean, icpmean, phipmean, &
                nneighbors, rmass, rphip, rpi)
        end do
    end if ! end Step 1b; nearestneighbor
#ifdef PERF
    tfinal = MPI_WTIME()
    PPICLF_TParticleParticleModels = PPICLF_TParticleParticleModels + (tfinal - tstart)
    ! TODO: Should this be tstart
    tfinal = MPI_WTIME()
#endif 

    !
    ! Step 2: Force component quasi-steady
    !
    if (qs_flag==1) then 
        call ppiclf_user_QS_Parmar(mp, phi, re, cd, beta, rmachp, rphip, rep, dp, rmu)
    else if (qs_flag==2) then 
        call ppiclf_user_QS_Osnes(mp, phi, re, cd, beta, rmachp, rphip, rep, dp, rmu)
    else if (qs_flag==3) then 
        call ppiclf_user_QS_ModifiedParmar(mp, phi, re, cd, beta, rmachp, rphip, rep, dp, rmu)
    else if (qs_flag==4) then 
        call ppiclf_user_QS_Gidaspow(phi, re, cd, beta, rphip, rep, dp, rmu, rphif, rhof, vmag)
    else
        print*, "***PPICLF: Error in QS Model Selection!"   
        call ppiclf_exittr('Wrong QS Model Choice$', 0.0d0, qs_flag)
    endif

    ! fqsx = beta*vx
    ! fqsy = beta*vy
    ! fqsz = beta*vz
    fqs = v * beta

    !
    ! Step 3: Force fluctuation for quasi-steady force
    !
    ! Note: QS fluctuations needs nearest neighbors,
    !   and is called above in Step 1b
    if (qs_fluct_flag==1) then
        call ppiclf_user_QS_fluct_Lattanzi(particle, interp, fqs_fluct, UnifRnd, upmean, u2pmean, icpmean, fac, rphip, vmag, rmachp, rmu, dp, rep, phipmean)
    elseif (qs_fluct_flag==2 .or. pseudoTurb_flag==1) then
        call ppiclf_user_QS_fluct_Osnes(particle, interp, fqs_fluct, UnifRnd, rsg, t_par, mp, phi, re, rem, upmean, u2pmean, &
                icpmean, rphip, vmag, rmachp, rmu, dp, rep, phipmean, fac, fqs, v, xi_par, xi_perp, xi_T)
    endif

    ! Add fluctuation part to quasi-steady force
    ! fqsx = fqsx + fqs_fluct(1)
    ! fqsy = fqsy + fqs_fluct(2)
    ! fqsz = fqsz + fqs_fluct(3)
    fqs = fqs + fqs_fluct

    ! Store quasi-steady fluctuating force
    ! @{USEPARTICLE(particle%rprop%FLUCTF%X)}@ = fqs_fluct(1)
    ! @{USEPARTICLE(particle%rprop%FLUCTF%Y)}@ = fqs_fluct(2)
    ! @{USEPARTICLE(particle%rprop%FLUCTF%Z)}@ = fqs_fluct(3)
    particle%rprop%FLUCTF%vec = fqs_fluct
    
    ! Store normally distributed random variables xi for PseudoTurbulence
    particle%rprop%XIPAR  = xi_par
    particle%rprop%XIPERP = xi_perp
    particle%rprop%XIT    = xi_T


    !
    ! Step 4: Force component added mass
    !
    if (am_flag == 1) then 
        call ppiclf_user_AM_Parmar(particle, interp, fam, rmass_add, rmachp, rphip, rhof, rphif, v)

        !-----------------------------------------------------------------------
        !Thierry - Added Mass code continues here
        
    elseif (am_flag == 2) then 

        ! Thierry - binary_model.f90 evaluates the terms
        ! in the folllowing order:
        !   (1) Unary Term
        !   (2) Evaluates Neighbor Acceleration
        !   (3) Binary Term
        ! We replicate that here by calling them in the same order
        ! Unary and Binary calculations are now under 
        !   two separate subroutines
        ! Thierry - need to make sure NearestNeighbor is called
        !    if fluct_flag = 0 (ie, no QS fluctuations)

        ! Binary subroutine only valid when number of neighbors .gt. 0
        if (nneighbors .gt. 0) then
            call ppiclf_user_AM_Briney_Binary(particle, interp, fam, rphip, Wdot_neighbor_mean, nneighbors)
            ! FamBinary(1) = famx - FamUnary(1)
            ! FamBinary(2) = famy - FamUnary(2)
            ! FamBinary(3) = famz - FamUnary(3)
            FamBinary = Fam - FamUnary
        else
            ! if particle has no neighbors, need to multiply added mass forces
            ! by volume, as this is taken care of in Binary subroutine
            ! famx = famx* particle%rprop%VOLP 
            ! famy = famy* particle%rprop%VOLP 
            ! famz = famz* particle%rprop%VOLP
            Fam = Fam * particle%rprop%VOLp
        endif
    endif


    !-----------------------------------------------------------------------

    !
    ! Step 5: Force component pressure gradient
    !
    if (pg_flag == 1) then
        ! fdpdx = -(particle%rprop%VOLP) * (interp%DPDX%vec(1))
        ! fdpdy = -(particle%rprop%VOLP) * (interp%DPDX%vec(2))
        ! fdpdz = -(particle%rprop%VOLP) * (interp%DPDX%vec(3))
        fdpd = interp%DPDX * (- particle%rprop%VOLp)

        if (flow_model == 1) then ! Navier-Stokes Flow Model
            ! fdpvdx = (particle%rprop%VOLP) * (interp%DPVDX%vec(1))
            ! fdpvdy = (particle%rprop%VOLP) * (interp%DPVDX%vec(2))
            ! fdpvdz = (particle%rprop%VOLP) * (interp%DPVDX%vec(3))
            fdpvd = interp%DPVDX * particle%rprop%VOLp
        endif ! flow_model

        ! fdpdx = fdpdx + fdpvdx
        ! fdpdy = fdpdy + fdpvdy
        ! fdpdz = fdpdz + fdpvdz
        fdpd = fdpd + fdpvd
    endif ! end pg_flag = 1

    !
    ! Step 6: Force component collisional force, ie, particle-particle
    !
    if (collisional_flag >= 1) then
        ! Collision force:
        !  A discrete numerical model for granular assemblies
        !  - Cundall and Strack (1979)
        !  - Geotechnique

        ! Sam - STILL NEED TO VALIDATE COLLISION FORCE
        
        ! fcx  = particle%ydotc%vel%vec(1)
        ! fcy  = particle%ydotc%vel%vec(2)
        ! fcz  = particle%ydotc%vel%vec(3)
        fc = particle%ydotc%vel


    endif ! collisional_flag >= 1


    !
    ! Step 7: Viscous unsteady force with history kernel
    !
    if (ViscousUnsteady_flag==1) then
        call ppiclf_user_VU_Rocflu(particle, interp,fvu, reyL, rnu, dp, rphip, vmag)
    endif


    !
    ! Step 8a: Combustion model for reactive particles
    !
    rmass_therm = rmass*rcp_part
    qq = 0.0d0

    if (burnrate_flag >= 1) then
        call ppiclf_user_BR_driver(particle,interp, iStage,burnrate_model,qq,mdot_me,mdot_ox, rpi, vmag, ierr)
        if (ierr .ne. 0) return
    endif

    !
    ! Step 8b: Heat transfer model
    !
    if (heattransfer_flag >= 1) then
        call ppiclf_user_HT_driver(particle, interp, qq, rkappa, dp, rep, rpr, rmachp, rphif, ierr)
    endif ! heattransfer_flag >= 1


    !
    ! Step 9a: Angular velocity model
    !
    rmass_omega = rmass*dp*dp/10.0d0

    if (collisional_flag >= 2) then
        ! taux  = particle%ydotc%ang_vel%vec(1)
        ! tauy  = particle%ydotc%ang_vel%vec(2)
        ! tauz  = particle%ydotc%ang_vel%vec(3)
        tau = particle%ydotc%ang_vel
        call ppiclf_user_Torque_driver(particle, interp, iStage, tau, tau_hydro, dp, rhop, rhof, rmass, rmu, rpi, ierr)
    endif ! collisional_flag >= 2

    !
    ! Step 9b: Saffman and Magnus Lift models
    !          Lift models requires gas-phase vorticity and
    !          particle angular velocity
    !
    if (collisional_flag == 4) then
        call ppiclf_user_Lift_driver(particle, interp, iStage,lift, v, vmag, dp, rep, rhof, rmu, rnu, ierr)
    endif ! collisional_flag == 4


    !
    ! Step 10: Set ydot for all PPICLF_SLN number of equations
    !
    ! particle%ydot%pos%vec(1)        = particle%y%vel%vec(1)
    ! particle%ydot%pos%vec(2)        = particle%y%vel%vec(2)
    ! particle%ydot%pos%vec(3)        = particle%y%vel%vec(3)
    particle%ydot%pos               = particle%y%vel
    ! particle%ydot%vel%vec(1)        = (fqsx+famx+fdpdx+fvux+liftx+fcx)/(rmass+rmass_add)
    ! particle%ydot%vel%vec(2)        = (fqsy+famy+fdpdy+fvuy+lifty+fcy)/(rmass+rmass_add)
    ! particle%ydot%vel%vec(3)        = (fqsz+famz+fdpdz+fvuz+liftz+fcz)/(rmass+rmass_add)
    particle%ydot%vel                = (fqs + fam + fdpd + fvu + lift + fc) / (rmass + rmass_add)
    particle%ydot%t                 = qq/rmass_therm
    ! particle%ydot%ang_vel%vec(1)    = taux/rmass_omega
    ! particle%ydot%ang_vel%vec(2)    = tauy/rmass_omega
    ! particle%ydot%ang_vel%vec(3)    = tauz/rmass_omega
    particle%ydot%ang_vel           = tau/rmass_omega
    particle%ydot%metal             = mdot_me
    particle%ydot%oxide             = mdot_ox

    !
    ! Update data for viscous unsteady case
    !
    if (ViscousUnsteady_flag>=1) then
        call ppiclf_user_UpdatePlag(particle, interp, v, rphif, rhof)
    endif

    !
    ! Step 11: Feed Back force to the gas phase
    !
    !    Note that Rocflu uses a negative of the RHS, and
    !    so ppiclf must respect this odd convention.
    !
    ! Project work done by hydrodynamic forces:
    !   Inter-phase heat transfer and energy coupling in turbulent 
    !   dispersed multiphase flows
    !   - Ling et al. (2016)
    !   - Physics of Fluids
    ! See also for more details
    !   Explosive dispersal of particles in high speed environments
    !   - Durant et al. (2022)
    !   - Journal of Applied Physics


    IF(feedback_flag==0) THEN
        feedback%FX = 0.0d0 
        feedback%FY = 0.0d0 
        feedback%FZ = 0.0d0 
        feedback%E  = 0.0d0
    END IF

    IF(feedback_flag==1) THEN
        ! Momentum equations feedback terms
        feedback%FX = (particle%rprop%JSPL) * ((particle%ydot%vel%vec(1))*rmass - fc%vec(1))
        feedback%FY = (particle%rprop%JSPL) * ((particle%ydot%vel%vec(2))*rmass - fc%vec(2))
        feedback%FZ = (particle%rprop%JSPL) * ((particle%ydot%vel%vec(3))*rmass - fc%vec(3))

        ! Energy equation feedback term
        ! 09/19/2025 - Thierry - Added Lift force
        ! Still need to add Torue \cdot angular velocity
        ! feedback%JE = particle%rprop%JSPL           &
        !     * ( (fqsx+fvux+liftx)*(@{USEPARTICLE(particle%y%vel%X)}@)    +               &
        !         (fqsy+fvuy+lifty)*(@{USEPARTICLE(particle%y%vel%Y)}@)    +               &
        !         (fqsz+fvuz+liftz)*(@{USEPARTICLE(particle%y%vel%Z)}@)    +               &
        !         famx*(@{USEPARTICLE(particle%rprop%U%X)}@)               +               &
        !         famy*(@{USEPARTICLE(particle%rprop%U%Y)}@)               +               &
        !         famz*(@{USEPARTICLE(particle%rprop%U%Z)}@)               +               &
        !         taux_hydro*(@{USEPARTICLE(particle%y%ang_vel%X)}@)       +               &
        !         tauy_hydro*(@{USEPARTICLE(particle%y%ang_vel%Y)}@)       +               &
        !         tauz_hydro*(@{USEPARTICLE(particle%y%ang_vel%Z)}@)       +               &
        !         qq )
        feedback%E = particle%rprop%JSPL * (                           &
            nvecComponentSum((fqs + fvu + lift) * particle%y%vel)   +   &
            nvecComponentSum(fam * interp%U)                        +   &
            nvecComponentSum(tau_hydro*particle%y%ang_vel)          +   &
            qq                                                          &
        )

        IF(pseudoTurb_flag==1) THEN
            ! 09/02/2025 -  Addition of PTKE to Rocflu's Energy Equation

            ! (fqsx+fvux+famx+liftx) * (@{USEPARTICLE(particle%y%Vel%X)}@)     +       &
            ! (fqsy+fvuy+famy+lifty) * (@{USEPARTICLE(particle%y%Vel%Y)}@)     +       &
            ! (fqsz+fvuz+famz+liftz) * (@{USEPARTICLE(particle%y%Vel%Z)}@)     +       &

            feedback%E = particle%rprop%JSPL * (nvecComponentSum((fqs + fvu + fam + lift) * particle%y%vel) + qq )
        ELSE
            Rsg   = 0.0D0
            T_par = 0.0D0
        END IF ! pseudoTurb_flag
        ! 07/21/2025 - Thierry - Added Reynolds Subgrid Stress Feedback
        feedback%RSG11 = Rsg(1,1) * (particle%rprop%JSPL)
        feedback%RSG12 = Rsg(1,2) * (particle%rprop%JSPL)
        feedback%RSG13 = Rsg(1,3) * (particle%rprop%JSPL)
        feedback%RSG21 = Rsg(2,1) * (particle%rprop%JSPL)
        feedback%RSG22 = Rsg(2,2) * (particle%rprop%JSPL)
        feedback%RSG23 = Rsg(2,3) * (particle%rprop%JSPL)
        feedback%RSG31 = Rsg(3,1) * (particle%rprop%JSPL)
        feedback%RSG32 = Rsg(3,2) * (particle%rprop%JSPL)
        feedback%RSG33 = Rsg(3,3) * (particle%rprop%JSPL)

        feedback%TSG1 = T_par(1) * particle%rprop%JSPL
        feedback%TSG2 = T_par(2) * particle%rprop%JSPL
        feedback%TSG3 = T_par(3) * particle%rprop%JSPL

    END IF ! Feedback flag

    ! Update volume fraction feedback quantities with feedback on or off
    feedback%PHIP    = (particle%rprop%VOLP) * (particle%rprop%JSPL)
    feedback%PHIPD   = (particle%rprop%VOLP) * (particle%rprop%RHOP)
    feedback%PHIPU   = (particle%rprop%VOLP) * (particle%y%Vel%vec(1))
    feedback%PHIPV   = (particle%rprop%VOLP) * (particle%y%Vel%vec(2))
    feedback%PHIPW   = (particle%rprop%VOLP) * (particle%y%Vel%vec(3))
    feedback%PHIPT   = (particle%rprop%VOLP) * (particle%y%T)


    ! Step 12: If stationary, don't move particles. Feedback can still be on
    ! though.
    !
    if (stationary .gt. 0) then
        if (stationary==1) then
            particle%ydot%pos%vec(1)         = 0.0d0
            particle%ydot%pos%vec(2)         = 0.0d0
            particle%ydot%pos%vec(3)         = 0.0d0
            particle%ydot%vel%vec(1)         = 0.0d0
            particle%ydot%vel%vec(2)         = 0.0d0
            particle%ydot%vel%vec(3)         = 0.0d0
            particle%ydot%T                  = 0.0d0
            particle%ydot%ang_vel%vec(1)     = 0.0d0
            particle%ydot%ang_vel%vec(2)     = 0.0d0
            particle%ydot%ang_vel%vec(3)     = 0.0d0
        else
            call ppiclf_exittr('Unknown stationary flag$', 0.0d0, 0)
        endif
    elseif(stationary .lt. 0) then
        call ppiclf_user_unit_tests(i,iStage,fam)
        ! particle%ydot%pos%vec(1)     = particle%y%vel%vec(1)
        ! particle%ydot%pos%vec(2)     = particle%y%vel%vec(2)
        ! particle%ydot%pos%vec(3)     = particle%y%vel%vec(3)
        particle%ydot%pos            = particle%y%vel
        ! particle%ydot%vel%vec(1)     = particle%ydot%vel%vec(1) + fam%vec(1)
        ! particle%ydot%vel%vec(2)     = particle%ydot%vel%vec(2) + fam%vec(2)
        ! particle%ydot%vel%vec(3)     = particle%ydot%vel%vec(3) + fam%vec(3)
        particle%ydot%vel            = particle%ydot%vel + fam
        particle%ydot%T              = 0.0d0
        particle%ydot%ang_vel%vec(1) = 0.0d0
        particle%ydot%ang_vel%vec(2) = 0.0d0
        particle%ydot%ang_vel%vec(3) = 0.0d0
    endif



    !
    ! Step 13: Store forces (DELETED)

    !
    ! Step 14: If debug mode is ON, calculate and print the max values.
    !          The user should not have this ON for production runs.
    ! This must be changed before GPU runs
    !
    if (ppiclf_debug .ge. 1) then
        phimax = max(phimax,abs(rphip))

        fqsx_max = max(fqsx_max,abs(fqs%vec(1)))
        fqsy_max = max(fqsy_max,abs(fqs%vec(2)))
        fqsz_max = max(fqsz_max,abs(fqs%vec(3)))
        fqs_mag  = max(fqs_mag,nvecMagnitude(fqs)) ! sqrt(fqsx*fqsx+fqsy*fqsy+fqsz*fqsz))

        fqsx_fluct_max = max(fqsx_fluct_max, abs(fqs_fluct(1)))
        fqsy_fluct_max = max(fqsy_fluct_max, abs(fqs_fluct(2)))
        fqsz_fluct_max = max(fqsz_fluct_max, abs(fqs_fluct(3)))

        ! Why are we doing this, it the exact same as above?
        fqsx_total_max = max(fqsx_total_max, abs(fqs%vec(1)))
        fqsy_total_max = max(fqsy_total_max, abs(fqs%vec(2)))
        fqsz_total_max = max(fqsz_total_max, abs(fqs%vec(3)))

        umean_max = max(umean_max, abs(upmean%vec(1)))
        vmean_max = max(vmean_max, abs(upmean%vec(2)))
        wmean_max = max(wmean_max, abs(upmean%vec(3)))

        famx_max = max(famx_max,abs(fam%vec(1)))
        famy_max = max(famy_max,abs(fam%vec(2)))
        famz_max = max(famz_max,abs(fam%vec(3)))
        fam_mag  = max(fam_mag,nvecMagnitude(fam)) ! sqrt(famx*famx+famy*famy+famz*famz))

        fdpdx_max = max(fdpdx_max,abs(fdpd%vec(1)))
        fdpdy_max = max(fdpdy_max,abs(fdpd%vec(2)))
        fdpdz_max = max(fdpdz_max,abs(fdpd%vec(3)))
        fdp_mag   = max(fdp_mag,nvecMagnitude(fdpd)) ! sqrt(fdpdx*fdpdx+fdpdy*fdpdy+fdpdz*fdpdz))

        fcx_max = max(fcx_max, abs(fc%vec(1)))
        fcy_max = max(fcy_max, abs(fc%vec(2)))
        fcz_max = max(fcz_max, abs(fc%vec(3)))
        fc_mag  = max(fc_mag,nvecMagnitude(fc)) ! sqrt(fcx*fcx+fcy*fcy+fcz*fcz))

        fvux_max = max(fvux_max, abs(fvu%vec(1)))
        fvuy_max = max(fvuy_max, abs(fvu%vec(2)))
        fvuz_max = max(fvuz_max, abs(fvu%vec(3)))

        qq_max = max(qq_max, abs(qq))

        tau_mag = nvecMagnitude(tau) ! sqrt(taux*taux + tauy*tauy + tauz*tauz)
        tau_max = max(tau_max, abs(tau_mag))

        lift_mag = nvecMagnitude(lift) ! sqrt(liftx**2 + lifty**2 + liftz**2)
        lift_max = max(lift_max,lift_mag)

        if (ppiclf_debug.eq.2 .and. ppiclf_nid.eq.0) then
            if (iStage==3) then
                if (i==1) then
                    write(7010,*) i,ppiclf_time,rmass,vmag,rhof,dp,rep,rphip,rphif,rmachp,rhop,rhoMixt,reyL,rmu,rnu,rkappa
                endif
                if (i==ppiclf_npart) then
                    write(7011,*) i,ppiclf_time,rmass,vmag,rhof,dp,rep,rphip,rphif,rmachp,rhop,rhoMixt,reyL,rmu,rnu,rkappa
                endif
            endif
        endif

    endif ! ppiclf_debug .ge. 1
        
    ! write out for debug
    if (ppiclf_debug==3) then
        if (ppiclf_nid==0 .and. iStage==1) then
            if (mod(idebug,1)==0) then
                if (i<=5) then
                    write(7020+i,*) i, ppiclf_time, rhof,   &
                        (interp%SDR),        &
                        (particle%ydot%vel),     &
                        (particle%y%Vel),        &
                        (particle%y%ang_vel)

                    write(7040+i,*) i, ppiclf_time,     &
                        interp%SDR,    &   ! Du/Dt
                        interp%SDO        ! DOmega/Dt

                    write(7050+i,*) i, ppiclf_time, fqs_mag,fam_mag,fdp_mag,fc_mag,tau_max

                    write(7060+i,*) i, ppiclf_time, fc, lift, tau !fcx,fcy,fcz,liftx,lifty,liftz,taux,tauy,tauz
                endif
            endif
        endif
    endif



    ! enddo ! do i=1,ppiclf_npart

    !
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    !06/05/2024 - Thierry - Store density-weighted acceleration
    !
    ! Briney Added Mass flag
    if (am_flag==2) then 
        ! do i=1,ppiclf_npart
     
        ! Substantial derivative of density - how rocflu does it  
        SDrho = nvecComponentSum(particle%y%vel * interp%RHOG) + interp%RHSR

        
        ! material derivative is phi weighted in Rocflu
        ! drho/dt
        SDrho = SDrho / (rphif)  
        ! vgradrhog = vx * (interpRHOG%vec(1)) +  &
        !             vy * (interpRHOG%vec(2)) +  &
        !             vz * (interpRHOG%vec(3))
        vgradrhog = nvecComponentSum(v * interp%RHOg)
    
        ! Fluid density
        rhof   = interp%RHOF
        ! this is just repeating work already done above
        ! vx = (@{USEPARTICLE(particle%rprop%U%X)}@) - (@{USEPARTICLE(particle%y%Vel%X)}@)
        ! vy = (@{USEPARTICLE(particle%rprop%U%Y)}@) - (@{USEPARTICLE(particle%y%Vel%Y)}@)
        ! vz = (@{USEPARTICLE(particle%rprop%U%Z)}@) - (@{USEPARTICLE(particle%y%Vel%Z)}@)
        
        ! ug = (@{USEPARTICLE(particle%rprop%U%X)}@)
        ! vg = (@{USEPARTICLE(particle%rprop%U%Y)}@)
        ! wg = (@{USEPARTICLE(particle%rprop%U%Z)}@)
        ug = interp%U
        ! Unary added mass solves rho^g d(u^p)/dt implicitly
        ! Binary added mass solves it explicitly and not implicitly
        ! WDOTX = D(rho^g u^g)/Dt - d(rho^g u^p)/dt)
        ! X-acceleration
        ! particle%rprop%WDOT%X = vx*SDrho + rhof*(particle%rprop%SDR%X) + ug*vgradrhog - rhof*(particle%ydot%Vel%X)
        
        ! ! Y-acceleration
        ! particle%rprop%WDOT%Y = vy*SDrho + rhof*(particle%rprop%SDR%Y) + vg*vgradrhog - rhof*(particle%ydot%Vel%Y)

        ! ! Z-acceleration
        ! particle%rprop%WDOT%Z = vz*SDrho + rhof*(particle%rprop%SDR%Z) + wg*vgradrhog - rhof*(particle%ydot%Vel%Z)
        particle%rprop%WDOT = v * SDrho + interp%SDR * rhof + ug * vgradrhog - particle%ydot%vel * rhof

        ! write out for debug
        if (ppiclf_debug==2) then
            if (ppiclf_nid==0 .and. iStage==1) then
                if (mod(idebug,10)==0) then
                    if (i<=3) then
                        write(7020+i,*) i, ppiclf_time, rhof,   &
                            interp%SDR,                         &
                            particle%ydot%vel,                  &
                            particle%y%vel

                        write(7030+i,*) i, ppiclf_time, particle%rprop%WDOT

                    endif
                endif
            endif
        endif

        ! enddo
    endif

    
#ifdef PERF
    tfinal = MPI_WTIME()
    PPICLF_TParticleParticleModels = PPICLF_TParticleParticleModels + (tfinal - tstart)
#endif
    ! ----------------------------------------------------------------------

    return
end procedure ppiclf_user_YdotParticle

module procedure ppiclf_user_SetYdotFinal
    !
    ! ----------------------------------------------------------------------
    !

    ! Use ppiclf ALLREDUCE to compute values across processors
    ! Note that ALLREDUCE uses MPI_BARRIER, which is cpu expensive
    ! Print out every 10th iStage=1 counts
    if (ppiclf_debug   .ge. 1) then
        if (iStage         .eq. 1) then
            if (mod(idebug,10) .eq. 0) then
                call ppiclf_user_debug
            endif
        endif
    endif

    !
    ! ----------------------------------------------------------------------
    !
    !

    !
    ! Reset arrays for Viscous Unsteady Force
    !
    if (ViscousUnsteady_flag>=1) then
        if (iStage==3) call ppiclf_user_ShiftUnsteadyData
        ! call ppiclf_user_plag2prop
    endif
end procedure ppiclf_user_SetYdotFinal
    
subroutine ppiclf_user_unit_tests(i,iStage,fam)
    integer*4 i, iStage
    type(PPICLF_t_realNVec) fam
end subroutine ppiclf_user_unit_tests

subroutine ppiclf_user_debug
end subroutine ppiclf_user_debug

end submodule ppiclf_user_SetYdot_imp