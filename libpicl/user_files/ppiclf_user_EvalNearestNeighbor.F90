#include "PPICLF_STD.h"
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine to find nearest neighbor for particle-particle and
!     particle-wall interactions. Subroutine also includes
!     the new added-mass binary terms, developed by Sam Briney.
!
! Added: if collisional_flag = 1  F = Fn
!                            = 2  F = Fn + Ft + Tt
!                            = 3  F = Fn + Ft + Tt + Th + Tr
! where Tt = collisional torque
!       Th = hydrodynamic torque
!       Tr = rolling torque
!
! Note that the collisional and rolling torques are due to
!     particle-particle interactions and thus are evaulated
!     in ppiclf_user_EvalNearestNeighbor. Therefore, only the
!     hydrodynamic torque is left to be calculated.
!-----------------------------------------------------------------------
!
module ppiclf_m_user_EvalNearestNeighbor
    ! particle data
    use ppiclf_data, only: ppiclf_npart
    use ppiclf_m_particledata, only: 
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
    use ppiclf_data, only: ppiclf_npart_gp
    ! wall support variables
    use ppiclf_data, only:
    ! AngularPeriodic variables (?)(SEE NOTE IN ppiclf_data)
    use ppiclf_data, only:


    use ppiclf_m_user_data
    use ppiclf_m_user_RFLUdata
    use ppiclf_user_particle

    use ppiclf_user_AM_functions, only: resistance_pair
    use ppiclf_op, only: ppiclf_exittr
    implicit none
    contains

subroutine ppiclf_user_EvalNearestNeighbor(i,particle, interp, j, neighbor, fam, Wdot_neighbor_mean, R_pair, &
    upmean, u2pmean, icpmean, phipmean, nneighbors, rmass, rphip, rpi)
    !
    ! Input:
    !
    integer*4 i ! this particle's index
    type(PPICLF_U_t_particle) particle ! this particle, access using this instead of ppiclf_parts(i)
    type(PPICLF_U_t_interp), intent(in) :: interp ! interpolation data for this particle
    integer*4 j ! neighbor's index, positive for real(local rank) particles, negative for ghost particles. 0 means neighbor isn't a particle, but a boundary
    type(PPICLF_U_t_ghostParticle) :: neighbor ! neighbor data, in a ppiclf_U_t_ghostParticle struct, regardless of if it is actually a ghost particle.
    ! input for variables from YdotParticle
    type(PPICLF_t_realNVec), intent(inout) :: fam, Wdot_neighbor_mean, upmean, u2pmean
    real*8, intent(inout) :: R_pair(6,6), phipmean, rmass, rphip, rpi
    integer*4, intent(inout) :: icpmean, nneighbors

    !
    ! Internal:
    !
    real*8 rdist ! center to center distance between particle and neighbor
    ! real*8 rxdiff, rydiff, rzdiff
    type(PPICLF_t_realNVec) rdiff ! vector between particle centers. rdist is the magnitude of this vector
    real*8 pi2, rthresh, rm1, rm2,   &
        rmult, eta_n, rbot,  rdelta12,       &
        rv12_mag, rv12_mage, rksp_max, rnmag, rksp_wall, rextra,    &
        JDP_i,JDP_j
    real*8 eps
    ! real*8 rn_12x, rn_12y, rn_12z
    type(PPICLF_t_realNVec) rn_12
    ! Sam - from TLJ for box filter
    integer*4 ifilt
    real*8 adptfilter, dpl, phip, dist2
    ! real*8 xdist2, ydist2, zdist2
    type(PPICLF_t_realNVec) dist2_vec
    real*8 dist, rsig
    real*8 sig2, gkern, pi
    ! TODO: Fix this
    parameter(pi = 3.1415926535)

    ! 06/06/2024 - Thierry - Added Mass code
    integer*4 k, l, kk, ll
    real*8 alpha_local, rad
    ! real*8 rxdiff1, rydiff1, rzdiff1
    type(PPICLF_t_realNVec) rdiff1

    ! 07/16/2024 - TLJ added tangential component
    ! - does not take into account angular velocity
    ! real*8 unx, uny, unz
    type(PPICLF_t_realNVec) un
    ! real*8 utx, uty, utz
    type(PPICLF_t_realNVec) ut
    real*8 ut_mag, rn_mag
    ! real*8 rt_12x, rt_12y, rt_12z
    type(PPICLF_t_realNVec) rt_12
    real*8 Fn_mag, Ftmin
    real*8 eta_t, mu_c
    ! real*8 A12x, A12y, A12z 
    type(PPICLF_t_realNVec) A12
    real*8 rad1, rad2
    ! real*8 u12x, u12y, u12z
    type(PPICLF_t_realNVec) u12
    ! real*8 tcx, tcy, tcz
    type(PPICLF_t_realNVec) tc
    ! real*8 trx, try, trz
    type(PPICLF_t_realNVec) tr
    real*8 thetar, dp1, dp2, r12
    ! real*8 omgrx, omgry, omgrz
    type(PPICLF_t_realNVec) omgr
    real*8 omgr_mag
    ! real*8 Ftx, Fty, Ftz
    type(PPICLF_t_realNVec) Ft
    
    ! 04/03/2025 - TLJ added for spring stiffness coefficient
    real*8 nu1, nu2
    real*8 E1, E2, Estar
    real*8 r1, r2, Rstar 
    real*8 ksp1, ksp2, ksp_min

    ! 01/29/2025 - Thierry - added for particle collision with conical
    !                         wall domain
    real*8 rp, yp, zp, vp, wp, rp_new, yp_new, zp_new, vp_new, wp_new, rbound, urp, thetap

    !
    ! Code:
    !
    pi2  = rpi*rpi

    ! other particles
    if (j .ne. 0) then
        !Added spload and radius factor

        ! Compute mean particle diameter between i and j; delta_{ij}
        rthresh  = 0.5d0*(particle%rprop%DP + neighbor%rprop%DP)

        ! Compute vector components and distance between 
        !    centers of particles i and j; D_{ij}
        ! rxdiff = (@{USEPARTICLE(neighbor%y%pos%X)}@) - (@{USEPARTICLE(particle%y%pos%X)}@)
        ! rydiff = (@{USEPARTICLE(neighbor%y%pos%Y)}@) - (@{USEPARTICLE(particle%y%pos%Y)}@)
        ! rzdiff = (@{USEPARTICLE(neighbor%y%pos%Z)}@) - (@{USEPARTICLE(particle%y%pos%Z)}@)
        rdiff = neighbor%y%pos - particle%y%pos
        
        rdist = nvecMagnitude(rdiff)

        !-----------------------------------------------------------------------
        !
        ! For binary added-mass for Briney model

        ! 06/06/2024 - Thierry - Added Mass code continues here
        ! 07/09/2024 - TLJ - Updated
        ! 07/14/2024 - Thierry - Updated overlapping particles if statement

        ! Filter widths are set to be equal to 2*cell length in x,y,z
        ! directions (1:3)
        ! ppiclf_nndist is neighbor width - user defined
        !   = max(NEIGHBORWIDTH,4*Dp)
    
        if (am_flag == 2 .and. rdist <= ppiclf_nndist) then
            ! Do not overwrite rxdiff, rydiff, rzdiff
            ! rxdiff1 = rxdiff
            ! rydiff1 = rydiff
            ! rzdiff1 = rzdiff
            rdiff1 = rdiff

            ! Check if particles are overlapping, replace value
            !   if yes, since we get crazy resistance values
            if (rdist < rthresh) then
                ! rxdiff1 = rthresh
                ! rydiff1 = rthresh
                ! rzdiff1 = rthresh
                rdiff1%vec = rthresh
            endif

            ! Model only valid for local volume fraction
            ! less than 0.4, so we limit it here without
            ! over riding rphip
            ! limit alpha to mitigate misuse
            alpha_local = min(0.4, rphip) 
            
            ! Compute the resistance matrix
            ! Only valid for monodispersed particles
            rad = 0.5d0*(particle%rprop%DP)

            call resistance_pair(rdiff1, alpha_local, rad, R_pair)
                
            ! accumulate number of neighbors
            nneighbors = nneighbors + 1
                
            ! do k=1,3
            !     do l=1,3
            !         ! added mass
            !         Fam(k) = Fam(k) + R_pair(k,l)   * @{USEPARTICLE(particle%rprop%WDOT, skipIndex)}@(l)
            !         ! induced added mass
            !         Fam(k) = Fam(k) + R_pair(k,l+3) * @{USEPARTICLE(neighbor%rprop%WDOT, skipIndex)}@(l)
            !     end do ! l-loop


            !     ! accumulate neighbor acceleration
            !     Wdot_neighbor_mean(k) = Wdot_neighbor_mean(k)  + @{USEPARTICLE(neighbor%rprop%WDOT, skipIndex)}@(k)
            ! end do ! k-loop

            do l=1,3
                ! added mass
                Fam = Fam + (R_pair(1:3,l)   * particle%rprop%WDOT%vec(l))
                ! induced added mass
                Fam = Fam + (R_pair(1:3,l+3) * neighbor%rprop%WDOT%vec(l))
            end do ! l-loop
            
            Wdot_neighbor_mean = Wdot_neighbor_mean + neighbor%rprop%wdot
        end if ! am_flag==2 .and. rdist <= ppiclf_nndist

        !-----------------------------------------------------------------------
        !
        ! For particle-particle collision

        ! Cycle if rdist > rthresh
        eps = 0.0d0
        if (rdist .lt. rthresh+eps) then

            ! Compute spring stiffness constant dynamically.
            ! The number of collision timesteps (ksp) is set by the user
            ! k1 = k_{n,limit}
            ksp1 = rmass*rpi*rpi/((ksp*ppiclf_dt)**2)
            ! k2 = k_{hertzian}
            E1  = 1.0d9  ! Assumed value for Young's modulus
            E2  = 1.0d9  ! Assumed value for Young's modulus
            nu1 = 0.35d0 ! Assumed value for Poisson's ratio
            nu2 = 0.35d0 ! Assumed value for Poisson's ratio
            Estar = (1.0d0-nu1*nu1)/E1 + (1.0d0-nu2*nu2)/E2
            Estar = 1.0d0/Estar
            r1 = 0.5d0*(particle%rprop%DP)
            r2 = 0.5d0*(neighbor%rprop%DP)
            Rstar = r1*r2/(r1+r2)
            ksp2 = (4.0d0/3.0d0)*Estar*sqrt(Rstar)
            ksp2 = ksp2*sqrt(abs(rdist-rthresh))
            ! kn = min(k1,k2)
            ksp_min = min(ksp1,ksp2)

            rm1 = (particle%rprop%RHOP) * (particle%rprop%VOLP)
            rm2 = (neighbor%rprop%RHOP) * (neighbor%rprop%VOLP)
            
            rmult = (rm1*rm2)/(rm1+rm2)
            eta_n = -2.0d0*sqrt(ksp_min)*log(erest)/sqrt(log(erest)**2+pi2)*sqrt(rmult)

            !            print*,'COLLS: ',i,j,ksp1,ksp2,ksp_min,
            !     >              eta_n,rdist-rthresh,vmag

            ! Compute unit normal vector along line of contact 
            !   pointing from particle i to particle j
            rbot = 1.0d0/rdist
            ! rn_12x = rxdiff*rbot
            ! rn_12y = rydiff*rbot
            ! rn_12z = rzdiff*rbot
            rn_12 = rdiff * rbot
            rn_mag = rdist
         
            ! Relative velocity in normal direction
            ! u12x = (@{USEPARTICLE(particle%y%Vel%X)}@) - (@{USEPARTICLE(neighbor%y%Vel%X)}@)
            ! u12y = (@{USEPARTICLE(particle%y%Vel%Y)}@) - (@{USEPARTICLE(neighbor%y%Vel%Y)}@)
            ! u12z = (@{USEPARTICLE(particle%y%Vel%Z)}@) - (@{USEPARTICLE(neighbor%y%Vel%Z)}@)
            u12 = particle%y%vel - neighbor%y%vel

            if (collisional_flag>=2) then
               ! Add contribution from angular velocity
               rad1 = 0.5d0*(particle%rprop%DP)
               rad2 = 0.5d0*(neighbor%rprop%DP)
            !    A12x = rad1 * (@{USEPARTICLE(particle%y%ang_vel%X)}@) + rad2 * (@{USEPARTICLE(neighbor%y%ang_vel%X)}@)
            !    A12y = rad1 * (@{USEPARTICLE(particle%y%ang_vel%Y)}@) + rad2 * (@{USEPARTICLE(neighbor%y%ang_vel%Y)}@)
            !    A12z = rad1 * (@{USEPARTICLE(particle%y%ang_vel%Z)}@) + rad2 * (@{USEPARTICLE(neighbor%y%ang_vel%Z)}@)
               A12 = particle%y%ang_vel * rad1 + neighbor%y%ang_vel * rad2

            !    u12x = u12x + (A12y*rn_12z - A12z*rn_12y)
            !    u12y = u12y + (A12z*rn_12x - A12x*rn_12z)
            !    u12z = u12z + (A12x*rn_12y - A12y*rn_12x)
               u12%vec(1) = u12%vec(1) + (A12%vec(2) * rn_12%vec(3) - A12%vec(3) * rn_12%vec(2))
               u12%vec(2) = u12%vec(2) + (A12%vec(3) * rn_12%vec(1) - A12%vec(1) * rn_12%vec(3))
               u12%vec(3) = u12%vec(3) + (A12%vec(1) * rn_12%vec(2) - A12%vec(2) * rn_12%vec(1))
            endif

            ! Compute (u_ij \cdot n_ij)
            rv12_mag = nvecComponentSum(u12 * rn_12) ! u12x*rn_12x + u12y*rn_12y + u12z*rn_12z
         
            ! Compute delta_12 and normal parameters
            rdelta12 = rthresh - rdist
            rksp_max  = ksp_min*rdelta12
            rv12_mage = rv12_mag*eta_n
            rnmag     = -rksp_max - rv12_mage

            ! Normal collision force Fn = -rnmag*n_{ij}
            ! Scalar magnitude |Fn| = abs(rnmag)
            Fn_mag = abs(rnmag)

            ! Compute tangential unit vector
            ! unx = rv12_mag*rn_12x
            ! uny = rv12_mag*rn_12y
            ! unz = rv12_mag*rn_12z
            un = rn_12 * rv12_mag
            ! utx = u12x - unx
            ! uty = u12y - uny
            ! utz = u12z - unz
            ut = u12 - un
            ut_mag = nvecMagnitude(ut) ! sqrt(utx*utx + uty*uty + utz*utz)
            ut_mag = max(ut_mag,1.0d-8)
            ! rt_12x = utx/ut_mag
            ! rt_12y = uty/ut_mag
            ! rt_12z = utz/ut_mag
            rt_12 = ut / ut_mag

            ! Compute tangential collision force
            Ftmin = 0.0d0
            if (collisional_flag>=2) then ! Tangential component
                if (ut_mag > 0) then
                    mu_c  = 0.4d0  ! Dimensionless; Coulomb
                    eta_t = eta_n  ! Set to normal; damping
                    Ftmin  = -min(mu_c*Fn_mag,eta_t*ut_mag)  
                endif
            endif

            ! Compute contributions to angular velocities
            ! tcx = 0.0d0; tcy = 0.0d0; tcz = 0.0d0;
            tc%vec = 0
            ! trx = 0.0d0; try = 0.0d0; trz = 0.0d0;
            tr%vec = 0

            if (collisional_flag>=2) then

                ! Tangential force and Collision torque contributions
                !    Ftx = Ftmin*rt_12x
                !    Fty = Ftmin*rt_12y
                !    Ftz = Ftmin*rt_12z
                Ft = rt_12 * Ftmin
                rad1 = 0.5d0* (particle%rprop%DP) ! rpropi(PPICLF_R_JDP)
                ! tcx = rad1*(rn_12y*Ftz - rn_12z*Fty)
                ! tcy = rad1*(rn_12z*Ftx - rn_12x*Ftz)
                ! tcz = rad1*(rn_12x*Fty - rn_12y*Ftx)
                tc%vec(1) = rad1*(rn_12%vec(2)*Ft%vec(3) - rn_12%vec(3)*Ft%vec(2))
                tc%vec(2) = rad1*(rn_12%vec(3)*Ft%vec(1) - rn_12%vec(1)*Ft%vec(3))
                tc%vec(3) = rad1*(rn_12%vec(1)*Ft%vec(2) - rn_12%vec(2)*Ft%vec(1))

                if (collisional_flag>=3) then
                    ! Add Rolling torque contribution
                    thetar = 0.06  ! Needs to be calibrated
                    dp1 = particle%rprop%DP
                    dp2 = neighbor%rprop%DP
                    r12 = 0.5d0*(dp1*dp2)/(dp1+dp2)
                    ! omgrx = @{USEPARTICLE(particle%y%ang_vel%X)}@ - @{USEPARTICLE(neighbor%y%ang_vel%X)}@
                    ! omgry = @{USEPARTICLE(particle%y%ang_vel%Y)}@ - @{USEPARTICLE(neighbor%y%ang_vel%Y)}@
                    ! omgrz = @{USEPARTICLE(particle%y%ang_vel%Z)}@ - @{USEPARTICLE(neighbor%y%ang_vel%Z)}@
                    omgr = particle%y%ang_vel - neighbor%y%ang_vel
                    omgr_mag = nvecMagnitude(omgr) ! sqrt(omgrx*omgrx+omgry*omgry+omgrz*omgrz)
                    omgr_mag = max(omgr_mag,1.d-8)
                    ! trx = -thetar*Fn_mag*r12*omgrx/omgr_mag
                    ! try = -thetar*Fn_mag*r12*omgry/omgr_mag
                    ! trz = -thetar*Fn_mag*r12*omgrz/omgr_mag
                    tr = (omgr / omgr_mag) * (r12 * Fn_mag * (-thetar))
                endif
            endif


            ! Now update that part of the RHS of equations 
            !   that involve nearest neighbors

            ! Particle velocities
            ! @{USEPARTICLE(particle%ydotc%vel%X)}@ = @{USEPARTICLE(particle%ydotc%vel%X)}@ + rnmag*rn_12x + Ftmin*rt_12x
            ! @{USEPARTICLE(particle%ydotc%vel%Y)}@ = @{USEPARTICLE(particle%ydotc%vel%Y)}@ + rnmag*rn_12y + Ftmin*rt_12y
            ! @{USEPARTICLE(particle%ydotc%vel%Z)}@ = @{USEPARTICLE(particle%ydotc%vel%Z)}@ + rnmag*rn_12z + Ftmin*rt_12z
            particle%ydotc%vel = particle%ydotc%vel + (rn_12 * rnmag) + (rt_12 * Ftmin)

            ! Particle angular velocities
            ! @{USEPARTICLE(particle%ydotc%ang_vel%X)}@ = @{USEPARTICLE(particle%ydotc%ang_vel%X)}@ + tcx + trx
            ! @{USEPARTICLE(particle%ydotc%ang_vel%Y)}@ = @{USEPARTICLE(particle%ydotc%ang_vel%Y)}@ + tcy + try
            ! @{USEPARTICLE(particle%ydotc%ang_vel%Z)}@ = @{USEPARTICLE(particle%ydotc%ang_vel%Z)}@ + tcz + trz
            particle%ydotc%ang_vel = particle%ydotc%ang_vel + tc + tr

        end if ! rdist lt rthresh

        !-----------------------------------------------------------------------
        !
        ! Feedback fluctuation mean

        dist2 = MAX(ppiclf_filter(1),ppiclf_filter(2),ppiclf_filter(3))

        ! Box filter half-width dist2
        IF(qs_fluct_filter_adapt_flag.NE.0) THEN
            ! Adaptive filter defined wrt particle i
            ! Used for adaptive box or gaussian
            dpl = particle%rprop%DP
            phip = interp%PHIP
            adptfilter = ( 10.*(dpl**3)/max(1.e-4,phip) )**(1./3.)
            adptfilter = adptfilter/2.0
            IF(adptfilter .GT. dist2) dist2 = adptfilter
        END IF

        ! Check if particle lies inside box or gaussian filter
        ! xdist2 = abs((@{USEPARTICLE(particle%y%pos%X)}@) - (@{USEPARTICLE(neighbor%y%pos%X)}@))
        ! if (xdist2 .gt. dist2) return

        ! ydist2 = abs((@{USEPARTICLE(particle%y%pos%Y)}@) - (@{USEPARTICLE(neighbor%y%pos%Y)}@))
        ! if (ydist2 .gt. dist2) return

        ! if (ppiclf_ndim .eq. 3) then
        !     zdist2 = abs((@{USEPARTICLE(particle%y%pos%Z)}@) - (@{USEPARTICLE(neighbor%y%pos%Z)}@))
        !     if (zdist2 .gt. dist2) return
        ! endif
        dist2_vec = particle%y%pos - neighbor%y%pos
        if (any(abs(dist2_vec%vec) .gt. dist2)) return

        !
        ! The mean is calcuated according to Lattanzi etal,
        !   Physical Review Fluids, 2022.
        !
        if (j.ne.0) then ! this shouldn't be needed, as j is already compared to zero above
            if (qs_fluct_filter_flag==0) then
                ! upmean   = upmean + (@{USEPARTICLE(neighbor%y%vel%X)}@)
                ! vpmean   = vpmean + (@{USEPARTICLE(neighbor%y%vel%Y)}@)
                ! wpmean   = wpmean + (@{USEPARTICLE(neighbor%y%vel%Z)}@)
                upmean = upmean + neighbor%y%vel
                ! u2pmean  = u2pmean + (@{USEPARTICLE(neighbor%y%vel%X)}@)**2
                ! v2pmean  = v2pmean + (@{USEPARTICLE(neighbor%y%vel%Y)}@)**2
                ! w2pmean  = w2pmean + (@{USEPARTICLE(neighbor%y%vel%Z)}@)**2
                u2pmean = u2pmean + (neighbor%y%vel * neighbor%y%vel)
                icpmean  = icpmean + 1
            else if (qs_fluct_filter_flag==1) then
                ! See https://dpzwick.github.io/ppiclF-doc/algorithms/overlap_mesh.html
                dist = nvecMagnitude(dist2_vec) ! sqrt(xdist2**2 + ydist2**2 + zdist2**2)
                gkern = sqrt(pi*dist2**2/(4.0d0*log(2.0d0)))**(-ppiclf_ndim) * exp(-dist**2/(dist2**2/(4.0d0*log(2.0d0))))

                phipmean = phipmean + gkern*(neighbor%rprop%VOLP)
                ! upmean   = upmean + gkern * (@{USEPARTICLE(neighbor%y%vel%X)}@) * (@{USEPARTICLE(neighbor%rprop%VOLP)}@)
                ! vpmean   = vpmean + gkern * (@{USEPARTICLE(neighbor%y%vel%Y)}@) * (@{USEPARTICLE(neighbor%rprop%VOLP)}@)
                ! wpmean   = wpmean + gkern * (@{USEPARTICLE(neighbor%y%vel%Z)}@) * (@{USEPARTICLE(neighbor%rprop%VOLP)}@)
                upmean = upmean + ((neighbor%y%vel) * gkern * neighbor%rprop%VOLp)
                ! u2pmean  = u2pmean + gkern * ((@{USEPARTICLE(neighbor%y%vel%X)}@)**2) * (@{USEPARTICLE(neighbor%rprop%VOLP)}@)
                ! v2pmean  = v2pmean + gkern * ((@{USEPARTICLE(neighbor%y%vel%Y)}@)**2) * (@{USEPARTICLE(neighbor%rprop%VOLP)}@)
                ! w2pmean  = w2pmean + gkern * ((@{USEPARTICLE(neighbor%y%vel%Z)}@)**2) * (@{USEPARTICLE(neighbor%rprop%VOLP)}@)
                u2pmean = u2pmean + ((neighbor%y%vel * neighbor%y%vel) * gkern * neighbor%rprop%VOLp)
                icpmean = icpmean + 1
            end if
        end if


    !-----------------------------------------------------------------------
    !
    ! boundaries
    elseif (j .eq. 0) then

        !rksp_wall = ksp
        !rksp_wall = 1000

        ! give a bit larger collision threshold for walls
        rextra   = 0.05d0 !
        ! add sploading and radius factor 
        rthresh  = (0.5d0+rextra)* (particle%rprop%DP)
        
        ! rxdiff = @{USEPARTICLE(neighbor%y%pos%X)}@ - @{USEPARTICLE(particle%y%pos%X)}@
        ! rydiff = @{USEPARTICLE(neighbor%y%pos%Y)}@ - @{USEPARTICLE(particle%y%pos%Y)}@
        ! rzdiff = @{USEPARTICLE(neighbor%y%pos%Z)}@ - @{USEPARTICLE(particle%y%pos%Z)}@
        rdiff = neighbor%y%pos - particle%y%pos
        
        ! rdist = rxdiff**2 + rydiff**2 + rzdiff**2
        ! rdist = sqrt(rdist)
        rdist = nvecMagnitude(rdiff)
        
        if (rdist .gt. rthresh) return

        rm1 = (particle%rprop%RHOP)*(particle%rprop%VOLP)

        ! Compute spring stiffness constant dynamically, 
        !   which overrides the user defined value
        ! Need to make sure this formula is valid for a wall
        ! k1 = k_{n,limit}
        ksp1 = rm1*rpi*rpi/((ksp*ppiclf_dt)**2)
        ! k2 = k_{hertzian}
        E1  = 1.0d9  ! Assumed value for Young's modulus
        nu1 = 0.35d0 ! Assumed value for Poisson's ratio
        Estar = E1/(1.0d0-nu1*nu1)
        r1 = 0.5d0*particle%rprop%DP
        r2 = r1
        Rstar = r1*r2/(r1+r2)
        ksp2 = (2.0d0/3.0d0)*Estar*sqrt(Rstar)
        ksp2 = ksp2*sqrt(abs(rdist-rthresh))
        ! kn = min(k1,k2)
        rksp_wall = min(ksp1,ksp2)
        
        rmult = sqrt(rm1)
        eta_n = 2.0d0*sqrt(rksp_wall)*log(erest)/sqrt(log(erest)**2+pi2)*rmult
         
        rbot = 1.0d0/rdist
        ! rn_12x = rxdiff*rbot
        ! rn_12y = rydiff*rbot
        ! rn_12z = rzdiff*rbot
        rn_12 = rdiff * rbot
    
        rdelta12 = rthresh - rdist
        
        ! rv12_mag = - (@{USEPARTICLE(particle%y%vel%X)}@) * rn_12x- (@{USEPARTICLE(particle%y%vel%Y)}@) * rn_12y- (@{USEPARTICLE(particle%y%vel%Z)}@) * rn_12z
        rv12_mag = nvecComponentSum((particle%y%vel * rn_12) * (-1.0d0))

        rv12_mage = rv12_mag*eta_n
        rksp_max  = rksp_wall*rdelta12
        rnmag     = -rksp_max - rv12_mage

        ! not sure why this is commented out?
        ! @{USEPARTICLE(particle%ydotc%Vel%X)}@ = @{USEPARTICLE(particle%ydotc%Vel%X)}@ + rnmag*rn_12x
        ! @{USEPARTICLE(particle%ydotc%Vel%Y)}@ = @{USEPARTICLE(particle%ydotc%Vel%Y)}@ + rnmag*rn_12y
        ! @{USEPARTICLE(particle%ydotc%Vel%Z)}@ = @{USEPARTICLE(particle%ydotc%Vel%Z)}@ + rnmag*rn_12z
        ! particle%ydotc%vel = particle%ydotc%vel + (rn_12*rnmag)


        ! Particles leavind the domain with wall collisions
        ! Simple fix for a conical geometry
        yp = particle%y%pos%vec(2) !yi(PPICLF_JY)
        zp = particle%y%pos%vec(3) !yi(PPICLF_JZ)
        vp = particle%y%vel%vec(2) !yi(PPICLF_JVY)
        wp = particle%y%vel%vec(3) !yi(PPICLF_JVZ)

        ! rbound = sqrt(yj(PPICLF_JY)**2 + yj(PPICLF_JZ)**2)
        rbound = sqrt((neighbor%y%pos%vec(2))**2 + (neighbor%y%pos%vec(3))**2)
        rp = sqrt(yp**2 + zp**2)
        urp = sqrt(vp**2 + wp**2)
        thetap = atan2(zp, yp)

        if(rp > rbound) then
        
            rp_new = rp - (rp - rbound)
            yp_new = rp_new * cos(thetap)
            zp_new = rp_new * sin(thetap)

            urp = - urp
            
            vp_new = urp * cos(thetap)
            wp_new = urp * sin(thetap)

            ! ppiclf_y(PPICLF_JY,i) = yp_new
            ! ppiclf_y(PPICLF_JZ,i) = zp_new

            ! ppiclf_y(PPICLF_JVY,i) = vp_new
            ! ppiclf_y(PPICLF_JVZ,i) = wp_new
            particle%y%pos%vec(2) = yp_new
            particle%y%pos%vec(3) = zp_new

            particle%y%Vel%vec(2) = vp_new
            particle%y%Vel%vec(3) = wp_new
        endif
        
        !write(*,*) "Wall NEAR",i,ppiclf_ydotc(PPICLF_JVY,i)  
    endif


    return
end subroutine ppiclf_user_EvalNearestNeighbor

end module ppiclf_m_user_EvalNearestNeighbor
