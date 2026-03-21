#include "PPICLF_STD.h"

PROGRAM main
    use mpi
    use UT_m_data
    use UT_m_setup
    use ppiclf_data, only: x_per_flag, x_per_min, x_per_max, y_per_flag, y_per_min, y_per_max, z_per_flag, z_per_min, z_per_max , ang_per_flag, ang_per_angle, ang_per_xangle, ang_per_rin, ang_per_rout
    use ppiclf_data, only: ppiclf_np, ppiclf_nndist, ppiclf_overlap, ppiclf_filter, ppiclf_n_bins, ppiclf_npart, ppiclf_binb, ppiclf_linperiodic
    use ppiclf_m_particledata, only: ppiclf_parts
    use ppiclf_data, only: ppiclf_equaldomain, PPICLF_TOTNNDIST, particle_nn

    use ppiclf_solve, only: ppiclf_solve_Initialize, ppiclf_solve_InitParticle, ppiclf_solve_Interpolate, ppiclf_solve_InterpFieldUser, ppiclf_solve_GetProFld, ppiclf_solve_InitSolve
    use ppiclf_m_comm, only: ppiclf_comm_CreateBin, ppiclf_comm_InitMPI, ppiclf_comm_InitOverlapGrid, ppiclf_comm_FindParticle, ppiclf_comm_CreateGhost,ppiclf_comm_MoveGhost, ppiclf_comm_MoveParticle
    use ppiclf_user, only: ppiclf_user_SetYdotInit
    use ppiclf_m_types
    IMPLICIT NONE

    INTEGER*4 i, j, k, l, ie

    ! Grid variables

    REAL*8    tpF(PPICLF_LEE), numBins, feedback1(PPICLF_LEE), feedback2(PPICLF_LEE), x_norm, y_norm, z_norm

    ! Particle variables
    REAL*8    xp, yp, zp, totErr, InteriorErr, xFaceErr, &
               yFaceErr, zFaceErr, xyEdgeErr, xzEdgeErr, yzEdgeErr, &
               xyzCornerErr, totCnt, InteriorCnt, &
               xFaceCnt, yFaceCnt, zFaceCnt, xyEdgeCnt, xzEdgeCnt, &
               yzEdgeCnt, xyzCornerCnt, weightCell(27), &
               weightTot, xcell, ycell, zcell, T_calc, &
               calcErr, eps, T_analytic, NNDistSQ

    INTEGER*4 ip, np(3), part_cell(3), xstart, xend, &
            ystart, yend, zstart, zend, icount, loopcount,  &
            NNCount, id1, id2, projCells

    LOGICAL   xFace, yFace, zFace, farAway, interpolation(8),  &
            projection(8), binGen, nnpart(8), interp_logical(8), &
            proj_logical(8), nn_logical(8)

    ! Projection variables
    REAL*8    wsum, dSQl, dSQi, dist, CellVol, GaussianConst, &
            w(PPICLF_LEE),part_feedbk1(PPICLF_LPART), &
            TrueFeedback1(PPICLF_LEE), TrueFeedback2(PPICLF_LEE), &
            part_feedbk2(PPICLF_LPART), error1, error2, e1avg, e2avg

    CHARACTER*50 filename, testcase, procString, npstring, par

    rootProc = 0
    DO i = 1,8
        nnpart(i)        = .TRUE.
    END DO
    gridBoundsDefined  = .FALSE.
    ! MPI Setup
    !**********************************************************************
    CALL MPI_Init(ierr)
    icomm = MPI_COMM_WORLD
    CALL MPI_Comm_size(icomm,nproc,ierr)
    CALL MPI_Comm_rank(icomm,nid,ierr)
    CALL ppiclf_comm_InitMPI(icomm, nid, nproc)
    WRITE(npstring, '(I0)') nproc
    WRITE(procString, '(I0)') nid
    IF(nid .LT. 10) THEN
        procString = '0' // TRIM(procString) //'_' // TRIM(npstring) 
    ELSE
        procString = TRIM(procString) // '_' // TRIM(npstring)
    END IF

    CALL UT_setup

    ! Particle to Particle Nearest Neighbor Test
    !**********************************************************************
    ! Start loop for varying periodicity cases
    DO test = 1,8
        ! Periodicity flag Setup
        CALL test_setperiodic(x_per_flag,y_per_flag,z_per_flag, test,testcase)
        ! Will handle angular periodicity separately
        ang_per_flag   = 0
        ang_per_angle  = 0.0D0
        ang_per_xangle = 0.0D0
        ang_per_rin    = 0.0D0
        ang_per_rout   = 0.0D0  
        CALL MPI_BARRIER(icomm,ierr)
        ! ppiclF Calls
        !**********************************************************************
        PPICLF_OVERLAP = .FALSE.
        CALL ppiclf_solve_InitParticle(2,3,0,npart_local, p_parts,filter,nndist)
        CALL ppiclf_solve_Initialize(.true.,                           &
                                    x_per_flag, x_per_min, x_per_max,  &
                                    y_per_flag, y_per_min, y_per_max,  &
                                    z_per_flag, z_per_min, z_per_max,  &
                                    ang_per_flag, ang_per_angle,       &
                                    ang_per_xangle, ang_per_rin,       &
                                    ang_per_rout)
        CALL ppiclf_comm_InitOverlapGrid(proc_ncells,p_grid)
        CALL ppiclf_comm_CreateBin
        CALL ppiclf_comm_FindParticle
        CALL ppiclf_comm_MoveParticle
        CALL ppiclf_comm_CreateGhost
        CALL ppiclf_comm_MoveGhost
        CALL MPI_BARRIER(icomm,ierr)
        IF(ppiclf_npart .GT. 0) THEN 
            ! SetYdot only runs with nearest neighbor subroutine
            ! due to preprocess compile flag.  PARTICLE_NN(i) saves
            ! the number of nearest neighbors per particle in ppiclF.
            CALL ppiclf_user_SetYdotInit

            ! Calculate number of nearest neighbors
            ! for all particles in domain, and compare
            ! with results from ppiclF
            DO i = 1,ppiclf_npart 
                ! id1 = ppiclf_iprop(1,i)
                ! id2 = ppiclf_iprop(2,i)
                id1 = ppiclf_parts(i)%tag%partNum
                id2 = ppiclf_parts(i)%tag%rankNum

                NNCount   = 0
                NNDistSQ  = 0.0D0
                ! Loop through all particles in domain
                ! parts is the array containing all particles
                DO j = 1, totalparticles
                    IF(j .EQ. (id1+id2*particlesPerProc)) CYCLE 
                    dSQi = 0.0D0
                    dSQl = 0.0D0
                    DO l = 1,3
                        IF(ppiclf_linperiodic(l) .AND. ppiclf_EqualDomain(l)) THEN
                            dSQl = MIN( (ppiclf_parts(i)%y%pos%vec(l) - parts(j)%y%pos%vec(l))**2, ( (gridDomain(2,l) - gridDomain(1,l)) - ABS(ppiclf_parts(i)%y%pos%vec(l) - parts(j)%y%pos%vec(l)) )**2 )
                        ELSE
                            dSQl = (ppiclf_parts(i)%y%pos%vec(l) - parts(j)%y%pos%vec(l))**2
                        END IF
                        dSQi = dSQi + dSQl
                    END DO !l
                    IF(dSQi .GT. ppiclf_nndist**2) CYCLE
                    NNCount = NNCount + 1
                    NNDistSQ = NNDistSQ + dSQi
                END DO
                IF(NNCount .NE. PARTICLE_NN(i) .OR. (ABS(NNDistSQ - PPICLF_TOTNNDIST(i)) /  ABS(NNDistSQ)) .GT. 1.0D-3) THEN
                    nnpart(test) = .FALSE.
                    !              PRINT*, 'Count Diff (ppiclf,ref):', PARTICLE_NN(i),NNCount
                    !              PRINT*, 'Dist SQ percent Error:', (ABS(NNDistSQ - 
                    !     >                 PPICLF_TOTNNDIST(i)) /  ABS(NNDistSQ))
                    IF(PARTICLE_NN(i) .GT. NNCount) THEN
                        PRINT*, 'More on ppiclf, i=',i,' diff/CountDiff:', (PPICLF_TOTNNDIST(i)-NNDistSQ)/(PARTICLE_NN(i)-NNCount), 'nndist^2:',ppiclf_nndist**2, PARTICLE_NN(i), "vs", NNCount
                        print*, "tag of incorrect particle: ", ppiclf_parts(i)%tag
                    ELSE
                        PRINT*, 'More on ref, i=',i,' diff/CountDiff:', (- PPICLF_TOTNNDIST(i)+NNDistSQ)/(-PARTICLE_NN(i)+NNCount)**2 , 'nndist:',ppiclf_nndist**2, NNCount, "vs",  PARTICLE_NN(i)
                        print*, "tag of incorrect particle: ", ppiclf_parts(i)%tag
                    END IF
                    EXIT
                END IF
            END DO
        END IF
        ! Ensure nnpart(test) is true across all processors
        ! CALL MPI_ALLREDUCE(nnpart(test), nn_logical(test), 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(icomm,ierr)
    END DO !test
    CALL MPI_ALLREDUCE(nnpart, nn_logical, 8, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
!********************************************************************** 
    IF(nproc .GT. 1) THEN
        par = ' PARALLEL'
    ELSE
        par = ''
    END IF
    IF(nid .EQ. rootProc) THEN
        PRINT*, ''
        PRINT*, 'NonPeriodic Results:'

        IF(nproc .EQ. 1) THEN
            IF(nn_logical(1)) THEN
                PRINT*, ' PARTICLE NEAREST NEIGHBOR SEARCH - PASSED'
            ELSE
                PRINT*, ' PARTICLE NEAREST NEIGHBOR SEARCH - FAILED'
            END IF
        ELSE
            IF(nn_logical(1)) THEN
                PRINT*, TRIM(par) // ' PARTILCE NN SEARCH ', 'WITH GHOST PARTICLES - PASSED'
            ELSE
                PRINT*, TRIM(par) // ' PARTILCE NN SEARCH ', 'WITH GHOST PARTICLES - FAILED'
            END IF

        END IF

        PRINT*, 'Linear Periodicity Results:'


        IF (all(nn_logical(2:8))) THEN
            PRINT*, TRIM(par) // ' PARTILCE NN SEARCH WITH ', 'GHOST PARTICLES - PASSED'
        ELSE
            PRINT*, TRIM(par) // ' PARTILCE NN SEARCH WITH ','GHOST PARTICLES - FAILED'
        END IF
    END IF

    CALL MPI_FINALIZE(ierr)
END PROGRAM

!----------------------------------------------------------------------
      SUBROUTINE test_setperiodic(xf,yf,zf,i,tc)

      ! Input/Output
      ! xf - x periodic flag
      ! yf - y periodic flag
      ! zf - z periodic flag
      ! i  - test iteration
      ! tc - test case name
      IMPLICIT NONE

      INTEGER*4 xf, yf, zf, i
      CHARACTER*50 tc

      xf = 0
      yf = 0
      zf = 0

      IF(i .EQ. 1) THEN
        tc = 'NonPeriodic' 
      ELSE IF(i .EQ. 2) THEN
        tc = 'Periodic_x' 
        xf = 1
      ELSE IF(i .EQ. 3) THEN
        tc = 'Periodic_y' 
        yf = 1
      ELSE IF(i .EQ. 4) THEN
        tc = 'Periodic_z' 
        zf = 1
      ELSE IF(i .EQ. 5) THEN
        tc = 'Periodic_xy' 
        xf = 1
        yf = 1
      ELSE IF(i .EQ. 6) THEN
        tc = 'Periodic_xz' 
        xf = 1
        zf = 1
      ELSE IF(i .EQ. 7) THEN
        tc = 'Periodic_yz'
        yf = 1
        zf = 1 
      ELSE IF(i .EQ. 8) THEN
        tc = 'Periodic_xyz' 
        xf = 1
        yf = 1
        zf = 1           
      END IF

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------

 
