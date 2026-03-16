#include "PPICLF_STD.h" 
! PPICLF_STD has been modified to include PPICLF_USER.h

module ppiclf_data
    use ppiclf_user_particle, only: PPICLF_U_t_interp, PPICLF_U_t_interp, PPICLF_U_t_feedback
    use ppiclf_m_wrapped_types, only: ppiclf_t_fluidCell_wrapped, ppiclf_t_interp_wrapped, ppiclf_t_feedback_wrapped
    implicit none

    save
    
    ! originally PPICLF_SOLN.h
    ! computational particl data
    ! REAL*8 PPICLF_Y     (PPICLF_LRS ,PPICLF_LPART)  &   ! Solution
    !     ,PPICLF_YDOT  (PPICLF_LRS ,PPICLF_LPART)    &   ! Total solution RHS
    !     ,PPICLF_YDOTC (PPICLF_LRS ,PPICLF_LPART)    &   ! Coupled solution RHS
    !     ,PPICLF_RPROP (PPICLF_LRP ,PPICLF_LPART)    &   ! Real particle properties
    !     ,PPICLF_RPROP2(PPICLF_LRP2,PPICLF_LPART)    &   ! Secondary real particle properties
    !     ,PPICLF_RPROP3(PPICLF_LRP3,PPICLF_LPART)    &   ! Third real particle properties (Was used for VU force model)
    !     ,PPICLF_RPROP4(PPICLF_LRP4,PPICLF_LPART)    &   ! Fourth real particle properties
    !     ,PPICLF_RPROP5(PPICLF_LRP5,PPICLF_LPART)    &   ! Fifth real particle properties
    !     ,PPICLF_FEEDBK(PPICLF_LRP_PRO,PPICLF_LPART)     !Feedback particle terms

    ! INTEGER*4 PPICLF_IPROP(PPICLF_LIP,PPICLF_LPART) ! Integer particle properties

    INTEGER*4 PPICLF_NPART 

    ! Previous time step solutions, may grow later
    ! REAL*8 PPICLF_Y1(PPICLF_LRS, PPICLF_LPART)

    ! Previous time step solutions, may grow later
    ! REAL*8 PPICLF_TIMEBH(PPICLF_VU)
    ! REAL*8 PPICLF_DRUDTPLAG(3,PPICLF_VU,PPICLF_LPART)
    ! REAL*8 PPICLF_DRUDTMIXT(3,PPICLF_VU,PPICLF_LPART)

    

    ! Originally PPICLF_GRID.h

    ! Grid
    REAL*8  PPICLF_XDRANGE(2,3)

    ! Projected info for the cells in the local rank,
    ! in the correct order for access by the fluid solver (via solve_GetProFld)
    ! ie same ordering as was originally given by InitOverlapGrid
    ! REAL*8 PPICLF_PRO_FLD(PPICLF_LEE,PPICLF_LRP_PRO)
    type(PPICLF_U_t_feedback) ppiclf_pro_fld(PPICLF_LEE)
    
    ! projection info for cells, used to transfer all projection info back to home rank
    ! before transfer, particle feed back is projected, into this array
    ! after transfer, the order is wrong, and there can be multiple indicies for a given cell,
    ! ppiclf_cell_map_proj is used to combine all the data and put it in the proper order in PPICLF_PRO_FLD
    type(ppiclf_t_feedback_wrapped) PPICLF_PRO_FLD_PICL(PPICLF_LEE)
    
    ! Holds the interpolation data passed in by the driver files until solve_InterpField is called by solve_InitSolve or solve_SetYdot
    ! cell indexing matches PPICLF_FLUID_GRID
    type(PPICLF_U_t_interp) PPICLF_INT_FLD_INPUT(PPICLF_LEE)
    

    ! Contains the mapped fluid cell location data
    ! set from received data in comm_MapOverlapGrid, ppiclf_cell_map maps the indicies in this back to their home rank/index
    ! used in ppiclf_solve, many places
    type(ppiclf_t_fluidCell_wrapped) PPICLF_PICL_GRID(PPICLF_LEE)
    
    ! precomputed distance from a particle to the center of each of its nnearest cells
    ! nnearest can be found in ppiclf_npart2cell, and the cell indicies can be found in PPICLF_PART2cell_map
    REAL*8 PPICLF_PART2CELL_DIST(PPICLF_LPART,27)
    
    
        
    ! this stores the location/bounds of the fluid grid from the local rank's fluid solver
    ! Indicies 1-3: Centroid x,y,z position
    ! Indicies 4-6: Max cell dx,dy,dz based on any vertex combination
    ! Index      7: Cell Volume 
    ! Set in comm_InitOverlapGrid
    REAL*8    PPICLF_FLUID_GRID (7,PPICLF_LEE)
    
    ! numbuer of fluid cells from the local rank's fluid solver (len of PPICLF_FLUID_GRID)
    INTEGER*4 PPICLF_NFVCELLS

    ! Stores interpolation data for use in solve_Interpolate
    ! before transfer, is filled from PPICLF_INT_FLD_INPUT using ppiclf_cell_map_interp in solve_InterpField
    ! it is transfered in solve_InterpTupleTransfer
    ! it is only used in solve_interpolate
    type(ppiclf_t_interp_wrapped) PPICLF_INT_FLD(PPICLF_LEE)

    ! Stores cell to rank mapping. From this rank's cells to the ranks that need them for computation.
    ! It is always sorted such that all of the cells going to a given rank are contiguous, and they are sorted in increasing order of rank num
    ! Index      1: Fluid solver cell id (index of the cell within home rank's PPICLF_FLUID_GRID)
    ! Index      2: Home rank # (ppiclf_nid of home rank)
    ! Index      3: Destination rank # (rank where this is being sent)
    ! Indicies 4-6: Bin indicies
    ! Filled in ppiclf_comm_MapOverlapGrid, sorted by ppiclf_CELL_MAP_GroupByDestRank, then used to fill ppiclf_picl_grid with real data, then sent out
    INTEGER*4 PPICLF_CELL_MAP(PPICLF_LRMAX,PPICLF_LEE)

    ! saved copies of the counts and displacements of the mapped cells within PPICLF_CELL_MAP for each of the ranks.
    ! allocated in comm_InitMPI with len=ppiclf_np + 1 (number of ranks indexed (0:ppiclf_np) ) and deallocated in comm_FinalizeMPI
    ! For use in MPI_ALLTOALLV
    INTEGER*4, allocatable :: PPICLF_CELL_MAP_SENDCOUNTS(:), PPICLF_CELL_MAP_SENDDISPS(:), PPICLF_CELL_MAP_RECVCOUNTS(:), PPICLF_CELL_MAP_RECVDISPS(:)

    ! backup copy of PPICLF_CELL_MAP from before transfer.
    ! used to create ppiclf_cell_map_interp
    ! INTEGER*4 PPICLF_CELL_MAP_ORIG(PPICLF_LRMAX,PPICLF_LEE)

    ! Indirect copy of PPICLF_CELL_MAP, used for transfering interpolation data
    ! needs to be the same as PPICLF_CELL_MAP was originally, so that interp data goes to the 
    ! same ranks as the coresponding cell location data, and importantly, has the same index.
    ! INTEGER*4 PPICLF_CELL_MAP_INTERP(PPICLF_LRMAX,PPICLF_LEE)

    ! Copied from PPICLF_CELL_MAP_INTERP in ppiclf_solve_ProjectParticleGrid.
    ! Not sure why it is coppied from map_interp instead of just PPICLF_CELL_MAP,
    !  I think they should be the same at this point in time
    ! Used to transfer projection data back to the rank where the cells came from,
    ! and after the transfer maps the projection data back to its original cell
    ! INTEGER*4 PPICLF_CELL_MAP_PROJ(PPICLF_LRMAX,PPICLF_LEE)

    ! # of mapped fluid cells recieved from other ranks
    INTEGER*4 PPICLF_NCELLS_FV2PICL
    ! # of mapped fluid cells sent to other ranks, set from PPICLF_NCELLS_FV2PICL before transfering in MapOverlapGrid
    ! used to set PPICLF_NCELLS_INTERP in solve_InitInterp
    INTEGER*4 PPICLF_NCELLS_FV2PICL_SENT !PPICLF_NCELLS_FV2PICL_ORIG

    ! = PPICLF_NCELLS_FV2PICL_ORIG
    ! used for array indexing when setting up interpolation data for transfer
    ! After transfer is equal to PPICLF_NCELLS_FV2PICL, and is used when looping over mapped cells in solve_SBParticleToCellMap
    ! INTEGER*4 PPICLF_NCELLS_INTERP

    ! = ppiclf_ncells_FV2PICL
    ! Used for array indexing when setting up projection data for transfer
    ! After transfer, used when setting ppiclf_pro_fld from recieved data
    ! INTEGER*4 PPICLF_NCELLS_PROJ

    ! number of cells in PART2CELL_MAP and PART2CELL_dist for a given particle
    INTEGER*4 PPICLF_NPART2CELL(PPICLF_LPART)

    ! cell index (into mapped data) the nearest cells to each particle
    INTEGER*4 PPICLF_PART2CELL_MAP(PPICLF_LPART,27)

    


    ! INTEGER*4 PPICLF_INT_ICNT                   ! number of fields that have had fluid data passed
    ! INTEGER*4 PPICLF_INT_MAP(PPICLF_LRP_INT)    ! map from the index of ppiclf_int_fluid_input to index of rprop


    ! Originally PPICLF_OPT.h
    ! Particle options
    LOGICAL PPICLF_RESTART, PPICLF_OVERLAP, PPICLF_LCOMM, PPICLF_LINIT, PPICLF_LINTP, PPICLF_LPROJ, PPICLF_LSUBSUBBIN,PPICLF_EQUALDOMAIN(3), PPICLF_LINPERIODIC(3), &
        PPICLF_REMOVE_PARTICLE, PPICLF_BINCHANGED, PPICLF_PRINTBINVTU, PPICLF_READYTOSOLVE

    DATA PPICLF_LCOMM /.false./
    DATA PPICLF_RESTART /.false./

    INTEGER*4 PPICLF_NDIM, PPICLF_IMETHOD, PPICLF_NGRIDS, PPICLF_CYCLE, PPICLF_IOSTEP, PPICLF_IENDIAN, PPICLF_IWALLM

    REAL*8 PPICLF_FILTER(3), PPICLF_RK3COEF(3,3), PPICLF_DT, PPICLF_TIME, PPICLF_NNDIST, PPICLF_INTERP_DCHK(3), PPICLF_TOTNNDIST(PPICLF_LPART)
    REAL*8 PPICLF_RK3ARK(3)


    ! Originally PPICLF_PARALLEL.h
    ! Communication

    INTEGER*4 PPICLF_COMM, PPICLF_NP, PPICLF_NID, PPICLF_GLNPART, PPICLF_CR_HNDL,PPICLF_FP_HNDL, PPICLF_COMM_NID
    DATA PPICLF_NID /0/

    ! Bins
    INTEGER*4 PPICLF_N_BINS(3),PPICLF_TOTALBINS, PPICLF_IMEDSLICE(2), PPICLF_ILARSLICE(2)

    REAL*8 PPICLF_BINS_DX(3), PPICLF_BINB(6), PPICLF_BIN_POS(2,3),PPICLF_PREVIOUSBINB(6), PPICLF_BINDOMLEN(3), PPICLF_RMEDSLICE(2)
      

    ! Ghost particles
    ! REAL*8 PPICLF_RPROP_GP(PPICLF_LRP_GP,PPICLF_LPART_GP),PPICLF_CP_MAP(PPICLF_LRP_GP,PPICLF_LPART)


    ! INTEGER*4 PPICLF_IPROP_GP(PPICLF_LIP_GP,PPICLF_LPART_GP),
    INTEGER*4 PPICLF_PARTICLEMOVED

    

    INTEGER*4  PPICLF_NPART_GP

    LOGICAL PPICLF_PPINTERACTIONS

    INTEGER*4 PARTICLE_NN(PPICLF_LPART)


    ! Originally PPICLF_GEOM.h

    ! Wall support
    REAL*8 PPICLF_WALL_C(9,PPICLF_LWALL),PPICLF_WALL_N(4,PPICLF_LWALL)

    INTEGER*4  PPICLF_NWALL

    ! Originally PPICLF_PERFORMACE.h

    ! Time per operation per stage
    REAL*8  PPICLF_TCreateBin               &
        ,PPICLF_TSendParticles              &
        ,PPICLF_TSendGridOverlap            &
        ,PPICLF_TSendFluidFields            &
        ,PPICLF_TParticleParticleModels     &
        ,PPICLF_TFluidParticleModels        &
        ,PPICLF_TSendGhostParticles         &
        ,PPICLF_TMapParticlesCells          &
        ,PPICLF_TInterpolation              &
        ,PPICLF_TProjection                 &
        ,PPICLF_TWriteSolution              &
        ,PPICLF_TIntegration                &
        ,PPICLF_TPeriodicity                &
        ,PPICLF_TDataTransfers              &
        ,PPICLF_TTotal                      

    ! Originally PPICLF_USER_COMMON.h
    ! this seems like really poor design, I dont think the user code should be defining variables used by 
    ! PPICLF core (configuring things like array sizes via defines are not the same)

    ! For ppiclf_solve_InitAngularPeriodic
    integer*4 x_per_flag, y_per_flag, z_per_flag, ang_per_flag,ang_case 
    real*8 ang_per_angle, ang_per_xangle,ang_per_rin, ang_per_rout,xrot(3) , vrot(3)
    real*8 x_per_min, x_per_max,y_per_min, y_per_max, z_per_min, z_per_max

end module ppiclf_data