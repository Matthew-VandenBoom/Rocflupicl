#include "PPICLF_STD.h"
module UT_m_data
    use ppiclf_user_particle
    ! General variables
    
    ! COMMON /Gen_Int/ 
    INTEGER*4 rootProc, nid, nproc, ierr, icomm, test


    ! COMMON /Gen_Real/ 
    REAL*8, parameter ::    PI       = 4.0D0*ATAN(1.0) ! pi
    
    real*8 randNum


    ! Grid variables
    ! COMMON /Grid_INT/ 
    INTEGER*4 nCells(3), proc_ncells, numCells, iCend, iCstart, cellsPerProc


    !COMMON /Grid_Real/ 
    REAL*8    p_grid(7,PPICLF_LEE), grid(7,PPICLF_LEE),gridDomain(2,3), gridDX(3), filter(3),nFilterCells, dx_min(3) 

    ! COMMON /Grid_Log/
    LOGICAL   gridBoundsDefined


    ! Particle variables
    ! COMMON /Part_Real/ 
    REAL*8 pdia, part_dx(3), nndist
    ! part_y(PPICLF_LRS,PPICLF_LPART), p_part_y(PPICLF_LRS,PPICLF_LPART), p_part_r(PPICLF_LRP,PPICLF_LPART),
    type(PPICLF_U_t_particle), dimension(PPICLF_LPART):: parts, p_parts


    !COMMON /Part_Int/ 
    INTEGER*4 particlesPerProc, npart_local, totalParticles, iPend, iPstart

    save
end module UT_m_data