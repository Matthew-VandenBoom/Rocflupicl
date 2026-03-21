#include <PPICLF_STD.h>

module ppiclf_m_transfers
    use mpi 

    ! used data/procedures
    use ppiclf_data, only: ppiclf_np, ppiclf_nid
    use ppiclf_op, only: ppiclf_exittr

    ! types we want to define an alltoallv for
    use ppiclf_user_particle, only: PPICLF_U_t_particle, PPICLF_U_t_ghostParticle
    use ppiclf_m_wrapped_types, only: ppiclf_t_fluidCell_wrapped, ppiclf_t_interp_wrapped, ppiclf_t_feedback_wrapped
    implicit none 
    type PPICLF_t_MPI_tags
        integer*4 ppiclf_alltoallv_sync_counts
        integer*4 ppiclf_alltoallv_transfer
    end type PPICLF_t_MPI_tags

    interface ppiclf_alltoallv
        procedure ppiclf_alltoallv_int, ppiclf_alltoallv_real, ppiclf_alltoallv_PPICLF_U_t_particle, ppiclf_alltoallv_PPICLF_U_t_ghostParticle, ppiclf_alltoallv_ppiclf_t_fluidCell_wrapped, ppiclf_alltoallv_ppiclf_t_interp_wrapped, ppiclf_alltoallv_ppiclf_t_feedback_wrapped
    end interface ppiclf_alltoallv

    type(PPICLF_t_MPI_tags), parameter :: PPICLF_MPI_tags = PPICLF_t_MPI_tags(  &
                                            ppiclf_alltoallv_sync_counts = 1,    &
                                            ppiclf_alltoallv_transfer = 2       )
    contains

! generate the definition of ppiclf_alltoallv for Integer4
#define GENERIC_TYPE_A_extention(function_name) CAT(function_name,_int)
#define GENERIC_TYPE_A_fullname INTEGER*4
#include "Generic_Procedures/ppiclf_alltoallv.subroutine.F90"
#undef GENERIC_TYPE_A_extention
#undef GENERIC_TYPE_A_fullname

! definition of ppiclf_alltoallv for real8
#define GENERIC_TYPE_A_extention(function_name) CAT(function_name,_real)
#define GENERIC_TYPE_A_fullname REAL*8
#include "Generic_Procedures/ppiclf_alltoallv.subroutine.F90"
#undef GENERIC_TYPE_A_extention
#undef GENERIC_TYPE_A_fullname


! generate ppiclf_alltoallv for PPICLF_U_t_particle
#define GENERIC_TYPE_A_extention(function_name) CAT(function_name,_PPICLF_U_t_particle)
#define GENERIC_TYPE_A_fullname PPICLF_U_t_particle
#include "Generic_Procedures/ppiclf_alltoallv.subroutine.F90"
#undef GENERIC_TYPE_A_extention
#undef GENERIC_TYPE_A_fullname

! generate ppiclf_alltoallv for PPICLF_U_t_ghostParticle
#define GENERIC_TYPE_A_extention(function_name) CAT(function_name,_PPICLF_U_t_ghostParticle)
#define GENERIC_TYPE_A_fullname PPICLF_U_t_ghostParticle
#include "Generic_Procedures/ppiclf_alltoallv.subroutine.F90"
#undef GENERIC_TYPE_A_extention
#undef GENERIC_TYPE_A_fullname

! generate ppiclf_alltoallv for ppiclf_t_fluidCell_wrapped
#define GENERIC_TYPE_A_extention(function_name) CAT(function_name,_ppiclf_t_fluidCell_wrapped)
#define GENERIC_TYPE_A_fullname ppiclf_t_fluidCell_wrapped
#include "Generic_Procedures/ppiclf_alltoallv.subroutine.F90"
#undef GENERIC_TYPE_A_extention
#undef GENERIC_TYPE_A_fullname

! generate ppiclf_alltoallv for ppiclf_t_interp_wrapped
#define GENERIC_TYPE_A_extention(function_name) CAT(function_name,_ppiclf_t_interp_wrapped)
#define GENERIC_TYPE_A_fullname ppiclf_t_interp_wrapped
#include "Generic_Procedures/ppiclf_alltoallv.subroutine.F90"
#undef GENERIC_TYPE_A_extention
#undef GENERIC_TYPE_A_fullname

! generate ppiclf_alltoallv for ppiclf_t_feedback_wrapped
#define GENERIC_TYPE_A_extention(function_name) CAT(function_name,_ppiclf_t_feedback_wrapped)
#define GENERIC_TYPE_A_fullname ppiclf_t_feedback_wrapped
#include "Generic_Procedures/ppiclf_alltoallv.subroutine.F90"
#undef GENERIC_TYPE_A_extention
#undef GENERIC_TYPE_A_fullname

end module ppiclf_m_transfers