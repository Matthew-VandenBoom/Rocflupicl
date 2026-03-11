submodule(ppiclf_m_wrapped_types) ppiclf_m_wrapped_types_imp
    use ppiclf_m_comm, only: user_particle_MPIH, user_ghost_MPIH, user_interp_MPIH, user_feedback_MPIH
    use mpi
    implicit none
    contains
    module procedure ppiclf_m_wrapped_Create_MPI_DerivedTypes
        !
        ! Internal:
        !
        integer :: extent_Real8
        integer :: i, ierr
        integer :: oldtypes(20), blockcounts(20), offsets(20)
        call MPI_TYPE_EXTENT(MPI_REAL8, extent_Real8, ierr)
        ! create and commit MPI derived type for ppiclf_t_fluidCell
        ! not actually a real type, just a type from mpi's perspective
        ! create and commit mpi derived type for PPICLF_t_realNVec
        call MPI_TYPE_CONTIGUOUS(7, MPI_REAL8, ppiclf_t_fluidCell_MPIH, ierr)
        call MPI_TYPE_COMMIT(ppiclf_t_fluidCell_MPIH, ierr)

        ! create and commit MPI derived type for ppiclf_t_fluidCell_wrapped
        i = 1
        ! 2 int4s (homeCellIndex, homeRank)
        oldtypes(i) = MPI_REAL8
        blockcounts(i) = 2
        offsets(i) = 0
        i = i + 1

        ! 1 fluidcell
        oldtypes(i) = ppiclf_t_fluidCell_MPIH
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_Real8

        call MPI_TYPE_STRUCT(i, blockcounts(1:i), offsets(1:i), oldtypes(1:i), ppiclf_t_fluidCell_wrapped_MPIH, ierr)
        call MPI_TYPE_COMMIT(ppiclf_t_fluidCell_wrapped_MPIH, ierr)

        ! create and commit MPI derived type for ppiclf_t_interp_wrapped
        i = 1
        ! 2 int4s (homeCellIndex, homeRank)
        oldtypes(i) = MPI_REAL8
        blockcounts(i) = 2
        offsets(i) = 0
        i = i + 1

        ! 1 t_interp
        oldtypes(i) = user_interp_MPIH
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_Real8

        call MPI_TYPE_STRUCT(i, blockcounts(1:i), offsets(1:i), oldtypes(1:i), ppiclf_t_interp_wrapped_MPIH, ierr)
        call MPI_TYPE_COMMIT(ppiclf_t_interp_wrapped_MPIH, ierr)

        ! create and commit MPI derived type for ppiclf_t_feedback_wrapped
        i = 1
        ! 2 int4s (homeCellIndex, homeRank)
        oldtypes(i) = MPI_REAL8
        blockcounts(i) = 2
        offsets(i) = 0
        i = i + 1

        ! 1 feedback
        oldtypes(i) = user_feedback_MPIH
        blockcounts(i) = 1
        offsets(i) = offsets(i - 1) + blockcounts(i - 1) * extent_Real8

        call MPI_TYPE_STRUCT(i, blockcounts(1:i), offsets(1:i), oldtypes(1:i), ppiclf_t_fluidCell_wrapped_MPIH, ierr)
        call MPI_TYPE_COMMIT(ppiclf_t_fluidCell_wrapped_MPIH, ierr)
    end procedure ppiclf_m_wrapped_Create_MPI_DerivedTypes
end submodule ppiclf_m_wrapped_types_imp