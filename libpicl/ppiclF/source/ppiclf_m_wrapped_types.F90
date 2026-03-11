#include <PPICLF_STD.h>

module ppiclf_m_wrapped_types
    use ppiclf_m_types
    use ppiclf_user_particle
    use mpi
    implicit none

    ! this stores the location/bounds of the fluid grid, along with the info to map the cells back to their home rank/cell index
    ! Indicies 1-3: Centroid x,y,z position
    ! Indicies 4-6: Max cell dx,dy,dz based on any vertex combination
    ! Index      7: Cell Volume 
    type ppiclf_t_fluidCell_wrapped
        integer*4 homeCellIndex ! index of the fluid cell on the home rank
        integer*4 homeRank      ! home rank for the fluid cell (when setting up for transfer, it will temporarily have the destination rank, but this is changed before it is transfered)
        REAL*8    FluidCell(7)  ! fluid cell location/bounds, see above
    end type ppiclf_t_fluidCell_wrapped


    ! this represents the interp data for a single fluid cell, with a bit of extra
    ! info to allow it to be associated with the correct fluid cell location data
    type ppiclf_t_interp_wrapped
        integer*4 homeCellIndex ! index of the fluid cell on the home rank
        integer*4 homeRank      ! home rank for the fluid cell (when setting up for transfer, it will temporarily have the destination rank, but this is changed before it is transfered)
        type(PPICLF_U_t_interp) interp
    end type ppiclf_t_interp_wrapped

    ! this represents the feedback data for a single fluid cell, with a bit of extra
    ! info to allow it to be associated with the correct fluid cell location data
    type ppiclf_t_feedback_wrapped
        integer*4 homeCellIndex ! index of the fluid cell on the home rank
        integer*4 homeRank      ! home rank for the fluid cell
        type(PPICLF_U_t_feedback) feedback
    end type ppiclf_t_feedback_wrapped


    integer, public, save :: ppiclf_t_fluidCell_wrapped_MPIH, ppiclf_t_interp_wrapped_MPIH, ppiclf_t_feedback_wrapped_MPIH
    integer, private, save :: ppiclf_t_fluidCell_MPIH
    interface
        module subroutine ppiclf_m_wrapped_Create_MPI_DerivedTypes
            
        end subroutine ppiclf_m_wrapped_Create_MPI_DerivedTypes
    end interface

    contains
        subroutine ppiclf_m_wrapped_Destroy_MPI_DerivedTypes
            integer ierr
            call MPI_TYPE_FREE(ppiclf_t_fluidCell_MPIH, ierr)
            call MPI_TYPE_FREE(ppiclf_t_fluidCell_wrapped_MPIH, ierr)
            call MPI_TYPE_FREE(ppiclf_t_interp_wrapped_MPIH, ierr)
            call MPI_TYPE_FREE(ppiclf_t_feedback_wrapped_MPIH, ierr)
        end subroutine ppiclf_m_wrapped_Destroy_MPI_DerivedTypes
end module ppiclf_m_wrapped_types