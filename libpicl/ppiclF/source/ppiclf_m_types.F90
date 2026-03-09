#include <PPICLF_STD.h>

module ppiclf_m_types

    
    type PPICLF_t_tag
        integer*4 :: partNum
        integer*4 :: rankNum
        integer*4 :: cycleNum
    end type PPICLF_t_tag

    type PPICLF_t_realNVec
        real*8 vec(NDIM)
    end type PPICLF_t_realNVec

    public
    
    integer :: PPICLF_t_tag_MPIH, PPICLF_t_realNVec_MPIH


    ! operator overloads for defined types.
    interface operator (+)
        module procedure nvecADDnVec
        module procedure nvecADDVec
        module procedure nvecADDScalar
    end interface
    interface operator(-)
        module procedure nvecSUBnVec
        module procedure nvecSUBVec
        module procedure nvecSUBScalar
    end interface
    interface operator(*)
        module procedure nvecMULTScalar
    end interface
    interface operator(/)
        module procedure nvecDIVScalar
    end interface

    contains

    ! nVec addition overloads
    FUNCTION nvecAddnVec(v1, v2) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1, v2
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec + v2%vec)
    end function nvecAddnVec

    FUNCTION nvecAddVec(v1, v2) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1
        real*8, intent(in)                  :: v2(NDIM)
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec + v2)
    end function nvecAddVec

    FUNCTION nvecADDScalar(v1, s1) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1
        real*8, intent(in)                  :: s1
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec + s1)
    end function nvecADDScalar

    ! nVec subtraction overloads
    FUNCTION nvecSUBnVec(v1, v2) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1, v2
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec - v2%vec)
    end function nvecSUBnVec

    FUNCTION nvecSUBVec(v1, v2) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1
        real*8, intent(in)                  :: v2(NDIM)
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec - v2)
    end function nvecSUBVec

    FUNCTION nvecSUBScalar(v1, s1) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1
        real*8, intent(in)                  :: s1
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec - s1)
    end function nvecSUBScalar

    ! nVec multiplication overloads
    FUNCTION nvecMULTScalar(v1, s1) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1
        real*8, intent(in)                  :: s1
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec * s1)
    end function nvecMULTScalar

    ! nVec division overloads
    FUNCTION nvecDIVScalar(v1, s1) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1
        real*8, intent(in)                  :: s1
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec / s1)
    end function nvecDIVScalar

    subroutine ppiclf_m_types_Create_MPI_Derivedtypes()
        integer ierr
        ! create and commit mpi derived type for PPICLF_t_realNVec
        call MPI_TYPE_CONTIGUOUS(NDIM, MPI_REAL8, PPICLF_t_realNVec_MPIH, ierr)
        call MPI_TYPE_COMMIT(PPICLF_t_realNVec_MPIH, ierr)

        ! create and commit mpi derived type for PPICLF_t_tag
        call MPI_TYPE_CONTIGUOUS(3, MPI_INTEGER4, PPICLF_t_tag_MPIH, ierr)
        call MPI_TYPE_COMMIT(PPICLF_t_tag_MPIH, ierr)
    end subroutine ppiclf_m_types_Create_MPI_Derivedtypes

    subroutine ppiclf_m_types_Destroy_MPI_Derivedtypes()
        integer ierr
        call MPI_TYPE_FREE(PPICLF_t_realNVec_MPIH, ierr)
        call MPI_TYPE_FREE(PPICLF_t_tag_MPIH, ierr)
    end subroutine ppiclf_m_types_Destroy_MPI_Derivedtypes

end module ppiclf_m_types




