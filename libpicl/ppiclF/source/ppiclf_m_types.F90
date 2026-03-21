#include <PPICLF_STD.h>

module ppiclf_m_types
    use mpi
    implicit none
    
    type PPICLF_t_tag
        integer*4 :: partNum
        integer*4 :: rankNum
        integer*4 :: cycleNum
    end type PPICLF_t_tag

    type PPICLF_t_realNVec
        real*8 vec(NDIM)
    end type PPICLF_t_realNVec

    ! holds the data for the nearest neighbor search for a given particle
    ! set up by ppiclf_solve_NearestNeighborSB_Start
    ! every time user calls ppiclf_solve_NearestNeighborSB, this is used to resume the
    ! iteration over all of the particles, ghost particles, and boundary points, to find the next
    ! neighbor for the particle to interact with
    type PPICLF_t_NNSB_Search_Data
        ! subbin info
        INTEGER*4 i ! particle index
        integer*4 SBt ! total # of subbins
        integer*4 SBn(3) ! number of subbins in each direction 
        integer*4 iB(3) ! this ranks bin number in each dimension  
        INTEGER*4, pointer :: SBc(:) ! number of particles in each subbin
        integer*4, pointer :: SBm(:,:) ! map of particles to subbins

        ! internal variables for the neighbor search
        integer*4 thisSB ! subbin number for the particle with index i
        logical doneParticles ! finished searching through particles, only need to search boundary points now
        logical done    ! the search is fully complete

        integer*4 iSB, jSB, kSB ! sb counters for looping over the subbins
        integer*4 k ! counter for looping over the particles in a given subbin, and for looping over the boundaries
        real*8 distSQ ! square of the NN dist. we compare against this instead of spending time doing a bunch of sqrts just to compare them
    end type PPICLF_t_NNSB_Search_Data
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
        module procedure nvecMULTnVec
        module procedure nvecMULTVec
    end interface
    interface operator(/)
        module procedure nvecDIVScalar
    end interface

    interface operator(.eq.)
        module procedure tagEQtag
    end interface

    contains

    ! nVec addition overloads
    pure FUNCTION nvecAddnVec(v1, v2) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1, v2
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec + v2%vec)
    end function nvecAddnVec

    pure FUNCTION nvecAddVec(v1, v2) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1
        real*8, intent(in)                  :: v2(NDIM)
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec + v2)
    end function nvecAddVec

    pure FUNCTION nvecADDScalar(v1, s1) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1
        real*8, intent(in)                  :: s1
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec + s1)
    end function nvecADDScalar

    ! nVec subtraction overloads
    pure FUNCTION nvecSUBnVec(v1, v2) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1, v2
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec - v2%vec)
    end function nvecSUBnVec

    pure FUNCTION nvecSUBVec(v1, v2) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1
        real*8, intent(in)                  :: v2(NDIM)
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec - v2)
    end function nvecSUBVec

    pure FUNCTION nvecSUBScalar(v1, s1) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1
        real*8, intent(in)                  :: s1
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec - s1)
    end function nvecSUBScalar

    ! nVec multiplication overloads
    pure FUNCTION nvecMULTScalar(v1, s1) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1
        real*8, intent(in)                  :: s1
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec * s1)
    end function nvecMULTScalar

    pure FUNCTION nvecMULTnVec(v1, v2) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1, v2
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec * v2%vec)
    end function nvecMULTnVec

    pure FUNCTION nvecMULTVec(v1, v2) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1
        real*8, intent(in)                  :: v2(NDIM)
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec * v2)
    end function nvecMULTVec

    ! nVec division overloads
    pure FUNCTION nvecDIVScalar(v1, s1) result(resVec)
        type(PPICLF_t_realNVec), intent(in) :: v1
        real*8, intent(in)                  :: s1
        type(PPICLF_t_realNVec)             :: resVec
        resVec = PPICLF_t_realNVec(v1%vec / s1)
    end function nvecDIVScalar

    ! other common operations for vectors
    pure function nvecMagnitude(v) result(mag)
        type(PPICLF_t_realNVec), intent(in) :: v
        real*8 mag
        mag = sqrt(nvecMagnitudeSQ(v))
    end function nvecMagnitude

    ! returns the square of the magnitude
    pure function nvecMagnitudeSQ(v) result(magSQ)
        type(PPICLF_t_realNVec), intent(in) :: v
        real*8 magSQ
        type(PPICLF_t_realNVec) temp
        temp = v * v
        magSQ = sum(temp%vec)
    end function nvecMagnitudeSQ

    pure function nvecComponentSum(v) result(compSum)
        type(PPICLF_t_realNVec), intent(in) :: v
        real*8 compSum
        compSum = sum(v%vec)
    end function nvecComponentSum

    pure function tagEQtag(t1, t2) result(equals)
        type(PPICLF_t_tag), intent(in) :: t1, t2
        logical equals
        equals = t1%partNum .eq. t2%partNum .and. &
                 t1%rankNum .eq. t2%rankNum .and. &
                 t1%cycleNum .eq. t2%cycleNum
    end function tagEQtag
    
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




