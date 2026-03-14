module ppiclf_user_random
    implicit none
    real*8, allocatable, private :: random_array(:)

    interface ppiclf_user_random_number
        ! pure module subroutine user_random_floatArray(state, targetArray)
        !     integer*8, intent(inout) :: state
        !     real*8, intent(out)      :: targetArray(:)    
        ! end subroutine user_random_floatArray
        procedure ppiclf_user_random_floatArray
    end interface

    contains
    subroutine ppiclf_user_random_setup(npart, nvals_per_part)
        integer*4, intent(in) :: npart, nvals_per_part
        if (allocated(random_array)) deallocate(random_array)
        allocate(random_array(0:npart*nvals_per_part - 1))
        call RANDOM_NUMBER(random_array)
    end subroutine ppiclf_user_random_setup
    pure subroutine ppiclf_user_random_floatArray(state, targetArray)
        integer*8, intent(inout) :: state
        real*8, intent(out)      :: targetArray(:)
        targetArray = random_array(state: state + size(targetArray) - 1)
        state = state + size(targetArray)
    end subroutine 

end module ppiclf_user_random