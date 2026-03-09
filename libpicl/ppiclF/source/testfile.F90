#include <PPICLF_STD.h>

program myTestProg
    use ppiclf_user_particle
    use ppiclf_m_types
    ! use ppiclf_m_overload
    implicit none
    type(PPICLF_t_solutionProps) solStructs(10)
    type(PPICLF_t_solutionProps) resultStruct
    integer i
    ! type(PPICLF_t_realNVec) testvec, testvec2, resultVec

    ! testvec = PPICLF_t_realNVec([1.0d0, 1.0d0, 1.0d0])
    ! testvec2 = PPICLF_t_realNVec([2.0d0, 2.0d0, 2.0d0])
    ! resultVec = testvec + testvec2
    ! print*, resultVec

    ! resultVec = testvec + 4.3d0
    ! print*, resultVec
    do i = 1, 10
        solStructs(i) = PPICLF_t_solutionProps(             &
            PPICLF_t_realNVec([i*1.0d0, i*1.0d0, i*1.0d0]), &   ! pos
            PPICLF_t_realNVec([i*2.0d0, i*2.0d0, i*2.0d0]), &   ! vel
            i*3.0d0,                                        &   ! T
            PPICLF_t_realNVec([i*4.0d0, i*4.0d0, i*4.0d0]), &   ! ang_vel
            i*5.0d0,                                        &   ! metal
            i*6.0d0                                         &   ! oxide
            )
        ! print*, solStructs(i)
    end do
    print*, solStructs(1)
    print*, solStructs(2)
    resultStruct = solStructs(1) + solStructs(2)
    print*, resultStruct


#include <PPICLF_USER_PART_VTU_MAP.h>
end program myTestProg