#include "PPICLF_STD.h"
module ppiclf_user

    use ppiclf_user_particle
    implicit none
    interface
        module subroutine ppiclf_user_SetYdot
        end subroutine ppiclf_user_SetYdot
        module subroutine ppiclf_user_MapProjPart
        ! module subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop,rprop4)
        !     real*8, intent(in) :: y    (PPICLF_LRS)
        !     real*8, intent(in) :: ydot (PPICLF_LRS)
        !     real*8, intent(in) :: ydotc(PPICLF_LRS)
        !     real*8, intent(in) :: rprop(PPICLF_LRP)
        !     real*8, intent(in) :: rprop4(PPICLF_LRP4)

        !     real*8, intent(out) :: map (PPICLF_LRP_PRO)
        end subroutine ppiclf_user_MapProjPart
        module subroutine ppiclf_user_NearestNeighbor(i)
            integer*4 i
        end subroutine ppiclf_user_NearestNeighbor
        module subroutine ppiclf_user_EvalNearestNeighbor(i,j, neighbor)!,yi,rpropi,yj,rpropj)
            integer*4 i
            integer*4 j
            type(PPICLF_U_t_ghostParticle) :: neighbor
            ! removed the extra parameters, because now i and j can be used to index into the particle array, regardless of ghost or not.
            ! real*8 yi    (PPICLF_LRS)    
            ! real*8 rpropi(PPICLF_LRP)
            ! real*8 yj    (PPICLF_LRS)    
            ! real*8 rpropj(PPICLF_LRP)
        end subroutine ppiclf_user_EvalNearestNeighbor

        module subroutine ppiclf_user_InitZero
        end subroutine ppiclf_user_InitZero

        ! this function zeros out all particle properties
        module pure subroutine ppiclf_user_ZeroParticle(particle)
            type(PPICLF_U_t_particle), intent(inout) :: particle
        end subroutine ppiclf_user_ZeroParticle
        module pure subroutine ppiclf_user_ZeroInterp(interp)
            type(PPICLF_U_t_interp), intent(inout) :: interp
        end subroutine ppiclf_user_ZeroInterp
    end interface
end module ppiclf_user


