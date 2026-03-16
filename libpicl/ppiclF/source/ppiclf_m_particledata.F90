#include "PPICLF_STD.h"

module ppiclf_m_particledata
    use ppiclf_m_types
    use ppiclf_user_particle
    implicit none

    type(PPICLF_U_t_particle), save :: ppiclf_parts(PPICLF_LPART)
    type(PPICLF_U_t_ghostParticle), save :: ppiclf_gparts(PPICLF_LPART_GP)

    type(PPICLF_U_t_feedback), save :: ppiclf_part_feedback(PPICLF_LPART)
    contains

    ! TODO: Move this to a user file
    subroutine CopyRealToGhost(particle, ghost)
        type(PPICLF_U_t_particle), intent(in) :: particle
        type(PPICLF_U_t_ghostParticle), intent(out) :: ghost

        ghost%tag   = particle%tag
        ghost%Y     = particle%Y
        ghost%rprop = particle%rprop
        ghost%iprop = particle%iprop
    end subroutine CopyRealToGhost
end module ppiclf_m_particledata