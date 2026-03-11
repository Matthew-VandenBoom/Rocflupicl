#include <PPICLF_STD.h>

module ppiclf_m_particle_ops

    use ppiclf_data, only: ppiclf_totalBins, ppiclf_npart, ppiclf_np

    use ppiclf_user_particle, only: PPICLF_U_t_particle, PPICLF_U_t_ghostParticle
    implicit none
    private

    public :: ppiclf_particles_GroupBy
    public :: ppiclf_particles_GroupByIntoArray
    public :: ppiclf_ghostParticles_GroupBy
    public :: ppiclf_ghostParticles_GroupByIntoArray

    type GroupByKeys_t
        integer BinNum, RankNum
    end type

    type(GroupByKeys_t), public, parameter :: GroupByKeys = GroupByKeys_t(0, 1)

    contains

    ! "sorts" the particles based on one of the elements of iprop. Since this isn't a full sort, it is done in O(n) time
    ! takes the particles as an argument, since it allows the subroutine to be declared as pure
    ! TODO: when merging with averys code, this needs to be changed to update the start/end indexes for the part to bin map
    pure subroutine ppiclf_particles_GroupBy(particles, keys, ierr)
        type(PPICLF_U_t_particle), intent(inout) :: particles(PPICLF_LPART)
        integer, intent(in) :: keys(:)
        integer, intent(out) :: ierr
        !
        ! Internal:
        !
        integer i
        type(PPICLF_U_t_particle), allocatable :: tempParticles(:)
        allocate(tempParticles(ppiclf_npart))

        do i = 1, size(keys)
            select case (keys(i))
            case (GroupByKeys%BinNum)
                call ppiclf_particles_GroupBy_Bin(particles, tempParticles)
            case (GroupByKeys%RankNum)
                call ppiclf_particles_GroupBy_Rank(particles, tempParticles)
            case default
                ierr = 1
            end select 
            particles(1:ppiclf_npart) = tempParticles(1:ppiclf_npart)
        end do
        ierr = 0
        deallocate(tempParticles)
    end subroutine ppiclf_particles_GroupBy

    ! this version of the function gets the particle data from ppiclf_parts, and stores the grouped output into the passed array
    pure subroutine ppiclf_particles_GroupByIntoArray(inParticles, outParticles, keys, ierr)
        type(PPICLF_U_t_particle), intent(in)  :: inParticles(PPICLF_LPART)

        type(PPICLF_U_t_particle), intent(out) :: outParticles(:)
        integer, intent(in) :: keys(:)
        integer, intent(out) :: ierr
        !
        ! Internal:
        !
        integer i
        type(PPICLF_U_t_particle), allocatable :: tempParticles(:)
        allocate(tempParticles(ppiclf_npart))

        if (size(outParticles) < ppiclf_npart) then
            ierr = 2
            return
        end if
        tempParticles(1:ppiclf_npart)  = inParticles(1:ppiclf_npart)
        do i = 1, size(keys)
            if (i .ne. 1) then
                tempParticles(1:ppiclf_npart)  = outParticles(1:ppiclf_npart)
            end if
            select case (keys(i))
            case (GroupByKeys%BinNum)
                call ppiclf_particles_GroupBy_Bin(tempParticles, outParticles)
            case (GroupByKeys%RankNum)
                call ppiclf_particles_GroupBy_Rank(tempParticles, outParticles)
            case default
                ierr = 1
            end select 
            
        end do
        ierr = 0

        deallocate(tempParticles)
    end subroutine ppiclf_particles_GroupByIntoArray


    pure subroutine ppiclf_particles_GroupBy_Bin(inParticles, outParticles)
        type(PPICLF_U_t_particle), intent(in)  :: inParticles(PPICLF_LPART)

        type(PPICLF_U_t_particle), intent(out) :: outParticles(:)

        
        ! (:,1) counts the number of particles in the bin
        ! (:,2) start index in the sorted array
        ! (:,3) counter when sorting
        integer :: binInfo(0:ppiclf_totalBins, 3)
        
        integer :: i, partBin

        binInfo = 0
        
        do i = 1, ppiclf_npart
            binInfo(inParticles(i)%iprop%binNum, 1) = binInfo(inParticles(i)%iprop%binNum, 1) + 1
        end do

        do i = 1, ppiclf_totalBins
            ! this start index = prev start index + n prev
            binInfo(i,2)  = binInfo(i-1, 2) + binInfo(i-1, 1)
        end do

        binInfo(:, 3) = binInfo(:, 2) + 1
        do i = 1, ppiclf_npart
            partBin = inParticles(i)%iprop%binNum
            outParticles(binInfo(partBin, 3)) = inParticles(i)
            binInfo(partBin, 3) = binInfo(partBin, 3) + 1
        end do

    end subroutine ppiclf_particles_GroupBy_Bin

    pure subroutine ppiclf_particles_GroupBy_Rank(inParticles, outParticles)
        type(PPICLF_U_t_particle), intent(in)  :: inParticles(PPICLF_LPART)

        type(PPICLF_U_t_particle), intent(out) :: outParticles(:)


        
        ! (:,1) counts the number of particles in the rank
        ! (:,2) start index in the sorted array
        ! (:,3) counter when sorting
        integer :: rankInfo(0:ppiclf_np, 3)
        
        integer :: i, partBin

        rankInfo = 0
        
        do i = 1, ppiclf_npart
            rankInfo(inParticles(i)%iprop%ParticleRank, 1) = rankInfo(inParticles(i)%iprop%ParticleRank, 1) + 1
        end do

        do i = 1, ppiclf_totalBins
            ! this start index = prev start index + n prev
            rankInfo(i,2)  = rankInfo(i-1, 2) + rankInfo(i-1, 1)
        end do

        rankInfo(:, 3) = rankInfo(:, 2) + 1
        do i = 1, ppiclf_npart
            partBin = inParticles(i)%iprop%ParticleRank
            outParticles(rankInfo(partBin, 3)) = inParticles(i)
            rankInfo(partBin, 3) = rankInfo(partBin, 3) + 1
        end do

    end subroutine ppiclf_particles_GroupBy_Rank


    ! "sorts" the particles based on one of the elements of iprop. Since this isn't a full sort, it is done in O(n) time
    ! takes the particles as an argument, since it allows the subroutine to be declared as pure
    ! TODO: when merging with averys code, this needs to be changed to update the start/end indexes for the part to bin map
    pure subroutine ppiclf_ghostParticles_GroupBy(particles, keys, ierr)
        type(PPICLF_U_t_ghostParticle), intent(inout) :: particles(PPICLF_LPART)
        integer, intent(in) :: keys(:)
        integer, intent(out) :: ierr
        !
        ! Internal:
        !
        integer i
        type(PPICLF_U_t_ghostParticle), allocatable :: tempParticles(:)
        allocate(tempParticles(ppiclf_npart))

        do i = 1, size(keys)
            select case (keys(i))
            case (GroupByKeys%BinNum)
                call ppiclf_ghostParticles_GroupBy_Bin(particles, tempParticles)
            case (GroupByKeys%RankNum)
                call ppiclf_ghostParticles_GroupBy_Rank(particles, tempParticles)
            case default
                ierr = 1
            end select 
            particles(1:ppiclf_npart) = tempParticles(1:ppiclf_npart)
        end do
        ierr = 0
        deallocate(tempParticles)
    end subroutine ppiclf_ghostParticles_GroupBy

    ! this version of the function gets the particle data from ppiclf_parts, and stores the grouped output into the passed array
    pure subroutine ppiclf_ghostParticles_GroupByIntoArray(inParticles, outParticles, keys, ierr)
        type(PPICLF_U_t_ghostParticle), intent(in)  :: inParticles(PPICLF_LPART)

        type(PPICLF_U_t_ghostParticle), intent(out) :: outParticles(:)
        integer, intent(in) :: keys(:)
        integer, intent(out) :: ierr
        !
        ! Internal:
        !
        integer i
        type(PPICLF_U_t_ghostParticle), allocatable :: tempParticles(:)
        allocate(tempParticles(ppiclf_npart))

        if (size(outParticles) < ppiclf_npart) then
            ierr = 2
            return
        end if
        tempParticles(1:ppiclf_npart)  = inParticles(1:ppiclf_npart)
        do i = 1, size(keys)
            if (i .ne. 1) then
                tempParticles(1:ppiclf_npart)  = outParticles(1:ppiclf_npart)
            end if
            select case (keys(i))
            case (GroupByKeys%BinNum)
                call ppiclf_ghostParticles_GroupBy_Bin(tempParticles, outParticles)
            case (GroupByKeys%RankNum)
                call ppiclf_ghostParticles_GroupBy_Rank(tempParticles, outParticles)
            case default
                ierr = 1
            end select 
            
        end do
        ierr = 0

        deallocate(tempParticles)
    end subroutine ppiclf_ghostParticles_GroupByIntoArray


    pure subroutine ppiclf_ghostParticles_GroupBy_Bin(inParticles, outParticles)
        type(PPICLF_U_t_ghostParticle), intent(in)  :: inParticles(PPICLF_LPART)

        type(PPICLF_U_t_ghostParticle), intent(out) :: outParticles(:)

        
        ! (:,1) counts the number of particles in the bin
        ! (:,2) start index in the sorted array
        ! (:,3) counter when sorting
        integer :: binInfo(0:ppiclf_totalBins, 3)
        
        integer :: i, partBin

        binInfo = 0
        
        do i = 1, ppiclf_npart
            binInfo(inParticles(i)%iprop%binNum, 1) = binInfo(inParticles(i)%iprop%binNum, 1) + 1
        end do

        do i = 1, ppiclf_totalBins
            ! this start index = prev start index + n prev
            binInfo(i,2)  = binInfo(i-1, 2) + binInfo(i-1, 1)
        end do

        binInfo(:, 3) = binInfo(:, 2) + 1
        do i = 1, ppiclf_npart
            partBin = inParticles(i)%iprop%binNum
            outParticles(binInfo(partBin, 3)) = inParticles(i)
            binInfo(partBin, 3) = binInfo(partBin, 3) + 1
        end do

    end subroutine ppiclf_ghostParticles_GroupBy_Bin

    pure subroutine ppiclf_ghostParticles_GroupBy_Rank(inParticles, outParticles)
        type(PPICLF_U_t_ghostParticle), intent(in)  :: inParticles(PPICLF_LPART)

        type(PPICLF_U_t_ghostParticle), intent(out) :: outParticles(:)


        
        ! (:,1) counts the number of particles in the rank
        ! (:,2) start index in the sorted array
        ! (:,3) counter when sorting
        integer :: rankInfo(0:ppiclf_np, 3)
        
        integer :: i, partBin

        rankInfo = 0
        
        do i = 1, ppiclf_npart
            rankInfo(inParticles(i)%iprop%ParticleRank, 1) = rankInfo(inParticles(i)%iprop%ParticleRank, 1) + 1
        end do

        do i = 1, ppiclf_totalBins
            ! this start index = prev start index + n prev
            rankInfo(i,2)  = rankInfo(i-1, 2) + rankInfo(i-1, 1)
        end do

        rankInfo(:, 3) = rankInfo(:, 2) + 1
        do i = 1, ppiclf_npart
            partBin = inParticles(i)%iprop%ParticleRank
            outParticles(rankInfo(partBin, 3)) = inParticles(i)
            rankInfo(partBin, 3) = rankInfo(partBin, 3) + 1
        end do

    end subroutine ppiclf_ghostParticles_GroupBy_Rank

end module ppiclf_m_particle_ops