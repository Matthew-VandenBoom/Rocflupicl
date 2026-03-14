#include "PPICLF_STD.h"
module ppiclf_m_user_data

    implicit none
    save

    
    

    ! COMMON Name: RFLU_user
    ! moved into ppiclf_user_YdotParticle, we need a separate copy of these for every 
    ! call of YdotParticle for GPU
    
    !
    ! For misc values
    !
    ! COMMON Name: ppiclf_misc01
    real*8, parameter ::  OneThird = 1.0d0/3.0d0
    real*8, parameter :: rpi        = acos(-1.0d0)

 

    !
    ! For ppiclf_user_Fluctuations.f
    !
    ! COMMON Name: user_flcut01, user_fluct02, user_fluct03
    


    !
    ! For ppiclf_user_debug.f
    !
    ! COMMON Name: user_debug
    real*8 phimax,                                      &
        fqsx_max,fqsy_max,fqsz_max,                     &
        famx_max,famy_max,famz_max,                     &
        fdpdx_max,fdpdy_max,fdpdz_max,                  &
        fcx_max,fcy_max,fcz_max,                        &
        umean_max,vmean_max,wmean_max,                  &
        fqs_mag,fam_mag,fdp_mag,fc_mag,                 &
        fqsx_fluct_max,fqsy_fluct_max,fqsz_fluct_max,   &
        fqsx_total_max,fqsy_total_max,fqsz_total_max,   &
        fvux_max,fvuy_max,fvuz_max,                     &
        qq_max,tau_max,lift_max


    !
    ! For ppiclf_user_AddedMass.f
    !
    ! COMMON Name: user_AddedMass01, user_AddedMass02

    !
    ! For ppiclf_solve_InitAngularPeriodic
    !
        
    ! MOVED TO ppiclf_data.f90

    ! From ppiclf_user_SetYdot
    INTEGER*4, allocatable, target :: SBin_map(:,:)
    INTEGER*4, allocatable, target :: SBin_counter(:)
    INTEGER*4 i_Bin(3), n_SBin(3), tot_SBin

    ! Finite Diff Material derivative Variables
    integer*4 nstage, istage
    integer*4 icallb
    save      icallb
    data      icallb /0/
    integer*4 idebug
    save      idebug
    data      idebug /0/

    ! Print Data to file
    LOGICAL I_EXIST 
    Character(LEN=25) str 
    integer*4 f_dump
    save      f_dump  
    data      f_dump /1/

    logical exist_file

end module ppiclf_m_user_data