#ifndef GENERIC_TYPE_A_extention
#error "GENERIC_TYPE_A_extention was not defined"
#endif
#ifndef GENERIC_TYPE_A_fullname
#error "GENERIC_TYPE_A_fullname was not defined"
#endif
subroutine GENERIC_TYPE_A_extention(ppiclf_alltoallv)(sendbuf, sendcounts, senddispls, recvbuf, recvcounts, recvdispls, type_handle, comm, maxrecvtotal)
    ! this function fills a similar role to mpi_alltoallv, but allows for unknown receive counts/displacements.
    ! The only restriction this function enforces that mpi_alltoallv doesn't is both the send and recieved data is the same type.

    ! If the first element of recvcounts is < 0, an extra alltoall communication will be initiated, using the 
    ! sendcounts to update every rank's recvcounts and recvdispls.
    ! If the first element of recvcounts is >= 0, it and recvdispls will be used as is for the transfer.

    ! NOTE: if any rank has recvcounts < 0, all ranks must, otherwise everything will deadlock. This is not checked by
    ! this function for efficiency
    ! maxrecvtotal is the maximum number of elements that can be received from all ranks. Typically dictated by the size of recvbuf.
    !
    ! Inputs:
    !
    type(GENERIC_TYPE_A_fullname), intent(in)     :: sendbuf(:)
    type(GENERIC_TYPE_A_fullname), intent(out)    :: recvbuf(:)
    integer, intent(in)     :: sendcounts(0:ppiclf_np-1), senddispls(0:ppiclf_np-1)
    integer, intent(inout)  :: recvcounts(0:ppiclf_np-1), recvdispls(0:ppiclf_np-1)
    integer, intent(in)     :: type_handle, comm, maxrecvtotal
    !
    ! Internal:
    !
    integer :: ialltoall_request
    integer :: send_request_handles(0:ppiclf_np-1), recv_request_handles(0:ppiclf_np-1)
    integer :: n_sends, n_recvs
    integer stat_local(MPI_STATUS_SIZE)
    integer stats(MPI_STATUS_SIZE, 0:ppiclf_np-1)
    logical :: syncing_counts
    integer :: i, ierr
    !
    ! Code:
    !
    syncing_counts = .false.
    ! if needed, start the transfer of the send counts
    if (recvcounts(0) .lt. 0) then
        syncing_counts = .true.
        ! if (rank .eq. 0) then
        !     print*, "need to sync send/recv counts"
        ! end if
        call mpi_ialltoall(sendcounts, 1, MPI_INTEGER4, recvcounts, 1, MPI_INTEGER4, MPI_COMM_WORLD, ialltoall_request, ierr)
    end if ! recvcounts(0) < 0

    ! initiate the sending of data
    n_sends = 0
    do i = 0, ppiclf_np - 1
        ! don't send if we dont need to
        if (sendcounts(i) .eq. 0) cycle

        call MPI_ISEND(sendbuf(senddispls(i) + 1), sendcounts(i), type_handle, i, PPICLF_MPI_tags%ppiclf_alltoallv_transfer, comm, send_request_handles(n_sends), ierr)
        n_sends = n_sends + 1
    end do ! data send

    ! recieve the transfer of the send counts from other ranks
    if (syncing_counts) then
        ! if (rank .eq. 0) then
        !     print*, "about to wait on completion of ialltoall"
        ! end if
        call MPI_WAIT(ialltoall_request, stat_local, ierr)
        ! if (rank .eq. 0) then
        !     print*, "recieved send counts from other ranks"
        ! end if
        ! recompute displacements
        recvdispls(0) = 0
        do i = 1, ppiclf_np - 1
            recvdispls(i) = recvdispls(i-1) + recvcounts(i-1)
        end do
        
        ! verify that we dont intend to recieve too much data
        if (sum(recvcounts) .gt. maxrecvtotal) then
            call ppiclf_exittr("tried to recieve too many values in ppiclf_alltoallv", real(sum(recvcounts), kind(0.0d0)), maxrecvtotal)
        end if
    end if ! recvcounts(0) < 0

    ! initiate the receiving of data
    n_recvs = 0
    do i = 0, ppiclf_np - 1
        ! don't receive if there is nothing being sent
        if (recvcounts(i) .eq. 0) cycle

        call MPI_IRECV(recvbuf(recvdispls(i) + 1), recvcounts(i), type_handle, i, PPICLF_MPI_tags%ppiclf_alltoallv_transfer, comm, recv_request_handles(n_recvs), ierr)
        n_recvs = n_recvs + 1
    end do ! data send
    call MPI_WAITALL(n_sends, send_request_handles, stats, ierr)
    call MPI_WAITALL(n_recvs, recv_request_handles, stats, ierr)
end subroutine GENERIC_TYPE_A_extention(ppiclf_alltoallv)