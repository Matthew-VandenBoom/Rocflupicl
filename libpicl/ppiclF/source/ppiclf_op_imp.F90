#include <PPICLF_STD.h>
submodule (ppiclf_op) ppiclf_op_imp
    use mpi
    ! comm variables
    use ppiclf_data, only: ppiclf_comm, ppiclf_nid, ppiclf_np
    implicit none
    
    contains

    ! group op ? gather op ?
    module procedure ppiclf_gop
        !
        ! Input:
        !
        ! real*8 x(n), w(n)
        ! character*3 op
        ! integer*4 n
        !
        ! Internal:
        !
        integer*4 i, ie
        !
        if (op.eq.'+  ') then
            call mpi_allreduce(x,w,n,MPI_DOUBLE_PRECISION,mpi_sum,ppiclf_comm,ie)
        elseif (op.EQ.'M  ') then
            call mpi_allreduce(x,w,n,MPI_DOUBLE_PRECISION,mpi_max,ppiclf_comm,ie)
        elseif (op.EQ.'m  ') then
            call mpi_allreduce(x,w,n,MPI_DOUBLE_PRECISION,mpi_min,ppiclf_comm,ie)
        elseif (op.EQ.'*  ') then
            call mpi_allreduce(x,w,n,MPI_DOUBLE_PRECISION,mpi_prod,ppiclf_comm,ie)
        endif

        do i=1,n
            x(i) = w(i)
        enddo

        return
    END procedure ppiclf_gop

    ! integer gop
    module procedure ppiclf_igop

        !
        ! Input:
        !
        ! integer*4 x(n), w(n)
        ! character*3 op
        ! integer*4 n
        !
        ! Internal:
        !
        integer*4 i, ierr
        !
        if     (op.eq.'+  ') then
            call MPI_Allreduce (x,w,n,mpi_integer,mpi_sum ,ppiclf_comm,ierr)
        elseif (op.EQ.'M  ') then
            call MPI_Allreduce (x,w,n,mpi_integer,mpi_max ,ppiclf_comm,ierr)
        elseif (op.EQ.'m  ') then
            call MPI_Allreduce (x,w,n,mpi_integer,mpi_min ,ppiclf_comm,ierr)
        elseif (op.EQ.'*  ') then
            call MPI_Allreduce (x,w,n,mpi_integer,mpi_prod,ppiclf_comm,ierr)
        endif

        do i=1,n
            x(i) = w(i)
        enddo

        return
    END procedure ppiclf_igop

    ! this method of having 
    module procedure ppiclf_iglsum
        ! 
        ! Input:
        ! 
        ! integer*4 a(n)
        ! integer*4 n
        ! 
        ! Internal:
        ! 
        integer*4 tsum
        integer*4 tmp(n),work(n)
        integer*4 i
        !
        tsum= 0
        do i=1,n
            tsum=tsum+a(i)
        enddo
        tmp(1)=tsum
        call ppiclf_igop(tmp,work,'+  ',1)
        ppiclf_iglsum=tmp(1)
        return
    END procedure ppiclf_iglsum

    module procedure ppiclf_glsum
        ! 
        ! Input:
        ! 
        ! TLJ changed dimension of A
        ! real*8 x(n)
        ! integer*4 n
        !
        ! Internal: 
        !
        real*8 tmp(1),work(1),tsum
        integer*4 i
        TSUM = 0.0d0
        DO I=1,N
            TSUM = TSUM+X(I)
        END DO
        TMP(1)=TSUM
        CALL ppiclf_GOP(TMP,WORK,'+  ',1)
        ppiclf_GLSUM = TMP(1)
        return
    END procedure ppiclf_glsum

    module procedure ppiclf_glmax
        ! 
        ! Input:
        ! 
        ! TLJ changed dimension of A
        ! REAL*8 A(n)
        ! integer*4 n
        !
        ! Internal:
        !
        integer*4 i
        real*8 TMP(1),WORK(1),tmax
        
        TMAX=-99.0e20
        DO I=1,N
            TMAX=MAX(TMAX,A(I))
        END DO
        TMP(1)=TMAX
        CALL ppiclf_GOP(TMP,WORK,'M  ',1)
        ppiclf_GLMAX=TMP(1)
        return
    END procedure ppiclf_glmax

    module procedure ppiclf_iglmax
        ! 
        ! Input:
        ! 
        ! TLJ changed dimension of A
        ! integer*4 a(n)
        ! integer*4 n
        !
        ! Internal:
        !
        integer*4 tmp(1),work(1),tmax
        integer*4 i
        tmax= -999999999
        do i=1,n
            tmax=max(tmax,a(i))
        enddo
        tmp(1)=tmax
        call ppiclf_igop(tmp,work,'M  ',1)
        ppiclf_iglmax=tmp(1)
        return
    END procedure ppiclf_iglmax

    module procedure ppiclf_glmin

        ! 
        ! Input:
        ! 
        ! TLJ changed dimension of A
        ! REAL*8 A(n)
        ! integer*4 n
        !
        ! Internal:
        !
        integer*4 i
        real*8 TMP(1),WORK(1),tmin
        TMIN=99.0e20
        DO I=1,N
            TMIN=MIN(TMIN,A(I))
        END DO
        TMP(1)=TMIN
        CALL ppiclf_GOP(TMP,WORK,'m  ',1)
        ppiclf_GLMIN = TMP(1)
        return
    END procedure ppiclf_glmin

    module procedure ppiclf_iglmin
        ! 
        ! Input:
        ! 
        ! TLJ changed dimension of a
        ! integer*4 a(n)
        ! integer*4 n
        !
        ! Internal:
        !
        integer*4 i
        integer*4 tmp(1),work(1),tmin
        tmin=  999999999
        do i=1,n
            tmin=min(tmin,a(i))
        enddo
        tmp(1)=tmin
        call ppiclf_igop(tmp,work,'m  ',1)
        ppiclf_iglmin=tmp(1)
        return
    END procedure ppiclf_iglmin

    module procedure ppiclf_vlmin
        ! 
        ! Input:
        ! 
        ! TLJ changed dimension of VEC
        ! REAL*8 VEC(n)
        ! integer*4 n
        !
        ! Internal:
        !
        real*8 tmin
        integer*4 i
        TMIN = 99.0E20
        DO I=1,N
            TMIN = MIN(TMIN,VEC(I))
        END DO
        ppiclf_VLMIN = TMIN
        return
    END procedure ppiclf_vlmin

    module procedure ppiclf_vlmax
        ! 
        ! Input:
        ! 
        ! TLJ changed dimension of VEC
        ! REAL*8 VEC(n)
        ! integer*4 n
        !
        ! Internal:
        !
        integer*4 i
        real*8 tmax
        TMAX =-99.0E20
        do i=1,n
            TMAX = MAX(TMAX,VEC(I))
        enddo
        ppiclf_VLMAX = TMAX
        return
    END procedure ppiclf_vlmax

    module procedure ppiclf_copy
        ! 
        ! Input:
        ! 
        ! TLJ changed dimension of A and B
        ! real*8 a(n),b(n)
        ! integer*4 n
        !
        ! Internal:
        !
        integer*4 i
        do i=1,n
            a(i)=b(i)
        enddo

        return
    END procedure ppiclf_copy

    module procedure ppiclf_icopy
        ! 
        ! Input:
        ! 
        ! TLJ changed dimension of A and B
        ! INTEGER*4 A(n), B(n)
        ! integer*4 n
        !
        ! Internal:
        !
        integer*4 i
        DO I = 1, N
            A(I) = B(I)
        END DO
        return
    END procedure ppiclf_icopy

    module procedure ppiclf_chcopy
        ! 
        ! Input:
        ! 
        ! TLJ changed A and B dimenions
        ! CHARACTER*1 A(n), B(n)
        ! integer*4 n
        !
        ! Internal:
        ! 
        integer*4 i
        DO I = 1, N
            A(I) = B(I)
        END DO
        return
    END procedure ppiclf_chcopy
    
    module procedure ppiclf_exittr
        ! 
        ! Input:
        ! 
        ! character(len=*) stringi
        ! real*8 rdata
        ! integer*4 idata
        !
        ! Internal:
        !
        character(len=fixed_str_len) stringo
        character*25 s25
        integer*4 ilen, ierr, k
        
        call ppiclf_blank(stringo,fixed_str_len)
        call ppiclf_chcopy(stringo,stringi,len(stringi))
        ilen = len(stringi) !ppiclf_indx1(stringo,'$')
        write(s25,25) rdata,idata
        25 format(1x,1p1e14.6,i10)
        call ppiclf_chcopy(stringo(ilen + 1:),s25,25)

        ! changed to always print the error message, not just on rank 0
        write(6,11) ppiclf_nid, stringo(1:ilen+24)
        11 format('PPICLF: ERROR on rank ',i2,":",a)

        !     call mpi_finalize (ierr)
        call mpi_abort(ppiclf_comm, 1, ierr)

        return
    END procedure ppiclf_exittr

    module procedure ppiclf_printsri
        !
        ! Input:
        !
        ! character(len=*) stringi
        ! real*8 rdata
        ! integer*4 idata
        !
        ! Internal:
        !
        character(len=fixed_str_len) stringo
        character*25 s25
        integer*4 ilen, k, ierr

#ifdef TEST
        return
#endif
        call ppiclf_blank(stringo,fixed_str_len)
        call ppiclf_chcopy(stringo,stringi,len(stringi))
        ilen = ppiclf_indx1(stringo,'$')
        write(s25,25) rdata,idata
        25 format(1x,1p1e14.6,i10)
        call ppiclf_chcopy(stringo(ilen:),s25,25)

        call mpi_barrier(ppiclf_comm,ierr)

        if (ppiclf_nid.eq.0) write(6,12) stringo(1:ilen+24)
        12 format('PPICLF: ',a)

        call mpi_barrier(ppiclf_comm,ierr)

        return
    END procedure ppiclf_printsri

    module procedure ppiclf_printsi
        !
        ! Input:
        !
        ! character(len=*) stringi
        ! integer*4 idata
        !
        ! Internal:
        !
        character(len=fixed_str_len) stringo
        character*10 s10
        integer*4 ilen, k, ierr

#ifdef TEST
        return
#endif
        call ppiclf_blank(stringo,fixed_str_len)
        call ppiclf_chcopy(stringo,stringi,len(stringi))
        ilen = ppiclf_indx1(stringo,'$')
        write(s10,10) idata
        10 format(1x,i9)
        call ppiclf_chcopy(stringo(ilen:),s10,10)

        call mpi_barrier(ppiclf_comm,ierr)

        if (ppiclf_nid.eq.0) write(6,13) stringo(1:ilen+9)
        13 format('PPICLF: ',a)

        call mpi_barrier(ppiclf_comm,ierr)

        return
    END procedure ppiclf_printsi

    module procedure ppiclf_printsr
        !
        ! Input:
        !
        ! character(len=*) stringi
        ! real*8 rdata
        !
        ! Internal:
        !
        character(len=fixed_str_len) stringo
        character*15 s15
        integer*4 ilen, k, ierr

#ifdef TEST
        return
#endif
        call ppiclf_blank(stringo,fixed_str_len)
        call ppiclf_chcopy(stringo,stringi,len(stringi))
        ilen = ppiclf_indx1(stringo,'$')
        write(s15,15) rdata
        15 format(1x,1p1e14.6)
        call ppiclf_chcopy(stringo(ilen:),s15,15)

        call mpi_barrier(ppiclf_comm,ierr)

        if (ppiclf_nid.eq.0) write(6,14) stringo(1:ilen+14)
        14 format('PPICLF: ',a)

        call mpi_barrier(ppiclf_comm,ierr)

        return
    END procedure ppiclf_printsr

    module procedure ppiclf_prints
        !
        ! Input:
        !
        ! character(len=*) stringi
        !
        ! Internal:
        !
        character(len=fixed_str_len) stringo
        integer*4 ilen, k, ierr

#ifdef TEST
        return
#endif
        call ppiclf_blank(stringo,fixed_str_len)
        call ppiclf_chcopy(stringo,stringi,len(stringi))
        ilen = ppiclf_indx1(stringo,'$')

        call mpi_barrier(ppiclf_comm,ierr)

        if (ppiclf_nid.eq.0) write(6,21) stringo(1:ilen-1)
        21 format('PPICLF: ',a)

        call mpi_barrier(ppiclf_comm,ierr)

        return
    END procedure ppiclf_prints

    module procedure PPICLF_BLANK
        ! 
        ! Input:
        !
        ! TLJ changed dimension of A
        ! CHARACTER(len=*), intent(out) :: A
        ! integer*4 n
        !
        ! Internal:
        !
        CHARACTER*1 BLNK
        SAVE        BLNK
        DATA        BLNK /' '/
        integer*4 i
        !
        !
        DO I=1,N
            A(I:I)=BLNK
        END DO
        RETURN
    END procedure PPICLF_BLANK

    ! finds the first index of a character(S2) in a string(S1)
    module procedure PPICLF_INDX1
        ! 
        ! Input:
        !
        ! CHARACTER(len=*) S1
        ! character(len=1) S2
        !
        ! Internal:
        !
        integer*4 n1, i
        !
        N1 = LEN(S1)
        PPICLF_INDX1=0
        ! IF (N1.LT.1) return
        !
        DO I=1,N1
            IF (S1(I:I).EQ.S2) THEN
                PPICLF_INDX1=I
                return
            END IF
        END DO
        !
        return
    END procedure PPICLF_INDX1

    module procedure PPICLF_CHSTR
        ! 
        ! Input:
        !
        ! TLJ modified, but not sure why I had to
        ! CHARACTER S1
        ! INTEGER indx1
        !
        PPICLF_CHSTR = S1(1:indx1)
        return
    END procedure PPICLF_CHSTR

    module procedure ppiclf_byte_open_mpi
        !
        ! Input:
        !
        ! character fnamei*(*)
        ! logical ifro
        ! integer*4 ierr, mpi_fh
        ! 
        ! Internal:
        !
        CHARACTER*1 BLNK
        DATA BLNK/' '/
        character*132 fname
        character*1   fname1(132)
        equivalence  (fname1,fname)
        integer*4 imode
        !
        imode = MPI_MODE_WRONLY+MPI_MODE_CREATE
        if(ifro) then
            imode = MPI_MODE_RDONLY 
        endif

        call MPI_file_open(ppiclf_comm,fnamei,imode, MPI_INFO_NULL,mpi_fh,ierr)

        return
    END procedure ppiclf_byte_open_mpi

    module procedure ppiclf_byte_read_mpi
        !
        ! Input:
        !
        ! real*4 buf(1)          ! buffer
        ! integer*4 icount, mpi_fh, ierr
        
        !
        ! Internal:
        !
        integer*4 iout

        iout = icount ! icount is in 4-byte words
        call MPI_file_read_all(mpi_fh,buf,iout,MPI_REAL, MPI_STATUS_IGNORE,ierr)

        return
    END procedure ppiclf_byte_read_mpi

    module procedure ppiclf_byte_write_mpi
        !
        ! Input:
        !
        ! real*4 buf(1)          ! buffer
        ! integer*4 icount, iorank, mpi_fh, ierr
        !
        ! Internal:
        !
        integer*4 iout

        iout = icount ! icount is in 4-byte words
        if(iorank.ge.0 .and. ppiclf_nid.ne.iorank) iout = 0
        call MPI_file_write_all(mpi_fh,buf,iout,MPI_REAL, MPI_STATUS_IGNORE,ierr)

        return
    END procedure ppiclf_byte_write_mpi

    module procedure ppiclf_byte_close_mpi
        !
        ! Input:
        !
        ! integer*4 mpi_fh, ierr
        !

        call MPI_file_close(mpi_fh,ierr)

        return
    END procedure ppiclf_byte_close_mpi

    module procedure ppiclf_byte_set_view
        !
        ! Input:
        !
        ! integer*8 ioff_in
        ! integer*4 mpi_fh
        !
        ! Internal:
        !
        integer*4 ierr
        call MPI_file_set_view(mpi_fh,ioff_in,MPI_BYTE,MPI_BYTE,'native',MPI_INFO_NULL,ierr)

        return
    END procedure ppiclf_byte_set_view

    module procedure ppiclf_bcast
        !
        ! Input:
        !
        ! type(*) buf(*)
        ! integer*4 len, 
        
        !
        ! Internal:
        !
        integer*4 ierr

        call mpi_bcast (buf,len,mpi_byte,0,ppiclf_comm,ierr)

        return
    END procedure ppiclf_bcast

    module procedure ppiclf_CELL_MAP_GroupByDestRank
        !
        ! Internal:
        !
        ! integer*4, intent(in)    :: nMappedCells, nranks
        ! INTEGER*4, intent(inout) :: cell_map(PPICLF_LRMAX,PPICLF_LEE), cell_map_counts(0:ppiclf_np), cell_map_disps(0:ppiclf_np)
        !
        ! Internal:
        !
        integer*4 :: temp_cell_map(PPICLF_LRMAX, PPICLF_LEE)
        integer*4 :: rankInfo(0:ppiclf_np)
        integer*4 i, cellRank

        cell_map_counts = 0

        do i = 1, nMappedCells
            cell_map_counts(cell_map(3, i)) = cell_map_counts(cell_map(3, i)) + 1
        end do

        cell_map_disps(0) = 0
        do i = 1, ppiclf_np
            ! this start index = prev start index + n prev
            cell_map_disps(i)  = cell_map_disps(i-1) + cell_map_counts(i-1)
        end do

        rankInfo(:) = cell_map_disps + 1
        do i = 1, nMappedCells
            cellRank = cell_map(3, i)
            temp_cell_map(:, rankInfo(cellRank)) = cell_map(:, i)
            rankInfo(cellRank) = rankInfo(cellRank) + 1
        end do
        cell_map = temp_cell_map
    end procedure ppiclf_CELL_MAP_GroupByDestRank

end submodule ppiclf_op_imp