#include <PPICLF_STD.h>
module ppiclf_op
    implicit none
    private

    public :: ppiclf_gop
    public :: ppiclf_igop
    public :: ppiclf_iglsum
    public :: ppiclf_glsum
    public :: ppiclf_glmax
    public :: ppiclf_iglmax
    public :: ppiclf_glmin
    public :: ppiclf_iglmin
    public :: ppiclf_vlmin
    public :: ppiclf_vlmax
    public :: ppiclf_copy
    public :: ppiclf_icopy
    public :: ppiclf_chcopy
    public :: ppiclf_exittr
    public :: ppiclf_printsri
    public :: ppiclf_printsi
    public :: ppiclf_printsr
    public :: ppiclf_prints
    public :: PPICLF_BLANK
    public :: PPICLF_INDX1
    public :: PPICLF_CHSTR
    public :: ppiclf_byte_open_mpi
    public :: ppiclf_byte_read_mpi
    public :: ppiclf_byte_write_mpi
    public :: ppiclf_byte_close_mpi
    public :: ppiclf_byte_set_view
    public :: ppiclf_bcast
    public :: ppiclf_CELL_MAP_GroupByDestRank

    public fixed_str_len
    ! internal module variables
    integer, parameter :: fixed_str_len = 132
    interface

        ! group op ? gather op ?
        module SUBROUTINE ppiclf_gop( x, w, op, n)
            real*8 x(n), w(n)
            character*3 op
            integer*4 n
        END SUBROUTINE ppiclf_gop

        ! integer gop
        module SUBROUTINE ppiclf_igop( x, w, op, n)
            integer*4 x(n), w(n)
            character*3 op
            integer*4 n
        END SUBROUTINE ppiclf_igop

        integer*4 module FUNCTION ppiclf_iglsum(a,n)

            integer*4 a(n)
            integer*4 n
        END FUNCTION ppiclf_iglsum

        real*8 module FUNCTION ppiclf_glsum (x,n)
            real*8 x(n)
            integer*4 n
        END FUNCTION ppiclf_glsum

        real*8 module FUNCTION ppiclf_glmax(a,n)
            REAL*8 A(n)
            integer*4 n
        END FUNCTION ppiclf_glmax

        integer*4 module FUNCTION ppiclf_iglmax(a,n)
            integer*4 a(n)
            integer*4 n
        END FUNCTION ppiclf_iglmax

        real*8 module FUNCTION ppiclf_glmin(a,n)
            REAL*8 A(n)
            integer*4 n
        END FUNCTION ppiclf_glmin

        integer*4 module FUNCTION ppiclf_iglmin(a,n)
            integer*4 a(n)
            integer*4 n
        END FUNCTION ppiclf_iglmin

        real*8 module FUNCTION ppiclf_vlmin(vec,n)
            REAL*8 VEC(n)
            integer*4 n
        END FUNCTION ppiclf_vlmin

        real*8 module FUNCTION ppiclf_vlmax(vec,n)
            REAL*8 VEC(n)
            integer*4 n
        END FUNCTION ppiclf_vlmax

        module SUBROUTINE ppiclf_copy(a,b,n)
            real*8 a(n),b(n)
            integer*4 n
        END SUBROUTINE ppiclf_copy

        module SUBROUTINE ppiclf_icopy(a,b,n)
            INTEGER*4 A(n), B(n)
            integer*4 n
        END SUBROUTINE ppiclf_icopy

        module SUBROUTINE ppiclf_chcopy(a,b,n)
            CHARACTER*1 A(n), B(n)
            integer*4 n
        END SUBROUTINE ppiclf_chcopy
        
        module SUBROUTINE ppiclf_exittr(stringi,rdata,idata)
            character(len=*) stringi
            real*8 rdata
            integer*4 idata
        END SUBROUTINE ppiclf_exittr

        module SUBROUTINE ppiclf_printsri(stringi,rdata,idata)
            character(len=*) stringi
            real*8 rdata
            integer*4 idata
        END SUBROUTINE ppiclf_printsri

        module SUBROUTINE ppiclf_printsi(stringi,idata)
            character(len=*) stringi
            integer*4 idata
        END SUBROUTINE ppiclf_printsi

        module SUBROUTINE ppiclf_printsr(stringi,rdata)
            character(len=*) stringi
            real*8 rdata
        END SUBROUTINE ppiclf_printsr

        module SUBROUTINE ppiclf_prints(stringi)
            character(len=*) stringi
        END SUBROUTINE ppiclf_prints

        module SUBROUTINE PPICLF_BLANK(A,N)
            CHARACTER(len=*), intent(out) :: A
            integer*4 n
        END SUBROUTINE PPICLF_BLANK

        ! finds the first index of a character(S2) in a string(S1)
        INTEGER*4 module FUNCTION PPICLF_INDX1(S1,S2)
            CHARACTER(len=*) S1
            character(len=1) S2
        END FUNCTION PPICLF_INDX1

        character*132 module FUNCTION PPICLF_CHSTR(S1,indx1)
            CHARACTER S1
            INTEGER indx1
        END FUNCTION PPICLF_CHSTR

        module SUBROUTINE ppiclf_byte_open_mpi(fnamei,mpi_fh,ifro,ierr)
            character fnamei*(*)
            logical ifro
            integer*4 ierr, mpi_fh
        END SUBROUTINE ppiclf_byte_open_mpi

        module SUBROUTINE ppiclf_byte_read_mpi(buf,icount,mpi_fh,ierr)
            real*4 buf(1)          ! buffer
            integer*4 icount, mpi_fh, ierr
        END SUBROUTINE ppiclf_byte_read_mpi

        module SUBROUTINE ppiclf_byte_write_mpi(buf,icount,iorank,mpi_fh,ierr)
            real*4 buf(1)          ! buffer
            integer*4 icount, iorank, mpi_fh, ierr
        END SUBROUTINE ppiclf_byte_write_mpi

        module SUBROUTINE ppiclf_byte_close_mpi(mpi_fh,ierr)
            integer*4 mpi_fh, ierr
        END SUBROUTINE ppiclf_byte_close_mpi

        module SUBROUTINE ppiclf_byte_set_view(ioff_in,mpi_fh)
            integer*8 ioff_in
            integer*4 mpi_fh
        END SUBROUTINE ppiclf_byte_set_view

        module SUBROUTINE ppiclf_bcast(buf,len)
            type(*) buf(*)
            integer*4 len
        END SUBROUTINE ppiclf_bcast

        module pure subroutine ppiclf_CELL_MAP_GroupByDestRank(nMappedCells, cell_map, cell_map_counts, cell_map_disps, nranks)
            integer*4, intent(in)    :: nMappedCells, nranks
            INTEGER*4, intent(inout) :: cell_map(PPICLF_LRMAX,PPICLF_LEE), cell_map_counts(0:nranks-1), cell_map_disps(0:nranks-1)
        end subroutine ppiclf_CELL_MAP_GroupByDestRank
    end interface
end module ppiclf_op
