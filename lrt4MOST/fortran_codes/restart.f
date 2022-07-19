        subroutine restart(ndone, fp_in, fp_out, fname_out)
        implicit real*8(a-h,o-z)
        integer fp_in, fp_out
        character*200 fname_out
        logical file_exists

c       First, check if the file exists.
        inquire(unit=fp_out, opened=file_exists)

        ndone = 0
        if(file_exists) then
c       If the file exists, count how many rows have been processed.
            do
                read(fp_out,*,end=100)
                ndone = ndone + 1
            enddo 
100         continue
        endif
        close(fp_out)

c     Close and reopen in append mode.
        close(fp_out)
        open(unit=fp_out, file=trim(fname_out), status='unknown', access='append')

c     Skip the ones already done. 
        do i=1,ndone
            read(fp_in,*)
        enddo

        return 
        end

        