        subroutine restart(fname, ndone)
        implicit real*8(a-h,o-z)
        logical file_exists

c       First, check if the file exists.
        inquire(unit=12, opened=file_exists)

        ndone = 0
        if(file_exists) then
c       If the file exists, count how many rows have been processed.
            do
                read(12,*,end=100)
                ndone = ndone + 1
            enddo 
100         continue
        endif
        close(12)

        return 


        