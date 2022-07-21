        program main
        implicit real*8(a-h,o-z)
        parameter (NCMAX=32)

        real*8 corr(NCMAX)
        character*200 fname_in

c     Open the files used.
        if(iargc().ne.1) then
                print*,"Correct use: ./run_fitzero fname_in"
                stop
        endif
        call getarg(1,fname_in)

c     Initialize the subroutines.
        call kcinit('bandmag.dat',1,1,1,1)

        call fitzero_omp(fname_in,10,0.8d0,corr,0)

        nchan = 13
        open(unit=50,file="channel.zpc",status='unknown')
        do j=1,nchan
                write(50,130)corr(j)
        enddo
        close(50)
130     format(1ES20.6)

        stop
        end
