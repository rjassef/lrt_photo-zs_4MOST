        program main
        implicit real*8 (a-h,o-z)
        parameter (NCMAX=40)

        integer*8 id_obj

        real*8 jy(NCMAX),ejy(NCMAX),comp(2)
        integer jyuse(NCMAX)

        character*200 line
        character*200 stype_str, fname_in, fname_out

c     Open the files used.
        if(iargc().ne.3) then
            print*,"Correct use: ./star_fit stype fname_in fname_out"
            stop
        endif
        call getarg(1, stype_str)
        read(stype_str,*)istype
        call getarg(2,fname_in)
        call getarg(3,fname_out)
        open(unit=11, file=trim(fname_in), status='old')
        open(unit=12, file=trim(fname_out), status='unknown')

c     Initialize the subroutines.
        call starinit('bandmag.dat',istype,0)

200     read(11,'(a)')line
            if(line(1:1).eq.'#') goto 200
        continue
        read(line,*)ntarg,nchan

c     Check if we are restarting from a mid point or not. 
        call restart(ndone, 11, 12, fname_out)

c     Start the main cycle.
        do i=ndone+1,ntarg

c     Read the 6'' aperture Photometry File
            read(11,*)id_obj, z, (jy(j),j=1,nchan),
     *        (ejy(j),j=1,nchan), (jyuse(j),j=1,nchan)

c     Calculate number of usable bands
            m = 0
            do j = 1,nchan
                if(jyuse(j).gt.0) m = m + 1
            enddo
            if(m.lt.2) then
                ns_best = -1
                chi2 = -1.d0
                do l=1,2
                    comp(l) = 0.d0
                enddo
                goto 100
            endif

c     Calculate the stellar fits.
            call stf(jy,ejy,jyuse,jymod,comp,ns_best,chi2,0)

100         continue

            write(12,120)id_obj,m,ns_best,chi2,(comp(l),l=1,2)
120             format(3i20,25E20.6)

        enddo

        close(11)
        close(12)

        stop
        end
