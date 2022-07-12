cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Program to calculate photometric redshifts for testing the LRT
c     libraries.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program main
      implicit real*8 (a-h,o-z)

      integer*8 id_obj

      real*8 jy(13),ejy(13),comp(4)
      integer jyuse(13)

      character*200 line
      character*200 fname_in, fname_out

c     Open the files used.
      if(iargc().ne.2) then
        print*,"Correct use: ./zphot fname_in fname_out"
        stop
      endif
      call getarg(1,fname_in)
      call getarg(2,fname_out)
      open(unit=11, file=trim(fname_in), status='old')
      open(unit=12, file=trim(fname_out), status='unknown')

c     Initialize the subroutines.
      call pzinit('bandmag.dat',1,1,0,0.d0,6.0d0,0.01d0,0)

 200  read(11,'(a)')line
         if(line(1:1).eq.'#') goto 200
      continue
      read(line,*)ntarg,nchan

c      print*
c      print*,'Calculating photometric redshifts for ',ntarg,' objects'


c     Start the main cycle.
      do i=1,ntarg

c     Read the 6'' aperture Photometry File
         read(11,*)id_obj, z, (jy(j),j=1,nchan),
     *   (ejy(j),j=1,nchan), (jyuse(j),j=1,nchan)

c     Calculate number of usable bands
         m = 0
         do j = 1,nchan
            if(jyuse(j).gt.0) m = m + 1
         enddo
         if(m.lt.4) then
           zp = -1.d0
           chigal = -1.d0
           chinop = -1.d0
           goto 100
         endif

c     Calculate the photometric redshifts.
         call pza(jy,ejy,jyuse,zp,chigal,chinop,0,0)

 100     continue

         write(12,120)id_obj,m,zp,z,chigal,chinop
 120     format(2i20,4E20.6)

      enddo

      close(11)
      close(12)

      stop
      end