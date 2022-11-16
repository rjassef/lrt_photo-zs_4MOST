cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Program to calculate photometric redshifts for testing the LRT
c     libraries.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program main
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=40,NGMAX=17000,NSMAX=4,NWMAX=350,NXMAX=5000)
      parameter (NEMAX=37,NIMAX=8)

      integer*8 id_obj

      real*8 jy(NCMAX),ejy(NCMAX),comp(4)
      integer jyuse(NCMAX)

      character*200 line
      character*200 fname_in, fname_out
      character*200 zmin_str, zmax_str, dz_str
      real*8 zmin, zmax, dz

      real*8 zgal(NXMAX),galjy(NXMAX,NSMAX,NCMAX,NEMAX,NIMAX)
      common /galtemp/zgal,galjy,ngalz,ngalt

c     Open the files used.
      if(iargc().ne.5) then
        print*,"Correct use: ./zphot_noagn fname_in fname_out zmin zmax dz"
        stop
      endif
      call getarg(1,fname_in)
      call getarg(2,fname_out)
      call getarg(3,zmin_str)
      call getarg(4,zmax_str)
      call getarg(5,dz_str)
      open(unit=11, file=trim(fname_in), status='old')
      open(unit=12, file=trim(fname_out), status='unknown')

c     Initialize the subroutines.
      read(zmin_str,*) zmin
      read(zmax_str,*) zmax
      read(dz_str  ,*) dz
      call pzinit('bandmag.dat',1,0,0,zmin,zmax,dz,0)
c      call pzinit('bandmag.dat',1,0,0,0.d0,6.0d0,0.01d0,0)

 200  read(11,'(a)')line
         if(line(1:1).eq.'#') goto 200
      continue
      read(line,*)ntarg,nchan

c     Do not use the AGN template.
      do i=1,ngalz
        do j=1,nchan
          galjy(i,1,j,1,1) = 0.d0
        enddo
      enddo

c      print*
c      print*,'Calculating photometric redshifts for ',ntarg,' objects'

c     Check if we are restarting from a mid point or not. 
      call restart(ndone, 11, 12, fname_out)

c     Start the main cycle.
      do i=ndone+1,ntarg

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
