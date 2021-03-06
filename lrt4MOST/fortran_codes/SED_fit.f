cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Program to calculate K corrections redshifts for testing the LRT
c     libraries.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      program main
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=40,NSMAX=4)

      integer*8 id_obj

      real*8 vec(NSMAX)
      real*8 jymod(NSMAX,NCMAX)
      real*8 jymodtot(NCMAX)
      common /models/jymod,jymodtot,vec

      real*8 jy(NCMAX),ejy(NCMAX)
      integer jyuse(NCMAX)
      real*8 jymodel(NCMAX),jycorr(NCMAX),comp(4)
      real*8 mag_kcorr(NCMAX),mag(NCMAX)

      real*8 jyzero(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,con,lbar

      real*8 ebv,igm
      real*8 cov(4,4)

      character*200 line
      character*200 fname_in, fname_out

c     Open the files used.
      if(iargc().ne.2) then
        print*,"Correct use: ./SED_fit fname_in fname_out"
        stop
      endif
      call getarg(1,fname_in)
      call getarg(2,fname_out)
      open(unit=11, file=trim(fname_in), status='old')
      open(unit=12, file=trim(fname_out), status='unknown')

c     Initialize the subroutines.
      call kcinit('bandmag.dat',1,1,1,0)

 200  read(11,'(a)')line
         if(line(1:1).eq.'#') goto 200
      continue
      read(line,*)ntarg,nchan

c     Check if we are restarting from a mid point or not. 
      call restart(ndone, 11, 12, fname_out)

c     Start the main cycle.
      do i=ndone+1,ntarg

c     Read the 6'' aperture Photometry File
         read(11,*)id_obj, z, (jy(j),j=1,nchan),
     *     (ejy(j),j=1,nchan), (jyuse(j),j=1,nchan)

c     Calculate number of usable bands
         m = 0
         do j = 1,nchan
            if(jyuse(j).gt.0) m = m + 1
         enddo
         if(m.lt.4) then
           chi2 = -1.d0
           ahat = -1.d0
           ebv  = -1.d0
           igm  = -1.d0
           do l=1,4
             comp(l) = 0.d0
             vec(l) = 0.d0
           enddo
           goto 100
        endif

c     Calculate the photometric redshifts
         z0 = 0.d0

         call kca(jy,ejy,jyuse,z,z0,jymodel,jycorr,comp,
     *        cov,ebv,igm,chi2,0)

         ahat = comp(1)/(comp(1)+comp(2)+comp(3)+comp(4))

 100     continue

         write(12,120)id_obj,m,z,chi2,ahat,ebv,igm,
     *   (comp(l),l=1,4),(vec(l),l=1,4)
 120     format(2i20,25E20.6)

      enddo

      close(11)
      close(12)

      stop
      end
