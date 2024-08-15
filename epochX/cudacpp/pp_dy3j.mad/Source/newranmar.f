      subroutine ntuple_new(x,a,b,ii,jconfig)
c-------------------------------------------------------
c     Front to ranmar which allows user to easily
c     choose the seed.
c------------------------------------------------------
      implicit none
c
c     Arguments
c
      double precision x,a,b
      integer ii,jconfig
c
c     Local
c
      integer init, ioffset, joffset
      integer     ij, kl, iseed1,iseed2

c
c     Global
c
c-------
c     18/6/2012 tjs promoted to integer*8 to avoid overflow for iseed > 60K
c------
      integer*8       iseed
      common /to_seed/iseed
c
c     Data
c
      data init /1/
      save ij, kl
c-----
c  Begin Code
c-----
      if (init .eq. 1) then
         init = 0
         call get_offset(ioffset)
         if (iseed .eq. 0) call get_base(iseed)
c
c     TJS 3/13/2008
c     Modified to allow for more sequences 
c     iseed can be between 0 and 30081*30081
c     before pattern repeats
c
c
c     TJS 12/3/2010
c     multipied iseed to give larger values more likely to make change
c     get offset for multiple runs of single process
c
c
c     TJS 18/6/2012
c     Updated to better divide iseed among ij and kl seeds
c     Note it may still be possible to get identical ij,kl for
c     different iseed, if have exactly compensating joffset, ioffset, jconfig
c
         call get_moffset(joffset)
         joffset = joffset * 3157
         iseed = iseed * 31300       
         ij=1802+jconfig + mod(iseed,30081)
         kl=9373+(iseed/30081)+ioffset + joffset     !Switched to 30081  20/6/12 to avoid dupes in range 30082-31328
         write(*,'(a,i6,a3,i6)') 'Using random seed offsets',jconfig," : ",ioffset
         write(*,*) ' with seed', iseed/31300
         do while (ij .gt. 31328)
            ij = ij - 31328
         enddo
         do while (kl .gt. 30081)
            kl = kl - 30081
         enddo
        call rmarin(ij,kl)         
      endif
      call ranmar(x)
      do while (x .lt. 1d-16)
         call ranmar(x)
      enddo
      x = a+x*(b-a)
      end
