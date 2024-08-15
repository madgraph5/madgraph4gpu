      subroutine ntuple_new(x)
c-------------------------------------------------------
c     Front to ranmar, assuming seeds have already been set
c------------------------------------------------------
      implicit none
c
c     Arguments
c
      double precision x
c-----
c  Begin Code
c-----
      call ranmar(x)
      do while (x .lt. 1d-16)
         call ranmar(x)
      enddo
      end
