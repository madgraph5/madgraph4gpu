c-------------------------------------------------------------------------------

      subroutine sample_get_x(wgt, x, j, ipole, xmin, xmax)
c************************************************************************
c     Returns maxdim random numbers between 0 and 1, and the wgt
c     associated with this set of points, and the iteration number
c     This routine chooses the point within the range specified by
c     xmin and xmax for dimension j in configuration ipole
c************************************************************************
c     INPUT  arguments: wgt, j, ipole, xmin, xmax
c     OUTPUT arguments: wgt, x
c************************************************************************
c     INPUT  common: ituple (only ituple=1 is supported)
c     INPUT  common: minvar
c     INPUT  common: grid
c     INPUT  common: tx, nzoom
c     INPUT  common: swidth, spole
c     OUTPUT common: lastbin
c************************************************************************
      implicit none
c
c     Constants
c
      include 'genps.inc' ! relevant parameters: ng, maxdim, maxinvar
      include 'maxconfigs.inc' ! relevant parameters: lmaxconfigs
c
c     Arguments
c
      double precision wgt, x, xmin, xmax
      integer j, ipole
c
c     Local
c
      integer im, ip, ij
      double precision xbin_min, xbin_max, ddum, xo
      data ddum/0d0/

      integer xbinarraydim
      parameter (xbinarraydim=maxdim*lmaxconfigs)
      double precision xbin_min0_array(maxdim, lmaxconfigs)
      double precision xbin_max1_array(maxdim, lmaxconfigs)
      logical xbin_min0_saved(maxdim, lmaxconfigs)
      logical xbin_max1_saved(maxdim, lmaxconfigs)
      save xbin_min0_array, xbin_max1_array
      save xbin_min0_saved, xbin_max1_saved
      data xbin_min0_saved/xbinarraydim*.false./
      data xbin_max1_saved/xbinarraydim*.false./
      
      integer icall
      save icall
      data icall/0/
c
c     External
c
      double precision xbin
      external xbin
c
c     Global
c
      double precision grid(2, ng, 0:maxinvar)
      common /data_grid/ grid

      integer Minvar(maxdim,lmaxconfigs)
      common /to_invar/ Minvar

      integer ituple
      common /to_random/ituple

      double precision spole(maxinvar), swidth(maxinvar), bwjac
      common/to_brietwigner/spole, swidth, bwjac

      integer nzoom
      double precision tx(1:3, maxinvar)
      common/to_xpoints/tx, nzoom

      integer lastbin(maxdim)
      common /to_lastbin/lastbin

c-----
c  Begin Code
c-----
c     CALL COUNTERS_START_COUNTER( 10, 1 ) ! 10=PROGRAM-SampleGetX
      icall = icall + 1
c     if (icall.le.100) then
c       write(6,*) 'sample_get_x', wgt, x, j, ipole, xmin, xmax
c     endif
      ij = Minvar(j, ipole)
c
c Fall back to the default old implementation under some circumstances
c - first call  => take care of ranmar initialization in the old ntuple
c - swidth > 0  => take care of calling transpole
c - ituple != 1 => take care of error message      
c - nzoom > 0   => take care of updating tx
c
      if (icall.eq.1) then ! first call
        call sample_get_x_old( wgt, x, j, ipole, xmin, xmax )
        goto 999
      endif
      if (ij.gt.0 .and. swidth(ij).gt.0d0) then
        call sample_get_x_old( wgt, x, j, ipole, xmin, xmax )
        goto 999
      endif
      if (ituple.ne.1) then
        call sample_get_x_old( wgt, x, j, ipole, xmin, xmax )
        goto 999
      endif
      if (nzoom.gt.0) then
        call sample_get_x_old( wgt, x, j, ipole, xmin, xmax )
        goto 999
      endif
c
c Proceed with the simplified new implementation
c (NB: ituple=1 and nzoom<=0 are expected from now on)
c
      if(xmax.ne.1 .or. .not.xbin_max1_saved(j,ipole)) then
        xbin_max = xbin(xmax, minvar(j,ipole))
        if(xmax.eq.1) then
          xbin_max1_array(j,ipole) = xbin_max
          xbin_max1_saved(j,ipole) = .true.
        endif
      else
        xbin_max = xbin_max1_array(j,ipole)
      endif

      if(xmin.ne.0 .or. .not.xbin_min0_saved(j,ipole)) then
        xbin_min = xbin(xmin, minvar(j,ipole))
        if (xbin_min .gt. xbin_max-1) then
          xbin_min = min(xbin_min, xbin_max)
        endif
        if(xmin.eq.0) then
          xbin_min0_array(j,ipole) = xbin_min
          xbin_min0_saved(j,ipole) = .true.
        endif
      else
        xbin_min = xbin_min0_array(j,ipole)
      endif

      call ntuple_new(ddum)
      ddum = xbin_min + ddum * (xbin_max - xbin_min)
      tx(1,j) = xbin_min
      tx(2,j) = ddum
      tx(3,j) = xbin_max
      im = ddum
      if (im.ge.ng)then
        im = ng -1
        ddum = ng
      endif
      if (im.lt.0) im = 0
      ip = im + 1
c------
c     tjs 3/5/2011  save bin used to avoid looking up when storing wgt
c------
      lastbin(j) = ip
c
c     New method of choosing x from bins
c
      if (ip .eq. 1) then ! first bin
        xo = grid(2,ip,ij) - xgmin
        x = grid(2,ip,ij) - xo * (dble(ip) - ddum)
      else           
        xo = grid(2,ip,ij) - grid(2,im,ij)
        x = grid(2,ip,ij) - xo * (dble(ip) - ddum)
      endif
      wgt = wgt * xo * dble(xbin_max - xbin_min)
 999  CONTINUE
c     CALL COUNTERS_STOP_COUNTER( 10 ) ! 10=PROGRAM-SampleGetX
      end

c-------------------------------------------------------------------------------
