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
      double precision xbin_min, xbin_max, ddum, xo, y
c
c     External
c
      double precision xbin
      external         xbin
c
c     Global
c

      double precision    grid(2, ng, 0:maxinvar)
      common /data_grid/ grid

      integer           Minvar(maxdim,lmaxconfigs)
      common /to_invar/ Minvar

      integer           ituple
      common /to_random/ituple

      double precision      spole(maxinvar),swidth(maxinvar),bwjac
      common/to_brietwigner/spole        ,swidth        ,bwjac

      integer nzoom
      double precision  tx(1:3,maxinvar)
      common/to_xpoints/tx, nzoom

      data ddum/0d0/

      integer            lastbin(maxdim)
      common /to_lastbin/lastbin

      integer icall
      save icall
      data icall/0/

c-----
c  Begin Code
c-----
      CALL COUNTERS_START_COUNTER( 10, 1 ) ! 10=PROGRAM-SampleGetX
      icall = icall + 1
c     if (icall.le.100) then
c       write(6,*) 'sample_get_x', wgt, x, j, ipole, xmin, xmax
c     endif
      ij = Minvar(j,ipole)
c
c Fall back to the default old implementation under some circumstances
c - first call  => take care of ranmar initialization in the old ntuple
c - swidth > 0  => take care of calling transpole
c - ituple != 1 => take care of error message      
c
      if (icall.eq.1) then ! first call
        call sample_get_x_old( wgt, x, j, ipole, xmin, xmax )
        return
      endif
      if (ij.gt.0 .and. swidth(ij).gt.0d0) then
        call sample_get_x_old( wgt, x, j, ipole, xmin, xmax )
        goto 999
      endif
      if (ituple.ne.1) then
        call sample_get_x_old( wgt, x, j, ipole, xmin, xmax )
        goto 999
      endif
c
c Proceed with the simplified new implementation
c (NB: ituple=1 is expected from now on)
c
      xbin_min = xbin(xmin,minvar(j,ipole))
      xbin_max = xbin(xmax,minvar(j,ipole))
      if (xbin_min .gt. xbin_max-1) then
        xbin_max = xbin(xmax,minvar(j,ipole))
        xbin_min = min(xbin(xmin,minvar(j,ipole)), xbin_max)
      endif
c
c     Line which allows us to keep choosing same x
c
      if (nzoom .le. 0) then
        call ntuple_new(ddum, xbin_min, xbin_max)
      else
        call ntuple_new(ddum, max(xbin_min, dble(int(tx(2,j)))),
     $    min(xbin_max,dble(int(tx(2,j))+1)))
      endif
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
      if (ip .eq. 1) then       !This is in the first bin
        xo = grid(2, ip, ij)-xgmin
        x = grid(2, ip, ij) - xo * (dble(ip) - ddum)
      else           
        xo = grid(2, ip, ij)-grid(2,im,ij)
        x = grid(2, ip, ij) - xo * (dble(ip) - ddum)
      endif
c
      wgt = wgt * xo * dble(xbin_max-xbin_min)
 999  CALL COUNTERS_STOP_COUNTER( 10 ) ! 10=PROGRAM-SampleGetX
      end

c-------------------------------------------------------------------------------
