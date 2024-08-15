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
      integer icount, it_warned
      integer im, ip, ij
      double precision xbin_min, xbin_max, ddum(maxdim), xo, y
c
c     External
c
      double precision xbin
      external         xbin
c
c     Global
c
      double precision tmean, trmean, tsigma
      integer             dim, events, itm, kn, cur_it, invar, configs
      common /sample_common/
     .  tmean, trmean, tsigma, dim, events, itm, kn, cur_it, invar, configs

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

      data ddum/maxdim*0d0/
      data icount/0/
      data it_warned/0/

      integer            lastbin(maxdim)
      common /to_lastbin/lastbin

c-----
c  Begin Code
c-----
      CALL COUNTERS_START_COUNTER( 10, 1 ) ! 10=PROGRAM-SampleGetX
      ij = Minvar(j,ipole)
c
c Fall back to the default old implementation under some circumstances
c
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
        call ntuple(ddum(j), xbin_min,xbin_max, j, ipole)
      else
        call ntuple(ddum(j),max(xbin_min,dble(int(tx(2,j)))),
     $    min(xbin_max,dble(int(tx(2,j))+1)),j,ipole)
      endif
      tx(1,j) = xbin_min
      tx(2,j) = ddum(j)
      tx(3,j) = xbin_max

      im = ddum(j)
      if (im.ge.ng)then
        im = ng -1
        ddum(j) = ng
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
        x = grid(2, ip, ij) - xo * (dble(ip) - ddum(j))
      else           
        xo = grid(2, ip, ij)-grid(2,im,ij)
        x = grid(2, ip, ij) - xo * (dble(ip) - ddum(j))
      endif
c
c     Simple checks to see if we got the right point note 1e-3 corresponds
c     to the fact that the grids are required to be separated by 1e-14. Since
c     double precision is about 18 digits, we expect things to agree to
c     3 digit accuracy.
c
      if (it_warned .ne. cur_it) then
        icount=0
        it_warned = cur_it
      endif
      if (abs(ddum(j)-xbin(x,ij))/(ddum(j)+1d-22) .gt. 1e-3) then
        if (icount .lt. 5) then
          write(*,'(a,i4,2e14.6,1e12.4)')
     &      'Warning xbin not returning correct x', ij,
     &      ddum(j),xbin(x,ij),xo
        elseif (icount .eq. 5) then
          write(*,'(a,a)')'Warning xbin still not working well. ',
     &      'Last message this iteration.'
        endif
        icount=icount+1
      endif
      wgt = wgt * xo * dble(xbin_max-xbin_min)
 999  CALL COUNTERS_STOP_COUNTER( 10 ) ! 10=PROGRAM-SampleGetX
      end

c-------------------------------------------------------------------------------
