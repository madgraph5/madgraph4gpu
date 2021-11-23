      program test_matrix

      integer m(5,6)

      do 20 x = 1,5
        do 10 y = 1,6
          m(x,y) = x*10 + y
  10    continue
  20  continue

C      do 40 x = 1,5
C        do 30 y = 1,6
C          write(*,*) m(x,y)
C  30    continue
C  40  continue

      do 50 x = 1,30
        write(*,*) m(x)
  50  continue

      stop
      end
