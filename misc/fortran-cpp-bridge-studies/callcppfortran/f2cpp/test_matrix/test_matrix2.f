      program test

      integer m(4,5,6), mo/4/, pa/5/, ev/6/
      integer i,j
      integer V(15)
      POINTER ( P, V )

 22   format(I4)

      do 5 k = 1,ev
      do 20 i = 1,mo
         do 10 j = 1,pa
              m(i,j,k) = k*100 + i*10 + j
 10      continue
 20   continue
 5         continue

      do 25 k = 1,ev
        do 40 i = 1,mo
          do 30 j = 1,pa
             write(6,22,advance="no") m(i,j,k)
 30      continue
         write(*,*)
 40   continue
      write(*,*)
 25        continue

      write(*,*)

      P = LOC(m)

      do 50 i = 1,mo*pa*ev
         write(6,22,advance="no") V(i)
 50   continue

      write(*,*)

      stop
      end
