      program test

      integer array_i(3)
      real array_f(3)
      double precision array_d(3)
      complex*8 array_cf(3)
      complex*16 array_cd(3)

      integer matrix_i(4,5), r/4/, c/5/
      integer i,j

 22   format(I4)

      data array_i  / 1,2,3 /
      data array_f  / 1,2,3 /
      data array_d  / 1,2,3 /
      data array_cf / (1,1),(2,2),(3,3) /
      data array_cd / (1,1),(2,2),(3,3) /

      write(*,*) ( array_i(i), i=1,3)
      call foo_i(array_i)
      write(*,*) ( array_i(i), i=1,3)

      write(*,*)

      write(*,*) ( array_f(i), i=1,3)
      call foo_f(array_f)
      write(*,*) ( array_f(i), i=1,3)

      write(*,*)

      write(*,*) ( array_d(i), i=1,3)
      call foo_d(array_d)
      write(*,*) ( array_d(i), i=1,3)

      write(*,*)

      write(*,*) ( array_cf(i), i=1,3)
      call foo_cf(array_cf)
      write(*,*) ( array_cf(i), i=1,3)

      write(*,*)

      write(*,*) ( array_cd(i), i=1,3)
      call foo_cd(array_cd)
      write(*,*) ( array_cd(i), i=1,3)

      do 20 i = 1,r
         do 10 j = 1,c
            matrix_i(i,j) = i*10 + j
 10      continue
 20   continue

      do 40 i = 1,r
        do 30 j = 1,c
          write(6,22,advance="no") matrix_i(i,j)
 30     continue
        write(*,*)
 40   continue
      call foo_mi(matrix_i, r, c)

      write(*,*)

      stop
      end
