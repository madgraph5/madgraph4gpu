program pointrectangle_main

use pointrectangle_module

type(pointrectangle_type) :: rect
type(fPoint) :: center
type(fRectangle) :: mybox
double precision :: area

mybox%corner%x  = 10.0d0
mybox%corner%y  =  0.0d0
mybox%width = 100.0d0
mybox%height = 200.0d0

call new(rect, mybox)

center = findCenter(rect)
area   = findArea(rect)

write(*,*) '-------------------------------------'
write(*,*) 'Print out from C++'
write(*,*) '-------------------------------------'

call printInfo(rect)

write(*,*) '-------------------------------------'
write(*,*) 'Print out from Fortran'
write(*,*) '-------------------------------------'
write(*,*) 'Center Point x: ', center%x
write(*,*) 'Center Point y: ', center%y
write(*,*) 'Area:', area

call delete(rect)

end program pointrectangle_main
