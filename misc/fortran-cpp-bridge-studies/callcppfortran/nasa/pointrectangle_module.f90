module pointrectangle_module

use, intrinsic :: ISO_C_Binding, only: C_int, C_double
use, intrinsic :: ISO_C_Binding, only: C_ptr, C_NULL_ptr

implicit none

private

TYPE, BIND(C) :: fPoint
REAL(C_double) :: x
REAL(C_double) :: y
END TYPE fPoint

TYPE, BIND(C) :: fRectangle
TYPE(fPoint) :: corner
REAL(C_double) :: width
REAL(C_double) :: height
END TYPE fRectangle

type pointrectangle_type
private
type(C_ptr) :: object = C_NULL_ptr
end type pointrectangle_type

!------------------------
! C function declarations
!------------------------

interface

function C_pointrectangle__new (a) result(this) bind(C,name="pointrectangle__new")
  import
  type(C_ptr) :: this
  TYPE(fRectangle), value :: a
end function C_pointrectangle__new

subroutine C_pointrectangle__delete (this) bind(C,name="pointrectangle__delete")
  import
  type(C_ptr), value :: this
end subroutine C_pointrectangle__delete

function C_pointrectangle__findCenter (this) result(center) bind(C,name="pointrectangle__findCenter")
  import
  TYPE(fPoint) :: center
  type(C_ptr), value :: this
end function C_pointrectangle__findCenter

function C_pointrectangle__findArea (this) result(area) bind(C,name="pointrectangle__findArea")
import
  REAL(C_double) :: area
  type(C_ptr), value :: this
end function C_pointrectangle__findArea

subroutine C_pointrectangle__printInfo (this) bind(C,name="pointrectangle__printInfo")
  import
  type(C_ptr), value :: this
end subroutine C_pointrectangle__printInfo

end interface

interface new
  module procedure pointrectangle__new
end interface new

interface delete
  module procedure pointrectangle__delete
end interface delete

interface findArea
  module procedure pointrectangle__findArea
end interface findArea

interface printInfo
  module procedure pointrectangle__printInfo
end interface printInfo

interface findCenter
  module procedure pointrectangle__findCenter
end interface findCenter

public :: new, delete, findCenter, findArea, printInfo
public :: pointrectangle_type, fPoint, fRectangle

!------------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------

!-------------------------------------------------
! Fortran wrapper routines to interface C wrappers
!-------------------------------------------------
subroutine pointrectangle__new(this,a)
  type(pointrectangle_type), intent(out) :: this
  TYPE(fRectangle) :: a
  this%object = C_pointrectangle__new(a)
end subroutine pointrectangle__new

subroutine pointrectangle__delete(this)
  type(pointrectangle_type), intent(inout) :: this
  call C_pointrectangle__delete(this%object)
  this%object = C_NULL_ptr
end subroutine pointrectangle__delete

subroutine pointrectangle__printInfo(this)
  type(pointrectangle_type), intent(in) :: this
  call C_pointrectangle__printInfo(this%object)
end subroutine pointrectangle__printInfo

function pointrectangle__findCenter(this) result(center)
  type(pointrectangle_type), intent(in) :: this
  TYPE(fPoint) :: center
  center = C_pointrectangle__findCenter(this%object)
end function pointrectangle__findCenter

function pointrectangle__findArea(this) result(area)
  type(pointrectangle_type), intent(in) :: this
  REAL(C_double) :: area
  area = C_pointrectangle__findArea(this%object)
end function pointrectangle__findArea
!------------------------------------------------------------------------------

end module pointrectangle_module
