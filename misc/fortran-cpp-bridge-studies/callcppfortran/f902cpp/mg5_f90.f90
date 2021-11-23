module test_module

use, intrinsic :: ISO_C_Binding, only: C_ptr, C_NULL_ptr

! implicit none

private

type test_type
  private
  type(C_ptr) :: object = C_NULL_ptr
end type test_type

! trampoline functions

interface

function f2c_test_new() result(this) bind(C,name="test_new_c")
  import
  type(C_ptr) :: this
end function f2c_test_new

subroutine f2c_test_delete(this) bind(C,name="test_delete_c")
  import
  type(C_ptr), value :: this
end subroutine f2c_test_delete

subroutine f2c_test_hello(this) bind(C,name="test_hello_c")
  import
  type(C_ptr), value :: this
end subroutine f2c_test_hello

end interface

! Declare interfaces

interface i_new
  module procedure w_test_new
end interface i_new

interface i_delete
  module procedure w_test_delete
end interface i_delete

interface i_hello
  module procedure w_test_hello
end interface i_hello

public :: i_new, i_delete, i_hello
public :: test_type

! Fortran wrapper routines

CONTAINS

  subroutine w_test_new(this)
    type(test_type), intent(out) :: this
    this%object = f2c_test_new()
  end subroutine w_test_new

  subroutine w_test_delete(this)
    type(test_type), intent(inout) :: this
    call f2c_test_delete(this%object)
    this%object = C_NULL_ptr
  end subroutine w_test_delete

  subroutine w_test_hello(this)
    type(test_type), intent(inout) :: this
    call f2c_test_hello(this%object)
  end subroutine w_test_hello


end module test_module
