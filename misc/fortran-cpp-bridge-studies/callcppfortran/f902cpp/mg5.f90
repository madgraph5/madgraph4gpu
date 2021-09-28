program test_main

use test_module

type(test_type) :: t

call i_new(t)

call i_hello(t)

call i_delete(t)

end program test_main
