#include <iostream>

class Test {

public:

  void hello() {
    std::cout << "hello world" << std::endl;
  }

};

// C wrapper

extern "C" {

  Test *test_new_c() {
    return new Test();
  }

  void test_hello_c(Test *This) {
      This->hello();
  }

  void test_delete_c(Test *This) {
    delete This;
  }

}
