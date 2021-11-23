#include <iostream>

using namespace std;

struct Point {
  double x;
  double y;
};

struct Rectangle {
  Point corner;
  double width;
  double height;
};

class pointrectangle {

    Rectangle r;

  public:

    pointrectangle (Rectangle);

    Point findCenter ()
    {
      double x = r.corner.x + r.width/2;
      double y = r.corner.y + r.height/2;
      Point result = {x, y};
      return result;
    }

    double findArea ()
    {
       return (r.width*r.height);
    }

    void printInfo()
    {
       Point center = findCenter();
       cout << "Center Point coordinates: " << endl;
       cout << "    x-value: " << center.x << endl;
       cout << "    y-value: " << center.y << endl;
       cout << "The area is: " << findArea() << endl;
    }
};

//------------------------------
// Constructor of pointrectangle
//------------------------------

pointrectangle::pointrectangle (Rectangle box) {
  r.corner.x = box.corner.x;
  r.corner.y = box.corner.y;
  r.width    = box.width;
  r.height   = box.height;
}

//-------------------------------------
// C wrapper interfaces to C++ routines
//-------------------------------------

extern "C" {
  pointrectangle *pointrectangle__new (Rectangle a) {
    return new pointrectangle(a);
  }

  Point pointrectangle__findCenter (pointrectangle *This) {
    return This->findCenter();
  }

  double  pointrectangle__findArea (pointrectangle *This) {
    return This->findArea();
  }

  void  pointrectangle__printInfo (pointrectangle *This) {
    return This->printInfo();
  }

  void pointrectangle__delete (pointrectangle *This) {
    delete This;
  }
}
