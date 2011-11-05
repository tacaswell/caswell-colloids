%module test_wrap
%{
#include "test.h"
%}
int t_fun(int j);

class dummy{
 public:
  int x;
  int y;
  int d_fun();
};



