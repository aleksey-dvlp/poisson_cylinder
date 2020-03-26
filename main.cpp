//
// 1   d       du
//--- --- ( r ---- ) = r^2
// r   dr      dr
// u(R) = 2
//     du |
// - r ---|    = 0
//     dr |
//        |r=0
// R = 1
//
// analitical solution u(r) = 31/16 + r^4/16
//
//

#include "pcylind.h"
#include "pcylind.cpp"

int main (int argc, char* argv[]) {

  const size_t N = 10;
  PCylinder<N> cyl;
  cyl.SolvePoisson();
  cyl.Analitical();
  cyl.write_in_file();
  cyl.inf_norm();

}
