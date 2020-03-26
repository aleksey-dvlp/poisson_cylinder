#if !defined(_PCYLIND_H)
#define      _PCYLIND_H

#include <fstream>


template <size_t N>
class PCylinder {
  
  double *pu; //numerical solution
  double *pua; //analitical solution
  const double R = 1.0;
  double h, inverseh2;
  const double bound = 2.0;

  std::ofstream file;

  public:
  PCylinder();
  ~PCylinder();
  void SolvePoisson ();
  void Analitical ();
  void write_in_file();
  inline double rhside(double r) {return r*r;}
  double inf_norm (); //compute infinity norm

};




#endif


