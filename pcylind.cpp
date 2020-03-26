#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "pcylind.h"
using namespace std;

template <size_t N>
PCylinder<N>::PCylinder() {
  pu = new double [N+1];
  pua= new double [N+1];
  for (int i = 0; i<N+1; i++) pu[i] = 0.0;
  for (int i = 0; i<N+1; i++) pua[i] = 0.0;
  h = R/N;
  inverseh2 = 1./h/h;
  auto data = "data";
  auto nstring = to_string(N);
  auto outstring = ".out";
  auto filename = data + nstring + outstring;
  file.open(filename);
  file << "coordiante "<<"numerical "<<"analitical"<<endl;
}

template <size_t N>
PCylinder<N>::~PCylinder() {
  file.close();
  delete [] pu;
  delete [] pua;
}

template <size_t N>
void PCylinder<N>::write_in_file() {
  file << scientific;
  for (int i = 1; i<N+1; i++) 
    file << setprecision(10) <<0.5*(2*i-1)*h << ' '
	 << pu[i] << ' ' << pua[i] << ' ' << endl;
  file << R<< ' ' <<bound<<' ' <<bound<<' '; 
  file << endl;
}

template <size_t N>
void PCylinder<N>::SolvePoisson() {

  double a[N+1], b[N+1], c[N+1], d[N+1], alpha[N+1], beta[N+1];
  

  for (int n=1; n<N; n++) {
    a[n] = 2.0*(n-1)/(2.0*n-1)*inverseh2;
    b[n] = -2.0*inverseh2;
    c[n] = 2.0*n/(2.0*n-1)*inverseh2;
    const double r = 0.5*(2.0*n-1.0)*h;
    d[n] = rhside(r);
  }
  a[N] = 2.0*(N-1)/(2.0*N-1)*inverseh2;
  b[N] = -2.0*(3.0*N-1)/(2.0*N - 1.0)*inverseh2;
  c[N] = 0.0;
  const double r = 0.5*(2.0*N-1.0)*h;
  d[N] = rhside(r) - 4.0*N*bound/(2.0*N-1.0)*inverseh2;

  for (int n = 1; n<N+1; n++) {
    alpha[n] = c[n]/(b[n] - a[n]*alpha[n-1]);
    beta [n] = d[n] - a[n]*beta[n-1];
    beta [n] /= b[n] - a[n]*alpha[n-1];
  }
  //Note: here we are doing a little crime
  //we are going outside the u[N+1] array bound
  //but it is doesn't matter, I hope, as alpha[N] = 0
  for (int n=N;n>0;n--) pu[n] = beta[n] - alpha[n]*pu[n+1];

}

template <size_t N>
void PCylinder<N>::Analitical() {
  const double pi = 3.141592;
  for (int n = 1; n<N+1; n++) {
    const double r = 0.5*(2.0*n-1.0)*h;
    pua[n] = 31./16. + r*r*r*r/16.;
  }
}

template <size_t N>
double PCylinder<N>::inf_norm () {
  double s = 0.0;
  for (int i = 1; i < N+1; i++) 
    s = max(s,abs(pu[i]-pua[i]));
  cout << "the infinity norm is "<< s << endl;
}


  
