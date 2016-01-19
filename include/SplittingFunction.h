#ifndef SPLITTING_FUNCTIONS_H_
#define SPLITTING_FUNCTIONS_H_

#include <cmath>
#include "PartonicCrossSection.h"
#include "SgTop.h"


class SplittingFunction {

public:
  enum KERNEL { PQQ0  = 1, PQG0  = 2, PGG0   = 3, PGQ0 = 4,
		PQQ1S = 5, PQQ1V = 6, PQQB1V = 7,
		PQG1  = 8, PGG1  = 9, PGQ1   = 10 };

  SplittingFunction(KERNEL kernel);
  SplittingFunction(KERNEL kernel, KERNEL kernel2);

  double convolve(double rho, double rnd, PartonicCrossSection & xs);
  double convolve(double rho, double r1, double r2, PartonicCrossSection & xs);

private:
  static const int MAX_CONVOLUTIONS = 2;
  
  double (*continuous_part[MAX_CONVOLUTIONS])(double);
  double (*alpha_part[MAX_CONVOLUTIONS])(double);
  double (*gamma_part[MAX_CONVOLUTIONS])(double);

  void setKernel(KERNEL kernel, double (*&c_part)(double),
		 double (*&a_part)(double), double (*&g_part)(double));

  static double Pqq0(double x);
  static double Pqg0(double x);
  static double Pgq0(double x);
  static double Pgg0(double x);
  
  static double Pqq1v(double x);
  static double Pqqb1v(double x);
  static double Pqq1s(double x);
  
  static double Pqg1(double x);
  static double Pgq1(double x);
  static double Pgg1(double x);
  
  static double zero(double x);	
  static double alpha_qq0(double x);
  static double alpha_qq1(double x);
  static double alpha_gg0(double x);
  static double alpha_gg1(double x);
  static double gamma_qq0(double x);
  static double gamma_qq1(double x);
  static double gamma_gg0(double x);
  static double gamma_gg1(double x);
  
};

#endif // SPLITTING_FUNCTIONS_H_
