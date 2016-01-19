#ifndef SGTOPSCHANNEL_H_
#define SGTOPSCHANNEL_H_

#include "SgTop.h"
#include "PartonicCrossSection.h"
#include "SplittingFunction.h"

class HathorSgTopS_LO_qqbar : public SgTopCrossSection {
 public:
  virtual double eval(double s);
};

class HathorSgTopS : public SgTop {
  
 public:
  HathorSgTopS(Pdf &pdf_);
  
 protected:
  virtual void   setPartonicEnergy(const double x[]);
  virtual void   evaluatePDFs(const double h1[], const double h2[],
			      const double h2left[] = 0, 
			      const double h2right[] = 0);
  virtual void   evaluateScalingFunctions();
  virtual double evaluateIntegral(double as, double wgt);
  virtual void   update();

private:
  double sparton, sqrts, rho;
  double jacobian;
  
  HathorSgTopS_LO_qqbar pxs_qqbar_00;
  SampledPartonicCrossSection pxs_qqbar_10,pxs_gq_10,pxs_gqbar_10;
  
  SplittingFunction xPqq0, xPgq0, xPgg0, xPqg0, xPqg1, xPqq1v, 
    xPqq1s, xPqqb1v,xPqqqg, xPggqg, xPqqqq, xPgqqg, xPqgqg, 
    xPqggg, xPqggq;
  
  double qqbar;
  double fqqbar00;

  double gqbar, gq;
  double nloqqbar, nlogqbar, nlogq;
  
  double qqbar2, gg, qq, qqv;
  double nnloqqbar, nnloqqbar2, nnlogq, nnlogg, nnloqq, nnloqqv;
  
  // approx NNLO
  double  rt, rs4;
  double  approxqqbar1, approxqqbar2;
  void approxNNLO(double t, double s4, double result[8]);
  void    evaluateApproxNNLO();
  double  xs_qqbar_00(double t);
  
  // analytically derived expressions
  double xs_qqbar_11();
  double xs_gq_11();
  
  double xs_qqbar_21();
  double xs_qqbar_22();
  double xs_qq_21();
  double xs_qq_22();
  double xs_gg_21();
  double xs_gg_22();
  double xs_gq_21();
  double xs_gq_22();
  double xs_qqv_21();
  double xs_qqv_22();
};

#endif // SGTOPSCHANNEL_H_
