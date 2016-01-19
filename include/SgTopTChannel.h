#ifndef SGTOPTCHANNEL_H_
#define SGTOPTCHANNEL_H_

#include "SgTop.h"
#include "PartonicCrossSection.h"
#include "SplittingFunction.h"

class HathorSgTopT_LO_qq : public SgTopCrossSection {
 public:
  virtual double eval(double s);
};

class HathorSgTopT_LO_qqbar : public SgTopCrossSection {
 public:
  virtual double eval(double s);
};


class HathorSgTopT : public SgTop {
	
 public:
  HathorSgTopT(Pdf &pdf_);

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

  HathorSgTopT_LO_qq    pxs_qq_00;
  HathorSgTopT_LO_qqbar pxs_qqb_00;
  
  SampledPartonicCrossSection pxs_qq_10,pxs_qqb_10,
    pxs_gu_10,pxs_gb_10,pxs_gqb_10;
  
  SplittingFunction xPqq0, xPqg0, xPgg0, xPgq0, xPqg1, xPqq1s, 
    xPqq1v, xPqqb1v, xPqggg, xPqggq, xPqgqg, xPqqqg, xPqqqq;
  
  double qq, qqbar,gqbar, gb, gu,gg, uq_ig, dq_ig, dbarq_ig, qqbv, qqv;
  double fqq00, fqqbar00;
  
  double nloqqbar, nloqq, nlogqbar, nlogb, nlogu;
  double nnloqqbar, nnloqq, nnlogqbar, nnlogb, nnlogu, nnlogg, 
    nnlouqig, nnlodqig, nnlodbarqig, nnloqqbv, nnloqqv;
  
  // approx (N)NLO
  double  rt, rs4;
  double  approxqq1, approxqqbar1,approxqq2, approxqqbar2;

  void approxNNLO(double t, double s4, double result[8]);
  void    evaluateApproxNNLO();
  double  ampsq_qq_00(double t);
  double  ampsq_qqbar_00(double t);
  
  // analytically derived expressions
  double xs_qq_11();
  double xs_gu_11();
  
  double xs_qqbar_11();
  double xs_gqbar_11();
  double xs_gb_11();
  
  double xs_qq_21();
  double xs_qq_22();
  double xs_qqbar_21();
  double xs_qqbar_22();
  double xs_gu_21();
  double xs_gu_22();
  double xs_gb_21();
  double xs_gb_22();
  double xs_gqbar_21();
  double xs_gqbar_22();
  double xs_gg_21();
  double xs_gg_22();
  double xs_uqig_21();
  double xs_uqig_22();
  double xs_dqig_21();
  double xs_dqig_22();
  double xs_dbarqig_21();
  double xs_dbarqig_22();
  double xs_qqbv_21();
  double xs_qqbv_22();
  double xs_qqv_21();
  double xs_qqv_22();
};

#endif // SGTOPTCHANNEL_H_
