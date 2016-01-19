#ifndef SGTOPWTCHANNEL_H_
#define SGTOPWTCHANNEL_H_

#include "SgTop.h"
#include "PartonicCrossSection.h"
#include "SplittingFunction.h"



class HathorSgTopWt_LO_gq : public SgTopCrossSection {
 public:
  virtual double eval(double s);
};



class HathorSgTopWt : public SgTop {
  
 public:
  HathorSgTopWt(Pdf &pdf_);

 protected:
  virtual void   setPartonicEnergy(const double x[]);
  virtual void   evaluatePDFs(const double h1[], const double h2[],
			      const double h2left[] = 0, 
			      const double h2right[] = 0);
  virtual void   evaluateScalingFunctions();
  virtual double evaluateIntegral(double as, double wgt);
  virtual void   update();
  
private:
  static const bool enableBB = false;
  
  double sparton, sqrts, rho;
  double jacobian;
  
  HathorSgTopWt_LO_gq pxs_gq_00;
  SampledPartonicCrossSection pxs_gq_10,pxs_gg_10,pxs_qq_10;
  
  SplittingFunction xPqq0,  xPqg0,  xPgg0,  xPgq0, xPgg1, xPgq1, xPqg1,  
    xPqq1s, xPqq1v, xPqqb1v,xPgggg, xPgqgg, xPqggg, xPqggq,  xPqqgg,  
    xPqqgq,xPqqqg, xPqqqq;
	
  double gb, gb_h;
  double fgq00;

  double gg, bq, bqbar;
  double nlogq, nlogg, nlobq;
  
  double gqig, gqbarig, gqbarig2, gqbarig2_h;
  double nnlogq, nnlogg, nnlobq, nnlogqig, nnlogqbig2, nnlogq_h,nnlogqbig2_h;
  
  // approx NNLO
  double  rt, rs4;
  double  approxgq1, approxgq2;
  void approxNNLO(double t, double s4, double result[]);
  void    evaluateApproxNNLO();
  double  xs_gq_00(double t, double u);
  
  // scale dependencies
  double xs_gq_11();
  double xs_gg_11();
  double xs_qq_11();
  
  double xs_gq_21();
  double xs_gq_21_h();
  double xs_gq_22();
  double xs_gq_22_h();
  double xs_gg_21();
  double xs_gg_22();
  double xs_qq_21();
  double xs_qq_22();
  double xs_gqig_21();
  double xs_gqig_22();
  double xs_gqbig2_21();
  double xs_gqbig2_21_h();
  double xs_gqbig2_22();
  double xs_gqbig2_22_h();
};

#endif // SGTOPWTCHANNEL_H_
