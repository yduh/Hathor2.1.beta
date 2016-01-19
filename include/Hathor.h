#ifndef HATHOR_H_
#define HATHOR_H_

#include "AbstractHathor.h"
#include "HathorFits.h"
#include "HathorPdf.h"
#include "SgTop.h"
#include "SgTopTChannel.h"
#include "SgTopSChannel.h"
#include "SgTopWtChannel.h"
#include "HathorWeakCorrections.h"


class Hathor : public AbstractHathor {
  
public:
  enum MASSRENORMALISATIONSCHEME { POLE_MASS=0, MS_MASS=1 };
  enum APPROXSCHEMES { LOG_ONLY=(1<<4), NOQG= (1<<8), NOHIGHENERGY= (1<<9), 
		       NNLOAPPROX=(1<<10), NOGG = (1<<11), NOQQB = (1<<12),
		       NOQQ = (1<<13), NOQQP = (1<<14), NOQQPB = (1<<15)};
  enum ABSCHEMES { ASCHEME=0, BSCHEME=1, ABSCHEME=2 };
  
  static const double z2,z3,ln2;
  static const double ca,ca2,cf;
  
  Hathor(Pdf & pdf_);
  void setMass(PARTONS parton, double mass);
  void setSwq(double newval){weak.setSwq(newval);};
  void setAlpha(double newval){weak.setAlpha(newval);};
  void setNf(int n);
  void setCqq(double tmp);
  void setCgg(double tmp);
  void setHighEnergyScheme(ABSCHEMES scheme);
  inline double dsigmagg(const double mt, const double shat, const double z){
    return(weak.dsigmagg(mt,shat,z));};
  inline double dsigmaqq(const double mt, const double shat, const double z){
    return(weak.dsigmaqq(mt, shat, z));
  };
  inline double dsigmaWeakgg(const double mt, const double shat, 
			     const double z){
    return(weak.dsigmaWeakgg(mt,shat,z));};
  inline void dsigmaWeakqq(const double mt, const double shat, const double z,
			   double & up, double & down){
    weak.dsigmaWeakqq(mt, shat, z, up, down);
  };

  virtual void setScheme(unsigned int newscheme);
  virtual void PrintOptions();

protected:
  virtual void update();
  
  virtual void setPartonicEnergy(const double x[]);
  virtual void evaluatePDFs(const double h1[], const double h2[],
			    const double h2left[] = 0, 
			    const double h2right[] = 0);
  virtual void evaluateScalingFunctions();
  virtual double evaluateIntegral(double as, double wgt);

private:
  void evaluateNLO(double beta);
  void evaluateNNLO(double beta);
  static inline double fluxgg(const double h1[], const double h2[]);
  static inline double fluxqqb(const double h1[], const double h2[]);
  static inline double fluxqq(const double h1[], const double h2[]);
  static inline double fluxqqp(const double h1[], const double h2[]);
  static inline double fluxqqpb(const double h1[], const double h2[]);
  static inline double fluxqg(const double h1[], const double h2[]);
  static inline double fluxgqb(const double h1[], const double h2[]);
  
  static const double inv4pi, inv4pi2, inv4pi4;

  double mw,mz,mh,mb,swq,alpha;
  WeakCorrections weak;
  double sparton, jacobian;
  double beta, rho, beta2, beta4;
  

  double gg, qqb, qq, qqp, qqpb, qg, gqb;
  double dgg, dqqb, dqg, dgqb;
  
  int nf;
  double as,as2,a;
  double b0, b1,d1,d2;
  double Cqq, Cgg;
  double fgg00,fqqb00;
  double fgg10, fgg11, fqqb10, fqqb11, fgq10, fgq11, nlogg,nloqqb,nlogq;
  //
  double fgg20, fgg21, fgg22;
  double fqqb20, fqqb21, fqqb22,fqq20,fqqp20,fqqpb20, fqq21,fqq22,
    fqqp21,fqqp22,fqqpb21,fqqpb22;
  double fgq20, fgq21,fgq22;
  double nnlogg,nnloqqb,nnloqq, nnloqqp, nnloqqpb, nnlogq;
  
};
#endif // HATHOR_H_
