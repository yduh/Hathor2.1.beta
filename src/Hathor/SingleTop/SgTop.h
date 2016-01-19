#ifndef SGTOP_H_
#define SGTOP_H_

#include "AbstractHathor.h"



// The following constants MUST NOT BE CHANGED for reasons of consistency!
const double pi          = M_PI;
const double pi2         = pi*pi;
const double pi4         = pi2*pi2;
const double twopi       = 2.*pi;
const double invpi       = 1. / pi;
const double invpi2      = 1. / pi2;
const double invpi4      = 1. / pi4;
const double fourpi      = 4. * pi;
const double fourpi2     = fourpi*fourpi;
const double mw          = 80.385;           // W mass (GeV)
const double mw2         = mw*mw;
const double mw4         = mw2*mw2;

const double MCFMalpha   = 1. / 132.2332298; // alpha(M_Z) as set in MCFM
const double MCFMswq    = 0.2228972;        // sin^2(theta_W)

// color
const double Nc = 3.;
const double CF = (Nc*Nc - 1.) / (2. * Nc);
const double CA = Nc;

// number of light quark flavors
const double nf = 5.;

// QCD beta function
const double beta0 = (11.*CA - 2.*nf)/3.;
const double beta1 = (34.*CA*CA - 10.*CA*nf)/3. - 2.*CF*nf;
// QCD beta function for use in scale-dependant terms
const double b0    = beta0 / (16.  * pi2);
const double b1    = beta1 / (256. * pi4);

const double zeta2 = pi2/6.;
const double zeta3 = 1.20205690315959428539973816151144999;
                     
const double K  = CA * (67./18. - pi2/6.) - 5.*nf/9.;
const double D2 = + CF * CA * ( -101./54. + 11./6.*zeta2 + 1.75*zeta3 ) 
  + CF * nf * ( 7./27. - zeta2/3. );
const double D2g = (CA/CF)*D2;


class SgTop : public AbstractHathor {
	
public:
  enum PARTICLE { TOPQUARK=0, ANTITOPQUARK=1, BOTH=2 };
  enum APPROXSCHEMES { NNLOAPPROX=(1<<10), NLOAPPROX=(1<<16)};
  enum LEADINGLOG { NLO_LL=(1<<0), NLO_NLL=(1<<1), NLO_VIRT=(1<<2),
		    NNLO_LL=(1<<3), NNLO_NLL=(1<<4), NNLO_NNLL=(1<<5),
		    NNLO_NNNLL=(1<<6), NLL_CALCULATION=(1<<7) };
  enum APPROXTEST { FAKE_NNLO=(1<<11) };
  
  SgTop(Pdf &pdf_);
  
  void setAlpha(const double value);
  void setSwq(const double value);
  void setParticle(PARTICLE particle_);
  void setCkmMatrix(const double ckm_[3][3]);
  void getCkmMatrix(double ckm_[3][3]);
  void PrintCkmMatrix();
  
  void setLLflags(const unsigned int flags);
  
  void PrintContributions();
  
 protected:
  virtual void update();
  
  double m4;
  long llflags;
  
  PARTICLE particle;
  double ckm[3][3]; 
  double &Vud, &Vus, &Vub, &Vcd, &Vcs, &Vcb, &Vtd, &Vts, &Vtb;
  double alpha, alpha2, swq;
  
  double contribution[6];
};

#endif // SGTOP_H_
