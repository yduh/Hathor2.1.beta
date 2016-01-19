#include "SplittingFunction.h"
#include "SgTop.h"
#include "AbstractHathor.h"

/*
 * 
 * Conventions from Ellis, Stirling, Webber: QCD and Collider Physics
 * Cambridge University Press 1996 [Ellis,Stirling,Webber '96]
 *
 */
const double TR=0.5;
const double TF=nf/2.0;
/*
 * To match different conventions define translation factors:
 * P0 = ESW_P0 * P0_ESW, P1 = ESW_P1 * P1_ESW
 *      ^^^^^^                ^^^^^^
 */
const double ESW_P0=1.0/(8.0*pi2);
const double ESW_P1=1.0/(64.0*pi4);

SplittingFunction::SplittingFunction(KERNEL kernel) {
  /*
   *  Constructor for convolution: P_ij \otimes f_kl
   */
  setKernel(kernel, continuous_part[0], alpha_part[0], gamma_part[0]);
}

SplittingFunction::SplittingFunction(KERNEL kernel1, KERNEL kernel2) {
  /*
   *  Constructor for convolution: P_ij \otimes P_kl \otimes f_mn
   */  
  setKernel(kernel1, continuous_part[0], alpha_part[0], gamma_part[0]);
  setKernel(kernel2, continuous_part[1], alpha_part[1], gamma_part[1]);
}

void SplittingFunction::setKernel(KERNEL kernel, double (*&c_part)(double),
				  double (*&a_part)(double), 
				  double (*&g_part)(double)) {
  switch (kernel) {
  case PQQ0:   c_part = &Pqq0;
    a_part = &alpha_qq0;
    g_part = &gamma_qq0;
    break;
  case PQG0:   c_part = &Pqg0;
    a_part = &zero;
    g_part = &zero;
    break;
  case PGG0:   c_part = &Pgg0;
    a_part = &alpha_gg0;
    g_part = &gamma_gg0;
    break;
  case PGQ0:   c_part = &Pgq0;
    a_part = &zero;
    g_part = &zero;
    break;
  case PQQ1V:  c_part = &Pqq1v;
    a_part = &alpha_qq1;
    g_part = &gamma_qq1;
    break;
  case PQQ1S:  c_part = &Pqq1s;
    a_part = &zero;
    g_part = &zero;
    break;
  case PQQB1V: c_part = &Pqqb1v;
    a_part = &zero;
    g_part = &zero;
    break;
  case PQG1:   c_part = &Pqg1;
    a_part = &zero;
    g_part = &zero;
    break;
  case PGG1:   c_part = &Pgg1;
    a_part = &alpha_gg1;
    g_part = &gamma_gg1;
    break;
  case PGQ1:   c_part = &Pgq1;
    a_part = &zero;
    g_part = &zero;
    break;
  }
}

double SplittingFunction::convolve(double rho, double rnd, 
				   PartonicCrossSection & xs) {

  double result = 0.;
  double newarg = rnd * (1. - rho) + rho;
  double xsrho  = xs.evalr(rho);
  double xsnew  = xs.evalr(rho/newarg);
  
  result = xsrho * (gamma_part[0](1.) * log(1. - rho) + alpha_part[0](1.));
  
  result += (1. - rho) 
    * (continuous_part[0](newarg)*xsnew
       + (gamma_part[0](newarg)*xsnew - gamma_part[0](1.)*xsrho)/(1.-newarg));
  
  return result;
}

double SplittingFunction::convolve(double rho, double r1, double r2,
				   PartonicCrossSection & xs) {

  double result = 0.;
  double newarg = r1 * (1. - rho) + rho;
  double xsrho  = convolve(rho,        r2, xs);
  double xsnew  = convolve(rho/newarg, r2, xs);
  
  result = xsrho * (gamma_part[1](1.) * log(1. - rho) + alpha_part[1](1.));
  
  result += (1. - rho) 
    * (continuous_part[1](newarg)*xsnew
       + (gamma_part[1](newarg)*xsnew - gamma_part[1](1.)*xsrho)/(1.-newarg));
	return result;
}

/*
 * Splitting functions in LO
 * Eq. (4.94) in [Ellis,Stirling,Webber '96], continuous part
 */
double SplittingFunction::Pqq0(double x) {
  return CF * ESW_P0 * (-1. - x);
}

double SplittingFunction::Pqg0(double x) {
  return TR * ESW_P0 * (1. - 2.*x + 2.*x*x);
}

double SplittingFunction::Pgq0(double x) {
  return CF * ESW_P0 * (2./x - 2. + x);
}

double SplittingFunction::Pgg0(double x) {
  return 2.0 * CA * ESW_P0 * (1./x - 2. + x - x*x);
}

// For splitting functions without a delta or plus distribution
double SplittingFunction::zero(double x) {
  return 0.;
}

// From Eq. (4.94) in [Ellis,Stirling,Webber '96]: alpha*delta() + gamma * Plus_Distribution
double SplittingFunction::alpha_qq0(double x) {
  return CF * ESW_P0 * 3.0/2.0;
}

double SplittingFunction::gamma_qq0(double x) {
  return CF * ESW_P0 * 2.0;
}

double SplittingFunction::alpha_gg0(double x) {
  return ESW_P0*(11.0*CA-4.0*nf*TR)/6.0;
}

double SplittingFunction::gamma_gg0(double x) {
  return 2.0 * CA * ESW_P0 * 1.0;
}

/*
 * Splitting functions in NLO
 *
 * qq case, eq. (4.107) in [Ellis,Stirling,Webber '96]:
 */
double SplittingFunction::Pqq1v(double x) {
  const double log_x=log(x);
  return(ESW_P1 * ( CF * CF 
		    * (-(2.0*log_x*log(1.0-x)+1.5*log_x)*(-1-x)
		       -(1.5+3.5*x)*log_x-0.5*(1.0+x)*log_x*log_x-5.0*(1.0-x))
		    + CF*CA *( (0.5*log_x*log_x+11.0/6.0*log_x+67.0/18.0
			       -pi2/6.0)*(-1.0-x)
			       +(1.0+x)*log_x+20.0/3.0*(1.0-x))
		    + CF * TF * (-(2.0/3.0*log_x+10.0/9.0)*(-1.0-x)
				 -4.0/3.0*(1.0-x))
		    ));
}

/*
 * Eq. (4.108) in [Ellis,Stirling,Webber '96]:
 */
double SplittingFunction::Pqqb1v(double x) {
  const double log_x=log(x);
  const double log_plus=log(1.0+x);
  return(
    CF*ESW_P1*(CF-CA/2.0)*(2.0*(2.0/(1.0+x)-1.0+x)*(-2.0*AbstractHathor::Li2(-x)
    +0.5*log_x*log_x-2.0*log_x*log_plus-pi2/6.0)+2.0*(1.0+x)*log_x+4.0*(1.0-x))
    );
}

/*
 * Eq. (4.109) in [Ellis,Stirling,Webber '96], 
 * Pqq1s=(Pqq1-PqqV1-PqqbarV1)/(2nf)
 * Pqq1s(x) -> 0 near threshold --> loss of digits 
 */
double SplittingFunction::Pqq1s(double x) {
  const double log_x=log(x);
  return(CF/4.0*ESW_P1
	 *(-4.0+12.0*x-112.0/9.0*x*x+40.0/(9.0*x)
	   +(10.0*x+16.0/3.0*x*x+2.0)*log_x-2.0*(1.0+x)*log_x*log_x)
	 );
}

/*
 * Eq. (4.119) in [Ellis,Stirling,Webber '96]
 */
double SplittingFunction::alpha_qq1(double x) {
  return (ESW_P1*(
		  + CF*CF*(3./8. - pi2/2.0 + 6.*zeta3)
		  + CA*CF*(17./24. + 11.0*pi2/18.0 - 3.*zeta3)
		  - CF*TF*(1.0/6.0 + 2.0*pi2/9.0)
		  ));
}

/*
 * Eq. (4.109) in [Ellis,Stirling,Webber '96], pqq() contains 2 * Plus_Distribution
 */
double SplittingFunction::gamma_qq1(double x) {
  if (x == 1.)
    return CF*(CA*(67./9. - pi2/3.)*2. + nf*(-20./9.))/(128.*pi4);
  const double log_x=log(x);
  return(ESW_P1*(
		 + CF * CF * (-(2.0*log_x*log(1.0-x)+1.5*log_x)*2.0)
		 + CF * CA * ((0.5*log_x*log_x+11.0/6.0*log_x+67.0/18.0
			       -pi2/6.0)*2.0)
		 + CF * TF * (-(2.0/3.0*log(x)+10.0/9.0)*2.0)
		 ));
}

/*
 * Eq. (4.110) in [Ellis,Stirling,Webber '96]
 * Note erratum for ESW, factor (2nf) included in eq. (4.110)
 */
double SplittingFunction::Pqg1(double x) {
  const double log_x = log(x);
  const double log_plus = log (1.0 + x);  
  const double log_minus = log (1.0 - x);  
  return(ESW_P1/4.0 *
	 (+ CF * ( + 4.0 -9.0*x-(1.0-4.0*x)*log_x-(1.0-2.0*x)*log_x*log_x
		 + 4.0*log_minus
		 + (2.0*(log_minus-log_x)*(log_minus-log_x)
		    -4.0*(log_minus-log_x)
		    -2.0/3.0*pi2+10.0)
		 *(1.0-2.0*x+2.0*x*x) )
	  + CA * (+ 182.0/9.0+14.0/9.0*x+40.0/(9.0*x)
		  + (136.0/3.0*x-38.0/3.0)*log_x
		  - 4.0*log_minus-(2.0+8.0*x)*log_x*log_x
		  + 2.0*(1.0+2.0*x+2.0*x*x)
		  * (-2.0*AbstractHathor::Li2(-x)
		     +0.5*log_x*log_x-2.0*log_x*log_plus-pi2/6.0)
		  + (-log_x*log_x+44.0/3.0*log_x
		     -2.0*log_minus*log_minus+4.0*log_minus
		     +pi2/3-218.0/9.0)*(1.0-2.0*x+2.0*x*x))
	  ));
}

/*
 * Eq. (4.111) in [Ellis,Stirling,Webber '96]
 */
double SplittingFunction::Pgq1(double x) {
  const double log_x=log(x);
  const double log_plus = log (1.0 + x);  
  const double log_minus = log (1.0 - x);  
  return( ESW_P1*(+ CF*CF*(-2.5-3.5*x+(2.0+3.5*x)*log_x-(1.0-0.5*x)*log_x*log_x
			   -2.0*x*log_minus
			   -(3.0*log_minus+log_minus*log_minus)*(2.0/x-2.0+x))
		  + CF*TF*(-4.0/3.0*x-(20.0/9.0+4.0/3.0*log_minus)
			   *(2.0/x-2.0+x) )
		  + CF*CA*(28.0/9.0+65.0/18.0*x+44.0/9.0*x*x
			   -(12.0+5.0*x+8.0/3.0*x*x)*log_x
			   +(4.0+x)*log_x*log_x+ 2.0*x*log_minus
			   +(-2.0*AbstractHathor::Li2(-x)+0.5*log_x
			     *log_x-2.0*log_x*log_plus-pi2/6)*(-2.0/x-2.0-x)
			   +(0.5-2.0*log_x
			     *log_minus+0.5*log_x*log_x+11.0/3.0*log_minus
			     +log_minus*log_minus
			     -pi2/6.0)*(2.0/x-2.0+x)))
	  ); 
}
/*
 * Eq.(4.112) in [Ellis,Stirling,Webber '96], uses Li2 with negative argument!
 */
double SplittingFunction::Pgg1(double x) {
  const double log_x=log(x);
  const double log_plus = log (1.0 + x);  
  const double log_minus = log (1.0 - x);  
  return(ESW_P1*(
		 CF*TF*(-16.0+8.0*x+20.0/3.0*x*x+4.0/(3.0*x)
			-(6.0+10.0*x)*log_x
			-(2.0+2.0*x)*log_x*log_x)
		 +CA*TF*(2.0-2.0*x+26.0/9.0*(x*x-1.0/x)-4.0/3.0*(1.0+x)*log_x
			 -20.0/9.0*(1.0/x-2.0+x*(1.0-x)))
		 +CA*CA*(27.0/2.0*(1.0-x)+67.0/9.0*(x*x-1.0/x)
			 -(25.0/3.0-11.0/3.0*x
			   +44.0/3.0*x*x)*log_x+4.0*(1.0+x)*log_x*log_x
			 +2.0*(1.0/(1.0+x) -1.0/x-2.0-x*(1.0+x))
			 *(-2.0*AbstractHathor::Li2(-x)
			   +0.5*log_x*log_x-2.0*log_x
			   *log_plus-pi2/6.0)
			 +(67.0/9.0-4.0*log_x*log_minus+log_x*log_x
			   -pi2/3)*(1.0/x-2.0+x*(1.0-x)))
		 ));
}

/*
 * Eq.(4.120) in [Ellis,Stirling,Webber '96]
 */
double SplittingFunction::alpha_gg1(double x) {
  return ESW_P1*(-4./3.*CA*TF + CA*CA*(8./3. + 3.*zeta3) - CF*TF);
}

/*
 * Eq.(4.112) in [Ellis,Stirling,Webber '96]
 */
double SplittingFunction::gamma_gg1(double x) {
  if (x == 1.)
    return (CA*TF*(-20./9.) + CA*CA*(67./9. - pi2/3.)) *ESW_P1;
  const double log_x=log(x);
  return(ESW_P1*(+ CA * TF * (-20.0/9.0)
		 + CA * CA *(67.0/9.0-4.0*log_x*log(1.0-x)+log_x*log_x-pi2/3.0))
	 );
}
