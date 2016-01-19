#include "Hathor.h"

using namespace std;

const double Hathor::z2  = 1.6449340668482264365;   /* Zeta(2) */
const double Hathor::z3  = 1.2020569031595942854;   /* Zeta(3) */
const double Hathor::ln2 = 0.69314718055994530942;  /* log(2)  */

const double Hathor::inv4pi  = 1./(4.*M_PI);
const double Hathor::inv4pi2 = inv4pi*inv4pi;
const double Hathor::inv4pi4 = inv4pi2*inv4pi2;

const double Hathor::ca  = 3.0;
const double Hathor::ca2 = ca*ca;
const double Hathor::cf  = 4.0/3.0;

Hathor::Hathor(Pdf & pdf_) 
  : AbstractHathor(pdf_),
    mw(80.385),mz(91.1876),mh(126.),mb(4.82),swq(1.0 - mw*mw/(mz*mz)),
    alpha(1./126.3), weak(mb,mz,mw,mh,alpha,swq,hcq)
{
  // weak.check();
  nf = 5;
  Cqq = 0;
  Cgg = 0;
  scheme = LO | NLO | NNLO | POLE_MASS;
  setab_(ABSCHEME);
  update();
  info("Dear user,");
  info("please note that per default the exact NNLO results [1-4] are used");
  info("for the evaluation of the inclusive cross section");
  info("The weak corrections are based on [5-7]"); 
  info("For the numerical integration Hathor uses an adopted version ");
  info("of the Vegas integrator [8,9]. The required random numbers are ");
  info("generated using Ranlux [10]. The one-loop integrals are calculated");
  info("using the FF library. [11]");
  info("[1]  Baernreuther, Czakon, Mitov: Percent level precision physics at ");
  info("     the Tevatron: first genuine NNLO QCD corrections to ");
  info("     q qbar -> t tbar + X,  Phys.Rev.Lett. 109 (2012) 132001.");
  info("[2]  Czakon, Mitov: NNLO corrections to top-pair production at ");
  info("     hadron colliders: the all-fermionic scattering channels,");
  info("     JHEP 1212 (2012) 054");
  info("[3]  Czakon, Mitov: NNLO corrections to top pair production at");
  info("     hadron colliders: the quark-gluon reaction,");
  info("     JHEP 1301 (2013) 080"); 
  info("[4]  Czakon, Fiedler, Mitov: The total top quark pair production ");
  info("     cross-section at hadron colliders through O(alpha_S^4),");
  info("     Phys.Rev.Lett. 110 (2013) 25, 252004"); 
  info("[5]  K¨uhn, Scharf and Uwer: Electroweak corrections to");
  info("     top-quark pair production in quark-antiquark annihilation");
  info("     Eur. Phys. J. C 45 (2006) 139 [hep-ph/0508092].");
  info("[6]  K¨uhn, Scharf and Uwer: Electroweak effects in");
  info("     top-quark pair production at hadron colliders,");  
  info("     Eur. Phys. J. C 51 (2007) 37 [hep-ph/0610335].");
  info("[7]  K¨uhn, Scharf and Uwer: Weak Interactions in"); 
  info("     Top-Quark Pair Production at Hadron Colliders: An Update,"); 
  info("     arXiv:1305.5773 [hep-ph]."); 
  info("[8]  Lepage:A New Algorithm for Adaptive Multidimensional Integration");
  info("     J.Comput.Phys.27 (1978) 192");
  info("[9] Lepage: VEGAS: AN ADAPTIVE MULTIDIMENSIONAL INTEGRATION PROGRAM");
  info("     CLNS-80/447 (1980).");
  info("[10] Luescher: A Portable high quality random number generator");
  info("     for lattice field theory simulations,");
  info("     Comput.Phys.Commun.79 (1994) 100 [hep-lat/930902]");
  info("[11] van Oldenborgh: FF: A Package to evaluate one loop");
  info("     Feynman diagrams, Comput.Phys.Commun 66 (1991) 1");
  info("");
  info("Please cite the corresponding articles if you quote Hathor results");
  info("in your article relying on the aformentioned calculations.");
}

void Hathor::update() {
  AbstractHathor::update();
  as = as_pdf[0];
  as2 = as*as;
  a = as/M_PI;
  
  /* Update QCD beta function in case that nf has changed */
  b0 = ( 11.0/3.0*ca - 4.0/3.0*nf/2.0 ) * inv4pi2;
  b1 = (+ 34.0/3.0*ca2 - 20.0/3.0*ca*nf/2.0 - 4.0*cf*nf/2.0 ) * inv4pi4;
  
  /* 
   * In the definition of d1,d2 the renormalisation scale mur 
   * is set to m, as a consequence no logs appear here.
   */
  d1 = 4./3.;
  d2 = 307./32. + 2.*z2 + 2./3.*z2*ln2 - 1./6.*z3 
    - ( 71./144. + 1./3.*z2 )*nf;
  /*
   * Update nf in the Fortran fit functions:
   */
  setnf_(nf);
}

#define PROPTS(_X_) if ( _X_ & scheme) std::cout << "Hathor::" << #_X_ <<" | ";

void Hathor::PrintOptions() {
  cout << " HATHOR: Options: ";
  PROPTS(LO);
  PROPTS(NLO);
  PROPTS(NNLO);
  PROPTS(NNLOAPPROX);
  PROPTS(NOHIGHENERGY);
  PROPTS(LOG_ONLY);
  PROPTS(PDF_SCAN);
  PROPTS(PDF_SYM_ERR);
  PROPTS(MS_MASS);
  PROPTS(NNPDF);
  PROPTS(NOQG);
  PROPTS(NOGG);
  PROPTS(NOQQB);
  PROPTS(NOQQ);
  PROPTS(NOQQP);
  PROPTS(NOQQPB); 
  std::cout << std::endl;
}

#undef PROPTS

void Hathor::setNf(int n) {
  nf = n;
  info("Number of massless flavours set to", nf);
  info("Note that this will not affect the scheme for alpha_s");
  info("provided by the pdf set.");
  setnf_(nf);
}

void Hathor::setScheme(unsigned int newscheme) {
  if (LOG_ONLY & newscheme) {
    if (MS_MASS & newscheme) {
      err("Incompatible use: In LOG_ONLY the MS_MASS scheme");
      err("is not supported. If your really believe you need it,");
      err("please contact the authors.");
      exit(1);
    }
    if (! (NNLO & newscheme)) {
      err("LOG_ONLY makes only sense in NNLO order");
      err("Please use in addition Hathor::NNLO as option");
    }
  }
  
  if (NNPDF & newscheme)
    newscheme |= PDF_SYM_ERR;
  
  if (PDF_SCAN & newscheme) {
    pdfmax = pdf.NumberPdf()+1;
    if (pdfmax <= 2) {
      err("Pdf does not seem to support error pdfs, at least ");
      err("the number of sets returned by the pdf is <= 1");
      exit(1);
    }
  } else
    pdfmax = 1;
  
  scheme = newscheme;
}

void Hathor::setCqq(double tmp) {
  Cqq = tmp;
  info("Cqq set to ", Cqq);
};

void Hathor::setCgg(double tmp) {
  Cgg = tmp;
  info("Cgg set to ", Cgg);
};

void Hathor::setMass(PARTONS parton, double mass){
  switch(parton){
  case ABOTTOM:
  case BOTTOM:
    mb = mass;
    info("b-quark mass is set to mb = ",mb);
    break;
  case WBOSON:
    mw = mass;
    info("W boson mass is set to mw = ",mw);
    break;
  case ZBOSON:
    mz = mass;
    info("Z boson mass is set to mz = ",mz);
    break;
  case HIGGS:
    mh = mass;
    info("Higgs boson mass is set to mh = ",mh);
    break;
  default:
    info("Only mw,mz,mw,mh entering the weak coorections for");
    info("top-quark pair production can be set using setmass()");
    info("The top mass is set dynamically when the cross sections");
    info("are calculated.");
  }
  weak.updateCouplings();
}

void Hathor::setHighEnergyScheme(ABSCHEMES s) {
  if (scheme & NOHIGHENERGY) {
    warn("The option Hathor::NOHIGHENERGY is used, the high energy");
    warn("approximation is thus not used. You should remove first the");
    warn("option Hathor::NOHIGHENERGY before changing the scheme");
  }
  if (s == ASCHEME) {
    info("A scheme for the high-energy approximation is used");
    info("as described in arXiv:1203.6282.");
    setab_(ASCHEME);
    return;
  }
  
  if (s == BSCHEME) {
    info("B scheme for the high-energy approximation is used");
    info("as described in arXiv:1203.6282.");
    setab_(BSCHEME);
    return;
  }
  
  if (s == ABSCHEME) {
    info("(A+B)/2 scheme for the high-energy approximation is used");
    info("as described in arXiv:1203.6282.");
    setab_(ABSCHEME);
    return;
  }
  
  err("Error in setHighEnergyScheme: unknown scheme");
  err("Use Hathor::ASCHEME, Hathor::BSCHEME or Hathor::AB2SCHEME");
  exit(1);
  return;
}

void Hathor::setPartonicEnergy(const double x[]) {
  /*
   *  Define kinematics
   */
  sparton = 4*m2*pow(shad/4./m2,x[0]);
  double r = sparton/shad;
  x1 = r *pow(r,-x[1]);
  x2 = r/x1;
  jacobian = -r*log(shad/4./m2)*log(r);
  
  beta = sqrt(1.- 4.*m2/sparton);
  rho =   4.*m2/sparton;
  beta2 = 1.- rho;
  beta4 = beta2*beta2;
  
  // delta for derivatives
  delta = x2*0.00001;
}

void Hathor::evaluatePDFs(const double h1[], const double h2[],
			  const double h2left[], const double h2right[]) {
  gg  = fluxgg(h1,h2);
  qqb = fluxqqb(h1,h2);
  qq = fluxqq(h1,h2);
  qqp = fluxqqp(h1,h2);
  qqpb = fluxqqpb(h1,h2);
  qg  = fluxqg(h1,h2);
  gqb = fluxgqb(h1,h2);
  
  dgg  = x2/(2.*delta) * ( fluxgg(h1,h2right) -  fluxgg(h1,h2left) );
  dqqb = x2/(2.*delta) * (fluxqqb(h1,h2right) - fluxqqb(h1,h2left) );
  dqg  = x2/(2.*delta) * ( fluxqg(h1,h2right) -  fluxqg(h1,h2left) );
  dgqb = x2/(2.*delta) * (fluxgqb(h1,h2right) - fluxgqb(h1,h2left) );

  if ( NOQG & scheme )
    qg = gqb = dqg = dgqb = 0.;
  if ( NOQQB & scheme )
    qqb = dqqb = 0.;
  if ( NOQQ & scheme )
    qq = 0;
  if ( NOGG & scheme )
    gg = dgg  = 0.;
  if ( NOQQP & scheme )
    qqp = 0;
  if ( NOQQPB & scheme )
    qqpb = 0.;
}

void Hathor::evaluateNLO(double beta) {
  fgg10 = tfgg10_(beta); fgg11 = tfgg11_(beta);  
  fqqb10 = tfqqb10_(beta); fqqb11 = tfqqb11_(beta);
  fgq10 = tfgq10_(beta); fgq11 = tfgq11_(beta);
  
  nlogg = 4.0 * M_PI * as * (fgg10 + Lmf * fgg11 + 2.0*b0 * Lmr * fgg00);
  nloqqb = 4.0 * M_PI * as * (fqqb10 + Lmf * fqqb11 + 2.0*b0 * Lmr * fqqb00);
  nlogq = 4.0 * M_PI * as * (fgq10 + Lmf*fgq11);
}

void Hathor::evaluateNNLO(double beta) {
  double xfqqb00, xfgg00;
  double xfqqb10, xfgg10, xfqqb11, xfgg11;
  double not_log_enhanced = 1.;

  double as2tmp = pow(4.0 * M_PI * as,2); 
  
  if ( LOG_ONLY & scheme ) {
    
    double lnbeta  = log(beta);
    double lnbeta2 = lnbeta*lnbeta;
    double lnbeta3 = lnbeta2*lnbeta;
    double lnbeta4 = lnbeta3*lnbeta;
    
    not_log_enhanced = 0.;
    //xfqqb00 = fqqb00; xfgg00 = fgg00;
    xfqqb00 = M_PI/9.*beta; xfgg00 = 7./192.*M_PI*beta;
    
    fqqb20 = xfqqb00 * inv4pi4 * 
      ( + 3.60774/beta/beta
	+ 1/beta*(-140.368*lnbeta2+32.106*lnbeta+3.95105)
	+ 910.222*lnbeta4-1315.53*lnbeta3+592.292*lnbeta2
	+ 528.557*lnbeta);
    
    fgg20 = xfgg00 * inv4pi4 * 
      ( + 68.5471/beta/beta
	+ 1/beta*(496.3*lnbeta2+321.137*lnbeta-8.62261)
	+ 4608*lnbeta4 - 1894.91*lnbeta3 
	- 912.349*lnbeta2+2456.74*lnbeta);
    
    fqqb21 = xfqqb00 * inv4pi4 *
      ( + 1/beta*(-49.3480 + 70.1839*lnbeta)
	- 910.222*lnbeta3+1358.99*lnbeta2
	- 476.012*lnbeta+163.524);
    
    fgg21 = xfgg00 * inv4pi4 * 
      ( + 1/beta*(39.6351-248.150*lnbeta)
	-4608*lnbeta3+2606.77*lnbeta2+368.035*lnbeta
	-171.000);
    
    fqqb22 = xfqqb00 * inv4pi4 * (227.555*lnbeta2-377.874*lnbeta+96.3456);
    
    fgg22 = xfgg00 * inv4pi4 * (1152*lnbeta2-890.989*lnbeta-104.291);
    
    xfqqb10 = xfqqb00 * inv4pi2 * (42.66667*lnbeta2
				 -20.6105*lnbeta+14.5312-3.28986/beta);
    xfgg10 = xfgg00 * inv4pi2 * (96*lnbeta2
				 -9.51647*lnbeta+25.8548+5.16979/beta);
    
    xfqqb11 = xfqqb00 * inv4pi2 * (-21.3333*lnbeta+13.8795);
    xfgg11 = xfgg00 * inv4pi2 * (-48*lnbeta+14.7289);
    
    nnlogq = 0;

    fqqpb20 = fqqp20 = fqq20 = 0.;

  } else { // LOG_ONLY else
    /* We use always the exact scheme dependence: */
    fgg21 = tfgg21_(beta); fgg22 = tfgg22_(beta);
    fqqb21 = tfqqb21_(beta); fqqb22 = tfqqb22_(beta);
    fgq21 = tfgq21_(beta); fgq22 = tfgq22_(beta);
    fqq21 = tfqq21_(beta); fqq22 = tfqq22_(beta);
    fqqp22 = fqqpb22 = fqq22;
    fqqp21 = fqqpb21 = tfqqp21_(beta);

    if ( NNLOAPPROX & scheme ){
      if ( NOHIGHENERGY & scheme ) {
	fqqb20 = tfqqb20approx_(beta);
	fgq20 = tfgq20approx_(beta);
	fgg20 = tfgg20approx_(beta); 
      } else {
	fqqb20 = tfqqb20_asym_(beta);
	fgq20 = tfgq20_asym_(beta); 
	fgg20 = tfgg20_asym_(beta); 
      }
    } else {
      /* Default: use exact results */
      fqqb20 = tfqqb20_exact_(beta);
      fgq20 = tfgq20_(beta); 
      fgg20 = tfgg20_(beta); 
    }
  
    /* Remaining channels: */
    fqq20 = tfqq20_(beta);
    fqqp20 = tfqqp20_(beta);
    fqqpb20 = tfqqpb20_(beta);
  
    xfqqb10 = fqqb10; xfgg10 = fgg10; xfqqb11 = fqqb11; xfgg11 = fgg11;
    xfqqb00 = fqqb00; xfgg00 = fgg00;
    
    nnlogq = as2tmp * 
      (
       + fgq20 + Lmf*fgq21 + Lmf*Lmf*fgq22
       + 3.0*b0*Lmr* fgq10 
       + 3.0*b0*Lmr*Lmf* fgq11
       );  
  } // LOG_ONLY else   
  
  nnlogg = as2tmp * 
    ( + fgg20 + Lmf*fgg21 + Lmf*Lmf*fgg22 
      + Lmr * ( 3.0*b0*xfgg10 + 2.0 * b1 * fgg00 * not_log_enhanced )
      + 3.0 * b0 * Lmr*Lmf* xfgg11  
      + 3.0 * b0*b0 * Lmr*Lmr* fgg00 * not_log_enhanced 
      + xfgg00*Cgg*inv4pi4 );  

  nnloqqb = as2tmp * 
    ( + fqqb20 + Lmf*fqqb21 + Lmf*Lmf*fqqb22 
      + Lmr * ( 3.0 * b0 * xfqqb10 + 2.0 * b1 * fqqb00 * not_log_enhanced )
      + 3.0 * b0 * Lmr*Lmf * xfqqb11 
      + 3.0 * b0*b0 * Lmr*Lmr * fqqb00 * not_log_enhanced 
      + xfqqb00*Cqq*inv4pi4 );  
  
  // additional qq, qqp, qqpb channels
  nnloqq = as2tmp * ( + fqq20 + Lmf*fqq21 + Lmf*Lmf*fqq22 );
  nnloqqp = as2tmp * ( + fqqp20 + Lmf*fqqp21 + Lmf*Lmf*fqqp22 );
  nnloqqpb = as2tmp * ( + fqqpb20 + Lmf*fqqpb21 + Lmf*Lmf*fqqpb22 );

}

void Hathor::evaluateScalingFunctions() {
  // leading-order scaling functions
  fgg00 = tfgg00_(beta);
  fqqb00 = tfqqb00_(beta);
  
  // next-to-leading-order scaling functions and cross sections
  if ((scheme & NNLO) || (scheme & NLO))
    evaluateNLO(beta);
  
  // next-to-next-to-leading order scaling functions and cross sections
  if (scheme & NNLO)
    evaluateNNLO(beta);
}

double Hathor::evaluateIntegral(double asi, double wgt) {
  double lo = 0, nlo = 0, nnlo = 0;
  
  // build integrand
  if ( scheme & LO ) {
    lo =  gg * fgg00 + qqb * fqqb00;
  }
  
  if ( scheme & NLO ) {
    nlo =  gg  * nlogg + qqb * nloqqb + ( qg + gqb ) * nlogq;
    if ( scheme & MS_MASS ){
      nlo += + 2.* d1 * a * ( dgg * fgg00 + dqqb * fqqb00 );
    }
  }
  
  if ( scheme & NNLO ) {
    nnlo = gg * nnlogg + qqb * nnloqqb + ( qg + gqb ) * nnlogq;
    nnlo += qq * nnloqq + qqp * nnloqqp + qqpb * nnloqqpb ;
    
    if ( scheme & MS_MASS ) {
      double mdiffgg = +1./192.*M_PI*rho/beta 
	* ( beta*(36.-40*beta2+4*beta4)*log((1.+beta)/(1.-beta))
	    -7.-116.*beta2+91.*beta4);
      double mdiffqqb = -1./9.*M_PI*pow(rho,3)/beta;
      
      nnlo += 
	+ dgg * a * ( + (2.*d2-d1*d1) * a * fgg00 + 2.*d1 * nlogg
		      + d1*d1 * a *mdiffgg ) 
	+ dqqb * a * ( + (2.*d2-d1*d1) * a * fqqb00 + 2.*d1 * nloqqb
		       + d1*d1 * a * mdiffqqb ) 
	+ ( dqg + dgqb ) * 2.*d1 * a * nlogq
	+ (-2.*d1*a*4.*M_PI*as) * ( gg*fgg11 + qqb*fqqb11 + (qg+gqb)*fgq11 )
	+ 2.*d1* a * 4.*M_PI*as* log(mur2/m2)*b0 
	* ( dgg * fgg00 + dqqb * fqqb00) ;
    }
  }
  
  //result *= as2/m2  * jacobian * hcq;
  return as2 * pow(asi/as,2) / m2  * jacobian * hcq 
    * ( lo + (asi/as) * nlo + pow(asi/as,2) * nnlo );
  
}

double Hathor::fluxgg(const double h1[], const double h2[]) {
  return(h1[GLUON] * h2[GLUON]);
}

double Hathor::fluxqqb(const double h1[], const double h2[]) {
  return( // Q QB
	 + h1[UP] * h2[AUP] 
	 + h1[DOWN] * h2[ADOWN]
	 + h1[STRANGE] * h2[ASTRANGE] 
	 + h1[CHARM] * h2[ACHARM] 
	 + h1[BOTTOM] * h2[ABOTTOM] 
	 // QB Q
	 + h1[AUP] * h2[UP] 
	 + h1[ADOWN] * h2[DOWN]
	 + h1[ASTRANGE] * h2[STRANGE] 
	 + h1[ACHARM] * h2[CHARM] 
	 + h1[ABOTTOM] * h2[BOTTOM] 
	  );
}

double Hathor::fluxqq(const double h1[], const double h2[]) {
  return(
	 + h1[UP] * h2[UP] 
	 + h1[DOWN] * h2[DOWN]
	 + h1[STRANGE] * h2[STRANGE] 
	 + h1[CHARM] * h2[CHARM] 
	 + h1[BOTTOM] * h2[BOTTOM] 

	 + h1[AUP] * h2[AUP] 
	 + h1[ADOWN] * h2[ADOWN]
	 + h1[ASTRANGE] * h2[ASTRANGE] 
	 + h1[ACHARM] * h2[ACHARM] 
	 + h1[ABOTTOM] * h2[ABOTTOM] 
	 );
}
double Hathor::fluxqqp(const double h1[], const double h2[]) {
  return( 
	 + h1[UP] * ( h2[DOWN] + h2[STRANGE] + h2[CHARM] + h2[BOTTOM])
	 + h1[DOWN] * ( h2[UP] + h2[STRANGE] + h2[CHARM] + h2[BOTTOM])
	 + h1[STRANGE] * ( h2[UP] + h2[DOWN] + h2[CHARM] + h2[BOTTOM])
	 + h1[CHARM] * ( h2[UP] + h2[DOWN] + h2[STRANGE] + h2[BOTTOM])
	 + h1[BOTTOM] * ( h2[UP] + h2[DOWN] + h2[STRANGE] + h2[CHARM])

	 + h1[AUP] * ( h2[ADOWN] + h2[ASTRANGE] + h2[ACHARM] + h2[ABOTTOM])
	 + h1[ADOWN] * ( h2[AUP] + h2[ASTRANGE] + h2[ACHARM] + h2[ABOTTOM])
	 + h1[ASTRANGE] * ( h2[AUP] + h2[ADOWN] + h2[ACHARM] + h2[ABOTTOM])
	 + h1[ACHARM] * ( h2[AUP] + h2[ADOWN] + h2[ASTRANGE] + h2[ABOTTOM])
	 + h1[ABOTTOM] * ( h2[AUP] + h2[ADOWN] + h2[ASTRANGE] + h2[ACHARM])

);
}
double Hathor::fluxqqpb(const double h1[], const double h2[]) {
  return( 
	 + h1[UP] * ( h2[ADOWN] + h2[ASTRANGE] + h2[ACHARM] + h2[ABOTTOM] ) 
	 + h1[DOWN] * ( h2[AUP] + h2[ASTRANGE] + h2[ACHARM] + h2[ABOTTOM] ) 
	 + h1[STRANGE] * ( h2[AUP] + h2[ADOWN] + h2[ACHARM] + h2[ABOTTOM] ) 
	 + h1[CHARM] * ( h2[AUP] + h2[ADOWN] + h2[ASTRANGE] + h2[ABOTTOM] ) 
	 + h1[BOTTOM] * ( h2[AUP] + h2[ADOWN] + h2[ASTRANGE] + h2[ACHARM] ) 
	 
	 + h1[AUP] * ( h2[DOWN] + h2[STRANGE] + h2[CHARM] + h2[BOTTOM] ) 
	 + h1[ADOWN] * ( h2[UP] + h2[STRANGE] + h2[CHARM] + h2[BOTTOM] ) 
	 + h1[ASTRANGE] * ( h2[UP] + h2[DOWN] + h2[CHARM] + h2[BOTTOM] ) 
	 + h1[ACHARM] * ( h2[UP] + h2[DOWN] + h2[STRANGE] + h2[BOTTOM] ) 
	 + h1[ABOTTOM] * ( h2[UP] + h2[DOWN] + h2[STRANGE] + h2[CHARM] ) 
	  );
}

double Hathor::fluxqg(const double h1[], const double h2[]) {
  return(+ h1[GLUON] 
	 * ( h2[UP] + h2[DOWN] + h2[STRANGE] + h2[CHARM] + h2[BOTTOM] )
	 + ( h1[UP] + h1[DOWN] + h1[STRANGE] + h1[CHARM] + h1[BOTTOM] ) 
	 * h2[GLUON]
	 );
}

double Hathor::fluxgqb(const double h1[], const double h2[]) {
  return(+ h1[GLUON] 
	 * ( h2[AUP] + h2[ADOWN] + h2[ASTRANGE] + h2[ACHARM] + h2[ABOTTOM] )
	 + ( h1[AUP] + h1[ADOWN] + h1[ASTRANGE] + h1[ACHARM] + h1[ABOTTOM] ) 
	 * h2[GLUON] );
}
