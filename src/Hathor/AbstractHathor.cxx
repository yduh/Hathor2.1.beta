#include "AbstractHathor.h"
#include "FF.h"
/*
 *  Output helper functions.
 */

double AbstractHathor::hcq = 0.389379323e9;
bool AbstractHathor::debug = false;

void AbstractHathor::info(std::string s) {
  if (debug) {
    std::cout << " HATHOR: " << s << std::endl;
  }
}

void AbstractHathor::info(std::string s, double x) {
  if (debug) {
    std::cout << " HATHOR: " << s << " " << x << std::endl;
  }
}

void AbstractHathor::warn(std::string s) {
  std::cerr << " HATHOR: WARNING: " << s << std::endl;
}

void AbstractHathor::err(std::string s) {
  std::cerr << " HATHOR: ERROR: " << s << std::endl;
}

/*
 *  Configuration functions.
 */

void AbstractHathor::setScheme(unsigned int newscheme) {
  if (NNPDF & newscheme)
    newscheme |= PDF_SYM_ERR;
  
  if (PDF_SCAN & newscheme) {
    pdfmax = pdf.NumberPdf() + 1;
    if (pdfmax <= 2) {
      err("PDF does not seem to support error pdfs, at least ");
      err("the number of sets returned by the pdf is <= 1");
      exit(1);
    }
  }
  else
    pdfmax = 1;
  scheme = newscheme;
}

void AbstractHathor::setPdfNumber(int npdf = 0) {
  if (npdf == 0) {
    // default depends on setting of PDF_SCAN
    if (PDF_SCAN & scheme) {
      pdfmax = pdf.NumberPdf() + 1;
    } else {
      pdfmax = 1;
    }
  } else {
    // limit evaluation of member PDFs
    if (npdf > (pdf.NumberPdf() + 1)) {
      err("Chosen number of pdfs is larger than then number of sets");
      err("returned by the pdf.");
      exit(1);
    }
    pdfmax = npdf;
    std::cout << " HATHOR: Evaluating " << pdfmax
	      << " / " << (pdf.NumberPdf() + 1)
	      << " member pdfs" << std::endl;
  }
}

void AbstractHathor::sethc2(double hcq_) {
  hcq = hcq_;
  info("(hc)^2 set to", hcq);
}

void AbstractHathor::setSqrtShad(double ecms) {
  info("Collider energy set to [GeV]: ", ecms);
  shad = ecms*ecms;
}

void AbstractHathor::setColliderType(COLLIDERTYPE type) {
  if (type == PP) {
    //shad = 8000.*8000.;
    shad = 13000.*13000.;
    collider = PP;
    info("Collider set to Proton-Proton");
    info("sqrt(shad) is set to [GeV]: ", sqrt(shad));
    info("Change if needed with setSqrtShad(newval)");
    return;
  }
  
  if (type == PPBAR) {
    shad = 1960.*1960.;
    collider = PPBAR;
    info("Collider set to Proton-AntiProton");
    info("sqrt(shad) is set to [GeV]: ", sqrt(shad));
    info("Change if needed with setSqrtShad(newval)");
    return;
  }
  
  // By construction we should never arrive here;
  err("Unknown collider used in setColliderType");
  exit(1);
}

void AbstractHathor::setPrecision(int n) {
  calls = n;
  info("Precision set to ", calls);
};

void AbstractHathor::setPDF(Pdf &newpdf) {
  pdf = newpdf;
}

/*
double AbstractHathor::getAlphas(double mur, int member) {
  pdf.InitMember(member);
  return pdf.GetAlphas(mur);
}
*/

#define PROPTS(_X_) if ( _X_ & scheme) std::cout << "Hathor::" << #_X_ <<" | ";
void AbstractHathor::PrintOptions() {
  std::cout << " HATHOR: Options: ";
  PROPTS(LO);
  PROPTS(NLO);
  PROPTS(NNLO);
  PROPTS(PDF_SCAN);
  PROPTS(PDF_SYM_ERR);
  std::cout << std::endl;
}
#undef PROPTS

AbstractHathor::AbstractHathor(Pdf &pdf_) : collider(PP), pdf(pdf_), 
					    calls(MEDIUM), pdfmax(1), 
					    shad(8000.*8000.) {
  debug = true;
  scheme = LO | NLO | NNLO;
  evalerr = 0;
  delta = 0.;
  dimension = 2;
  mur=muf=170.;
  update();
}

/*
 *  Access to results
 */

void AbstractHathor::getResult(int pdfset, double &integral_, 
                               double &err_, double & chi2a_) {
  
  if ( ! (scheme & PDF_SCAN) && (pdfset > 0) ){
    err("Use the option PDF_SCAN to evaluate pdf uncertainties");
    err("pdfset > 0 requires PDF_SCAN set first to work, otherwise");
    err("only the `central pdf' (pdfset=0) is evaluated.");
    exit(1);
  }
  
  integral_ = integral[pdfset]; 
  err_ = error[pdfset];
  chi2a_ = chi2a[pdfset];
  
  if (fabs(err_/integral_) > 0.02) {
    warn("Numerical integration error larger than 2 %");
    warn("Consider reruning with setPrecision set to a higher value");
  }
  
  if (fabs(err_/integral_) > 0.1) {
    err("Numerical integration returned error larger 10 %");
    err("Result is unreliable rerun with setPrecision set to higher value");
  }
  
}

void AbstractHathor::update() {
  /*
   *  Update constants in case that the setting has changed,
   *  i.e. alphas value, nf, or the scales, etc...
   */
  pdf.InitMember(0);
  as_pdf[0] = pdf.GetAlphas(mur);
  
  if (PDF_SCAN & scheme) {
    for (int ipdf = 1; ipdf < pdfmax; ipdf++) {
      pdf.InitMember(ipdf);
      as_pdf[ipdf] = pdf.GetAlphas(mur);
    }
    pdf.InitMember(0);
  }
  
  m2 = m*m;
  muf2 = muf*muf;
  mur2 = mur*mur;
  Lmf = log(muf2/m2);
  Lmr = log(mur2/muf2);
  
  /*
   * Clear result from last integration
   */
  for(int i = 0; i < 100; i++) {
    integral[i] = 0; error[i] = 0; chi2a[i] = 0;
  }
}


/*
 *  Execute Vegas.
 */
double AbstractHathor::getXsection(double m_, double mur_, double muf_) {
  
  m = m_; mur = mur_; muf = muf_;
  update();
  
  Vegas<AbstractHathor> integrator(*this, pdfmax);
  //integrator.setnprn(1);
  integrator.vegas();
  
  for (int i = 0; i < pdfmax; i++) {
    integral[i] = integrator.avgi_[i];
    error[i] = integrator.sd_[i];
    chi2a[i] = integrator.chi2a_[i];
  }
  
  // reset pdf member
  if (pdfmax > 1)
    pdf.InitMember(0);
  
  return integrator.avgi;
}

/*
 *  Internal routines for integration.
 */

double AbstractHathor::f(const double x[], const double wgt, double res[]) {
  /*
   *  Differential cross section including the parton luminosities.
   */
  
  // calculate kinematic variables
  setPartonicEnergy(x);
  
  if ((delta > 0) && (x2 + delta > 1.))
    return 0.;
  
  // evaluate partonic cross sections
  evaluateScalingFunctions();
  
  double hi[13], hj[13], hjleft[13], hjright[13];
  // loop over pdf set members
  for(int ipdf = 0; ipdf < pdfmax; ipdf++) {
    if (pdfmax > 1)
      pdf.InitMember(ipdf);

    pdf.GetPdf(x1, muf, hi);
    pdf.GetPdf(x2, muf, hj);

    if (delta > 0) {
      pdf.GetPdf(x2-delta, muf, hjleft);
      pdf.GetPdf(x2+delta, muf, hjright);
    }
    
    if (collider == PPBAR) {
      ChargeConjugation(hj);
      if (delta > 0) {
	ChargeConjugation(hjleft);
	ChargeConjugation(hjright);
      }
    }
    
    // evaluate parton fluxes
    if (delta > 0) {
      evaluatePDFs(hi, hj, hjleft, hjright);
    } else {
      evaluatePDFs(hi, hj);
    }
    
    // calculate full integral
    res[ipdf] = evaluateIntegral(as_pdf[ipdf], wgt);
  }
  
  return res[0];
}

void AbstractHathor::getPdfErr(double &up, double &down) {
  /*
   *  Calculate PDF uncertainties.
   */
  
  if (! (PDF_SCAN & scheme)) {
    err("Option PDF_SCAN is required to evaluate pdf uncertainties");
    exit(1);
  }
  
  if (evalerr != 0 ) {
    /* Use error function provided by user to evaluate pdf uncertainties */
    evalerr(pdfmax-1, integral, up, down);
    return;
  }
  
  up = 0; down = 0;
  
  unsigned int pdferr = 0;
  std::string pdfname = pdf.GetName();
  
  /* Determine how the pdf sets are to be interpreted: Symmetric
     or asymmetric uncertainties: */
  
  /* NNPDF scheme corresponds to symmetric errors */  
  if (pdfname.find("NNPDF") != pdfname.npos )
    pdferr = PDF_SYM_ERR | NNPDF; 
  
  /* ABKM/ABM gives symmetric errors */  
  if (pdfname.find("abkm") != pdfname.npos 
      || pdfname.find("abm") != pdfname.npos )
    pdferr = PDF_SYM_ERR; 
  
  /* Allow user to force a setting via scheme */
  pdferr = pdferr | (scheme & PDF_SYM_ERR) | (scheme & NNPDF);
  
  /* Symmetric error case */
  if (PDF_SYM_ERR & pdferr) {
    for (int i = 1; i < pdfmax; i++) {
      up += pow(integral[i]-integral[0], 2);
    }
    up = sqrt(up);
    
    if (NNPDF & pdferr) {
      up /= sqrt(pdfmax-1);
    }
    down = up;
    return;
  }
  
  /* asymmetric error case */
  double tmp;
  for (int i = 0; i < pdfmax/2; i++) {
    tmp = std::max(integral[2*i+1]-integral[0], integral[2*i+2]-integral[0]); 
    tmp = std::max(tmp,0.);
    up += tmp*tmp;
    tmp = std::max(integral[0]-integral[2*i+1], integral[0]-integral[2*i+2]); 
    tmp = std::max(tmp,0.);
    down += tmp*tmp;
  }
  up = sqrt(up); down = sqrt(down);
  
}

double  AbstractHathor::Li2(double x){
  /*
   * Calculate Li2 using the ff library. Note that only the real part
   * is calculated.
   */
  int ierr = 0;
  double xxli2,xxli1;
  if ( x >= 2 ){
    double res = M_PI*M_PI/3 - 0.5 * log(x)*log(x);
    x=1/x;
    ffxli2_(xxli2,xxli1,x,ierr);
    return(res-xxli2);
  }
  if ( x > 1. ){
    double res = M_PI*M_PI/6 - 0.5*log(x)*log(x) + log(1/x)*log(1.-1./x);
    x=1.-1./x;
    ffxli2_(xxli2,xxli1,x,ierr);
    return(res+xxli2);
  }
  if ( x > 0.5 ){
    double res = M_PI*M_PI/6 - log(x)*log(1.-x);
    x=1.-x;
    ffxli2_(xxli2,xxli1,x,ierr);
    return(res-xxli2);
  }
  if ( x < -1. ){
    double res = -M_PI*M_PI/6 - 0.5 * log(-x)*log(-x);
    x=1./x;
    ffxli2_(xxli2,xxli1,x,ierr);
    return(res-xxli2);
  }
  ffxli2_(xxli2,xxli1,x,ierr);
  return(xxli2);

}
