#include "SgTop.h"
#include <cmath>
#include <iomanip>

using namespace std;

/* CKM matrix as used in MCFM */ 
static const double ckmMCFM[3][3] =  {{0.97500, 0.22200, 0.00000},
				      {0.22200, 0.97500, 0.00000},
				      {0.00000, 0.00000, 1.00000}};

/* CKM as in PDG 12 */ 
static const double ckmPDG12[3][3] = {{0.97427, 0.22534, 0.00351},
				      {0.22520, 0.97344, 0.0412},
				      {0.00867, 0.0404,  0.999146}};

SgTop::SgTop(Pdf &pdf_) 
  : AbstractHathor(pdf_), 
    Vud(ckm[0][0]),Vus(ckm[0][1]),Vub(ckm[0][2]),
    Vcd(ckm[1][0]),Vcs(ckm[1][1]),Vcb(ckm[1][2]),
    Vtd(ckm[2][0]),Vts(ckm[2][1]),Vtb(ckm[2][2]) {

  // for single-top default to NLO for the time being
  scheme = LO | NLO;
  
  // default to top-quark production
  particle = TOPQUARK;

  setCkmMatrix(ckmPDG12);

  // default to NLL calculation, until c1 is known
  llflags = NLO_LL | NLO_NLL | NNLO_LL | NNLO_NLL | NNLO_NNLL | NNLO_NNNLL
    | NLL_CALCULATION;
  
  /*
   * Couplings in the G-scheme using mw = 80.385 GeV, mz=91.1876
   * and GF = 1.1663787e-5 as input, corresponds to MCFM using the
   * ewscheme = 1 and adjusting GF, mw,mz to the aforementioned values.
   */
  setAlpha(1. / 132.2332298);
  setSwq(0.2228972);
  update();

 info("Dear user,");
 info("the QCD corrections for single top-quark production used in this");
 info("code are based on MCFM [1-3] and a private implementation of the");
 info("formulae presented in [4].");
 info("In addition leading and sub-leading logarithms as given in [5-8]");
 info("are implemented."); 
 info("For the numerical integration Hathor uses an adopted version ");
 info("of the Vegas integrator [9,10]. The required random numbers are ");
 info("generated using Ranlux [11]. The dilogarithm is calculated using");
 info("the FF library. [12]");
 info("[1]  Campbell,Frederix, Maltoni and Tramontano,: t-channel single");
 info("     top production at hadron colliders,");
 info("     Phys. Rev. Lett. 102 (2009) 182003, [arXiv:0903.0005 [hep-ph]].");
 info("[2]  Campbell, Tramontano: Next-to-leading order corrections ");
 info("     to Wt production and decay,");
 info("     Nucl. Phys. B 726, 109 (2005) [arXiv:hep-ph/0506289");
 info("[3]  Campbell, Ellis and Tramontano: Single top production and ");
 info("     decay at next-to-leading order,");
 info("     Phys. Rev. D 70, 094012 (2004) [arXiv:hep-ph/0408158].");
 info("[4]  Harris, Laenen, Phaf, Sullivan, Weinzierl: The Fully ");
 info("     differential single top quark cross-section in next");
 info("     to leading order QCD,Phys.Rev. D66 (2002) 054024 [hep-ph/0207055]");
#if 0
 info("[5]  Kidonakis: Single top production at the Tevatron: Threshold");
 info("     resummation and finite-order soft gluon corrections,"); 
 info("     Phys.Rev. D74 (2006) 114012, [hep-ph/0609287]");
 info("[6]  Kidnoakis: NNLL resummation for s-channel single top quark");
 info("     production, Phys.Rev. D81 (2010) 054028, [arXiv:1001.5034]");
 info("[7]  Kidonakis: Two-loop soft anomalous dimensions for single top");
 info("     quark associated production with a W- or H-,");
 info("     Phys.Rev. D82 (2010) 054018 arXiv:1005.4451"); 
 info("[8]  Kidonakis: Next-to-next-to-leading-order collinear and soft");
 info("     gluon corrections for t-channel single top quark production,");
 info("     Phys.Rev. D83 (2011) 091503, [arXiv:1103.2792].");
#endif
 info("[5]  Lepage:A New Algorithm for Adaptive Multidimensional Integration");
 info("     J.Comput.Phys.27 (1978) 192");
 info("[6]  Lepage: VEGAS: AN ADAPTIVE MULTIDIMENSIONAL INTEGRATION PROGRAM");
 info("     CLNS-80/447 (1980).");
 info("[7]  Luescher: A Portable high quality random number generator");
 info("     for lattice field theory simulations,");
 info("     Comput.Phys.Commun.79 (1994) 100 [hep-lat/930902]");
 info("[8]  van Oldenborgh: FF: A Package to evaluate one loop");
 info("     Feynman diagrams, Comput.Phys.Commun 66 (1991) 1");
 info("");
 info("Please cite the corresponding articles if you quote Hathor results");
 info("in your article relying on the aformentioned calculations.");
}

void SgTop::update() {
  /*
   *  Update constants and clear results from previous integration.
   */
  AbstractHathor::update();
  
  m4 = m2*m2;
  
  if ( (scheme & NNLO) && ( (mur != m) || (muf != m) ) ) {
    // NNLO scales need two additional integration dimensions
    dimension = 4;
  } else if (scheme & (NLOAPPROX | NNLOAPPROX)) {
    // approx (N)NLO always needs two additional integration dimensions
    dimension = 4;
  } else if ((scheme & NLO) && ((mur != m) || (muf != m))) {
    // (N)NLO scale dependencies need one additional dimension
    dimension = 3;
  } else {
    // for (N)LO and mu=m, run the usual two dimensional integration
    dimension = 2;
  }
  
  // check settings
  if ((scheme & NLO) && (scheme & NLOAPPROX)) {
    err("The NLO and NLOAPPROX calculations are mutually exclusive.");
    err("Please choose only one NLO contribution.");
    exit(1);
  }

  if ((scheme & NNLO) && (scheme & NNLOAPPROX)) {
    err("The NNLO and NNLOAPPROX calculations are mutually exclusive.");
    err("Please choose only one NNLO contribution.");
    exit(1);
  }
  
  // reset contribution counters
  for (int i = 0; i < 6; i++)
    contribution[i] = 0.;
}

void SgTop::setLLflags(const unsigned int flags) {
  llflags = flags;
  if (debug) {
    cout << " HATHOR: Chosen leading logs: ";
    if (llflags & NLO_LL)   cout << "NLO-LL ";
    if (llflags & NLO_NLL)  cout << "NLO-NLL ";
    if (llflags & NNLO_LL)  cout << "NNLO-LL ";
    if (llflags & NNLO_NLL) cout << "NNLO-NLL ";
    if (llflags & NNLO_NNLL)  cout << "NNLO-NNLL ";
    if (llflags & NNLO_NNNLL) cout << "NNLO-NNNLL ";
    if (llflags & NLL_CALCULATION) cout << " (NLL from hep-ph/0609287) ";
    cout << endl;
  }
}

void SgTop::setAlpha(const double alpha_) {
  /*
   * Set electro-weak alpha(M_Z) (e.g. 1/132)
   */
  alpha  = alpha_;
  alpha2 = alpha * alpha;
}

void SgTop::setSwq(const double value) {
  /*
   * Set weak coupling, sin(theta_W)^2.
   */
  swq  = value;
}

void SgTop::setParticle(PARTICLE particle_) {
  /*
   *  Set desired result: top, antitop, top+antitop.
   */
  particle = particle_;
  if (particle_ == TOPQUARK) {
    info("Calculating top-quark production");
    return;
  }
  if (particle_ == ANTITOPQUARK) {
    info("Calculating antitop-quark production");
    return;
  }
  if (particle_ == BOTH) {
    info("Calculating sum of top-quark and antitop-quark production");
    return;
  }
  
  err("Unknown particle given to setParticle().");
  exit(1);
}

void SgTop::setCkmMatrix(const double ckm_[3][3]) {
  /* we dont test the unitarity */
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      ckm[i][j] = ckm_[i][j];
}

void SgTop::getCkmMatrix(double ckm_[3][3]) {
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      ckm_[i][j] = ckm[i][j];
}

void SgTop::PrintCkmMatrix() {
  /*
   *  Print the CKM matrix.
   */
  info("The following CKM matrix is used:");
  cout << setw(9) << Vud << setw(9) << Vus << setw(9) << Vub << endl;
  cout << setw(9) << Vcd << setw(9) << Vcs << setw(9) << Vcb << endl;
  cout << setw(9) << Vtd << setw(9) << Vts << setw(9) << Vtb << endl;
}

void SgTop::PrintContributions() {
  /* 
   *  Print contributions from different initial states.
   */
  
  if (particle == BOTH) {
    warn("The contribution from different initial states can not be");
    warn("computed with setParticle(SgTop::BOTH).");
    warn("Please specify either top or antitop production to use this feature.");
    return;
  }
  
  double contrib_norm = 0.;
  for (int i = 0; i < 6; i++)
    contrib_norm += contribution[i];
  
  cout << "    QQ: " << setw(10) 
       << (contribution[0] / contrib_norm * integral[0]) << " pb " 
       << "(" << setw(7) << setprecision(2) << fixed 
       << (contribution[0] / contrib_norm * 100.) << " %)" << endl;
  cout << " QQBAR: " << setw(10) 
       << (contribution[1] / contrib_norm * integral[0]) << " pb "
       << "(" << setw(7) << setprecision(2) << fixed 
       << (contribution[1] / contrib_norm * 100.) << " %)" << endl;
  cout << "    GQ: " << setw(10) 
       << (contribution[2] / contrib_norm * integral[0]) << " pb "
       << "(" << setw(7) << setprecision(2) << fixed 
       << (contribution[2] / contrib_norm * 100.) << " %)" << endl;
  cout << " GQBAR: " << setw(10) 
       << (contribution[3] / contrib_norm * integral[0]) << " pb "
       << "(" << setw(7) << setprecision(2) << fixed 
       << (contribution[3] / contrib_norm * 100.) << " %)" << endl;
  cout << "    GG: " << setw(10) 
       << (contribution[4] / contrib_norm * integral[0]) << " pb "
       << "(" << setw(7) << setprecision(2) << fixed 
       << (contribution[4] / contrib_norm * 100.) << " %)" << endl;
  cout << "  QBQB: " << setw(10) 
       << (contribution[5] / contrib_norm * integral[0]) << " pb "
       << "(" << setw(7) << setprecision(2) << fixed 
       << (contribution[5] / contrib_norm * 100.) << " %)" << endl;
}
