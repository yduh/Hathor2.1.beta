#include "Hathor.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>

using namespace std;

/*
 * This example code shows, how to get a cross section.
 */

void SgTopExample() {
  // set mu_r and mu_f to desired scales
  //-double mt = 173.5, mur = mt*1., muf = mt*1.;
  double mt = 173.5, mur = mt*1.7, muf = mt*1.5;
  // enable LO and NLO contributions
  unsigned int scheme = Hathor::LO | Hathor::NLO;
  // choose pdf string
  string pdfname = "MSTW2008nnlo68cl";

  // initialize LHApdf and Hathor

  double val, err;
  Lhapdf pdf(pdfname);
  //-  HathorSgTopS XS(pdf);
  HathorSgTopT XS(pdf);
  XS.sethc2(0.38937911e9);	    // value used by MCFM
  XS.setColliderType(Hathor::PP);   // choose PP or PPBAR
  XS.setSqrtShad(8000.);	    // energy in GeV
  XS.setParticle(SgTop::TOPQUARK);  // choose TOPQUARK, ANTITOPQUARK or BOTH
  XS.setScheme(scheme);		    // see above
  XS.setPrecision(Hathor::HIGH);    // choose LOW, MEDIUM or HIGH



  // EXAMPLE: modify CKM matrix by setting V_tb = 0.9999
  XS.PrintCkmMatrix();
  double ckm[3][3];
  XS.getCkmMatrix(ckm);
  ckm[2][2] = 0.9999;
  XS.setCkmMatrix(ckm);
  XS.PrintCkmMatrix();

  // calculate cross section
  XS.getXsection(mt, mur, muf);

  // fetch/print result
  XS.getResult(0, val, err);	     // evaluate cross section!

  cout << endl;
  cout << "      alpha_s = " << XS.getAlphas(mur) << endl;
  cout << "cross section = " << val << " +/- " << err << " pb" << endl;
  cout << endl;
}

int main() {
  SgTopExample();
}

