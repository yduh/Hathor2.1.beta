// $Modified: Thu Jan 30 18:16:07 2014 by puwer $
#include "SgTopTChannel.h"
#include "SingleTopGrids.h"

using namespace std;

HathorSgTopT::HathorSgTopT(Pdf &pdf_) 
  : SgTop(pdf_),
    pxs_qq_10(t_sqrts_ub,t_xs_ub, massval),
    pxs_qqb_10(t_sqrts_dbarb,t_xs_dbarb, massval),
    pxs_gu_10(t_sqrts_gu,t_xs_gu, massval),
    pxs_gb_10(t_sqrts_gb,t_xs_gb, massval),
    pxs_gqb_10(t_sqrts_gdbar,t_xs_gdbar, massval),
    xPqq0(SplittingFunction::PQQ0),
    xPqg0(SplittingFunction::PQG0),
    xPgg0(SplittingFunction::PGG0),
    xPgq0(SplittingFunction::PGQ0),
    xPqg1(SplittingFunction::PQG1),
    xPqq1s(SplittingFunction::PQQ1S),
    xPqq1v(SplittingFunction::PQQ1V),
  xPqqb1v(SplittingFunction::PQQB1V),
  xPqggg(SplittingFunction::PQG0,SplittingFunction::PGG0),
  xPqggq(SplittingFunction::PQG0,SplittingFunction::PGQ0),
  xPqgqg(SplittingFunction::PQG0,SplittingFunction::PQG0),
  xPqqqg(SplittingFunction::PQQ0,SplittingFunction::PQG0),
  xPqqqq(SplittingFunction::PQQ0,SplittingFunction::PQQ0){

  double scalefac = 1. / fourpi * (MCFMswq*MCFMswq) / (MCFMalpha*MCFMalpha);

  pxs_qq_10.setScalefactor(scalefac);
  pxs_qqb_10.setScalefactor(scalefac);
  pxs_gu_10.setScalefactor(scalefac);
  pxs_gb_10.setScalefactor(scalefac/4.);
  pxs_gqb_10.setScalefactor(scalefac);

}

void HathorSgTopT::update() {
  /*
   * The user has changed the setting for the top mass.
   * We have to update the mass parameter and the production threshold
   * in interpolation routine. 
   */
  SgTop::update();

  pxs_qq_00.setMass(m, m);
  pxs_qqb_00.setMass(m, m);
  pxs_qq_10.setMass(m, m);
  pxs_qqb_10.setMass(m, m);
  pxs_gu_10.setMass(m, m);
  pxs_gb_10.setMass(m, m);
  pxs_gqb_10.setMass(m, m);
  
  fqq00  = fqqbar00  = 0.;
  nloqq  = nloqqbar  = nlogb  = nlogu  = nlogqbar  = 0.;
  nnloqq = nnloqqbar = nnlogb = nnlogu = nnlogqbar = nnlogg =
    nnlouqig = nnlodqig = nnlodbarqig = nnloqqbv = nnloqqv = 0.;
}

void HathorSgTopT::setPartonicEnergy(const double x[]) {
  /*
   *  Calculate kinematic variables.
   */
  
  // exclude energies below kinematic threshold
  sparton = m2*pow(shad/m2,x[0]);
  double r = sparton/shad;
  x1 = r*pow(r,-x[1]);
  x2 = r/x1;
  jacobian = -r*log(shad/m2)*log(r); 
  
  sqrts = sqrt(sparton);
  rho   = m2/sparton;
  
  rt  = (dimension > 2) ? x[2] : 0.5;
  rs4 = (dimension > 3) ? x[3] : 0.5;
}

void HathorSgTopT::evaluatePDFs(const double h1[], const double h2[],
				const double h2left[], 
				const double h2right[]) {
  /*
   *  Evaluate parton fluxes.
   */

  double ckm_u = Vud*Vud + Vus*Vus + Vub*Vub;
  double ckm_c = Vcd*Vcd + Vcs*Vcs + Vcb*Vcb;
  double ckm_t = Vtd*Vtd + Vts*Vts + Vtb*Vtb;
  double ckm_d = Vud*Vud + Vcd*Vcd;
  double ckm_s = Vus*Vus + Vcs*Vcs;
  double ckm_b = Vub*Vub + Vcb*Vcb;
  
  double t_bu, t_bdbar, t_bg, t_gdbar, t_gu;
  double tbar_bu, tbar_bdbar, tbar_bg, tbar_gdbar, tbar_gu;
  double t_gg, t_uqig, t_dqig, t_dbarqig, t_qqbv, t_qqv;
  double tbar_gg, tbar_uqig, tbar_dqig, tbar_dbarqig, tbar_qqbv, tbar_qqv;
  
  // required for leading order
  t_bu = 
    + (h1[DOWN] * Vtd*Vtd + h1[STRANGE] * Vts*Vts + h1[BOTTOM] * Vtb*Vtb) 
    * (h2[UP] * ckm_u + h2[CHARM] * ckm_c)
    + (h2[DOWN] * Vtd*Vtd + h2[STRANGE] * Vts*Vts + h2[BOTTOM] * Vtb*Vtb) 
    * (h1[UP] * ckm_u + h1[CHARM] * ckm_c);

  tbar_bu = 
    + (h1[ADOWN] * Vtd*Vtd + h1[ASTRANGE] * Vts*Vts + h1[ABOTTOM] * Vtb*Vtb) 
    * (h2[AUP] * ckm_u + h2[ACHARM] * ckm_c)
    + (h2[ADOWN] * Vtd*Vtd + h2[ASTRANGE] * Vts*Vts + h2[ABOTTOM] * Vtb*Vtb) 
    * (h1[AUP] * ckm_u + h1[ACHARM] * ckm_c);
  
  t_bdbar = 
    + (h1[DOWN] * Vtd*Vtd + h1[STRANGE] * Vts*Vts + h1[BOTTOM] * Vtb*Vtb) 
    * (h2[ADOWN] * ckm_d + h2[ASTRANGE] * ckm_s + h2[ABOTTOM] * ckm_b)
    + (h2[DOWN] * Vtd*Vtd + h2[STRANGE] * Vts*Vts + h2[BOTTOM] * Vtb*Vtb) 
    * (h1[ADOWN] * ckm_d + h1[ASTRANGE] * ckm_s + h1[ABOTTOM] * ckm_b);

  tbar_bdbar = 
    + (h1[ADOWN] * Vtd*Vtd + h1[ASTRANGE] * Vts*Vts + h1[ABOTTOM] * Vtb*Vtb) 
    * (h2[DOWN] * ckm_d + h2[STRANGE] * ckm_s + h2[BOTTOM] * ckm_b)
    + (h2[ADOWN] * Vtd*Vtd + h2[ASTRANGE] * Vts*Vts + h2[ABOTTOM] * Vtb*Vtb) 
    * (h1[DOWN] * ckm_d + h1[STRANGE] * ckm_s + h1[BOTTOM] * ckm_b);
  
  // required for next-to-leading order
  t_bg = 
    + (h1[GLUON] * (h2[BOTTOM] * Vtb*Vtb + h2[STRANGE] * Vts*Vts 
		    + h2[DOWN] * Vtd*Vtd) 
       * (ckm_u + ckm_c + ckm_d + ckm_s + ckm_b))
    + (h2[GLUON] * (h1[BOTTOM] * Vtb*Vtb + h1[STRANGE] * Vts*Vts 
		    + h1[DOWN] * Vtd*Vtd) 
       * (ckm_u + ckm_c + ckm_d + ckm_s + ckm_b));

  tbar_bg = 
    + (h1[GLUON] * (h2[ABOTTOM] * Vtb*Vtb + h2[ASTRANGE] * Vts*Vts 
		    + h2[ADOWN] * Vtd*Vtd) 
       * (ckm_u + ckm_c + ckm_d + ckm_s + ckm_b))
    + (h2[GLUON] * (h1[ABOTTOM] * Vtb*Vtb + h1[ASTRANGE] * Vts*Vts 
		    + h1[ADOWN] * Vtd*Vtd) 
       * (ckm_u + ckm_c + ckm_d + ckm_s + ckm_b));

  t_gdbar = (h1[GLUON] * (h2[ADOWN] * ckm_d + h2[ASTRANGE] * ckm_s 
			  + h2[ABOTTOM] * ckm_b)) * ckm_t
    + (h2[GLUON] * (h1[ADOWN] * ckm_d + h1[ASTRANGE] * ckm_s 
		    + h1[ABOTTOM] * ckm_b)) * ckm_t;
  
  tbar_gdbar = 
    + (h1[GLUON] * (h2[DOWN] * ckm_d + h2[STRANGE] * ckm_s 
		    + h2[BOTTOM] * ckm_b)) * ckm_t
    + (h2[GLUON] * (h1[DOWN] * ckm_d + h1[STRANGE] * ckm_s 
		    + h1[BOTTOM] * ckm_b)) * ckm_t;

  t_gu  = 
    + (h1[GLUON] * (h2[UP] * ckm_u + h2[CHARM] * ckm_c)) * ckm_t 
    + (h2[GLUON] * (h1[UP] * ckm_u + h1[CHARM] * ckm_c)) * ckm_t;
  
  tbar_gu    = 
    + (h1[GLUON] * (h2[AUP] * ckm_u + h2[ACHARM] * ckm_c)) * ckm_t 
    + (h2[GLUON] * (h1[AUP] * ckm_u + h1[ACHARM] * ckm_c)) * ckm_t;
  
  // required for next-to-next-to-leading order
  t_gg    = 
    + (h1[GLUON] * h2[GLUON]) * ckm_t * (ckm_c + ckm_u); //(ckm_d + ckm_s + ckm_b + ckm_u + ckm_c);

  tbar_gg = t_gg;

  t_uqig  = 
    + (h1[UP] * ckm_u + h1[CHARM] * ckm_c) 
    * (h2[UP] + h2[CHARM] + h2[DOWN] + h2[STRANGE] + h2[BOTTOM] 
       + h2[AUP] + h2[ACHARM] + h2[ADOWN] + h2[ASTRANGE] + h2[ABOTTOM]) * ckm_t
    + (h2[UP] * ckm_u + h2[CHARM] * ckm_c) 
    * (h1[UP] + h1[CHARM] + h1[DOWN] + h1[STRANGE] + h1[BOTTOM] + h1[AUP] 
       + h1[ACHARM] + h1[ADOWN] + h1[ASTRANGE] + h1[ABOTTOM]) * ckm_t;
  
  tbar_uqig = 
    + (h1[AUP] * ckm_u + h1[ACHARM] * ckm_c) 
    * (h2[UP] + h2[CHARM] + h2[DOWN] + h2[STRANGE] + h2[BOTTOM] + h2[AUP] 
       + h2[ACHARM] + h2[ADOWN] + h2[ASTRANGE] + h2[ABOTTOM]) * ckm_t
    + (h2[AUP] * ckm_u + h2[ACHARM] * ckm_c) 
    * (h1[UP] + h1[CHARM] + h1[DOWN] + h1[STRANGE] + h1[BOTTOM] 
       + h1[AUP] + h1[ACHARM] + h1[ADOWN] + h1[ASTRANGE] + h1[ABOTTOM]) * ckm_t;

  t_dbarqig = 
    + (h1[ADOWN] * ckm_d + h1[ASTRANGE] * ckm_s + h1[ABOTTOM] * ckm_b) 
    * (h2[UP] + h2[CHARM] + h2[DOWN] + h2[STRANGE] + h2[BOTTOM] + h2[AUP] 
       + h2[ACHARM] + h2[ADOWN] + h2[ASTRANGE] + h2[ABOTTOM]) * ckm_t
    + (h2[ADOWN] * ckm_d + h2[ASTRANGE] * ckm_s + h2[ABOTTOM] * ckm_b) 
    * (h1[UP] + h1[CHARM] + h1[DOWN] + h1[STRANGE] + h1[BOTTOM] + h1[AUP] 
       + h1[ACHARM] + h1[ADOWN] + h1[ASTRANGE] + h1[ABOTTOM]) * ckm_t;

  tbar_dbarqig = 
    + (h1[DOWN] * ckm_d + h1[STRANGE] * ckm_s + h1[BOTTOM] * ckm_b) 
    * (h2[UP] + h2[CHARM] + h2[DOWN] + h2[STRANGE] + h2[BOTTOM] 
       + h2[AUP] + h2[ACHARM] + h2[ADOWN] + h2[ASTRANGE] + h2[ABOTTOM]) * ckm_t
    + (h2[DOWN] * ckm_d + h2[STRANGE] * ckm_s + h2[BOTTOM] * ckm_b) 
    * (h1[UP] + h1[CHARM] + h1[DOWN] + h1[STRANGE] + h1[BOTTOM] + h1[AUP] 
       + h1[ACHARM] + h1[ADOWN] + h1[ASTRANGE] + h1[ABOTTOM]) * ckm_t;

  t_dqig = 
    + (h1[DOWN] * Vtd*Vtd + h1[STRANGE] * Vts*Vts + h1[BOTTOM] * Vtb*Vtb) 
    * (h2[UP] + h2[CHARM] + h2[DOWN] + h2[STRANGE] + h2[BOTTOM] + h2[AUP] 
       + h2[ACHARM] + h2[ADOWN] + h2[ASTRANGE] + h2[ABOTTOM]) 
    * (ckm_u + ckm_c) * 2.
    + (h2[DOWN] * Vtd*Vtd + h2[STRANGE] * Vts*Vts + h2[BOTTOM] * Vtb*Vtb) 
    * (h1[UP] + h1[CHARM] + h1[DOWN] + h1[STRANGE] + h1[BOTTOM] + h1[AUP] 
       + h1[ACHARM] + h1[ADOWN] + h1[ASTRANGE] + h1[ABOTTOM]) 
    * (ckm_u + ckm_c) * 2.;

  tbar_dqig = 
    + (h1[ADOWN] * Vtd*Vtd + h1[ASTRANGE] * Vts*Vts + h1[ABOTTOM] * Vtb*Vtb) 
    * (h2[UP] + h2[CHARM] + h2[DOWN] + h2[STRANGE] + h2[BOTTOM] + h2[AUP] 
       + h2[ACHARM] + h2[ADOWN] + h2[ASTRANGE] + h2[ABOTTOM]) 
    * (ckm_u + ckm_c) * 2.
    + (h2[ADOWN] * Vtd*Vtd + h2[ASTRANGE] * Vts*Vts + h2[ABOTTOM] * Vtb*Vtb) 
    * (h1[UP] + h1[CHARM] + h1[DOWN] + h1[STRANGE] + h1[BOTTOM] + h1[AUP] 
       + h1[ACHARM] + h1[ADOWN] + h1[ASTRANGE] + h1[ABOTTOM]) 
    * (ckm_u + ckm_c) * 2.;

  t_qqbv = 
    + (h1[UP] * ckm_u + h1[CHARM] * ckm_c) 
    * (h2[ADOWN]*Vtd*Vtd + h2[ASTRANGE]*Vts*Vts + h2[ABOTTOM]*Vtb*Vtb)
    + (h1[DOWN]*Vtd*Vtd + h1[STRANGE]*Vts*Vts + h1[BOTTOM]*Vtb*Vtb) 
    * (h2[AUP]*ckm_u + h2[ACHARM]*ckm_c)
    + (h2[UP] * ckm_u + h2[CHARM] * ckm_c) 
    * (h1[ADOWN]*Vtd*Vtd + h1[ASTRANGE]*Vts*Vts + h1[ABOTTOM]*Vtb*Vtb)
    + (h2[DOWN]*Vtd*Vtd + h2[STRANGE]*Vts*Vts + h2[BOTTOM]*Vtb*Vtb) 
    * (h1[AUP]*ckm_u + h1[ACHARM]*ckm_c);

  tbar_qqbv = 
    + (h1[UP] * ckm_u + h1[CHARM] * ckm_c) 
    * (h2[ADOWN]*Vtd*Vtd + h2[ASTRANGE]*Vts*Vts + h2[ABOTTOM]*Vtb*Vtb)
    + (h1[DOWN]*Vtd*Vtd + h1[STRANGE]*Vts*Vts + h1[BOTTOM]*Vtb*Vtb) 
    * (h2[AUP]*ckm_u + h2[ACHARM]*ckm_c)
    + (h2[UP] * ckm_u + h2[CHARM] * ckm_c) 
    * (h1[ADOWN]*Vtd*Vtd + h1[ASTRANGE]*Vts*Vts + h1[ABOTTOM]*Vtb*Vtb)
    + (h2[DOWN]*Vtd*Vtd + h2[STRANGE]*Vts*Vts + h2[BOTTOM]*Vtb*Vtb) 
    * (h1[AUP]*ckm_u + h1[ACHARM]*ckm_c);

  t_qqv =  
    + (h1[DOWN]*ckm_d  + h1[STRANGE]*ckm_s  + h1[BOTTOM]*ckm_b)  
    * (h2[DOWN]*Vtd*Vtd  + h2[STRANGE]*Vts*Vts  + h2[BOTTOM]*Vtb*Vtb)
    + (h1[ADOWN]*ckm_d + h1[ASTRANGE]*ckm_s + h1[ABOTTOM]*ckm_b) 
    * (h2[ADOWN]*Vtd*Vtd + h2[ASTRANGE]*Vts*Vts + h2[ABOTTOM]*Vtb*Vtb)
    + (h2[ADOWN]*ckm_d + h2[ASTRANGE]*ckm_s + h2[ABOTTOM]*ckm_b) 
    * (h1[ADOWN]*Vtd*Vtd + h1[ASTRANGE]*Vts*Vts + h1[ABOTTOM]*Vtb*Vtb)
    + (h2[DOWN]*ckm_d  + h2[STRANGE]*ckm_s  + h2[BOTTOM]*ckm_b)  
    * (h1[DOWN]*Vtd*Vtd  + h1[STRANGE]*Vts*Vts  + h1[BOTTOM]*Vtb*Vtb);

  tbar_qqv =  
    + (h1[DOWN]*ckm_d  + h1[STRANGE]*ckm_s  + h1[BOTTOM]*ckm_b)  
    * (h2[DOWN]*Vtd*Vtd  + h2[STRANGE]*Vts*Vts  + h2[BOTTOM]*Vtb*Vtb)
    + (h1[ADOWN]*ckm_d + h1[ASTRANGE]*ckm_s + h1[ABOTTOM]*ckm_b) 
    * (h2[ADOWN]*Vtd*Vtd + h2[ASTRANGE]*Vts*Vts + h2[ABOTTOM]*Vtb*Vtb)
    + (h2[ADOWN]*ckm_d + h2[ASTRANGE]*ckm_s + h2[ABOTTOM]*ckm_b) 
    * (h1[ADOWN]*Vtd*Vtd + h1[ASTRANGE]*Vts*Vts + h1[ABOTTOM]*Vtb*Vtb)
    + (h2[DOWN]*ckm_d  + h2[STRANGE]*ckm_s  + h2[BOTTOM]*ckm_b)  
    * (h1[DOWN]*Vtd*Vtd  + h1[STRANGE]*Vts*Vts  + h1[BOTTOM]*Vtb*Vtb);
  
  // set fluxes according to chosen particle
  if (particle == TOPQUARK) {
    qq = t_bu; qqbar = t_bdbar;
    gb = t_bg; gqbar = t_gdbar; gu = t_gu;
    gg = t_gg;
    uq_ig = t_uqig;
    dq_ig = t_dqig;
    dbarq_ig = t_dbarqig;
    qqbv = t_qqbv; qqv = t_qqv;
    return;
  } 
  
  if (particle == ANTITOPQUARK) {
    qq = tbar_bu; qqbar = tbar_bdbar;
    gb = tbar_bg; gqbar = tbar_gdbar; gu = tbar_gu;
    gg = tbar_gg;
    uq_ig = tbar_uqig;
    dq_ig = tbar_dqig; dbarq_ig = tbar_dbarqig; 
    qqbv = tbar_qqbv; qqv = tbar_qqv;
    return;
  }
  
  if (particle == BOTH) {
    qq = t_bu + tbar_bu; qqbar = t_bdbar + tbar_bdbar;
    gb = t_bg + tbar_bg; gqbar = t_gdbar + tbar_gdbar; gu = t_gu + tbar_gu;
    gg = t_gg + tbar_gg;
    uq_ig = t_uqig + tbar_uqig;
    dq_ig = t_dqig + tbar_dqig; dbarq_ig = t_dbarqig + tbar_dbarqig;
    qqbv = t_qqbv + tbar_qqbv; qqv = t_qqv + tbar_qqv;
    return;
  }

  cout << "# set particle to TOPQUARK, ANTITOPQUAR or BOTH \n";
  exit(1);

}

void HathorSgTopT::evaluateScalingFunctions() {
  /*
   *  Evaluate scaling functions for further use.
   */

  // Leading order
  fqq00    = pxs_qq_00.eval(sparton);
  fqqbar00 = pxs_qqb_00.eval(sparton);
  
  if ( scheme & ( NLOAPPROX | NNLOAPPROX ) )
    evaluateApproxNNLO();
  
  // Next-to-leading order
  double fqqbar10 = 0., fqq10 = 0., fgb10 = 0., fgqbar10 = 0., fgu10 = 0.;
  double fqqbar11 = 0., fqq11 = 0., fgb11 = 0., fgqbar11 = 0., fgu11 = 0.;

  if ( scheme & ( NLO | NNLO | NLOAPPROX | NNLOAPPROX ) ) {
    // NLO scaling functions to calculate mu-dependence
    fqq10    = pxs_qq_10.eval(sparton);
    fqqbar10 = pxs_qqb_10.eval(sparton);
    fgu10    = pxs_gu_10.eval(sparton);
    fgb10    = pxs_gb_10.eval(sparton);
    fgqbar10 = pxs_gqb_10.eval(sparton);
    fqq11    = xs_qq_11();
    fqqbar11 = xs_qqbar_11();
    fgb11    = xs_gb_11();
    fgqbar11 = xs_gqbar_11();
    fgu11    = xs_gu_11();
  }

  /* 
   * Combine the scale dependent contributions with the scale
   * independent contributions for the two different approximations:
   */
  if ( scheme & NLO ) {
    nloqqbar = fourpi * (fqqbar10 + Lmf * fqqbar11);
    nloqq    = fourpi * (fqq10    + Lmf * fqq11   );
    nlogb    = fourpi * (fgb10    + Lmf * fgb11   );
    nlogqbar = fourpi * (fgqbar10 + Lmf * fgqbar11);
    nlogu    = fourpi * (fgu10    + Lmf * fgu11   );
  } else if (scheme & NLOAPPROX) {
    nloqqbar = fourpi * (         + Lmf * fqqbar11) + approxqqbar1;
    nloqq    = fourpi * (         + Lmf * fqq11   ) + approxqq1;
    nlogb    = fourpi * (         + Lmf * fgb11   );
    nlogqbar = fourpi * (         + Lmf * fgqbar11);
    nlogu    = fourpi * (         + Lmf * fgu11   );
  }
  
  
  if ( scheme & FAKE_NNLO ) {
    /* This option is used to estimate what can be expected from the 
     * full NNLO calcu;ation as far as the scaledependence is concerned.
     * We estimate the size of the NNLO corrections for mu=mt using
     * a simple ansatz:
     */
    // - scale NLO result with ratio (NLO/LO)
    approxqq2    = fourpi * fqq10    * (fqq10    / fqq00);
    approxqqbar2 = fourpi * fqqbar10 * (fqqbar10 / fqqbar00);
    // - simply use LO result
    //approxqq2    = fqq00;
    //approxqqbar2 = fqqbar00;
    // - model with sqrt(mt^2 / s)
    //approxqq2    = sqrt(rho / hcq) / (twopi*twopi);
    //approxqqbar2 = sqrt(rho / hcq) / (twopi*twopi);
    
  }
  
  /*
   * Calculate the full scale dependence
   * Whatever approximation is used we always include the complete
   * scale (in)dependence!
   */
  if ( scheme & ( NNLO | NNLOAPPROX | FAKE_NNLO ) ) {

    double fqq21, fqqbar21, fgb21, fgqbar21, fgu21, fgg21, fuqig21, 
      fdqig21, fdbarqig21, fqqbv21, fqqv21, fqq22, fqqbar22, fgb22, 
      fgqbar22, fgu22, fgg22, fuqig22,fdqig22, fdbarqig22, fqqbv22, fqqv22;

    fqq21    = xs_qq_21();
    fqqbar21 = xs_qqbar_21();
    fgb21    = xs_gb_21();
    fgqbar21 = xs_gqbar_21();
    fgu21    = xs_gu_21();
    fgg21    = xs_gg_21();
    fgg22    = xs_gg_22();
    fuqig21  = xs_uqig_21();
    fdqig21  = xs_dqig_21();
    fdbarqig21 = xs_dbarqig_21();
    fdbarqig22 = xs_dbarqig_22();
    fqq22    = xs_qq_22();    // these function calls are ordered 
    fqqbar22 = xs_qqbar_22(); // to increase the performance
    fgu22    = xs_gu_22();
    fgb22    = xs_gb_22();
    fgqbar22 = xs_gqbar_22();
    fuqig22  = xs_uqig_22();
    fdqig22  = xs_dqig_22();
    fqqbv21  = xs_qqbv_21();
    fqqbv22  = xs_qqbv_22();
    fqqv21   = xs_qqv_21();
    fqqv22   = xs_qqv_22();

    nnloqq       = fourpi2 * Lmf * (fqq21 + Lmf * fqq22)
      + fourpi2 * Lmr * (fqq10 + Lmf * fqq11) * b0;
    nnloqqbar    = fourpi2 * Lmf * (fqqbar21 + Lmf * fqqbar22)
      + fourpi2 * Lmr * (fqqbar10 + Lmf * fqqbar11) * b0;
    nnlogb       = fourpi2 * Lmf * (fgb21 + Lmf * fgb22)
      + fourpi2 * Lmr * (fgb10 + Lmf * fgb11) * b0;
    nnlogqbar    = fourpi2 * Lmf * (fgqbar21 + Lmf * fgqbar22)
      + fourpi2 * Lmr * (fgqbar10 + Lmf * fgqbar11) * b0;
    nnlogu       = fourpi2 * Lmf * (fgu21 + Lmf * fgu22)
      + fourpi2 * Lmr * (fgu10 + Lmf * fgu11) * b0;
    nnlogg       = fourpi2 * Lmf * (fgg21 + Lmf * fgg22);
    nnlouqig     = fourpi2 * Lmf * (fuqig21 + Lmf * fuqig22);
    nnlodqig     = fourpi2 * Lmf * (fdqig21 + Lmf * fdqig22);
    nnlodbarqig  = fourpi2 * Lmf * (fdbarqig21 + Lmf * fdbarqig22);
    nnloqqbv     = fourpi2 * Lmf * (fqqbv21 + Lmf * fqqbv22);
    nnloqqv      = fourpi2 * Lmf * (fqqv21 + Lmf * fqqv22);
    // Add scale-independent approx. NNLO terms
    if (scheme & (NNLOAPPROX | FAKE_NNLO)) {
      nnloqq    += approxqq2;
      nnloqqbar += approxqqbar2;
    }
  }
}

double HathorSgTopT::evaluateIntegral(double as, double wgt) {
  /*
   *  Calculate the full integrand, using the given alpha_s.
   */
  
  double lo = 0, nlo = 0, nnlo = 0;
  const double as2 = as * as;
  
  double fac = jacobian * hcq * wgt;

  if (scheme & LO) {
    lo = qq * fqq00 + qqbar * fqqbar00;
    if (particle == TOPQUARK) {
      contribution[0] += fac * qq * fqq00;       // q q
      contribution[1] += fac * qqbar * fqqbar00; // q qbar
    } else if (particle == ANTITOPQUARK) {
      contribution[5] += fac * qq * fqq00;       // qbar qbar
      contribution[1] += fac * qqbar * fqqbar00; // qbar q
    }
  }
  
  fac *= as;

  if (scheme & (NLO | NLOAPPROX)) {
    nlo = qq * nloqq + qqbar * nloqqbar + gb * nlogb 
      + gqbar * nlogqbar + gu * nlogu;
    if (particle == TOPQUARK) {
      contribution[0] += fac * qq * nloqq;                // q q
      contribution[1] += fac * qqbar * nloqqbar;          // qbar q
      contribution[2] += fac * (gu * nlogu + gb * nlogb); // g q
      contribution[3] += fac * gqbar * nlogqbar;          // g qbar
    } else if (particle == ANTITOPQUARK) {
      contribution[5] += fac * qq * nloqq;                // qbar qbar
      contribution[1] += fac * qqbar * nloqqbar;          // q qbar
      contribution[3] += fac * (gu * nlogu + gb * nlogb); // g qbar
      contribution[2] += fac * gqbar * nlogqbar;          // g q
    }
  }
  
  fac *= as;

  if (scheme & (NNLO | NNLOAPPROX | FAKE_NNLO)) {
    nnlo = qq * nnloqq + qqbar * nnloqqbar + gb * nnlogb 
      + gqbar * nnlogqbar + gu * nnlogu
      + gg * nnlogg + uq_ig * nnlouqig + dq_ig * nnlodqig 
      + dbarq_ig * nnlodbarqig
      + qqbv * nnloqqbv + qqv * nnloqqv;
    if (particle == TOPQUARK) {
      contribution[0] += fac * qq * nnloqq;		// q q
      contribution[1] += fac * qqbar * nnloqqbar;	// q qbar
      contribution[2] += fac * (gb * nnlogb + gu * nnlogu); // g q
      contribution[3] += fac * gqbar * nnlogqbar;	// g qbar
      contribution[4] += fac * gg * nnlogg;		// g g
      contribution[5] += fac * 
	(uq_ig * nnlouqig + dq_ig * nnlodqig 
	 + dbarq_ig * nnlodbarqig); // qbar qbar (nto really)
    } else if (particle == ANTITOPQUARK) {
      contribution[5] += fac * qq * nnloqq;		// qbar qbar
      contribution[1] += fac * qqbar * nnloqqbar;	// q qbar
    }
  }
  
  return( jacobian * hcq * alpha2 / swq /swq * ( lo + as * nlo + as2 * nnlo ) );
}

double HathorSgTopT::ampsq_qq_00(double t) {
  /*
   * t-channel amplitude squared for u b -> d t
   *
   * Taken from hep-ph/0609289, Eq.(A.10)
   */
  
  return( 4. * pi2 * sparton * (sparton - m2) / pow(t - mw2, 2.));
}

double HathorSgTopT::ampsq_qqbar_00(double t) {
  /*
   * t-channel amplitude squared for dbar b -> ubar t
   *
   * Taken from hep-ph/0609289, Eq.(A.11)
   */
  
  const double spt = sparton + t;
  
  return( 4. * pi2 / pow(t - mw2, 2.) * spt * (spt - m2));
}

void  HathorSgTopT::approxNNLO(double t, double s4, double result[8]) {
  /*
   * Evaluate the coefficents of the different plus-distributions 
   * for a specific phase space point determined through s,t
   *
   * The NNLO terms are based on various soft gluon approximations,
   * for the references see below.
   *
   * All scale dependent terms are disabled here, as the full 
   * scale dependency is available with the 'NNLO' option.
   *
   * result array -> 0,1,2 = NLO, 3,4,5,6,7 = NNLO
   *
   */
  // kinematics
  const double s = sparton;
  const double m32 = 0.;
  double u = s4 - s - t + m32 + m2;
  
  for (int i = 0; i < 8; i++)
    result[i] = 0.;
  
  if (llflags & NLL_CALCULATION) {
    // Formulae taken from hep-ph/0609287
    double c3 = 0., c2 = 0., c2mu = 0., T2 = 0., c1mu = 0., G1S11 = 0.;
    // 1-loop soft anomalous dimension (in axial gauge)
    G1S11 = CF * (log(-t/s) + log((m2-t)/(m*sqrts)) + 1.); // Eq.(3.6)
    // everything else (in axial gauge)
    c3 = 3. * CF;
    T2 = 2.*G1S11 - 15./4.*CF - 2.*CF*log((t-m2)*(u-m2)/m4) 
      - 3.*CF*log(m2/s); // Eq.(3.3)
    //c2mu = -2.*CF*Lmf;
    //c1mu = (CF*log((t-m2)*(u-m2)/m4) - 1.5*CF)*Lmf;
    c2 = T2 + c2mu;
    
    if (scheme & NLOAPPROX) {
      // NO PLUS-DISTRIBUTION
      if (llflags & NLO_VIRT)
	result[0] += invpi * c1mu;
      // Coefficient of [1/s4]_+
      if (llflags & NLO_NLL)
	result[1] += invpi * c2; // Eq.(3.1)
      // Coefficient of [log(s4/M^2)/s4]_+
      if (llflags & NLO_LL)
	result[2] += invpi * c3; // Eq.(3.1)
    }
    
    if (scheme & NNLOAPPROX) {
      // Coefficient of [1/s4]_+
      // NNLO-NNNLL: not all terms are known yet
      if (llflags & NNLO_NNNLL)
	result[4] += invpi2 * ( - zeta2*c3*c2 + zeta3*c3*c3); // Eq.(3.11)
      //    + c2mu*c1mu + beta0/4.*c2mu*log(mur2/mt2) 
      //    + CF*beta0/4.*pow(log(muf2/mt2),2.) 
      // Coefficient of [log(s4/M^2)/s4]_+
      if (llflags & NNLO_NNLL)
	result[5] += invpi2 * ( - zeta2*c3*c3 ); // Eq.(3.11)
      //     + c2mu*c1mu + beta0/4.*c2mu*log(mur2/m2)
      //     + CF*beta0/4.*pow(log(muf2/m2),2.)
      // Coefficient of [log^2(s4/M^2)/s4]_+
      if (llflags & NNLO_NLL)
	result[6] += invpi2 * (1.5*c3*c2 - beta0/4.*c3 + CF*beta0/8.); 
      // Eq.(3.11)
      // Coefficient of [log^3(s4/M^2)/s4]_+
      if (llflags & NNLO_LL)
	result[7] += invpi2 * 0.5*c3*c3; // Eq.(3.11)
    }
  } else {
    /*
     * The results below are based on [arXiv:1103.2792]
     */
    double c3 = 0., c2 = 0., c2mu = 0., T2 = 0., c1 = 0., c1mu = 0., 
      G1S11 = 0., G1S12 = 0., G1S21 = 0.;
    // 1-loop soft anomalous dimension (in Feynman gauge)
    G1S11 = CF * (log(-t/s) + log((m2-t)/(m*sqrts)) - 0.5); //Eq.(2.9)
    G1S21 = log(u*(u-m2)/(s*(s-m2))); //Eq.(2.11)
    G1S12 = CF * G1S21 / (2. * Nc); //Eq.(2.11)
    // everything else (in Feynman gauge)
    c3 = 3. * CF;
    T2 = 2.*G1S11 - 0.75*CF - 2.*CF*log((t-m2)*(u-m2)/m4) - 3.*CF*log(m2/s);
    //c2mu = -2.*CF*Lmf; // = 0 for mu = mt, Eq.(2.8)
    //c1mu = (CF*log((t-m2)*(u-m2)/m4) - 1.5*CF)*Lmf; 
    // = 0 for mu = mt, Eq.(2.10)
    c2 = T2 + c2mu; 
    
    if (scheme & NLOAPPROX) {
      // "constant term":
      if (llflags & NLO_VIRT)
	result[0] += invpi * c1mu;
      // Coefficient of [1/s4]_+
      if (llflags & NLO_NLL)
	result[1] += invpi * c2; //Eq.(2.6)
      // Coefficient of [log(s4/M^2)/s4]_+
      if (llflags & NLO_LL)
	result[2] += invpi * c3; //Eq.(2.6)
    }
    
    if (scheme & NNLOAPPROX) {
      // Coefficient of [1/s4]_+
      // NNLO-NNNLL: not all terms are known yet
      if (llflags & NNLO_NNNLL) {
	const double G2S11 = 0.5*K*G1S11 + CF*CA*(1. - zeta3)/4.; //Eq.(2.12)
	
	const double B2 = CF*CF*(-3./32. + 0.75*zeta2 - 1.5*zeta3)
	  + CF*CA*(-1539./864. - 11./12.*zeta2 + 0.75*zeta3)
	  + nf*CF*(135./432. + zeta2/6.); //Eq.(2.5)
	
	result[4] += invpi2 *
	  (c2*c1 - zeta2*c3*c2 + zeta3*c3*c3 + 0.25*beta0*c2*log(m2/s)
	   - 0.5*beta0*CF*pow(log((m2-t)/m2),2.) 
	   - 0.5*beta0*CF*pow(log((m2-u)/m2),2.)
	   - CF*K*log((m2-u)*(m2-t)/m4) + B2 + 3.*D2 
	   + 0.25*CF*beta0*pow(log(m2/s),2.) - CF*K*log(m2/s)
	   + 0.375*beta0*CF*pow(log(m2/s),2.) 
	   - CF*(0.5*K - 0.1875*beta0)*log(m2/s) + 2.*G2S11
	   + (4.*G1S12*G1S21 + 4.*G1S11*G1S11) * log(m2/s)); //Eq.(2.13)
      }
      // Coefficient of [log(s4/M^2)/s4]_+ Eq.(2.13)
      if (llflags & NNLO_NNLL)
	result[5] += invpi2 *
	  (c3*c1 + c2*c2 - zeta2*c3*c3 - 0.5*beta0*T2 //+ 0.25*beta0*c3*Lmrm
	   + 1.5*CF*K - 3./16.*CF*beta0 + 4.*G1S12*G1S21); 
      // Coefficient of [log^2(s4/M^2)/s4]_+ Eq.(2.13)
      if (llflags & NNLO_NLL)
	result[6] += invpi2 * (1.5*c3*c2 - beta0/4.*c3 + CF*beta0/8.);
      // Coefficent of [log^3(s4/M^2)/s4]_+ Eq.(2.13)
      if (llflags & NNLO_LL)
	result[7] += invpi2 * 0.5*c3*c3;
    }
  }
  
}

void HathorSgTopT::evaluateApproxNNLO() {
  // calculate kinematics
  const double s = sparton;
  double tmin  = m2 - s;
  double t     = (-tmin) * rt + tmin;
  double s4max = s + t - m2;
  double s4    = s4max * rs4;
  // jacobian from transformation min..max -> 0..1
  double jac = (-tmin) * s4max;
  
  double f_s4[8],f_0[8];
  
  approxNNLO(t, s4, f_s4);
  approxNNLO(t, 0., f_0);

  // Evaluate Born term
  double FBqq = 1. / (16. * pi * s*s)*ampsq_qq_00(t);
  double FBqqbar = 1. / (16. * pi * s*s)*ampsq_qqbar_00(t);;
  
  // Evaluate plus-distributions
  double logs1 = 0., logs2 = 0.;
  //// NLO
  // no plus-distribution
  logs1 += f_0[0] / s4max;
  // [1/s4]_+
  logs1 += (f_s4[1] - f_0[1]) *           1. / s4     
    +         log(s4max/m2)         * f_0[1] / s4max;
  // [log(s4/M^2)/s4]_+
  logs1 += (f_s4[2] - f_0[2]) *     log(s4/m2)/s4     
    + 0.5 *   pow(log(s4max/m2),2.) * f_0[2] / s4max;
  //// NNLO
  // no plus-distribution
  logs2 += f_0[3] / s4max;
  // [1/s4]_+
  logs2 += (f_s4[4] - f_0[4]) *           1. / s4     
    +         log(s4max/m2)         * f_0[4] / s4max;
  // [log(s4/M^2)/s4]_+
  logs2 += (f_s4[5] - f_0[5]) *     log(s4/m2)/s4     
    + 0.5 *   pow(log(s4max/m2),2.) * f_0[5] / s4max;
  // [log^2(s4/M^2)/s4]_+
  logs2 += (f_s4[6] - f_0[6]) * pow(log(s4/m2),2.)/s4 
    + (1./3.)*pow(log(s4max/m2),3.) * f_0[6] / s4max;
  // [log^3(s4/M^2)/s4]_+
  logs2 += (f_s4[7] - f_0[7]) * pow(log(s4/m2),3.)/s4 
    + 0.25 *  pow(log(s4max/m2),4.) * f_0[7] / s4max;
  
  approxqq1    = FBqq    * jac * logs1;
  approxqqbar1 = FBqqbar * jac * logs1;
  approxqq2    = FBqq    * jac * logs2;
  approxqqbar2 = FBqqbar * jac * logs2;
}

double HathorSgTopT::xs_qq_11() {
  return( - 2. * xPqq0.convolve(rho, rt, pxs_qq_00) );
}

double HathorSgTopT::xs_gu_11() {
  return( - 1. * xPqg0.convolve(rho, rt, pxs_qq_00) );
}

double HathorSgTopT::xs_qqbar_11() {
  return( - 2. * xPqq0.convolve(rho, rt, pxs_qqb_00) );
}

double HathorSgTopT::xs_gqbar_11() {
  return( - 1. * xPqg0.convolve(rho, rt, pxs_qqb_00) );
}

double HathorSgTopT::xs_gb_11() {
  return( - 0.5 * xPqg0.convolve(rho, rt, pxs_qq_00)
	  - 0.5 * xPqg0.convolve(rho, rt, pxs_qqb_00) );
}

double HathorSgTopT::xs_qq_21() {
  return( - 2. * xPqq1v.convolve(rho, rt, pxs_qq_00)
	  - 2. * xPqq0.convolve(rho, rt, pxs_qq_10)
	  + 1. * b0 * pxs_qq_10.eval(sparton) );
}

double HathorSgTopT::xs_qq_22() {
  return( + 2. * xPqqqq.convolve(rho, rt, rs4, pxs_qq_00)
	  - 1. * b0 * xPqq0.convolve(rho, rt, pxs_qq_00) );
}

double HathorSgTopT::xs_qqbar_21() {
  return( - 2. * xPqq1v.convolve(rho, rt, pxs_qqb_00)
	  - 2. * xPqq0.convolve(rho, rt, pxs_qqb_10)
	  + 1. * b0 * pxs_qqb_10.eval(sparton) );
}

double HathorSgTopT::xs_qqbar_22() {
  return( + 2. * xPqqqq.convolve(rho, rt, rs4, pxs_qqb_00)
	  - 1. * b0 * xPqq0.convolve(rho, rt, pxs_qqb_00) );
}

double HathorSgTopT::xs_gu_21() {
  return( - 1. * xPqg1.convolve(rho, rt, pxs_qq_00)
	  - 1. * xPqg0.convolve(rho, rt, pxs_qq_10)
	  - 1. * xPqq0.convolve(rho, rt, pxs_gu_10)
	  - 1. * xPgg0.convolve(rho, rt, pxs_gu_10)
	  + 1. * b0 * pxs_gu_10.eval(sparton) );
}

double HathorSgTopT::xs_gu_22() {
  return( + 1.5 * xPqqqg.convolve(rho, rt, rs4, pxs_qq_00)
	  + 0.5 * xPqggg.convolve(rho, rt, rs4, pxs_qq_00)
	  - 0.5 * b0 * xPqg0.convolve(rho, rt, pxs_qq_00));
}

double HathorSgTopT::xs_gb_21() {
  return( - 0.5 * xPqg1.convolve(rho, rt, pxs_qq_00)
	  - 0.5 * xPqg0.convolve(rho, rt, pxs_qq_10)
	  - 0.5 * xPqg1.convolve(rho, rt, pxs_qqb_00)
	  - 0.5 * xPqg0.convolve(rho, rt, pxs_qqb_10)
	  - 0.5 * xPqq0.convolve(rho, rt, pxs_gb_10)
	  - 0.5 * xPgg0.convolve(rho, rt, pxs_gb_10)
	  + 0.5 * b0 * pxs_gb_10.eval(sparton) );
}

double HathorSgTopT::xs_gb_22() {
  return( + 0.75 * xPqqqg.convolve(rho, rt, rs4, pxs_qq_00)
	  + 0.25 * xPqggg.convolve(rho, rt, rs4, pxs_qq_00)
	  - 0.25 * b0 * xPqg0.convolve(rho, rt, pxs_qq_00)
	  + 0.75 * xPqqqg.convolve(rho, rt, rs4, pxs_qqb_00)
	  + 0.25 * xPqggg.convolve(rho, rt, rs4, pxs_qqb_00)
	  - 0.25 * b0 * xPqg0.convolve(rho, rt, pxs_qqb_00) );
}

double HathorSgTopT::xs_gqbar_21() {
  return( - 1. * xPqg1.convolve(rho, rt, pxs_qqb_00)
	  - 1. * xPqg0.convolve(rho, rt, pxs_qqb_10)
	  - 1. * xPqq0.convolve(rho, rt, pxs_gqb_10)
	  - 1. * xPgg0.convolve(rho, rt, pxs_gqb_10)
	  + 1. * b0 * pxs_gqb_10.eval(sparton) );
}

double HathorSgTopT::xs_gqbar_22() {
  return( + 1.5 * xPqqqg.convolve(rho, rt, rs4, pxs_qqb_00)
	  + 0.5 * xPqggg.convolve(rho, rt, rs4, pxs_qqb_00)
	  - 0.5 * b0 * xPqg0.convolve(rho, rt, pxs_qqb_00) );
}

double HathorSgTopT::xs_gg_21() {
  return( - 2. * xPqg0.convolve(rho, rt, pxs_gu_10)
	  - 2. * xPqg0.convolve(rho, rt, pxs_gb_10)
	  - 2. * xPqg0.convolve(rho, rt, pxs_gqb_10) );
}

double HathorSgTopT::xs_gg_22() {
  return( + 2. * xPqgqg.convolve(rho, rt, rs4, pxs_qq_00)
	  + 2. * xPqgqg.convolve(rho, rt, rs4, pxs_qqb_00) );
}

double HathorSgTopT::xs_uqig_21() {
  return( - 1. * xPqq1s.convolve(rho, rt, pxs_qq_00)
	  - 1. * xPgq0.convolve(rho, rt, pxs_gu_10) );
}

double HathorSgTopT::xs_uqig_22() {
  return( + 0.5 * xPqggq.convolve(rho, rt, rs4, pxs_qq_00) );
}

double HathorSgTopT::xs_dqig_21() {
  return( - 0.5 * xPqq1s.convolve(rho, rt, pxs_qq_00)
	  - 0.5 * xPqq1s.convolve(rho, rt, pxs_qqb_00)
	  - 0.5 * xPgq0.convolve(rho, rt, pxs_gb_10) );
}

double HathorSgTopT::xs_dqig_22() {
  return( + 0.25 * xPqggq.convolve(rho, rt, rs4, pxs_qq_00)
	  + 0.25 * xPqggq.convolve(rho, rt, rs4, pxs_qqb_00) );
}

double HathorSgTopT::xs_dbarqig_21() {
  return( - 1. * xPqq1s.convolve(rho, rt, pxs_qqb_00)
	  - 1. * xPgq0.convolve(rho, rt, pxs_gqb_10) );
}

double HathorSgTopT::xs_dbarqig_22() {
  return( + 0.5 * xPqggq.convolve(rho, rt, rs4, pxs_qqb_00) );
}

double HathorSgTopT::xs_qqbv_21() {
  return( - 1. * xPqqb1v.convolve(rho, rt, pxs_qq_00) );
}

double HathorSgTopT::xs_qqbv_22() {
  return( 0. );
}

double HathorSgTopT::xs_qqv_21() {
  return( - 1. * xPqqb1v.convolve(rho, rt, pxs_qqb_00) );
}

double HathorSgTopT::xs_qqv_22() {
  return( 0. );
}



double HathorSgTopT_LO_qq::eval(double s) {
  /*
   * partonic cross section for u b -> d t
   */
  
  return( pi / (4.) * pow(s - m2, 2.) / (mw2 * s * (s - m2 + mw2)) );
}

double HathorSgTopT_LO_qqbar::eval(double s) {
  /*
   * partonic cross section for dbar b -> ubar t
   */  
  const double s2 = s*s;
  const double a = (s - m2) / 2.;
  
  return( pi * (s - m2) / ( 4. * s2)
	  * (a * (2.*a*mw2 + 2.*mw4 + 2.*mw2*s + s2 - m2*(mw2+s))
	     + mw2 * (2.*a + mw2) * (m2 - 2.*(mw2 + s)) * 0.5 
	     * log(1. + 2.*a/mw2))
	  / (a * mw2 * (2.*a + mw2))
	  );
}

