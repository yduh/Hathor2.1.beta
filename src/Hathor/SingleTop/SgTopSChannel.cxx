#include "SgTopSChannel.h"
#include "SingleTopGrids.h"


using namespace std;


HathorSgTopS::HathorSgTopS(Pdf &pdf_) 
  : SgTop(pdf_),
    pxs_qqbar_10(s_sqrts_udbar,s_xs_udbar, massval),
    pxs_gq_10(s_sqrts_gu,s_xs_gu,massval),
    pxs_gqbar_10(s_sqrts_gdbar,s_xs_gdbar, massval),
    xPqq0(SplittingFunction::PQQ0),
    xPgq0(SplittingFunction::PGQ0),
    xPgg0(SplittingFunction::PGG0),
    xPqg0(SplittingFunction::PQG0),
    xPqg1(SplittingFunction::PQG1),
    xPqq1v(SplittingFunction::PQQ1V),
    xPqq1s(SplittingFunction::PQQ1S),
    xPqqb1v(SplittingFunction::PQQB1V),
    xPqqqg(SplittingFunction::PQQ0,SplittingFunction::PQG0),
  xPggqg(SplittingFunction::PGG0,SplittingFunction::PQG0),
  xPqqqq(SplittingFunction::PQQ0,SplittingFunction::PQQ0),
  xPgqqg(SplittingFunction::PGQ0,SplittingFunction::PQG0),
  xPqgqg(SplittingFunction::PQG0,SplittingFunction::PQG0),
  xPqggg(SplittingFunction::PQG0,SplittingFunction::PGG0),
  xPqggq(SplittingFunction::PQG0,SplittingFunction::PGQ0){

  double scalefac = 1. / fourpi * (MCFMswq*MCFMswq) / (MCFMalpha*MCFMalpha);

  pxs_qqbar_10.setScalefactor(scalefac);
  pxs_gq_10.setScalefactor(scalefac);
  pxs_gqbar_10.setScalefactor(scalefac);
  
}

void HathorSgTopS::update() {
  /*
   * The user has changed the setting for the top mass.
   * We have to update the mass parameter and the production threshold
   * in interpolation routines. 
   */
  SgTop::update();
  
  pxs_qqbar_00.setMass(m, m);
  pxs_qqbar_10.setMass(m, m);
  pxs_gq_10.setMass(m, m);
  pxs_gqbar_10.setMass(m, m);
  
  fqqbar00  = 0.;
  nloqqbar  = nlogq  = nlogqbar = 0.;
  nnloqqbar = nnlogq = nnloqqv  = nnlogg = nnloqqbar2 = nnloqq = 0.;
}

void HathorSgTopS::setPartonicEnergy(const double x[]) {
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
  rho = m2/sparton;
  
  rt  = (dimension > 2) ? x[2] : 0.5;
  rs4 = (dimension > 3) ? x[3] : 0.5;
}

void HathorSgTopS::evaluatePDFs(const double h1[], const double h2[], 
                                const double h2left[], 
				const double h2right[]) {
  /*
   *  Evaluate parton fluxes.
   */
  
  const double ckm_u = Vud*Vud + Vus*Vus + Vub*Vub;
  const double ckm_c = Vcd*Vcd + Vcs*Vcs + Vcb*Vcb;
  const double ckm_t = Vtd*Vtd + Vts*Vts + Vtb*Vtb;
  const double ckm_d = Vud*Vud + Vcd*Vcd;
  const double ckm_s = Vus*Vus + Vcs*Vcs;
  const double ckm_b = Vub*Vub + Vcb*Vcb;
  
  double t_qqbar,    t_gq,    t_gqbar,    t_gg,    t_qqbar2,    t_qq,    t_qqv;
  double tbar_qqbar, tbar_gq, tbar_gqbar, tbar_gg, tbar_qqbar2, 
    tbar_qq, tbar_qqv;
  
  // required for leading order
  t_qqbar    = 
    + (h1[UP] * h2[ADOWN]    * Vud*Vud + h1[CHARM] * h2[ADOWN]    * Vcd*Vcd + 
       h1[UP] * h2[ASTRANGE] * Vus*Vus + h1[CHARM] * h2[ASTRANGE] * Vcs*Vcs +
       h1[UP] * h2[ABOTTOM]  * Vub*Vub + h1[CHARM] * h2[ABOTTOM]  * Vcb*Vcb) 
    * ckm_t 
    + (h2[UP] * h1[ADOWN]    * Vud*Vud + h2[CHARM] * h1[ADOWN]    * Vcd*Vcd + 
       h2[UP] * h1[ASTRANGE] * Vus*Vus + h2[CHARM] * h1[ASTRANGE] * Vcs*Vcs +
       h2[UP] * h1[ABOTTOM]  * Vub*Vub + h2[CHARM] * h1[ABOTTOM]  * Vcb*Vcb) 
    * ckm_t;

  tbar_qqbar = 
    + (h1[DOWN]    * h2[AUP] * Vud*Vud + h1[DOWN]    * h2[ACHARM] * Vcd*Vcd + 
       h1[STRANGE] * h2[AUP] * Vus*Vus + h1[STRANGE] * h2[ACHARM] * Vcs*Vcs +
       h1[BOTTOM]  * h2[AUP] * Vub*Vub + h1[BOTTOM]  * h2[ACHARM] * Vcb*Vcb) 
    * ckm_t 
    + (h2[DOWN]    * h1[AUP] * Vud*Vud + h2[DOWN]    * h1[ACHARM] * Vcd*Vcd + 
       h2[STRANGE] * h1[AUP] * Vus*Vus + h2[STRANGE] * h1[ACHARM] * Vcs*Vcs +
       h2[BOTTOM]  * h1[AUP] * Vub*Vub + h2[BOTTOM]  * h1[ACHARM] * Vcb*Vcb) 
    * ckm_t;
  
  // required for next-to-leading order
  t_gq    = 
    (h1[GLUON] * h2[UP] * ckm_u + h1[GLUON] * h2[CHARM] * ckm_c) * ckm_t +
    (h2[GLUON] * h1[UP] * ckm_u + h2[GLUON] * h1[CHARM] * ckm_c) * ckm_t;
  
  tbar_gq = 
    (h1[GLUON] * h2[DOWN]   * ckm_d + h1[GLUON] * h2[STRANGE] * ckm_s +
     h1[GLUON] * h2[BOTTOM] * ckm_b) * ckm_t +
    (h2[GLUON] * h1[DOWN]   * ckm_d + h2[GLUON] * h1[STRANGE] * ckm_s +
     h2[GLUON] * h1[BOTTOM] * ckm_b) * ckm_t;
 
  t_gqbar    = 
    (h1[GLUON] * h2[ADOWN]   * ckm_d + h1[GLUON] * h2[ASTRANGE] * ckm_s +
     h1[GLUON] * h2[ABOTTOM] * ckm_b) * ckm_t + 
    (h2[GLUON] * h1[ADOWN]   * ckm_d + h2[GLUON] * h1[ASTRANGE] * ckm_s +
     h2[GLUON] * h1[ABOTTOM] * ckm_b) * ckm_t;

  tbar_gqbar = 
    (h1[GLUON] * h2[AUP] * ckm_u + h1[GLUON] * h2[ACHARM] * ckm_c) * ckm_t + 
    (h2[GLUON] * h1[AUP] * ckm_u + h2[GLUON] * h1[ACHARM] * ckm_c) * ckm_t;
  
  // required for next-to-next-to-leading order
  t_gg    = (h1[GLUON] * h2[GLUON]) * ckm_t * (ckm_u + ckm_c);
  tbar_gg = t_gg;
  
  t_qqbar2   =  
    (h1[UP] * h2[ADOWN]    * (ckm_u + ckm_d) 
     + h1[CHARM] * h2[ADOWN]    * (ckm_c + ckm_d) +
     h1[UP] * h2[ASTRANGE] * (ckm_u + ckm_s) 
     + h1[CHARM] * h2[ASTRANGE] * (ckm_c + ckm_s) +
     h1[UP] * h2[ABOTTOM]  * (ckm_u + ckm_b) 
     + h1[CHARM] * h2[ABOTTOM]  * (ckm_c + ckm_b)) * ckm_t 
    + (h2[UP] * h1[ADOWN]    * (ckm_u + ckm_d) 
       + h2[CHARM] * h1[ADOWN]    * (ckm_c + ckm_d) + 
       h2[UP] * h1[ASTRANGE] * (ckm_u + ckm_s) 
       + h2[CHARM] * h1[ASTRANGE] * (ckm_c + ckm_s) +
       h2[UP] * h1[ABOTTOM]  * (ckm_u + ckm_b) 
       + h2[CHARM] * h1[ABOTTOM]  * (ckm_c + ckm_b)) * ckm_t;

  tbar_qqbar2 = (h1[DOWN]    * h2[AUP] * (ckm_u + ckm_d) 
		 + h1[DOWN]    * h2[ACHARM] * (ckm_c + ckm_d) + 
		 h1[STRANGE] * h2[AUP] * (ckm_u + ckm_s) 
		 + h1[STRANGE] * h2[ACHARM] * (ckm_c + ckm_s) +
		 h1[BOTTOM]  * h2[AUP] * (ckm_u + ckm_b) 
		 + h1[BOTTOM]  * h2[ACHARM] * (ckm_c + ckm_b)) * ckm_t +
    (h2[DOWN]    * h1[AUP] * (ckm_u + ckm_d) 
     + h2[DOWN]    * h1[ACHARM] * (ckm_c + ckm_d) + 
     h2[STRANGE] * h1[AUP] * (ckm_u + ckm_s) 
     + h2[STRANGE] * h1[ACHARM] * (ckm_c + ckm_s) +
     h2[BOTTOM]  * h1[AUP] * (ckm_u + ckm_b) 
     + h2[BOTTOM]  * h1[ACHARM] * (ckm_c + ckm_b)) * ckm_t;
  
  t_qq = (h1[UP] * ckm_u + h1[CHARM] * ckm_c) * ckm_t
    * (h2[UP] + h2[CHARM] + h2[AUP] + h2[ACHARM] 
       + h2[DOWN] + h2[STRANGE] + h2[BOTTOM])
    + (h1[ADOWN] * ckm_d + h1[ASTRANGE] * ckm_s + h1[ABOTTOM] * ckm_b) * ckm_t
    * (h2[AUP] + h2[ACHARM] + h2[DOWN] + h2[STRANGE] + h2[BOTTOM] + h2[ADOWN] 
       + h2[ASTRANGE] + h2[ABOTTOM])
    + (h2[UP] * ckm_u + h2[CHARM] * ckm_c) * ckm_t
    * (h1[UP] + h1[CHARM] + h1[AUP] + h1[ACHARM] + h1[DOWN] + h1[STRANGE] 
       + h1[BOTTOM])
    + (h2[ADOWN] * ckm_d + h2[ASTRANGE] * ckm_s + h2[ABOTTOM] * ckm_b) * ckm_t
    * (h1[AUP] + h1[ACHARM] + h1[DOWN] + h1[STRANGE] + h1[BOTTOM] + h1[ADOWN] 
       + h1[ASTRANGE] + h1[ABOTTOM]);

  tbar_qq = (h1[DOWN] * ckm_d + h1[STRANGE] * ckm_s + h1[BOTTOM] * ckm_b) 
    * ckm_t * (h2[UP] + h2[CHARM] + h2[ADOWN] + h2[ASTRANGE] + h2[ABOTTOM] 
	       + h2[DOWN] + h2[STRANGE] + h2[BOTTOM])
    + (h1[AUP] * ckm_u + h1[ACHARM] * ckm_c) * ckm_t
    * (h2[AUP] + h2[ACHARM] + h2[UP] + h2[CHARM] + h2[ADOWN] 
       + h2[ASTRANGE] + h2[ABOTTOM])
    + (h2[DOWN] * ckm_d + h2[STRANGE] * ckm_s + h2[BOTTOM] * ckm_b) * ckm_t
    * (h1[UP] + h1[CHARM] + h1[ADOWN] + h1[ASTRANGE] + h1[ABOTTOM] 
       + h1[DOWN] + h1[STRANGE] + h1[BOTTOM])
    + (h2[AUP] * ckm_u + h2[ACHARM] * ckm_c) * ckm_t
    * (h1[AUP] + h1[ACHARM] + h1[UP] + h1[CHARM] + h1[ADOWN] + h1[ASTRANGE] 
       + h1[ABOTTOM]);
  
  t_qqv = (h1[UP] * h2[DOWN] * Vud*Vud      + h1[UP] * h2[STRANGE] * Vus*Vus 
	   + h1[UP] * h2[BOTTOM] * Vub*Vub)
    + (h1[CHARM] * h2[DOWN] * Vcd*Vcd   + h1[CHARM] * h2[STRANGE] * Vcs*Vcs   
       + h1[CHARM] * h2[BOTTOM] * Vcb*Vcb)
    + (h1[AUP] * h2[ADOWN] * Vud*Vud    + h1[AUP] * h2[ASTRANGE] * Vus*Vus    
       + h1[AUP] * h2[ABOTTOM] * Vub*Vub)
    + (h1[ACHARM] * h2[ADOWN] * Vcd*Vcd + h1[ACHARM] * h2[ASTRANGE] * Vcs*Vcs 
       + h1[ACHARM] * h2[ABOTTOM] * Vcb*Vcb)
    + (h2[UP] * h1[DOWN] * Vud*Vud      + h2[UP] * h1[STRANGE] * Vus*Vus  
       + h2[UP] * h1[BOTTOM] * Vub*Vub)
    + (h2[CHARM] * h1[DOWN] * Vcd*Vcd   + h2[CHARM] * h1[STRANGE] * Vcs*Vcs 
       + h2[CHARM] * h1[BOTTOM] * Vcb*Vcb)
    + (h2[AUP] * h1[ADOWN] * Vud*Vud    + h2[AUP] * h1[ASTRANGE] * Vus*Vus 
       + h2[AUP] * h1[ABOTTOM] * Vub*Vub)
    + (h2[ACHARM] * h1[ADOWN] * Vcd*Vcd + h2[ACHARM] * h1[ASTRANGE] * Vcs*Vcs 
       + h2[ACHARM] * h1[ABOTTOM] * Vcb*Vcb);

  t_qqv *= ckm_t;

  tbar_qqv = t_qqv;
  
  // set fluxes according to chosen particle
  if (particle == TOPQUARK) {
    qqbar  = t_qqbar;
    gqbar  = t_gqbar;  gq = t_gq;
    qqbar2 = t_qqbar2; gg = t_gg; qq = t_qq;
    qqv    = t_qqv;
    return;
  } 
  
  if (particle == ANTITOPQUARK) {
    qqbar  = tbar_qqbar;
    gqbar  = tbar_gqbar;  gq = tbar_gq;
    qqbar2 = tbar_qqbar2; gg = tbar_gg; qq = tbar_qq;
    qqv    = tbar_qqv;
    return;
  } 

  if (particle == BOTH) {
    qqbar  = t_qqbar  + tbar_qqbar;
    gqbar  = t_gqbar  + tbar_gqbar;  gq = t_gq + tbar_gq;
    qqbar2 = t_qqbar2 + tbar_qqbar2; gg = t_gg + tbar_gg; qq = t_qq + tbar_qq;
    qqv    = t_qqv    + tbar_qqv;
    return;
  }

  cout << "# Hathor error: set particle to TOPQUARK, ANTITOPQUAR or BOTH \n";
  exit(1);

}

void HathorSgTopS::evaluateScalingFunctions() {
  /*
   *  Evaluate scaling functions for further use.
   */

  // Leading order
  fqqbar00 = pxs_qqbar_00.eval(sparton);
  
  if ( scheme & ( NLOAPPROX | NNLOAPPROX ) )
    evaluateApproxNNLO();
  
  // Next-to-leading order
  double fqqbar10 = 0., fgqbar10 = 0., fgq10 = 0.;
  double fqqbar11 = 0., fgq11 = 0.;

  if ( scheme & ( NLO | NNLO | NLOAPPROX | NNLOAPPROX ) ) {
    // NLO scaling functions to calculate mu-dependence
    fqqbar10 = pxs_qqbar_10.eval(sparton);
    fgqbar10 = pxs_gqbar_10.eval(sparton);
    fgq10    = pxs_gq_10.eval(sparton);
    fqqbar11 = xs_qqbar_11();
    fgq11    = xs_gq_11();
  }

  /* 
   * Combine the scale dependent contributions with the scale
   * independent contributions for the two different approximations:
   */
  if ( scheme & NLO ) {
    nloqqbar = fourpi * ( fqqbar10 + Lmf * fqqbar11 );
    nlogqbar = fourpi * ( fgqbar10 + Lmf * fgq11    );
    nlogq    = fourpi * ( fgq10    + Lmf * fgq11    );
  } else if (scheme & NLOAPPROX) {
    // always use full scale-dependency
    nloqqbar = fourpi * (         + Lmf * fqqbar11) + approxqqbar1;
    nlogqbar = fourpi * (         + Lmf * fgq11   );
    nlogq    = fourpi * (         + Lmf * fgq11   );
  }
  
  if ( scheme & FAKE_NNLO ) {
    /* This option is used to estimate what can be expected from the 
     * full NNLO calcu;ation as far as the scaledependence is concerned.
     * We estimate the size of the NNLO corrections for mu=mt using
     * a simple ansatz:
     */
    // - scale NLO result with ratio (NLO/LO)
    approxqqbar2 = fourpi * fqqbar10 * (fqqbar10 / fqqbar00);
    // - simply use LO result
    //approxqqbar2 = fqqbar00;
    // - model with sqrt(mt^2 / s)
    // approxqqbar2 = sqrt(rho / hcq) / (twopi*twopi);
  }
  
  /*
   * Calculate the full scale dependence
   * Whatever approximation is used we always include the complete
   * scale (in)dependence!
   */ 
  if ( scheme & ( NNLO | NNLOAPPROX | FAKE_NNLO ) ) {

    double fqqbar21, fgg21, fgq21, fqq21, fqqv21,fqqbar22, fgg22, 
      fgq22, fqq22, fqqv22;
    
    fqqbar21   = xs_qqbar_21();
    fgg21      = xs_gg_21();
    fgq21      = xs_gq_21();
    fqq21      = xs_qq_21();
    fqqv21     = xs_qqv_21();
    fqqbar22   = xs_qqbar_22();
    fgg22      = xs_gg_22();
    fgq22      = xs_gq_22();
    fqq22      = xs_qq_22();
    fqqv22     = xs_qqv_22();
    

    nnloqqbar  = fourpi2 * Lmf * (fqqbar21 + Lmf * fqqbar22)
      + fourpi2 * Lmr * (fqqbar10 + Lmf * fqqbar11) * b0;
    nnloqqbar2 = fourpi2 * Lmf * (fqq21    + Lmf * fqq22);
    nnlogg     = fourpi2 * Lmf * (fgg21    + Lmf * fgg22);
    nnlogq     = fourpi2 * Lmf * (fgq21    + Lmf * fgq22)
      + fourpi2 * Lmr * (fgq10    + Lmf * fgq11) * b0;
    nnloqq     = fourpi2 * Lmf * (fqq21    + Lmf * fqq22);
    nnloqqv    = fourpi2 * Lmf * (fqqv21   + Lmf * fqqv22);
    // Add scale-independent approx. NNLO terms
    if ( scheme & ( NNLOAPPROX | FAKE_NNLO ) )
      nnloqqbar  += approxqqbar2;
  }
}

double HathorSgTopS::evaluateIntegral(double as, double wgt) {
  /*
   *  Calculate the full integrand, using the given alpha_s.
   */

  double lo = 0, nlo = 0, nnlo = 0;
  const double as2 = as * as;
  
  double fac = jacobian * hcq;

  if (scheme & LO) {
    lo = qqbar * fqqbar00;
    contribution[1] +=  fac * qqbar * fqqbar00;	// q qbar
  }

  fac *= as;
 
  if (scheme & (NLO | NLOAPPROX)) {
    nlo = qqbar * nloqqbar + gqbar * nlogqbar + gq * nlogq;
    contribution[1] += fac * qqbar * nloqqbar; // q qbar
    contribution[2] += fac * gq    * nlogq;    // g q
    contribution[3] += fac * gqbar * nlogqbar; // g qbar
  }
  
  fac *= as;

  if (scheme & (NNLO | NNLOAPPROX | FAKE_NNLO)) {
    nnlo = qqbar  * nnloqqbar  + qqbar2 * nnloqqbar2 
      + gg * nnlogg   + (gq + gqbar) * nnlogq  + qq * nnloqq
      + qqv * nnloqqv;
    contribution[1] += fac * (qqbar * nnloqqbar + qqbar2 * nnloqqbar2);//q qbar
    contribution[2] += fac * gq    * nnlogq; // g q
    contribution[3] += fac * gqbar * nnlogq; // g qbar
    contribution[4] += fac * gg    * nnlogg; // g g
    contribution[5] += fac * qq    * nnloqq; // "qbar qbar"
  }

  return( jacobian * hcq * alpha2 /swq/swq * ( lo + as * nlo + as2 * nnlo ) );
}

double HathorSgTopS::xs_qqbar_00(double t) {
  /*
   * s-channel amplitude squared for u dbar -> bbar t
   *
   * Taken from hep-ph/0609289, Eq.(A.12)
   */
  return(4. * pi2 * t * (t - m2) / pow(sparton - mw2, 2.));
}

void HathorSgTopS::approxNNLO(double t, double s4, double result[8]) {
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
   */
  // kinematics
  const double s = sparton;
  double u = s4 - s - t + m2;
  
  for (int i = 0; i < 8; i++)
    result[i] = 0.;

  if ( llflags & NLL_CALCULATION ) {
    /*
     * NLL results from hep-ph/0609287
     */ 
    double c3 = 0., c2 = 0., c2mu = 0., T2 = 0., c1mu = 0.;
    // 1-loop soft anomalous dimension (in axial gauge)
    double G1S11 = CF * (log((s-m2)/(m*sqrts)) + 1.); // Eq.(3.7)
    // Everything else (in axial gauge)
    c3 = 3. * CF;
    T2 = + 2.*G1S11 - 15./4.*CF - 2.*CF*log((t-m2)*(u-m2)/m4) 
      - 3.*CF*log(m2/s) ; // Eq.(3.2)
    c2mu = 0.; // -2.*CF*Lmf = 0 for mu=m;
    c1mu = 0.; // (CF*log((m2-t)*(m2-u)/m4) - 1.5*CF)*Lmf = 0 for mu = m;
    c2 = T2 + c2mu;

    if (scheme & NLOAPPROX) {
      // NO PLUS-DISTRIBUTION
      if (llflags & NLO_VIRT)
	result[0] += invpi * c1mu;
      // INTEGRAND IN FRONT OF [1/s4]_+
      if (llflags & NLO_NLL)
	result[1] += invpi * c2; // Eq.(3.1)
      // INTEGRAND IN FRONT OF [log(s4/M^2)/s4]_+
      if (llflags & NLO_LL)
	result[2] += invpi * c3; // Eq.(3.1)
    }
    
    if (scheme & NNLOAPPROX) {
      // Coefficient of [1/s4]_+
      // NNLO-NNNLL: not all terms are known yet
      if (llflags & NNLO_NNNLL)
	result[4] += invpi2 * ( - zeta2*c3*c2 + zeta3*c3*c3); // Eq.(3.12)
      // + c2mu*c1mu + beta0/4.*c2mu*log(mur2/m2) = 0 for mu = m
      // + CF*beta0/4.*pow(log(muf2/m2),2.) = 0 for mu = m;
      // Coefficient of [log(s4/M^2)/s4]_+
      if (llflags & NNLO_NNLL)
	result[5] += invpi2 * ( - zeta2*c3*c3); // Eq.(3.12)
      // + c3*c1mu + c2mu*c2mu + 2.*c2mu*T2 + beta0/4.*c3*log(mur2/mt2)
      // Coefficient of [log^2(s4/M^2)/s4]_+
      if (llflags & NNLO_NLL)
	result[6] += invpi2 * (1.5*c3*c2 - beta0/4.*c3 + CF*beta0/8.);
      // Coefficient of [log^3(s4/M^2)/s4]_+
      if (llflags & NNLO_LL)
	result[7] += invpi2 * 0.5*c3*c3; // Eq.(3.12)
    }
  } else {
    /*
     * NNLL from arXiv:1001.5034
     */
    double c3 = 0., c2 = 0., c2mu = 0., T2 = 0., c1 = 0., c1mu = 0., T1 = 0.;
    double G1S11 = 0., G1S12 = 0., G1S21 = 0.;
    // 1-loop soft anomalous dimension (in Feynman gauge)
    G1S11 = CF * ( log( (s-m2)/(m*sqrts) ) - 0.5); // Eq.(2.6)
    G1S21 = log( u*(m2-u) / (t*(m2-t)) ); // Eq.(2.7)
    G1S12 = CF * G1S21 / (2. * Nc); // Eq.(2.7)
    // Everything else (in Feynman gauge)
    c3 = 3. * CF; // Eq.(2.10)
    T2 = -2.*CF*log( (m2-t)*(m2-u)/m4 )  + 0 
      - 3.*CF*log(m2/s) - 0.75*CF + 2.*G1S11; // Eq.(2.12)
    c2mu = 0.;// Eq.(2.11) -2.*CF*Lmf = 0 for mu=mt;
    c1mu = 0.;// Eq.(2.13) (CF*log((m2-t)*(m2-u)/m4) - 1.5*CF)*Lmf =0 for mu=mt;
    c2 = T2 + c2mu;
    T1 = 0.;   // not yet known
    c1 = T1 + c1mu;
    
    if (scheme & NLOAPPROX) {
      // NO PLUS-DISTRIBUTION
      result[0] += invpi * c1;
      // Coefficient of {\cal D}_0(s4) = [1/s4]_+
      if (llflags & NLO_NLL)
	result[1] += invpi * c2; // Eq.(2.9)
      // Coefficient of {\cal D}_0(s4) = [log(s4/M^2)/s4]_+
      if (llflags & NLO_LL)
	result[2] += invpi * c3; // Eq.(2.9)
    }
    
    if (scheme & NNLOAPPROX) {
      // Coefficient of {\cal D}_0(s_4):
      // NNLO-NNNLL: not all terms are known yet
      if (llflags & NNLO_NNNLL) {
	const double G2S11 = 0.5*K*G1S11 + CF*CA*(1. - zeta3)/4.; //Eq.(2.8)
	const double B2 = CF*CF*(-3./32. + 0.75*zeta2 - 1.5*zeta3)
	  + CF*CA*(77./864. - 11./4.*zeta2 - zeta3)
	  + nf*CF*(23./432. + zeta2/2.); //Eq.(2.5)
	// Eq.(2.14):
	result[4] += invpi2 *
	  (c2*c1 - zeta2*c3*c2 + zeta3*c3*c3 + 0.25*beta0*c2*log(m2/s) 
	   - 0.5*beta0*CF*pow(log((m2-t)/m2),2.)
	   - 0.5*beta0*CF*pow(log((m2-u)/m2),2.)
	   - CF*K*log((m2-u)*(m2-t)/m4) + B2 + 3.*D2 
	   + 0.25*CF*beta0*pow(log(m2/s),2.) - CF*K*log(m2/s)
	   + 0.375*beta0*CF*pow(log(m2/s),2.) 
	   - CF*(0.5*K - 0.1875*beta0)*log(m2/s) + 2.*G2S11
	   + (4.*G1S12*G1S21) * log(m2/s));
      }
      // Coefficient of {\cal D}_1(s4) = [log(s4/M^2)/s4]_+
      if (llflags & NNLO_NNLL){
	// Eq.(2.14):
	result[5] += invpi2 *
	  (c3*c1 + c2*c2 - zeta2*c3*c3 - 0.5*beta0*T2 //+ 0.25*beta0*c3*Lmrm
	   + 1.5*CF*K - 3./16.*CF*beta0 + 4.*G1S12*G1S21);
      }
      // Coefficient of {\cal D}_2(s4) = [log^2(s4/M^2)/s4]_+
      if (llflags & NNLO_NLL){
	// Eq.(2.14):
	result[6] += invpi2 * (1.5*c3*c2 - beta0/4.*c3 + CF*beta0/8.);
      }
      // Coefficient of {\cal D}_3(s4) = [log^3(s4/M^2)/s4]_+
      if (llflags & NNLO_LL){
	// Eq.(2.14):
	result[7] += invpi2 * 0.5*c3*c3;
      }
    }
  }
}

void HathorSgTopS::evaluateApproxNNLO() {
  // calculate kinematics
  const double s = sparton;
  double t0   = 0.5*(-s + m2);
  double tpm  = 0.5*sqrt(s*s - 2.*s*(m2) + pow(- m2,2.));
  double tmax = t0 + tpm;
  double tmin = t0 - tpm;
  double t  = (tmax - tmin) * rt + tmin;
  double s4max = s + t - m2;
  double s4 = s4max * rs4;
  // jacobian from transformation min..max -> 0..1
  double jac = (tmax - tmin) * s4max;
  
  double f_s4[8], f_0[8];


  approxNNLO(t, s4, f_s4);
  approxNNLO(t, 0., f_0);

  // evaluate Born term
  double FBqqbar = 1. / (16. * pi * s*s) * xs_qqbar_00(t);

  // evaluate plus-distributions
  double logs1 = 0., logs2 = 0.;
  // NLO accuracy
  // no plus-distribution
  logs1 += f_0[0] / s4max;
  // [1/s4]_+
  logs1 += (f_s4[1] - f_0[1]) * 1. / s4   
    + log(s4max/m2)  * f_0[1] / s4max;
  // [log(s4/M^2)/s4]_+
  logs1 += (f_s4[2] - f_0[2]) * log(s4/m2)/s4 
    + 0.5 *   pow(log(s4max/m2),2.) * f_0[2] / s4max;
  // NNLO accuracy
  // no plus-distribution
  logs2 += f_0[3] / s4max;
  // [1/s4]_+
  logs2 += (f_s4[4] - f_0[4]) * 1. / s4 
    + log(s4max/m2) * f_0[4] / s4max;
  // [log(s4/M^2)/s4]_+
  logs2 += (f_s4[5] - f_0[5]) * log(s4/m2)/s4 
    + 0.5 *   pow(log(s4max/m2),2.) * f_0[5] / s4max;
  // [log^2(s4/M^2)/s4]_+
  logs2 += (f_s4[6] - f_0[6]) * pow(log(s4/m2),2.)/s4 
    + (1./3.)*pow(log(s4max/m2),3.) * f_0[6] / s4max;
  // [log^3(s4/M^2)/s4]_+
  logs2 += (f_s4[7] - f_0[7]) * pow(log(s4/m2),3.)/s4 
    + 0.25 *  pow(log(s4max/m2),4.) * f_0[7] / s4max;


  approxqqbar1 = FBqqbar * jac * logs1;
  approxqqbar2 = FBqqbar * jac * logs2;
}

double HathorSgTopS::xs_qqbar_11() {
  return -2. * xPqq0.convolve(rho, rt, pxs_qqbar_00);
}

double HathorSgTopS::xs_qqbar_21() {
  return( - 2. * xPqq1v.convolve(rho, rt, pxs_qqbar_00)
	  - 2. * xPqq0.convolve(rho, rt, pxs_qqbar_10)
	  + 1. * b0 * pxs_qqbar_10.eval(sparton));
}

double HathorSgTopS::xs_qqbar_22() {
  return + 2. * xPqqqq.convolve(rho, rt, rs4, pxs_qqbar_00)
    - 1. * b0 * xPqq0.convolve(rho, rt, pxs_qqbar_00);
}

double HathorSgTopS::xs_gq_11() {
  return - xPqg0.convolve(rho, rt, pxs_qqbar_00);
}

double HathorSgTopS::xs_gq_21() {
  return - xPqq0.convolve(rho, rt, pxs_gq_10)
    - xPgg0.convolve(rho, rt, pxs_gq_10)
    - xPqg0.convolve(rho, rt, pxs_qqbar_10)
    - xPqg1.convolve(rho, rt, pxs_qqbar_00)
    + b0 * pxs_gq_10.eval(sparton);
}

double HathorSgTopS::xs_gq_22() {
  return 1.5 * xPqqqg.convolve(rho, rt, rs4, pxs_qqbar_00)
    + 0.5 * xPqggg.convolve(rho, rt, rs4, pxs_qqbar_00)
    - 0.5 * b0 * xPqg0.convolve(rho, rt, pxs_qqbar_00);
}

double HathorSgTopS::xs_qq_21() {
  return - 1. * xPgq0.convolve(rho, rt, pxs_gq_10) 
    - 1. * xPqq1s.convolve(rho, rt, pxs_qqbar_00);
}

double HathorSgTopS::xs_qq_22() {
  return + 0.5 * xPqggq.convolve(rho, rt, rs4, pxs_qqbar_00);
}

double HathorSgTopS::xs_gg_21() {
  return - 4. * xPqg0.convolve(rho, rt, pxs_gq_10);
}

double HathorSgTopS::xs_gg_22() {
  return + 2. * xPqgqg.convolve(rho, rt, rs4, pxs_qqbar_00);
}

double HathorSgTopS::xs_qqv_21() {
  return - 1. * xPqqb1v.convolve(rho, rt, pxs_qqbar_00);
}

double HathorSgTopS::xs_qqv_22() {
  return 0.;
}



double HathorSgTopS_LO_qqbar::eval(double s) {
  /*
   * partonic cross section for dbar b -> ubar t
   */
  
  const double s2 = s*s;
  
  return pi / (24.0) * (s - m2) * (2.*s2 - s*m2 - m2*m2)
    / (pow(s - mw2, 2.) * s2);
}

