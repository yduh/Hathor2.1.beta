#include "SgTopWtChannel.h"
#include "SingleTopGrids.h"

using namespace std;



HathorSgTopWt::HathorSgTopWt(Pdf &pdf_) 
  : SgTop(pdf_), 
    pxs_gq_10(wt_sqrts_gb, wt_xs_gb, massval),
    pxs_gg_10(wt_sqrts_gg, wt_xs_gg, massval),
    pxs_qq_10(wt_sqrts_ub, wt_xs_ub, massval),
    xPqq0(SplittingFunction::PQQ0),
    xPqg0(SplittingFunction::PQG0),
    xPgg0(SplittingFunction::PGG0),
    xPgq0(SplittingFunction::PGQ0),
    xPgg1(SplittingFunction::PGG1),
    xPgq1(SplittingFunction::PGQ1),
    xPqg1(SplittingFunction::PQG1),
    xPqq1s(SplittingFunction::PQQ1S),
    xPqq1v(SplittingFunction::PQQ1V),
  xPqqb1v(SplittingFunction::PQQB1V),
  xPgggg(SplittingFunction::PGG0,SplittingFunction::PGG0),
  xPgqgg(SplittingFunction::PGQ0,SplittingFunction::PGG0),
  xPqggg(SplittingFunction::PQG0,SplittingFunction::PGG0),
  xPqggq(SplittingFunction::PQG0,SplittingFunction::PGQ0),
  xPqqgg(SplittingFunction::PQQ0,SplittingFunction::PGG0),
  xPqqgq(SplittingFunction::PQQ0,SplittingFunction::PGQ0),
  xPqqqg(SplittingFunction::PQQ0,SplittingFunction::PQG0),
  xPqqqq(SplittingFunction::PQQ0,SplittingFunction::PQQ0){
  
  pxs_gq_10.setScalefactor(1./fourpi * MCFMswq / MCFMalpha);
  pxs_gg_10.setScalefactor(1./fourpi * MCFMswq / MCFMalpha);
  pxs_qq_10.setScalefactor(1./fourpi * MCFMswq / MCFMalpha);

}

void HathorSgTopWt::update() {
  /*
   *  The user has changed the setting for the top mass.
   * We have to update the mass parameter and the production threshold
   * in interpolation routines. 
   */
  
  SgTop::update();
  
  pxs_gq_00.setMass(m, m+mw);
  pxs_gq_10.setMass(m, m+mw);
  pxs_gg_10.setMass(m, m+mw);
  pxs_qq_10.setMass(m, m+mw);
  
  fgq00  = 0.;
  nlogq  = nlogg  = nlobq  = 0.;
  nnlogq = nnlogg = nnlobq = nnlogqig = nnlogqbig2 = 0.;
}

void HathorSgTopWt::setPartonicEnergy(const double x[]) {
  /*
   *  Calculate kinematic variables.
   */
  
  // exclude energies below kinematic threshold
  double mthresh2 = pow(m + mw, 2.);
  sparton = mthresh2*pow(shad/mthresh2,x[0]);
  double r = sparton/shad;
  x1 = r*pow(r,-x[1]);
  x2 = r/x1;
  jacobian = -r*log(shad/mthresh2)*log(r);
  
  sqrts = sqrt(sparton);
  rho   = mthresh2/sparton;
  
  rt  = (dimension > 2) ? x[2] : 0.5;
  rs4 = (dimension > 3) ? x[3] : 0.5;
}

void HathorSgTopWt::evaluatePDFs(const double h1[], const double h2[], 
				 const double h2left[], 
				 const double h2right[]) {
  /*
   *  Evaluate parton fluxes.
   */
  
  double ckm_t = Vtd*Vtd + Vts*Vts + Vtb*Vtb;
  
  double t_gb,t_gb_h, t_bq, t_gg, t_bqbar, t_gqig, t_gqbarig, t_gqbarig2,
    t_gqbarig2_h;
  double tbar_gb,tbar_gb_h, tbar_bq, tbar_gg, tbar_bqbar;
  double tbar_gqig, tbar_gqbarig, tbar_gqbarig2,tbar_gqbarig2_h;
  
  // required for leading order
  t_gb = h1[GLUON] * 
    (h2[BOTTOM] * Vtb*Vtb + h2[STRANGE] * Vts*Vts + h2[DOWN] * Vtd*Vtd) +
    h2[GLUON] * (h1[BOTTOM] * Vtb*Vtb + h1[STRANGE] * Vts*Vts 
		 + h1[DOWN] * Vtd*Vtd);

  t_gb_h = 0.0;
  

  tbar_gb = h1[GLUON] * (h2[ABOTTOM] * Vtb*Vtb + h2[ASTRANGE] * Vts*Vts 
			 + h2[ADOWN] * Vtd*Vtd) 
    + h2[GLUON] * (h1[ABOTTOM] * Vtb*Vtb + h1[ASTRANGE] * Vts*Vts 
		   + h1[ADOWN] * Vtd*Vtd);

  tbar_gb_h = 0.0;

  // required for next-to-leading order
  t_gg    = h1[GLUON] * h2[GLUON] * ckm_t;
  tbar_gg = h1[GLUON] * h2[GLUON] * ckm_t;

  t_bq    = 
    (h1[BOTTOM] * Vtb*Vtb + h1[STRANGE] * Vts*Vts + h1[DOWN] * Vtd*Vtd) * 
    (h2[UP] + h2[DOWN] + h2[CHARM] + h2[STRANGE] + h2[BOTTOM]) +
    (h2[BOTTOM] * Vtb*Vtb + h2[STRANGE] * Vts*Vts + h2[DOWN] * Vtd*Vtd) *
    (h1[UP] + h1[DOWN] + h1[CHARM] + h1[STRANGE] + h1[BOTTOM]);
  
  tbar_bq = 
    (h1[ABOTTOM] * Vtb*Vtb + h1[ASTRANGE] * Vts*Vts + h1[ADOWN] * Vtd*Vtd) *
    (h2[UP] + h2[DOWN] + h2[CHARM] + h2[STRANGE] + h2[BOTTOM] ) +
    (h2[ABOTTOM] * Vtb*Vtb + h2[ASTRANGE] * Vts*Vts + h2[ADOWN] * Vtd*Vtd) *
    (h1[UP] + h1[DOWN] + h1[CHARM] + h1[STRANGE] + h1[BOTTOM]);

  t_bqbar    = 
    (h1[BOTTOM] * Vtb*Vtb + h1[STRANGE] * Vts*Vts + h1[DOWN] * Vtd*Vtd) *
    (h2[AUP] + h2[ADOWN] + h2[ACHARM] + h2[ASTRANGE]) +
    (h2[BOTTOM] * Vtb*Vtb + h2[STRANGE] * Vts*Vts + h2[DOWN] * Vtd*Vtd) *
    (h1[AUP] + h1[ADOWN] + h1[ACHARM] + h1[ASTRANGE]);

  tbar_bqbar = 
    (h1[ABOTTOM] * Vtb*Vtb + h1[ASTRANGE] * Vts*Vts + h1[ADOWN] * Vtd*Vtd) *
    (h2[AUP] + h2[ADOWN] + h2[ACHARM] + h2[ASTRANGE]) +
    (h2[ABOTTOM] * Vtb*Vtb + h2[ASTRANGE] * Vts*Vts + h2[ADOWN] * Vtd*Vtd) *
    (h1[AUP] + h1[ADOWN] + h1[ACHARM] + h1[ASTRANGE]);
  
  // required for next-to-next-to-leading order
  t_gqig = h1[GLUON] 
    * (h2[UP] + h2[CHARM] + h2[DOWN] + h2[STRANGE] + h2[BOTTOM]) * ckm_t
    + h2[GLUON] * (h1[UP] + h1[CHARM] + h1[DOWN] + h1[STRANGE] + h1[BOTTOM]) 
    * ckm_t;

  tbar_gqig = t_gqig;

  t_gqbarig = h1[GLUON] * (h2[AUP] + h2[ACHARM] + h2[ADOWN] 
			   + h2[ASTRANGE] + h2[ABOTTOM]) * ckm_t
    + h2[GLUON] * (h1[AUP] + h1[ACHARM] + h1[ADOWN] + h1[ASTRANGE] 
		   + h1[ABOTTOM]) * ckm_t;
  
  tbar_gqbarig = t_gqbarig;
  
  t_gqbarig2 = h1[GLUON] * (h2[ADOWN]*Vtd*Vtd + h2[ASTRANGE]*Vts*Vts 
			    + h2[ABOTTOM]*Vtb*Vtb)
    + h2[GLUON] * (h1[ADOWN]*Vtd*Vtd + h1[ASTRANGE]*Vts*Vts 
		   + h1[ABOTTOM]*Vtb*Vtb);

  t_gqbarig2_h = 0.0;
  
  tbar_gqbarig2 = h1[GLUON] * (h2[DOWN]*Vtd*Vtd + h2[STRANGE]*Vts*Vts 
			       + h2[BOTTOM]*Vtb*Vtb)
    + h2[GLUON] * (h1[DOWN]*Vtd*Vtd + h1[STRANGE]*Vts*Vts + h1[BOTTOM]*Vtb*Vtb);

  tbar_gqbarig2_h = 0.0;
  
  if (!enableBB) {
    // if set, include b b and b bbar initial states
    t_bq    -= h1[BOTTOM] * Vtb*Vtb  * h2[BOTTOM] * 2.;
    tbar_bq -= (h1[ABOTTOM] * h2[BOTTOM] + h1[BOTTOM] * h2[ABOTTOM]) * Vtb*Vtb;
    t_bqbar -= (h1[BOTTOM] * h2[ABOTTOM]  + h1[ABOTTOM] * h2[BOTTOM]) * Vtb*Vtb;
    tbar_bqbar -= h1[ABOTTOM] * h2[ABOTTOM] * Vtb*Vtb * 2.;
    // the channels gb and gbbar change in NNLO, if bb and b bbar don't exist in
    // NLO
    t_gb_h  = (h1[GLUON] * h2[BOTTOM] + h1[BOTTOM] * h2[GLUON]) * Vtb*Vtb;
    tbar_gb_h  = (h1[GLUON] * h2[ABOTTOM] + h1[ABOTTOM] * h2[GLUON]) * Vtb*Vtb;
    t_gqbarig2_h =  tbar_gb_h;
    tbar_gqbarig2_h = t_gb_h;
  }

  // set fluxes according to chosen particle
  if (particle == TOPQUARK) {
    gb = t_gb; gb_h = t_gb_h;
    gg = t_gg; bq = t_bq; bqbar = t_bqbar;
    gqig = t_gqig; gqbarig = t_gqbarig; 
    gqbarig2 = t_gqbarig2; 
    gqbarig2_h = t_gqbarig2_h;
    return;
  }

  if (particle == ANTITOPQUARK) {
    gb = tbar_gb; gb_h = tbar_gb_h;
    gg = tbar_gg; bq = tbar_bq; bqbar = tbar_bqbar;
    gqig = tbar_gqig; gqbarig = tbar_gqbarig; 
    gqbarig2 = tbar_gqbarig2; gqbarig2_h = tbar_gqbarig2_h;
    return;
  }

  if (particle == BOTH) {
    gb = t_gb + tbar_gb; gb_h = t_gb_h + tbar_gb_h;
    gg = t_gg + tbar_gg; bq = t_bq + tbar_bq; bqbar = t_bqbar + tbar_bqbar;
    gqig = t_gqig + tbar_gqig; 
    gqbarig = t_gqbarig + tbar_gqbarig;
    gqbarig2 = t_gqbarig2 + tbar_gqbarig2; 
    gqbarig2_h = t_gqbarig2_h + tbar_gqbarig2_h;
    return;
  }
  cout << "# Hathor error: set particle to TOPQUARK, ANTITOPQUAR or BOTH \n";
  exit(1);
}

void HathorSgTopWt::evaluateScalingFunctions() {
  /*
   *  Evaluate scaling functions for further use.
   */
  
  if (scheme & (NLOAPPROX | NNLOAPPROX))
    evaluateApproxNNLO();
  
  // leading order
  fgq00 = pxs_gq_00.eval(sparton);
  
  // next-to-leading order
  double fqq10 = 0., fgg10 = 0., fgq10 = 0.;
  double fqq11 = 0., fgg11 = 0., fgq11 = 0.;

  if (scheme & (NLO | NNLO | NLOAPPROX | NNLOAPPROX)) {
    // NLO scaling functions to calculate mu-dependence
    fgq10    = pxs_gq_10.eval(sparton);
    fgg10    = pxs_gg_10.eval(sparton);
    fqq10    = pxs_qq_10.eval(sparton);
    fgq11    = xs_gq_11();
    fgg11    = xs_gg_11();
    fqq11    = xs_qq_11();
  }

  /* 
   * Combine the scale dependent contributions with the scale
   * independent contributions for the two different approximations:
   */
  if (scheme & NLO) {
    nlobq = fourpi * ( fqq10 + Lmf * fqq11 );
    nlogg = fourpi * ( fgg10 + Lmf * fgg11 );
    nlogq = fourpi * ( fgq10 + Lmf * fgq11 + Lmr * fgq00 * b0 );
  } else if (scheme & NLOAPPROX) {
    // always use full scale dependency
    nlobq = fourpi * (       + Lmf * fqq11 );
    nlogg = fourpi * (       + Lmf * fgg11 );
    nlogq = fourpi * (         Lmf * fgq11 + Lmr * fgq00 * b0 )
      + approxgq1;
  }
  
  if (scheme & FAKE_NNLO) {
    /* This option is used to estimate what can be expected from the 
     * full NNLO calcu;ation as far as the scaledependence is concerned.
     * We estimate the size of the NNLO corrections for mu=mt using
     * a simple ansatz:
     */
    // - scale NLO result with ratio (NLO/LO)
    approxgq2 = fourpi * fgq10 * (fgq10 / fgq00);
    // - simply use LO result
    //approxgq2 = fgq00;
    // - model with sqrt(mt^2 / s)
    //approxgq2    = sqrt(rho / hcq) / (twopi*twopi);
  }
  
  /*
   * Calculate the full scale dependence
   * Whatever approximation is used we always include the complete
   * scale (in)dependence!
   */ 
  if (scheme & (NNLO | NNLOAPPROX | FAKE_NNLO)) {
    double fgg21, fbq21, fgq21, fgq21_h, fgqig21, fgqbig21, fgqbig21_h;
    double fgg22, fbq22, fgq22, fgq22_h, fgqig22, fgqbig22, fgqbig22_h;
    fgq21    = xs_gq_21();
    fgq22    = xs_gq_22();
    fgq21_h  = xs_gq_21_h();
    fgq22_h  = xs_gq_22_h();
    fgg21    = xs_gg_21();
    fgg22    = xs_gg_22();
    fbq21    = xs_qq_21();
    fbq22    = xs_qq_22();
    fgqig21  = xs_gqig_21();
    fgqig22  = xs_gqig_22();
    fgqbig21 = xs_gqbig2_21();
    fgqbig22 = xs_gqbig2_22();
    fgqbig21_h = xs_gqbig2_21_h();
    fgqbig22_h = xs_gqbig2_22_h();

    // always use full scale dependency
    nnlogq       = fourpi2 * Lmf * (fgq21      + Lmf * fgq22)
      + fourpi2 * Lmr * (fgq00 * b1
			 + 2. * b0 * (fgq10 + Lmf * fgq11)
			 + b0 * b0 * Lmr * fgq00);
    nnlogq_h     = fourpi2 * Lmf * (fgq21_h    + Lmf * fgq22_h)
      + fourpi2 * Lmr * (fgq00 * b1
			 + 2. * b0 * (fgq10 + Lmf * fgq11)
			 + b0 * b0 * Lmr * fgq00);
    nnlogg       = fourpi2 * Lmf * (fgg21      + Lmf * fgg22)
      + fourpi2 * Lmr * (2. * b0 * (fgg10 + Lmf * fgg11));
    nnlobq       = fourpi2 * Lmf * (fbq21      + Lmf * fbq22)
      + fourpi2 * Lmr * (2. * b0 * (fqq10 + Lmf * fqq11));
    nnlogqig     = fourpi2 * Lmf * (fgqig21    + Lmf * fgqig22);
    nnlogqbig2   = fourpi2 * Lmf * (fgqbig21   + Lmf * fgqbig22);
    nnlogqbig2_h = fourpi2 * Lmf * (fgqbig21_h + Lmf * fgqbig22_h);
    // add scale-independent approx. NNLO terms
    if (scheme & (NNLOAPPROX | FAKE_NNLO))
      nnlogq  += approxgq2;
  }
}

double HathorSgTopWt::evaluateIntegral(double as, double wgt) {
  /*
   *  Calculate the full integrand, using the given alpha_s.
   */
  
  double lo = 0, nlo = 0, nnlo = 0;
  const double as2 = as * as;
  const double as3 = as * as2;

  double fac = jacobian * hcq * as * wgt;

  if (scheme & LO) {
    lo = gb * fgq00;
    if (particle == TOPQUARK) {
      contribution[2] += fac * gb * fgq00; // g q
    } else if (particle == ANTITOPQUARK) {
      contribution[3] += fac * gb * fgq00; // g qbar
    }
  }
  
  fac *= as;

  if (scheme & (NLO | NLOAPPROX)) {
    nlo = gb * nlogq + gg * nlogg + (bq + bqbar) * nlobq;
    if (particle == TOPQUARK) {
      contribution[2] += fac * gb * nlogq;       // g q
      contribution[4] += fac * gg * nlogg;       // g g
      contribution[1] += fac * bqbar * nlobq;    // q qbar
      contribution[0] += fac * bq * nlobq;       // q q
    } else if (particle == ANTITOPQUARK) {
      contribution[3] += fac * gb * nlogq;       // g qbar
      contribution[4] += fac * gg * nlogg;       // g g
      contribution[5] += fac * bqbar * nlobq;    // qbar qbar
      contribution[1] += fac * bq * nlobq;       // q qbar
    }			
  }
  
  if (scheme & (NNLO | NNLOAPPROX | FAKE_NNLO)) {
    nnlo = gb * nnlogq -gb_h * nnlogq_h 
      + gg * nnlogg + (bq + bqbar) * nnlobq 
      + (gqig + gqbarig) * nnlogqig 
      + gqbarig2 * nnlogqbig2 - gqbarig2_h * nnlogqbig2_h ;
    if (particle == TOPQUARK) {
      contribution[2] += fac * gb * nnlogq; // g q
    } else if (particle == ANTITOPQUARK) {
      contribution[3] += fac * gb * nnlogq; // g qbar
    }
  }

  return jacobian * hcq * alpha /swq * ( as * lo + as2 * nlo + as3 * nnlo );
}

double HathorSgTopWt::xs_gq_00(double t, double u) {
  /*
   * Wt-channel amplitude squared for b g -> t W-
   *
   * Taken from hep-ph/0609289, Eq.(A.13)
   */
  const double s = sparton;
  double u1 = u - m2, u2 = u - mw2, t1 = -m2 + mw2 - s - u1;
  
  double A1 = -u2 * (s - m2 - mw2) * (mw2 + m2/2.) 
    -(t1/2.) * (-2.*mw2*mw2 + mw2*m2 + m4)
    -2.*u2*m2*(2.*mw2 + m2);
  
  double A2 = -t1*(-mw2+m2)*mw2 - (u2/2.)*t1*m2 - u1*(u2/2.)*m2
    -s*mw2*m2 - (s/2.)*m4;

  double A3 = -u1*(s/4.)*(2.*mw2 + m2);
  
  return( 4. * pi2 / (3. * mw2) * (A1/(u1*u1) - 2.*A2/(u1*s) + 2.*A3/(s*s)) );
}

void HathorSgTopWt::approxNNLO(double t, double s4, double result[8]) {
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
  double u = s4 - s - t + mw2 + m2;
  
  for (int i = 0; i < 8; i++)
    result[i] = 0.;
  
  // calculate Born term
  double FB = xs_gq_00(t, u) / (16. * pi * s*s);
  
  if (llflags & NLL_CALCULATION) {
    /*
     * NLL results from hep-ph/0609287
     */ 
    double c3 = 0., c2 = 0., c2mu = 0., T2 = 0., c1mu = 0., G1S11 = 0.;
    // 1-loop soft anomalous dimension (in axial gauge)
    G1S11 = CF*log((m2-t)/(m*sqrts)) + 0.5*CA*log((m2-u)/(m2-t)) + CA/2.;
    // everything else
    c3 = 2.*(CF + CA);
    T2 = 2.*G1S11 - CF - CA - 2.*CF*log((-u+mw2)/m2)
      - 2.*CA*log((-t+mw2)/m2) - (CF+CA)*log(m2/s);
    //c2mu = -(CF + CA)*Lmf;
    //c1mu = (CF*log((-u+mw2)/m2) + CA*log((-t+mw2)/m2) - 0.75*CF 
    //	   - beta0/4.)*Lmf + 0.25*beta0*Lmrm;
    c2 = T2 + c2mu;
    
    if (scheme & NLOAPPROX) {
      // "constant term" 
      result[0] += FB*invpi * c1mu;
      // Coefficient of [1/s4]_+ Eq.(3.1)
      if (llflags & NLO_NLL)
	result[1] += FB*invpi * c2;
      // Coefficient of [log(s4/M^2)/s4]_+ Eq.(3.1)
      if (llflags & NLO_LL)
	result[2] += FB*invpi * c3;
    }
    
    if (scheme & NNLOAPPROX) {
      // Coefficient of [1/s4]_+  (Eq.(3.12))
      // NNLO-NNNLL: not all terms are known yet
      if (llflags & NNLO_NNNLL)
	result[4] += FB*invpi2 * ( - zeta2*c3*c2 + zeta3*c3*c3);
      //    (c2mu*c1mu + beta0/4.*c2mu*log(mur2/m2) + 
      //    (CF+CA)*beta0/8.*pow(log(muf2/m2),2.)
      // Coefficient [log(s4/M^2)/s4]_+ (Eq.(3.12))
      if (llflags & NNLO_NNLL)
	result[5] += FB*invpi2 * ( - zeta2*c3*c3);
      //    c3*c1mu + c2mu*c2mu + 2.*c2mu*T2 + beta0/4.*c3*log(mur2/m2)
      // Coefficient [log^2(s4/M^2)/s4]_+ (Eq.(3.12))
      if (llflags & NNLO_NLL)
	result[6] += FB*invpi2 * (1.5*c3*c2 - beta0/4.*c3);
      // Coefficient [log^3(s4/M^2)/s4]_+ (Eq.(3.12))
      if (llflags & NNLO_LL)
	result[7] += FB*invpi2 * 0.5*c3*c3;
    }
  } else {
    /*
     * Results from arXiv 1005.4451
     */ 
    double c3 = 0., c2 = 0., c2mu = 0., T2 = 0., c1 = 0., c1mu = 0., G1S11 = 0.;
    // 1-loop soft anomalous dimension (in Feynman gauge)
    G1S11 = CF*(log((m2-t)/(m*sqrts)) - 0.5) + 0.5*CA*log((m2-u)/(m2-t));
    // everything else
    c3 = 2.*(CF + CA);
    T2 = -2.*CF*log((mw2-u)/m2) - 2.*CA*log((mw2-t)/m2) + 2.*G1S11 
      - (CF+CA)*log(m2/s);
    //c2mu = -(CF + CA)*Lmf;
    //c1mu = (CF*log((-u+mw2)/m2) + CA*log((-t+mw2)/m2) - 0.75*CF 
    //	   - beta0/4.)*Lmf + 0.25*beta0*Lmrm;
    c2 = T2 + c2mu;
    c1 = c1mu;
    
    if (scheme & NLOAPPROX) {
      // "constant term" 
      result[0] += FB*invpi * c1mu;
      // Coefficient of [1/s4]_+ (Eq.(3.7))
      if (llflags & NLO_NLL)
	result[1] += FB*invpi * c2;
      // Coefficient of [log(s4/M^2)/s4]_+ (Eq.(3.7))
      if (llflags & NLO_LL)
	result[2] += FB*invpi * c3;
    }
    
    if (scheme & NNLOAPPROX) {
      // Coefficient of [1/s4]_+ (Eq.(3.12))
      // NNLO-NNNLL: not all terms are known yet
      if (llflags & NNLO_NNNLL) {
	const double G2S11 = 0.5*K*G1S11 + CF*CA*(1. - zeta3)/4.;
	result[4] += FB*invpi2 * (c2*c1 - zeta2*c3*c2 + zeta3*c3*c3 
				  + 0.25*beta0*c2*log(mur2/s)
				  - 0.5*beta0*CF*pow(log((mw2-u)/m2),2.) 
				  - 0.5*beta0*CA*pow(log((mw2-t)/m2),2.)
				  - CF*K*log((mw2-u)/m2) 
				  - CA*K*log((mw2-t)/m2) + D2 + D2g
				  + 0.125*beta0*(CF+CA)*pow(log(muf2/s),2.) 
				  - 0.5*(CF+CA)*K*log(muf2/s)
				  + 2.*G2S11);
      }
      // Coefficient of [log(s4/M^2)/s4]_+ (Eq.(3.12))
      if (llflags & NNLO_NNLL)
	result[5] += FB*invpi2 * (c3*c1 + c2*c2 - zeta2*c3*c3 - 0.5*beta0*T2
				  // + 0.25*beta0*c3*Lmrm 
				  + CF*K + CA*K);
      // Coefficient of [log^2(s4/M^2)/s4]_+ (Eq.(3.12))
      if (llflags & NNLO_NLL)
	result[6] += FB*invpi2 * (1.5*c3*c2 - beta0/4.*c3);
      // Coefficient of [log^3(s4/M^2)/s4]_+ (Eq.(3.12))
      if (llflags & NNLO_LL)
	result[7] += FB*invpi2 * 0.5*c3*c3;
    }
  }
}

void HathorSgTopWt::evaluateApproxNNLO() {
  // calculate kinematics
  const double s = sparton;
  const double m32 = mw2;
  double t0   = 0.5*(-s + m32 + m2);
  double tpm  = 0.5*sqrt(s*s - 2.*s*(m32 + m2) + pow(m32 - m2,2.));
  double tmax = t0 + tpm;
  double tmin = t0 - tpm;
  double t  = (tmax - tmin) * rt + tmin;
  double s4max = s + t - m2 + m32*s/(t-m32);
  double s4 = s4max * rs4;
  // jacobian from transformation min..max -> 0..1
  double jac = (tmax - tmin) * s4max;
  
  double f_s4[8],f_0[8];

  approxNNLO(t, s4,f_s4);
  approxNNLO(t, 0.,f_0);
  
  // evaluate plus-distributions
  double logs1 = 0., logs2 = 0.;
  if (scheme & NLOAPPROX) {
    // no plus-distribution
    logs1 += f_0[0] / s4max;
    // [1/s4]_+
    logs1 += (f_s4[1] - f_0[1]) *1. / s4     
      +         log(s4max/m2)         * f_0[1] / s4max;
    // [log(s4/M^2)/s4]_+
    logs1 += (f_s4[2] - f_0[2]) * log(s4/m2)/s4     
      + 0.5 *   pow(log(s4max/m2),2.) * f_0[2] / s4max;
  }
  if (scheme & NNLOAPPROX) {
    // no plus-distribution
    logs2 += f_0[3] / s4max;
    // [1/s4]_+
    logs2 += (f_s4[4] - f_0[4]) * 1. / s4     
      +  log(s4max/m2)         * f_0[4] / s4max;
    // [log(s4/M^2)/s4]_+
    logs2 += (f_s4[5] - f_0[5]) * log(s4/m2)/s4 
      + 0.5 *   pow(log(s4max/m2),2.) * f_0[5] / s4max;
    // [log^2(s4/M^2)/s4]_+
    logs2 += (f_s4[6] - f_0[6]) * pow(log(s4/m2),2.)/s4 
      + (1./3.)*pow(log(s4max/m2),3.) * f_0[6] / s4max;
    // [log^3(s4/M^2)/s4]_+
    logs2 += (f_s4[7] - f_0[7]) * pow(log(s4/m2),3.)/s4 
      + 0.25 *  pow(log(s4max/m2),4.) * f_0[7] / s4max;
  }
  
  
  approxgq1 = jac * logs1;
  approxgq2 = jac * logs2;
}

double HathorSgTopWt::xs_gq_11() {
  return( - 1. * xPqq0.convolve(rho, rt, pxs_gq_00)
	  - 1. * xPgg0.convolve(rho, rt, pxs_gq_00)
	  + 1. * b0 * pxs_gq_00.eval(sparton));
}

double HathorSgTopWt::xs_gg_11() {
  return( - 2. * xPqg0.convolve(rho, rt, pxs_gq_00) );
}

double HathorSgTopWt::xs_qq_11() {
  return( - 1. * xPgq0.convolve(rho, rt, pxs_gq_00) );
}

double HathorSgTopWt::xs_gq_21() {
  return( - 2.*nf * xPqg0.convolve(rho, rt, pxs_qq_10)
	  -  1. * xPgg1.convolve(rho, rt, pxs_gq_00)
	  -  1. * xPqq1v.convolve(rho, rt, pxs_gq_00)
	  -  1. * xPqq0.convolve(rho, rt, pxs_gq_10)
	  -  1. * xPgg0.convolve(rho, rt, pxs_gq_10)
	  +  2. * b0 * pxs_gq_10.eval(sparton)
	  +  1. * b1 * pxs_gq_00.eval(sparton));
}

double HathorSgTopWt::xs_gq_21_h() {
  return(- 3. * xPqg0.convolve(rho, rt, pxs_qq_10));
}


double HathorSgTopWt::xs_gq_22() {
  return( + 0.5 * xPqqqq.convolve(rho, rt, rs4, pxs_gq_00)
	  +  1. * xPqqgg.convolve(rho, rt, rs4, pxs_gq_00)
	  - 1.5 * b0 * xPqq0.convolve(rho, rt, pxs_gq_00)
	  +  nf * xPqggq.convolve(rho, rt, rs4, pxs_gq_00)
	  + 0.5 * xPgggg.convolve(rho, rt, rs4, pxs_gq_00)
	  - 1.5 * b0 * xPgg0.convolve(rho, rt, pxs_gq_00)
	  +  1. * b0 * b0 * pxs_gq_00.eval(sparton));
}

double HathorSgTopWt::xs_gq_22_h() {
  return( 1.5 * xPqggq.convolve(rho, rt, rs4, pxs_gq_00));
}

double HathorSgTopWt::xs_gg_21() {
  return( - 2. * xPqg1.convolve(rho, rt, pxs_gq_00)
	  - 2. * xPqg0.convolve(rho, rt, pxs_gq_10)
	  - 2. * xPgg0.convolve(rho, rt, pxs_gg_10)
	  + 2. * b0 * pxs_gg_10.eval(sparton));
}

double HathorSgTopWt::xs_gg_22() {
  return( + 1. * xPqqqg.convolve(rho, rt, rs4, pxs_gq_00)
	  + 3. * xPqggg.convolve(rho, rt, rs4, pxs_gq_00)
	  - 3. * b0 * xPqg0.convolve(rho, rt, pxs_gq_00));
}
    
double HathorSgTopWt::xs_qq_21() {
  return( - 2. * xPqq0.convolve(rho, rt, pxs_qq_10)
	  + 2. * b0 * pxs_qq_10.eval(sparton)
	  - 1. * xPgq1.convolve(rho, rt, pxs_gq_00)
	  - 1. * xPgq0.convolve(rho, rt, pxs_gq_10));
}

double HathorSgTopWt::xs_qq_22() {
  return( + 1.5 * xPqqgq.convolve(rho, rt, rs4, pxs_gq_00)
	  + 0.5 * xPgqgg.convolve(rho, rt, rs4, pxs_gq_00)
	  - 1.5 * b0 * xPgq0.convolve(rho, rt, pxs_gq_00));
}

double HathorSgTopWt::xs_gqig_21() {
  return( - 1. * xPqg0.convolve(rho, rt, pxs_qq_10)
	  - 1. * xPqq1s.convolve(rho, rt, pxs_gq_00)
	  - 1. * xPgq0.convolve(rho, rt, pxs_gg_10));
}

double HathorSgTopWt::xs_gqig_22() {
  return( 1.5 * xPqggq.convolve(rho, rt, rs4, pxs_gq_00));
}

double HathorSgTopWt::xs_gqbig2_21() {
  return( - 1. * xPqqb1v.convolve(rho, rt, pxs_gq_00));
}

double HathorSgTopWt::xs_gqbig2_21_h() {
  return( -1. * xPqg0.convolve(rho, rt, pxs_qq_10));
}


double HathorSgTopWt::xs_gqbig2_22() {
  return( 0. );
}

double HathorSgTopWt::xs_gqbig2_22_h() {
  return( + 0.5 * xPqggq.convolve(rho, rt, rs4, pxs_gq_00));
}



double HathorSgTopWt_LO_gq::eval(double s) {
  /*
   * partonic cross section for g b -> W- t
   */
  
  const double sqrts = sqrt(s);
  double E2 = (m2*m2 - 2.*m2*mw2 + mw4 - 2.*m2*s - 2.*mw2*s + s*s) / (4. * s);
  double E = sqrt(E2);
  
  double prefactor = pi * E / (3. * mw2)/ ( 2. * s * sqrts);
  double a = E * sqrts;
  double b = sqrt(E2 + m2)  * sqrts;
  double c = sqrt(E2 + mw2) * sqrts;
  
  double A1 = 0.5 *(m2 + 2.*mw2) 
    * ((m2 - mw2) * (- b - c + 2.*m2 + s) 
       / ((a + c + m2 - mw2) * (a - c - m2 + mw2))
       + (2.*m2 + s) / a * 0.5 
       * log((a + c + m2 - mw2) / (c - a + m2 - mw2)));
  double A2 = 1. / (2. * s) 
    * (m2 * (b + c + m2) + m2 * mw2 - 2.*mw4 
       - 1./a * (m2 + 2.*mw2)
       * ((m2 - mw2) * (b + c + m2 - mw2) - m2*s)
       * 0.5 * log((a + c + m2 - mw2) / (c - a + m2 - mw2)));
  double A3 = (c + m2 - mw2) * (m2 + 2.*mw2) / (2. * s);

  
  return(prefactor * ( A1 - 2. * A2 + A3 ) );
}

