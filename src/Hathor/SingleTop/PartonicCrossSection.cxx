#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

#include "AbstractHathor.h"
#include "PartonicCrossSection.h"

using namespace std;
/*
 *  The SgTopCrossSection class is intended to be subclassed.
 *  The subclasses should then implement the eval() method, which
 *  analytically calculates the partonic cross section.
 */

SgTopCrossSection::SgTopCrossSection() {
  /*
   *  Simple constructor setting default masses.
   */
  mw  = 80.385;
  mw2 = mw*mw;
}

void SgTopCrossSection::setMass(double _m, double _threshold) {
  /*
   *  Adjust the top mass.
   */
  m  = _m;
  m2 = m*m;
  threshold  = _threshold;
  threshold2 = threshold*threshold;
}



/*
 *  The SampledPartonicCrossSection class is intended to be instanciated.
 *  It can then be used to handle partonic cross sections, which have been
 *  sampled for different values of sqrt(s) and m_t.
 */


unsigned int SampledPartonicCrossSection::lastip = -1;



int SampledPartonicCrossSection::find(const int len, const double xvals[], 
				      const double x){

  // Do a simple binary search...
  int right,left,middle;

  left  = 0;
  right = len-1;
  while ( right-left > 1) {
    middle = (right+left) / 2;
    if ( x > xvals[middle] )
      left=middle;
    else
      right=middle;
  }
  return(left);
}


SampledPartonicCrossSection::SampledPartonicCrossSection(const double (&sqrts_)[NMASS][NE],
							 const double (&xs_)[NMASS][NE],
							 const double (&masses_)[NMASS])
  : sqrts(sqrts_), xs(xs_), masses(masses_) {
 
  scaleFactor = 1;
  nXSPoly = defaultXSPoly;
  
  // Don't rescale the partonic cross section by default.
  scaleFactor = 1.;
  
  // clear cache
  icache = 0;
  c_s[0]  = c_s[1]  = 0.;
  c_xs[0] = c_xs[1] = 0.;
}


void SampledPartonicCrossSection::setScalefactor(const double _scaleFactor) {
  /*
   *  Multiply the partonic cross section with the given factor.
   */
  
  scaleFactor = _scaleFactor;
}

void SampledPartonicCrossSection::setMass(const double _m, const double _thr) {
  /*
   *  Load sampled cross section for given mass.
   *  The _thr parameter should be set to the production threshold of the
   *  process (e.g. _thr=_m for s- and t-channel, _thr=_m+mw for Wt-channel).
   */
  
  // store configuration
  threshold  = _thr;
  threshold2 = _thr*_thr;
  m  = _m;
  m3 = threshold - m;
  
  // determine, which samples will be used
  findNearestMasses();
}

double SampledPartonicCrossSection::eval(double s) {
  /*
   *  Interpolate between the partonic cross sections
   *  to the top mass m and the given energy.
   */
  
  // check for cache hit
  // (huge performance boost when evaluating scale-dependent functions)
  if (c_s[0] == s) return c_xs[0];
  if (c_s[1] == s) return c_xs[1];

  const double sqx = sqrt(threshold2 / s);
  
  double result = 0.;
  switch(nMassPoly) {
  case 0: // exact sample available for chosen top mass
    result = getXS(imass, (masses[imass]+m3)/sqx);
    break;
  case 1: // perform interpolation
    if ( m < 165 ) {
      result = evalpoly1(m, 
			 masses[imass], 
			 getXS(imass, (masses[imass]+m3)/sqx),
			 masses[imass+1], 
			 getXS(imass+1, (masses[imass+1]+m3)/sqx));
    } else {
      result = evaldlog(m, 
			masses[imass], 
			getXS(imass, (masses[imass]+m3)/sqx),
			masses[imass+1], 
			getXS(imass+1, (masses[imass+1]+m3)/sqx));
    }
    break;
  default:
    AbstractHathor::err("Error in PartonicCrossSection.cxx: nMassPoly != 0,1 not implemented");
  }
  
  // update cache
  icache = 1 - icache;
  c_s[icache]  = s;
  c_xs[icache] = result * scaleFactor;
  
  return(result * scaleFactor);
}


double SampledPartonicCrossSection::getXS(const int imass, const double x) {
  /*
   *  Return the cross section for a given energy x.
   */
  
  //  cout <<"search for energy " << x << endl;;
  int isqrts = find(NE,&sqrts[imass][0],x);

  if ( nXSPoly == 1 ) {
    return evalpoly1(x, sqrts[imass][isqrts], xs[imass][isqrts],
		        sqrts[imass][isqrts+1], xs[imass][isqrts+1]);
  } else if (nXSPoly == 2) {
    return evalpoly2(x, sqrts[imass][isqrts], xs[imass][isqrts],
		     sqrts[imass][isqrts+1], xs[imass][isqrts+1],
		     sqrts[imass][isqrts+2], xs[imass][isqrts+2]);	
  } else {
    cout << "Unsupported degree of interpolation!" << endl;
    exit(1);
  }
  
}

void SampledPartonicCrossSection::findNearestMasses() {
  /*
   *  Find the required set of top masses closest to the chosen mass.
   */
  imass = find(NMASS,masses,m);

  //  cout << "imass = " << imass << " mass value = " << masses[imass] << endl;
  if ( masses[imass] == m ){
    nMassPoly = 0;
    return;
  } 
  
  // generic case
  nMassPoly = defaultMassPoly;
}



inline double SampledPartonicCrossSection::evalpoly1(double x, 
						     double x1,double y1, double x2, double y2) {
  /*
   *  Linear interpolation between two points.
   */
  const double m = (y2 - y1) / (x2 - x1);
  const double n = y1 - m * x1;
  return m * x + n;
}

inline double SampledPartonicCrossSection::evaldlog(double x, 
						    double x1, double y1, double x2, double y2) {
  /*
   *  Double-logarithmic interpolation between two points.
   */
  
  const double base = y1 / y2;
  const double expo = log(x / x1) / log(x1 / x2);
  if (base < 0.) {
    return evalpoly1(x, x1, y1, x2, y2);
  }
  return y1 * pow(base, expo);
}

inline double SampledPartonicCrossSection::evalpoly2(double x, 
	double x1, double y1, double x2, double y2, double x3, double y3) {
  /*
   *  Quadratic interpolation between three points.
   */
  return (x-x2)*(x-x3)/((x1-x2)*(x1-x3))*y1 +
    (x-x1)*(x-x3)/((x2-x1)*(x2-x3))*y2 +
    (x-x1)*(x-x2)/((x3-x1)*(x3-x2))*y3;
}

inline double SampledPartonicCrossSection::evalpoly3(double x, 
						     double x1, double y1, double x2, double y2, 
						     double x3, double y3, double x4, double y4) {
  /*
   *  Cubic interpolation between four points.
   */
  return (x-x2)*(x-x3)*(x-x4)/((x1-x2)*(x1-x3)*(x1-x4))*y1 +
    (x-x1)*(x-x3)*(x-x4)/((x2-x1)*(x2-x3)*(x2-x4))*y2 +
    (x-x1)*(x-x2)*(x-x4)/((x3-x1)*(x3-x2)*(x3-x4))*y3 +
    (x-x1)*(x-x2)*(x-x3)/((x4-x1)*(x4-x2)*(x4-x3))*y4;
}
