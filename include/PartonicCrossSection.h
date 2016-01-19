#ifndef PARTONIC_CROSS_SECTION_H_
#define PARTONIC_CROSS_SECTION_H_

#include <string>
#include <vector>

#define NMASS 18
#define NE 64

class PartonicCrossSection {
  
 public:
  virtual double eval(double s) = 0;
  double evalr(double rho) { return eval(threshold2/rho); };
  
 protected:
  double threshold, threshold2;
};


class SgTopCrossSection : public PartonicCrossSection {
  
 public:
  SgTopCrossSection();
  void setMass(const double _m, const double _threshold);
  virtual double eval(double s) = 0;
  
 protected:
  double m,  m2;
  double mw, mw2;
};


class SampledPartonicCrossSection : public PartonicCrossSection {
  
 public:
  SampledPartonicCrossSection(const double (&sqrts_)[NMASS][NE],
			      const double (&xs_)[NMASS][NE],
			      const double (&masses_)[NMASS]);
  
  void setScalefactor(const double scale);
  void setMass(const double _m, const double _threshold);
  virtual double eval(double s);
  
 protected:
  double m, m3;
  
 private:
  const double (&sqrts)[NMASS][NE];
  const double (&xs)[NMASS][NE];
  const double (&masses)[NMASS];
  // degree of sqrt(s) interpolation
  const static unsigned int defaultXSPoly   = 2;
  unsigned int nXSPoly;
  // degree of m_t interpolation
  const static unsigned int defaultMassPoly = 1;
  unsigned int nMassPoly;
  
  int imass;
  // factor to multiply the partonic cross section with
  double scaleFactor;
  
  static unsigned int lastip;

  // cache for the two most recent results
  int icache;
  double c_s[2], c_xs[2];
  
  void    findNearestMasses();
  static int find(const int len, const double xvals[],const double x);
  double  getXS(const int imass, const double sqrts);
  
  
  // interpolation routines
  static inline double evaldlog(double x, 
			 double x1, double y1, double x2, double y2);
  static inline double evalpoly1(double x, 
			  double x1, double y1, double x2, double y2);
  static inline double evalpoly2(double x,
			  double x1, double y1, double x2, double y2,
			  double x3, double y3);
  static inline double evalpoly3(double x,
			  double x1, double y1, double x2, double y2, 
			  double x3, double y3, double x4, double y4);
  
};


#endif // PARTONIC_CROSS_SECTION_H_
