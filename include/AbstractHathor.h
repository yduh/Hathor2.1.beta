#ifndef ABSTRACTPROCESS_H_
#define ABSTRACTPROCESS_H_

#include "HathorVegas.h"
#include <iostream>
#include <cstdlib>


class Pdf {
public:
  Pdf() { };
  virtual void GetPdf(double x, double muf, double h[13]) = 0;
  virtual double GetAlphas(double ) = 0;
  virtual int NumberPdf(void) = 0;
  virtual void InitMember(int ) = 0;
  virtual std::string GetName(void) = 0;
  // virtual int InitPdf(std::string) = 0;
};


class AbstractHathor {
	
public:
  enum APPROXSCHEMES { NNLOAPPROX=(1<<10), NLOAPPROX=(1<<16)};
  enum PERTURBATIVEORDER { LO=2, NLO=4, NNLO=8 };
  enum COLLIDERTYPE { PPBAR=0, PP=1 };
  enum ACCURACY { LOW=1000, MEDIUM=10000, HIGH=100000 };  
  enum SETTINGS { PDF_SCAN=32, PDF_SYM_ERR=64, NNPDF=128 };
  enum PARTONS  { GLUON = 0, 
		  DOWN=1, UP=2, STRANGE=3, CHARM=4, BOTTOM=5, TOP=6,
		  ADOWN=7, AUP=8, ASTRANGE=9, ACHARM=10, ABOTTOM=11, ATOP=12,
		  ZBOSON=20,WBOSON=21,HIGGS=22};

  virtual void setScheme(unsigned int newscheme);
  void sethc2(double hcq_);
  void setSqrtShad(double sqrts);
  void setColliderType(COLLIDERTYPE type);
  void setPDF(Pdf &newpdf);
  void setPrecision(int n);
  virtual void PrintOptions();
  
  double getAlphas(double mur) { return pdf.GetAlphas(mur); }
  double getXsection(double m, double mur, double muf);
  void getResult(int pdfset, double &integral, double &err){
    double dummy; getResult(pdfset,integral,err,dummy);}; 
  void getResult(int pdfset, double &integral, double &err, double & chi2a); 
  void getPdfErr(double &up, double &down);
  int getPdfNumber() { return pdfmax; }
  void setPdfNumber(int npdf);

  void setErrorFun(void (*fun)(int n, double res[], double &up, double &down)) 
  { evalerr = fun; }
  void clearErrorFun() { evalerr = 0; }

  AbstractHathor(Pdf &pdf_);

  static void info(std::string s);
  static void info(std::string s, double x);
  static void warn(std::string s);
  static void err(std::string s);

 protected:
  static double hcq;
  static bool debug;
  COLLIDERTYPE collider;
  Pdf &pdf;
  int calls;
  int pdfmax;
  unsigned int scheme;
  
  double delta;
  
  double shad, x1, x2;
  double m, m2;
  double mur, muf, mur2, muf2, Lmr, Lmf;
  
  virtual void update();
  
  
  virtual void setPartonicEnergy(const double x[]) = 0;
  virtual void evaluatePDFs(const double h1[], const double h2[], 
			    const double h2left[] = 0, 
			    const double h2right[] = 0) = 0;
  virtual void evaluateScalingFunctions() = 0;
  virtual double evaluateIntegral(double as, double wgt) = 0;
  
  int dimension;
  double integral[1024], error[1024], chi2a[1024];
  double as_pdf[1024];
  
 private:
  void (*evalerr)(int n, double res[], double &up, double &down); 
  
  // Functions used by Vegas
  double f(const double x[], const double wgt, double res[]);
  int getDimension() { return dimension; }
  int VegasIterations() { return 10; }
  int VegasPoints() { return calls; }
  double VegasAcc() { return 1.e-10; }
  
  static inline void ChargeConjugation(double h[13]) {
    // Charge conjugate proton pdf's to obtain anti-proton pdf's 
    double tmp;
    for (int i = DOWN; i < ADOWN; i++) {
      tmp = h[i];
      h[i] = h[i + ADOWN-DOWN];
      h[i + ADOWN-DOWN] = tmp;
    }
  }

  static double Li2(double x);
  friend class Vegas<AbstractHathor>;
  friend class WeakCorrections;
  friend class SplittingFunction; 
};

#endif // ABSTRACTPROCESS_H_
