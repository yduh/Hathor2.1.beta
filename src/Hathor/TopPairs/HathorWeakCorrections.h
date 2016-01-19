/* $Modified: Thu Jan 23 13:33:40 2014 by uwer $ */
#ifndef HATHORWEAKCORRECTIONS_H_
#define  HATHORWEAKCORRECTIONS_H_
#include <complex>
#include <string>

class WeakCorrections {
 public:
  WeakCorrections(double & mb_, double & mz_,double & mw_,double & mh_,
		  double & alpha_,double & swq_,double & hcq_);
  void printParameters();
  void setAlphas(double value){alphas = value;}
  void setAlpha(double value){alpha = value;}
  void setHcq(double value){hcq = value;}
  void setLambdat(double value);
  void setSwq(double value);
  void updateCouplings(void);
  double dsigmagg(const double mt, const double sparton, const double z);
  double dsigmaqq(const double mt, const double sparton, const double z);
  double dsigmaWeakgg(const double mt, const double shat, const double z);
  void dsigmaWeakqq(const double mt, const double shat, const double z,
		    double & up, double & down);
  void check();

 private:
  int ggtriangle,ggself,ggvertex,ggbox;
  double N;
  double GammaH,GammaHq;
  double &mb,&mz,&mw,&mh,&alpha,&swq,&hcq;
  double mt,mtq,mu2,alphas,cwq,gvt,gat,gab,gvb,lambdat,lambdat2;
  double gvt2,gat2,gvb2,gab2,gw,gw2;

  void info(std::string s); 
  void setScale(double value){mu2 = value*value;}
  double ReC0m(const double shat, const double mq);
  double ImC0m(const double shat, const double mq);
  double A(double mq);
  double B(double p1q, double m0q, double m1q);
  double C(double p1q, double p2q, double p5q,
	   double m0q,double m1q, double m2q);
  double D(double p1q, double p2q, double p3q, double p4q, 
	   double p5q, double p7q, 
	   double m0q,double m1q, double m2q, double m3q);
  double diffB0(double pq, double m0, double m1);
  double ggtriangles(const double shat, double z);
  double Higgs_s_channel(const double shat, double z);
  double Z_Chi_s_channel(const double shat, double z);
  double ggselfenergies(const double shat, double z);
  double ggvertices(const double shat, double z);
  double ggboxes(const double shat, double z);

  double f1(double x);
  double F2xs(const double s, const double z, const double Fe, 
	      const double Fm, const double Fb);
  inline double mycast(std::complex<double> cint) {return cint.real();}

};
#endif
