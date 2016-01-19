/* $Modified: Tue May 18 14:44:01 2010 by uwer $ */
#ifndef PDF_H_
#define PDF_H_
#include "AbstractHathor.h"
#include "LHAPDF/LHAPDF.h"

class  MSTW : public Pdf {
  char fullname[255];  
  std::string prefix;
  std::string PDFname;
  int imember;
 public:
  MSTW();
  void GetPdf(double x, double muf, double h[13]);
  double GetAlphas(double mu);
  void InitMember(int i);
  int NumberPdf(void);
  std::string GetName(void);
};

class  CTEQ : public Pdf {
  char pdfname[255];
  int imember;
 public:
  CTEQ();
  void GetPdf(double x, double muf, double h[13]);
  double GetAlphas(double mu);
  void InitMember(int i);
  int NumberPdf(void);
};

#if LHAPDF_MAJOR_VERSION > 5
class Lhapdf : public Pdf {
  LHAPDF::PDFSet pdfset;
  std::vector<LHAPDF::PDF*> pdfs;
  const LHAPDF::PDF* pdf;
 public:
  Lhapdf();
  Lhapdf(std::string s);
  ~Lhapdf();
  void GetPdf(double x, double muf, double h[13]);
  double GetAlphas(double mu);
  void InitMember(int i); 
  int NumberPdf(void);
  std::string GetName(void);
};
#else
class Lhapdf : public Pdf {
  std::string PDFname;
 public:
  Lhapdf();
  Lhapdf(std::string s);
  void GetPdf(double x, double muf, double h[13]);
  double GetAlphas(double mu);
  void InitMember(int i); 
  int NumberPdf(void);
  std::string GetName(void);
};
#endif

#endif 
