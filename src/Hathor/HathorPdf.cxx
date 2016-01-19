// $Modified: Tue Ma 25 14:43:23 20: Mon May 31 19:29:55 2010 by uwer 10 by uwer $
#include "AbstractHathor.h"
#include "HathorPdf.h"
#include <string>
#include <string.h>
#include "LHAPDF/LHAPDF.h"
#include <cstdlib>


using namespace std;
using namespace LHAPDF;


#if LHAPDF_MAJOR_VERSION > 5

Lhapdf::Lhapdf(){
 try {
  pdfset = getPDFSet("MSTW2008nnlo68cl");
 } catch (LHAPDF::ReadError) {
  cerr << "HathorPdf: Cannot load PDF set from LHAPDF." << endl;
  cerr << "HathorPdf: Please check the PDF name 'MSTW2008nnlo68cl' and" << endl;
  cerr << "HathorPdf: ensure its presence in your LHAPDF setup." << endl;
  exit(1);
 }
 pdfs = pdfset.mkPDFs();
 pdf = pdfs[0];
}

Lhapdf::Lhapdf(string s){
 try {
  pdfset = getPDFSet(s);
 } catch (LHAPDF::ReadError) {
  cerr << "HathorPdf: Cannot load PDF set from LHAPDF." << endl;
  cerr << "HathorPdf: Please check the PDF name '" << s << "' and" << endl;
  cerr << "HathorPdf: ensure its presence in your LHAPDF setup." << endl;
  exit(1);
 }
 pdfs = pdfset.mkPDFs();
 pdf = pdfs[0];
}

Lhapdf::~Lhapdf(){
 for(unsigned int i = 0; i < pdfs.size(); i++)
  delete pdfs.at(i);
 pdfs.clear();
}

void Lhapdf::GetPdf(double x, double muf, double h[13]){
 vector<double> res;
 res.resize(13);
 pdf->xfxQ(x, muf, res);
 h[AbstractHathor::GLUON] = res[6+ 0];
 h[AbstractHathor::ABOTTOM] = res[6 -5];
 h[AbstractHathor::ACHARM] = res[6 -4];
 h[AbstractHathor::ASTRANGE] = res[6 -3];
 h[AbstractHathor::AUP] = res[6 -2];
 h[AbstractHathor::ADOWN] = res[6 -1];
 h[AbstractHathor::DOWN] = res[6+ 1];
 h[AbstractHathor::UP] = res[6+ 2];
 h[AbstractHathor::STRANGE] = res[6+ 3];
 h[AbstractHathor::CHARM] = res[6+ 4];
 h[AbstractHathor::BOTTOM] = res[6+ 5];
 
 for (int i=0; i<13; i++){
  h[i] /= x;
 }
}

double Lhapdf::GetAlphas(double mu){
 return(pdf->alphasQ(mu));
}

void Lhapdf::InitMember(int i){ 
 pdf = pdfs[i];
}

int Lhapdf::NumberPdf(void){ 
 return(pdfset.size()-1);
}

string Lhapdf::GetName(void){
 return(pdfset.name());
}

#else

Lhapdf::Lhapdf(){
  PDFname  = "MSTW2008nnlo68cl";
  initPDFSet(PDFname, LHGRID);
}

Lhapdf::Lhapdf(string s){
  if ( s.find(".") == string::npos ){
    PDFname  = s;
    initPDFSet(PDFname, LHGRID);
    return;
  } else {
    int n = s.find('.');
    string hlp(s,0,n);
    string extension(s,n+1);
    PDFname = hlp;
    if ( "LHpdf" == extension ){       
      initPDFSet(PDFname, LHPDF);
      return;
    }
    if ( "LHgrid" == extension ) {
      initPDFSet(PDFname, LHGRID);
      return;
    }
  }

  cerr << "HathorPdf: Cannot load LHA pdf, unknown extension\n";
  cerr << "HathorPdf: use .LHpdf or .LHgrid to specify the respective files\n";
  cerr << "HathorPdf: Without extension .LHgrid files are used.\n"; 
  exit(1);
}

void Lhapdf::GetPdf(double x, double muf, double h[13]){
  

  double res[13];
  xfx(x, muf, res);
  h[AbstractHathor::GLUON] = res[6+ 0];
  h[AbstractHathor::ABOTTOM] = res[6 -5];
  h[AbstractHathor::ACHARM] = res[6 -4];
  h[AbstractHathor::ASTRANGE] = res[6 -3];
  h[AbstractHathor::AUP] = res[6 -2];
  h[AbstractHathor::ADOWN] = res[6 -1];
  h[AbstractHathor::DOWN] = res[6+ 1];
  h[AbstractHathor::UP] = res[6+ 2];
  h[AbstractHathor::STRANGE] = res[6+ 3];
  h[AbstractHathor::CHARM] = res[6+ 4];
  h[AbstractHathor::BOTTOM] = res[6+ 5];

  for (int i=0; i<13; i++){
    h[i] /= x;
  }

}


double Lhapdf::GetAlphas(double mu){
  return(alphasPDF(mu) );
}

void Lhapdf::InitMember(int i){ 
  usePDFMember(i);
}

int Lhapdf::NumberPdf(void){ 
  return(numberPDF());
}

string Lhapdf::GetName(void){
  return(PDFname);
}

#endif // LHAPDF_MAJOR_VERSION > 5


/******************************************************************/
extern "C" {
  double mstwgetonepdf_(char *,const int & iset,const double & x,
			const double &q, const int & f);
}

MSTW::MSTW() {
  PDFname = "mstw2008nnlo.68cl";
  if (getenv("MSTWPATH")!=NULL) {
    prefix = std::string(getenv("MSTWPATH"));
  } else {
    prefix = std::string("./");
  }
  strcpy(fullname,(prefix+PDFname).c_str());
  imember = 0;
}

void MSTW::GetPdf(double x, double muf, double h[13]){

  h[AbstractHathor::GLUON] = mstwgetonepdf_(fullname,imember,x,muf,0);
  h[AbstractHathor::ABOTTOM] = mstwgetonepdf_(fullname,imember,x,muf,-5);
  h[AbstractHathor::ACHARM] = mstwgetonepdf_(fullname,imember,x,muf,-4);
  h[AbstractHathor::ASTRANGE] = mstwgetonepdf_(fullname,imember,x,muf,-3);
  h[AbstractHathor::AUP] = mstwgetonepdf_(fullname,imember,x,muf,-2);
  h[AbstractHathor::ADOWN] = mstwgetonepdf_(fullname,imember,x,muf,-1);
  h[AbstractHathor::DOWN] = mstwgetonepdf_(fullname,imember,x,muf,1);
  h[AbstractHathor::UP] = mstwgetonepdf_(fullname,imember,x,muf,2);
  h[AbstractHathor::STRANGE] = mstwgetonepdf_(fullname,imember,x,muf,3);
  h[AbstractHathor::CHARM] = mstwgetonepdf_(fullname,imember,x,muf,4);
  h[AbstractHathor::BOTTOM] = mstwgetonepdf_(fullname,imember,x,muf,5);

  for (int i=0; i<13; i++){
    h[i] /= x;
  }
}

double MSTW::GetAlphas(double mu){ return(1.);}

void MSTW::InitMember(int i){ imember=i;}

int MSTW::NumberPdf(void){return(40);}

string MSTW::GetName(void){
  return(PDFname);
}

