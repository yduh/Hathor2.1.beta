// $Modified: Thu Jun 24 09:56:51 2010 by uwer $
#include "Hathor.h"
#include "HathorPdf.h"
#include "HathorWeakCorrections.h"
#include <iostream>
#include <fstream>


using namespace std;

void QQtest(){

  double mt=173.3;

  unsigned int scheme = Hathor::NNLO | Hathor::NOHIGHENERGY;
  unsigned int channels[] = { 
    Hathor::NOQG | Hathor::NOQQB | Hathor::NOQQ | Hathor::NOQQP | Hathor::NOQQPB,
    Hathor::NOGG | Hathor::NOQQB | Hathor::NOQQ | Hathor::NOQQP | Hathor::NOQQPB,
    Hathor::NOQG | Hathor::NOGG | Hathor::NOQQ | Hathor::NOQQP | Hathor::NOQQPB,
    Hathor::NOQG | Hathor::NOGG | Hathor::NOQQB | Hathor::NOQQP| Hathor::NOQQPB,
    Hathor::NOQG | Hathor::NOGG | Hathor::NOQQB | Hathor::NOQQ | Hathor::NOQQPB,      
    Hathor::NOQG | Hathor::NOGG | Hathor::NOQQB | Hathor::NOQQ | Hathor::NOQQP};

  Hathor::COLLIDERTYPE collider[]={Hathor::PP, Hathor::PPBAR};
  //string pdfname[] = {"MSTW2008nnlo68cl"};
  string pdfname[] = {"NNPDF30_nlo_as_0118"};


  for (int icollider=0; icollider<1; icollider++){

    if ( Hathor::PP == collider[icollider] )
      cout 
	<< "---------------------------- LHC ------------------------------"
	<< endl;
    else
      cout 
	<< "------------------------- Tevatron ----------------------------"
	<< endl;

    for (int ipdf =0; ipdf<1;ipdf++){
      cout << "Using " << pdfname[ipdf] << ":" << endl;
      double val,err,chi;
      Lhapdf pdf(pdfname[ipdf]);
      Hathor XS(pdf);
      
      XS.setColliderType(collider[icollider]);
      if ( Hathor::PP == collider[icollider] )
	XS.setSqrtShad(7000.);
      for(int is=0; is<5; is++){
	cout << "--------------------------------------------------------\n";
	XS.setScheme(  scheme | channels[is] );      
	XS.PrintOptions();
	for (int imur=-1;imur<2;imur++){
	  for (int imuf=-1;imuf<2;imuf++){
	    double mur = mt * pow(2,-imur);
	    double muf = mt * pow(2,-imuf);
	    XS.getXsection(mt,mur,muf);
	    XS.getResult(0,val,err,chi);
	    cout << mt << "\t" 
		 << XS.getAlphas(mur) << "\t" 
		 << mur << "\t"
		 << muf << "\t"
		 << val << "\t" << err << " "
		 << endl;
	  }
	}
	cout << "\n\n\n";
      }
    }
  }
}

void CFMNR(){

  /* 
   * Reproduce fixed order results
   * from Cacciari, Frixione, Mangano, Nason and Ridolfi [JHEP 0809:127,2008]
   */
  double mt=171., mur=171., muf=171.;

  unsigned int scheme = Hathor::LO  | Hathor::NLO;
  Hathor::COLLIDERTYPE collider[]={Hathor::PP, Hathor::PPBAR};
  string pdfname[] = {"cteq65","MRST2006nnlo"};

  ofstream out("CFMNRcheck.log");
  out << "---------------------------------------------------------------"
      << endl
      << "In the following we reproduce the values quoted in" << endl
      << "Cacciari, Frixione, Mangano, Nason, Ridolfi, [JHEP 0809:127,2008]" 
      << endl
      << "The quoted results are: " << endl
      << "alphas(mur), crosssection, uncertainty from numerical integration,"
      << endl 
      << "asymmetric pdf uncertainty" << endl
      << "We use mt=mur=muf=171. GeV" << endl
      << "---------------------------------------------------------------"
      << endl;

  for (int icollider=0; icollider<2; icollider++){
    if ( Hathor::PP == collider[icollider] )
      out 
	<< "---------------------------- LHC ------------------------------"
	<< endl;
    else
      out 
	<< "------------------------- Tevatron ----------------------------"
	<< endl;

    for (int ipdf =0; ipdf<2;ipdf++){
      out << "Using " << pdfname[ipdf] << ":" << endl;
      double val,err,chi,pdfup,pdfdown;      
      Lhapdf pdf(pdfname[ipdf]);
      Hathor XS(pdf);
      
      XS.setColliderType(collider[icollider]);
      if ( Hathor::PP == collider[icollider] )
	XS.setSqrtShad(14000.);
      XS.setScheme(  scheme );      
      XS.getXsection(mt,mur,muf);
      XS.getResult(0,val,err,chi);

      XS.setScheme(  scheme | Hathor::PDF_SCAN );
      XS.setPrecision(Hathor::LOW);
      XS.getXsection(mt,mur,muf);
      XS.getPdfErr(pdfup,pdfdown);
      
      out << XS.getAlphas(mur) << " " << val << " " << err << " "
	   << pdfup << " " << -pdfdown << " (pdf)"
	   << endl;
    }
  }
}

void scaledep(Hathor & XS, const double m, ofstream &out){
  double muf, mur,val,err,chi;
  for (int imur=-1;imur<2;imur++){
    for (int imuf=-1;imuf<2;imuf++){
      muf = m * pow(2.,imuf);
      mur = m * pow(2.,imur);
      XS.getXsection(m,mur,muf);
      XS.getResult(0,val,err,chi);
      out  << mur <<" " << muf  << " " << XS.getAlphas(mur) << " "
	   << val << " " << err << endl;
    }
  }
}


void LMU09(){
  /*
   * This example reproduces the central curve (NNLO) of the right plot in Fig. 5
   * of Ref. [Langenfeld, Moch, Uwer, PRD 80, 054009 (2009)].
   */
  double val,err,chi;

  //Lhapdf lhapdf("MSTW2008nnlo68cl");
  Lhapdf lhapdf("NNPDF30_nlo_as_0118");
  Hathor XS(lhapdf);

  XS.setColliderType(Hathor::PPBAR);
  XS.setScheme(Hathor::LO | Hathor::NLO | Hathor::NNLO | Hathor::MS_MASS );
  //  XS.setPrecision(Hathor::LOW);

  for(double mt = 140; mt < 181.; mt++ ){
    XS.getXsection(mt,mt,mt);
    XS.getResult(0,val,err,chi);
    cout << mt << " " << XS.getAlphas(mt) <<" "<< val << " " << err << endl;
  }
}

void PdfErrDemo(string pdfname, AbstractHathor::COLLIDERTYPE collider){


  string logname;
  logname = pdfname;
  if ( Hathor::PP == collider ) 
    logname += "_LHC.log";
  else
    logname += "_TEV.log";
  ofstream out(logname.c_str());
  double val,err,chi,mt=173.3,ecms;
  double mur = 2*mt, muf = mt/2.;
  double pdfup,pdfdown;
  Lhapdf lhapdf(pdfname);

  out << "Pdf set: " << pdfname << endl;
  
  Hathor XS(lhapdf);
  XS.setColliderType(collider);
  if ( Hathor::PP == collider ) 
    ecms = 7000.;
  else
    ecms = 1960.;
  out << "Centre of mass energy: " << ecms << endl;
  XS.setSqrtShad(ecms);
  XS.setScheme(Hathor::LO | Hathor::NLO | Hathor::NNLO 
	       | Hathor::PDF_SCAN 
	       //	       | Hathor::NOQG | Hathor::NOGG
	       //	       | Hathor::NOGG | Hathor::NOQQ | Hathor::NOQQP | Hathor::NOQQPB | Hathor::NOQQB
	       | Hathor::NOQG | Hathor::NOQQ | Hathor::NOQQP | Hathor::NOQQPB | Hathor::NOQQB
	       //| Hathor::PDF_SYM_ERR 
	       );
  XS.setPrecision(Hathor::HIGH);
  XS.getXsection(mt,mur,muf);

  for (int i=0; i< XS.getPdfNumber(); i++) {
    XS.getResult(i,val,err,chi);
    out << i << " " << val << " " << err << endl;
  } 

  XS.getPdfErr(pdfup,pdfdown);
  XS.getResult(0,val,err,chi); 
    out << val  
       << " + " << pdfup 
       << " - " << pdfdown
       << endl;
}

void CFM13(){
  double mt=173.3; 
  double val,err,chi;

  double ecms[4]={1960.,7000.,8000.,14000.};
  Hathor::COLLIDERTYPE ctype[]={Hathor::PPBAR,Hathor::PP,Hathor::PP,Hathor::PP};
  //Lhapdf lhapdf("MSTW2008nnlo68cl");
  Lhapdf lhapdf("NNPDF30_nlo_as_0118");
  Hathor XS(lhapdf);

  for(int i=0;i<4;i++){
    XS.setColliderType(ctype[i]);
    XS.setSqrtShad(ecms[i]);
    XS.setScheme(Hathor::LO | Hathor::NLO | Hathor::NNLO);
    XS.setPrecision(Hathor::HIGH);
    
    XS.getXsection(mt,mt,mt);
    XS.getResult(0,val,err,chi);
    cout << mt << " " << XS.getAlphas(mt) <<" "<< val << " " << err << endl;
  }
}

void nonSM(){

  //Lhapdf lhapdf("NNPDF30_nlo_as_0118");
  //WeakCorrections(double & mb_, double & mz_,double & mw_,double & mh_,
  //  double & alpha_,double & swq_,double & hcq_);	
  double Mb = 4.82, MZ = 91.1876, MW = 80.385, Mh = 126. ;
  double Alpha = 1/126.3, Sin2 = 1-(MW/MZ)*(MW/MZ), Hcq;
  WeakCorrections weak(Mb, MZ, MW, Mh, Alpha, Sin2, Hcq);
  
  double factor = 100.;
  weak.setLambdat(factor);

  //cout << "yukuwa const " << yukuwaC << endl;
   double mt=173.3;
  //double mt=500.;
  double val,err,chi;

  // double ecms[4]={1960.,7000.,8000.,14000.};
  double ecms[1] = {13000.};
  // Hathor::COLLIDERTYPE ctype[]={Hathor::PPBAR,Hathor::PP,Hathor::PP,Hathor::PP};
  Hathor::COLLIDERTYPE ctype[]={Hathor::PP};
  //Lhapdf lhapdf("MSTW2008nnlo68cl");
  Lhapdf lhapdf("NNPDF30_nlo_as_0118");
  Hathor XS(lhapdf);

  for(int i=0;i<1;i++){
    XS.setColliderType(ctype[i]);
    XS.setSqrtShad(ecms[i]);
    XS.setScheme(Hathor::LO | Hathor::NLO | Hathor::NNLO);
    XS.setPrecision(Hathor::HIGH);
    
    XS.getXsection(mt,mt,mt);
    XS.getResult(0,val,err,chi);

	double kin[4] = {300., 400., 500., 600.};
	double weight = 1;
	double results[4];

    cout << mt << " " << XS.getAlphas(mt) <<" "<< val << " " << err << endl;
	
	
	//AbstractHathor *dXS(lhapdf);
	//dXS = &XS;
	//cout << dXS->f(kin, weight, results) << endl;
  }

}
int main(){
  //  CFMNR();
  //QQtest();
  // PdfErrDemo("MSTW2008nnlo68cl", Hathor::PP);
  //PdfErrDemo("abm11_5n_nnlo", Hathor::PP);
  //PdfErrDemo("MSTW2008nnlo68cl", Hathor::PPBAR);
  //PdfErrDemo("abm11_5n_nnlo", Hathor::PPBAR);
  //CFM13();
  nonSM();
}
