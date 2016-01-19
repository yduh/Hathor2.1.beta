// $Modified: Thu Jun 24 09:56:51 2010 by uwer $
#include "Hathor.h"
#include "HathorPdf.h"
#include "HathorWeakCorrections.h"
#include <iostream>
#include <fstream>

#include <TFile.h>
//#include <TH1D.h>
#include <TH1F.h>
//#include <TH2D.h>
#include <TH2F.h>
#include <TGraph.h>


using namespace std;

double s_step = 1;
double z_step = 0.05;
double pt_step = 5;
double y_step = 0.1;
double beta_step = 0.1;

int s_low = 350;
int s_high = 1500;
int z_low = -1;
int z_high = 1;
int pt_low = 30;
int pt_high = 1000;
double y_low = -5;
double y_high = 5;
int beta_low = 0.0;
int beta_high = 2.0;

int binX = int((s_high-s_low)/s_step);
int binY = int((z_high-z_low)/z_step);
int binX_pt = int((pt_high-pt_low)/pt_step);
int binX_y = int((y_high-y_low)/y_step);
int binX_beta = int((beta_high-beta_low)/beta_step);


double mt=173;
double Mb = 4.82, MZ = 91.1876, MW = 80.385, Mh = 126. ;
double Alpha = 1/126.3, Sin2 = 1-(MW/MZ)*(MW/MZ), Hcq = 0.389379323e9;
double factor = 1.;


TFile f("yukawa.root", "recreate");

TH2F* hWXS = new TH2F("WXS", "WXS", binX, s_low, s_high, binY, z_low, z_high);
TH2F* hLOXS = new TH2F("LOXS", "LOXS", binX, s_low, s_high, binY, z_low, z_high);
TH2F* hXS = new TH2F("XS", "XS", binX, s_low, s_high, binY, z_low, z_high);
TH2F* hresult = new TH2F("result", "result", binX, s_low, s_high, binY, z_low, z_high);

TH1F* hWXS_X = new TH1F("WXS_X", "WXS_X", binX, s_low, s_high);
TH1F* hLOXS_X = new TH1F("LOXS_X", "LOXS_X", binX, s_low, s_high);
TH1F* hXS_X = new TH1F("XS_X", "XS_X", binX, s_low, s_high);
TH1F* hresult_X = new TH1F("result_X", "result_X", binX, s_low, s_high);

TH1F* hWXS_Y = new TH1F("WXS_Y", "WXS_Y", binY, z_low, z_high);
TH1F* hLOXS_Y = new TH1F("LOXS_Y", "LOXS_Y", binY, z_low, z_high);
TH1F* hXS_Y = new TH1F("XS_Y", "XS_Y", binY, z_low, z_high);
TH1F* hresult_Y = new TH1F("result_Y", "result_Y", binY, z_low, z_high);

//more kinematic variables
TH1F* hWXS_pt = new TH1F("WXS_pt", "WXS_pt", binX_pt, pt_low, pt_high);
TH1F* hLOXS_pt = new TH1F("LOXS_pt", "LOXS_pt", binX_pt, pt_low, pt_high);
TH1F* hXS_pt = new TH1F("XS_pt", "XS_pt", binX_pt, pt_low, pt_high);
TH1F* hXSR_pt = new TH1F("XSR_pt", "XSR_pt", binX_pt, pt_low, pt_high);

TH1F* hWXS_y = new TH1F("WXS_y", "WXS_y", binX_y, y_low, y_high);
TH1F* hLOXS_y = new TH1F("LOXS_y", "LOXS_y", binX_y, y_low, y_high);
TH1F* hXS_y = new TH1F("XS_y", "XS_y", binX_y, y_low, y_high);
TH1F* hXSR_y = new TH1F("XSR_y", "XSR_y", binX_y, y_low, y_high);
TH1F* hdensity = new TH1F("density", "density", binX_y, y_low, y_high);
TH2F* hLOXS_s_y = new TH2F("LOXS_s_y", "LOXS_s_y", binX, s_low, s_high, binX_y, y_low, y_high);
TH2F* hXS_s_y = new TH2F("XS_s_y", "XS_s_y", binX, s_low, s_high, binX_y, y_low, y_high);
TH2F* hXSR_s_y = new TH2F("XSR_s_y", "XSR_s_y", binX, s_low, s_high, binX_y, y_low, y_high);

TH1F* hWXS_beta = new TH1F("WXS_beta", "WXS_beta", binX_beta, beta_low, beta_high);
TH1F* hLOXS_beta = new TH1F("LOXS_beta", "LOXS_beta", binX_beta, beta_low, beta_high);
TH1F* hXS_beta = new TH1F("XS_beta", "XS_beta", binX_beta, beta_low, beta_high);
TH1F* hXSR_beta = new TH1F("XSR_beta", "XSR_beta", binX_beta, beta_low, beta_high);



void diffMtt(){

  Lhapdf pdf("NNPDF30_nlo_as_0118");
  Hathor weak(pdf);
  
  double ggW, qqW, ggW_test;
  double ggLO, qqLO, ggLO_test;
  double upflavor, downflavor;

  //double ecms[1] = {13000.};
 // Hathor::COLLIDERTYPE ctype[]={Hathor::PP};

/*  for(int i=0;i<1;i++){
    XS.setColliderType(ctype[i]);
    XS.setSqrtShad(ecms[i]);
    XS.setScheme(Hathor::LO | Hathor::NLO | Hathor::NNLO);
    XS.setPrecision(Hathor::HIGH);
    
    XS.getXsection(mt,mt,mt);
    XS.getResult(0,val,err,chi);

    cout << mt << " " << XS.getAlphas(mt) <<" "<< val << " " << err << endl;		
  }
*/
  //tune the yukawa coupling const.
  WeakCorrections weak2(Mb, MZ, MW, Mh, Alpha, Sin2, Hcq);
  weak2.setLambdat(factor);

  for(int z_bin = 0; z_bin < binY; z_bin++){
	for(int s_bin = 0; s_bin < binX; s_bin++){
	      double s = s_bin* s_step + s_low;//[350,800]
		  double z = z_bin* z_step + z_low;//[-1,1]
		  ggW = weak2.dsigmaWeakgg(mt, s*s, z);
		  //weak2.dsigmaWeakqq(mt, s*s, z, upflavor, downflavor);
		  //cout<<upflavor<<" "<<downflavor<<endl;
		  ggLO = weak2.dsigmagg(mt, s*s, z);
		  //qqLO = weak2.dsigmaqq(mt, s*s, z);
	      //cout<< mt <<", "<< z <<", "<< s <<", "<< ggLO <<", "<< ggW << endl;
	      //cout<< mt <<", "<< z <<", "<< s <<", "<< qqLO <<", "<< qqW << endl;

		  hWXS->SetBinContent(s_bin+1, z_bin+1, ggW);
		  hWXS_X->SetBinContent(s_bin+1, ggW);
		  hWXS_Y->SetBinContent(z_bin+1, ggW);
		  
		  hLOXS->SetBinContent(s_bin+1, z_bin+1, ggLO);
		  hLOXS_X->SetBinContent(s_bin+1, ggLO);
		  hLOXS_Y->SetBinContent(z_bin+1, ggLO);


		  hXS->SetBinContent(s_bin+1, z_bin+1, ggLO+ggW);
		  hXS_X->SetBinContent(s_bin+1, ggLO+ggW);
		  hXS_Y->SetBinContent(z_bin+1, ggLO+ggW);

		  hresult->SetBinContent(s_bin+1, z_bin+1, ggW/ggLO);
		  hresult_X->SetBinContent(s_bin+1, ggW/ggLO);
		  hresult_Y->SetBinContent(z_bin+1, ggW/ggLO);
	  }
  }

  //hXS_X->Clone();
  //hresult_X->Scale(1/hresult_X->GetEntries());
  //for(int i=1; i<=binX; i++){
  //cout<<i<<" = "<<hresult_X->GetBinContent(i)<<endl;
  //}

  hWXS->Draw();
  hWXS_X->Draw();
  hWXS_Y->Draw();
  hLOXS->Draw();
  hLOXS_X->Draw();
  hLOXS_Y->Draw();
  hXS->Draw();
  hXS_X->Draw();
  hXS_Y->Draw();
  hresult->Draw();
  hresult_X->Draw();
  hresult_Y->Draw();

  hWXS->Write();
  hWXS_X->Write();
  hWXS_Y->Write();
  hLOXS->Write();
  hLOXS_X->Write();
  hLOXS_Y->Write();
  hXS->Write();
  hXS_X->Write();
  hXS_Y->Write();
  hresult->Write();
  hresult_X->Write();
  hresult_Y->Write();

  f.Close();
}

void diffPt(){

  Lhapdf pdf("NNPDF30_nlo_as_0118");
  Hathor weak(pdf);
  
  double ggW, ggW_test;
  double ggLO, qqLO, ggLO_test;
  double upflavor, downflavor;

  //tune the yukawa coupling const.
  WeakCorrections weak2(Mb, MZ, MW, Mh, Alpha, Sin2, Hcq);
  weak2.setLambdat(factor);

  for(int z_bin = 0; z_bin < binY; z_bin++){
	for(int pt_bin = 0; pt_bin < binX_pt; pt_bin++){
		  double z = z_bin* z_step + z_low;//[-1,1]
		  double s = 2*std::sqrt(std::pow(pt_bin*pt_step+pt_low, 2)/(1-z*z)+mt*mt);
		  
		  ggW = weak2.dsigmaWeakgg(mt, s*s, z);
		  weak2.dsigmaWeakqq(mt, s*s, z, upflavor, downflavor);
		  ggLO = weak2.dsigmagg(mt, s*s, z);
		  qqLO = weak2.dsigmaqq(mt, s*s, z);

		  //other kinematic variables:
		  //double pt = std::sqrt((s*s/4 - mt*mt)*(1 - z*z));
	      //x[z_bin] = pt;
	      //y[z_bin] = ggLO_test;
	      cout<< mt <<", "<< z <<", "<< s <<", "<< pt_bin <<", "<< ggLO <<", "<< ggW << endl;

		  hWXS_pt->SetBinContent(pt_bin+1, ggW);
		  hLOXS_pt->SetBinContent(pt_bin+1, ggLO);
		  hXS_pt->SetBinContent(pt_bin+1, ggW+ggLO);
		  hXSR_pt->SetBinContent(pt_bin+1, ggW/ggLO);
	  }
  }
		
  hWXS_pt->Draw();
  hWXS_pt->Write();
  hLOXS_pt->Draw();
  hLOXS_pt->Write();
  hXS_pt->Draw();
  hXS_pt->Write();
  hXSR_pt->Draw();
  hXSR_pt->Write();

  f.Close();
}

void diffDeltaY(){

  Lhapdf pdf("NNPDF30_nlo_as_0118");
  Hathor weak(pdf);
  
  double ggW, ggW_test;
  double ggLO, qqLO, ggLO_test;
  double upflavor, downflavor;

  //tune the yukawa coupling const.
  WeakCorrections weak2(Mb, MZ, MW, Mh, Alpha, Sin2, Hcq);
  weak2.setLambdat(factor);

  for(int y_bin = 0; y_bin < binX_y; y_bin++){
	for(int s_bin = 0; s_bin < binX; s_bin++){
		  double s = s_bin*s_step + s_low;
		  double z = 1/std::sqrt(1-4*(mt*mt)/(s*s)) * (std::exp(y_bin*y_step + y_low)-1)/(std::exp(y_bin*y_step + y_low)+1);
		  //cout<< y_bin*y_step+y_low <<endl;
		  //double density = 2*std::exp(y_bin*y_step + y_low)/std::pow(std::exp(y_bin*y_step + y_low) +1, 2);
		  double density = s*std::sqrt(s*s/4-mt*mt) / ((1-z*z)*s*s/4 + mt*mt*z*z);
		  density = 1/density;
		  /*double temp = std::pow(   (std::exp(y_bin*y_step + y_low) -1)   /   (std::exp(y_bin*y_step + y_low) +1)   ,2);
		  double missing_factor = 1 + 1/(2*mt)* std::pow(4*mt*mt/(s*s)-1, -1.5)* pow(z/ (z*z-temp), 2);//your missing part
		  density = density * missing_factor;
*/
		  if(std::abs(z)>1){
			  cout<<" >>> "<< s <<", "<< z <<", "<<density<<endl;
			  continue;
		  }

  /*for(int z_bin = 0; z_bin < binY; z_bin++){
	for(int y_bin = 0; y_bin < binX_y; y_bin++){
		  double z = z_bin* z_step + z_low;//[-1,1]
		  double y = y_bin* y_step + y_low;
		  double s =  1 - (std::exp(y)-1)/(std::exp(y)+1) * (1/z) ;
		  cout <<s<< ", "<<z<<endl;
		  s = 2*mt/s;
*/	 
		  ggW = weak2.dsigmaWeakgg(mt, s*s, z);
		  //cout << ggW << ", "<< ggW*density <<endl;
		  weak2.dsigmaWeakqq(mt, s*s, z, upflavor, downflavor);
		  ggLO = weak2.dsigmagg(mt, s*s, z);
		  qqLO = weak2.dsigmaqq(mt, s*s, z);
		  cout << s <<", "<< z <<", "<< ggLO <<", "<< ggW <<", "<<density<<endl;

		  hWXS_y->SetBinContent(y_bin+1, ggW*density);
		  hLOXS_y->SetBinContent(y_bin+1, ggLO*density);
		  hXS_y->SetBinContent(y_bin+1, (ggW+ggLO)*density);
		  hXSR_y->SetBinContent(y_bin+1, ggW/ggLO);
		  hdensity->SetBinContent(y_bin+1, density);
		  hLOXS_s_y->SetBinContent(s_bin+1, y_bin+1, ggLO*density);
		  hXS_s_y->SetBinContent(s_bin+1, y_bin+1, (ggW+ggLO)*density);
		  if(ggLO<0 || ggW ==0) continue;
		  hXSR_s_y->SetBinContent(s_bin+1, y_bin+1, ggW/ggLO);
	  }
  }
		
  hWXS_y->Draw();
  hWXS_y->Write();
  hLOXS_y->Draw();
  hLOXS_y->Write();
  hXS_y->Draw();
  hXS_y->Write();
  hXSR_y->Draw();
  hXSR_y->Write();

  hdensity->Draw();
  hdensity->Write();
  hLOXS_s_y->Draw();
  hLOXS_s_y->Write();
  hXS_s_y->Draw();
  hXS_s_y->Write();
  hXSR_s_y->Draw();
  hXSR_s_y->Write();


  f.Close();
}

void diffBeta(){

  Lhapdf pdf("NNPDF30_nlo_as_0118");
  Hathor weak(pdf);
  
  double ggW, ggW_test;
  double ggLO, qqLO, ggLO_test;
  double upflavor, downflavor;

  //tune the yukawa coupling const.
  WeakCorrections weak2(Mb, MZ, MW, Mh, Alpha, Sin2, Hcq);
  weak2.setLambdat(factor);

  for(int z_bin = 0; z_bin < binY; z_bin++){
	for(int beta_bin = 0; beta_bin < binX_beta; beta_bin++){
		  double z = z_bin* z_step + z_low;//[-1,1]
		  double beta = beta_bin* beta_step + beta_low;
		  double s = 2*mt* std::sqrt( 1/(1-std::pow(beta/4,2)) );
		  
		  ggW = weak2.dsigmaWeakgg(mt, s*s, z);
		  weak2.dsigmaWeakqq(mt, s*s, z, upflavor, downflavor);
		  ggLO = weak2.dsigmagg(mt, s*s, z);
		  qqLO = weak2.dsigmaqq(mt, s*s, z);
	      cout<< mt <<", "<< z <<", "<< s <<", "<< ggLO <<", "<< ggW << endl;

		  hWXS_beta->SetBinContent(beta_bin+1, ggW);
		  hLOXS_beta->SetBinContent(beta_bin+1, ggLO);
		  hXS_beta->SetBinContent(beta_bin+1, ggW+ggLO);
		  hXSR_beta->SetBinContent(beta_bin+1, ggW/ggLO);
	  }
  }
		
  hWXS_beta->Draw();
  hWXS_beta->Write();
  hLOXS_beta->Draw();
  hLOXS_beta->Write();
  hXS_beta->Draw();
  hXS_beta->Write();
  hXSR_beta->Draw();
  hXSR_beta->Write();

  f.Close();
}


int main(){

  //diffMtt();
  //diffPt();
  diffDeltaY();
  //diffBeta();

}
