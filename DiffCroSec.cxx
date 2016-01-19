// $Modified: Thu Jun 24 09:56:51 2010 by uwer $
#include "Hathor.h"
#include "HathorPdf.h"
#include "HathorWeakCorrections.h"
#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TH1D.h>
//#include <TH1F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TProfile.h>


using namespace std;

double s_step = 1;
double z_step = 0.02;
double pt_step = 5;
double y_step = 0.01;
double beta_step = 0.05;

int s_low = 346;
int s_high = 3000;
double z_low = -1;
double z_high = 1;
int pt_low = 30;
int pt_high = 1000;
double y_low = -5.69; //-4.88;
double y_high = 5.69; //4.88;
double beta_low = 0.3;
double beta_high = 2;

int binX = int((s_high-s_low)/s_step +0.5);
int binY = int((z_high-z_low)/z_step +0.5);
int binX_pt = int((pt_high-pt_low)/pt_step +0.5);
int binX_y = int((y_high-y_low)/y_step +0.5);
int binX_beta = int((beta_high-beta_low)/beta_step +0.5);


double mt=173;
double Mb = 4.82, MZ = 91.1876, MW = 80.385, Mh = 126. ;
double Alpha = 1/126.3, Sin2 = 1-(MW/MZ)*(MW/MZ), Hcq = 0.389379323e9;
double gt = 0;


TFile f("DiffCroSec.root", "recreate");
/*
TH2D* hWXS = new TH2D("WXS", "WXS", binX, s_low, s_high, binY, z_low, z_high);
TH2D* hLOXS = new TH2D("LOXS", "LOXS", binX, s_low, s_high, binY, z_low, z_high);
TH2D* hXSR = new TH2D("XSR", "XSR", binX, s_low, s_high, binY, z_low, z_high);
TH1D* hWXS_x = new TH1D("WXS_x", "WXS_x", binX, s_low, s_high);
TH1D* hLOXS_x = new TH1D("LOXS_x", "LOXS_x", binX, s_low, s_high);
TH1D* hXSR_x = new TH1D("XSR_x", "XSR_x", binX, s_low, s_high);
TH1D* hWXS_y = new TH1D("WXS_y", "WXS_y", binY, z_low, z_high);
TH1D* hLOXS_y = new TH1D("LOXS_y", "LOXS_y", binY, z_low, z_high);
TH1D* hXSR_y = new TH1D("XSR_y", "XSR_y", binY, z_low, z_high);
*/

/*
TH1F* hWXS_pt = new TH1F("WXS_pt", "WXS_pt", binX_pt, pt_low, pt_high);
TH1F* hLOXS_pt = new TH1F("LOXS_pt", "LOXS_pt", binX_pt, pt_low, pt_high);
TH1F* hXS_pt = new TH1F("XS_pt", "XS_pt", binX_pt, pt_low, pt_high);
TH1F* hXSR_pt = new TH1F("XSR_pt", "XSR_pt", binX_pt, pt_low, pt_high);
*/

TH2D* hdensity = new TH2D("density", "Jacobian matrix dz/d#Deltay; M_{t#bar{t}} (GeV); #Deltay_{t#bar{t}}", binX, s_low, s_high, binX_y, y_low, y_high);
TH2D* hWXS_gg = new TH2D("WXS_gg", "EW d#sigma(s, #Deltay)/d#Deltay (p^{-1}b); M_{t#bar{t}} (GeV); #Deltay_{t#bar{t}}", binX, s_low, s_high, binX_y, y_low, y_high);
TH2D* hLOXS_gg = new TH2D("LOXS_gg", "LO #Delta#sigma(s, #Deltay) (p^{-1}b); M_{t#bar{t}} (GeV); #Deltay_{t#bar{t}}", binX, s_low, s_high, binX_y, y_low, y_high);
TH2D* hXSR_gg = new TH2D("XSR_gg", "ratio matrix of EW/LO; M_{t#bar{t}} (GeV); #Deltay_{t#bar{t}}", binX, s_low, s_high, binX_y, y_low, y_high);

TH1D* hWXS_gg_x = new TH1D("WXS_gg_x", "EW ", binX, s_low, s_high);
TH1D* hWXS_gg_y = new TH1D("WXS_gg_y", "WXS_y", binX_y, y_low, y_high);
TH1D* hLOXS_gg_x = new TH1D("LOXS_gg_x", "LO ", binX, s_low, s_high);
TH1D* hLOXS_gg_y = new TH1D("LOXS_gg_y", "LOXS_y", binX_y, y_low, y_high);
TH1D* hXSR_gg_x = new TH1D("XSR_gg_x", "#frac{d#sigma_{EW}}{dM_{t#bar{t}}}/#frac{d#sigma_{LO}}{dM_{t#bar{t}}}; M_{t#bar{t}}", binX, s_low, s_high);
TH1D* hXSR_gg_y = new TH1D("XSR_gg_y", "#frac{d#sigma_{EW}}{dM_{t#bar{t}}}/#frac{d#sigma_{LO}}{dM_{t#bar{t}}}; #Deltay_{t#bar{t}}", binX_y, y_low, y_high);


TH2D* hWXSup_qq = new TH2D("WXSup_qq", "EW d#sigma(s, #Deltay)/d#Deltay (p^{-1}b); M_{t#bar{t}} (GeV); #Deltay_{t#bar{t}}", binX, s_low, s_high, binX_y, y_low, y_high);
TH2D* hWXSdown_qq = new TH2D("WXSdown_qq", "EW d#sigma(s, #Deltay)/d#Deltay (p^{-1}b); M_{t#bar{t}} (GeV); #Deltay_{t#bar{t}}", binX, s_low, s_high, binX_y, y_low, y_high);
TH2D* hLOXS_qq = new TH2D("LOXS_qq", "LO #Delta#sigma(s, #Deltay) (p^{-1}b); M_{t#bar{t}} (GeV); #Deltay_{t#bar{t}}", binX, s_low, s_high, binX_y, y_low, y_high);
TH2D* hXSRup_qq = new TH2D("XSRup_qq", "ratio matrix of EW/LO; M_{t#bar{t}} (GeV); #Deltay_{t#bar{t}}", binX, s_low, s_high, binX_y, y_low, y_high);
TH2D* hXSRdown_qq = new TH2D("XSRdown_qq", "ratio matrix of EW/LO; M_{t#bar{t}} (GeV); #Deltay_{t#bar{t}}", binX, s_low, s_high, binX_y, y_low, y_high);

TH1D* hWXSup_qq_x = new TH1D("WXSup_qq_x", "EW ", binX, s_low, s_high);
TH1D* hWXSup_qq_y = new TH1D("WXSup_qq_y", "WXS_y", binX_y, y_low, y_high);
TH1D* hWXSdown_qq_x = new TH1D("WXSdown_qq_x", "EW ", binX, s_low, s_high);
TH1D* hWXSdown_qq_y = new TH1D("WXSdown_qq_y", "WXS_y", binX_y, y_low, y_high);
TH1D* hXSRdown_qq_x = new TH1D("XSRdown_qq_x", "#frac{d#sigma_{EW}}{dM_{t#bar{t}}}/#frac{d#sigma_{LO}}{dM_{t#bar{t}}}; M_{t#bar{t}}", binX, s_low, s_high);
TH1D* hXSRdown_qq_y = new TH1D("XSRdown_qq_y", "#frac{d#sigma_{EW}}{dM_{t#bar{t}}}/#frac{d#sigma_{LO}}{dM_{t#bar{t}}}; #Deltay_{t#bar{t}}", binX_y, y_low, y_high);
TH1D* hLOXS_qq_x = new TH1D("LOXS_qq_x", "LO ", binX, s_low, s_high);
TH1D* hLOXS_qq_y = new TH1D("LOXS_qq_y", "LOXS_y", binX_y, y_low, y_high);
TH1D* hXSRup_qq_x = new TH1D("XSRup_qq_x", "#frac{d#sigma_{EW}}{dM_{t#bar{t}}}/#frac{d#sigma_{LO}}{dM_{t#bar{t}}}; M_{t#bar{t}}", binX, s_low, s_high);
TH1D* hXSRup_qq_y = new TH1D("XSRup_qq_y", "#frac{d#sigma_{EW}}{dM_{t#bar{t}}}/#frac{d#sigma_{LO}}{dM_{t#bar{t}}}; #Deltay_{t#bar{t}}", binX_y, y_low, y_high);


/*
TH2D* hdensity = new TH2D("density", "density_beta", binX_beta, beta_low, beta_high, 100, -1, 1);
TH2D* hLOXS = new TH2D("LOXS", "LOXS_beta_z", binX_beta, beta_low, beta_high, 100, -1, 1);
TH2D* hWXS = new TH2D("WXS", "WXS_beta_z", binX_beta, beta_low, beta_high, 100, -1, 1);
TH2D* hXSR = new TH2D("XSR", "XSR_beta_z", binX_beta, beta_low, beta_high, 100, -1, 1);
TH1D* hWXS_x = new TH1D("WXS_x", "WXS_x", binX_beta, beta_low, beta_high);
TH1D* hLOXS_x = new TH1D("LOXS_x", "LOXS_x", binX_beta, beta_low, beta_high);
TH1D* hXSR_x = new TH1D("XSR_x", "XSR_x", binX_beta, beta_low, beta_high);
TH1D* hWXS_y = new TH1D("WXS_y", "WXS_y", 100, -1, 1);
TH1D* hLOXS_y = new TH1D("LOXS_y", "LOXS_y", 100, -1, 1);
TH1D* hXSR_y = new TH1D("XSR_y", "XSR_y", 100, -1, 1);
*/

/*
void diffMtt(){

  Lhapdf pdf("NNPDF30_nlo_as_0118");
  Hathor weak(pdf);
  
  double ggW, qqW, ggW_test;
  double ggLO, qqLO, ggLO_test;
  double upflavor, downflavor;

  WeakCorrections weak2(Mb, MZ, MW, Mh, Alpha, Sin2, Hcq);
  weak2.setLambdat(gt);

  for(int z_bin = 0; z_bin < binY; z_bin++){
	for(int s_bin = 0; s_bin < binX; s_bin++){
		  double s = hLOXS -> GetXaxis() -> GetBinCenter(s_bin+1);
		  double z = hLOXS -> GetYaxis() -> GetBinCenter(z_bin+1);
	      //double s = s_bin* s_step + s_low;
		  //double z = z_bin* z_step + z_low;
		  ggW = weak2.dsigmaWeakgg(mt, s*s, z);
		  //weak2.dsigmaWeakqq(mt, s*s, z, upflavor, downflavor);
		  //cout<<upflavor<<" "<<downflavor<<endl;
		  ggLO = weak2.dsigmagg(mt, s*s, z);
		  //qqLO = weak2.dsigmaqq(mt, s*s, z);
	      //cout<< mt <<", "<< z <<", "<< s <<", "<< ggLO <<", "<< ggW << endl;
	      //cout<< mt <<", "<< z <<", "<< s <<", "<< qqLO <<", "<< qqW << endl;

		  hWXS->SetBinContent(s_bin+1, z_bin+1, ggW);
		  hLOXS->SetBinContent(s_bin+1, z_bin+1, ggLO);
		  hXSR->SetBinContent(s_bin+1, z_bin+1, ggW/ggLO);
	  }
  }

  hWXS_x = hWXS->ProjectionX();
  hWXS_y = hWXS->ProjectionY();
  hWXS_x ->Scale(z_step); 
  hWXS_y ->Scale(s_step);

  hLOXS_x = hLOXS->ProjectionX();
  hLOXS_y = hLOXS->ProjectionY();
  hLOXS_x ->Scale(z_step);
  hLOXS_y ->Scale(s_step);
  
  hXSR_x ->Divide(hWXS_x, hLOXS_x, 1, 1);
  hXSR_y ->Divide(hWXS_y, hLOXS_y, 1, 1);


  hWXS->Draw();
  hWXS_x->Draw();
  hWXS_y->Draw();
  hLOXS->Draw();
  hLOXS_x->Draw();
  hLOXS_y->Draw();
  hXSR->Draw();
  hXSR_x->Draw();
  hXSR_y->Draw();

  hWXS->Write();
  hWXS_x->Write();
  hWXS_y->Write();
  hLOXS->Write();
  hLOXS_x->Write();
  hLOXS_y->Write();
  hXSR->Write();
  hXSR_x->Write();
  hXSR_y->Write();

  f.Close();
}
*/
/*
void diffPt(){

  Lhapdf pdf("NNPDF30_nlo_as_0118");
  Hathor weak(pdf);
  
  double ggW, ggW_test;
  double ggLO, qqLO, ggLO_test;
  double upflavor, downflavor;

  WeakCorrections weak2(Mb, MZ, MW, Mh, Alpha, Sin2, Hcq);
  weak2.setLambdat(gt);

  for(int z_bin = 0; z_bin < binY; z_bin++){
	for(int pt_bin = 0; pt_bin < binX_pt; pt_bin++){
		  double z = z_bin* z_step + z_low;//[-1,1]
		  double s = 2*std::sqrt(std::pow(pt_bin*pt_step+pt_low, 2)/(1-z*z)+mt*mt);
		  
		  ggW = weak2.dsigmaWeakgg(mt, s*s, z);
		  weak2.dsigmaWeakqq(mt, s*s, z, upflavor, downflavor);
		  ggLO = weak2.dsigmagg(mt, s*s, z);
		  qqLO = weak2.dsigmaqq(mt, s*s, z);

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
*/

void diffDeltaY(){

  Lhapdf pdf("NNPDF30_nlo_as_0118");
  Hathor weak(pdf);
  
  double ggW, ggW_test;
  double ggLO, qqLO, ggLO_test;
  double upflavor, downflavor;

  WeakCorrections weak2(Mb, MZ, MW, Mh, Alpha, Sin2, Hcq);
  weak2.setLambdat(gt);

  for(int y_bin = 0; y_bin < binX_y; y_bin++){
	for(int s_bin = 0; s_bin < binX; s_bin++){
		  double s = hLOXS_gg -> GetXaxis() -> GetBinCenter(s_bin+1);
		  double y = hLOXS_gg -> GetYaxis() -> GetBinCenter(y_bin+1);
		  //double s = s_bin*s_step + s_low;
		  //double z = 1/std::sqrt(1-4*(mt*mt)/(s*s)) * (std::exp(y_bin*y_step + y_low)-1)/(std::exp(y_bin*y_step + y_low)+1);
		  double z = 1/std::sqrt(1-4*(mt*mt)/(s*s)) * (std::exp(y)-1)/(std::exp(y)+1);
		  //cout<< y_bin*y_step+y_low <<endl;
		  double density = s*std::sqrt(s*s/4-mt*mt) / ((1-z*z)*s*s/4 + mt*mt*z*z);
		  density = 1/density;
		  if(std::abs(z)>1){
			  cout<<" >>> "<< s <<", "<< z <<", "<<density<<endl;
			  continue;
		  }
		  

		  ggW = weak2.dsigmaWeakgg(mt, s*s, z);
		  ggLO = weak2.dsigmagg(mt, s*s, z);
		  qqLO = weak2.dsigmaqq(mt, s*s, z);
		  weak2.dsigmaWeakqq(mt, s*s, z, upflavor, downflavor);
		  //qqLO = weak2.dsigmaqq(mt, s*s, z);
		  //cout << s <<", "<< z <<", "<< ggLO <<", "<< ggW <<", "<<density<<endl;

		  hdensity->SetBinContent(s_bin+1, y_bin+1, density);

		  hLOXS_gg->SetBinContent(s_bin+1, y_bin+1, ggLO*density);
		  hWXS_gg->SetBinContent(s_bin+1, y_bin+1, ggW*density);
		  hXSR_gg->SetBinContent(s_bin+1, y_bin+1, ggW/ggLO);
		  
		  hLOXS_qq->SetBinContent(s_bin+1, y_bin+1, qqLO*density);
		  hWXSup_qq->SetBinContent(s_bin+1, y_bin+1, upflavor*density);
		  hWXSdown_qq->SetBinContent(s_bin+1, y_bin+1, downflavor*density);
		  hXSRup_qq->SetBinContent(s_bin+1, y_bin+1, upflavor/qqLO);
		  hXSRdown_qq->SetBinContent(s_bin+1, y_bin+1, downflavor/qqLO);
	  

	  }
  }

  // the gg
  hWXS_gg_y = hWXS_gg->ProjectionY();
  hWXS_gg_x = hWXS_gg->ProjectionX();
  hWXS_gg_y ->Scale(s_step); 
  hWXS_gg_x ->Scale(y_step);

  hLOXS_gg_y = hLOXS_gg->ProjectionY();
  hLOXS_gg_x = hLOXS_gg->ProjectionX();
  hLOXS_gg_y ->Scale(s_step);
  hLOXS_gg_x ->Scale(y_step);

  hXSR_gg_y ->Divide(hWXS_gg_y, hLOXS_gg_y, 1, 1);
  hXSR_gg_x ->Divide(hWXS_gg_x, hLOXS_gg_x, 1, 1);
	
  hWXS_gg->Draw();
  hWXS_gg_x->Draw();
  hWXS_gg_y->Draw();
  hLOXS_gg->Draw();
  hLOXS_gg_x->Draw();
  hLOXS_gg_y->Draw();
  hXSR_gg->Draw();
  hXSR_gg_x->Draw();
  hXSR_gg_y->Draw();

  hWXS_gg->Write();
  hWXS_gg_x->Write();
  hWXS_gg_y->Write();
  hLOXS_gg->Write();
  hLOXS_gg_x->Write();
  hLOXS_gg_y->Write();
  hXSR_gg->Write();
  hXSR_gg_x->Write();
  hXSR_gg_y->Write();
  
  // the qq 
  hWXSup_qq_y = hWXSup_qq->ProjectionY();
  hWXSup_qq_x = hWXSup_qq->ProjectionX();
  hWXSup_qq_y ->Scale(s_step); 
  hWXSup_qq_x ->Scale(y_step);  
  hWXSdown_qq_y = hWXSdown_qq->ProjectionY();
  hWXSdown_qq_x = hWXSdown_qq->ProjectionX();
  hWXSdown_qq_y ->Scale(s_step); 
  hWXSdown_qq_x ->Scale(y_step);


  hLOXS_qq_y = hLOXS_qq->ProjectionY();
  hLOXS_qq_x = hLOXS_qq->ProjectionX();
  hLOXS_qq_y ->Scale(s_step);
  hLOXS_qq_x ->Scale(y_step);

  hXSRup_qq_y ->Divide(hWXSup_qq_y, hLOXS_qq_y, 1, 1);
  hXSRup_qq_x ->Divide(hWXSup_qq_x, hLOXS_qq_x, 1, 1);
  hXSRdown_qq_y ->Divide(hWXSdown_qq_y, hLOXS_qq_y, 1, 1);
  hXSRdown_qq_x ->Divide(hWXSdown_qq_x, hLOXS_qq_x, 1, 1);
	
  hWXSup_qq->Draw();
  hWXSdown_qq->Draw();
  hWXSup_qq_x->Draw();
  hWXSup_qq_y->Draw();
  hWXSdown_qq_x->Draw();
  hWXSdown_qq_y->Draw();
  hLOXS_qq->Draw();
  hLOXS_qq_x->Draw();
  hLOXS_qq_y->Draw();
  hXSRup_qq->Draw();
  hXSRdown_qq->Draw();
  hXSRup_qq_x->Draw();
  hXSRup_qq_y->Draw();
  hXSRdown_qq_x->Draw();
  hXSRdown_qq_y->Draw();

  hWXSup_qq->Write();
  hWXSdown_qq->Write();
  hWXSup_qq_x->Write();
  hWXSup_qq_y->Write();
  hWXSdown_qq_x->Write();
  hWXSdown_qq_y->Write();
  hLOXS_qq->Write();
  hLOXS_qq_x->Write();
  hLOXS_qq_y->Write();
  hXSRup_qq->Write();
  hXSRdown_qq->Write();
  hXSRup_qq_x->Write();
  hXSRup_qq_y->Write();
  hXSRdown_qq_x->Write();
  hXSRdown_qq_y->Write();

  hdensity->Draw();
  hdensity->Write();

  f.Close();
}

/*
void diffBeta(){

  Lhapdf pdf("NNPDF30_nlo_as_0118");
  Hathor weak(pdf);
  
  double ggW, ggW_test;
  double ggLO, qqLO, ggLO_test;
  double upflavor, downflavor;

  WeakCorrections weak2(Mb, MZ, MW, Mh, Alpha, Sin2, Hcq);
  weak2.setLambdat(gt);

  for(int z_bin = 0; z_bin < binY; z_bin++){
	for(int beta_bin = 0; beta_bin < binX_beta; beta_bin++){
		  //double z = z_bin* z_step + z_low;//[-1,1]
		  //double beta = beta_bin* beta_step + beta_low;
		  double z = hLOXS -> GetYaxis() -> GetBinCenter(z_bin+1);
		  double beta = hLOXS -> GetXaxis() -> GetBinCenter(beta_bin+1);
		  //cout<< z_bin <<",  "<< beta_bin <<endl;
		  //double s = 2*mt* std::sqrt( 1/(1-std::pow(beta/4,2)) );
		  double s = 2*mt/std::sqrt(1-beta*beta/4);
		  double density = 8*mt*mt/( std::sqrt(1-4*mt*mt/(s*s)) * std::pow(s,3));
		  density = 1/density;
		  //double density = 1;
		  //double s = std::sqrt(4-beta*beta)/(mt*mt);
		  //if(s<350*350) { 
			//  cout <<">>>> "<< z <<", "<< std::sqrt(s) <<", "<< ggLO <<endl;
			  //continue;}

		  ggW = weak2.dsigmaWeakgg(mt, s*s, z);
		  //weak2.dsigmaWeakqq(mt, s*s, z, upflavor, downflavor);
		  ggLO = weak2.dsigmagg(mt, s*s, z);
		  //qqLO = weak2.dsigmaqq(mt, s*s, z);
	      cout<< mt <<", "<< z <<", "<< beta <<", "<< s <<", "<< ggLO <<", "<< ggW << endl;

		  hLOXS->SetBinContent(beta_bin+1, z_bin+1, ggLO*density);
		  hWXS->SetBinContent(beta_bin+1, z_bin+1, ggW*density);
		  hXSR->SetBinContent(beta_bin+1, z_bin+1, ggW/ggLO);
		  hdensity->SetBinContent(beta_bin+1, z_bin+1, density);
	  }
  }

  hWXS_x = hWXS->ProjectionX();
  hWXS_y = hWXS->ProjectionY();
  hWXS_x ->Scale(z_step); 
  hWXS_y ->Scale(beta_step); 

  hLOXS_x = hLOXS->ProjectionX();
  hLOXS_y = hLOXS->ProjectionY();
  hLOXS_x ->Scale(z_step); 
  hLOXS_y ->Scale(beta_step);

  hXSR_x ->Divide(hWXS_x, hLOXS_x, 1, 1);
  hXSR_y ->Divide(hWXS_y, hLOXS_y, 1, 1);

  hWXS->Draw();
  hWXS_x->Draw();
  hWXS_y->Draw();
  hLOXS->Draw();
  hLOXS_x->Draw();
  hLOXS_y->Draw();
  hXSR->Draw();
  hXSR_x->Draw();
  hXSR_y->Draw();

  hWXS->Write();
  hWXS_x->Write();
  hWXS_y->Write();
  hLOXS->Write();
  hLOXS_x->Write();
  hLOXS_y->Write();
  hXSR->Write();
  hXSR_x->Write();
  hXSR_y->Write();

  hdensity->Draw();
  hdensity->Write();

  f.Close();
}
*/

int main(){

  //diffMtt();
  //diffPt();
  diffDeltaY();
  //diffBeta();

}
