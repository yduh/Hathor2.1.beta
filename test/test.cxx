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


using namespace std;

int binY = 18;
int binX = 8;
double s_step = 50;
double z_step = 0.1;

TFile f("test.root", "recreate");

TH2F* hWXS = new TH2F("WXS", "WXS", binX, 400, 800, binY, -0.9, 0.9);
TH2F* hLOXS = new TH2F("LOXS", "LOXS", binX, 400, 800, binY, -0.9, 0.9);
TH2F* hXS = new TH2F("XS", "XS", binX, 400, 800, binY, -0.9, 0.9);

TH1F* hWXS_X = new TH1F("WXS_X", "WXS_X", binX, 400, 800);
TH1F* hLOXS_X = new TH1F("LOXS_X", "LOXS_X", binX, 400, 800);
TH1F* hXS_X = new TH1F("XS_X", "XS_X", binX, 400, 800);

TH1F* hWXS_Y = new TH1F("WXS_Y", "WXS_Y", binY, -0.9, 0.9);
TH1F* hLOXS_Y = new TH1F("LOXS_Y", "LOXS_Y", binY, -0.9, 0.9);
TH1F* hXS_Y = new TH1F("XS_Y", "XS_Y", binY, -0.9, 0.9);

void nonSM(){

  Lhapdf pdf("NNPDF30_nlo_as_0118");
  Hathor weak(pdf);
  
//  double factor = 100.;
//  weak.setLambdat(factor);

  //cout << "yukuwa const " << yukuwaC << endl;
   double mt=173;//, s=400, z=0.5;
  //double mt=500.;
  //double val,err,chi;
  double ggW;
  double ggLO, qqLO;
  double upflavor, downflavor;

 // double ecms[1] = {13000.};
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
  double Mb = 4.82, MZ = 91.1876, MW = 80.385, Mh = 126. ;
  double Alpha = 1/126.3, Sin2 = 1-(MW/MZ)*(MW/MZ), Hcq = 0.389379323e9;
  WeakCorrections weak2(Mb, MZ, MW, Mh, Alpha, Sin2, Hcq);

  double factor = 1.;
  weak2.setLambdat(factor);

 /* 
  int Zstep = 10;
  for(float i=0; i<=Zstep; i++){
	gg = weak2.dsigmaWeakgg(mt, s*s, i/Zstep);
	cout << mt <<", "<< s <<", "<< i/Zstep <<", "<< gg << endl;
	//weak.printParameters();
  }

  //int Sstep = 30;
  for(float i=400; i<=800; i=i+50){
	gg = weak2.dsigmaWeakgg(mt, i*i, z);
	cout << mt <<", "<< i*i <<", "<< z <<", "<< gg << endl;
  }
*/

  //for(double z=-0.9; z<0.9; z=z+0.1){
	//for(double s=400; s<=800; s=s+50){
	for(int z_bin = 0; z_bin < binY; z_bin++){
	  for(int s_bin = 0; s_bin < binX; s_bin++){
		  double z = z_bin* z_step - 0.9;
		  double s = s_bin* s_step + 400;

		  ggW = weak2.dsigmaWeakgg(mt, s*s, z);
		  weak2.dsigmaWeakqq(mt, s*s, z, upflavor, downflavor);
		  cout << mt <<", "<< s <<", "<< z <<", "<< ggW << endl;
		  hWXS->SetBinContent(s_bin+1, z_bin+1, ggW);
		  hWXS_X->SetBinContent(s_bin+1, ggW);
		  hWXS_Y->SetBinContent(z_bin+1, ggW);
		  
		  ggLO = weak2.dsigmagg(mt, s*s, z);
		  qqLO = weak2.dsigmaqq(mt, s*s, z);
		  cout << mt <<", "<< s <<", "<< z <<", "<< ggLO <<" "<< qqLO << endl;
		  hLOXS->SetBinContent(s_bin+1, z_bin+1, ggLO);
		  hLOXS_X->SetBinContent(s_bin+1, ggLO);
		  hLOXS_Y->SetBinContent(z_bin+1, ggLO);


		  hXS->SetBinContent(s_bin+1, z_bin+1, ggLO-ggW);
		  hXS_X->SetBinContent(s_bin+1, ggLO-ggW);
		  hXS_Y->SetBinContent(z_bin+1, ggLO-ggW);
		  
	  }
  }

  hWXS->Draw();
  hWXS_X->Draw();
  hWXS_Y->Draw();
  hLOXS->Draw();
  hLOXS_X->Draw();
  hLOXS_Y->Draw();
  hXS->Draw();
  hXS_X->Draw();
  hXS_Y->Draw();

  hWXS->Write();
  hWXS_X->Write();
  hWXS_Y->Write();
  hLOXS->Write();
  hLOXS_X->Write();
  hLOXS_Y->Write();
  hXS->Write();
  hXS_X->Write();
  hXS_Y->Write();

  f.Close();
}
int main(){

  nonSM();
}
