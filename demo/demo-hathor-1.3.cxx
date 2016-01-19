// $Modified: Thu Apr 26 15:27:28 2012 by uwer $
#include "Hathor.h"
#include "HathorPdf.h"
#include <iostream>
#include <fstream>


using namespace std;



void MochUwerVogt2012(int ipdf, unsigned int what, double mt){


  string pdfname[]={"abm11_5n_nnlo","MSTW2008nnlo68cl"};

  Lhapdf lhapdf(pdfname[ipdf]);
  Hathor XS(lhapdf);
  XS.setNf(5);

  unsigned int scheme = Hathor::LO  | Hathor::NLO | Hathor::NNLO | what;

  scheme |= Hathor::PDF_SCAN 
    | Hathor::NNLOAPPROX  
    | Hathor::NOQQ | Hathor::NOQQP| Hathor::NOQQPB;

  XS.setScheme(scheme);
  XS.setPrecision(Hathor::HIGH);

  Hathor::COLLIDERTYPE collider[]={Hathor::PPBAR, 
				   Hathor::PP,Hathor::PP,Hathor::PP };

  double ecms[] = {1960.,7000.,8000.,14000};

  
  double val,err,chi,up,down;

  for(int i=0; i<4;i++){
    XS.setColliderType(collider[i]);
    XS.setSqrtShad(ecms[i]);
    XS.getXsection(mt,mt,mt);
    XS.getResult(0,val,err,chi);
    cout <<ecms[i]<< " " << val << " " << err << endl;
    XS.getPdfErr(up,down);
    cout << "PDF uncertainty: +" << up << " -" << down << endl; 
  }
}

void BaernreutherCzakonMitov2012(unsigned int what){

  double val,err,chi,up,down;
  double mt = 173.3;
  string pdfname="MSTW2008nnlo68cl";

  Lhapdf lhapdf(pdfname);
  Hathor XS(lhapdf);
  XS.setNf(5);

  unsigned int scheme = Hathor::LO  | Hathor::NLO | Hathor::NNLO 
    | what;
  XS.setScheme(scheme);
  XS.setColliderType(Hathor::PPBAR);
  XS.setSqrtShad(1960.);
  XS.setPrecision(Hathor::MEDIUM);
  XS.getXsection(mt,mt,mt);
  XS.getResult(0,val,err,chi);
  cout  << val << " " << err << endl;
  
}

int main(int argc, char* argv[]){

  cout << "---------------------------------------------------------\n";
  cout << "-- Results Moch, Uwer, Vogt arXiv:1203.6282  ------------\n";
  cout << "-- to appear in Physics Letters B -----------------------\n";
  cout << "---------------------------------------------------------\n";
  cout << "------------- ABM - POLE MASS 173 -----------------------\n"; 
  MochUwerVogt2012(0,Hathor::NOHIGHENERGY,173.);
  cout << "------------- ABM - POLE MASS 173, HIGH ENERGY ----------\n"; 
  MochUwerVogt2012(0,0,173.);
  cout << "------------- ABM - MS MASS 164 -------------------------\n"; 
  MochUwerVogt2012(0,Hathor::MS_MASS|Hathor::NOHIGHENERGY,164.);
  cout << "------------- ABM - MS MASS 164, HIGH ENERGY ------------\n"; 
  MochUwerVogt2012(0,Hathor::MS_MASS,164.);

  cout << "------------- MSTW - POLE MASS 173 ----------------------\n"; 
  MochUwerVogt2012(1,Hathor::NOHIGHENERGY,173.);  
  cout << "------------- MSTW - POLE MASS 173, HIGH ENERGY ---------\n"; 
  MochUwerVogt2012(1,0,173.);

  cout << "---------------------------------------------------------\n";
  cout << "-- Comparison with Baernreuther, Czakon, Mitov ----------\n";
  cout << "-- arXiv:1206.0621 --------------------------------------\n";
  cout << "---------------------------------------------------------\n";
  cout << "NNLO exact for qq, NOHIGHENERGY" << endl;
  BaernreutherCzakonMitov2012(Hathor::NOHIGHENERGY);
  cout << "NNLO exact for qq, + HIGHENERGY approximation for qg, gg" << endl;
  BaernreutherCzakonMitov2012(0);
  cout << "NNLOAPPROX for qq, NOHIGHENERGY (old results)" << endl;
  BaernreutherCzakonMitov2012(Hathor::NNLOAPPROX | Hathor::NOHIGHENERGY);
  cout << "NNLOAPPROX for qq,+ HIGHENERGY approximation" << endl;
  BaernreutherCzakonMitov2012(Hathor::NNLOAPPROX );

  return(0);
}
