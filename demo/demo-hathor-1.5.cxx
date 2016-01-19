// $Modified: Thu Jun 24 09:56:51 2010 by uwer $
#include "Hathor.h"
#include "HathorPdf.h"
#include <iostream>
#include <fstream>


using namespace std;


void CFM13(){
  double mt=173.3; 
  double val,err,chi;

  double ecms[4]={1960.,7000.,8000.,14000.};
  Hathor::COLLIDERTYPE ctype[]={Hathor::PPBAR,Hathor::PP
				,Hathor::PP,Hathor::PP};
  Lhapdf lhapdf("MSTW2008nnlo68cl");
  Hathor XS(lhapdf);

  cout << "Reproduce fixed order results from arXiv:1303.6254\n";
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

int main(){
  CFM13();
}
