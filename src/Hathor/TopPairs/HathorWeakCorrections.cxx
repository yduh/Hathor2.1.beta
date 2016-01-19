// $Modified: Mon Feb  3 15:17:55 2014 by puwer $
#include <cmath>
#include <iostream>
#include <string> 
#include <sstream> 
#include <stdlib.h>
#include "AbstractHathor.h"
#include "HathorWeakCorrections.h"
#include "FF.h"

using namespace std;

#define INFO(_X_) AbstractHathor::info(_X_)
#define PRINT(_X_) cout <<"# "<< #_X_ <<" = " << _X_ << endl;

class CheckData {
public:
  double ecms,z,res;
};

class InitLoopIntegrals {
public:
  InitLoopIntegrals(){
#ifdef QCDLOOP
    qlinit_();
    INFO("FF [1] and QCDLoop [2] are used to calculate the");
    INFO("scalar one-loop integrals");
    INFO("[1] van Oldenborgh,");
    INFO("    \"FF: A Package To Evaluate One Loop Feynman Diagrams\",");
    INFO("    Comput.Phys.Commun.66:1-15,1991");
    INFO("[2] R.Keith Ellis, Giulia Zanderighi,");
    INFO("    \"Scalar one-loop integrals for QCD\", ");
    INFO("    JHEP 0802:002,2008");
    INFO("Please cite this properly when you use this program");
#else
    INFO("FF [1] is used to calculate the scalar one-loop integrals");
    INFO("[1] van Oldenborgh,"); 
    INFO("    \"FF: A Package To Evaluate One Loop Feynman Diagrams\",");
    INFO("    Comput.Phys.Commun.66:1-15,1991");
    INFO("Please cite this properly when you use this program");
    ffini_();
    ffcut_.delta = 0.0;
#endif
  }
};

static InitLoopIntegrals dummy;

WeakCorrections::WeakCorrections(double & mb_, double & mz_,double & mw_,
				 double & mh_,
				 double & alpha_,double & swq_,double & hcq_) :
  mb(mb_), mz(mz_), mw(mw_),mh(mh_), alpha(alpha_),swq(swq_), hcq(hcq_){

  /*
   * Set default values:
   */
  N = 3.0; alphas = 1.;
  ggtriangle = ggself = ggbox = ggvertex = 1;
  GammaH = 0; GammaHq = GammaH*GammaH;
  setScale(170.);
  lambdat = 1; lambdat2 = lambdat*lambdat;
  updateCouplings();

}

void WeakCorrections::setSwq(double value){
  swq = value;
  updateCouplings();
}

void WeakCorrections::setLambdat(double value){
  lambdat = value; lambdat2 = lambdat*lambdat;
  AbstractHathor::info("Using non-SM value for top-quark Yukawa coupling lambda_t = ", lambdat );
}

void WeakCorrections::check(void){
  /* 
   * To make sure that the formulae are correct we compare 
   * with old results produced with the code used for
   * the respective publications.
   */
  
  CheckData gg[] = {{400.,0.4711,-0.0901199730013604},
		    {1000.,0.8,-0.216713029203266},
		    {1400,0.2,-0.040614655540509}};

  CheckData qqup[] = {{400.,0.4711,0.0307124934100407},
		    {1000.,0.8,-0.210832061451692},
		    {1400,0.2,-0.0949539670877666}};

  CheckData qqdown[] = {{400.,0.4711,0.0295240331116438},
		    {1000.,0.8,-0.214186485748532},
		    {1400,0.2,-0.0967086144237238}};

  double alphas_data = 0.106823089396409;
  double alpha_data = 1./126.3;
  double hcq_data = 0.389379323e9;
  double mw_data = 80.385, mt_data = 173.2, mb_data = 4.82,
    mz_data = 91.1876, mh_data = 126.;
  
  double hcq_tmp = hcq, alpha_tmp=alpha, swq_tmp=swq,
    mw_tmp = mw, mb_tmp=mb, mz_tmp=mz, mh_tmp=mh,
    lambdat_tmp=lambdat;
  
  // Set all parameters to the values used to produce the data above:
  mb = mb_data;
  mz = mz_data;
  mw = mw_data;
  mh = mh_data;
  hcq = hcq_data;
  swq = 1.0 - (mw*mw)/(mz*mz);
  alpha = alpha_data;
  lambdat = 1; lambdat2 = lambdat*lambdat;
  updateCouplings();

  stringstream s;
  s.precision(15);

  INFO("Check implementation of weak corrections,");
  INFO("compare with reference values...");
  INFO("------------------ Gluon channel -------------------");
  for(unsigned int i=0;i<sizeof(gg)/3/sizeof(double);i++){
    double res = dsigmaWeakgg(mt_data,gg[i].ecms*gg[i].ecms, gg[i].z); 
    s.str("");
    s << "Ratio this_prg/ref: "
		   << alphas_data*alphas_data*res/gg[i].res;
    INFO(s.str());
  }

  INFO("----------------- Quark channel --------------------");
  double up,down;
  for(unsigned int i=0;i<sizeof(gg)/3/sizeof(double);i++){
    dsigmaWeakqq(mt_data,gg[i].ecms*gg[i].ecms, gg[i].z,up,down); 
    s.str("");
    s << "Ratio this_prg/ref: "
		   << alphas_data*alphas_data*up/qqup[i].res << " " 
		   << alphas_data*alphas_data*down/qqdown[i].res;
    INFO(s.str());
  }
  // Restore the parameter to saved values:
  mb = mb_tmp;
  mz = mz_tmp;
  mw = mw_tmp;
  mh = mh_tmp;
  hcq = hcq_tmp;
  alpha = alpha_tmp;
  swq = swq_tmp;
  lambdat = lambdat_tmp; lambdat2 = lambdat*lambdat;
  updateCouplings();
}


void WeakCorrections::updateCouplings(void){
  cwq = 1.-swq;
  gat = 0.5 / sqrt(swq*cwq) * (0.5);
  gab = 0.5 / sqrt(swq*cwq) * (-0.5);
  gvt = 0.5 / sqrt(swq*cwq ) * ( + 0.5 - 2.0 * swq * 2./3. );  
  gvb = 0.5 / sqrt(swq*cwq ) * ( - 0.5 - 2.0 * swq * (-1./3.) ); 
  gvt2 = gvt*gvt; gat2=gat*gat;
  gvb2 = gvb*gvb; gab2=gab*gab;
  gw2 = 1./8./swq;
  gw  = sqrt(gw2);
}

void WeakCorrections::printParameters(){
  cout << "# The following parameters are currently used:" << endl;
  PRINT(hcq);
  PRINT(alphas);
  PRINT(1/alpha);
  PRINT(swq);
  PRINT(lambdat);
  PRINT(mt);
  PRINT(mb);
  PRINT(mw);
  PRINT(mz);
}

double WeakCorrections::dsigmagg(const double mt_, const double sparton,
				 const double z){

  mt = mt_; mtq = mt*mt;

  double z2=z*z,z4=z2*z2;
  double beta2 = 1.-4*mtq/sparton, beta4=beta2*beta2, beta=sqrt(beta2);
  double kappa = M_PI*M_PI*alphas*alphas/(1.-z2*beta2)/(1.-z2*beta2)*
    (N*N-2+N*N*z2*beta2)/N/(N*N-1.);

  return(beta / 4. / M_PI / sparton * kappa * 
	 ( 1.+2*beta2-2*z2*beta2-2*beta4+2*z2*beta4-z4*beta4) * hcq );
}

double WeakCorrections::dsigmaqq(const double mt_, const double sparton,
				 const double z){
  mt = mt_; mtq = mt*mt;
  double beta2 = 1.-4.*mtq/sparton;
  double beta = sqrt(beta2);
  return( 1./8. * (N*N-1.)/N/N * M_PI * alphas * alphas
	  * beta / sparton * ( 2 + ( z*z - 1 ) *beta2 ) *hcq );
}
				  
double WeakCorrections::ReC0m(const double shat, const double mq){
  const double betam = sqrt(1.0-4.0*mq/shat);
  return(1.0/2.0/shat 
	 * ( pow(log( ( 1.0 + betam ) / ( 1.0 - betam ) ),2) 
	     - pow(M_PI,2) ) );
}
 
double WeakCorrections::ImC0m(const double shat, const double mq){
  const double betam = sqrt(1.0-4.0*mq/shat);
  return(- M_PI/shat * log( ( 1.0 + betam ) / ( 1.0 - betam ) ) );
} 

double WeakCorrections::ggtriangles(const double shat, double z){
  return( Higgs_s_channel(shat,z) + Z_Chi_s_channel(shat,z) );
}

double WeakCorrections::Higgs_s_channel(const double shat, double z){
  const double mbq=mb*mb,mhq=mh*mh,mwq=mw*mw;
  const double beta2 = 1.0-4.0*mtq/shat;
  const double beta = sqrt(beta2);

  double C0mtRealPart = ReC0m(shat,mtq); 
  double C0mbRealPart = ReC0m(shat,mbq);
  double C0mtImPart = ImC0m(shat,mtq); 
  double C0mbImPart = ImC0m(shat,mbq);

  if ( mhq > 4. * mtq ) {
    if ( mhq == 1000.*1000. ){
      GammaH = 492.471324;
      GammaHq = GammaH * GammaH;
    } else {
      cout << "Required Higgs width not known" << endl;
      exit(1);
    } 
  }
 
  double sigma0 =  M_PI * alphas * alphas 
    / 4.0 / (N*N-1.0) * beta/shat; 

  return( sigma0 * alpha / M_PI 
	  * mtq / ( swq * mwq ) 
	   / ( (shat-mhq)*(shat-mhq) + mhq*GammaHq )
	  * beta2 / (1.0 - beta2*z*z)
	  * ( 
	     ( + mtq * ( shat - 4.0 * mtq ) * C0mtRealPart * lambdat2
	       + mbq * ( shat - 4.0 * mbq ) * C0mbRealPart * lambdat
	       - 2.0 * ( mtq * lambdat2 + mbq * lambdat )
	       ) * (shat-mhq) 
	     +
	     ( + mtq * ( shat - 4.0 * mtq ) * C0mtImPart * lambdat2
	       + mbq * ( shat - 4.0 * mbq ) * C0mbImPart * lambdat
	       ) * mh * GammaH
	     )
	  );

}

double WeakCorrections::Z_Chi_s_channel(const double shat, double z){

  const double mbq = mb*mb, mzq = mz*mz;
  const double beta2 = 1.0-4.0*mtq/shat;
  const double beta = sqrt(beta2);

  double C0mt = ReC0m(shat,mtq); 
  double C0mb = ReC0m(shat,mbq);
    
  /*
   * Note that in the following we used the on-shell relation between
   * mwq, mzq and swq. Using the MSbar value for swq leads than to an
   * 1 % effect with respect to the version where the on-shell relation
   * is not used.
   */
  double sigma0 =  M_PI * alphas * alphas 
    / 4.0 / (N*N-1.0) * beta/shat; 

  return(  16.0 * sigma0 * alpha / M_PI 
	   * mtq / (  mzq * ( 1.0 - beta2 * z * z ) ) 
	   * ( gat*gat*mtq*C0mt + gat*gab*mbq*C0mb ) );

}

double WeakCorrections::ggselfenergies(const double shat, double z){

  const double mwq=mw*mw, mzq=mz*mz,mhq=mh*mh, mbq=mb*mb;
  double selfenergy_Z = 0, selfenergy_H = 0, selfenergy_W = 0,
    selfenergy_Phi = 0, selfenergy_Chi = 0;

  double t;
  const double s = shat;
  const double N2 = N*N;

  double beta2 = 1.0-4.0*mtq/shat;
  double beta = sqrt(beta2);
  double beta3 = beta2*beta;
  double beta4 = beta2*beta2;
  double beta5 = beta4*beta;
  double beta6 = beta4*beta2;
  double beta7 = beta6*beta;
  double beta8 = beta6*beta2;
#include "auto/self.dec"

  double z2 = z*z;
  double z3 = z2 *z;
  double z4 = z2 * z2;
  double z5  = z4*z;
  double z6  = z4*z2;

  t = mtq - s/2.0 * ( 1.0 - beta*z );

  double sigma0 =  M_PI * alphas * alphas 
    / 4.0 / (N*N-1.0) * beta/shat; 

  for(int i = 0; i < 2; i++ ) { 
#include "auto/self.cpp"
    z = -z;
    z3 = z2 *z;
    z5  = z4*z;
    t = mtq - s/2.0*(1.0-beta*z);
  }
  return( 
	 + selfenergy_H * lambdat2  
	 + selfenergy_Z 
	 + selfenergy_Phi 
	 + selfenergy_W
	 + selfenergy_Chi 
	 );
}

double WeakCorrections::ggvertices(const double shat, double z){

  const double mzq=mz*mz,mwq=mw*mw,mhq=mh*mh,mbq=mb*mb;
  double Z_vertex, H_vertex, W_vertex,Phi_vertex,Chi_vertex;
  double Z_svertex, H_svertex, W_svertex,Phi_svertex,Chi_svertex;
  
  double vertices = 0.0, svertices =0.0;
  double t;
  const double s = shat;

  const double N2 = N*N;

  const double beta2 = 1.0-4.0*mtq/shat;
  const double beta = sqrt(beta2);


  const double beta3 = beta2*beta;
  const double beta4 = beta2*beta2;
  const double beta5 = beta4*beta;
  const double beta6 = beta4*beta2;
  const double beta7 = beta6*beta;
  const double beta8 = beta6*beta2;
  const double z2 = z*z;
  double z3 = z2 *z;
  const double z4 = z2 * z2;
  double z5  = z4*z;
  const double z6  = z4*z2;

#include "auto/vertices.dec"

  t = mtq - s/2.0*(1.0-beta*z);

  double sigma0 =  M_PI * alphas * alphas 
    / 4.0 / (N*N-1.0) * beta/shat; 

  for(int i = 0; i < 2; i++ ) {
#include "auto/vertices.cpp"
    vertices += 
      + Z_vertex 
      + H_vertex * lambdat2
      + W_vertex 
      + Phi_vertex 
      + Chi_vertex
      ;
    // The s-vertex is taken only once !!!
    svertices = 
      + Z_svertex 
      + H_svertex * lambdat2
      + W_svertex 
      + Phi_svertex 
      + Chi_svertex
      ;

    z  = -z;
    z3 = -z3;
    z5 = -z5;
    t = mtq - s/2.0*(1.0-beta*z);
  }
  
  return( vertices + svertices );
}

double WeakCorrections::ggboxes(const double shat, double z){

  const double mzq=mz*mz,mwq=mw*mw,mhq=mh*mh,mbq=mb*mb;
  double Z_box, H_box, W_box,Phi_box,Chi_box;
  
  double boxes_sum = 0.0;
  double t;
  const double s = shat;

  const double beta = sqrt(1.0-4.0*mtq/shat);

#include "auto/boxes.dec"

  double sigma0 =  M_PI * alphas * alphas 
    / 4.0 / (N*N-1.0) * beta/shat; 

  t = mtq - s/2.0*(1.0-beta*z);
  for(int i = 0; i < 2; i++ ) {
#include "auto/boxes.cpp"
    boxes_sum += Z_box + H_box * lambdat2 + W_box + Phi_box + Chi_box;
    z = -z;
    t = mtq - s/2.0*(1.0-beta*z);
  }
  
  return( boxes_sum );
}

double WeakCorrections::dsigmaWeakgg(const double mt_, const double shat, 
				     const double z){

  mt = mt_; mtq=mt*mt;
  double total = 0;

  if (ggtriangle)
    total += ggtriangles(shat,z);
  
  if (ggself)
    total += ggselfenergies(shat,z);

  if (ggbox)
    total += ggboxes(shat,z);
  
  if (ggvertex) 
    total += ggvertices(shat,z);
  
  return( hcq * total );

}

double WeakCorrections::A(double mq){
  complex<double> cint;
  int ierr=0;    
  ffxa0_(cint,0.,mu2,mq,ierr);
  return mycast(cint);
}

double WeakCorrections::B(double p1q, double m0q, double m1q){
  complex<double> cint;
  int ierr=0;    
  ffxb0_(cint, 0.,mu2,p1q , m0q, m1q, ierr);
  return mycast(cint);
}

double WeakCorrections::C(double p1q, double p2q, double p5q, 
          double m0q,double m1q, double m2q){
  /*
   * Use FF lib to calculate the finite integral:
   */

  complex <double> cint;
  double xxpi[6];
  int ierr=0;
  xxpi[0] = m0q;
  xxpi[1] = m1q;
  xxpi[2] = m2q;
  xxpi[3] = p1q; 
  xxpi[4] = p2q; 
  xxpi[5] = p5q; 
  ffxc0_(cint,xxpi,ierr);
  return mycast( cint);
}

double WeakCorrections::D(double p1q, double p2q, double p3q, double p4q, 
          double p5q, double p7q, 
          double m0q,double m1q, double m2q, double m3q){
 /*
   * Use FF lib to calculate the finite integral:
   */
  complex <double> cint;
  double xxpi[13];
  int ierr = 0;
  xxpi[0] = m0q;
  xxpi[1] = m1q;
  xxpi[2] = m2q;
  xxpi[3] = m3q;
  xxpi[4] = p1q; 
  xxpi[5] = p2q; 
  xxpi[6] = p3q; 
  xxpi[7] = p4q; 
  xxpi[8] = p5q;
  xxpi[9] = p7q;
  xxpi[10] = 0.0;
  xxpi[11] = 0.0;
  xxpi[12] = 0.0;  
  ffxd0_(cint,xxpi,ierr);
  return(mycast(cint));
}


double WeakCorrections::diffB0(double pq, double m0, double m1){
  complex<double> r,cint;
  
  if ( m0*m1 == 0.0 ) {
    cout << "wrong arguments in RediffB0, m0, m1 must be non-zero\n";
    exit(1);
  }
  double m0q = m0*m0;
  double m1q = m1*m1;
  r = ( m0q+m1q-pq
	+ sqrt( complex<double>( (-pq+m1q-2.0*m0*m1+m0q) * 
				 (-pq+m1q+2.0*m0*m1+m0q) ) ) ) /m0/m1/2.0;
  
  if ( r == complex<double>(1.0,0.0) ) 
    cint = - (m0q-m1q)/pq/pq*log(m1/m0) - 1.0/pq*2.0;
  else   
    cint = - (m0q-m1q)/pq/pq*log(m1/m0)
      + m0*m1/pq/pq*(1.0/r-r)*log(r)
      - 1.0/pq*(1.0+(r*r+1.0)/(r*r-1.0)*log(r));
  return(mycast(cint));
}

double WeakCorrections::f1(double x) {
  return( 
	 1.0 + 2.0 * ( + ( 1.0 + log(x) ) * ( 2.0*x + 3. ) 
		       - 2*pow(1.0+x,2)*(AbstractHathor::Li2(1.0+1./x)
					 - M_PI*M_PI/6. ) 
		       )
	  );
}

double WeakCorrections::F2xs(const double s, const double z, const double Fe, 
			     const double Fm, const double Fb) {
  double beta2 = 1.0-4.0*mtq/s, beta=sqrt(beta2);
  double sigma0 = 1./8.*M_PI*alphas*alphas*(N*N-1.)/N/N*beta/s;
  double dsigmadz = sigma0*(2.-beta2+beta2*z*z);

  return(sigma0*alpha/M_PI*((1.-beta2)/beta2*Fe*(1.-z*z) + Fm*(1.+z*z))
			    + dsigmadz*Fb);
    
}


void WeakCorrections::dsigmaWeakqq(const double mt_, const double s, 
				   const double z,double & up, double & down){


  /* 
   * Formulae are taken from Kühn,Scharf,Uwer, EPJ C45 (2005) 139
   *
   * see also below for a mistake in Eq.(35) and Eq.(36) 
   */
  mt = mt_; mtq=mt*mt;
  const double mzq=mz*mz,mwq=mw*mw,mhq=mh*mh,mbq=mb*mb;
  double rz = mzq/s, rw=mwq/s,rb=mbq/s,rh=mhq/s;
  double rz2 = rz*rz,rw2=rw*rw, rh2=rh*rh, rb2=rb*rb;
  double beta2 = 1.0-4.0*mtq/s, beta=sqrt(beta2), 
    beta4=beta2*beta2, beta6=beta4*beta2;

  double Amz = A(mzq);
  double Amt = A(mtq);
  double Amw = A(mwq);
  double Amb = A(mbq);
  double Amh = A(mhq);

  double B01_34 = B(mtq,mzq,mtq);
  double B03_13 = B(s,mtq,mtq);
  double B04_12 = B(mtq,mbq,mwq);
  double B04_13 = B(s,mbq,mbq);
  double B05_12 = B(mtq,mtq,mhq);

  double C03 = C(mtq,mtq,s,mtq,mzq,mtq);
  double C04 = C(mtq,mtq,s,mbq,mwq,mbq);
  double C05 = C(mtq,mtq,s,mtq,mhq,mtq);

  double F_e_Z = (2.0*rz*(gvt2+gat2) - beta2*(gvt2-3*gat2))*(B03_13-B01_34) 
    + 0.5*((4*rz2-beta2)*(gvt2+gat2)-beta4*(gvt2-3*gat2)
	   + 8*rz*beta2*gat2)*s*C03;  
    
  double F_m_Z = (gvt2+gat2)*(-0.5 + 2./( s*(1-beta2))*(Amz-Amt)
		       -(3./2.+rz/beta2)*B03_13
		       +(2+(rz*(1.-3.*beta2))/(beta2*(1-beta2)))*B01_34
		       -0.5/beta2*(beta2*(1+beta2)
				       + 4*rz*beta2+2*rz2)*s*C03); 

  double F_B_Z = alpha/2/M_PI*(2.*rz*(gvt2+gat2)+(1.-beta2)*(gvt2-3*gat2))*s
    * diffB0(mtq,mz,mt);



  double F_e_W = gw2*( (1+beta2+4*(rw-rb))*(B04_13-B04_12)
		       + 0.25*(pow(1+beta2+4*(rw-rb),2)-4*beta2)*s*C04); 

  double F_m_W = gw2*(-1.+4./( s*(1.-beta2))*(Amw-Amb)
		      + 0.5/(beta2*(1-beta2)) 
		      *(1.+4*beta2-5*beta4+4*(rw-rb)*(1-3*beta2))*B04_12 
		      -0.5/beta2*(1.+5*beta2+4*(rw-rb))*B04_13
		      -0.125/beta2*(1.+10*beta2+5*beta4+8*(rw-rb)
				    +8*beta2*(3*rw-rb)+16*pow(rw-rb,2))*s*C04); 

  double F_B_W = -alpha/2/M_PI*gw2*(1.-beta2-4*(rw-rb))*s*diffB0(mtq,mb,mw);

  double F_e_H = gw2*mtq/mwq*(2*(beta2+rh)*(B03_13-B05_12) 
			      -(beta2*(1-beta2)-3*rh*beta2-2.*rh2)*s*C05);

  double F_m_H = gw2*mtq/mwq*( -0.5 + 2./( s*(1-beta2))*(Amh-Amt)
			 + (0.5-rh/beta2)*B03_13 
			 + rh*(1.-3*beta2)/(beta2*(1-beta2))*B05_12 
			-1./beta2*(beta2*(1.-beta2)+rh2)*s*C05);

  double F_B_H = -alpha/2/M_PI*gw2*mtq/mwq*2*(1-beta2-rh)*s*diffB0(mtq,mt,mh);


  double F_e_chi = gat2*mtq/mzq*2*(2*rz*(B03_13-B01_34)+(2*rz+beta2)*rz*s*C03);
    
  double F_m_chi = gat2*mtq/mzq*(-1.+4/( s*(1.-beta2))*(Amz-Amt)
			  +(1-2*rz/beta2)*B03_13 
			  +2*rz*(1-3*beta2)/(beta2*(1-beta2))*B01_34
			  -2*rz2/beta2*s*C03); 

  double F_B_chi = alpha/2/M_PI*gat2*mtq/mzq*4*rz*s *diffB0(mtq,mt,mz);


  /*
   * In the journal version there is an error in F_e_phi and F_m_phi 
   * affecting the coefficients of (B04_13-B04_12) in F_e_phi 
   * and B04_13 in F_m_phi
   */
  double F_e_phi = gw2/8.*s/mwq*( -(beta4+8*beta2*rb+4*rw*beta2-4*rw-1.
				    -16*rw*rb+16*rb2)*(B04_13-B04_12)
				  +0.25*(1-beta6-3*beta2*(1-beta2)
					 -4*beta4*(rb+2*rw) 
					 + 8*beta2*rb+16*(1-beta2)*(rw2-rb2) 
					 + 64*rb*pow(rb-rw,2)
					 +4*(2*rw-rb))*s*C04);
  
  double F_m_phi = gw2/2*s/mwq* (-0.25*(1.-beta2+4*rb)
				 + (1-beta2+4*rb) / ( s*(1-beta2) )*(Amw-Amb) 
				 +0.125/(beta2*(1-beta2))*(1-3*beta2)
				 *(1-beta2+4*rb)*(1-beta2+4*(rw-rb))*B04_12
				 -0.125/beta2*(-1+beta2-4*rb)
				 *(3*beta2-4*rw-1+4*rb)*B04_13
				 -1./(32*beta2)*(1+beta2-5*beta4+3*beta6
						 +4*(2*rw-rb)
						 +64*rb*pow(rw-rb,2) 
						 +16*rw2*(1-beta2)
						 -16*rb2*(1-5*beta2)
						 +4*beta4*(2*rw+7*rb)
						 -8*beta2*(2*rw+3*rb))*s*C04);

  double F_B_phi = -alpha/16./M_PI*gw2*s/mwq*
    (pow(1-beta2,2)-4*(1-beta2)*(2*rb+rw)-16*rb*(rw-rb))*s*diffB0(mtq,mb,mw);
 

  double Fe = F_e_Z + F_e_W + F_e_H + F_e_chi + F_e_phi;
  double Fm = F_m_Z + F_m_W + F_m_H + F_m_chi + F_m_phi;
  double Fb = F_B_Z + F_B_W + F_B_H + F_B_chi + F_B_phi;

  up=down=F2xs(s,z,Fe,Fm,Fb);

  double sigma0 = 1./8.*M_PI*alphas*alphas*(N*N-1.)/N/N*beta/s;
  double dsigmadz = sigma0*(2.-beta2+beta2*z*z);

  up   += -1./8.*alpha/M_PI*((gvt2+gat2)*f1(rz)+2*gw2*f1(rw)) * 2 *dsigmadz;
  down += -1./8.*alpha/M_PI*((gvb2+gab2)*f1(rz)+2*gw2*f1(rw)) * 2 *dsigmadz;

  // Convert to pb:
  up *= hcq;
  down *= hcq;

  
}
