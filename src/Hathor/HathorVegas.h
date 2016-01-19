// $Modified: Mon Jun 28 16:09:47 2010 by uwer $
#ifndef HATHORVEGAS_H_
#define HATHORVEGAS_H_

#include <iostream>
#include <cmath>

extern "C" void ranlxd(double r[],int n); 

//#define FDIMMAX 64
#define FDIMMAX 1024
#define MAXDIM 4
#define NDMX 50

/*****************************************************************************
 * Vegas by G.P.Lepage
 *
 * --- Modified version to allow the integration of a vector ---
 *
 * References:
 * 
 * [1] G.P.Lepage, "Vegas: An Adaptive Multidimensional Integration Program", 
 *     CLNS-80/447
 *
 * [2] G.P.Lepage, "A New Algorithm for Adaptive Multidimensional Integration",
 *     J.Comput.Phys.27:192,1978
 *
 *****************************************************************************/

template <class Integrand> 
class Vegas {
  
 private:
  
  const int ndim;
  double alpha;
  
  int fdim;
  int ncall, itmx, nprn ;  
  int npg, ndo, ndm,  mds, it, nd, ng;
  int nxi[MAXDIM][NDMX], ia[MAXDIM], kg[MAXDIM];
  
  double si, si2, swgt, schi, scalls;
  double ti2, calls, ti, tsi, acc;
  
  double si_[FDIMMAX], si2_[FDIMMAX], swgt_[FDIMMAX], schi_[FDIMMAX];
  double ti2_[FDIMMAX], ti_[FDIMMAX], tsi_[FDIMMAX];
  double res[FDIMMAX];
  
  double dt[MAXDIM], qran[MAXDIM], x[MAXDIM], xl[MAXDIM], xu[MAXDIM], 
    r[NDMX], xin[NDMX], xi[MAXDIM][NDMX], 
    d[MAXDIM][NDMX], di[MAXDIM][NDMX], dx[MAXDIM];
  
  Integrand & integ;
  
 public:  
  double avgi, sd, chi2a;
  double avgi_[FDIMMAX], sd_[FDIMMAX], chi2a_[FDIMMAX];
  
 Vegas(Integrand & f, int flen=0) : 
  ndim(f.getDimension()), alpha(0.5), itmx(integ.VegasIterations()),
    nprn(0), ndo(1), mds(1), integ(f) {
    acc=integ.VegasAcc();
    ncall=integ.VegasPoints();
    /*
     * Set default for integration boundaries 
     */
    for (int i = 0; i < ndim; i++) {
      xi[i][0] = 1.0;
      xl[i] = 0.0; 
      xu[i] = 1.0;
    }
    fdim = flen;
  }
    
  void setnprn(int i){nprn = i;}
  
  void vegas() { vegas1();}		
  
  void vegas1(){
    it = 0;
    si = si2 = swgt = schi = scalls = 0.;
    for(int i=0; i<fdim; i++){
      si_[i] = 0; si2_[i] = 0; swgt_[i] = 0; schi_[i] = 0.;
    }
    vegas2();
  }
  
  void vegas2(){
    double xn, xo, wgt,xjac, rc, xnd, dxg, dv2g;
    int k, iaj, iaj1;
    nd = NDMX;
    ng = 1;
    
    if (mds != 0) {
      ng = static_cast<int>( pow(ncall * 0.5, 1.0 / ndim) );
      mds = 1;
      if ( (2 * ng - NDMX) >= 0) {
	mds = -1;
	npg = ng / NDMX + 1;
	nd = ng / npg;
	ng = npg * nd;
      }
    }
    
    k = static_cast<int>( pow(static_cast<double>(ng), ndim) );
    npg = ncall / k;
    if (npg < 2) npg = 2;
    calls = npg * k;
    dxg = 1.0 / ng;
    dv2g = pow(dxg, 2 * ndim) / ( npg * npg * (npg - 1.0) );
    xnd = nd;
    ndm = nd - 1;
    dxg = dxg * xnd;
    xjac = 1.0;

    for (int j = 0; j < ndim; j++) {
      dx[j] = xu[j] - xl[j];
      xjac *= dx[j];
    }

    if (nd != ndo) {//8
      rc = ndo / xnd;
      for (int j = 0; j < ndim; j++) {
	k = 0;
	xo = 0.;
	double dr = 1.0;
	xn = xi[j][k];

	for ( int i=0; i < ndm; i++) {
	  while ( rc > dr ) {
	    k++;
	    dr += 1.0;
	    xo = xn;
	    xn = xi[j][k];
	  }
	  dr -= rc;
	  xin[i] = xn - (xn - xo) * dr;
	}

	for (int i = 0; i < ndm; i++)
	  xi[j][i] = xin[i];
	xi[j][ndm] = 1.0;
      }
      ndo = nd;
    }
    
    if (nprn != 0) {
      std::cout << "INPUT PARAMETERS FOR VEGAS   NDIM=  " 
		<< ndim << "  NCALL=  " << calls << std::endl;
      std::cout << "                             IT=    " << it   
		<< "  ITMX =  " << itmx << std::endl;
      std::cout << "                             ACC=   " << acc << std::endl;
      std::cout << "                             MDS=   " << mds  
		<< "  ND= " << nd << std::endl;
    }
    
    do {//9 Main integration loop
      it++;
      ti = tsi = 0.;
      for (int i=0; i<fdim; i++){
	ti_[i] = 0.; tsi_[i] = 0.;
      }
      for (int j = 0; j < ndim; j++) {
	kg[j] = 1;
	for (int i = 0; i < nd; i++) {
	  nxi[j][i] = 0;
	  d[j][i] = 0.;
	  di[j][i] = 0.;
	}
      }
      
      do { //11
	for(kg[ndim-1]=1; kg[ndim-1]<=ng; kg[ndim-1]++) {
	  double fb=0, f2b = 0.;
	  double fb_[FDIMMAX]={0.}, f2b_[FDIMMAX] = {0.};
	  for ( int kk=1; kk <= npg; kk++) {
	    ranlxd(qran,ndim);
	    wgt = xjac;
	    for (int j = 0; j < ndim; j++) {
	      xn = (kg[j] - qran[j]) * dxg;
	      ia[j] = static_cast<int>(xn);
	      iaj = ia[j];
	      iaj1 = iaj - 1;
	      if (iaj > 0) {
		xo = xi[j][iaj] - xi[j][iaj1];
		rc = xi[j][iaj1] + (xn - iaj) * xo;
	      } else {
		xo = xi[j][iaj];
		rc = (xn - iaj) * xo;
	      }
	      x[j] = xl[j] + rc * dx[j];
	      wgt = wgt * xo * xnd;
	    }

	    const double f = integ.f(x, wgt,res) * wgt;

	    const double f1 = f / calls;
	    const double f2 = f * f;
	    fb +=  f;
	    f2b +=  f2;
	    //START----------------
	    for (int i=0;i<fdim;i++){
	      fb_[i] += res[i]*wgt;
	      f2b_[i] += res[i]*res[i]*wgt*wgt; 
	    }
	    //END------------------
	    for (int j = 0; j < ndim; j++) {
	      iaj = ia[j];
	      nxi[j][iaj]++;
	      di[j][iaj] += f1;
	      if (mds >= 0)
		d[j][iaj] +=  f2;
	    }
	  } // for(int kk=1;...

	  f2b = sqrt(f2b*npg);
	  f2b = (f2b - fb) * (f2b + fb);
	  ti  += fb;
	  tsi += f2b;
	  //START----------------
	  for (int i=0;i<fdim;i++){
	    f2b_[i] = sqrt(f2b_[i]*npg);
	    f2b_[i] = (f2b_[i] - fb_[i]) * (f2b_[i] + fb_[i]);
	    ti_[i]  += fb_[i];
	    tsi_[i] += f2b_[i];
	  }
	  //END------------------
	  if ( mds < 0 ) {
	    for (int j = 0; j < ndim; j++) {
	      d[j][ia[j]] += f2b;
	    }
	  }
	} //for(; kg[dim]
	
	for (k=ndim-2; k>=0; k--){
	  kg[k] = (kg[k] % ng) + 1;
	  if ( kg[k]!=1 )
	    break;
	}
	
      } while ( (kg[k] != 1) && (k >= 0) ); //11
      
      ti /= calls;
      tsi *= dv2g;
      ti2 = ti * ti;
      wgt = ti2 / tsi;
      si += ti * wgt;
      si2 +=  ti2;
      swgt += wgt;
      schi += ti2 * wgt;
      scalls += calls;
      avgi = si / swgt;
      sd = swgt * it / si2;
      chi2a = (it > 1 ) 
	? sd * (schi / swgt - avgi * avgi) / (it - 1) : 0. ;
      sd = 1.0 / sd;
      sd = sqrt(sd);
      tsi = sqrt(tsi);
      //>>>>
      //START----------------
      for (int i=0;i<fdim;i++){
	ti_[i] /= calls;
	tsi_[i] *= dv2g;
	ti2_[i] = ti_[i] * ti_[i];
	double wgti = ti2_[i] / tsi_[i];
	si_[i] += ti_[i] * wgti;
	si2_[i] +=  ti2_[i];
	swgt_[i] += wgti;
	schi_[i] += ti2_[i] * wgti;
	avgi_[i] = si_[i] / swgt_[i];
	sd_[i] = swgt_[i] * it / si2_[i];
	chi2a_[i] = (it > 1 ) 
	  ? sd_[i] * (schi_[i] / swgt_[i] - avgi_[i] * avgi_[i]) 
	  / (it - 1) : 0. ;
	sd_[i] = 1.0 / sd_[i];
	sd_[i] = sqrt(sd_[i]);
	tsi_[i] = sqrt(tsi_[i]);

      }
      //END------------------
      if (nprn != 0) {
	std::cout << std::endl ;
	std::cout << "INTEGRATION BY VEGAS" << std::endl;
	std::cout << "ITERATION NO  " << it << "   INTEGRAL =" 
		  << ti << std::endl;    
	std::cout << "                    STD DEV  =" << tsi << std::endl;
	std::cout << "ACCUMULATED RESULTS.   INTEGRAL =" <<  avgi << std::endl;
	std::cout << "                       STD DEV  =" << sd << std::endl;
	std::cout << "                       CHI**2 PER ITN   =" << chi2a 
		  << std::endl;
	if (fdim>0){
	  std::cout << std::endl << std::endl;
	  std::cout <<
	    "f[i]    |    this iter.     this err.    acc.   acc. err.   chi2 "
		    << std::endl;
	}
	for (int i=0;i<fdim;i++){
	std::cout << i << "  " 
		  << ti_[i] << "  " << tsi_[i] << "  "
		  << avgi_[i] << "  " << sd_[i] << "  "
		  << chi2a_[i] << std::endl;
	} 	
      }
      
      if ( 0 == ti ){
	/* Integrand is zero?, don't start second iteration */
	avgi = 0.; sd = 0;
	std::cout << "Vegas: Integral zero?" << std::endl;
	return;
      }
      for (int j = 0; j < ndim; j++) {
	dt[j] = 0.;
	for (int i = 0; i < nd; i++) {
	  if (  nxi[j][i] > 0 )
	    d[j][i] = d[j][i] / nxi[j][i];
	  dt[j] = dt[j] + d[j][i];
	}
      }

      for (int j = 0; j < ndim; j++) {
	rc = 0.;
	for (int i = 0; i < nd; i++) {
	  r[i] = 0.;
	  if ( d[j][i] > 0. ) {
	    xo = dt[j] / d[j][i];
	    r[i] = pow( (xo - 1.0) /  xo / log(xo) , alpha);
	  }
	  rc += r[i];
	}
	rc /=  xnd;
	k = 0;
	xo = 0.;
	double dr = 0.;
	dr = dr + r[k];
	xn = xi[j][k];
	
	for ( int i=0; i < ndm; i++) {
	  while (rc > dr) {
	    k++;
	    dr +=  r[k];
	    xo = xn;
	    xn = xi[j][k];
	  }
	  dr -= rc;
	  xin[i] = xn - (xn - xo) * dr / r[k];
	} 
	
	for (int i = 0; i < ndm; i++) {
	  xi[j][i] = xin[i];
	}
	xi[j][ndm] = 1.0;
      }
    } while ( ( it < itmx ) && ( acc < fabs(sd / avgi) ) ); //9
  }
  
};

#endif // HATHORVEGAS_H_
