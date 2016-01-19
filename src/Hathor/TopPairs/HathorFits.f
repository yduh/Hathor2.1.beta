C$Modified: Fri Jun  1 10:07:10 2012 by uwer $
! =====================================================================
      subroutine setnf(val)
      implicit none
      integer val
      integer nf
      common /nfnumber/ nf
      nf = val
      end      
! =====================================================================
      subroutine setAB(ischeme)
      implicit none
      integer ischeme
      real*8 Arxgg,Arxqq,Arxgq,Acqq20ab,Acgg20ab
      real*8 Brxgg,Brxqq,Brxgq,Bcqq20ab,Bcgg20ab
      real*8 rxgg,rxqq,rxgq,cqq20ab,cgg20ab
      common/rvals/rxgg,rxqq,rxgq,cqq20ab,cgg20ab
      real*8 ln2,ln22,z2
      parameter (ln2=0.693147180559945d0,ln22=0.480453013918201d0
     .     ,z2=1.6449340668482264365d0)
      integer nf
      common /nfnumber/ nf
C     Scheme A:
      Acqq20ab = 0.d0
      Acgg20ab = 0.d0
      Arxqq = -3.0d0
      Arxgq = -5.6d0
      Arxgg = -4.8d0
C     Scheme B:
      Bcqq20ab = (1276.d0/9.d0 -172.d0*ln2 + 256.d0/3.d0*ln22
     .     -86.d0/3.d0*z2 - 20.d0/9.d0*nf+8.d0/3.d0*nf*ln2)**2
      Bcgg20ab = (4444.d0/21.d0 - 2136.d0/7.d0*ln2+192.d0*ln22
     .     - 283.d0/7.d0*z2)**2
      Brxqq = -5.1d0
      Brxgq = -7.7d0
      Brxgg = -8.2d0
      if (0.eq.ischeme) then    
C     A-Scheme
         cqq20ab = Acqq20ab
         cgg20ab = Acgg20ab
         rxqq = Arxqq
         rxgq = Arxgq
         rxgg = Arxgg
      endif
      if (1.eq.ischeme) then    
C     B-Scheme
         cqq20ab = Bcqq20ab
         cgg20ab = Bcgg20ab
         rxqq = Brxqq
         rxgq = Brxgq
         rxgg = Brxgg
      endif
      if (2.eq.ischeme) then    
C     AB-Scheme
         cqq20ab = (Acqq20ab+Bcqq20ab)/2.d0
         cgg20ab = (Acgg20ab+Bcgg20ab)/2.d0
         rxqq = (Arxqq+Brxqq)/2.d0
         rxgq = (Arxgq+Brxgq)/2.d0
         rxgg = (Arxgg+Brxgg)/2.d0
      endif
      end
! ---------------------------------------------------------------------
      real*8 function tfqqb20_exact(beta)
      implicit none
C     Fits to exact 2-loop results taken from:
C     [P.Bärnreuther, M.Czakon, A.Mitov, arXiv:1204.5201]
      real*8 beta,rho,lnbeta,pi,z2
      real*8 tfqqb00,Lr,Fbeta0,Fbeta1,Ffit0,Ffit1,F2
      real*8 csx,dxfqqbns,xfqqpb
      integer nf
      common /nfnumber/ nf

      rho = 1d0 - beta**2

      pi = 2d0*dasin(1d0)
      lnbeta = dlog(beta)
      Lr = dlog(rho)

*----1204.5201------------
      F2 = tfqqb00(beta)/(108d0*pi**2) * (25d0 - 3d0*pi**2 
     ~+ 30d0*dlog(rho/4d0) + 9d0*dlog(rho/4d0)**2 ) 

      Fbeta1 = tfqqb00(beta)*( (0.0116822d0 - 0.0277778d0*lnbeta)/beta     
     ~+ 0.353207d0*lnbeta - 0.575239d0*lnbeta**2 + 0.240169d0*lnbeta**3)

      Fbeta0 = tfqqb00(beta)*( 0.0228463d0/beta**2 
     ~+ ( -0.0333905d0 + 0.342203d0*lnbeta - 0.888889d0*lnbeta**2)/beta
     ~+ 1.58109d0*lnbeta + 6.62693d0*lnbeta**2 - 9.53153d0*lnbeta**3
     ~+ 5.76405d0*lnbeta**4 )

      Ffit1 = (0.90756d0 - 6.75556d0*beta + 9.18183d0*beta**2)*beta*rho    
     ~+ ( -0.99749d0 + 27.7454d0*beta - 12.9055d0*beta**2)*beta*rho**2   
     ~+ ( -0.0077383d0 - 4.49375d0*beta + 3.86854d0*beta**2)*beta*rho**3   
     ~- 0.380386d0*beta**4*rho**4 + Lr*rho*(1.3894d0 + 6.13693d0*rho   
     ~+ 8.78276d0*rho**2 - 0.0504095d0*rho**3) + Lr**2*0.137881d0*rho 

      Ffit0 = (-2.32235d0 + 44.3522d0*beta - 24.6747d0*beta**2)*beta*rho   
     ~+ (2.92101d0 + 224.311d0*beta + 21.5307d0*beta**2)*beta*rho**2   
     ~+ (2.05531d0 + 945.506d0*beta + 36.1101d0*beta**2   
     ~- 176.632d0*beta**3)*beta*rho**3 + 7.68918d0*beta**4*rho**4 
     ~+ Lr*rho*(3.11129d0 + 100.125d0*rho + 563.1d0*rho**2 
     ~+ 568.023d0*rho**3)

      tfqqb20_exact = (Ffit0+Fbeta0+(Ffit1+Fbeta1)*nf+F2*nf**2)
     .     /(4.d0*pi)**2

*----1207.0236------------
      csx = - 0.4768323995789214d0

      xfqqpb = - 0.740572d0 - 31.2117d0*beta**2 - 0.31495d0*beta**3
     ~+ 15.8601d0*beta**4 - 1.64639d0*beta**5 + 18.9767d0*beta**6
     ~+ Lr**2*(-3.16565d0*rho+12.3828d0*rho**2) + Lr*(-19.6977d0*rho
     ~- 16.1386d0*rho**2 + 4.17707d0*rho**3) 

      dxfqqbns = (1.53647d0*beta**3 + 10.7411d0*beta**4)*rho 
     ~- 24.3298d0*beta**4*rho**2 
     ~+ (-4.50719d0*beta**3 + 15.4975d0*beta**4)*rho**3 
     ~+ (2.90068d0*beta**3 - 4.98808d0*beta**4)*rho**4 
     ~- 1.26644d0*beta**20*lnbeta 
     ~+ Lr**2*(0.327143d0*rho - 10.7669d0*rho**2) 
     ~+ Lr*(3.86236d0*rho - 21.332d0*rho**2 + 17.4705d0*rho**3)

      tfqqb20_exact = tfqqb20_exact 
     ~     + (  csx*Lr - beta**2*dexp(xfqqpb) + dxfqqbns)/(4.d0*pi)**2

*----------------

      return
      end
! =====================================================================
      real*8 function tfgg20_asym(beta)
      implicit none
      real*8 z2,vc,nc,pi,inv4pi2
      parameter (z2=1.6449340668482264365d0,nc=3.d0,Vc=8.d0,
     .     pi=3.1415926535897932385d0,
     .     inv4pi2 = 0.63325739776461107150d-2)
      real*8 beta,rho,Lb
      real*8 C,gamma,eta
      real*8 tfgg20approx,tfgg00
      real*8 rxgg,rxqq,rxgq,cqq20ab,cgg20ab
      common/rvals/rxgg,rxqq,rxgq,cqq20ab,cgg20ab
      C = 12.6d0
      gamma = 0.84d0 ! Sven used: gamma = 0.844d0
      rho = 1d0 - beta**2
      eta = 1.d0/rho - 1.d0
      tfgg20_asym = beta**3 / pi * (
     .     + 3089.d0/2250.d0 * nc/Vc + 19696.d0/3375.d0 * nc 
     .     - ( 59.d0/90.d0 * nc/Vc + 4.d0/15.d0*nc ) * z2 )*inv4pi2
      tfgg20_asym = tfgg20_asym 
     .     *  ( -dlog(rho) + rxgg * eta**gamma/(C+eta**gamma) )
      tfgg20_asym = tfgg20_asym 
     .     + tfgg00(beta) * inv4pi2*inv4pi2 * cgg20ab
      tfgg20_asym = tfgg20_asym + tfgg20approx(beta)
      return
      end
!------------------------------------------------------------------------
      real*8 function tfqqb20_asym(beta)
      implicit none
      real*8 z2,vc,nc,pi,inv4pi2
      parameter (z2=1.6449340668482264365d0,nc=3.d0,Vc=8.d0,
     .     pi=3.1415926535897932385d0,
     .     inv4pi2 = 0.63325739776461107150e-2)
      real*8 beta,rho,Lb
      real*8 C,gamma,eta
      real*8 tfqqb20approx,tfqqb00
      real*8 rxgg,rxqq,rxgq,cqq20ab,cgg20ab
      common/rvals/rxgg,rxqq,rxgq,cqq20ab,cgg20ab
      C = 47.9d0
      gamma = 1.37d0 ! Sven used: gamma = 1.374d0
      rho = 1d0 - beta**2
      eta = 1.d0/rho - 1.d0

      tfqqb20_asym = beta**3  / pi * (
     .     2462.d0/3375.d0*nc-88463d0/81000d0/nc+235d0/648d0/nc**3
     .     -(1d0/15d0*nc+11d0/360d0/nc-7d0/72d0/nc**3) * z2)*inv4pi2
      tfqqb20_asym = tfqqb20_asym
     .     *  ( - dlog(rho) + rxqq * eta**gamma/(C+eta**gamma) )
      tfqqb20_asym = tfqqb20_asym 
     .     + tfqqb00(beta) * inv4pi2*inv4pi2 * cqq20ab
      tfqqb20_asym = tfqqb20_asym + tfqqb20approx(beta)
      return
      end
!------------------------------------------------------------------------
      real*8 function tfgq20_asym(beta)
      implicit none
      real*8 z2,vc,nc,pi,inv4pi2
      parameter (z2=1.6449340668482264365d0,nc=3.d0,Vc=8.d0,
     .     pi=3.1415926535897932385d0,
     .     inv4pi2 = 0.63325739776461107150e-2)
      real*8 beta,rho,Lb
      real*8 C,gamma,eta
      real*8 tfgq20approx
      real*8 rxgg,rxqq,rxgq,cqq20ab,cgg20ab
      common/rvals/rxgg,rxqq,rxgq,cqq20ab,cgg20ab
      C = 16.4d0
      gamma = 0.90d0 ! Sven used: gamma = 0.898d0
      rho = 1d0 - beta**2
      eta = 1.d0/rho - 1.d0
      tfgq20_asym = beta**5  / pi * (
     .     2462.d0/1125d0*nc-479d0/324d0/nc
     .     - (2d0/15d0*nc+7d0/36d0/nc) * z2 ) * inv4pi2
      tfgq20_asym = tfgq20_asym
     .     *  ( -dlog(rho) + rxgq * eta**gamma/(C+eta**gamma) )
      tfgq20_asym = tfgq20_asym + rho*tfgq20approx(beta)
      return
      end
!------------------------------------------------------------------------
!------------------------------------------------------------------------
      real*8 function tfgg00(beta)
!     Ref: NPB303 (1988) 607-633, Eq. (15)
      implicit none
      integer nf
      real*8 beta,rho,pi,Lb
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)
      Lb = dlog((1d0 + beta)/(1d0 - beta))

      rho = 1d0 - beta**2
      tfgg00 = pi*beta*rho/192d0*((rho**2 + 16d0*rho + 16d0) * Lb/beta
     &                               -28d0 - 31d0 * rho)
      return
      end
!------------------------------------------------------------------------
      real*8 function tfgg10(beta)
      implicit none
!     Fit to Czakon and Mitov, arXiv:0811:4119 [hep-ph]
!     changed 10.02.09 UL
      real*8 f(17),c(17),beta,rho,pi,ln2,b1,b2,b3,coul,g1,g2,g3,a1
      integer nf,i
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)
      ln2 = dlog(2d0)
      rho = 1d0 - beta**2

      call fitfunctions(f,beta)

cc    Fit coefficients
      c(1) = -0.89256322246229405215D1
      c(2) = 0.14990572831107905388D3
      c(3) = -0.14055601422883656703D3
      c(4) = -0.34115615355928585466D0
      c(5) = -0.24104983305594559225D1
      c(6) = 0.54733818886364359100D2
      c(7) = 0.90915480146693527772D2
      c(8) = -0.48840100786375938770D1
      c(9) = -0.17466779250042875647D0
      c(10) = 0.13470336282052909270D2
      c(11) = 0.22664827097020649681D2
      c(12) = 0.46072668228382836681D1
      c(13) = -0.67623423283924284200D2
      c(14) = -0.97039142697975481608D1
      c(15) = 0.65080508879278166115D2
      c(16) = 0.50966326006659876509D1
      c(17) = -0.20122253412990102013D2


!     exact part
      g1 = beta * dlog(8d0*beta**2) ** 2
      g2 = beta * dlog(8d0*beta**2)
      g3 = rho**2*(dlog((1d0 + beta)/(1d0 - beta)) -2d0*beta)

      b1 = 7d0/pi/128d0
      b2 = -61d0/pi/256d0
      b3 = (nf - 4)/1024d0/Pi
      coul = 11d0*Pi/6d0/1536d0
      a1 =   1111d0/2304d0/Pi - 283d0*pi/18432d0  + 5d0/256d0*ln2/Pi
     &    - 7d0/128d0*ln2**2/Pi

      tfgg10 = b1*g1 + b2* g2 + b3*g3 + coul + a1 * beta

      do i =1,17
            tfgg10 = tfgg10 + c(i) * f(i)
      enddo

      return
      end

!------------------------------------------------------------------------
      real*8 function tfgg11(beta)
      implicit none
!     Ref: NPB303 (1988) 607-633, Eq. (19)
!     13. Feb. 2009 UL
      real*8 beta,rho,pi,h1,h2,tfgg00,ddilog,Lb
      integer nf
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)
      rho = 1d0 - beta**2

      h1 = dlog((1d0 + beta)/2d0)**2 - dlog((1d0 - beta)/2d0)**2
     &     + 2d0 * ddilog((1d0 + beta)/2d0)
     &     - 2d0 * ddilog((1d0 - beta)/2d0)

      h2 = ddilog(2d0*beta/(1d0+beta)) - ddilog(-2d0*beta/(1d0-beta))

      Lb = dlog((1d0 + beta)/(1d0 - beta))

      tfgg11 = (pi/192d0*(2*rho * (59*rho**2 + 198*rho - 288)*Lb
     &                   + 12*rho * (rho**2 + 16*rho + 16) * h2
     &                   -  6*rho * (rho**2 - 16*rho + 32) * h1
     &                   -(4d0/15d0)*beta*(7449*rho**2-3328*rho+724))
     &          + 12*tfgg00(beta) * dlog(rho/4d0/beta**2))/8d0/pi**2

      return
      end

! ---------------------------------------------------------------------
      real*8 function tfgg20approx(beta)
      implicit none
      real*8 beta,lnbeta,pi
      real*8 tfgg00
      integer nf
      common /nfnumber/ nf
!     Ref. arXiv:0911.5166, Eq. (5) - (6)

      pi = 2d0*dasin(1d0)
      lnbeta = dlog(beta)

      tfgg20approx = tfgg00(beta)/(4.d0*pi)**4*(
     ~  + 6.854713813d+1/beta**2
     ~  - 3.791058422d+0/beta - 9.663111475d-1/beta*nf
     ~  + lnbeta * ( + 2.346888886d+3 + 2.196952908d+1*nf
     ~    + 2.866713181d+2/beta + 6.893057042d+0/beta*nf )
     ~  + lnbeta**2 * ( - 3.155721841d+2 - 1.193552921d+2*nf
     ~    + 4.963001070d+2/beta )
     ~  + lnbeta**3 * ( - 2.321581037d+3 + 8.533333333d+1*nf )
     ~  + lnbeta**4 * (+ 4.608d+3) )

      return
      end

! ---------------------------------------------------------------------
      real*8 function tfgg20(beta)
      implicit none
      real*8 beta,rho,lnbeta,lnrho,pi
      real*8 tfgg00
      integer nf
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)

      rho = 1d0 - beta**2
      lnrho = dlog(rho)
      lnbeta = dlog(beta)

      tfgg20 = tfgg00(beta)/(4.d0*pi)**4*(
     ~  + 6.854713813d+1/beta**2
     ~  - 3.791058422d+0/beta - 9.663111475d-1/beta*nf
     ~  + lnbeta * ( + 2.346888886d+3 + 2.196952908d+1*nf
     ~    + 2.866713181d+2/beta + 6.893057042d+0/beta*nf )
     ~  + lnbeta**2 * ( - 3.155721841d+2 - 1.193552921d+2*nf
     ~    + 4.963001070d+2/beta )
     ~  + lnbeta**3 * ( - 2.321581037d+3 + 8.533333333d+1*nf )
     ~  + lnbeta**4 * (+ 4.608d+3) )

* exact gg-result from 1303.6254
      tfgg20 = tfgg20 + 
     ~ nf**2*(
     ~ (6.44022d-4*beta - 4.8664d-4*beta**2 - 0.0324653d-4*lnrho**2)*rho     
     ~ +(-13.8424d-4*beta + 4.7366d-4*beta**2 - 2.91398d-4*lnrho)*rho**2    
     ~ + (8.43828d-4*beta - 2.78748d-4*beta**2 + 2.38971d-4*beta**3)*
     ~     rho**3 )/(4.d0*pi)**2

      tfgg20 = tfgg20 + 
     ~ nf*(
     ~ - 0.0195046d0*beta - 1.4717d0*beta**2 - 0.223616d0*beta**3    
     ~ + 0.499196d0*beta**5 + 1.32756d0*beta**7 
     ~ + 0.00466872d0*beta**3*lnbeta + 0.0321469d0*beta**6*lnbeta**2 
     ~ + (0.579781d0*lnrho**2 + 0.166646d0*lnrho**3)*rho    
     ~ + (-1.36644d0*lnrho + 2.24909d0*lnrho**2)*rho**2 )/(4.d0*pi)**2

      tfgg20 = tfgg20 + (
     ~ 581.27542d0*beta + 1251.4057d0*beta**2 - 60.478096d0*beta**3 
     ~ + 1101.2272d0*beta**4 - 2905.3858d0*beta**5 
     ~ + 629.9128d0*beta**4*lnbeta - 5.1891075d0*lnrho 
     ~ + (1200.741d0*lnrho + 162.50333d0*lnrho**2)*rho     
     ~ + (36.074524d0*lnrho - 1192.8918d0*lnrho**2 
     ~     - 1810.2849d0*beta)*rho**2     
     ~ + 1568.7591d0*beta*rho**3 - 461.21326d0*beta*rho**4 
     ~     + 121.6379d0*beta*rho**5 )/(4.d0*pi)**2

      return
      end

! ---------------------------------------------------------------------
      real*8 function tfgg21(beta)
!     13. Feb. 2009 UL
!     Ref. PRD 80, 054009 (arXiv:0906.5273), Eq. (A.5)
      implicit none
      real*8 beta,rho,pi,ln2,approx00,a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,
     &       a(17),b(17),c(17),f(17)
      integer nf,i
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)
      ln2 = dlog(2d0)
      rho = 1d0 - beta**2
      a = 0d0

      call fitfunctions(f,beta)

      b(1) = -0.41893146444814722238D1
      b(2) = 0.82350664056295042675D2
      b(3) = -0.87873119691146622380D2
      b(4) = 0.98025932783251437483D1
      b(5) = -0.11226855014431170876D1
      b(6) = 0.29518302253901308394D2
      b(7) = 0.48361106944274381146D2
      b(8) = -0.70626177024572682864D1
      b(9) = -0.80252261565496867482D-1
      b(10) = 0.70149377905789963124D1
      b(11) = 0.15005881396623939798D2
      b(12) = 0.38414244102933078200D1
      b(13) = -0.47021617893150508163D2
      b(14) = -0.80558337917689825702D1
      b(15) = 0.47027405354396281189D2
      b(16) = 0.42143805176352359112D1
      b(17) = -0.14995997323716623358D2

      c(1) = 0.12306772109269179547D0
      c(2) = -0.27580880638729862191D1
      c(3) = 0.31973927251127193202D1
      c(4) = -0.56233044586307285123D0
      c(5) = 0.32400482270261797861D-1
      c(6) = -0.92541787970821446847D0
      c(7) = -0.15715471195386525654D1
      c(8) = 0.35109759921606889967D0
      c(9) = 0.22793615804648390035D-2
      c(10) = -0.21030153083493005391D0
      c(11) = -0.63688407214636296543D0
      c(12) = -0.12959775586040827193D0
      c(13) = 0.19169021640898521897D1
      c(14) = 0.26755747230533596919D0
      c(15) = -0.18615442306802792378D1
      c(16) = -0.13795864739683204589D0
      c(17) = 0.58155055639729381126D0

      approx00 = 7d0/192d0 * pi * beta

      b3 = -4608d0
      b2 = 109920d0/7d0 - 18432d0*ln2
      b1 = 69.64718 - 248.1500536/beta
      b0 = 56.86772061/beta + 17.010070

      c2 = -64d0
      c1 = 4048d0/21d0 - 192d0*ln2
      c0 = -3.446528522/beta - 37.60200419

      a3 = b3
      a2 = b2 + nf * c2
      a1 = b1 + nf * c1
      a0 = b0 + nf * c0

      tfgg21 = approx00 * (a3 * dlog(beta)**3 + a2 * dlog(beta)**2
     &                     + a1 * dlog(beta) + a0)/(16d0*pi**2)**2

      do i = 1,17
            a(i) = b(i) + nf * c(i)
            tfgg21 = tfgg21 + f(i) * a(i)
      enddo

      return
      end

! ---------------------------------------------------------------------
      real*8 function tfgg22(beta)
!     13. Feb. 2009 UL
!     Ref. PRD 80, 054009 (arXiv:0906.5273), Eq. (A.6)
      implicit none
      real*8 beta,rho,pi,ln2,approx00,a0,a1,a2,b0,b1,b2,c0,c1,
     &       a(17),b(17),c(17),f(17)
      integer nf,i
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)
      ln2 = dlog(2d0)
      rho = 1d0 - beta**2
      a = 0d0

      call fitfunctions(f,beta)

      b(1) = 0.12227826357971659934D-1
      b(2) = -0.77856183592455153899D0
      b(3) = 0.13395569795681659581D1
      b(4) = -0.59108409449884185759D0
      b(5) = 0.24833334390921618940D-2
      b(6) = -0.23827213128265742371D0
      b(7) = -0.38868909605647599698D0
      b(8) = 0.28342152835139977545D0
      b(9) = 0.10876067670746644266D-3
      b(10) = -0.33838620324281485002D-1
      b(11) = -0.29071015969066110884D0
      b(12) = -0.11473653816963170572D0
      b(13) = 0.98929369111933271629D0
      b(14) = 0.24899069031778152525D0
      b(15) = -0.10609632095432079805D1
      b(16) = -0.13425338490565164468D0
      b(17) = 0.35935659955161491823D0


      c(1) = -0.38038632524673614310D-2
      c(2) = 0.87577660922313445824D-1
      c(3) = -0.10742267112088965482D0
      c(4) = 0.23827064628577118431D-1
      c(5) = -0.99760082284744132732D-3
      c(6) = 0.29329412090521597518D-1
      c(7) = 0.49061467754082112120D-1
      c(8) = -0.13737341195966558972D-1
      c(9) = -0.69861280959753105896D-4
      c(10) = 0.65837107808845763385D-2
      c(11) = 0.20893206587541243495D-1
      c(12) = 0.49541414933149772185D-2
      c(13) = -0.65534585247274909019D-1
      c(14) = -0.10466353469288984754D-1
      c(15) = 0.65591296235852678153D-1
      c(16) = 0.55121822531690348902D-2
      c(17) = -0.20950586659315014504D-1

      approx00 = 7d0/192d0 * pi * beta

      b2 =  1152d0
      b1 = -2568d0 + 2304d0*ln2
      b0 = -79.743121401032

      c1 =  16d0
      c0 = -16d0 + 16d0*ln2

      a2 = b2
      a1 = b1 + nf * c1
      a0 = b0 + nf * c0

      tfgg22 = approx00 * (a2 * dlog(beta)**2
     &                     + a1 * dlog(beta) + a0)/(16d0*pi**2)**2

      do i = 1,17
            a(i) = b(i) + nf * c(i)
            tfgg22 = tfgg22 + f(i) * a(i)
      enddo

      return
      end

! ----------------------------------------------------------------------
      real*8 function tfqqb00(beta)
!     Ref: NPB303 (1988) 607-633, Eq. (14)
      implicit none
      integer nf
      real*8 beta,rho,pi
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)

      rho = 1d0 - beta**2
      tfqqb00 = pi * beta * rho * (2d0 + rho)/27d0
      return
      end


! -----------------------------------------------------------------------

      real*8 function tfqqb10(beta)
      implicit none
!     Fit to Czakon and Mitov, arXiv:0811:4119 [hep-ph]
!     changed 10. Feb. 2009 UL
      real*8 beta,rho,pi,ln2,approx00,b1,b2,b3,coul,a1,c(17),f(17)
      integer nf,i
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)
      ln2 = dlog(2d0)
      rho = 1d0 - beta**2

      call fitfunctions(f,beta)

!     Fit coefficients
      c(1) = 0.71206030042166749076D-1
      c(2) = -0.12716999877731435587D1
      c(3) = 0.12409953645972994269D1
      c(4) = -0.40504433345222059745D-1
      c(5) = 0.20537372584649030226D-1
      c(6) = -0.31763337032897441141D0
      c(7) = -0.71439686025892158319D0
      c(8) = 0.11700022807041564582D-1
      c(9) = 0.14891794812794497077D-2
      c(10) = -0.14451496658703046734D0
      c(11) = -0.13906364471709057710D0
      c(12) = 0.10767564560584650012D-1
      c(13) = 0.49397844716572649140D0
      c(14) = -0.56738105040293545974D-2
      c(15) = -0.53741901335282913682D0
      c(16) = -0.50937802879072537262D-2
      c(17) = 0.18250365495575112787D0


!     exact part
      approx00 = rho*beta/72d0/pi
      b3 = (nf - 4)*beta*rho*(2d0 + rho)/27d0*
     &      (2d0/3d0*dlog(4d0/rho) - 10d0/9d0)/8d0/pi
      b2 = approx00 * 16d0/3d0*dlog(8d0 * beta**2)**2
      b1 = approx00 * (-82d0/3d0)*dlog(8d0 * beta**2)
      coul = -approx00*pi**2/6d0/beta
      a1 = 299d0/324d0/Pi - 43d0/1296d0*Pi + ln2/54d0/Pi 
     &     - 2d0/27d0*ln2**2/Pi

      tfqqb10 = b3 + b2 + b1 + coul + a1 * beta * rho 

      do i = 1,17
            tfqqb10 = tfqqb10 + f(i) * c(i)
      enddo

      return
      end

!------------------------------------------------------------------------
      real*8 function tfqqb11(beta)
      implicit none
!     Ref: NPB303 (1988) 607-633, Eq. 18
!     13. Feb. 2009 UL
      real*8 beta,rho,pi,tfqqb00,ddilog,Lb,La
      integer nf
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)
      rho = 1d0 - beta**2

      Lb = dlog((1d0 + beta)/(1d0 - beta))
      La = dlog(rho/4d0/beta**2)

      tfqqb11 = (16d0*pi/81d0 * rho * Lb
     &         + (127 - 6*nf + 48 * La)*tfqqb00(beta)/9d0)/8d0/pi**2 

      return
      end

! ---------------------------------------------------------------------
      real*8 function tfqqb20approx(beta)
      implicit none
!     Ref. arXiv:0911.5166, Eq. (5)
      real*8 beta,lnbeta,pi
      real*8 tfqqb00
      integer nf
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)
      lnbeta = dlog(beta)

      tfqqb20approx = tfqqb00(beta)/(4.d0*pi)**4*(
     ~  + 3.607744112d+0/beta**2
     ~  - 5.272824211d+0/beta + 1.844775827d+0/beta*nf
     ~  + lnbeta * ( + 2.496754699d+2 + 5.577627475d+1*nf
     ~     + 5.403845377d+1/beta - 4.386490844d+0/beta*nf )
     ~  + lnbeta**2 * ( + 1.046483122d+3 - 9.083813530d+1*nf
     ~     - 1.403677070d+2/beta )
     ~  + lnbeta**3 * ( - 1.505158864d+3 + 3.792592592d+1*nf )
     ~  + lnbeta**4 * ( + 9.102222222d+2 ) )

      return
      end

! ---------------------------------------------------------------------
      real*8 function tfqqb21(beta)
      implicit none
      real*8 beta,rho,pi,ln2,approx00,b0,b1,b2,b3,c0,c1,c2,a(17),
     &       b(17),c(17),f(17),tfqqb00,a0,a1,a2,a3
      integer nf,i
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)
      ln2 = dlog(2d0)
      rho = 1d0 - beta**2
      a = 0d0

      call fitfunctions(f,beta)

! nf - independent coeffs
      b(1) = -0.3375529691516198d0
      b(2) = 6.164600455759970d0
      b(3) = -7.593493073912857d0
      b(4) = 1.768059646904276d0
      b(5) = -0.08762489912696413d0
      b(6) = 1.845938857400107d0
      b(7) = 3.523088585321746d0
      b(8) = -1.037199323293990d0
      b(9) = -0.006036427734527996d0
      b(10) = 0.3519758077975276d0
      b(11) = 1.661772750000683d0
      b(12) = 0.3546028381444250d0
      b(13) = -5.079079442873409d0
      b(14) = -0.7358243011649564d0
      b(15) = 5.000835768715067d0
      b(16) = 0.3812187492866887d0
      b(17) = -1.581977962829636d0

! nf - part
      c(1) = -0.003049222256391878d0
      c(2) = 0.04834626907759019d0
      c(3) = -0.03815845749760503d0
      c(4) = -0.006326122842925337d0
      c(5) = -0.0008178303428226445d0
      c(6) = 0.01467247159742959d0
      c(7) = 0.02996648869425292d0
      c(8) = 0.001675197532242481d0
      c(9) = -0.00005955543416873546d0
      c(10) = 0.006032742013617652d0
      c(11) = 0.003500283928605217d0
      c(12) = -0.00008169988675434426d0
      c(13) = -0.01412202098087665d0
      c(14) = -0.0001396481624019000d0
      c(15) = 0.01603850585097242d0
      c(16) = 0.0002213488752513140d0
      c(17) = -0.005416751281030128d0

! The leading dlogs ...
      approx00 = pi*beta/9d0

      b3 = -8192d0/9d0
      b2 = 12928d0/3d0 - 32768d0/9d0*ln2
      b1 = -840.510648d0 + 70.18385354d0/beta
      b0 = -82.24670336d0/beta + 467.9040166d0

      c2 = -256d0/3d0
      c1 = 2608d0/9d0-2816d0/9d0*ln2
      c0 = 6.579736270d0/beta - 64.6142758d0

      a3 = b3
      a2 = b2 + nf * c2
      a1 = b1 + nf * c1
      a0 = b0 + nf * c0

      tfqqb21 = approx00 * (a3 * dlog(beta)**3 + a2 * dlog(beta)**2
     &                     + a1 * dlog(beta) + a0)/(16d0*pi**2)**2

      do i = 1,17
         tfqqb21 = tfqqb21 + f(i) * (b(i) + nf * c(i))
      enddo

      tfqqb21 = tfqqb21 
     &      + nf**2 * 6d0/(16d0*pi**2)**2 * (-2d0/3d0) 
     &        * tfqqb00(beta) * (2d0/3d0 * dlog(4d0/rho)-10d0/9d0)


      return
      end

! ---------------------------------------------------------------------
      real*8 function tfqqb22(beta)
!     12. Feb. 2009
!     Ref. PRD 80, 054009 (arXiv:0906.5273), Eq. (A.2)
      implicit none
      real*8 beta,rho,pi,ln2,approx00,b0,b1,b2,c0,c1,a(17),
     &       b(17),c(17),f(17),tfqqb00,a0,a1,a2
      integer nf,i
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)
      ln2 = dlog(2d0)
      rho = 1d0 - beta**2
      a = 0d0

      call fitfunctions(f,beta)

      b(1) = 0.37947055850913903857D0
      b(2) = -0.42513804092471165186D1
      b(3) = 0.29171609398058425708D1
      b(4) = 0.94994469664640028618D0
      b(5) = 0.10537528686368018953D0
      b(6) = -0.16968987425222360344D1
      b(7) = -0.26097718103497107484D1
      b(8) = -0.27215566700172012463D0
      b(9) = 0.78785461740931727269D-2
      b(10) = -0.47933826922430962866D0
      b(11) = -0.18217132030886794008D0
      b(12) = -0.40679715559499412066D-1
      b(13) = 0.54147193643441573596D0
      b(14) = 0.84044060231810111539D-1
      b(15) = -0.51918414221346771405D0
      b(16) = -0.43364515724043944100D-1
      b(17) = 0.15957987999816306239D0

      c(1) = -0.22411359606195057267D-2
      c(2) = 0.26855760938516506711D-1
      c(3) = -0.17771261333389585834D-1
      c(4) = -0.62612116927097226015D-2
      c(5) = -0.62062291008162207036D-3
      c(6) = 0.98099888998031007537D-2
      c(7) = 0.16311753506096732526D-1
      c(8) = 0.18249960003284008512D-2
      c(9) = -0.46265779861705765304D-4
      c(10) = 0.28617550235673534383D-2
      c(11) = 0.11145907414448274122D-2
      c(12) = 0.17424723698612763258D-3
      c(13) = -0.35959286749381003025D-2
      c(14) = -0.35339444745074357479D-3
      c(15) = 0.36330037916478041857D-2
      c(16) = 0.17914863228005325235D-3
      c(17) = -0.11516365864674351892D-2


!     The leading dlogs...
      approx00 = pi*beta/9d0

      b2 =  2048d0/9d0
      b1 = -7840d0/9d0 + 4096d0/9d0 * ln2
      b0 =  270.8972388

      c1 =  320d0/9d0
      c0 = -596d0/9d0 + 320d0/9d0 * ln2

      a2 = b2
      a1 = b1 + nf * c1
      a0 = b0 + nf * c0

!     putting everything together...
      tfqqb22 = approx00 * (a2 * dlog(beta)**2
     &                     + a1 * dlog(beta) + a0)/(16d0*pi**2)**2

      do i = 1,17
         tfqqb22 = tfqqb22 + f(i) * (b(i) + nf * c(i))
      enddo

      tfqqb22 = tfqqb22
     &  + 3d0*nf**2*(-2d0/3d0)**2/(16d0*pi**2)**2 * tfqqb00(beta)

      return
      end
!----------------------------------------------------------------------

      real*8 function tfgq10(beta)
      implicit none
!     Fit to Czakon and Mitov, arXiv:0811:4119 [hep-ph]
      real*8 f(15),c(15),eta,beta,rho,pi
      integer i

      pi = 2d0 * dasin(1d0)
      rho = 1d0 - beta**2

!     Basis functions
      f(1) = beta ** 5 * dlog(beta)
      f(2) = beta ** 5
      f(3) = beta ** 2 * rho * dlog(rho)
      f(4) = beta ** 2 * rho * dlog(rho) ** 2
      f(5) = beta ** 3 * rho * dlog(rho)
      f(6) = beta ** 3 * rho * dlog(rho) ** 2
      f(7) = beta ** 4 * rho * dlog(rho)
      f(8) = beta ** 4 * rho * dlog(rho) ** 2
      f(9) = beta ** 2 * rho * dlog(rho) ** 3
      f(10) = beta ** 2 * rho * dlog(rho) ** 4
      f(11) = beta ** 2 * rho * dlog(rho) ** 5
      f(12) = beta ** 4 * dlog(beta)
      f(13) = beta ** 4
      f(14) = beta ** 6 * dlog(beta)
      f(15) = beta ** 6

!     Fit coefficients
      c(1) = -0.46606991381972084042D-1
      c(2) = 0.30192672424987652733D0
      c(3) = -0.25397760661778883776D0
      c(4) = -0.99912922440706274255D-2
      c(5) = 0.39878716751613268232D0
      c(6) = -0.24441723594649294420D-1
      c(7) = -0.14178345977173127890D0
      c(8) = 0.18672866178156529517D-1
      c(9) = 0.23865595016548048119D-2
      c(10) = -0.33988404124914615341D-4
      c(11) = -0.88807954159926924735D-6
      c(12) = -0.14214989118958294437D-2
      c(13) = -0.26103969870391417974D0
      c(14) = -0.15089037700587308619D0
      c(15) = -0.15054872592708121318D-1

!     exact part
      tfgq10 =   5d0/144d0/Pi * beta**3 * dlog(beta)
     &        + (-73d0/108d0 + 5d0/6d0*dlog(2d0))/16d0/pi * beta**3

      do i = 1,15
         tfgq10 = tfgq10 + f(i) * c(i)
      enddo

      return
      end

!------------------------------------------------------------------------
      real*8 function tfgq11(beta)
      implicit none
!     Ref: NPB303 (1988) 607-633, Eq. (20)
!     13. Feb. 2009 UL
      real*8 beta,rho,pi,h1,ddilog,Lb
      integer nf
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)
      rho = 1d0 - beta**2

      h1 = dlog((1d0 + beta)/2d0)**2 - dlog((1d0 - beta)/2d0)**2
     &     + 2d0 * ddilog((1d0 + beta)/2d0)
     &     - 2d0 * ddilog((1d0 - beta)/2d0)

      Lb = dlog((1d0 + beta)/(1d0 - beta))

      tfgq11 = (4d0/9d0*rho*(14*rho**2 + 27*rho - 136) * Lb
     &            - 32d0/3d0 * rho * (2 - rho) * h1
     &            - 8d0/135d0 * beta * (1319*rho**2 - 3468*rho + 724))
     &             /(1536d0*pi)

      return
      end
! ---------------------------------------------------------------------
      real*8 function tfgq20(beta)
      implicit none
      real*8 beta,rho,lnbeta,pi
      real*8 Lr,Ffit0,Ffit1
      integer nf
      common /nfnumber/ nf

      rho = 1d0 - beta**2

      pi = 2d0*dasin(1d0)
      lnbeta = dlog(beta)
      Lr = dlog(rho)

*----1210.6832------------
      Ffit1 = 0.363838d0*beta**2 - 1.44391d0*beta**3 + 1.1146d0*beta**7 
     ~- 0.309165d0*beta**3*lnbeta + 0.990057d0*beta**4*lnbeta**2 
     ~+ 0.362183d0*rho**2*Lr + (0.194867d0*rho + 1.57274d0*rho**2)*Lr**2 
     ~+ 0.0401411d0*rho*Lr**3 

      Ffit0 = 28.0998d0*beta**2 + 24.1753d0*beta**3 - 12.3211d0*beta**5 
     ~- 49.909d0*beta**7 + 11.7853d0*beta**3*lnbeta 
     ~+ 28.6697d0*beta**6*lnbeta**2 + (-1.68957d0 + 30.6335d0*rho**2)*Lr 
     ~+ (-9.80339d0*rho - 76.7407d0*rho**2)*Lr**2 - 3.82993d0*rho*Lr**3 

      tfgq20 = (Ffit0 + Ffit1*nf)/(4.d0*pi)**2
*----------------

      return
      end

! ---------------------------------------------------------------------
      real*8 function tfgq20approx(beta)
      implicit none
! csm     /date here/
!     Ref. PRD 80, 054009 (arXiv:0906.5273), Eq. (13)
      real*8 beta, tfgg00,a3,pi
      integer nf
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)

      a3 = pi*65d0/54d0

      tfgq20approx = beta**3 * (
     .     a3 * dlog(8.d0*beta**2)**3)/(16d0*pi**2)**2

      end

!------------------------------------------------------------------------
      real*8 function tfgq21(beta)
      implicit none
!     10. March 2009 UL
!     Ref. PRD 80, 054009 (arXiv:0906.5273), Eq. (A.3)
      real*8 beta,rho,pi,a(18),b(18),c(18),f(18),b2,b1,b0,c1,c0,
     &       a2,a1,a0 ,ln2
      integer nf,i
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)
      ln2 = dlog(2d0)
      rho = 1d0 - beta**2

      call fitfunctionsgq(f,beta)

      b(1) = -0.12053186893285318348733D-2
      b(2) = -0.49063528679563403272462D-1
      b(3) = -0.20885724666183118845068D0
      b(4) = -0.13731372244977537180718D2
      b(5) = 0.14018188404386073507504D2
      b(6) = -0.93048801637120825828700D-2
      b(7) = -0.52223667615322527066930D0
      b(8) = -0.46844051532118554322090D1
      b(9) = -0.76104616573715383957303D1
      b(10) = 0.13668774273125213200139D1
      b(11) = 0.18469829055221048908740D1
      b(12) = -0.72626598778080121804946D1
      b(13) = -0.48936402632047117522337D1
      b(14) = 0.11045667843237368689545D2
      b(15) = 0.41366018969437059321786D1
      b(16) = -0.63347705124519230661008D1
      b(17) = -0.10899543974745331643300D1
      b(18) = 0.11901056056132948056907D1

      c(1) = 0.32569675086132001407558D-4
      c(2) = 0.14275541264869577559384D-3
      c(3) = -0.40201696322625403789000D-2
      c(4) = 0.63298312610007168431600D-1
      c(5) = -0.59528252250125844431054D-1
      c(6) = 0.26940784931906550876628D-4
      c(7) = 0.15980360721538197059145D-2
      c(8) = 0.15226724258513526648356D-1
      c(9) = 0.28694384898798230551075D-1
      c(10) = -0.87558890588123621791996D-2
      c(11) = -0.80027092174209965701398D-2
      c(12) = 0.40434790805899279071914D-1
      c(13) = 0.19658783154071089450474D-1
      c(14) = -0.52622926106998094031341D-1
      c(15) = -0.14573949991897106947442D-1
      c(16) = 0.23146161288878520618164D-1
      c(17) = 0.29179200053640971870784D-2
      c(18) = -0.22011541343835690032074D-2

      b2 = 770d0/27d0
      b1 = -6805d0/81d0 + 6160d0/81d0*ln2
      b0 = 0.22068868d0 + 0.1370778390d0/beta

      c1 = 46d0/81d0
      c0 = -163d0/243d0 + 76d0/81d0*ln2

      a2 = - b2
      a1 = -(b1 + nf*c1)
      a0 = -(b0 + nf*c0)

      tfgq21 = 0d0
      tfgq21 = tfgq21 + beta**3/256d0/Pi**3 * (
     &            a2 * dlog(beta)**2 + a1 * dlog(beta) + a0)

      do i = 1,18
         tfgq21 = tfgq21 + f(i) * (b(i) + nf * c(i))
      enddo


      return
      end

!------------------------------------------------------------------------
      real*8 function tfgq22(beta)
      implicit none
!     10. March 2009 UL
!     Ref. PRD 80, 054009 (arXiv:0906.5273), Eq. (A.4)
      real*8 beta,rho,pi,a(18),b(18),c(18),f(18),a1,a0,ln2
      integer nf,i
      common /nfnumber/ nf

      pi = 2d0*dasin(1d0)
      ln2 = dlog(2d0)
      rho = 1d0 - beta**2

      call fitfunctionsgq(f,beta)

      b(1) = -0.22246839799351461611D-3
      b(2) = 0.50421878843683227151D-3
      b(3) = -0.29455040751936927466D-1
      b(4) = 0.34340411533975352371D0
      b(5) = -0.31894916719919378518D0
      b(6) = 0.92131181120549456986D-4
      b(7) = 0.69040164066403043799D-2
      b(8) = 0.78472327728269188833D-1
      b(9) = 0.16042051461979639693D0
      b(10) = -0.51869743021192863507D-1
      b(11) = -0.38610214184274521455D-1
      b(12) = 0.21650362064161570252D0
      b(13) = 0.10137656016574391131D0
      b(14) = -0.28056264283796080921D0
      b(15) = -0.80904694657086655177D-1
      b(16) = 0.13077889249832628608D0
      b(17) = 0.18138617319482943600D-1
      b(18) = -0.15857574767408567661D-1

      c(1) = 0.17886394829838176307D-4
      c(2) = 0.71275007337739231392D-6
      c(3) = -0.20581202660824972748D-3
      c(4) = 0.10875906288790151140D-2
      c(5) = -0.86284454134216735528D-3
      c(6) = 0.10164051081032675698D-6
      c(7) = 0.16380910241604766689D-4
      c(8) = 0.22729893902561754776D-3
      c(9) = 0.45698418149457103189D-3
      c(10) = -0.25619590149383899676D-3
      c(11) = -0.16025625676917503360D-3
      c(12) = 0.70713254738599103434D-3
      c(13) = 0.34937478554858063282D-3
      c(14) = -0.72546951880236132068D-3
      c(15) = -0.25524991490893572525D-3
      c(16) = 0.34015324829640885491D-3
      c(17) = 0.66131364068768656960D-4
      c(18) = -0.65620823765804253336D-4


      a1 = 385d0/81d0
      a0 = -1540d0/243d0 + 385d0/81d0*ln2

      tfgq22 = 0d0
      tfgq22 = tfgq22 + beta**3/256d0/Pi**3 * (
     &            a1 * dlog(beta) + a0)

      do i = 1,18
         tfgq22 = tfgq22 + f(i) * (b(i) + nf * c(i))
      enddo


      return
      end
! =====================================================================

      subroutine fitfunctions(f,beta)
      real*8 f(17),beta,rho

      rho = 1d0 - beta**2

      f(1) = beta ** 2
      f(2) = beta ** 3
      f(3) = beta ** 4
      f(4) = beta ** 5
      f(5) = beta ** 2 * dlog(beta)
      f(6) = beta ** 3 * dlog(beta)
      f(7) = beta ** 4 * dlog(beta)
      f(8) = beta ** 5 * dlog(beta)
      f(9) = beta ** 2 * dlog(beta) ** 2
      f(10) = beta ** 3 * dlog(beta) ** 2
      f(11) = beta * dlog(rho)
      f(12) = beta * dlog(rho) ** 2
      f(13) = beta ** 2 * dlog(rho)
      f(14) = beta ** 2 * dlog(rho) ** 2
      f(15) = beta ** 3 * dlog(rho)
      f(16) = beta ** 3 * dlog(rho) ** 2
      f(17) = beta ** 4 * dlog(rho)
      end
! =====================================================================

      subroutine fitfunctionsgq(f,beta)
      real*8 f(18),beta,rho

      rho = 1d0 - beta**2

      f(1) = beta ** 3
      f(2) = beta ** 4
      f(3) = beta ** 5
      f(4) = beta ** 6
      f(5) = beta ** 7
      f(6) = beta ** 4 * dlog(beta)
      f(7) = beta ** 5 * dlog(beta)
      f(8) = beta ** 6 * dlog(beta)
      f(9) = beta ** 7 * dlog(beta)
      f(10) = beta ** 3 * dlog(rho)
      f(11) = beta ** 3 * dlog(rho) ** 2
      f(12) = beta ** 4 * dlog(rho)
      f(13) = beta ** 4 * dlog(rho) ** 2
      f(14) = beta ** 5 * dlog(rho)
      f(15) = beta ** 5 * dlog(rho) ** 2
      f(16) = beta ** 6 * dlog(rho)
      f(17) = beta ** 6 * dlog(rho) ** 2
      f(18) = beta ** 7 * dlog(rho)

      end
! =====================================================================
! The dilog function
! Ref: ask Tord

      DOUBLE PRECISION FUNCTION DDILOG(X)
*     ====== ========= ======== =========
C
      DOUBLE PRECISION X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
      DOUBLE PRECISION C(0:18),H,ALFA,B0,B1,B2
C
      DATA ZERO /0.0D0/, ONE /1.0D0/
      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/
      DATA PI3 /3.28986 81336 96453D0/, PI6 /1.64493 40668 48226D0/
C
      DATA C( 0) / 0.42996 69356 08137 0D0/
      DATA C( 1) / 0.40975 98753 30771 1D0/
      DATA C( 2) /-0.01858 84366 50146 0D0/
      DATA C( 3) / 0.00145 75108 40622 7D0/
      DATA C( 4) /-0.00014 30418 44423 4D0/
      DATA C( 5) / 0.00001 58841 55418 8D0/
      DATA C( 6) /-0.00000 19078 49593 9D0/
      DATA C( 7) / 0.00000 02419 51808 5D0/
      DATA C( 8) /-0.00000 00319 33412 7D0/
      DATA C( 9) / 0.00000 00043 45450 6D0/
      DATA C(10) /-0.00000 00006 05784 8D0/
      DATA C(11) / 0.00000 00000 86121 0D0/
      DATA C(12) /-0.00000 00000 12443 3D0/
      DATA C(13) / 0.00000 00000 01822 6D0/
      DATA C(14) /-0.00000 00000 00270 1D0/
      DATA C(15) / 0.00000 00000 00040 4D0/
      DATA C(16) /-0.00000 00000 00006 1D0/
      DATA C(17) / 0.00000 00000 00000 9D0/
      DATA C(18) /-0.00000 00000 00000 1D0/
 
      IF(X .EQ. ONE) THEN
       DDILOG=PI6
       RETURN
      ELSE IF(X .EQ. MONE) THEN
       DDILOG=MALF*PI6
       RETURN
      END IF
      T=-X
      IF(T .LE. MTWO) THEN
       Y=MONE/(ONE+T)
       S=ONE
       A=-PI3+HALF*(DLOG(-T)**2-DLOG(ONE+ONE/T)**2)
      ELSE IF(T .LT. MONE) THEN
       Y=MONE-T
       S=MONE
       A=DLOG(-T)
       A=-PI6+A*(A+DLOG(ONE+ONE/T))
      ELSE IF(T .LE. MALF) THEN
       Y=(MONE-T)/T
       S=ONE
       A=DLOG(-T)
       A=-PI6+A*(MALF*A+DLOG(ONE+T))
      ELSE IF(T .LT. ZERO) THEN
       Y=-T/(ONE+T)
       S=MONE
       A=HALF*DLOG(ONE+T)**2
      ELSE IF(T .LE. ONE) THEN
       Y=T
       S=ONE
       A=ZERO
      ELSE
       Y=ONE/T
       S=MONE
       A=PI6+HALF*DLOG(T)**2
      END IF
 
      H=Y+Y-ONE
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO 1 I = 18,0,-1
      B0=C(I)+ALFA*B1-B2
      B2=B1
    1 B1=B0
      DDILOG=-(S*(B0-H*B2)+A)
      RETURN
      END


! ---------------------------------------------------------------------

      real*8 function tfqq20(beta)
      implicit none
      real*8 beta,rho,lnbeta,pi
      real*8 Lr,csx,xfqq
      integer nf
      common /nfnumber/ nf

      rho = 1d0 - beta**2

      pi = 2d0*dasin(1d0)
      lnbeta = dlog(beta)
      Lr = dlog(rho)

*----1207.0236------------
      csx = - 0.4768323995789214d0

      xfqq = - 0.740558d0 - 22.8129d0*beta**2 - 0.191648d0*beta**3 
     ~- 6.58031d0*beta**4 - 0.537669d0*beta**5 + 31.7872d0*beta**6 
     ~+ Lr**2*(-3.25313d0*rho + 15.8988d0*rho**2) + Lr*(-21.0783d0*rho 
     ~- 10.8176d0*rho**2 + 8.64557d0*rho**3) 

      tfqq20 = ( csx*Lr - beta**2*exp(xfqq) )/(4.d0*pi)**2
*----------------

      return
      end
!----------------------------------------------------------------------

      real*8 function tfqqp20(beta)
      implicit none
      real*8 beta,rho,lnbeta,pi
      real*8 Lr,csx,xfqqp
      integer nf
      common /nfnumber/ nf

      rho = 1d0 - beta**2

      pi = 2d0*dasin(1d0)
      lnbeta = dlog(beta)
      Lr = dlog(rho)

*----1207.0236------------
      csx = - 0.4768323995789214d0

      xfqqp = - 0.740558d0 - 23.4518d0*beta**2 - 0.193073d0*beta**3
     ~- 5.97215d0*beta**4 - 0.541402d0*beta**5 + 31.8227d0*beta**6
     ~+ Lr**2*(-3.29162d0*rho + 15.9932d0*rho**2) + Lr*(-21.3725d0*rho
     ~- 11.1642d0*rho**2 + 8.64746d0*rho**3) 

      tfqqp20 = ( csx*Lr - beta**2*dexp(xfqqp) )/(4.d0*pi)**2
*----------------

      return
      end
!----------------------------------------------------------------------

      real*8 function tfqqpb20(beta)
      implicit none
      real*8 beta,rho,lnbeta,pi
      real*8 Lr,csx,xfqqpb
      integer nf
      common /nfnumber/ nf

      rho = 1d0 - beta**2

      pi = 2d0*dasin(1d0)
      lnbeta = dlog(beta)
      Lr = dlog(rho)

*----1207.0236------------
      csx = - 0.4768323995789214d0

      xfqqpb = - 0.740572d0 - 31.2117d0*beta**2 - 0.31495d0*beta**3
     ~+ 15.8601d0*beta**4 - 1.64639d0*beta**5 + 18.9767d0*beta**6
     ~+ Lr**2*(-3.16565d0*rho+12.3828d0*rho**2) + Lr*(-19.6977d0*rho
     ~- 16.1386d0*rho**2 + 4.17707d0*rho**3) 

      tfqqpb20 = ( csx*Lr - beta**2*dexp(xfqqpb) )/(4.d0*pi)**2
*----------------

      return
      end

! ---------------------------------------------------------------------
      real*8 function tfqq21(beta)
      implicit none
      real*8 beta,rho,lnbeta,pi
      real*8 a(18),f(18)
      integer i

csm      pi = 2d0*dasin(1d0)

      a(1) = 0.00001383272300617876d0
      a(2) = 0.005596389340552027d0
      a(3) = 0.08106875480083188d0
      a(4) = -0.2781118398761645d0
      a(5) = 0.1977046777683356d0
      a(6) = 0.001279728034360501d0
      a(7) = 0.02997677540598500d0
      a(8) = 0.08266102969883754d0
      a(9) = -0.04561706606044799d0
      a(10) = 0.01802546474047812d0
      a(11) = 0.04273912867464031d0
      a(12) = -0.2543139108404375d0
      a(13) = -0.08420511868735086d0
      a(14) = 0.3573814574329435d0
      a(15) = 0.02316467394477266d0
      a(16) = -0.09071309320285033d0
      a(17) = 0.01828902864539101d0
      a(18) = -0.02899311678117691d0

      call fitfunctionsgq(f,beta)
      tfqq21 = 0d0
      do i = 1,18
            tfqq21 = tfqq21 + f(i) * a(i)
      enddo

csm      tfqq21 = tfqq21

      return
      end

!------------------------------------------------------------------------
      real*8 function tfqqp21(beta)
      implicit none
      real*8 scms,m2,beta,rho,lnbeta,pi
      real*8 a(18),f(18)
      integer i

csm      pi = 2d0*dasin(1d0)

      a(1) = 0.00001141914739815969d0
      a(2) = 0.004805421114065240d0
      a(3) = 0.06354433647150350d0
      a(4) = -0.2145055945769053d0
      a(5) = 0.1523638815296002d0
      a(6) = 0.001090030057141765d0
      a(7) = 0.02669797190799809d0
      a(8) = 0.08272125404332966d0
      a(9) = -0.02159096935671626d0
      a(10) = 0.005726223748346590d0
      a(11) = 0.03577974601059812d0
      a(12) = -0.2107243003496479d0
      a(13) = -0.06709366815637684d0
      a(14) = 0.3062136012969388d0
      a(15) = 0.01010509267456901d0
      a(16) = -0.06796885727994953d0
      a(17) = 0.02119576569747843d0
      a(18) = -0.03187223188420510d0

      call fitfunctionsgq(f,beta)
      tfqqp21 = 0d0
      do i = 1,18
         tfqqp21 = tfqqp21 + f(i) * a(i)
      enddo

csm      tfqqp21 = tfqqp21

      return
      end

!------------------------------------------------------------------------
      real*8 function tfqq22(beta)
      implicit none
      real*8 beta,rho,lnbeta,pi
      real*8 a(18),f(18)
      integer i

csm      pi = 2d0*dasin(1d0)

      a(1) = -0.00008017385357243393d0
      a(2) = 0.006148272192159677d0
      a(3) = 0.00008853309139700135d0
      a(4) = -0.007058857025945600d0
      a(5) = -0.00001961617963363922d0
      a(6) = 0.001248079789216109d0
      a(7) = 0.003736251736647367d0
      a(8) = 0.002677010535814378d0
      a(9) = -0.1268628200145811d-5
      a(10) = 0.0002178767656461540d0
      a(11) = 0.003186176605532504d0
      a(12) = -0.0008762216054656224d0
      a(13) = -0.001234451190039539d0
      a(14) = 0.002783757233927430d0
      a(15) = -0.005937067609989587d0
      a(16) = -0.001906132255281418d0
      a(17) = 0.003707035671613175d0

      call fitfunctions(f,beta)
      tfqq22 = 0d0
      do i = 1,17
         tfqq22 = tfqq22 + f(i) * a(i)
      enddo

csm      tfqq22 = tfqq22

      return
      end

!------------------------------------------------------------------------
