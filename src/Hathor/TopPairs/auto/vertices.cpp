      t1 = C(mtq,0.0,t,mzq,mtq,mtq);
      t2 = gvt*gvt;
      t4 = beta-1.0;
      t5 = t4*t4;
      t6 = beta+1.0;
      t7 = t6*t6;
      t8 = t5*t7;
      t9 = t8*s;
      t13 = N2+N2*beta*z-2.0;
      t14 = alpha*t13;
      t15 = sigma0;
      t16 = t14*t15;
      t17 = 1./M_PI;
      t18 = beta*z;
      t19 = -1.0+t18;
      t20 = t19*t19;
      t21 = 1/t20;
      t23 = 1/N;
      t24 = t17*t21*t23;
      t25 = t16*t24;
      t28 = beta4*mzq;
      t29 = 6.0*t28;
      t30 = mzq*beta6;
      t31 = 8.0*t30;
      t32 = s*beta2;
      t33 = 3.0*t32;
      t34 = s*beta6;
      t35 = s*beta4;
      t37 = beta2*mzq;
      t38 = 8.0*t37;
      t39 = s*beta;
      t40 = t39*z;
      t41 = 3.0*t40;
      t42 = s*beta8;
      t43 = t42*z4;
      t44 = 2.0*t43;
      t45 = t42*z2;
      t46 = 6.0*t45;
      t47 = 2.0*mzq;
      t48 = 2.0*s;
      t50 = beta5*mzq*z;
      t51 = 10.0*t50;
      t52 = t28*z2;
      t53 = 8.0*t52;
      t54 = t32*z2;
      t56 = t37*z2;
      t57 = 4.0*t56;
      t58 = s*beta5;
      t59 = t58*z;
      t60 = 13.0*t59;
      t61 = s*beta7;
      t62 = t61*z;
      t64 = t35*z4;
      t66 = -t29+t31-t33+t34+6.0*t35-t38+t41-t44+t46+t47-t48-t51+t53+7.0*t54+
t57-t60+7.0*t62+10.0*t64;
      t67 = t34*z4;
      t68 = 6.0*t67;
      t71 = 2.0*mzq*beta*z;
      t72 = t35*z2;
      t74 = 4.0*t42;
      t75 = z5*beta5;
      t77 = 4.0*t75*mzq;
      t78 = beta3*z3;
      t80 = 12.0*t78*mzq;
      t83 = 12.0*z3*mzq*beta5;
      t84 = t34*z6;
      t85 = 2.0*t84;
      t86 = t61*z5;
      t88 = t61*z3;
      t90 = t58*z5;
      t91 = 8.0*t90;
      t92 = s*beta3;
      t93 = t92*z3;
      t95 = t58*z3;
      t96 = 26.0*t95;
      t97 = t34*z2;
      t99 = t92*z;
      t100 = 7.0*t99;
      t102 = 4.0*t30*z4;
      t104 = 12.0*t30*z2;
      t106 = beta3*mzq*z;
      t107 = 16.0*t106;
      t108 = -t68-t71-18.0*t72-t74-t77-t80+t83+t85+4.0*t86-12.0*t88-t91-14.0*
t93+t96+3.0*t97+t100+t102-t104+t107;
      t110 = B(t,mtq,mzq);
      t114 = 2.0*t18;
      t115 = -1.0-beta2+t114;
      t116 = 1/t115;
      t117 = 1/s;
      t120 = 1/(1.0+t18);
      t123 = 1/t20/t19;
      t126 = t116*t117*t120*t17*t123*t23;
      t129 = B(mtq,mtq,mzq);
      t131 = 22.0*t28;
      t134 = 16.0*t37;
      t135 = 6.0*t50;
      t136 = 32.0*t52;
      t137 = 2.0*t54;
      t138 = 20.0*t56;
      t139 = -t131+t31-t32+9.0*t34-5.0*t35+t134+t40-t44+t46+t47+s+t135+t136-
t137-t138;
      t140 = 3.0*t59;
      t141 = 2.0*t64;
      t143 = 10.0*t72;
      t145 = 8.0*t28*z4;
      t147 = 3.0*t99;
      t148 = 8.0*t106;
      t149 = t140-t62-t141+4.0*t67-t71+t143-t145-t74+t77+t80-t83-14.0*t97-t147+
t102-t104-t148;
      t153 = 1/t6;
      t154 = 1/t4;
      t156 = t153*t154*t117;
      t157 = t120*t17;
      t158 = t123*t23;
      t160 = t156*t157*t158;
      t163 = A(mtq);
      t165 = beta5*z;
      t167 = beta3*z;
      t169 = beta4*z4;
      t171 = 2.0*beta2;
      t173 = beta4*z2;
      t176 = beta2*z2;
      t178 = beta5*z3;
      t180 = t18-3.0*t165+4.0*t167-1.0+2.0*t169+beta4-t171-6.0*t78-4.0*t173-2.0
*t75+4.0*t176+6.0*t178;
      t183 = t14*t15*t116;
      t185 = t21*t23;
      t187 = t156*t157*t185;
      t190 = A(mzq);
      t196 = t15*t13;
      t197 = -s+t32-t47;
      t202 = -3.0*t176+t171-t167+t18+1.0-t169-2.0*beta4+3.0*t173;
      t205 = diffB0(mtq,mt,mz);
      t212 = t6*alpha;
      t214 = t196*t24;
      t217 = gat*gat;
      t223 = 7.0*t34;
      t225 = 6.0*t43;
      t226 = 18.0*t45;
      t229 = -t29+t31+5.0*t32-t223-10.0*t35-t38+t41+t225-t226+t47-t48-t51+t53-
t54+t57+35.0*t59-25.0*t62+t141;
      t230 = 6.0*t72;
      t231 = 12.0*t42;
      t234 = 2.0*t93;
      t236 = 11.0*t97;
      t237 = 9.0*t99;
      t238 = -t68-t71+t230+t231-t77-t80+t83+t85-12.0*t86+36.0*t88+t91+t234-38.0
*t95+t236-t237+t102-t104+t107;
      t249 = 6.0*t54;
      t250 = -t131+t31-9.0*t32-31.0*t34+27.0*t35+t134+t40+t225-t226+t47+s+t135+
t136+t249-t138;
      t255 = t140-t62+6.0*t64-12.0*t67-t71-30.0*t72-t145+t231+t77+t80-t83+42.0*
t97-t147+t102-t104-t148;
      t271 = 3.0*s;
      t272 = -t271+t33+t47;
      Z_vertex = -t1*t2*t9*t25/4.0-t16*(t66+t108)*t110*t2*t126/2.0-t129*t2*(
t139+t149)*t16*t160/2.0+2.0*t2*t163*t180*t183*t187-2.0*t2*t190*t180*t183*t187-
t196*alpha*t197*t202*t2*t205*t120*t24/2.0+t2*t4*t212*t214/2.0-t1*t217*t9*t25/
4.0-t16*(t229+t238)*t110*t217*t126/2.0-t129*t217*(t250+t255)*t16*t160/2.0+2.0*
t217*t163*t180*t183*t187-2.0*t217*t190*t180*t183*t187+t196*alpha*t272*t202*t217
*t205*t120*t24/2.0+t217*t4*t212*t214/2.0;
      t285 = 1/swq;
      t286 = t185*t285;
      t287 = t16*t286;
      t288 = mw*mw;
      t289 = 1/t288;
      t290 = t289*t17;
      t291 = t290*t5;
      t294 = s*s;
      t295 = C(mtq,0.0,t,mhq,mtq,mtq);
      t301 = t4*t6;
      t304 = 4.0*t43;
      t305 = 12.0*t45;
      t308 = 4.0*t64;
      t309 = 8.0*t42;
      t312 = 4.0*t90;
      t313 = -t33+5.0*t34+8.0*t35+t40-t304+t305+t54-23.0*t59+17.0*t62-t308+t68-
t309-t85+8.0*t86-24.0*t88-t312-t234;
      t314 = mhq*beta3;
      t316 = 12.0*t314*z3;
      t317 = mhq*beta2;
      t318 = t317*z2;
      t320 = mhq*beta6;
      t322 = 4.0*t320*z4;
      t323 = t314*z;
      t326 = 12.0*t320*z2;
      t329 = 2.0*mhq*beta*z;
      t330 = mhq*beta4;
      t331 = t330*z2;
      t333 = mhq*beta5;
      t334 = t333*z;
      t337 = 12.0*t333*z3;
      t339 = 4.0*t333*z5;
      t342 = 8.0*t320;
      t343 = 2.0*mhq;
      t344 = t96-t236+t99+t316-4.0*t318-t322-16.0*t323+t326+t329-8.0*t331+10.0*
t334-t337+t339+6.0*t330+8.0*t317-t342-t343;
      t346 = B(t,mtq,mhq);
      t350 = t158*t285;
      t353 = t350*t290*t116*t120;
      t356 = t15*t23;
      t358 = t14*t356*t123;
      t359 = t285*t289;
      t360 = t359*t17;
      t366 = -t33-21.0*t34+17.0*t35+t40+t304-t305-s+t249+t140-t62+t308-8.0*t67
-24.0*t72+t309+30.0*t97;
      t375 = -t147+t316-20.0*t318+t322-8.0*t323-t326-t329+32.0*t331+6.0*t334-
t337+t339-22.0*t330+16.0*t317+t342+t343-8.0*t330*z4;
      t378 = B(mtq,mtq,mhq);
      t384 = mhq-s+t32;
      t388 = diffB0(mtq,mt,mh);
      t391 = t359*t185;
      t397 = t196*t120;
      t401 = t116*t17*t289*t286;
      t403 = t163*t180*alpha*t397*t401/16.0;
      t404 = A(mhq);
      t411 = t14*t356*t21;
      t414 = t411*t360*t9/64.0;
      H_vertex = t287*t291*t7*(beta2+t114-3.0)*t294*t295/128.0-t16*t301*(t313+
t344)*t346*t353/64.0+t358*t360*(t366+t375)*t120*t378/64.0-t16*s*t4*t6*t384*t202
*t388*t157*t391/32.0-t403+t404*t180*alpha*t397*t401/16.0-t414;
      t415 = t285*t17;
      t418 = C(mtq,0.0,t,mwq,mbq,mbq);
      t424 = 3.0*t35;
      t431 = -t33+t34-t424+5.0*t40+t44-t46-t271+8.0*t54+9.0*t59-13.0*t62+12.0*
t64-10.0*t67;
      t435 = 2.0*t90;
      t438 = mbq*beta2;
      t440 = beta6*mwq;
      t441 = 16.0*t440;
      t442 = mbq*beta6;
      t443 = 16.0*t442;
      t444 = -t143+t74+4.0*t84-6.0*t86+18.0*t88-t435-18.0*t93+4.0*t97+t100+16.0
*t438+t441-t443;
      t446 = beta2*mwq;
      t448 = mwq*beta4;
      t450 = mbq*beta4;
      t452 = 4.0*mwq;
      t453 = 4.0*mbq;
      t454 = beta5*mwq;
      t456 = 24.0*t454*z3;
      t458 = 8.0*t440*z4;
      t460 = 24.0*t440*z2;
      t463 = 4.0*mwq*beta*z;
      t464 = beta3*mwq;
      t465 = t464*z;
      t467 = t454*z;
      t469 = t448*z2;
      t471 = -16.0*t446-12.0*t448+12.0*t450+t452-t453+t456+t458-t460-t463+32.0*
t465-20.0*t467+16.0*t469;
      t472 = t446*z2;
      t474 = mbq*beta5;
      t476 = 24.0*t474*z3;
      t478 = 8.0*t474*z5;
      t479 = mbq*beta3;
      t481 = 24.0*t479*z3;
      t483 = 8.0*t442*z4;
      t485 = 24.0*t442*z2;
      t486 = t479*z;
      t488 = t438*z2;
      t490 = t450*z2;
      t492 = t474*z;
      t494 = mbq*beta;
      t496 = 4.0*t494*z;
      t498 = 8.0*t454*z5;
      t500 = 24.0*t464*z3;
      t501 = 8.0*t472-t476+t478+t481-t483+t485-32.0*t486-8.0*t488-16.0*t490+
20.0*t492+t496-t498-t500;
      t506 = t117*t120;
      t507 = B(t,mbq,mwq);
      t512 = B(mtq,mbq,mwq);
      t515 = -t33-t223+t424+t40+t44-t46+t271-t249+t60-5.0*t62-2.0*t67+t230;
      t518 = 6.0*t93;
      t524 = t74-2.0*t86+6.0*t88+t435+t518-12.0*t95+6.0*t97-t237-32.0*t438+t441
-t443+32.0*t446-44.0*t448;
      t531 = 44.0*t450+t452-t453-t456+t458-t460-t463-16.0*t465+12.0*t467+64.0*
t469-40.0*t472+t476;
      t540 = -t478-t481-t483+t485+16.0*t486+40.0*t488-64.0*t490-12.0*t492+t496+
t498+t500-16.0*t169*mwq+16.0*t169*mbq;
      t547 = t154*t117;
      t548 = t547*t120;
      t553 = t452-t453-s+t32;
      t557 = diffB0(mtq,mb,mw);
      t563 = A(mbq);
      t566 = t116*t153;
      t567 = t196*t566;
      t570 = t548*t415*t185;
      t573 = A(mwq);
      W_vertex = t411*t415*t4*t6*mbq*t418/4.0-t16*t350*t17*(t431+t444+t471+t501
)*t116*t506*t507/16.0-t512*(t515+t524+t531+t540)*alpha*t196*t153*t548*t415*t158
/16.0+t196*alpha*t553*t202*t557*t120*t17*t286/16.0+t563*t180*alpha*t567*t570/
2.0-t573*t180*alpha*t567*t570/2.0+t14*t356*t21*t285*t17*t4*t6/8.0;
      t595 = mb*mb;
      t596 = t595*t595;
      t599 = 16.0*t596*beta*z;
      t600 = beta*beta;
      t601 = t600*t600;
      t602 = t601*t601;
      t604 = t294*t602*beta;
      t605 = t604*z;
      t607 = mbq*s;
      t609 = 8.0*t607*t18;
      t611 = t35*z4*mbq;
      t613 = s*mwq;
      t615 = 4.0*t613*t18;
      t616 = mwq*z;
      t617 = t92*t616;
      t619 = mwq*z2;
      t620 = t35*t619;
      t622 = t613*t176;
      t624 = t607*t167;
      t628 = t607*t173;
      t630 = t294+t599+13.0*t605-t609+80.0*t611-t615+36.0*t617+8.0*t620+8.0*
t622+88.0*t624+40.0*t607*t176-144.0*t628;
      t631 = t607*t165;
      t633 = beta6*z2;
      t634 = t607*t633;
      t636 = beta8*z2;
      t638 = 48.0*t607*t636;
      t639 = beta8*z4;
      t641 = 16.0*t607*t639;
      t642 = beta6*z4;
      t643 = t607*t642;
      t646 = t607*beta7*z;
      t649 = 48.0*t178*t613;
      t650 = t58*t616;
      t652 = t633*t613;
      t654 = t61*t616;
      t656 = beta7*z3;
      t658 = 24.0*t656*t613;
      t659 = beta7*z5;
      t661 = 8.0*t659*t613;
      t662 = -104.0*t631+40.0*t634+t638-t641-48.0*t643+56.0*t646+t649-52.0*t650
-40.0*t652+20.0*t654-t658+t661;
      t665 = 24.0*t636*t613;
      t667 = 8.0*t639*t613;
      t669 = 24.0*t78*t613;
      t671 = 8.0*t75*t613;
      t672 = t642*t613;
      t676 = 32.0*t474*z5*mwq;
      t677 = t438*t619;
      t679 = t479*t616;
      t681 = z3*mwq;
      t683 = 96.0*t479*t681;
      t685 = 96.0*t474*t681;
      t686 = t450*t619;
      t688 = t474*t616;
      t690 = t665-t667-t669-t671+8.0*t672-t676+32.0*t677+128.0*t679-t683+t685+
64.0*t686-80.0*t688;
      t692 = 16.0*t494*t616;
      t693 = z4*mwq;
      t695 = 32.0*t442*t693;
      t697 = 96.0*t442*t619;
      t712 = t602*t600;
      t713 = t294*t712;
      t714 = 4.0*t713;
      t715 = t596*beta6;
      t716 = 64.0*t715;
      t717 = t294*beta8;
      t719 = t294*beta6;
      t721 = -t692+t695-t697+32.0*t607*t659-96.0*t607*t656-64.0*t607*t75-80.0*
t92*z3*mbq+16.0*t607*beta6*z6+176.0*t607*t178-t714-t716+3.0*t717+8.0*t719;
      t724 = t596*beta2;
      t726 = t294*beta4;
      t728 = t294*beta2;
      t730 = t613*beta4;
      t732 = t607*beta4;
      t734 = t728*z2;
      t738 = t294*beta3;
      t739 = t738*z;
      t744 = t294*beta*z;
      t746 = t32*mwq;
      t748 = t607*beta2;
      t749 = 24.0*t748;
      t750 = 64.0*t724-4.0*t726-4.0*t728+4.0*t730+32.0*t732+4.0*t734-14.0*t726*
z2+18.0*t739+12.0*t726*z4-3.0*t744-20.0*t746-t749;
      t751 = 16.0*t596;
      t752 = 4.0*t613;
      t754 = 16.0*mbq*mwq;
      t755 = t596*beta4;
      t757 = t596*beta3;
      t758 = t757*z;
      t760 = t607*beta6;
      t763 = 32.0*t607*beta8;
      t765 = 96.0*t757*z3;
      t766 = t596*beta5;
      t768 = 32.0*t766*z5;
      t770 = 96.0*t766*z3;
      t772 = 32.0*t715*z4;
      t773 = t724*z2;
      t775 = -t751+t752+t754+48.0*t755-128.0*t758+8.0*t760-t763+t765+t768-t770-
t772-32.0*t773;
      t777 = t755*z2;
      t780 = 96.0*t715*z2;
      t781 = t34*mwq;
      t784 = 16.0*t42*mwq;
      t785 = t294*z4;
      t786 = t785*beta8;
      t788 = t294*beta5;
      t789 = t788*z;
      t791 = t766*z;
      t793 = t294*beta7;
      t795 = 4.0*t793*z5;
      t797 = 2.0*t788*z5;
      t798 = t793*z;
      t800 = t793*z3;
      t802 = t604*z3;
      t804 = -64.0*t777+t780+28.0*t781-t784+12.0*t786-6.0*t789+80.0*t791-t795-
t797-22.0*t798+26.0*t800-18.0*t802;
      t806 = 2.0*t785*t712;
      t808 = 6.0*t713*z2;
      t809 = t450*mwq;
      t811 = t438*mwq;
      t814 = 64.0*t442*mwq;
      t815 = t717*z2;
      t817 = t788*z3;
      t819 = t738*z3;
      t821 = t719*z2;
      t823 = t785*beta6;
      t829 = t604*z5;
      t831 = -t806+t808-48.0*t809-64.0*t811+t814-14.0*t815+2.0*t817-10.0*t819+
18.0*t821-22.0*t823+4.0*t719*z6-4.0*t717*z6+6.0*t829;
      t853 = -t294+t599+5.0*t605+t609-16.0*t611-t615-12.0*t617+104.0*t620-40.0*
t622-24.0*t624+48.0*t628;
      t861 = 24.0*t631-96.0*t634+t638-t641+32.0*t643-8.0*t646-t649+28.0*t650
-88.0*t652-12.0*t654+t658;
      t866 = -t661+t665-t667+t669+t671+24.0*t672+t676-160.0*t677-64.0*t679+t683
-t685;
      t871 = 8.0*t607;
      t874 = 256.0*t686+48.0*t688-t692+t695-t697-t714-t716+11.0*t717-6.0*t719-
t871-128.0*t724-6.0*t726;
      t883 = 6.0*t728-76.0*t730-56.0*t732-2.0*t734-10.0*t739+t744+28.0*t746+
t749-t751+t752+t754;
      t890 = 176.0*t755+64.0*t758+72.0*t760-t763-t765-t768+t770-t772+160.0*t773
-256.0*t777+t780+60.0*t781;
      t898 = -t784+4.0*t786+22.0*t789-48.0*t791-t795+t797-18.0*t798+18.0*t800
-6.0*t802-t806+t808;
      t913 = -176.0*t809+128.0*t811+t814-16.0*t815-18.0*t817+6.0*t819+12.0*t821
-2.0*t823+2.0*t829-16.0*t35*t693-64.0*t450*t693+64.0*t755*z4;
      t924 = 8.0*t748;
      t925 = 4.0*t746;
      t927 = -t754+t751-t871+t924-t752+t925-2.0*t728+t294+t726;
      t935 = t32-t453-s;
      t941 = t547*t157*t391;
      Phi_vertex = t287*t290*t4*t6*(-t32+4.0*t40+t453-t271)*mbq*t418/32.0-t16*(
t630+t662+t690+t721+t750+t775+t804+t831)*t507*t116*t506*t23*t123*t285*t290/
128.0-t16*t158*t359*t17*(t853+t861+t866+t874+t883+t890+t898+t913)*t120*t547*
t153*t512/128.0-t16*t927*t202*t557*t157*t289*t286/128.0-t563*t935*t180*alpha*
t567*t941/16.0+t573*t935*t180*t15*t14*t566*t941/16.0-t411*t360*t301*t935/64.0;
      t963 = t29-t31+t32+t34+t38+t40-t47+t51-t53-3.0*t54-t57+t59+t62-8.0*t64;
      t967 = t68+t71+12.0*t72+t77+t80-t83-t85+t312+t518-6.0*t95-7.0*t97-t100-
t102+t104-t107;
      t974 = -t131+t31+t32-t34+t35+t134+t40+t47-s+t135+t136+t137-t138;
      t977 = t140-t62-t71-4.0*t72-t145+t77+t80-t83+2.0*t97-t147+t102-t104-t148;
      t990 = s*mzq;
      Chi_vertex = -t287*t291*t7*t115*t294*t1/128.0-t16*t301*(t963+t967)*t110*
t353/64.0+t358*t360*(t974+t977)*t120*t129/64.0-t403+t190*t180*t15*t14*t120*t401
/16.0-t16*t990*t301*t202*t205*t157*t391/32.0-t414;
      t999 = N*t4*t6;
      t1000 = z2*alpha;
      t1001 = t1000*t15;
      t1002 = t999*t1001;
      t1003 = mz*mz;
      t1004 = t1003*t1003;
      t1005 = 4.0*t1004;
      t1007 = 1/t19;
      t1010 = t117*t17;
      t1011 = C(mtq,mtq,s,mtq,mzq,mtq);
      t1018 = (-t47+t32)*t2;
      t1020 = t1007*t120;
      t1021 = t1020*t1010;
      t1024 = B(s,mtq,mtq);
      t1045 = (t33+t47)*t217;
      Z_svertex = t1002*(-t1005+t728+t726)*t1007*t120*t1010*t2*t1011/2.0-t1002*
t1018*t129*t1021+t1002*t1018*t1024*t1021+t1002*beta2*t197*t205*t2*t1007*t157/
2.0-t1002*(-t728+t1005+8.0*t990*beta2+3.0*t726)*t217*t1011*t1021/2.0+t1002*
t1045*t129*t1021-t1002*t1045*t1024*t1021-t1002*beta2*t272*t205*t217*t1007*t157/
2.0;
      t1060 = t8*N*t1001;
      t1061 = t359*t1007;
      t1065 = mh*mh;
      t1066 = t1065*t1065;
      t1069 = C(mtq,mtq,s,mtq,mhq,mtq);
      t1075 = t32+mhq;
      t1086 = N*z2;
      t1087 = alpha*t15;
      t1090 = t359*t1020;
      t1091 = t17*beta2;
      H_svertex = t1060*t1061*t157*(t726-t728+3.0*mhq*s*beta2+2.0*t1066)*t1069/
32.0-t1060*t1061*t157*t1075*t378/16.0+t1060*t1061*t157*t1075*t1024/16.0+t9*
t1086*t1087*t1090*t1091*t384*t388/32.0;
      t1097 = 2.0*t39;
      t1103 = C(mtq,mtq,s,mbq,mwq,mbq);
      t1110 = t1020*t117;
      t1111 = t32+t452-t453+s;
      t1117 = B(s,mbq,mbq);
      W_svertex = t999*t1000*t15*(-t32+t1097-s-t452+t453)*(t32+t1097+s+t452-
t453)*t1103*t1007*t506*t415/32.0+t1002*t1110*t415*t1111*t512/8.0-t1002*t1110*
t415*t1111*t1117/8.0-t1002*t1020*t285*t1091*t553*t557/16.0;
      t1129 = t294*s;
      t1131 = t294*mwq;
      t1136 = t294*mbq;
      t1139 = t288*t288;
      t1144 = s*t596;
      t1160 = t1129*beta6+8.0*t1131*beta4-3.0*t1129*beta4+4.0*t1136*beta4+16.0*
t32*t1139-8.0*t1136*beta2-16.0*t1144*beta2+3.0*t1129*beta2-64.0*t596*t595+128.0
*t596*mwq-64.0*mbq*t1139-8.0*t1131-16.0*s*t1139+16.0*t1144-t1129+4.0*t1136;
      t1172 = t999*t1000*t15*(t726+t924+t925-t294+t751-t754-t752);
      t1186 = t289*t1007;
      Phi_svertex = t999*t1000*t15*t1160*t1103*t285*t289*t1021/256.0-t1172*t512
*t285*t289*t1021/64.0+t1172*t1117*t285*t289*t1021/64.0+t999*t1000*t15*t285*
t1186*t120*t1091*t927*t557/128.0;
      t1193 = t1086*t5;
      t1194 = t7*mzq;
      t1196 = t1193*t1194*t1087;
      t1205 = t1193*t1194*alpha;
      t1207 = t15*t289*t285;
      Chi_svertex = t1196*(t32+t47)*t1011*t285*t1186*t157/32.0-t1205*t1207*
t1020*t17*t129/16.0+t1205*t1207*t1020*t17*t1024/16.0+t1196*t1090*t17*s*beta2*
t205/32.0;
