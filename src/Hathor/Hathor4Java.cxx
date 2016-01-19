// $MSTWed: Mon Jun 28 21:18:34 2010 by uwer $
#include "Hathor.h"
#include "HathorPdf.h"
#include <jni.h>
#include <iostream>
#include "LHAPDF/LHAPDF.h"

using namespace std;

/*
 * Java-Interface
 */

static double alphas;
static double pdfup;
static double pdfdown;

extern "C" JNIEXPORT jdouble JNICALL Java_JHathor_xs(JNIEnv *env, jobject, 
						     jint scheme, jint acc, 
						     jint collider, 
						     jdouble energy,
						     jstring pdfname,
						     jdouble Cqq, jdouble Cgg,
						     jdouble mt, jdouble mur, 
						     jdouble muf){

  const char *str;
  str = env->GetStringUTFChars(pdfname, NULL);

  cout << " HATHOR: use " << str << endl;

  Lhapdf pdf(str);
  Hathor xs(pdf);
  
  if ( 0==collider )
    xs.setColliderType(Hathor::PPBAR);
  else
    xs.setColliderType(Hathor::PP);
  
  xs.setSqrtShad(energy);
  xs.setScheme(scheme);
  xs.setPrecision(acc);
  xs.setCqq(Cqq);
  xs.setCgg(Cgg);
  
  double res = xs.getXsection(mt,mur,muf);
  alphas = xs.getAlphas(mur);
  if (scheme & Hathor::PDF_SCAN)
    xs.getPdfErr(pdfup,pdfdown);
  else {
    pdfup = 0;
    pdfdown = 0;
  }
  double result,err,chi;
  xs.getResult(0, result, err, chi);
  xs.PrintOptions();
  cout << " HATHOR: mt  / mur / muf / alphas / sigma / integration-err" << endl;
  cout << " HATHOR: " << mt << " " 
       << mur << " " 
       << muf << " " 
       << alphas << " "
       << res  << " "
       << err
       << endl;

  env->ReleaseStringUTFChars(pdfname, str);
  return(res);
}

extern "C" JNIEXPORT jdouble JNICALL Java_JHathor_alphas
(JNIEnv *env, jobject){
  return(alphas);
}

extern "C" JNIEXPORT jdouble JNICALL Java_JHathor_pdfup
(JNIEnv *env, jobject){
  return(pdfup);
}

extern "C" JNIEXPORT jdouble JNICALL Java_JHathor_pdfdown
(JNIEnv *env, jobject){
  return(pdfdown);
}

extern "C" JNIEXPORT jstring JNICALL Java_JHathor_pdfpath
(JNIEnv *env, jobject obj){
  string path = LHAPDF::pdfsetsPath();
  jstring comment = env->NewStringUTF(path.c_str());
  return comment;
}

extern "C" JNIEXPORT jint JNICALL Java_JHathor_lhapdfversion
(JNIEnv *env, jobject obj){
#ifdef LHAPDF_MAJOR_VERSION
  int version = LHAPDF_MAJOR_VERSION;
#else
  int version = 5;
#endif
  return version;
}

extern "C" JNIEXPORT jobjectArray JNICALL Java_JHathor_getDefaultCkmMatrix
(JNIEnv *env, jobject obj){
  MSTW mstw;
  HathorSgTopT xs(mstw);
  double ckm[3][3];
  xs.getCkmMatrix(ckm);

  jclass doubleArrayClass = env->FindClass("[D");
  jobjectArray jckm = env->NewObjectArray((jsize)3, doubleArrayClass, NULL);
  for (int i = 0; i < 3; i++) {
    jdoubleArray doubleArray = env->NewDoubleArray(3);
    env->SetDoubleArrayRegion(doubleArray, (jsize)0, (jsize)3, (jdouble*)ckm[i]);
    env->SetObjectArrayElement(jckm, (jsize)i, doubleArray);
    env->DeleteLocalRef(doubleArray);
  }
  return jckm;
}

unsigned int scheme_=Hathor::LO|Hathor::NLO|Hathor::NNLO;
int acc_=10000;;
int collider_=1;
double energy_=7000;
double Cqq_=0, Cgg_=0;
char pdfname[255] = "MSTW2008nnlo68cl";

extern "C" double xs(double mt, double mur, double muf){

  cout << " HATHOR: use " << pdfname << endl;

  Lhapdf pdf(pdfname);
  Hathor xs(pdf);
  
  if ( 0==collider_ )
    xs.setColliderType(Hathor::PPBAR);
  else
    xs.setColliderType(Hathor::PP);

  xs.setSqrtShad(energy_);
  xs.setScheme(scheme_);
  xs.setPrecision(acc_);
  xs.setCqq(Cqq_);
  xs.setCgg(Cgg_);

  double res = xs.getXsection(mt,mur,muf);
  alphas = xs.getAlphas(mur);
  if (scheme_ & Hathor::PDF_SCAN)
    xs.getPdfErr(pdfup,pdfdown);
  else {
    pdfup = 0;
    pdfdown = 0;
  }

  xs.PrintOptions();
  cout << " HATHOR: mt  / mur / muf / alphas / sigma " << endl;
  cout << " HATHOR: " << mt << " " 
       << mur << " " 
       << muf << " " 
       << alphas << " "
       << res << endl;

  return(res);
}

double xs_sgtop(SgTop &xs, int particle,
 int scheme, int acc, int collider, double energy,
 double ckm[3][3], double mt, double mur, double muf) {
  if ( 0==collider )
    xs.setColliderType(Hathor::PPBAR);
  else
    xs.setColliderType(Hathor::PP);

  xs.setSqrtShad(energy);
  xs.setScheme(scheme);
  xs.setPrecision(acc);
  switch (particle) {
    case 0: xs.setParticle(SgTop::TOPQUARK);     break;
    case 1: xs.setParticle(SgTop::ANTITOPQUARK); break;
    case 2: xs.setParticle(SgTop::BOTH);         break;
  }
  xs.setCkmMatrix(ckm);

  double res = xs.getXsection(mt,mur,muf);
  alphas = xs.getAlphas(mur);
  if (scheme & Hathor::PDF_SCAN)
    xs.getPdfErr(pdfup,pdfdown);
  else {
    pdfup = 0;
    pdfdown = 0;
  }
  double result,err,chi;
  xs.getResult(0, result, err, chi);
  xs.PrintOptions();
  cout << " HATHOR: mt  / mur / muf / alphas / sigma / integration-err" << endl;
  cout << " HATHOR: " << mt << " " 
       << mur << " " 
       << muf << " " 
       << alphas << " "
       << res  << " "
       << err
       << endl;

  return(res);
}

extern "C" JNIEXPORT jdouble JNICALL Java_JHathor_xssgtop
(JNIEnv *env, jobject, jint channel, jint particle,
 jint scheme, jint acc, jint collider, jdouble energy, jstring pdfname,
 jobjectArray jckm, jdouble mt, jdouble mur, jdouble muf){

  const char *str;
  str = env->GetStringUTFChars(pdfname, NULL);

  cout << " HATHOR: use " << str << endl;

  Lhapdf pdf(str);
  HathorSgTopT  xs_t  = HathorSgTopT(pdf);
  HathorSgTopS  xs_s  = HathorSgTopS(pdf);
  HathorSgTopWt xs_wt = HathorSgTopWt(pdf);

  double ckm[3][3];
  for (int i = 0; i < 3; i++) {
    jdoubleArray ar = (jdoubleArray)env->GetObjectArrayElement(jckm, i);
    jdouble *el = env->GetDoubleArrayElements(ar, 0);
    for (int j = 0; j < 3; j++)
      ckm[i][j] = el[j];
  }

  double res = 0.;
  switch (channel) {
    case 0: res = xs_sgtop(xs_t,  particle, scheme, acc, collider, energy,
                           ckm, mt, mur, muf);
            break;
    case 1: res = xs_sgtop(xs_s,  particle, scheme, acc, collider, energy,
                           ckm, mt, mur, muf);
            break;
    case 2: res = xs_sgtop(xs_wt, particle, scheme, acc, collider, energy,
                           ckm, mt, mur, muf);
            break;
    default: cout << "Invalid channel." << endl;
             exit(1);
  }

  env->ReleaseStringUTFChars(pdfname, str);

  return(res);
}

