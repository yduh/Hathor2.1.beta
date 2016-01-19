/* $Modified: Tue Apr 27 16:11:21 2010 by uwer $ */
#ifndef FITS_H_
#define FITS_H_

#define FUN(NAME) double NAME(const double & beta)

extern "C" {
  void setnf_(const int & nf); 
  double tfgg00_(const double & beta);
  double tfgg10_(const double & beta);
  double tfgg11_(const double & beta);
  double tfgg20approx_(const double & beta);
  double tfgg21_(const double & beta);
  double tfgg22_(const double & beta);
  double tfqqb00_(const double & beta);
  double tfqqb10_(const double & beta);
  double tfqqb11_(const double & beta);
  double tfqqb20approx_(const double & beta);
  double tfqqb21_(const double & beta);
  double tfqqb22_(const double & beta);
  double tfgq10_(const double & beta);
  double tfgq11_(const double & beta);
  double tfgq20approx_(const double & beta);
  double tfgq21_(const double & beta);
  double tfgq22_(const double & beta);
  double tfgg20_asym_(const double & beta);
  double tfqqb20_asym_(const double & beta);
  double tfgq20_asym_(const double & beta);
  double tfqqb20_exact_(const double & beta);
  FUN(tfgg20_);
  FUN(tfgq20_);
  FUN(tfqq20_);
  FUN(tfqqp20_);
  FUN(tfqqpb20_);
  FUN(tfqq21_);
  FUN(tfqqp21_);
  FUN(tfqq22_);
  void setab_(const int & i); 

#undef FUN

}
#endif
