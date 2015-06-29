
#ifndef __CLUSTERS_TCC9A_H
#define __CLUSTERS_TCC9A_H

#include "tcc_base.h"

class tcc9A : public tcc_base
{
private:
  int nsp4b;    // number of sqaure pyramid clusters in configuration
  int msp4b;    // maximum number of square pyramids
  std::vector<std::vector<int> > sp4b;
                // msp4b x 5 array listing the particles within a square pyramid cluster

private:
  void Rings_gSP3( int );
  void Rings_gSP4( int, int, int );
  void Rings_aSP4( int, int, int, int );
  void Clusters_Get9A_D3h();

public:
  tcc9A( int );
  virtual int Find();
};

#endif
