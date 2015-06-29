
#ifndef __CLUSTERS_TCCBCC9_H
#define __CLUSTERS_TCCBCC9_H

#include "tcc_base.h"

class tccBCC9 : public tcc_base
{
private:
  int nsp4b, nsp4c;    // number of sqaure pyramid and Octahedral 6A clusters in configuration
  int msp4b, msp4c;    // maximum number of square pyramid and 6A
  std::vector<std::vector<int> > sp4b;
                // msp4b x 5 array listing the particles within a square pyramid cluster
  std::vector<std::vector<int> > sp4c;
                // msp4c x 6 array listing the particles within a 6A cluster

private:
  void Rings_gSP3( int );
  void Rings_gSP4( int, int, int );
  void Rings_aSP4( int, int, int, int );
  void Clusters_GetBCC9();

public:
  tccBCC9( int );
  virtual int Find();
};

#endif
