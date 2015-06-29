
#ifndef __CLUSTERS_TCC11A_H
#define __CLUSTERS_TCC11A_H

#include "tcc_base.h"

class tcc11A : public tcc_base
{
private:
  int nsp4c;    // number of Octahedral 6A clusters in configuration
  int msp4c;    // maximum number of 6A
  std::vector<std::vector<int> > sp4c;
                // msp4c x 6 array listing the particles within a 6A cluster

private:
  void Rings_gSP3( int );
  void Rings_gSP4( int, int, int );
  void Rings_aSP4( int, int, int, int );
  void Clusters_Get11A_D4d();

public:
  tcc11A( int );
  virtual int Find();
};

#endif
