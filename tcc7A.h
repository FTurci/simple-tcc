
#ifndef __CLUSTERS_TCC7A_H
#define __CLUSTERS_TCC7A_H

#include "tcc_base.h"

class tcc7A : public tcc_base
{
private:
  int nsp5c;    // number of pentagonal bipyramid 7A clusters in configuration
  int msp5c;    // maximum number of 7A
  std::vector<std::vector<int> > sp5c;
                // msp5c x 7 array listing the particles within a 7A cluster

private:
  void Rings_gSP3( int );
  void Rings_gSP4( int, int, int );
  void Rings_gSP5( int, int, int, int );
  void Rings_aSP5( int, int, int, int, int );

public:
  tcc7A( int );
  virtual int Find();
};

#endif
