
#ifndef __TCC_BASE_H
#define __TCC_BASE_H

#include "voronoi.h"

/**
 * Base class for the TCC cluster detection.
 */
class tcc_base : public voronoi
{
public:
  int nClu;               // number of clusters
  std::vector<clusterinfo> clusters;
  std::vector<int> sClu;  // N-length array saying if particle i is a member of a cluster
  std::vector<int> sSpe;  // N-length array saying if particle i is a special member of a cluster

public:
  tcc_base( int _N ) : voronoi( _N ) {}
  virtual ~tcc_base() {}

  virtual int Find() = 0;
};

#endif
