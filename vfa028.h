
#ifndef __CLUSTERS_VFA028_H
#define __CLUSTERS_VFA028_H

#include "voronoi.h"

struct vfa028 : public voronoi
{
  int n_0_2_8;    // number of (0,2,8) clusters
  std::vector<int> s_0_2_8;
                  // N-length array saying if particle i is a membered of a (0,2,8) cluster

  vfa028( simbox&, particle_vector&, double, bool );

  void Clusters_Get_0_2_8();
  int Find_N_0_2_8_N();
};

#endif
