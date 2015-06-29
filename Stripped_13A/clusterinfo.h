
#ifndef __CLUSTERINFO_H
#define __CLUSTERINFO_H

struct clusterinfo
{
  // list of particles that are part of this cluster
  std::vector<int> particles;
  // index of central particle (if any)
  int c;
  // axis of the cluster
  // vector axis;

  size_t size() { return particles.size(); }

  clusterinfo() : c(-1) {}
};

#endif
