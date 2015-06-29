
#ifndef __CLUSTERS_VORONOI_H
#define __CLUSTERS_VORONOI_H

/**
 * Base class for Voronoi tessellation.
 */
class voronoi
{
protected:
  int N;
  std::vector<double> x, y, z;
  double side;
  double halfSide;
  int PBCs;

  double fc;
  double rcut;
  double rcut2;
  int nB;
  std::vector<int> cnb;
  std::vector<std::vector<int> > bNums;

protected:
  voronoi( int );

  double Bonds_GetR2( int, int );
  double Bonds_GetR2_PBCs( int, int );
  void Bonds_GetDisp_PBCs( int, int, vector& );
  int Bonds_BondCheck( int, int );

public:
  void set_particles( simbox&, particle_vector&, bool );
  void Bonds_GetBondsV();
};

#endif
