
#include "global.h"
#include "supercool.h"
#include "voronoi.h"

#define EPS 0.000001  // for when comparing difference between two floats

voronoi::voronoi( int _N ) : N( _N )
{
  // global parameters
  fc = g_ctx.props.get( "tcc.fc", 1.0 );
                  // parameter for Williams modified Voronoi method, 1.0 is standard Voronoi method
  rcut = g_ctx.props.get( "tcc.rc", 2.0 );
                  // maximum permissible bond-length, setting this to something like 2.0 speeds up the calculation
                  // but don't set too low otherwise bonds will be missed
  rcut2 = rcut*rcut;
                  // squared cut-off length
  nB = 30;        // maximum number of bonds per particle

  // prepare arrays
  x.resize( N );
  y.resize( N );
  z.resize( N );
  cnb.resize( N );
  bNums.resize( N );
  for( int i = 0; i < N; i++ )
    bNums[i].resize( nB );
}

void voronoi::set_particles( simbox &box, particle_vector &pv, bool is )
{
  assert( pv.size() == N );

  // box dimensions
  // XXX: only works for cubic box!
  side = box.get_Lx();
  halfSide = .5*side;
  PBCs = 1;

  // copy particle data
  for( int n = 0; n < N; n++ ) {
    if( is ) {
      x[n] = pv[n].is_r[0];
      y[n] = pv[n].is_r[1];
#if NDIM == 3
      z[n] = pv[n].is_r[2];
#else
      z[n] = 0.0;
#endif
    } else {
      x[n] = pv[n].r[0];
      y[n] = pv[n].r[1];
#if NDIM == 3
      z[n] = pv[n].r[2];
#else
      z[n] = 0.0;
#endif
    }
  }
}

double voronoi::Bonds_GetR2(int i, int j)
{  // get squared separation between particles i and j
  double dx, dy, dz;

  dx = x[i] - x[j];
  dy = y[i] - y[j];
  dz = z[i] - z[j];
  return dx * dx + dy * dy + dz * dz;
}

double voronoi::Bonds_GetR2_PBCs(int i, int j)
{ // get PBC wrapped squared separation between particles i and j
  double dx, dy, dz;

  dx = x[i] - x[j];
  dy = y[i] - y[j];
  dz = z[i] - z[j];
  if (dx<-halfSide) dx+=side;
  else if (dx>halfSide) dx-=side;
  if (dy<-halfSide) dy+=side;
  else if (dy>halfSide) dy-=side;
  if (dz<-halfSide) dz+=side;
  else if (dz>halfSide) dz-=side;
  return dx * dx + dy * dy + dz * dz;
}

void voronoi::Bonds_GetDisp_PBCs(int i, int j,vector &axis)
{ // get PBC wrapped squared separation between particles i and j
  double dx, dy, dz;

  dx = x[i] - x[j];
  dy = y[i] - y[j];
  dz = z[i] - z[j];
  if (dx<-halfSide) dx+=side;
  else if (dx>halfSide) dx-=side;
  if (dy<-halfSide) dy+=side;
  else if (dy>halfSide) dy-=side;
  if (dz<-halfSide) dz+=side;
  else if (dz>halfSide) dz-=side;

  double r = sqrt(dx*dx+dy*dy+dz*dz);
  axis(0) = dx/r;
  axis(1) = dy/r;
  axis(2) = dz/r;
}

int voronoi::Bonds_BondCheck(int i, int j)
{ // Returns 1 if i & j are bonded; 0 otherwise
  int k;

  for (k=0; k<cnb[i]; ++k) {
    if (bNums[i][k] == j) return 1;
  }
  return 0;
}

void voronoi::Bonds_GetBondsV()
{ // Get bonds using Voronoi
  int rank = 0;
  int i, j, k, l, m;
  const int nBs = 4 * nB;
  int cnbs, cnbs2;
  int S[nBs], S2[nBs];
  double Sr[nBs], Sr2[nBs];
  double x1, x2, dr2;
  double rijx, rijy, rijz, rikx, riky, rikz, rjkx, rjky, rjkz;
  int Sb[nBs];

  for (i=0; i<N; ++i) { // reset bond lists
    cnb[i] = 0;
    for (j=0; j<nB; j++) bNums[i][j] = 0;
  }

  for (i=0; i<N; ++i) { // loop over all particles i
    cnbs = 0;
    for (j=0; j<N; ++j) { // loop over all particles j
      if (i==j) continue; // i cannot be bonded to itself
      if (PBCs == 1) dr2 = Bonds_GetR2_PBCs(i,j); // get squared separation of particles
      else dr2 = Bonds_GetR2(i,j);
      if (dr2 < rcut2) {
        if (cnbs < nBs) {  // max number of bonds, do ith particle
          k = cnbs++;
          S[k] = j;
          Sb[k] = 1;
          Sr[k] = dr2;
        }
        else {    // list is now full
          printf("d%d Bonds_GetBondsV(): nBs %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",rank,nBs,i,j);
          exit(1);
        }
      }
    } // We've now filled up S with all particles which could possibly be bonded to i

    cnbs2 = 0;  // Sort S so runs from closest particle to furthest particle from i
    for (j=0; j<cnbs; ++j) {  // maybe could use quickSort algorithm to make quicker here
      for(k=0; k<cnbs2; ++k) { // find spot to insert S[j]
        if (Sr[j] < Sr2[k]){
          for (l=cnbs2; l>k; --l) {
            S2[l] = S2[l-1];
            Sr2[l] = Sr2[l-1];
          }
          S2[k] = S[j];
          Sr2[k] = Sr[j];
          break;
        }
      }
      if (k==cnbs2){
        S2[cnbs2] = S[j];
        Sr2[cnbs2] = Sr[j];
      }
      ++cnbs2;
    } // Now sorted the list in order of distance from i

    for (j=0; j<cnbs2; ++j) Sb[j] = 1;

    for (l=0; l<cnbs2-1; ++l) { // Now eliminate particles from list which do not satisfy equation (4) from original TCC paper
      k = S2[l];
      for (m=l+1; m<cnbs2; ++m) {
        j = S2[m];
        rijx = x[i] - x[j];
        rijy = y[i] - y[j];
        rijz = z[i] - z[j];
        rikx = x[i] - x[k];
        riky = y[i] - y[k];
        rikz = z[i] - z[k];
        rjkx = x[j] - x[k];
        rjky = y[j] - y[k];
        rjkz = z[j] - z[k];
        if (PBCs==1) { // if PBCs are being used
          if (rijx>halfSide) rijx-=side;
          else if (rijx<-halfSide) rijx+=side;
          if (rijy>halfSide) rijy-=side;
          else if (rijy<-halfSide) rijy+=side;
          if (rijz>halfSide) rijz-=side;
          else if (rijz<-halfSide) rijz+=side;
          if (rikx>halfSide) rikx-=side;
          else if (rikx<-halfSide) rikx+=side;
          if (riky>halfSide) riky-=side;
          else if (riky<-halfSide) riky+=side;
          if (rikz>halfSide) rikz-=side;
          else if (rikz<-halfSide) rikz+=side;
          if (rjkx>halfSide) rjkx-=side;
          else if (rjkx<-halfSide) rjkx+=side;
          if (rjky>halfSide) rjky-=side;
          else if (rjky<-halfSide) rjky+=side;
          if (rjkz>halfSide) rjkz-=side;
          else if (rjkz<-halfSide) rjkz+=side;
        }
        x1 = rijx * rikx + rijy * riky + rijz * rikz;
        x1 -= rijx * rjkx + rijy * rjky + rijz * rjkz;
        x2 = rikx * rikx + riky * riky + rikz * rikz;
        x2 += rjkx * rjkx + rjky * rjky + rjkz * rjkz;
        x1 = x1 / x2;
        if (x1-fc > EPS) { // Eliminate particle m from S
          Sb[m] = 0;
        }
      }
    } // bond detection complete

    for (l=0; l<cnbs2; ++l) { // put bonds into bond lists
      if (Sb[l]) {
        j = S2[l];
        if (cnb[i] < nB && cnb[j] < nB) {  // do not max number of bonds for ith and jth particle
          k = cnb[i]++;
          bNums[i][k] = j;
        }
        else {    // list is now full
          printf("d%d Bonds_GetBondsV(): nB %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",rank,nB,i,j);
          exit(1);
        }
      }
    }
  } // End i loop
}
