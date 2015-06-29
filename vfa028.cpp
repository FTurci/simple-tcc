
#include "global.h"
#include "vfa028.h"

vfa028::vfa028( simbox &box, particle_vector &pv, double rc, bool is )
  : voronoi( box, pv, rc, is )
{
  n_0_2_8 = 0;
  s_0_2_8.resize( N );
}

void vfa028::Clusters_Get_0_2_8()
{ // Find (0,2,8) clusters
  int i, j, k, cnt, n4, n5;
  int flg;

  for (i=0; i<N; i++) s_0_2_8[i] = 0;

  for (i=0; i<N; i++) { // loop over all particles, could perhaps improve by just looping over B particles if saying (0,2,8) has to be centered on a B particle
    n4=n5=flg=0;

    for (j=0; j<cnb[i]; j++) {  // loop j over all neighbours of partilce i
      cnt=0;
      for (k=0; k<cnb[i]; k++) {  // loop k over all neighbours of partilce i
        if (k==j) continue; // don't check if particle is bonded to itself
        if (Bonds_BondCheck(bNums[i][j],bNums[i][k])==1) cnt++; // check if two neighbours of particle i are bonded
        if (cnt>5) {  // too many neighbours of i are bonded to the jth neighbour of i
          flg=1;
          break;
        }
      }
      if (flg==1) break;  // i cannot be (0,2,8)

      if (cnt==2) { // too few neighbours of i are bonded to the jth neighbour of i
        flg=1;
        break;
      }
      else if (cnt==3) {  // too few neighbours of i are bonded to the jth neighbour of i
        flg=1;
        break;
      }

      else if (cnt==4) {
        n4++;
        if (n4>2) { // too many neighbours of i have four-bonds to the other neighbours of i
          flg=1;
          break;
        }
      }

      else if (cnt==5) {
        n5++;
        if (n5>8) { // too many neighbours of i have five-bonds to the other neighbours of i
          flg=1;
          break;
        }
      }
    }

    if (flg==1) continue; // i cannot be (0,2,8)

    if (n4==2 && n5==8) {
      // have found (0,2,8) cluster
      n_0_2_8++;

      s_0_2_8[i]=1;   // identify (0,2,8) particles in s_0_2_8 array
      for (j=0; j<cnb[i]; j++) s_0_2_8[bNums[i][j]]=1;
    }
  }
}

int vfa028::Find_N_0_2_8_N()
{ // Find fraction of particles within (0,2,8) clusters
  int i;
  int total_0_2_8_parts;

  Clusters_Get_0_2_8(); // find the (0,2,8)

  total_0_2_8_parts=0;
  for(i=0; i<N; i++) {
    total_0_2_8_parts+=s_0_2_8[i];
  }

  return total_0_2_8_parts;
}
