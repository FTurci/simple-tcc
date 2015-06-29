
#include "global.h"
#include "supercool.h"
#include "tcc11A.h"

tcc11A::tcc11A( int _N ) : tcc_base( _N )
{
  try { fc = g_ctx.props.get<double>( "tcc.fc_11A" ); }
  catch(...) {}

  nsp4c = 0;
  msp4c = 6*N;
  sp4c.resize( msp4c );
  for( int n = 0; n < msp4c; n++ )
    sp4c[n].resize( 6 );

  nClu = 0;
  sClu.resize( N );
  sSpe.resize( N );
}

void tcc11A::Rings_gSP3(int n0)
{ // get SP3/4/5 rings including particle n0
  // ROUTINES FOR SP3 AND SP4 RINGS REMOVED
  int i,j;
  int n1, n2;

  for (i=0; i<cnb[n0]-1; i++){
    n1=bNums[n0][i];
    if (n1 < n0) continue;  // don:t find previously detected ring
    for (j=i+1; j<cnb[n0]; ++j){
      n2=bNums[n0][j];
      if (n2<n0) continue;  // don:t find previously detected ring
      if (!Bonds_BondCheck(n1,n2)) { // if n1 bonded to n2 have three membered ring
        if (n1<n2) Rings_gSP4(n0,n1,n2);
        else Rings_gSP4(n0,n2,n1);
      }
    }
  }
}

void tcc11A::Rings_gSP4(int n0, int n1, int n2)
{ // {n0,n1,n2} is not an SP3 ring, is it an SP4 ring?
  int i;
  int n3;

  for (i=0; i<cnb[n1]; ++i) {
    n3=bNums[n1][i];
    if (n3 <= n0) continue; // don:t find previously detected ring
    if (!Bonds_BondCheck(n0,n3)) {  // n1 not bonded to n2 & n0 not bonded to n3
      if (Bonds_BondCheck(n2,n3)) { // 4 membered ring found
        Rings_aSP4(n0,n1,n3,n2); // check SP4 type and store
      }
    }
  }
}

void tcc11A::Rings_aSP4(int n0, int n1, int n2, int n3)
{ // Take {n0,n1,n2,n3}, check SP4 ring and if so detect SP4a/b/c cluster
  // ROUTINES SP4a/b clusters REMOVED (i.e. squares and square pyramids - only detecting octahedra)
  int i, j;
  int type = 0;
  int cp[2];  // common spindles - particles bonded to all members of four membered ring
  int bcheck;
  char errMsg[1000];

  cp[0]=cp[1]=-1;
  for (i=0; i<cnb[n0]; ++i) { // find out how many spindle particles bonded to four membered ring
    j = bNums[n0][i];
    bcheck = j == n1 || j == n2;
    if (bcheck) continue;
    bcheck = Bonds_BondCheck(n1,j)==1 && Bonds_BondCheck(n2,j)==1 && Bonds_BondCheck(n3,j)==1;
    if (bcheck) {
      if (type<2) {
        cp[type] = j;
        type++;
      }
      else type++;
    }
  }

  if (type==2) {
    if (nsp4c == msp4c) { printf( "Rings_aSP4(): msp4c too small\n" ); }
    sp4c[nsp4c][0] = n0;
    sp4c[nsp4c][1] = n1;
    sp4c[nsp4c][2] = n2;
    sp4c[nsp4c][3] = n3;
    if (cp[0]<cp[1]) {
      sp4c[nsp4c][4] = cp[0];
      sp4c[nsp4c][5] = cp[1];
    }
    else {
      sp4c[nsp4c][4] = cp[1];
      sp4c[nsp4c][5] = cp[0];
    }
    ++nsp4c;
  }
}

void tcc11A::Clusters_Get11A_D4d()
{ // Detect 11A D4d clusters
  //  Difficult to be desisive about. Made from 2 sp4c clusters with a common sp4 spindle
  // particle. Big gaps in the two 4 membered rings. Does work if the bond length is large
  // enough.
  int i, j, k, l, m;
  int c, s0, s1;        // central and spindle particles

  for (i=0; i<N; ++i) sClu[i] = 0;

  for(i=0; i<nsp4c-1; ++i) {  // loop over all 6A
    for(j=i+1; j<nsp4c; ++j) {  // loop over all 6A
      c = -1;
      // check at least one common spindle between two 6A clusters
      if( sp4c[i][4] == sp4c[j][4] ) {
        c = sp4c[i][4];
        s0 = sp4c[i][5];
        s1 = sp4c[j][5];
      }
      if( sp4c[i][4] == sp4c[j][5] ) {
        c = sp4c[i][4];
        s0 = sp4c[i][5];
        s1 = sp4c[j][4];
      }
      if( sp4c[i][5] == sp4c[j][4] ) {
        c = sp4c[i][5];
        s0 = sp4c[i][4];
        s1 = sp4c[j][5];
      }
      if( sp4c[i][5] == sp4c[j][5] ) {
        c = sp4c[i][5];
        s0 = sp4c[i][4];
        s1 = sp4c[j][4];
      }
      if( c == -1 )
        continue;

      for(k=0; k<4; ++k) {  // check no overlap of particles within 4-rings of the two 6A clusters
        for(l=0; l<4; ++l) {
          if(sp4c[i][k] == sp4c[j][l]) break;
        }
        if(l<4) break;
      }
      if(k<4) continue;

      for(k=0; k<4; ++k) { // check each member of four-membered ring of 6A_i is bonded to two members of the four-membered ring in 6A_j
        m = 0;
        for(l=0; l<4; ++l) if(Bonds_BondCheck(sp4c[i][k], sp4c[j][l])) ++m;
        if(m!=2) break;
      }

      if (k!=4) continue;

      // now have found 11A cluster
      ++nClu;

      for (k=0; k<6; k++) { // identify 11A particles in s11A array
        sClu[sp4c[i][k]]=1;
        sClu[sp4c[j][k]]=1;
      }

      clusterinfo info;
      info.c = c;
      Bonds_GetDisp_PBCs(s0,s1,info.axis);
      clusters.push_back(info);

      // central particle
      sSpe[c] = 1;
    }
  }
}

int tcc11A::Find()
{
  int i;
  int total_11A_parts;

  for(i=0; i<N; i++) {  //  four membered shortest-path rings, then 6A
    Rings_gSP3(i);
  }
  Clusters_Get11A_D4d();  // find the 11A

  total_11A_parts=0;
  for(i=0; i<N; i++) {
    total_11A_parts+=sClu[i];
  }

  return total_11A_parts;
}
