
#include "global.h"
#include "supercool.h"
#include "tcc7A.h"

tcc7A::tcc7A( int _N ) : tcc_base( _N )
{
  try { fc = g_ctx.props.get<double>( "tcc.fc_7A" ); }
  catch(...) {}

  nsp5c = 0;
  msp5c = 12*N;
  sp5c.resize( msp5c );
  for( int n = 0; n < msp5c; n++ )
    sp5c[n].resize( 12 );

  nClu = 0;
  sClu.resize( N );
  sSpe.resize( N );
}

void tcc7A::Rings_gSP3(int n0)
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

void tcc7A::Rings_gSP4(int n0, int n1, int n2)
{ // {n0,n1,n2} is not an SP3 ring, is it an SP4 ring?
  int i;
  int n3;

  for (i=0; i<cnb[n1]; ++i) {
    n3=bNums[n1][i];
    if (n3 <= n0) continue; // don:t find previously detected ring
    if (!Bonds_BondCheck(n0,n3)) {  // n1 not bonded to n2 & n0 not bonded to n3
      if (!Bonds_BondCheck(n2,n3)) { // n2 not bonded to n3
        Rings_gSP5(n0,n1,n3,n2);
      }
    }
  }
}

void tcc7A::Rings_gSP5(int n0, int n1, int n2, int n3) 
{ // {n0,n1,n2,n3} is not an SP4 ring, is it an SP5 ring?
  int i,j;
  int n4,n5;
  int bond4_1;
	
  for (i=0; i<cnb[n2]; ++i){
    n4=bNums[n2][i];
    if(n4 < n0 || n4 == n3) continue; // Now: is n4 bonded to n2 and not to n1 or n0
    bond4_1 = 0;
    for (j=0; j<cnb[n4]; ++j){
      n5=bNums[n4][j];
      if (n5==n3) bond4_1 = 1;
      if (n5==n1 || n5==n0) break; // Not SP ring
    }
    if (j==cnb[n4] && bond4_1==1) {
      Rings_aSP5(n0, n1, n2, n4, n3); // check SP5 type and store 
    }
  }
}

void tcc7A::Rings_aSP5(int n0, int n1, int n2, int n3, int n4)
{ // Take {n0,n1,n2,n3,n4}, check SP5 ring and if so detect SP5a/b/c cluster
  // ROUTINES SP5a/b clusters REMOVED (i.e. pentagons and pentagonal pyramids - only detecting pentagonal bipyramids)
  int i, j;
  int type = 0;
  int cp[2];  // common spindles - particles bonded to all members of four membered ring
  int bcheck;
  char errMsg[1000];

  cp[0]=cp[1]=-1;
  for (i=0; i<cnb[n0]; ++i) { // find out how many spindle particles bonded to five membered ring
    j = bNums[n0][i];
    bcheck = j == n1 || j == n4;
    if (bcheck) continue;
    bcheck = Bonds_BondCheck(n1,j)==1 && Bonds_BondCheck(n2,j)==1 && Bonds_BondCheck(n3,j)==1 && Bonds_BondCheck(n4,j)==1;
    if (bcheck) {
      if (type<2) {
        cp[type] = j;
        type++;
      }
      else type++;
    }
  }

  if (type==2) {
    if (nsp5c == msp5c) { printf( "Rings_aSP5(): msp5c too small\n" ); }
    sp5c[nsp5c][0] = n0;
    sp5c[nsp5c][1] = n1;
    sp5c[nsp5c][2] = n2;
    sp5c[nsp5c][3] = n3;
    sp5c[nsp5c][4] = n4;
    if (cp[0]<cp[1]) {
      sp5c[nsp5c][5] = cp[0];
      sp5c[nsp5c][6] = cp[1];
    }
    else {
      sp5c[nsp5c][5] = cp[1];
      sp5c[nsp5c][6] = cp[0];
    }

    for (i=0; i<7; i++) {
      sClu[sp5c[nsp5c][i]] = 1;
    }
    ++nsp5c;
  }
}

int tcc7A::Find()
{
  int i;
  int total_7A_parts;

  for (i=0; i<N; ++i) sClu[i] = 0;

  for(i=0; i<N; i++) {  //  five membered shortest-path rings, then 7A
    Rings_gSP3(i);
  }

  total_7A_parts=0;
  for(i=0; i<N; i++) {
    total_7A_parts+=sClu[i];
  }

  return total_7A_parts;
}
