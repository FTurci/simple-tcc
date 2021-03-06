
#include "global.h"
#include "supercool.h"
#include "tcc13A.h"

tcc13A::tcc13A( int _N ) : tcc_base( _N )
{
  try { fc = g_ctx.props.get<double>( "tcc.fc_13A" ); }
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

void tcc13A::Rings_gSP3(int n0)
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

void tcc13A::Rings_gSP4(int n0, int n1, int n2)
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

void tcc13A::Rings_gSP5(int n0, int n1, int n2, int n3) 
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

void tcc13A::Rings_aSP5(int n0, int n1, int n2, int n3, int n4)
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
    ++nsp5c;
  }
}

void tcc13A::Clusters_Get13A_Ih()
{ // Detect 13A D4d clusters
  //  Difficult to be desisive about. Made from 2 sp4c clusters with a common sp4 spindle
  // particle. Big gaps in the two 4 membered rings. Does work if the bond length is large
  // enough.
  int i, j, k, l, m;
  int sp1, sp2, nSB1, nSB2;
  int flg;

  for (i=0; i<N; ++i) sClu[i] = 0;
	
  for (i=0; i<nsp5c; ++i) { //first 7A
    sp1 = sp5c[i][5];
    sp2 = sp5c[i][6];
    nSB1 = nSB2 = 0; // count up spindle bonds

    for (j=i+1; j<nsp5c; ++j) { // second 7A
      flg = sp1 == sp5c[j][5] && Bonds_BondCheck(sp2,sp5c[j][6]);
      flg = flg || (sp1 == sp5c[j][6] && Bonds_BondCheck(sp2,sp5c[j][5]));
      if (flg==1) {
        if (nSB1>=5) {
	  nSB1++;
	  break;
	}
	nSB1++;
      }
    }

    if(nSB1 == 5) {	 // possibly found 13A, definately found 12B, now establish status
      for (j=i+1; j<nsp5c; ++j) {
        if (sp1 == sp5c[j][5] || sp1 == sp5c[j][6]) {
	  for (k=0; k<5; ++k) {
	    for (l=0; l<5; ++l) {
	      if (sp5c[i][k] == sp5c[j][l]) break;
	    }
	    if (l<5) break;
	  }
	  if (k==5) { // got 13A, make sure not detected previously Check all sp5c[j][ring] - sp1 sp5c are less than i
	    for  (k=0; k<i; ++k) {
	      for (l=0; l<5; ++l) {
	        if (sp5c[j][l] == sp5c[k][5] && sp1 == sp5c[k][6]) break;
		if (sp5c[j][l] == sp5c[k][6] && sp1 == sp5c[k][5]) break;
	      }
	      if(l<5) break; // index k < i present
	    }
	    if(k==i) break; // no index k < i present
	  }
	}
      }
      if (j<nsp5c) { // 13A found
        ++nClu;
	
	sClu[sp5c[i][5]]=1;
	sClu[sp5c[i][6]]=1;
	sClu[sp5c[j][5]]=1;
	sClu[sp5c[j][6]]=1;
	sClu[sp5c[i][0]]=1;
        sClu[sp5c[i][1]]=1;
	sClu[sp5c[i][2]]=1;
	sClu[sp5c[i][3]]=1;
	sClu[sp5c[i][4]]=1;
	sClu[sp5c[j][0]]=1;
	sClu[sp5c[j][1]]=1;
	sClu[sp5c[j][2]]=1;
	sClu[sp5c[j][3]]=1;
	sClu[sp5c[j][4]]=1;
      }
    }
		
    for (j=i+1; j<nsp5c; ++j) { // second 7A
      flg = sp2 == sp5c[j][5] && Bonds_BondCheck(sp1,sp5c[j][6]);
      flg = flg || (sp2 == sp5c[j][6] && Bonds_BondCheck(sp1,sp5c[j][5]));
      if (flg==1) {
        if (nSB2>=5) {
	  nSB2++;
	  break;
	}
	nSB2++;
      }
    }

    if(nSB2 == 5) {	 // possibly found 13A, definately found 12B, now establish status
      for (j=i+1; j<nsp5c; ++j) {
        if (sp2 == sp5c[j][5] || sp2 == sp5c[j][6]) {
	  for (k=0; k<5; ++k) {
	    for (l=0; l<5; ++l) {
	      if (sp5c[i][k] == sp5c[j][l]) break;
	    }
	    if (l<5) break;
	  }
	  if (k==5) { // got 13A, make sure not detected previously Check all sp5c[j][ring] - sp1 sp5c are less than i
	    for  (k=0; k<i; ++k) {
	      for (l=0; l<5; ++l) {
	        if (sp5c[j][l] == sp5c[k][5] && sp2 == sp5c[k][6]) break;
		if (sp5c[j][l] == sp5c[k][6] && sp2 == sp5c[k][5]) break;
	      }
	      if(l<5) break; // index k < i present
	    }
	    if(k==i) break; // no index k < i present
	  }
	}
      }
      if (j<nsp5c) { // 13A found
        ++nClu;
	
	sClu[sp5c[i][5]]=1;
	sClu[sp5c[i][6]]=1;
	sClu[sp5c[j][5]]=1;
	sClu[sp5c[j][6]]=1;
	sClu[sp5c[i][0]]=1;
        sClu[sp5c[i][1]]=1;
	sClu[sp5c[i][2]]=1;
	sClu[sp5c[i][3]]=1;
	sClu[sp5c[i][4]]=1;
	sClu[sp5c[j][0]]=1;
	sClu[sp5c[j][1]]=1;
	sClu[sp5c[j][2]]=1;
	sClu[sp5c[j][3]]=1;
	sClu[sp5c[j][4]]=1;
      }
    }
  }
}

int tcc13A::Find()
{
  int i;
  int total_13A_parts;

  for(i=0; i<N; i++) {  //  five membered shortest-path rings, then 7A
    Rings_gSP3(i);
  }
  Clusters_Get13A_Ih();  // find the 13A

  total_13A_parts=0;
  for(i=0; i<N; i++) {
    total_13A_parts+=sClu[i];
  }

  return total_13A_parts;
}
