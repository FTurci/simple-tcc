
#include "global.h"
#include "supercool.h"
#include "tcc9A.h"

tcc9A::tcc9A( int _N ) : tcc_base( _N )
{
  try { fc = g_ctx.props.get<double>( "tcc.fc_9A" ); }
  catch(...) {}

  nsp4b = 0;
  msp4b = 5*N;
  sp4b.resize( msp4b );
  for( int n = 0; n < msp4b; n++ )
    sp4b[n].resize( 5 );

  nClu = 0;
  sClu.resize( N );
  sSpe.resize( N );
}

void tcc9A::Rings_gSP3(int n0)
{ // get SP3/4 rings including particle n0
  // ROUTINES FOR sp3 AND sp5 clusters REMOVED
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

void tcc9A::Rings_gSP4(int n0, int n1, int n2)
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

void tcc9A::Rings_aSP4(int n0, int n1, int n2, int n3)
{ // Take {n0,n1,n2,n3}, check SP4 ring and if so detect sp4a/b/c cluster
  // ROUTINES sp4a/sp4c clusters REMOVED (i.e. squares - only detecting square pyramids and octahedra)
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
  
  if (type==1) {
    if (nsp4b == msp4b) { printf( "Rings_aSP4(): msp4b too small\n" ); }
    sp4b[nsp4b][0] = n0;
    sp4b[nsp4b][1] = n1;
    sp4b[nsp4b][2] = n2;
    sp4b[nsp4b][3] = n3;
    sp4b[nsp4b][4] = cp[0];
    ++nsp4b;
  }
}

void tcc9A::Clusters_Get9A_D3h()
{ // Detect 9A clusters
   // Detected as 3 x sp4b, all three spindles unique and unbonded
   // 2x sp4b SP4 membered ring particles shared with each other cluster
	
	int i, j, k, l, m, n;
	int sp[3], db[2], ob[4];
	int flg;
	int s_com=-1;
	char errMsg[1000];
	
	for (i=0; i<nsp4b-2; ++i) {	// loop over all sp4b_i
		for (j=i+1; j<nsp4b-1; ++j) { // loop over all sp4b_j
			
			if (sp4b[i][4] == sp4b[j][4]) continue;	// spindles must be seperate
			if (Bonds_BondCheck(sp4b[i][4], sp4b[j][4])) continue;	// spindles must not be bonded
			sp[0]=sp4b[i][4];
			sp[1]=sp4b[j][4];

			m = 0; // two common particles between SP4 rings of sp4b_i and sp4b_j
			for(k=0; k<4; ++k) {
				for(l=0; l<4; ++l) {
					if(sp4b[i][k] == sp4b[j][l]){ 
						if(m<2) db[m] = sp4b[i][k];
						++m;
					}
				}
			}
			if (m!=2) continue;	 // checked	
				
			m = 0; // find SP4 particles of third sp4b cluster
			for (k=0; k<4; ++k) {
				if(sp4b[i][k] == db[0] || sp4b[i][k] == db[1]) continue;  // find particles in SP4 ring of sp4b_i that are not common to SP4 ring of sp4b_j
				for(l=0; l<4; ++l){
					if(sp4b[j][l] == db[0] || sp4b[j][l] == db[1]) continue;	  // find particles in SP4 ring of sp4b_j that are not common to SP4 ring of sp4b_i
					if(Bonds_BondCheck(sp4b[i][k], sp4b[j][l])) {	// check non-common SP4 ring particles from sp4b_i and sp4b_j are bonded
						if(m<4) ob[m] = sp4b[i][k];
						++m;
						if(m<4) ob[m] = sp4b[j][l];
						++m;
					}	
				}
			}
			if (m!=4) continue;  // done
			
			for(k=j+1; k<nsp4b; k++) {	// loop over all sp4b_k
				n = 0;  // check all four ring particles are in ob[.]
				for(l=0; l<4; ++l){
					for(m=0; m<4; ++m){
						if(sp4b[k][l] == ob[m]){
							++n;
							break;
						}	
					}
				}
				if (n != 4) continue; // done
				
				if (sp4b[k][4] == db[0]) continue;  // check sp4b_k spindle not in sp4b_i or sp4b_j
				if (sp4b[k][4] == db[1]) continue;
				if (sp4b[k][4] == sp[0]) continue;
				if (sp4b[k][4] == sp[1]) continue;
				
				// Now we have found the 9A D3h cluster
				++nClu;
				
				sp[2]=sp4b[k][4];
				
				// identify 9A particles in sClu array
				for (l=0; l<2; l++) sClu[db[l]]=1;
				for (l=0; l<3; l++) sClu[sp[l]]=1;
				for (l=0; l<4; l++) sClu[ob[l]]=1;
			}
		}
	}
}

int tcc9A::Find()
{
  int i;
  int total_9A_parts;

  for(i=0; i<N; i++) {  //  four membered shortest-path rings, then sp4b and 6A
    Rings_gSP3(i);
  }
  Clusters_Get9A_D3h();  // find the 9A

  total_9A_parts=0;
  for(i=0; i<N; i++) {
    total_9A_parts+=sClu[i];
  }

  return total_9A_parts;
}
