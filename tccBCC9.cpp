
#include "global.h"
#include "supercool.h"
#include "tccBCC9.h"

tccBCC9::tccBCC9( int _N ) : tcc_base( _N )
{
  try { fc = g_ctx.props.get<double>( "tcc.fc_BCC9" ); }
  catch(...) {}

  nsp4b = 0;
  msp4b = 5*N;
  sp4b.resize( msp4b );
  for( int n = 0; n < msp4b; n++ )
    sp4b[n].resize( 5 );
    
  nsp4c = 0;
  msp4c = 6*N;
  sp4c.resize( msp4c );
  for( int n = 0; n < msp4c; n++ )
    sp4c[n].resize( 6 );

  nClu = 0;
  sClu.resize( N );
  sSpe.resize( N );
}

void tccBCC9::Rings_gSP3(int n0)
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

void tccBCC9::Rings_gSP4(int n0, int n1, int n2)
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

void tccBCC9::Rings_aSP4(int n0, int n1, int n2, int n3)
{ // Take {n0,n1,n2,n3}, check SP4 ring and if so detect SP4a/b/c cluster
  // ROUTINES SP4a clusters REMOVED (i.e. squares - only detecting square pyramids and octahedra)
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

  else if (type==2) {
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

void tccBCC9::Clusters_GetBCC9()
{ // Detect BCC9 clusters
  // WARNING: as not storing all BCC9 clusters cannot check for multiple dection of same cluster
  //          therefore nBCC_9 is probably wrong! But NBB9N and total_BCC9_parts will be ok...
  // BCC9 is detected in three ways
  // 1. As 2 x sp4b square pyramids: the spindles are common and each SP4 ring particle is bonded to 
  //    exactly one particle from the SP4 ring of the other cluster
  // 2. As 2 x sp4c: exaclty one spindle is common and each SP4 ring particle is bonded to 
  //    exactly one particle from the SP4 ring of the other cluster. The non-common spindles are not 
  //    part of the cluster
  // 3. As 1 x sp4b and 1 x sp4c: exaclty one spindle is common and each SP4 ring particle is bonded to 
  //    exactly one particle from the SP4 ring of the other cluster. The non-common spindle of the 
  //    sp4c is not part of the cluster 
	
	int i, j, k, l, m;
	int flg;
	int s_com=-1;
	char errMsg[1000];
	
	for (i=0; i<nsp4b-1; i++) { // loop over all sp4b_i
		for (j=i+1; j<nsp4b; ++j) { // loop over all sp4b_j
			if (sp4b[i][4]!=sp4b[j][4]) continue; // sp4b spindles must be common
			s_com=sp4b[i][4]; // common spindle
			
			flg=0; // check number of common particles between SP4 rings
			for (k=0; k<4; k++) {
				for (l=0; l<4; l++) {
					if (sp4b[i][k]==sp4b[j][l]) {
						flg=1;
						break;
					}
				}
				if (flg==1) break;
			}
			if (flg==1) continue; // zero common particles between SP4 rings
			
			flg=0; // check each SP4 ring particle of sp4b_i is bonded to exactly one in SP4_j
			for (k=0; k<4; k++) {
				m=0;
				for (l=0; l<4; l++) {
					if (Bonds_BondCheck(sp4b[i][k],sp4b[j][l])) {
						m++;
						if (m==2) {
							flg=1;
							break;
						}
					}
				}
				if (flg==1) break;
			}
			if (flg==1) continue; // checked
			
			flg=0; // check each SP4 ring particle of sp4b_j is bonded to exactly one in SP4_i
			for (k=0; k<4; k++) {
				m=0;
				for (l=0; l<4; l++) {
					if (Bonds_BondCheck(sp4b[j][k],sp4b[i][l])) {
						m++;
						if (m==2) {
							flg=1;
							break;
						}
					}
				}
				if (flg==1) break;
			}
			if (flg==1) continue; // checked

			// now have found BCC9 cluster
      		++nClu;

      		for (k=0; k<4; k++) { // identify BCC9 particles in sBCC9 array
        		sClu[sp4b[i][k]]=1;
        		sClu[sp4b[j][k]]=1;
      		}
      		sClu[s_com]=1;
		}
	}
	
	for (i=0; i<nsp4c-1; i++) { // loop over all sp4c_i
		for (j=i+1; j<nsp4c; ++j) { // loop over all sp4c_j

			m=0; // check exactly one common spindle particle between sp4c_i and sp4c_j
			for (k=4; k<6; k++) {
				for (l=4; l<6; l++) {
					if (sp4c[i][k]==sp4c[j][l]) {
						s_com=sp4c[i][k];
						m++;
					}
				}
			}
			if (m==0 || m>1) continue; // exactly one commond spindle
			
			flg=0;  // check number of common particles (not s_com) between sp4c clusters is 1
			for (k=0; k<6; k++) {
				if (sp4c[i][k]==s_com) continue;
				for (l=0; l<6; l++) {
					if (sp4c[j][l]==s_com) continue;
					if (sp4c[i][k]==sp4c[j][l]) {
						flg=1;
						break;
					}
				}
				if (flg==1) break;
			}
			if (flg==1) continue; // exactly one common particle (s_com)
			
			flg=0;  // check each SP4 ring particle of sp4c_i is bonded to exactly one in SP4_j
			for (k=0; k<4; k++) {
				m=0;
				for (l=0; l<4; l++) {
					if (Bonds_BondCheck(sp4c[i][k],sp4c[j][l])) {
						m++;
						if (m==2) {
							flg=1;
							break;
						}
					}
				}
				if (flg==1) break;
			}
			if (flg==1) continue; // checked
			
			flg=0;  // check each SP4 ring particle of sp4c_j is bonded to exactly one in SP4_i
			for (k=0; k<4; k++) {
				m=0;
				for (l=0; l<4; l++) {
					if (Bonds_BondCheck(sp4c[j][k],sp4c[i][l])) {
						m++;
						if (m==2) {
							flg=1;
							break;
						}
					}
				}
				if (flg==1) break;
			}
			if (flg==1) continue; // checked
			
			// now have found BCC9 cluster
      		++nClu;

      		for (k=0; k<4; k++) { // identify BCC9 particles in sBCC9 array
      			sClu[sp4c[i][k]]=1;
        		sClu[sp4c[j][k]]=1;
      		}
      		sClu[s_com]=1;
      	}
	}
	
	for (i=0; i<nsp4b; i++) { // loop over all sp4b_i
		for (j=0; j<nsp4c; ++j) { // loop over all sp4c_j
			
			m=0; // check exactly one common spindle particle between sp4b_i and sp4c_j
			for (k=4; k<6; k++) {
				if (sp4b[i][4]==sp4c[j][k]) {
					s_com=sp4b[i][4];
					m++;
				}
			}
			if (m==0 || m>1) continue; // exactly one commond spindle
			
			flg=0;  // check number of common particles (not s_com) between sp4c clusters is 1
			for (k=0; k<4; k++) {
				for (l=0; l<6; l++) {
					if (sp4c[j][l]==s_com) continue;
					if (sp4b[i][k]==sp4c[j][l]) {
						flg=1;
						break;
					}
				}
				if (flg==1) break;
			}
			if (flg==1) continue; // exactly one common particle (s_com)
			
			flg=0;  // check each SP4 ring particle of sp4b_i is bonded to exactly one in SP4_j
			for (k=0; k<4; k++) {
				m=0;
				for (l=0; l<4; l++) {
					if (Bonds_BondCheck(sp4b[i][k],sp4c[j][l])) {
						m++;
						if (m==2) {
							flg=1;
							break;
						}
					}
				}
				if (flg==1) break;
			}
			if (flg==1) continue; // checked
			
			flg=0;  // check each SP4 ring particle of sp4c_j is bonded to exactly one in SP4_i
			for (k=0; k<4; k++) {
				m=0;
				for (l=0; l<4; l++) {
					if (Bonds_BondCheck(sp4c[j][k],sp4b[i][l])) {
						m++;
						if (m==2) {
							flg=1;
							break;
						}
					}
				}
				if (flg==1) break;
			}
			if (flg==1) continue; // checked
			
			// now have found BCC9 cluster
      		++nClu;

      		for (k=0; k<4; k++) { // identify BCC9 particles in sBCC9 array
      			sClu[sp4b[i][k]]=1;
        		sClu[sp4c[j][k]]=1;
      		}
      		sClu[s_com]=1;
		}
	}
}

int tccBCC9::Find()
{
  int i;
  int total_BCC9_parts;

  for(i=0; i<N; i++) {  //  four membered shortest-path rings, then sp4b and 6A
    Rings_gSP3(i);
  }
  Clusters_GetBCC9();  // find the BCC9

  total_BCC9_parts=0;
  for(i=0; i<N; i++) {
    total_BCC9_parts+=sClu[i];
  }

  return total_BCC9_parts;
}
