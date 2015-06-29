
#include "global.h"
#include "supercool.h"
#include "tcc11W.h"

tcc11W::tcc11W( int _N ) : tcc_base( _N )
{
  try { fc = g_ctx.props.get<double>( "tcc.fc_11W" ); }
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

void tcc11W::Rings_gSP3(int n0)
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

void tcc11W::Rings_gSP4(int n0, int n1, int n2)
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

void tcc11W::Rings_gSP5(int n0, int n1, int n2, int n3) 
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

void tcc11W::Rings_aSP5(int n0, int n1, int n2, int n3, int n4)
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

void tcc11W::Clusters_Get11W_Cs()
{ 
  int i, j, k, l, m, n;
  int sp1, sp2i, sp2j, sp5com[2], ep;
  int trial9B[9], trial11W[11];
  int flg, flg1, flg2, break_out;

  for (i=0; i<N; ++i) sClu[i] = 0;
  
  for (i=0; i<nsp5c; ++i) { //first 7A
    for (j=i+1; j<nsp5c; ++j) { // second 7A
      flg = 0;
      if (sp5c[i][5] == sp5c[j][5] && sp5c[i][6] != sp5c[j][6]) {  // check one common spindle
        if (Bonds_BondCheck(sp5c[i][6], sp5c[j][6])) {  // check uncommon spindles are bonded
          flg = 1;
          sp1 = sp5c[i][5];  // s_com common spindle
          sp2i = sp5c[i][6];  // 2nd spindle particle of cluster 7A_i
          sp2j = sp5c[j][6];  // 2nd spindle particle of cluster 7A_j
        }
      }
      if (sp5c[i][6] == sp5c[j][6] && sp5c[i][5] != sp5c[j][5]) {  // check one common spindle
        if (Bonds_BondCheck(sp5c[i][5], sp5c[j][5])) {  // check uncommon spindles are bonded
          flg = 1;
          sp1 = sp5c[i][6];  // s_com common spindle
          sp2i = sp5c[i][5];  // 2nd spindle particle of cluster 7A_i
          sp2j = sp5c[j][5];  // 2nd spindle particle of cluster 7A_j
        }
      }
      if (sp5c[i][5] == sp5c[j][6] && sp5c[i][6] != sp5c[j][5]) {  // check one common spindle
        if (Bonds_BondCheck(sp5c[i][6], sp5c[j][5])) {  // check uncommon spindles are bonded
          flg = 1;
          sp1 = sp5c[i][5];  // s_com common spindle
          sp2i = sp5c[i][6];  // 2nd spindle particle of cluster 7A_i
          sp2j = sp5c[j][5];  // 2nd spindle particle of cluster 7A_j
        }          
      }
      if (sp5c[i][6] == sp5c[j][5] && sp5c[i][5] != sp5c[j][6]) {  // check one common spindle
        if (Bonds_BondCheck(sp5c[i][5], sp5c[j][6])) {  // check uncommon spindles are bonded
          flg = 1;
          sp1 = sp5c[i][6];  // s_com common spindle
          sp2i = sp5c[i][5];   // 2nd spindle particle of cluster 7A_i
          sp2j = sp5c[j][6];  // 2nd spindle particle of cluster 7A_j
        }          
      }
      if (flg==0) continue;  // one common spindle, other spindles are bonded
       
      flg1 = flg2 = 0;  // ensure the two distinct spindle particles are part of the SP5 ring of other sp5c cluster
      for (k=0; k<5; ++k) {
        if (sp2i == sp5c[j][k]) flg1 = 1;
        if (sp2j == sp5c[i][k]) flg2 = 1;
      }
      if (flg1==0 || flg2==0) continue;
      
      m = 0;  // check for two common particles between SP5 rings of sp5c_i and sp5c_j
      for (k=0; k<5; ++k) {
        for (l=0; l<5; ++l) {
          if (sp5c[i][k] == sp5c[j][l]) {
            if (m==2) {
              m++; 
              break; 
            }
            sp5com[m]=sp5c[i][k];
            ++m;
          }
        }
      }
      if (m!=2) continue;
      
      if (sp5com[0]<sp5com[1]) {
        trial9B[4]=sp5com[0];
        trial9B[5]=sp5com[1];
      }
      else {
        trial9B[4]=sp5com[1];
        trial9B[5]=sp5com[0];
      }
      
      if (sp2i<sp2j) {
        trial9B[6]=sp2i;
        trial9B[7]=sp2j;

        for (k=0; k<5; ++k) {
          if (Bonds_BondCheck(sp5c[i][k],trial9B[4]) && sp5c[i][k]!=trial9B[7] && sp5c[i][k]!=trial9B[4]) {
            trial9B[0]=sp5c[i][k];
          }
          if (Bonds_BondCheck(sp5c[i][k],trial9B[5]) && sp5c[i][k]!=trial9B[7] && sp5c[i][k]!=trial9B[5]) {
            trial9B[1]=sp5c[i][k];
          }
          if (Bonds_BondCheck(sp5c[j][k],trial9B[4]) && sp5c[j][k]!=trial9B[6] && sp5c[j][k]!=trial9B[4]) {
            trial9B[2]=sp5c[j][k];
          }
          if (Bonds_BondCheck(sp5c[j][k],trial9B[5]) && sp5c[j][k]!=trial9B[6] && sp5c[j][k]!=trial9B[5]) {
            trial9B[3]=sp5c[j][k];
          }
        }
      }
      else {
        trial9B[6]=sp2j;
        trial9B[7]=sp2i;
        
        for (k=0; k<5; ++k) {
          if (Bonds_BondCheck(sp5c[j][k],trial9B[4]) && sp5c[j][k]!=trial9B[7] && sp5c[j][k]!=trial9B[4]) {
            trial9B[0]=sp5c[j][k];
          }
          if (Bonds_BondCheck(sp5c[j][k],trial9B[5]) && sp5c[j][k]!=trial9B[7] && sp5c[j][k]!=trial9B[5]) {
            trial9B[1]=sp5c[j][k];
          }
          if (Bonds_BondCheck(sp5c[i][k],trial9B[4]) && sp5c[i][k]!=trial9B[6] && sp5c[i][k]!=trial9B[4]) {
            trial9B[2]=sp5c[i][k];
          }
          if (Bonds_BondCheck(sp5c[i][k],trial9B[5]) && sp5c[i][k]!=trial9B[6] && sp5c[i][k]!=trial9B[5]) {
            trial9B[3]=sp5c[i][k];
          }
        }
      }
      trial9B[8]=sp1;
      
      // Now we have found the 9B C2v cluster need to check for 10B C3v
      for (k=j+1; k<nsp5c; ++k) {
        if (sp5c[k][5] == sp1) {  // check one spindle of sp5c_k is the common spindle of 9B
          if (Bonds_BondCheck(sp5c[k][6], sp2i)==0) continue;  // check other spindle of sp5c_k is bonded to spindle sp2i of 9B
          if (Bonds_BondCheck(sp5c[k][6], sp2j)==0) continue;  // check other spindle of sp5c_k is bonded to spindle sp2j of 9B
          
          flg1=0;  // check sp2i and sp2j are in SP5 ring of sp5c_k
          flg2=0;
          for (l=0;l<5;l++) {
            if (sp5c[k][l]==sp2i) {
              flg1=1;
              continue;
            }
            if (sp5c[k][l]==sp2j) {
              flg2=1;
              continue;
            }
          }
          if (flg1==0 || flg2==0) continue;
          trial11W[6]=trial9B[6];
          trial11W[7]=trial9B[7];
          trial11W[8]=sp5c[k][6];
          trial11W[9]=trial9B[8];
          
          m=0;  // check other spindle of sp5c_k is in non-spindle SP5 ring particle in 9B cluster
          break_out=0;
          for (l=0;l<6;l++) {
            if (trial9B[l]==sp5c[k][6]) continue;
            if (m==5) {
              m++;
              break_out=1;
              break;
            }
            trial11W[m]=trial9B[l];
            m++;
          }
          if (break_out==1 || m!=5) continue;
          
          break_out=0;  // check sp2i and sp2j are in SP5 ring of sp5c_k, and that one particle in sp5c_k is distinct from 9B
          for (l=0;l<5;l++) {
            if (sp5c[k][l]==trial9B[6]) continue;
            if (sp5c[k][l]==trial9B[7]) continue;
            for (m=0;m<5;m++) {
              if (sp5c[k][l]==trial11W[m]) break;
            }
            if (m==5) {
              trial11W[5]=sp5c[k][l];
              break_out++;
            }
          }
          if (break_out!=1) continue;
          
          // Now we have found the 10B C3v cluster, need to check if also 11W
          
          if(cnb[sp1]!= 10) continue;        // sp1 has 10 bonds in total (all forming the shell)

          n = 0;        // find extra particle
          break_out=0;
          for (l=0; l<10; ++l) {
                  for (m=0; m<9; ++m) {
                          if(bNums[sp1][l] == trial11W[m]) break;
                  }
                  if(m==9){
                          if(n==1) { // two or more particles
                                  break_out=1;
                                  break;
                          }        
                          ep= bNums[sp1][l];
                          n++;
                  } 
          }
          if (break_out==1 || n<1) continue;
          
          if (Bonds_BondCheck(ep, trial11W[6])==1) continue;        // extra particles must not be bonded to three 7A spindles in shell of 10B
          if (Bonds_BondCheck(ep, trial11W[7])==1) continue;        // extra particles must not be bonded to three 7A spindles in shell of 10B
          if (Bonds_BondCheck(ep, trial11W[8])==1) continue;        // extra particles must not be bonded to three 7A spindles in shell of 10B
          
          // Now we have found the 11W Cs cluster
          ++nClu;
          
          sClu[trial11W[0]]=1;
          sClu[trial11W[1]]=1;
          sClu[trial11W[2]]=1;
          sClu[trial11W[3]]=1;
          sClu[trial11W[4]]=1;
          sClu[trial11W[5]]=1;
          sClu[trial11W[6]]=1;
          sClu[trial11W[7]]=1;
          sClu[trial11W[8]]=1;
          sClu[trial11W[9]]=1;
	  sClu[ep]=1;
      }
      if (sp5c[k][6] == sp1) {  // check one spindle of sp5c_k is the common spindle of 9B
          if (Bonds_BondCheck(sp5c[k][5], sp2i)==0) continue;  // check other spindle of sp5c_k is bonded to spindle sp2i of 9B
          if (Bonds_BondCheck(sp5c[k][5], sp2j)==0) continue;  // check other spindle of sp5c_k is bonded to spindle sp2j of 9B
          
          flg1=0;  // check sp2i and sp2j are in SP5 ring of sp5c_k
          flg2=0;
          for (l=0;l<5;l++) {
            if (sp5c[k][l]==sp2i) {
              flg1=1;
              continue;
            }
            if (sp5c[k][l]==sp2j) {
              flg2=1;
              continue;
            }
          }
          if (flg1==0 || flg2==0) continue;
          trial11W[6]=trial9B[6];
          trial11W[7]=trial9B[7];
          trial11W[8]=sp5c[k][5];
          trial11W[9]=trial9B[8];
          
          m=0;  // check other spindle of sp5c_k is in non-spindle SP5 ring particle in 9B cluster
          break_out=0;
          for (l=0;l<6;l++) {
            if (trial9B[l]==sp5c[k][5]) continue;
            if (m==5) {
              m++;
              break_out=1;
              break;
            }
            trial11W[m]=trial9B[l];
            m++;
          }
          if (break_out==1 || m!=5) continue;
          
          break_out=0;  // check sp2i and sp2j are in SP5 ring of sp5c_k, and that one particle in sp5c_k is distinct from 9B
          for (l=0;l<5;l++) {
            if (sp5c[k][l]==trial9B[6]) continue;
            if (sp5c[k][l]==trial9B[7]) continue;
            for (m=0;m<5;m++) {
              if (sp5c[k][l]==trial11W[m]) break;
            }
            if (m==5) {
              trial11W[5]=sp5c[k][l];
              break_out++;
            }
          }
          if (break_out!=1) continue;
          
          // Now we have found the 10B C3v cluster, need to check if also 11W
          
          if(cnb[sp1]!= 10) continue;        // sp1 has 10 bonds in total (all forming the shell)

          n = 0;        // find extra particle
          break_out=0;
          for (l=0; l<10; ++l) {
                  for (m=0; m<9; ++m) {
                          if(bNums[sp1][l] == trial11W[m]) break;
                  }
                  if(m==9){
                          if(n==1) { // two or more particles
                                  break_out=1;
                                  break;
                          }        
                          ep= bNums[sp1][l];
                          n++;
                  } 
          }
          if (break_out==1 || n<1) continue;
          
          if (Bonds_BondCheck(ep, trial11W[6])==1) continue;        // extra particles must not be bonded to three 7A spindles in shell of 10B
          if (Bonds_BondCheck(ep, trial11W[7])==1) continue;        // extra particles must not be bonded to three 7A spindles in shell of 10B
          if (Bonds_BondCheck(ep, trial11W[8])==1) continue;        // extra particles must not be bonded to three 7A spindles in shell of 10B
          
          // Now we have found the 11W Cs cluster
          ++nClu;
          
          sClu[trial11W[0]]=1;
          sClu[trial11W[1]]=1;
          sClu[trial11W[2]]=1;
          sClu[trial11W[3]]=1;
          sClu[trial11W[4]]=1;
          sClu[trial11W[5]]=1;
          sClu[trial11W[6]]=1;
          sClu[trial11W[7]]=1;
          sClu[trial11W[8]]=1;
          sClu[trial11W[9]]=1;
	  sClu[ep]=1;
        }
      }
    }
  }
}

int tcc11W::Find()
{
  int i;
  int total_11W_parts;

  for(i=0; i<N; i++) {  //  five membered shortest-path rings, then 7A
    Rings_gSP3(i);
  }
  Clusters_Get11W_Cs();  // find the 11W

  total_11W_parts=0;
  for(i=0; i<N; i++) {
    total_11W_parts+=sClu[i];
  }

  return total_11W_parts;
}
