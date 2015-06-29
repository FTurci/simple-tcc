
#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include "clusterinfo.h"
#include "tcc13A.h"
#include <string>


#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define BLUE    "\x1b[34m"
#define MAGENTA "\x1b[35m"
#define CYAN    "\x1b[36m"
#define RESET   "\x1b[0m"

// Declare a bunch of global variables, extracted from the many sparse headers of the TCC
int N;                  //number of particles
int nClu;               // number of clusters
std::vector<clusterinfo> clusters;
std::vector<int> sClu;  // N-length array saying if particle i is a member of a cluster
std::vector<int> sSpe;  // N-length array saying if particle i is a special member of a cluster
std::vector<int> cnb;   // N-length array storing the number of bonds 
std::vector<std::vector<int> > bNums; //NxNumBonds table of bonds for every particles

inline bool exists (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
}
//check bonds from external file
bool Bonds_BondCheck(int i, int j)
{ // Returns 1 if i & j are bonded; 0 otherwise
  int k;

  for (k=0; k<cnb[i]; ++k) {
    if (bNums[i][k] == j) return 1;
  }
  return 0;
}

tcc13A::tcc13A( int _N ) 
{

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

int sum(std::vector<int> v){
  int s=0;
  for (int i = 0; i < v.size(); ++i)
  {
    s+=v[i];
  }
  return s;
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


int main(int argc, char const *argv[])
{
  if(argc<3){
    std::cout<<"!!! Input error:\nPlease input the number of particles and the bond file name.\n";
    exit(0);
  }

  if(!exists(argv[2])){
        std::cout<<"!!! Input error:\nCannot find the bond file.\n";
    exit(1);
  }
  // set the number of particles
  N=atoi(argv[1]); 
  bNums.resize(N);
  // read the bonding information
  // the bonding file has to contain (one line per particle) :
  // first:
  //      the number of bonds
  // then:
  //      the ids of the bonded particles
  // Example: particle 0 has 8 bonds to particles (1,2,3,4,5,6,7,8) so the beginning of the file will be
  // 8 1 2 3 4 5 6 7 8
  std::cout<<RESET<<"\nReading the bonding information for "<<CYAN<<N<<RESET<<" particles from input file "<<CYAN<<argv[2]<<RESET<<" ...\n\n";
  std::ifstream fin(argv[2], std::ios::in);
  int dummy_nbonds;
  for (int i = 0; i < N; ++i)
  { 
    fin>>dummy_nbonds;
    cnb.push_back(dummy_nbonds);
    bNums[i].resize(dummy_nbonds);
    for (int k = 0; k < dummy_nbonds; ++k)
    {
      fin>>bNums[i][k];
    }
    
  }
  // Find the Icosahedra
  tcc13A icosahedra(N);
  icosahedra.Find();

  // output the results

  std::cout<<"The number of icosahedral clusters is "<<RED<<nClu<<RESET<<std::endl;
  std::cout<<"The number of particles in icosahedra is "<<RED<<sum(sClu)<<RESET<<std::endl; 
  char bufname[256];
  sprintf(bufname,"%s.icosahedra",argv[2]);
  std::cout<<"Writing the cluster file "<<RED<<bufname<<RESET<<" ...\n";


  std::ofstream fout(bufname,std::ios::out);
  for (int i = 0; i < N; ++i)
  {
    fout<<sClu[i]<<std::endl;
  }

  std::cout<<std::endl<<"Done.\n";
  return 0;
}
