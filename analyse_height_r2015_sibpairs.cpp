#include <string.h>
#include <math.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <random>

#define PACK_DENSITY 4
#define MASK0   3 // 3 << 2 * 0
#define MASK1  12 // 3 << 2 * 1
#define MASK2  48 // 3 << 2 * 2
#define MASK3 192 // 3 << 2 * 3

using namespace std;

void decode_plink(char *output, const char *input, const int lengthInput){
  int i, k;
  char tmp, geno;
  int a1, a2;
  
  for(i=0;i<lengthInput;++i){
    tmp = input[i];
    k   = PACK_DENSITY * i;
    geno      = (tmp & MASK0);
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK1) >> 2; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK2) >> 4; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK3) >> 6; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
  }
}

// Main funtion used on cov_files.txt with a sequential access
int main(int argc, char *argv[]){
  // Input arguments
  string prefix   = "ht_sample_selected_for_analysis";
  double sigma2_e = 0.5334163;
  double sigma2_f = 0.4720280;
  double theta    = sigma2_f / sigma2_e;
  
  // Indices
  int i,j,k;

  string bedfile = prefix+".bed";
  string bimfile = prefix+".bim";
  string famfile = prefix+".fam.phen";

  // Few tools
  string line = "";
  string tok  = ""; 
  ifstream tmpStream;
  
  // Get number of SNPs
  int M = -1;
  tmpStream.open(bimfile.c_str());
  while(tmpStream){
    getline(tmpStream,line);
    M++;
  }
  tmpStream.close();
  cout<<"# Found "<<M<<" SNPs."<<endl;

  // Get sample size
  int N = -1; 
  tmpStream.open(famfile.c_str());
  while(tmpStream){
    getline(tmpStream,line);
    N++;
  }
  tmpStream.close();  
  cout<<"# Found "<<N<<" samples."<<endl;
  
  int    *FID   = new int[N];
  int    *IID   = new int[N];
  double *Y_tmp = new double[N];
  string fid, iid, pid, mid, sex, phn; // here - if fam file has changed dimensions
  tmpStream.open(famfile.c_str());
  for(i=0;i<N;i++){
    tmpStream >> fid;
    tmpStream >> iid;
    tmpStream >> pid;
    tmpStream >> mid;
    tmpStream >> sex;
    tmpStream >> phn;
    
    FID[i]   = atoi(fid.c_str());
    IID[i]   = atoi(iid.c_str());
    Y_tmp[i] = atof(phn.c_str());
  }
  tmpStream.close();
  
  int N_F = 8348; // here - maybe as an input parameter
  int *n  = new int[N_F];
  for(k=0;k<N_F;k++) n[k] = 0;
  for(i=0;i<N;i++){
    n[FID[i]]++;
  }
  
  double **X = new double*[N_F];
  double **Y = new double*[N_F];
  for(k=0;k<N_F;k++){
    X[k] = new double[n[k]];
    Y[k] = new double[n[k]];
  }

  for(i=0;i<N;i++){
    Y[FID[i]][IID[i]] = Y_tmp[i];
  }
  
  string *CHROM   = new string[M];
  string *SNPS    = new string[M];
  string *POS     = new string[M];
  string *A1      = new string[M];
  string *A2      = new string[M];

  tmpStream.open(bimfile.c_str());
  for(j=0;j<M;j++){
    tmpStream >> CHROM[j];
    tmpStream >> SNPS[j];
    tmpStream >> tok;
    tmpStream >> POS[j];
    tmpStream >> A1[j];
    tmpStream >> A2[j];
  }
  tmpStream.close();

  int numBytes   = (int)ceil((double)N / PACK_DENSITY);
  char* packed   = new char[numBytes];
  char* unpacked = new char[numBytes * PACK_DENSITY];

  // run GWAS
  ifstream influx;
  influx.open(bedfile.c_str(), std::ios::in | std::ios::binary);
  if(!influx){
    cerr << "[readGenotypes] Error reading file "<<bedfile<<endl;
    exit(1);
  }
  influx.seekg(0, ifstream::end);
  influx.seekg(3, ifstream::beg);
  
  double x;
  double A_11, A_12, A_22, B_11, B_21;
  double s_x, s_x2, s_y, s_xy;
  double n_k, s_k;
  double Determinant;
  double BETA_REML, SE_REML;
  double BETA_OLS,  SE_OLS;
  double numOLS;
  double denOLS;
  double RSS, e_ik;
  
  int *nEff   = new int[N_F];
  double *mx  = new double[N_F];
  double *my  = new double[N_F];
  int nEff_sample, nEff_family, nEff_sibpairs;
  
  string outfile = prefix+".gwas.assoc";
  ofstream fileOut(outfile.c_str());
  fileOut<<"CHROM\tSNP\tBP\tA1\tA2\tN_FAM\tN_SAMPLES\tN_SIBPAIRS\tBETA_REML\tSE_REML\tBETA_OLS\tSE_OLS\n";
  
  cout<<"#   0%[          ]\n"; //  0%
  bool completed10  = false;
  bool completed20  = false;
  bool completed30  = false;
  bool completed40  = false;
  bool completed50  = false;
  bool completed60  = false;
  bool completed70  = false;
  bool completed80  = false;
  bool completed90  = false;
  bool completed100 = false;
  
  int dp = M / 10;

  for(j=0;j<M;j++){
    if(j>=1*dp-1 and !completed10){
      completed10 = true;
      cout<<"#  10%[=         ]\n"; // 10%
    }
    if(j>=2*dp-1 and !completed20){
      completed20 = true;
      cout<<"#  20%[==        ]\n"; // 20%
    }
    if(j>=3*dp-1 and !completed30){
      completed30 = true;
      cout<<"#  30%[===       ]\n"; // 30%
    }
    if(j>=4*dp-1 and !completed40){
      completed40 = true;
      cout<<"#  40%[====      ]\n"; // 40%
    }
    if(j>=5*dp-1 and !completed50){
      completed50 = true;
      cout<<"#  50%[=====     ]\n"; // 50%
    }
    if(j>=6*dp-1 and !completed60){
      completed60 = true;
      cout<<"#  60%[======    ]\n"; // 60%
    }
    if(j>=7*dp-1 and !completed70){
      completed70 = true;
      cout<<"#  70%[=======   ]\n"; // 70%
    }
    if(j>=8*dp-1 and !completed80){
      completed80 = true;
      cout<<"#  80%[========  ]\n"; // 80%
    }
    if(j>=9*dp-1 and !completed90){
      completed90 = true;
      cout<<"#  90%[========= ]\n"; // 90%
    }
    if(j>=10*dp-1 and !completed100){
      completed100 = true;
      cout<<"# 100%[==========]\n"; //100%
    }
    
    // Initialize
    for(k=0;k<N_F;k++){
      mx[k]   = 0.;
      my[k]   = 0.;
      nEff[k] =  0;
    }
    
    influx.read((char*)packed, sizeof(char) * numBytes);
    decode_plink(unpacked, packed, numBytes);
    for(i=0;i<N;i++){
      x = (double) ((int) unpacked[i]);
      X[FID[i]][IID[i]] = x;
      if(x!=3.){
        nEff[FID[i]]++;
        mx[FID[i]] += x;
        my[FID[i]] += Y[FID[i]][IID[i]];
      }
    }
    
    // Scale mx
    for(k=0;k<N_F;k++){
      mx[k]   = mx[k] / nEff[k];
      my[k]   = my[k] / nEff[k];
    }
    
    for(i=0;i<N;i++){
      if(X[FID[i]][IID[i]]==3.){
        X[FID[i]][IID[i]] = mx[FID[i]];
      }
    }
    
    A_11 = 0.;
    A_12 = 0.;
    A_22 = 0.;
    B_11 = 0.;
    B_21 = 0.;
    
    nEff_sample   = 0;
    nEff_family   = 0;
    nEff_sibpairs = 0;
    numOLS = 0.;
    denOLS = 0.;
    for(k=0;k<N_F;k++){
      if(nEff[k]>=2){
        nEff_family++;
        nEff_sample += nEff[k];
        nEff_sibpairs += nEff[k]*(nEff[k]-1)/2;
        s_x   = 0.;
        s_x2  = 0.;
        s_y   = 0.;
        s_xy  = 0.;
        for(i=0;i<nEff[k];i++){
          s_x  += X[k][i];
          s_y  += Y[k][i];
          s_xy += X[k][i] * Y[k][i];
          s_x2 += X[k][i] * X[k][i];
          
          // OLS
          numOLS += (X[k][i]-mx[k])*(Y[k][i]-my[k]);
          denOLS += (X[k][i]-mx[k])*(X[k][i]-mx[k]);
        }
        // Update A's and B's
        n_k   = (double) nEff[k];
        s_k   = 1.0 / (1.0 + n_k * theta);
        A_11 += n_k * s_k;
        A_12 += s_x * s_k;
        A_22 += s_x2 - theta * s_x * s_x * s_k;
        B_11 += s_y * s_k;
        B_21 += s_xy - theta * s_x * s_y * s_k;
      }
    }

    // REML 
    Determinant = A_11 * A_22 - A_12 * A_12;
    BETA_REML   = (A_11*B_21-A_12*B_11) / Determinant;
    SE_REML     = sqrt( sigma2_e * A_11 / Determinant );
    
    // OLS
    BETA_OLS = numOLS / denOLS;
    RSS = 0.;
    for(k=0;k<N_F;k++){
      if(nEff[k]>=2){
        for(i=0;i<nEff[k];i++){
          e_ik = (Y[k][i]-my[k])-BETA_OLS*(X[k][i]-mx[k]);
          RSS += e_ik * e_ik;
        }
      }
    }
    SE_OLS = sqrt( (RSS/(nEff_sample-nEff_family)) / denOLS );

    // Write Results
    fileOut<<CHROM[j]<<"\t"<<SNPS[j]<<"\t"<<POS[j]<<"\t"<<A1[j]<<"\t"<<A2[j]<<"\t"<<nEff_family<<"\t"<<nEff_sample<<"\t"<<nEff_sibpairs<<"\t"<<BETA_REML<<"\t"<<SE_REML<<"\t"<<BETA_OLS<<"\t"<<SE_OLS<<"\n";
  }
  influx.close();
  fileOut.close();

  delete [] packed;
  delete [] unpacked;    
  delete [] CHROM;
  delete [] SNPS;
  delete [] POS;
  delete [] A1;
  delete [] A2;
  delete [] Y;
  delete [] X;
  delete [] FID;
  delete [] IID;
  delete [] Y_tmp;
  delete [] mx;
  delete [] my;
  delete [] nEff;

  return EXIT_SUCCESS;
}


