/*
  pf.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 3/2006, Justin Bois 1/2007
            Asif Khan 8/2009 Brian Wolfe 10/2009
  This file moves the partition function algorithm to a function, so
  that it can be more readily used as in a library.  In addition,
  scaling is incorporated to allow for longer sequences.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "pfuncUtilsHeader.h" //contains functions and structures
#include "DNAExternals.h"


/* ************************************************ */
// This is the main function for computing partition functions.

DEV
int getNodesFirstEntry(int L, int rank, int seqLen, int N) {
    int split, rem, last;
      last = seqLen - L + 1;
        if (rank >= N || rank < 0)
              return last; 
          split = (int)(last/N);
            rem = last % N;
              return rank*split + MIN(rank, rem);
}

DEV
int IsProcessorInUse(int rank, int L, int seqLen, int numCPUs ) {
    return (getNodesFirstEntry(L,rank,seqLen,numCPUs) != seqLen - L + 1);
}

__constant__ energy_model_t *energies;

GLB
void pfuncFullWithSymHelper(DBL_TYPE *pf, int ** inputSeqs, int * seqlengths,
    int * nStrands_arr, int * permSymmetries) {

  int * inputSeq = inputSeqs[blockIdx.x];
  int seqlength = seqlengths[blockIdx.x];
  int nStrands = nStrands_arr[blockIdx.x];
  int permSymmetry = permSymmetries[blockIdx.x];
  //complexity: 3 = N^3, 4 = N^4, 5 = N^5, 8 = N^8
  //naType: DNA = 0, RNA = 1
  //dangles: 0 = none, 1 = normal, 2 = add both

  __shared__ int *seq;
  if(threadIdx.x == 0) {
    seq = (int*)malloc((seqlength + 1) * sizeof(int));
  }

  __shared__ DBL_TYPE *Q;
  __shared__ DBL_TYPE *Qb;
  __shared__ DBL_TYPE *Qm; //O(N^2)


  //N^3 arrays
  __shared__ DBL_TYPE *Qx;
  __shared__ DBL_TYPE *Qx_1;
  __shared__ DBL_TYPE *Qx_2;
  __shared__ DBL_TYPE *Qs;
  __shared__ DBL_TYPE *Qms;

  /*
  The above matrices are dynamically allocated matrices that
  contain partition functions restricted to a subsequence of the
  strand.  Each of the above should be accessed by the call
  to Q[ pf_index(i, j)] indicate the partition function between
  i and j, inclusive. 

  They are described in the paper mentioned above.
  */

  int i, j; // the beginning and end bases for Q;
  int L; //This the length of the current subsequence 
  int pf_ij; //index for O(N^2) matrixes; used to reduce calls to pf_index;
  DBL_TYPE returnValue;


  int iMin;
  int iMax;
  

  __shared__ int nicks[ MAXSTRANDS];  //the entries must be strictly increasing
  if (threadIdx.x == 0) {
    for (i=0;i<MAXSTRANDS;i++){
      nicks[i]=-1;
    }
  }
  //nicks[i] = N means a strand ends with base N, and a new one starts at N+1
  // isNicked[n] is 0 if no nick at n, 1 otherwise

  __shared__ int **etaN;
  __shared__ int *etaN_space;

  if (threadIdx.x == 0) {
    processMultiSequence( inputSeq, seqlength, nStrands, seq, nicks);

    // Allocate and Initialize Matrices
    int arraySize = seqlength*(seqlength+1)/2+(seqlength+1);

    InitLDoublesMatrix( &Q, arraySize, "Q");
    InitLDoublesMatrix( &Qb, arraySize, "Qb");
    InitLDoublesMatrix( &Qm, arraySize, "Qm");

    etaN = (int**)malloc(arraySize * sizeof(int*));
    etaN_space = (int*)malloc(arraySize * 2 * sizeof(int));
    for (int i = 0; i < arraySize; ++i) {
      etaN[i] = etaN_space + (2 * i);
    }
    InitEtaN( etaN, nicks, seqlength);
    nonZeroInit( Q, seq, seqlength, energies);

    InitLDoublesMatrix( &Qs, arraySize, "Qs");
    InitLDoublesMatrix( &Qms, arraySize, "Qms");

    InitLDoublesMatrix( &Qx, arraySize/2, "Qx");
    InitLDoublesMatrix( &Qx_1, arraySize/2, "Qx_1");
    InitLDoublesMatrix( &Qx_2, arraySize/2, "Qx_2");

  }

  for( L = 1; L <= seqlength; L++) {
    /* Calculate all sub partition functions for
    distance = 0, then 1, then 2.... */

    if (threadIdx.x == 0) {
      manageQx( &Qx, &Qx_1, &Qx_2, L-1, seqlength);   
    }
    __syncthreads();

    iMin = getNodesFirstEntry(L, threadIdx.x, seqlength, blockDim.x);
    iMax = getNodesFirstEntry(L, threadIdx.x + 1, seqlength, blockDim.x) - 1;
    /*
    //Default without parallelization
    iMin = 0;
    iMax = seqlength - L; 
    */

  
    int active = IsProcessorInUse(threadIdx.x, L, seqlength, blockDim.x);
    for( i = iMin; i <= iMax && active; i++) {
      j = i + L - 1;
      pf_ij = pf_index( i, j, seqlength);
      /* Recursions for Qb.  See figure 13 of paper */
      /* bp = base pairs, pk = pseudoknots */
      if( CanPair( seq[ i], seq[ j]) == FALSE) {
        Qb[ pf_ij] = 0.0; //scaling still gives 0
      }
      else {
        Qb[ pf_ij] = ExplHairpin( i, j, seq, seqlength, etaN, energies);

        //no nicked haripins allowed in previous function
        if( etaN[ EtaNIndex(i+0.5, i+0.5, seqlength)][0] == 0 &&
           etaN[ EtaNIndex(j-0.5, j-0.5, seqlength)][0] == 0) {
             //regular multiloop.  No top-level nicks
             
              Qb[ pf_ij] += SumExpMultiloops(i, j, seq, Qms, Qm,
                                            seqlength, etaN, energies);
        }

        if( etaN[ EtaNIndex(i+0.5, j-0.5, seqlength)][0] >= 1) {
          //Exterior loop (created by nick)
          Qb[ pf_ij] += SumExpExteriorLoop( i, j, seq, seqlength, 
                                           Q, nicks, etaN, energies); 
        }

      }

      fastILoops( i, j, L, seqlength, seq, etaN, Qb, Qx, Qx_2, energies);


      /* Recursions for Qms, Qs */
      MakeQs_Qms( i, j, seq, seqlength, Qs, Qms, Qb, nicks, etaN, energies);
      
      /* Recursions for Q, Qm, Qz */
      MakeQ_Qm_N3( i, j, seq, seqlength, Q, Qs, Qms, Qm,
                  nicks,etaN, energies);


    }
    __syncthreads();
  }

  //adjust this for nStrands, symmetry at rank == 0 node
  if (threadIdx.x == 0) {
    returnValue = EXP_FUNC(
        -1*(energies->bimolecular + energies->salt_correction)*(nStrands-1)/
        (kB*energies->temp_k)
      ) * Q[ pf_index(0,seqlength-1, seqlength)]/((DBL_TYPE) permSymmetry);

    free( Q);
    free( Qb);
    free( Qm);

    Q = Qb = Qm = NULL;

    free( Qs);
    free( Qms);
    
    free( Qx);
    free( Qx_1);
    free( Qx_2);
    
    Qs = Qms = Qx = Qx_1 = Qx_2 = NULL;

    free( seq);

    for( i = 0; i <= seqlength-1; i++) {
      for( j = i-1; j <= seqlength-1; j++) {
        pf_ij = pf_index(i,j,seqlength);
        free( etaN[pf_ij]);
      }
    }

    free( etaN);
    pf[blockIdx.x] = returnValue;
  }
}
/* ****** */
//void pfuncInitialize(DBL_TYPE temp_lo, DBL_TYPE temp_hi, DBL_TYPE temp_step,

DBL_TYPE * pf;
int ** intSeqs;
int * seqlengths;
int * nStrands_arr;
int * permSymmetries;
int nblocks;

#define cudaCheck(call) {\
  cudaError_t e = (call); \
  if(e != cudaSuccess) {\
    fprintf(stderr, #call ": %s", cudaGetErrorString(e)); \
    abort(); \
  }\
}\


extern "C" void pfuncInitialize(int nblocks_in, DBL_TYPE temp_k,
    DBL_TYPE sodium_conc, DBL_TYPE magnesium_conc,
    int long_helix, int dangletype, int dnarnacount) {

  nblocks = nblocks_in;
  SODIUM_CONC = sodium_conc;
  MAGNESIUM_CONC = magnesium_conc;
  USE_LONG_HELIX_FOR_SALT_CORRECTION = long_helix;
  DANGLETYPE = dangletype;
  DNARNACOUNT = dnarnacount;

  // Load Energy Model
  energy_model_t em;
  LoadEnergies(&em, temp_k);

  // Allocate Device Memory
  energy_model_t *em_d;
  cudaMalloc(&em_d, sizeof(energy_model_t));
  cudaMemcpy(em_d, &em, sizeof(energy_model_t), cudaMemcpyHostToDevice);

  // Send Device Ptr to Device
  cudaMemcpyToSymbol(energies, &em_d, sizeof(energy_model_t*));

  //cudaMemcpyToSymbol(ENERGIES, &energies, sizeof(energy_model_t));

  cudaCheck(
    cudaMallocManaged(&pf, nblocks * sizeof(DBL_TYPE))
  );
  cudaCheck(
    cudaMallocManaged(&seqlengths, nblocks * sizeof(int))
  );
  cudaCheck(
    cudaMallocManaged(&nStrands_arr, nblocks * sizeof(int))
  );
  cudaCheck(
    cudaMallocManaged(&permSymmetries, nblocks * sizeof(int))
  );
  cudaCheck(
    cudaMallocManaged(&intSeqs, nblocks * sizeof(int*))
  );

  for(int i = 0; i < nblocks; ++i) {
    cudaCheck(
      cudaMallocManaged(&(intSeqs[i]), (MAXSEQLENGTH + 1) * sizeof(int))
    );
  }
  cudaCheck(
    cudaDeviceSetLimit(cudaLimitMallocHeapSize, 1 << 30)
  );
}

extern "C" void pfuncMulti(char ** inputSeqs, int nseqs, int * permSym, DBL_TYPE
    * result) {

  for (int s = 0; s < nseqs; s += nblocks) {
    for(int i = 0; i < nblocks; ++i) {
      int len = strlen(inputSeqs[s + i]);
      convertSeq(inputSeqs[s + i], intSeqs[i], len);
      seqlengths[i] = getSequenceLengthInt(intSeqs[i], nStrands_arr + i);
      permSymmetries[i] = permSym[s + i];
    }

    pfuncFullWithSymHelper<<<nblocks, 256>>>(
        pf, intSeqs, seqlengths, nStrands_arr, permSymmetries
    );
    cudaDeviceSynchronize();

    for(int i = 0; i < nblocks; ++i) {
      result[s + i] = pf[i];
    }
  }
}
