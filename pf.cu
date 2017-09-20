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
#include <assert.h>

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
__constant__ DBL_TYPE t_lo;
__constant__ DBL_TYPE t_hi;
__constant__ DBL_TYPE t_step;
DEV
energy_model_t* get_energy_model(DBL_TYPE temp_k) {
  assert(t_lo <= temp_k && temp_k <= t_hi);
  int index = round((temp_k - t_lo) / t_step);
  return energies + index;
}

GLB
void pfuncFullWithSymHelper(DBL_TYPE *pf, int ** inputSeqs, int * seqlengths,
    int * nStrands_arr, int * permSymmetries, DBL_TYPE * temps) {

  int * inputSeq = inputSeqs[blockIdx.x];
  int seqlength = seqlengths[blockIdx.x];
  int nStrands = nStrands_arr[blockIdx.x];
  int permSymmetry = permSymmetries[blockIdx.x];
  energy_model_t *em = get_energy_model(temps[blockIdx.x]);

  //complexity: 3 = N^3, 4 = N^4, 5 = N^5, 8 = N^8
  //naType: DNA = 0, RNA = 1
  //dangles: 0 = none, 1 = normal, 2 = add both

  __shared__ int *seq;

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

  __shared__ int nicks[ MAXSTRANDS];
  __shared__ int **etaN;
  __shared__ int *etaN_space;

  int arraySize = seqlength*(seqlength+1)/2+(seqlength+1);

  assert(blockDim.x >= 11);
  switch (threadIdx.x) {
    case 0: etaN = (int**) malloc(arraySize * sizeof(int*));         break;
    case 1: etaN_space = (int*) malloc(arraySize * 2 * sizeof(int)); break;
    case 2: InitLDoublesMatrix( &Q, arraySize, "Q");                 break;
    case 3: InitLDoublesMatrix( &Qb, arraySize, "Qb");               break;
    case 4: InitLDoublesMatrix( &Qm, arraySize, "Qm");               break;
    case 5: InitLDoublesMatrix( &Qs, arraySize, "Qs");               break;
    case 6: InitLDoublesMatrix( &Qms, arraySize, "Qms");             break;
    case 7: InitLDoublesMatrix( &Qx, arraySize/2, "Qx");             break;
    case 8: InitLDoublesMatrix( &Qx_1, arraySize/2, "Qx_1");         break;
    case 9: InitLDoublesMatrix( &Qx_2, arraySize/2, "Qx_2");         break;
    case 10: seq = (int*) malloc((seqlength + 1) * sizeof(int));     break;
  }
  for (int i = threadIdx.x; i < MAXSTRANDS; i += blockDim.x) {
    nicks[i]=-1;
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    processMultiSequence( inputSeq, seqlength, nStrands, seq, nicks);
    nonZeroInit( Q, seq, seqlength, em);
  }
  for (int i = threadIdx.x; i < arraySize; i += blockDim.x) {
    etaN[i] = etaN_space + (2 * i);
  }
  __syncthreads();

  if (threadIdx.x == 0) {
    InitEtaN( etaN, nicks, seqlength);
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
        Qb[ pf_ij] = ExplHairpin( i, j, seq, seqlength, etaN, em);

        //no nicked haripins allowed in previous function
        if( etaN[ EtaNIndex(i+0.5, i+0.5, seqlength)][0] == 0 &&
           etaN[ EtaNIndex(j-0.5, j-0.5, seqlength)][0] == 0) {
             //regular multiloop.  No top-level nicks
             
              Qb[ pf_ij] += SumExpMultiloops(i, j, seq, Qms, Qm,
                                            seqlength, etaN, em);
        }

        if( etaN[ EtaNIndex(i+0.5, j-0.5, seqlength)][0] >= 1) {
          //Exterior loop (created by nick)
          Qb[ pf_ij] += SumExpExteriorLoop( i, j, seq, seqlength, 
                                           Q, nicks, etaN, em); 
        }

      }

      fastILoops( i, j, L, seqlength, seq, etaN, Qb, Qx, Qx_2, em);


      /* Recursions for Qms, Qs */
      MakeQs_Qms( i, j, seq, seqlength, Qs, Qms, Qb, nicks, etaN, em);
      
      /* Recursions for Q, Qm, Qz */
      MakeQ_Qm_N3( i, j, seq, seqlength, Q, Qs, Qms, Qm,
                  nicks,etaN, em);


    }
    __syncthreads();
  }

  //adjust this for nStrands, symmetry at rank == 0 node
  if (threadIdx.x == 0) {
    returnValue = EXP_FUNC(
        -1*(em->bimolecular + em->salt_correction)*(nStrands-1)/
        (kB*em->temp_k)
      ) * Q[ pf_index(0,seqlength-1, seqlength)]/((DBL_TYPE) permSymmetry);
    pf[blockIdx.x] = returnValue;
  }
  switch(threadIdx.x) {
    case 0: free( Q); break;
    case 1: free( Qb); break;
    case 2: free( Qm); break;
    case 3: free( Qs); break;
    case 4: free( Qms); break;
    case 5: free( Qx); break;
    case 6: free( Qx_1); break;
    case 7: free( Qx_2); break;
    case 8: free( seq); break;
    case 9: free( etaN); break;
    case 10: free( etaN_space); break;
  }
}
/* ****** */

DBL_TYPE * pf;
int ** intSeqs;
int * seqlengths;
int * nStrands_arr;
int * permSymmetries;
DBL_TYPE * temps;
int nblocks;

#define cudaCheck(call) {\
  cudaError_t e = (call); \
  if(e != cudaSuccess) {\
    fprintf(stderr, #call ": %s", cudaGetErrorString(e)); \
    abort(); \
  }\
}\


extern "C" void pfuncInitialize(
    int nblocks_in,
    DBL_TYPE temp_lo, DBL_TYPE temp_hi, DBL_TYPE temp_step,
    DBL_TYPE sodium_conc, DBL_TYPE magnesium_conc,
    int long_helix, int dangletype, int dnarnacount) {

  nblocks = nblocks_in;
  SODIUM_CONC = sodium_conc;
  MAGNESIUM_CONC = magnesium_conc;
  USE_LONG_HELIX_FOR_SALT_CORRECTION = long_helix;
  DANGLETYPE = dangletype;
  DNARNACOUNT = dnarnacount;

  // Load Energy Models
  int ntemps = ceil((temp_hi - temp_lo) / temp_step) + 1;
  temp_lo += ZERO_C_IN_KELVIN;
  temp_hi += ZERO_C_IN_KELVIN;
  energy_model_t *em;
  em = (energy_model_t*)malloc(ntemps * sizeof(energy_model_t));
  DBL_TYPE temp_k = temp_lo;
  for(int i = 0; i < ntemps; ++i) {
    LoadEnergies(em + i, temp_k);
    temp_k += temp_step;
  }
  // Allocate device memory and copy
  energy_model_t *em_d;
  cudaMalloc(&em_d, ntemps * sizeof(energy_model_t));
  cudaMemcpy(em_d, em, ntemps * sizeof(energy_model_t), cudaMemcpyHostToDevice);

  // Set device constants
  cudaMemcpyToSymbol(energies, &em_d, sizeof(energy_model_t*));
  cudaMemcpyToSymbol(t_hi, &temp_hi, sizeof(DBL_TYPE));
  cudaMemcpyToSymbol(t_lo, &temp_lo, sizeof(DBL_TYPE));
  cudaMemcpyToSymbol(t_step, &temp_step, sizeof(DBL_TYPE));

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
    cudaMallocManaged(&temps, nblocks * sizeof(DBL_TYPE))
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

extern "C" void pfuncMulti(char ** inputSeqs, int nseqs, int * permSym,
    DBL_TYPE * temp_Cs, DBL_TYPE * result) {

  for (int s = 0; s < nseqs; s += nblocks) {
    for(int i = 0; i < nblocks; ++i) {
      int len = strlen(inputSeqs[s + i]);
      convertSeq(inputSeqs[s + i], intSeqs[i], len);
      seqlengths[i] = getSequenceLengthInt(intSeqs[i], nStrands_arr + i);
      permSymmetries[i] = permSym[s + i];
      temps[i] = temp_Cs[s + i] + ZERO_C_IN_KELVIN;
    }

    pfuncFullWithSymHelper<<<nblocks, 256>>>(
        pf, intSeqs, seqlengths, nStrands_arr, permSymmetries, temps
    );
    cudaDeviceSynchronize();

    for(int i = 0; i < nblocks; ++i) {
      result[s + i] = pf[i];
    }
  }
}
