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

#define cudaCheck(call) {\
  cudaError_t e = (call); \
  if(e != cudaSuccess) {\
    fprintf(stderr, #call ": %s", cudaGetErrorString(e)); \
    abort(); \
  }\
}\


class PFMemory {
  public:
    DBL_TYPE *Q;
    DBL_TYPE *Qb;
    DBL_TYPE *Qm;

    DBL_TYPE *Qs;
    DBL_TYPE *Qms;

    DBL_TYPE *Qx;
    DBL_TYPE *Qx_1;
    DBL_TYPE *Qx_2;

    int *seq;
    int **etaN;
    int *etaN_space;

    int arraySize;
    int seqlength;

    PFMemory() {}
    ~PFMemory() {}
    void init(int);
    void free();
    __device__ void clear();
};
void PFMemory::init(int seqlength) {
  this->seqlength = seqlength;
  arraySize = seqlength*(seqlength+1)/2+(seqlength+1);

  cudaCheck(cudaMalloc(&Q, arraySize * sizeof(DBL_TYPE)));
  cudaCheck(cudaMalloc(&Qb, arraySize * sizeof(DBL_TYPE)));
  cudaCheck(cudaMalloc(&Qm, arraySize * sizeof(DBL_TYPE)));

  cudaCheck(cudaMalloc(&Qs, arraySize * sizeof(DBL_TYPE)));
  cudaCheck(cudaMalloc(&Qms, arraySize * sizeof(DBL_TYPE)));

  cudaCheck(cudaMalloc(&Qx, arraySize/2 * sizeof(DBL_TYPE)));
  cudaCheck(cudaMalloc(&Qx_1, arraySize/2 * sizeof(DBL_TYPE)));
  cudaCheck(cudaMalloc(&Qx_2, arraySize/2 * sizeof(DBL_TYPE)));

  cudaCheck(cudaMalloc(&seq, (seqlength + 1) * sizeof(int)));
  cudaCheck(cudaMalloc(&etaN, arraySize * sizeof(int*)));
  cudaCheck(cudaMalloc(&etaN_space, arraySize * 2 * sizeof(int)));

}
void PFMemory::free() {
  cudaCheck(cudaFree(Q));
  cudaCheck(cudaFree(Qb));
  cudaCheck(cudaFree(Qm));

  cudaCheck(cudaFree(Qs));
  cudaCheck(cudaFree(Qms));

  cudaCheck(cudaFree(Qx));
  cudaCheck(cudaFree(Qx_1));
  cudaCheck(cudaFree(Qx_2));

  cudaCheck(cudaFree(seq));
  cudaCheck(cudaFree(etaN));
  cudaCheck(cudaFree(etaN_space));
}

__device__ void PFMemory::clear() {
  memset(Q, 0, arraySize * sizeof(DBL_TYPE));
  memset(Qb, 0, arraySize * sizeof(DBL_TYPE));
  memset(Qm, 0, arraySize * sizeof(DBL_TYPE));

  memset(Qs, 0, arraySize * sizeof(DBL_TYPE));
  memset(Qms, 0, arraySize * sizeof(DBL_TYPE));

  memset(Qx, 0, arraySize/2 * sizeof(DBL_TYPE));
  memset(Qx_1, 0, arraySize/2 * sizeof(DBL_TYPE));
  memset(Qx_2, 0, arraySize/2 * sizeof(DBL_TYPE));

  memset(seq, 0, (seqlength + 1) * sizeof(int));
  memset(etaN, 0, arraySize * sizeof(int*));
  memset(etaN_space, 0, arraySize * 2 * sizeof(int));
}

__constant__ PFMemory *pf_mem;

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
    PFMemory *mem = pf_mem + blockIdx.x;
    mem->clear();
    Q    = mem->Q;
    Qb   = mem->Qb;
    Qm   = mem->Qm;
    Qx   = mem->Qx;
    Qx_1 = mem->Qx_1;
    Qx_2 = mem->Qx_2;
    Qs   = mem->Qs;
    Qms  = mem->Qms;
    seq  = mem->seq;
    etaN = mem->etaN;
    etaN_space = mem->etaN_space;

    processMultiSequence( inputSeq, seqlength, nStrands, seq, nicks);

    // Allocate and Initialize Matrices
    int arraySize = seqlength*(seqlength+1)/2+(seqlength+1);

    for (int i = 0; i < arraySize; ++i) {
      etaN[i] = etaN_space + (2 * i);
    }
    InitEtaN( etaN, nicks, seqlength);
    nonZeroInit( Q, seq, seqlength, em);

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
}
/* ****** */

DBL_TYPE * pf;
int ** intSeqs;
int * seqlengths;
int * nStrands_arr;
int * permSymmetries;
DBL_TYPE * temps;
int nblocks;
int nthreads;

energy_model_t *em_host;
energy_model_t *em_dev;

PFMemory *pfm_host;
PFMemory *pfm_dev;

extern "C" void pfuncInitialize(
    int nblocks_in,
    int nthreads_in,
    int max_seqlen,
    DBL_TYPE temp_lo, DBL_TYPE temp_hi, DBL_TYPE temp_step,
    DBL_TYPE sodium_conc, DBL_TYPE magnesium_conc,
    int long_helix, int dangletype, int dnarnacount) {

  nblocks = nblocks_in;
  nthreads = nthreads_in;
  SODIUM_CONC = sodium_conc;
  MAGNESIUM_CONC = magnesium_conc;
  USE_LONG_HELIX_FOR_SALT_CORRECTION = long_helix;
  DANGLETYPE = dangletype;
  DNARNACOUNT = dnarnacount;

  // Load Energy Models
  int ntemps = ceil((temp_hi - temp_lo) / temp_step) + 1;
  temp_lo += ZERO_C_IN_KELVIN;
  temp_hi += ZERO_C_IN_KELVIN;

  em_host = (energy_model_t*)malloc(ntemps * sizeof(energy_model_t));

  DBL_TYPE temp_k = temp_lo;
  for(int i = 0; i < ntemps; ++i) {
    LoadEnergies(em_host + i, temp_k);
    temp_k += temp_step;
  }
  // Allocate device memory and copy
  cudaMalloc(&em_dev, ntemps * sizeof(energy_model_t));
  cudaMemcpy(em_dev, em_host, ntemps * sizeof(energy_model_t), cudaMemcpyHostToDevice);

  // Set device constants
  cudaMemcpyToSymbol(energies, &em_dev, sizeof(energy_model_t*));
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

  // perform allocations
  pfm_host = new PFMemory[nblocks];
  for(int i = 0; i < nblocks; ++i) { pfm_host[i].init(max_seqlen); }

  cudaMalloc(&pfm_dev, nblocks * sizeof(PFMemory));
  cudaMemcpy(pfm_dev, pfm_host, nblocks * sizeof(PFMemory), cudaMemcpyHostToDevice);

  // copy symbol
  cudaMemcpyToSymbol(pf_mem, &pfm_dev, sizeof(PFMemory*));
}

extern "C" void pfuncShutdown() {
  // free device memory for energy models
  cudaCheck(cudaFree(em_dev));
  // free host memory for energy models
  free(em_host);

  // free unified memory for parameters
  cudaCheck(cudaFree(pf));
  cudaCheck(cudaFree(seqlengths));
  cudaCheck(cudaFree(nStrands_arr));
  cudaCheck(cudaFree(permSymmetries));
  cudaCheck(cudaFree(temps));
  for (int i = 0; i < nblocks; ++i) {
    cudaCheck(cudaFree(intSeqs[i]));
  }
  cudaCheck(cudaFree(intSeqs));

  // free device memory for PF arrays
  for(int i = 0; i < nblocks; ++i) {
    pfm_host[i].free();
  }
  // free device memory for PF array pointers
  cudaCheck(cudaFree(pfm_dev));
  // free host memory for PF array pointers
  delete [] pfm_host;

}

extern "C" void pfuncMulti(char ** inputSeqs, int nseqs, int * permSym,
    DBL_TYPE * temp_Cs, DBL_TYPE * result) {

  for (int s = 0; s < nseqs; s += nblocks) {
    int njobs = MIN(nblocks, nseqs - s);
    for(int i = 0; i < njobs; ++i) {
      int len = strlen(inputSeqs[s + i]);
      convertSeq(inputSeqs[s + i], intSeqs[i], len);
      seqlengths[i] = getSequenceLengthInt(intSeqs[i], nStrands_arr + i);
      permSymmetries[i] = permSym[s + i];
      temps[i] = temp_Cs[s + i] + ZERO_C_IN_KELVIN;
    }

    pfuncFullWithSymHelper<<<njobs, nthreads>>>(
        pf, intSeqs, seqlengths, nStrands_arr, permSymmetries, temps
    );
    cudaDeviceSynchronize();

    for(int i = 0; i < njobs; ++i) {
      result[s + i] = pf[i];
    }
  }
}
