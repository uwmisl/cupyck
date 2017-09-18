/*
  DNAExternals.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 3/2006, Justin Bois 1/2007

  External variables describing fundamental properties of DNA, such as
  base types, dangles, etc.
*/

#ifndef DNAEXTERNALS_H
#define DNAEXTERNALS_H

#include "pfuncUtilsConstants.h"

extern int NUPACK_VALIDATE;
extern int mfe_sort_method; // A constant to allow forced sort-by-structure
extern int NupackShowHelp;

typedef struct {
  DBL_TYPE Stack[36];
  DBL_TYPE loop37[90];
  int tloops[6*4096];//has tetra loop sequences+cp, (6 extern ints per tetra loop)
  DBL_TYPE tloop_energy[ 4096]; //energies of tetraloops
  int triloops[5*1024]; //triloops equences + closing pairs
  DBL_TYPE triloop_energy[ 2048]; //number of triloops

  //Mismatch energies  (see functions.h)
  DBL_TYPE MMEnergiesHP[6*16];
  DBL_TYPE MMEnergiesIL[256];
  DBL_TYPE IL_SInt2[16*36]; //Symmetric extern interior Loops, size 2
  DBL_TYPE IL_SInt4[256*36]; // Symmetric extern interior Loops, size 4
  DBL_TYPE IL_AsInt1x2[64*36]; // Asymmetric extern interior Loop, size 3
  DBL_TYPE dangle_energy[48]; // Dangle Energies
  DBL_TYPE asymmetry_penalty[4]; // Asymmetric loop penalties
  DBL_TYPE max_asymmetry;
  DBL_TYPE *sizeTerm;
  DBL_TYPE *pairPr;
  DBL_TYPE bimolecular;

  DBL_TYPE at_penalty;
  DBL_TYPE polyc3;
  DBL_TYPE polycslope;
  DBL_TYPE polycint;
  DBL_TYPE alpha_1; //multiloop penalties
  DBL_TYPE alpha_2;
  DBL_TYPE alpha_3;
  DBL_TYPE beta_1; //pseudoknot penalties
  DBL_TYPE beta_2;
  DBL_TYPE beta_3;
  DBL_TYPE beta_1m;
  DBL_TYPE beta_1p;

  DBL_TYPE salt_correction;

  DBL_TYPE temp_k;
  int dangletype;
  int dnarnacount;

} energy_model_t;

/*
extern DBL_TYPE Stack[36];
extern DBL_TYPE loop37[90];
extern int tloops[6*4096];//has tetra loop sequences+cp, (6 extern ints per tetra loop)
extern DBL_TYPE tloop_energy[ 4096]; //energies of tetraloops
extern int triloops[5*1024]; //triloops equences + closing pairs
extern DBL_TYPE triloop_energy[ 2048]; //number of triloops

//Mismatch energies  (see functions.h)
extern DBL_TYPE MMEnergiesHP[6*16];
extern DBL_TYPE MMEnergiesIL[256];
extern DBL_TYPE IL_SInt2[16*36]; //Symmetric extern interior Loops, size 2
extern DBL_TYPE IL_SInt4[256*36]; // Symmetric extern interior Loops, size 4
extern DBL_TYPE IL_AsInt1x2[64*36]; // Asymmetric extern interior Loop, size 3
extern DBL_TYPE dangle_energy[48]; // Dangle Energies
extern DBL_TYPE asymmetry_penalty[4]; // Asymmetric loop penalties
extern DBL_TYPE max_asymmetry;
extern long int maxGapIndex;
extern DBL_TYPE *sizeTerm;
extern DBL_TYPE *pairPr;
extern DBL_TYPE * EXTERN_Q ;
extern DBL_TYPE * EXTERN_QB ;
extern DBL_TYPE BIMOLECULAR;

extern DBL_TYPE AT_PENALTY;
extern DBL_TYPE POLYC3;
extern DBL_TYPE POLYCSLOPE;
extern DBL_TYPE POLYCINT;
extern DBL_TYPE ALPHA_1; //multiloop penalties
extern DBL_TYPE ALPHA_2;
extern DBL_TYPE ALPHA_3;
extern DBL_TYPE BETA_1; //pseudoknot penalties
extern DBL_TYPE BETA_2;
extern DBL_TYPE BETA_3;
extern DBL_TYPE BETA_1M;
extern DBL_TYPE BETA_1P;
*/
extern DBL_TYPE SODIUM_CONC;
extern DBL_TYPE MAGNESIUM_CONC;
extern int USE_LONG_HELIX_FOR_SALT_CORRECTION;
//extern DBL_TYPE SALT_CORRECTION;
extern DBL_TYPE TEMP_K;
extern int DANGLETYPE;
extern int DNARNACOUNT;
extern int DO_PSEUDOKNOTS;
extern int ONLY_ONE_MFE;
extern int USE_MFE;
extern char PARAM_FILE[100];

#endif

