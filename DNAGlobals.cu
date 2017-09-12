/*
  DNAGlobals.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 3/2006, Justin Bois 1/2007

  Global variables describing fundamental properties of DNA, such as
  base types, dangles, etc.
*/
#include "pfuncUtilsHeader.h"

int mfe_sort_method = 0;
int NUPACK_VALIDATE = 0;
int NupackShowHelp = 0;

DEV energy_model_t ENERGIES;

DBL_TYPE SODIUM_CONC;
DBL_TYPE MAGNESIUM_CONC;
int USE_LONG_HELIX_FOR_SALT_CORRECTION;
DBL_TYPE TEMP_K;
int DANGLETYPE;
int DNARNACOUNT;
int DO_PSEUDOKNOTS;
int ONLY_ONE_MFE;
int USE_MFE;

DBL_TYPE *pairPr;
DBL_TYPE *pairPrPbg;
DBL_TYPE *pairPrPb;

DBL_TYPE * EXTERN_QB = NULL;
DBL_TYPE * EXTERN_Q = NULL;

char PARAM_FILE[100]="";

int nupack_sample = 0;
int nupack_num_samples = 0;
