/*
  DNAGlobals.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 3/2006, Justin Bois 1/2007

  Global variables describing fundamental properties of DNA, such as
  base types, dangles, etc.
*/
#include "pfuncUtilsHeader.h"

DBL_TYPE SODIUM_CONC;
DBL_TYPE MAGNESIUM_CONC;
int USE_LONG_HELIX_FOR_SALT_CORRECTION;
int DANGLETYPE;
int DNARNACOUNT;
int DO_PSEUDOKNOTS;

char PARAM_FILE[100]="";

