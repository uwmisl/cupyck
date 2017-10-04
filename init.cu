/*
  init.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 3/2006, Justin Bois 1/2007

  Functions to be run once at the beginning of the
  partition function algorithm
*/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<assert.h>

#include "pfuncUtilsHeader.h"
#include "DNAExternals.h"

/* **************************** */
int getSequenceLength( char *seq, int *nStrands /*, seq2, nicks*/) {
  
  int i; //position in sequence
  int done = FALSE;
  int seqlength = 0;
  char tmpC;

  i = 0;
  *nStrands = 1;
  while( done == FALSE) {
    tmpC = toupper( seq[i]);
    if( tmpC == '+') {
      (*nStrands)++;
      //nicks[ seqlength] = 1;
    }
    else if( tmpC != 'A' && tmpC != 'T' && tmpC != 'C'
            && tmpC != 'G' && tmpC != 'U') {
              done = TRUE; 
            }
    else {
      //seq2[ seqlength] = tmpC;
      seqlength++;
    }
    i++;
  }

  if( seqlength > MAXSEQLENGTH) {
    printf("Sequences longer than maximum of %d\n", MAXSEQLENGTH);
    assert(0);
  }

  return seqlength;
}

/* **************************** */
int getSequenceLengthInt( int seq[], int *nStrands /*, seq2, nicks*/) {

  int i; //position in sequence
  int done = FALSE;
  int seqlength = 0;
  int tmpC;

  i = 0;
  *nStrands = 1;
  while( done == FALSE) {
    tmpC = seq[i];
    if( tmpC == STRAND_PLUS) {
      (*nStrands)++;
      //nicks[ seqlength] = 1;
    }
    else if( tmpC != BASE_A && tmpC != BASE_T && tmpC != BASE_C
            && tmpC != BASE_G && tmpC != BASE_U) {
              done = TRUE; 
            }
    else {
      //seq2[ seqlength] = tmpC;
      seqlength++;
    }
    i++;
  }

  if( seqlength > MAXSEQLENGTH) {
    printf("Sequences longer than maximum of %d\n", MAXSEQLENGTH);
    assert(0);
  }

  return seqlength;
}

/* *********** */
DEV
void processMultiSequence( int inputSeq[], int seqlength, int nStrands,
                           int seq[], int nicks[]) {
  
  int i, j;
  int nNick = 0;
  int done;

  j = 0;
  for( i = 0; i < seqlength; i++) {
    done = FALSE;
    while( !done) {
      if( inputSeq[j] == BASE_A || inputSeq[j] == BASE_C || inputSeq[j] == BASE_G ||
        inputSeq[j] == BASE_T || inputSeq[j] == BASE_U) {
        done = TRUE;
        seq[i] = inputSeq[j];
      }
      else if( inputSeq[j] == STRAND_PLUS) {
        nicks[nNick++] = i-1;
      }
      j++;

      if( j >= seqlength + nStrands && !done) {
        printf("\nError in processing sequence:\n%d\n", inputSeq[0]);
        printf("seqlength = %d, nStrands = %d\n", seqlength, nStrands);
        assert(0);
      }
    }
  }
  seq[ seqlength] = -1;
}


/* ********************************************* */
DEV
void InitLDoublesMatrix( DBL_TYPE **Q, int size, char name[]) {
  // Allocate cleared memory for a DBL_TYPEs matrix.
  *Q =  (DBL_TYPE *) malloc( size * sizeof( DBL_TYPE));
  memset(*Q, 0, size * sizeof(DBL_TYPE));
  if( *Q == NULL) {
    printf("InitLDoublesMatrix: unable to allocate %lu bytes for %s!\n", size * sizeof( DBL_TYPE),  name);
    assert(0);
  }
}

void ClearLDoublesMatrix(DBL_TYPE **Q, int size, char name[]) {
  memset(*Q, 0,size * sizeof(DBL_TYPE));
}

/* ******************************************** */
DEV
void nonZeroInit( DBL_TYPE Q[], int seq[], int seqlength, energy_model_t *em) {
  // Set Q[i, i-1] = 1.
  int i;

  for( i = 0; i <= seqlength; i++) {
    Q[ pf_index(i, i-1, seqlength)] = ExplDangle(i,i-1,seq,seqlength, em);
  }
}


/* *************************************************** */
DEV
void manageQx( DBL_TYPE **Qx, DBL_TYPE **Qx_1, DBL_TYPE **Qx_2, int len, int seqlength) {
  // Allocate and deallocate QbIx matrices

  int i;
  DBL_TYPE *temp;
  int arraySize = seqlength*(seqlength+1)/2 + (seqlength+1);

  if( len > 11) {

    temp = *Qx;
    *Qx = *Qx_1;
    *Qx_1 = *Qx_2;
    *Qx_2 = temp;

    for(i = 0; i < arraySize/2; ++i) {
      (*Qx_2)[i] = 0;
    }

  }
}


/* ************************************** */
DBL_TYPE computeSaltCorrection(DBL_TYPE sodiumConc, DBL_TYPE magnesiumConc,
			       int useLongHelix, DBL_TYPE temp_k) {

  // No correction for RNA since we don't have parameters
  if (DNARNACOUNT != DNA || (sodiumConc == 1.0 && magnesiumConc == 0.0)) { 
    return 0.0;
  }

  // Ignore magnesium for long helix mode (not cited why, for consistency with Mfold)
  if (useLongHelix) { 
    return -(0.2 + 0.175*log(sodiumConc)) * temp_k / 310.15;
  }

  return -0.114*log(sodiumConc + 3.3*sqrt(magnesiumConc)) * temp_k / 310.15;
}


/* ************************************** */
void LoadEnergies(energy_model_t *em, DBL_TYPE temp_k) {
  
  const char *default_param_files[] = { "dna1998", "rna1995", "rna1999"};
  
  DBL_TYPE H_loop37[90];  
  DBL_TYPE H_tloop_energy[4096];
  DBL_TYPE H_triloop_energy[2048];
  DBL_TYPE H_MMEnergiesHP[6*16];
  DBL_TYPE H_MMEnergiesIL[256];
  DBL_TYPE H_IL_SInt2[16*36];
  DBL_TYPE H_IL_SInt4[256*36];
  DBL_TYPE H_IL_AsInt1x2[64*36];
  DBL_TYPE H_dangle_energy[48];
  DBL_TYPE H_asymmetry_penalty[4];
  DBL_TYPE H_max_asymmetry;
  DBL_TYPE H_Stack[36];
  DBL_TYPE H_ALPHA_1, H_ALPHA_2, H_ALPHA_3, H_BETA_1, H_BETA_2, 
  H_BETA_3, H_BETA_1M, H_BETA_1P;
  DBL_TYPE H_POLYC3, H_POLYCINT, H_POLYCSLOPE;
  DBL_TYPE H_AT_PENALTY;
  DBL_TYPE H_BIMOLECULAR;
  
  DBL_TYPE G_loop37[90];  
  DBL_TYPE G_tloop_energy[4096];
  DBL_TYPE G_triloop_energy[2048];
  DBL_TYPE G_MMEnergiesHP[6*16];
  DBL_TYPE G_MMEnergiesIL[256];
  DBL_TYPE G_IL_SInt2[16*36];
  DBL_TYPE G_IL_SInt4[256*36];
  DBL_TYPE G_IL_AsInt1x2[64*36];
  DBL_TYPE G_dangle_energy[48];
  DBL_TYPE G_asymmetry_penalty[4];
  DBL_TYPE G_max_asymmetry;
  DBL_TYPE G_Stack[36];
  DBL_TYPE G_ALPHA_1, G_ALPHA_2, G_ALPHA_3, G_BETA_1, G_BETA_2,
  G_BETA_3, G_BETA_1M, G_BETA_1P;
  DBL_TYPE G_POLYC3, G_POLYCINT, G_POLYCSLOPE;
  DBL_TYPE G_AT_PENALTY;
  DBL_TYPE G_BIMOLECULAR;
  DBL_TYPE water_conc;


  //temporary storage of data
  int nRead;
  int array[MAXLINE];
  char *token;
  char tetraloop[6] = "\0";
  char triloop[5] = "\0";
  int indexL, tmpIndex, index4;
  
  FILE *fp;	
  char line[MAXLINE];
  int i, j, k;
  char fileG[300] = "\0";
  char fileH[300] = "\0";
  char *nupackhome = NULL;
  char fileNameRoot[300] = "\0";
  
  static char parameterFileName[300] = "\0";
  

  em->temp_k = temp_k;
  em->dnarnacount = DNARNACOUNT;
  em->dangletype = DANGLETYPE;

  if( DNARNACOUNT == COUNT) {
    setParametersToZero(em);
    return;
  }
  
  // Get density of water
  water_conc = (DBL_TYPE) WaterDensity(temp_k - ZERO_C_IN_KELVIN);

  // Compute the salt correction
  em->salt_correction = computeSaltCorrection(
      SODIUM_CONC, MAGNESIUM_CONC,
      USE_LONG_HELIX_FOR_SALT_CORRECTION, temp_k);
  
  
  //Parameter file input.  If a path was given as a command line parameter,
  //follow this path and return an error if the parameter files cannot be found.
  //Otherwise, check the local directory first, and if the files are not there,
  //check the NUPACKHOME/parameters directory.
  
  if( DNARNACOUNT < USE_SPECIFIED_PARAMETERS_FILE) {
    strcpy( fileNameRoot, default_param_files[ DNARNACOUNT]);
    strcpy( parameterFileName, ""); //Used to check if parameter files need to be reloaded.
  }
  else if( DNARNACOUNT == USE_SPECIFIED_PARAMETERS_FILE) {
    strcpy( fileNameRoot, PARAM_FILE);
    strcpy( parameterFileName, PARAM_FILE); //store this to check if parameter reload is needed.
  }
  
  //check first for .dG parameter file using current directory as home
  strcpy( fileG, fileNameRoot);
  strcat( fileG, ".dG");
  
  if( !fileExists( fileG) ) {
    //if files not found, use environment variable NUPACKHOME as root
    nupackhome = getenv("NUPACKHOME");
    if( nupackhome != NULL) {
      strcpy( fileG, nupackhome);
      strcat( fileG, "/parameters/");
      strcat( fileG, fileNameRoot);
      strcat( fileG, ".dG");
    }
    else {
      fprintf(stderr, "Unable to find %s.dG and NUPACKHOME environment variable is not set\n",
             fileNameRoot);
      exit(1);
    }
  }
  
  if( ! fileExists( fileG)) {
    fprintf(stderr, "Unable to find file %s.dG locally or in NUPACKHOME = %s/parameters\n", 
           fileNameRoot, 
           nupackhome);
    fprintf(stderr, "%s\n", fileG);
    exit(1);
  }
  
  fp = fopen( fileG, "r");
  
  if( fp == NULL) {  // Make sure input file exits 
    fprintf(stderr, "Error opening loop data file: %s\n", fileG);
    exit(1);  
  }
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //Read in Stacking data
  for( i = 0; i < 6; i++) {
    
    nRead = 0;
    token = strtok( line, " ");
    int tmp_array;
    while( token != NULL) {
      tmp_array=0;
      if( sscanf( token, "%d", &tmp_array ) == 1) {
        array[nRead]=tmp_array;
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 6) {
      fprintf(stderr, "Error in stacking data format\n");
      exit(1);
    }
    for( j = 0; j < 6; j++) {
      em->Stack[i*6+j] = G_Stack[i*6+j] = (DBL_TYPE) array[j]/100.0;
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  for( i = 0; i < 3; i++) {
    while( line[0] == '>') {
      fgets( line, MAXLINE, fp);
    }
    
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead > 30) {
      fprintf(stderr, "Error in Loop energies data\n");
      exit(1);
    }
    for( j = 0; j < 30; j++) {
      if( nRead-1 >= j) {
        em->loop37[30*(2-i)+j] = G_loop37[30*(2-i)+j] = (DBL_TYPE) array[j]/100.0;
      }
      else {
        em->loop37[30*(2-i)+j] = G_loop37[30*(2-i)+j] = G_loop37[30*(2-i)+nRead-1]+
          1.75*kB*temp_k*LOG_FUNC( (j+1)/(1.0*nRead));
      }
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 5) {
    fprintf(stderr, "Error in asymmetry terms!\n");
    exit(1);
  }
  
  for( j = 0; j < 4; j++) {
    em->asymmetry_penalty[j] = G_asymmetry_penalty[j] = (DBL_TYPE) array[j]/100.0;
  }
  em->max_asymmetry = G_max_asymmetry = (DBL_TYPE) array[4]/100.0;
  
  //Triloops
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  for( i = 0; i < 2048; i++) {
    em->triloop_energy[i] = G_triloop_energy[i] = 0;
  }
  while( line[0] != '>') {
    
    if( sscanf( line, "%s %d", triloop, &(array[0]) ) == 2) {
      indexL = 0;
      for( i = 0 ; i < 5; i++) {
        tmpIndex = 1;
        for( j = 0; j < i; j++) {
          tmpIndex *= 4;
        }
        if( triloop[4-i] == 'C') {
          indexL += tmpIndex;
        }
        else if( triloop[4-i] == 'G') {
          indexL += tmpIndex*2;
        }
        else if( triloop[4-i] == 'U' || triloop[4-i] == 'T') {
          indexL += tmpIndex*3;
        }
        else if( triloop[4-i] != 'A') {
          fprintf(stderr, "Error in triloop indexing %s\n", triloop);
        }
      }
      em->triloop_energy[ indexL] = G_triloop_energy[ indexL] = (DBL_TYPE) array[0]/100.0;
    }
    else {
      fprintf(stderr, "Error in triloop data\n%s\n",line);
    }
    fgets( line, MAXLINE, fp);
  }
  
  //Tetraloops
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  for( i = 0; i < 4096; i++) {
    em->tloop_energy[i] = G_tloop_energy[i] = 0;
  }
  while( line[0] != '>') {
    if( sscanf( line, "%s %d", tetraloop, &(array[0]) ) == 2) {
      indexL = 0;
      for( i = 0 ; i < 6; i++) {
        tmpIndex = 1;
        for( j = 0; j < i; j++) {
          tmpIndex *= 4;
        }
        if( tetraloop[5-i] == 'C') {
          indexL += tmpIndex;
        }
        else if( tetraloop[5-i] == 'G') {
          indexL += tmpIndex*2;
        }
        else if( tetraloop[5-i] == 'U' || tetraloop[5-i] == 'T') {
          indexL += tmpIndex*3;
        }
        else if( tetraloop[5-i] != 'A') {
          fprintf(stderr, "Error in tetraloop indexing %s\n", tetraloop);
        }
      }
      em->tloop_energy[ indexL] = G_tloop_energy[ indexL] = (DBL_TYPE) array[0]/100.0;
    }
    else {
      fprintf(stderr, "Error in tetraloop data\n%s\n",line);
    }
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //fprintf(fp, "Mismatch Hairpin: \n");
  for( i = 0; i < 16; i++) {
    
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 6) {
      fprintf(stderr, "Error in mismatch hairpin format! %d\n", nRead);
      exit(1);
    }
    
    for( j = 0; j < 6; j++) {
      em->MMEnergiesHP[ 6*i + j] = G_MMEnergiesHP[ 6*i + j] = (DBL_TYPE) array[j]/100.0;
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //fprintf(fp, "Mismatch Interior: \n");
  for( i = 0; i < 16; i++) {
    
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 6) {
      fprintf(stderr, "Error in mismatch Interior format!\n");
      exit(1);
    }
    
    for( j = 0; j < 6; j++) {
      em->MMEnergiesIL[ 6*i + j] = G_MMEnergiesIL[ 6*i + j] = (DBL_TYPE) array[j]/100.0;
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //Read in Dangles
  for( i = 0; i < 6; i++) {
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 4) {
      fprintf(stderr, "1. Error in dangle data format!\n");
      exit(1);
    }
    for( j = 0; j < 4; j++) {
      em->dangle_energy[i*4+j]  = G_dangle_energy[i*4+j] = (DBL_TYPE) array[j]/100.0;
      
      if( DANGLETYPE == 0)  em->dangle_energy[i*4+j] = G_dangle_energy[i*4+j] = 0.0;

    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //Read in Dangles
  for( i = 0; i < 6; i++) {
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 4) {
      fprintf(stderr, "2. Error in dangle data format!\n");
      exit(1);
    }
    for( j = 0; j < 4; j++) {
      em->dangle_energy[24+ i*4+j] = G_dangle_energy[24+ i*4+j] = (DBL_TYPE) array[j]/100.0;
      if( DANGLETYPE == 0) em->dangle_energy[24+ i*4+j] = G_dangle_energy[24+ i*4+j] = 0;
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }


  //Multiloop parameters
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) ) ==1) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  
  if( nRead != 3) {
    fprintf(stderr, "3. Error in dangle data format!\n");
    exit(1);
  }
  
  em->alpha_1 = G_ALPHA_1 = (DBL_TYPE) array[0]/100.0;
  em->alpha_2 = G_ALPHA_2 = (DBL_TYPE) array[1]/100.0;
  em->alpha_3 = G_ALPHA_3 = (DBL_TYPE) array[2]/100.0;
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //AT PENALTY
  if( sscanf(line, "%d", &(array[0]) ) == 1) {
    em->at_penalty = G_AT_PENALTY = (DBL_TYPE) array[0]/100.0;
  }
  else {
    fprintf(stderr, "Error in AT PENALTY data\n");
    exit(1);
  }
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //1x1 interior loop
  for( i = 0; i < 36; i++) {
    fgets( line, MAXLINE, fp); //read in label
    for( j = 0; j < 4; j++) {
      
      nRead = 0;
      token = strtok( line, " ");
      while( token != NULL) {
        if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
          
          nRead++;
        }
        token = strtok( NULL, " ");
      }
      
      if( nRead != 4) {
        fprintf(stderr, "Error in 1x1 format!\n");
        exit(1);
      }
      
      for( k = 0; k < 4; k++) {
        em->IL_SInt2[ i*16 + j*4 + k] = G_IL_SInt2[ i*16 + j*4 + k] = (DBL_TYPE) array[ k]/100.0;
      }
      
      fgets( line, MAXLINE, fp);
    }
  }
  
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //2x2 interior loop
  for( i = 0; i < 36*16; i++) {
    fgets( line, MAXLINE, fp); //read in label
    for( j = 0; j < 4; j++) {
      
      nRead = 0;
      token = strtok( line, " ");
      while( token != NULL) {
        if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
          
          nRead++;
        }
        token = strtok( NULL, " ");
      }
      
      if( nRead != 4) {
        fprintf(stderr, "Error in 1x1 format!\n");
        exit(1);
      }
      
      for( k = 0; k < 4; k++) {
        em->IL_SInt4[ 1536*(i/96)+256*((i%96)/16) +
                 64*((i%16)/4) + 4*(i%4) + k*16 + j] = 
          G_IL_SInt4[ 1536*(i/96)+256*((i%96)/16) +
                     64*((i%16)/4) + 4*(i%4) + k*16 + j] = 
          (DBL_TYPE) array[ k]/100.0;
      }
      
      fgets( line, MAXLINE, fp);
    }
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //1x2 interior loop
  for( i = 0; i < 144; i++) {
    fgets( line, MAXLINE, fp); //read in label
    for( j = 0; j < 4; j++) {
      
      nRead = 0;
      token = strtok( line, " ");
      while( token != NULL) {
        if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
          
          nRead++;
        }
        token = strtok( NULL, " ");
      }
      
      if( nRead != 4) {
        fprintf(stderr, "Error in 1x1 format!\n");
        exit(1);
      }
      
      for( k = 0; k < 4; k++) {
        em->IL_AsInt1x2[ 384*(i/24) + 4*((i%24)/4) + 24*(i%4) +
                    96*j + k] = G_IL_AsInt1x2[ 384*(i/24) + 4*((i%24)/4) + 24*(i%4) +
                                              96*j + k] = (DBL_TYPE) array[ k]/100.0;
      }
      
      fgets( line, MAXLINE, fp);
    }
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //polyC hairpin parameters
  nRead = 0;

  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 3) {
    fprintf(stderr, "4. Error in polyC hairpin parameters!\n");
    exit(1);
  }
  
  em->polyc3 = G_POLYC3 = (DBL_TYPE) array[0]/100.0;
  em->polycslope = G_POLYCSLOPE = (DBL_TYPE) array[1]/100.0;
  em->polycint = G_POLYCINT = (DBL_TYPE) array[2]/100.0;
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  
  //Pseudoknot parameters
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 5) {
    fprintf(stderr, "5. Error in dangle data format!\n");
    exit(1);
  }
  
  em->beta_1 = G_BETA_1 = (DBL_TYPE) array[0]/100.0;
  em->beta_2 = G_BETA_2 = (DBL_TYPE) array[1]/100.0;
  em->beta_3 = G_BETA_3 = (DBL_TYPE) array[2]/100.0;
  em->beta_1m = G_BETA_1M = (DBL_TYPE) array[3]/100.0;
  em->beta_1p = G_BETA_1P = (DBL_TYPE) array[4]/100.0;
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //BIMOLECULAR TERM
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 1) {
    fprintf(stderr, "Error in BIMOLECULAR format!\n");
    exit(1);
  }
  
  G_BIMOLECULAR = (DBL_TYPE) array[0]/100.0;
  em->bimolecular = G_BIMOLECULAR - kB*temp_k*LOG_FUNC( water_conc);
  
  fclose( fp);
  
  /* ****************************** */
  //Load Enthalpies and calculate modified G
  
  //If Temperature == 37 C, add the salt correction andskip he dH file.
  //  Else, make sure it is present
  if( temp_k > 37.0 + ZERO_C_IN_KELVIN - 0.001 && temp_k < 37.0 + ZERO_C_IN_KELVIN + 0.001) {

    // Make the salt corrections
    // Stacked bases
    for (i = 0; i < 36; i++) {
      em->Stack[i] += em->salt_correction;
    }
    
    // Loop correction.  Covers all hairpins, and bulges, but overcounts for
    // 1-base bulges.  This is corrected in the function InteriorEnergyFull.
    // This also covers large (non-tabulated) interior loops.
    for (i = 0; i < 90; i++) {
      em->loop37[i] += em->salt_correction;
    }
    
    // Corrections for tabulated interior loops
    for (i = 0; i < 16*36; i++) {
      em->IL_SInt2[i] += em->salt_correction;
    }
    for (i = 0; i < 256*36; i++) {
      em->IL_SInt4[i] += em->salt_correction;
    }
    for (i = 0; i < 64*36; i++) {
      em->IL_AsInt1x2[i] += em->salt_correction;
    }
    
    // Multiloop
    em->alpha_1 += em->salt_correction;
    
    return;
  }
  
  //check first for parameter files using current directory as home
  strcpy( fileH, fileNameRoot);
  strcat( fileH, ".dH");
  
  if( !fileExists( fileH) ) {
    //if files not found, use environment variable NUPACKHOME as root
    nupackhome = getenv("NUPACKHOME");
    if( nupackhome != NULL) {
      
      strcpy( fileH, nupackhome);
      strcat( fileH, "/parameters/");
      strcat( fileH, fileNameRoot);
      strcat( fileH, ".dH");
    }
    else {
      fprintf(stderr, "Unable to find %s.dH locally, and NUPACKHOME environment variable is not set.\n",
             fileNameRoot);
      fprintf(stderr, "Consequently, your temperature must be set to 37.0 C, not %.1f.  Job Aborted.\n", 
             (float) (temp_k - ZERO_C_IN_KELVIN));
      exit(1);
    }
  }
  
  if( ! fileExists( fileH)) {
    fprintf(stderr, "Unable to find file %s.dH locally or in NUPACKHOME = %s\n", fileNameRoot,
           nupackhome);
    fprintf(stderr, "Consequently, your temperature must be set to 37.0 C, not %.1f.  Job Aborted.\n",
           (float) (temp_k - ZERO_C_IN_KELVIN) );
    exit(1);
  }
  
  
  fp = fopen( fileH, "r");
  if( fp == NULL) {  // Make sure input file exits 
    fprintf(stderr, "Error opening loop data file: %s\n", fileH);
    exit(1);  
  }
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //Read in Stacking data
  for( i = 0; i < 6; i++) {
    
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) ) == 1) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 6) {
      fprintf(stderr, "Error in stacking data format\n");
      exit(1);
    }
    for( j = 0; j < 6; j++) {
      H_Stack[i*6+j] = (DBL_TYPE) array[j]/100.0;
      em->Stack[i*6+j] = (G_Stack[i*6+j] - H_Stack[i*6+j])*temp_k/310.15
        + H_Stack[i*6+j];
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  for( i = 0; i < 3; i++) {
    while( line[0] == '>') {
      fgets( line, MAXLINE, fp);
    }
    
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead > 30) {
      fprintf(stderr, "Error in Loop energies data\n");
      exit(1);
    }
    for( j = 0; j < 30; j++) {
      if( nRead-1 >= j) {
        H_loop37[30*(2-i)+j] = (DBL_TYPE) array[j]/100.0;
      }
      else {
        H_loop37[30*(2-i)+j] = H_loop37[30*(2-i)+nRead-1]+
          1.75*kB*temp_k*LOG_FUNC( (j+1)/(1.0*nRead));
      }
      
      em->loop37[30*(2-i)+j] = (G_loop37[30*(2-i)+j] - H_loop37[30*(2-i)+j])*
        temp_k/310.15 + H_loop37[30*(2-i)+j];
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 5) {
    fprintf(stderr, "Error in asymmetry terms!\n");
    exit(1);
  }
  
  for( j = 0; j < 4; j++) {
    H_asymmetry_penalty[j] = (DBL_TYPE) array[j]/100.0;
    em->asymmetry_penalty[j] = (G_asymmetry_penalty[j] - 
                            H_asymmetry_penalty[j])*
      temp_k/310.15 + H_asymmetry_penalty[j];
  }
  H_max_asymmetry = (DBL_TYPE) array[4]/100.0;
  em->max_asymmetry = (G_max_asymmetry - 
                   H_max_asymmetry)* temp_k/310.15 + H_max_asymmetry;
  
  
  //Triloops
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  for( i = 0; i < 2048; i++) {
    H_triloop_energy[i] = 0;
    em->triloop_energy[i] = (G_triloop_energy[i] - 
                         H_triloop_energy[i])* temp_k/310.15 + 
      H_triloop_energy[i];
  }
  while( line[0] != '>') {
    
    if( sscanf( line, "%s %d", triloop, &(array[0]) ) == 2) {
      indexL = 0;
      for( i = 0 ; i < 5; i++) {
        tmpIndex = 1;
        for( j = 0; j < i; j++) {
          tmpIndex *= 4;
        }
        if( triloop[4-i] == 'C') {
          indexL += tmpIndex;
        }
        else if( triloop[4-i] == 'G') {
          indexL += tmpIndex*2;
        }
        else if( triloop[4-i] == 'U' || triloop[4-i] == 'T') {
          indexL += tmpIndex*3;
        }
        else if( triloop[4-i] != 'A') {
          fprintf(stderr, "Error in triloop indexing %s\n", triloop);
        }
      }
      H_triloop_energy[ indexL] = (DBL_TYPE) array[0]/100.0;
      em->triloop_energy[indexL] = (G_triloop_energy[indexL] - 
                                H_triloop_energy[indexL])* temp_k/310.15 + 
        H_triloop_energy[indexL];
    }
    else {
      fprintf(stderr, "Error in triloop data\n%s\n",line);
    }
    fgets( line, MAXLINE, fp);
  }
  
  //Tetraloops
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  for( i = 0; i < 4096; i++) {
    H_tloop_energy[i] = 0;
    em->tloop_energy[i] = (G_tloop_energy[i] - 
                       H_tloop_energy[i])* temp_k/310.15 + 
      H_tloop_energy[i];
    //printf("%f ", temp_k);
  }
  while( line[0] != '>') {
    
    if( sscanf( line, "%s %d", tetraloop, &(array[0]) ) == 2) {
      indexL = 0;
      for( i = 0 ; i < 6; i++) {
        tmpIndex = 1;
        for( j = 0; j < i; j++) {
          tmpIndex *= 4;
        }
        if( tetraloop[5-i] == 'C') {
          indexL += tmpIndex;
        }
        else if( tetraloop[5-i] == 'G') {
          indexL += tmpIndex*2;
        }
        else if( tetraloop[5-i] == 'U' || tetraloop[5-i] == 'T') {
          indexL += tmpIndex*3;
        }
        else if( tetraloop[5-i] != 'A') {
          fprintf(stderr, "Error in tetraloop indexing %s\n", tetraloop);
        }
      }
      
      H_tloop_energy[ indexL] = (DBL_TYPE) array[0]/100.0;
      em->tloop_energy[indexL] = (G_tloop_energy[indexL] - 
                              H_tloop_energy[indexL])* temp_k/310.15 + 
        H_tloop_energy[indexL];
    }
    else {
      fprintf(stderr, "Error in tetraloop data\n%s\n",line);
    }
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //fprintf(fp, "Mismatch Hairpin: \n");
  for( i = 0; i < 16; i++) {
    
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 6) {
      fprintf(stderr, "Error in mismatch hairpin format! %d\n", nRead);
      exit(1);
    }
    
    for( j = 0; j < 6; j++) {
      H_MMEnergiesHP[ 6*i + j] = (DBL_TYPE) array[j]/100.0;
      em->MMEnergiesHP[6*i+j] = (G_MMEnergiesHP[6*i+j] - 
                             H_MMEnergiesHP[6*i+j])* temp_k/310.15 + 
        H_MMEnergiesHP[6*i+j];
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //fprintf(fp, "Mismatch Interior: \n");
  for( i = 0; i < 16; i++) {
    
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 6) {
      fprintf(stderr, "Error in mismatch Interior format!\n");
      exit(1);
    }
    
    for( j = 0; j < 6; j++) {
      H_MMEnergiesIL[ 6*i + j] = (DBL_TYPE) array[j]/100.0;
      em->MMEnergiesIL[6*i+j] = (G_MMEnergiesIL[6*i+j] - 
                             H_MMEnergiesIL[6*i+j])* temp_k/310.15 + 
        H_MMEnergiesIL[6*i+j];
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //Read in Dangles
  for( i = 0; i < 6; i++) {
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 4) {
      fprintf(stderr, "6. Error in dangle data format!\n");
      exit(1);
    }
    for( j = 0; j < 4; j++) {
      H_dangle_energy[i*4+j] = (DBL_TYPE) array[j]/100.0;
      if( DANGLETYPE == 0) H_dangle_energy[i*4+j] = 0;
      
      em->dangle_energy[i*4+j] = (G_dangle_energy[i*4+j] - 
                              H_dangle_energy[i*4+j])* temp_k/310.15 + 
        H_dangle_energy[i*4+j];
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //Read in Dangles
  for( i = 0; i < 6; i++) {
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 4) {
      fprintf(stderr, "7. Error in dangle data format!\n");
      exit(1);
    }
    for( j = 0; j < 4; j++) {
      H_dangle_energy[24+ i*4+j] = (DBL_TYPE) array[j]/100.0;
      if( DANGLETYPE == 0) 
        H_dangle_energy[24+ i*4+j] = 0;
      
      em->dangle_energy[24+i*4+j] = (G_dangle_energy[24+i*4+j] - 
                                 H_dangle_energy[24+i*4+j])* temp_k/310.15 + 
        H_dangle_energy[24+i*4+j];
      
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  //Multiloop parameters
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) ) ==1) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  
  if( nRead != 3) {
    fprintf(stderr, "8. Error in dangle data format!\n");
    exit(1);
  }
  
  H_ALPHA_1 = (DBL_TYPE) array[0]/100.0;
  H_ALPHA_2 = (DBL_TYPE) array[1]/100.0;
  H_ALPHA_3 = (DBL_TYPE) array[2]/100.0;
  
  em->alpha_1 = (G_ALPHA_1 - H_ALPHA_1)* temp_k/310.15 + 
    H_ALPHA_1;
  em->alpha_2 = (G_ALPHA_2 - H_ALPHA_2)* temp_k/310.15 + 
    H_ALPHA_2;
  em->alpha_3 = (G_ALPHA_3 - H_ALPHA_3)* temp_k/310.15 + 
    H_ALPHA_3;
  
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //AT PENALTY
  if( sscanf(line, "%d", &(array[0]) ) == 1) {
    H_AT_PENALTY = (DBL_TYPE) array[0]/100.0;
    em->at_penalty = (G_AT_PENALTY - H_AT_PENALTY)* temp_k/310.15 + 
      H_AT_PENALTY;
    
  }
  else {
    fprintf(stderr, "Error in AT PENALTY data\n");
    exit(1);
  }
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //1x1 interior loop
  for( i = 0; i < 36; i++) {
    fgets( line, MAXLINE, fp); //read in label
    for( j = 0; j < 4; j++) {
      
      nRead = 0;
      token = strtok( line, " ");
      while( token != NULL) {
        if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
          
          nRead++;
        }
        token = strtok( NULL, " ");
      }
      
      if( nRead != 4) {
        fprintf(stderr, "Error in 1x1 format!\n");
        exit(1);
      }
      
      for( k = 0; k < 4; k++) {
        H_IL_SInt2[ i*16 + j*4 + k] = (DBL_TYPE) array[ k]/100.0;
        em->IL_SInt2[i*16+j*4+k] = (G_IL_SInt2[i*16+j*4+k] - 
                                H_IL_SInt2[i*16+j*4+k])* temp_k/310.15 + 
          H_IL_SInt2[i*16+j*4+k];
        
      }
      
      fgets( line, MAXLINE, fp);
    }
  }
  
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //2x2 interior loop
  for( i = 0; i < 36*16; i++) {
    fgets( line, MAXLINE, fp); //read in label
    for( j = 0; j < 4; j++) {
      
      nRead = 0;
      token = strtok( line, " ");
      while( token != NULL) {
        if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
          
          nRead++;
        }
        token = strtok( NULL, " ");
      }
      
      if( nRead != 4) {
        fprintf(stderr, "Error in 1x1 format!\n");
        exit(1);
      }
      
      for( k = 0; k < 4; k++) {
        index4 = 1536*(i/96)+256*((i%96)/16) +
          64*((i%16)/4) + 4*(i%4) + k*16 + j;
        H_IL_SInt4[index4] = 
          (DBL_TYPE) array[ k]/100.0;
        em->IL_SInt4[index4] = 
          (G_IL_SInt4[index4] - 
           H_IL_SInt4[index4])* temp_k/310.15 + 
          H_IL_SInt4[index4];
        
      }
      
      fgets( line, MAXLINE, fp);
    }
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //1x2 interior loop
  for( i = 0; i < 144; i++) {
    fgets( line, MAXLINE, fp); //read in label
    for( j = 0; j < 4; j++) {
      
      nRead = 0;
      token = strtok( line, " ");
      while( token != NULL) {
        if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
          
          nRead++;
        }
        token = strtok( NULL, " ");
      }
      
      if( nRead != 4) {
        fprintf(stderr, "Error in 1x1 format!\n");
        exit(1);
      }
      
      for( k = 0; k < 4; k++) {
        index4 = 384*(i/24) + 4*((i%24)/4) + 24*(i%4) +
          96*j + k;
        H_IL_AsInt1x2[index4] = (DBL_TYPE) array[ k]/100.0;
        em->IL_AsInt1x2[index4] = 
          (G_IL_AsInt1x2[index4] - 
           H_IL_AsInt1x2[index4])* temp_k/310.15 + 
          H_IL_AsInt1x2[index4];
      }
      
      fgets( line, MAXLINE, fp);
    }
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //polyC hairpin parameters
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 3) {
    fprintf(stderr, "9. Error in dangle data format!\n");
    exit(1);
  }
  
  H_POLYC3 = (DBL_TYPE) array[0]/100.0;
  H_POLYCSLOPE = (DBL_TYPE) array[1]/100.0;
  H_POLYCINT = (DBL_TYPE) array[2]/100.0;
  
  em->polyc3 = (G_POLYC3 - H_POLYC3)* temp_k/310.15 + 
    H_POLYC3;
  em->polycslope = (G_POLYCSLOPE - H_POLYCSLOPE)* temp_k/310.15 + 
    H_POLYCSLOPE;
  em->polycint = (G_POLYCINT - H_POLYCINT)* temp_k/310.15 + H_POLYCINT;
  
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //Pseudoknot parameters
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 5) {
    fprintf(stderr, "Error in pseudoknot data format!\n");
    exit(1);
  }
  
  H_BETA_1 = (DBL_TYPE) array[0]/100.0;
  H_BETA_2 = (DBL_TYPE) array[1]/100.0;
  H_BETA_3 = (DBL_TYPE) array[2]/100.0;
  H_BETA_1M = (DBL_TYPE) array[3]/100.0;
  H_BETA_1P = (DBL_TYPE) array[4]/100.0;
  
  em->beta_1 = (G_BETA_1 - H_BETA_1)* temp_k/310.15 + 
    H_BETA_1;
  em->beta_2 = (G_BETA_2 - H_BETA_2)* temp_k/310.15 + 
    H_BETA_2;
  em->beta_3 = (G_BETA_3 - H_BETA_3)* temp_k/310.15 + 
    H_BETA_3;
  em->beta_1m = (G_BETA_1M - H_BETA_1M)* temp_k/310.15 + 
    H_BETA_1M;
  em->beta_1p = (G_BETA_1P - H_BETA_1P)* temp_k/310.15 + 
    H_BETA_1P;
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //BIMOLECULAR TERM
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 1) {
    fprintf(stderr, "Error in bimolecular data format!\n");
    exit(1);
  }
  
  H_BIMOLECULAR = (DBL_TYPE) array[0]/100.0;
  
  em->bimolecular = (G_BIMOLECULAR - H_BIMOLECULAR)* temp_k/310.15 + 
    H_BIMOLECULAR - kB*temp_k*LOG_FUNC( water_conc);
  
  fclose( fp);


  // Make the salt corrections
  // Stacked bases
  for (i = 0; i < 36; i++) {
    em->Stack[i] += em->salt_correction;
  }

  // Loop correction.  Covers all hairpins, and bulges, but overcounts for
  // 1-base bulges.  This is corrected in the function InteriorEnergyFull.
  // This also covers large (non-tabulated) interior loops.
  for (i = 0; i < 90; i++) {
    em->loop37[i] += em->salt_correction;
  }

  // Corrections for tabulated interior loops
  for (i = 0; i < 16*36; i++) {
    em->IL_SInt2[i] += em->salt_correction;
  }
  for (i = 0; i < 256*36; i++) {
    em->IL_SInt4[i] += em->salt_correction;
  }
  for (i = 0; i < 64*36; i++) {
    em->IL_AsInt1x2[i] += em->salt_correction;
  }

  // Multiloop
  em->alpha_1 += em->salt_correction;

}

/* ************** */

void setParametersToZero(energy_model_t *em) {
  
  int i;
	
  for( i = 0; i < 90; i++) {
    em->loop37[i] = 0;
  }
  for ( i = 0; i < 4096; i++) {
    em->tloop_energy[i] = 0;
  }
  for( i = 0; i < 2048; i++) {
    em->triloop_energy[i] = 0;
  }
  for( i = 0; i < 16*6; i++) {
    em->MMEnergiesHP[i] = 0;
  }
  for( i = 0; i < 256; i++) {
    em->MMEnergiesIL[i] = 0;
  }
  for( i = 0; i < 16*36; i++) {
    em->IL_SInt2[i] = 0;
  }
  for( i = 0; i < 256*36; i++) {
    em->IL_SInt4[ i] = 0;
  }
  for( i = 0; i < 64*36; i++) {
    em->IL_AsInt1x2[i] = 0;
  }
  for( i = 0; i < 48; i++) {
    em->dangle_energy[i] = 0;
  }
  for( i = 0; i < 4; i++) {
    em->asymmetry_penalty[i] = 0;
  }
  em->max_asymmetry = 0;
	
  for( i = 0; i < 36; i++) {
    em->Stack[i] = 0;
  }

  em->alpha_1 = em->alpha_2 = em->alpha_3 = em->beta_1 = em->beta_2 =
    em->beta_3 = em->beta_1m = em->beta_1p = em->polyc3 =
    em->polycint = em->polycslope = em->at_penalty = em-> bimolecular = 0;
}

/* *************** */
DEV
void InitEtaN( int **etaN, const int *nicks, int seqlength) {
  
  int i,j,k, nick;
  int indexE;
  
  for( i = 0; i <= seqlength-1; i++) {
    for( j = i-1; j <= seqlength-1; j++) {
      indexE = pf_index( i, j, seqlength);
      /*
      etaN[ indexE] = (int *) malloc( 2*sizeof( int));
      if(!etaN[indexE]) {
        printf("etaN allocation failed\n");
        assert(0);
      }
      */
      etaN[ indexE][0] = 0;
      etaN[ indexE][1] = -1;
      
    }
  }
  
  k = 0;
  nick = nicks[k];
  while( nick != -1) {
    for( i = 0; i <= nick; i++) {
      for( j = nick; j <= seqlength-1; j++) { 
        indexE =  pf_index(i,j,seqlength);
        etaN[ indexE][0]++;
        if( etaN[ indexE][1] == -1) { 
          //assume nicks are assigned in increasing order 
          etaN[indexE][1] = k;
        }
      }
    }
    nick = nicks[++k];
  }
}

/* *************************** */
int EtaNIndex_old( float i, float j, int seqlength) { 
  return pf_index( (int) i, (int) j, seqlength);
}

/* ******************* */


