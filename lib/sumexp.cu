/*
  sumexp.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 7/2006, Justin Bois 1/2007

  This file contains the functions that calculate the sum of exponentials
  in order to calculate a given partition function matrix. 

  See pfuncUtilsHeader.h for more specific descriptions of each function.
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<string.h>

#include "pfuncUtilsHeader.h" //contains functions and structures
/* ******************************************* */
DEV
DBL_TYPE ExplHairpin( int i, int j, int seq[], int seqlength, int **etaN,
    energy_model_t *em) {
  //this version disallows nicks here

  DBL_TYPE energy = 0;
  int index;
  int nNicks = 0;

  index = EtaNIndex(i+0.5, j-0.5, seqlength);
  nNicks = etaN[ index][0];

  if( nNicks >= 1) return 0;

  if( nNicks == 0 && j-i <= 3) {
    return 0;
  }

  if( nNicks == 0) {
    energy = HairpinEnergy( i, j, seq, em);
  }

  if( energy == NAD_INFINITY) {
    return 0.0;
  }

  return EXP_FUNC( -energy/( kB*em->temp_k) );
}

/* ********************* */
DEV
DBL_TYPE SumExpMultiloops( int i, int j, int seq[], 
                          DBL_TYPE *Qms, DBL_TYPE *Qm, int seqlength,
                          int **etaN, energy_model_t *em){
  // Decomposes the region inside pair i,j into multiloops, i.e.
  // and excludes the possibility of "top level" nicks

  DBL_TYPE sum_exp = 0.0;
  DBL_TYPE bp_penalty = 0.0;
  DBL_TYPE extraTerms;

  int d; // d is the left base of a rightmost paired base between i, j.
 
  if( ( seq[i]) + ( seq[j]) == 5) {
    for( d = i+3; d <= j - 2; d++) {
      //reset loop parameters
      bp_penalty = 0.0;
      if( etaN[EtaNIndex_same( d-0.5, seqlength)][0] == 0 ) {
        
        if( seq[i] != BASE_C  && seq[j] != BASE_C) {
          bp_penalty += em->at_penalty;
        }

        extraTerms = EXP_FUNC( -( em->alpha_1 + em->alpha_2 + bp_penalty) 
                          / (kB*em->temp_k) );
        if( em->dnarnacount == COUNT) 
          extraTerms = 1;

        sum_exp += Qm[ pf_index( i+1, d-1, seqlength)] *
          Qms[ pf_index(d, j-1, seqlength)] * extraTerms;
      }
      
    }
  }

  return sum_exp;  
}

/* *********************************************** */
DEV
DBL_TYPE SumExpExteriorLoop( int i,int j, int seq[], int seqlength, 
                            DBL_TYPE *Q, int *nicks, int **etaN, energy_model_t
                            *em) {

  DBL_TYPE sumExp = 0.0;
  DBL_TYPE bp_penalty = 0.0;
  int multiNick = -1;
  int index_ij;
  int leftIndex;
  int nNicks;
  int n;
  int iNicked, jNicked;
  DBL_TYPE extraTerms;

  index_ij = EtaNIndex(i+0.5, j-0.5, seqlength);
  iNicked = jNicked = FALSE;

  if( etaN[ EtaNIndex_same(j-0.5, seqlength)][0] != 0) {
    jNicked = TRUE;
  }

  if( etaN[ EtaNIndex_same(i+0.5, seqlength)][0] != 0) {
    iNicked = TRUE;
  }

  if( ( seq[i]) + ( seq[j]) == 5) {
    bp_penalty = 0.0;

    if( seq[i] != BASE_C  && seq[j] != BASE_C) {
      bp_penalty = em->at_penalty;
    }

    nNicks = etaN[ index_ij][0];
    leftIndex = etaN[ index_ij ][1];

    //treat each nick as rightmost nick in paired interval
    //All non-nicked cases handled outside this loop
    for( n = 0; n <= nNicks-1; n++) {
      multiNick = nicks[ leftIndex + n];

      extraTerms = 
        EXP_FUNC( -1*(bp_penalty)/(kB*em->temp_k));

      if( em->dnarnacount == COUNT) 
        if( extraTerms != 0) extraTerms = 1;

      if( (iNicked == FALSE && jNicked == FALSE) ||
         (i == j - 1) || 
         (multiNick == i && jNicked == FALSE) ||
         (multiNick == j-1 && iNicked == FALSE ) ) {

           sumExp +=
             Q[ pf_index(i+1, multiNick, seqlength)]*
             Q[ pf_index( multiNick+1, j-1, seqlength)] *
             extraTerms;
      }
    }
  }
  return sumExp;  
}

/* *********************************************** */
DEV
void fastILoops( int i, int j, int L, int seqlength, int seq[],
                 int **etaN, DBL_TYPE *Qb, DBL_TYPE *Qx, DBL_TYPE *Qx_2,
                 energy_model_t *em) {

  int size;
  int pf_ij = pf_index( i, j, seqlength);
  DBL_TYPE extraTerms;

  int isEndNicked = FALSE;
  if( etaN[ EtaNIndex( i-0.5,i-0.5, seqlength)][0] == 1 || 
     etaN[ EtaNIndex( j+0.5,j+0.5, seqlength)][0] == 1) 
    isEndNicked = TRUE;
  if( L >= 12) {
    makeNewQx( i, j, seq, seqlength, etaN, Qb, Qx, em);
  }

  //Use extensible cases              
  if( CanPair( seq[ i], seq[j]) == TRUE) {
    for( size = 8; size <= L - 4; size++) {
      
      extraTerms = EXP_FUNC( -InteriorMM( seq[i], seq[j], seq[i+1], 
                                     seq[j-1], em)/(kB*em->temp_k));
      if( em->dnarnacount == COUNT) 
        extraTerms = 1;
      
      Qb[ pf_ij] += 
        Qx[ fbixIndex( j-i, i, size, seqlength)] * 
        extraTerms;
    }
  }

  if( L >= 12 && i != 0 && j != seqlength -1 && isEndNicked == FALSE) {
    extendOldQx( i, j, seqlength, Qx,Qx_2, em );
  }
  
  /* Add in inextensible cases */  
  if( CanPair( seq[ i], seq[j]) == TRUE) {
    //first check inextensible cases

    Qb[ pf_ij] += SumExpInextensibleIL( i,j, seq, seqlength, Qb,  etaN, em);

  } 
}

/* *************** */

/* Qs, Qms  Recursion */
DEV
void MakeQs_Qms( int i, int j, int seq[], int seqlength, 
                DBL_TYPE *Qs, DBL_TYPE *Qms, DBL_TYPE *Qb,
                int *nicks, int **etaN, energy_model_t *em) {

  int d; //base pair is i,d
  DBL_TYPE bp_penalty = 0.0;
  int pf_ij = pf_index( i, j, seqlength);

  DBL_TYPE extraTerms;
  int nNicks;
  int index_ij = EtaNIndex( i+0.5, j-0.5, seqlength);
  int start;

  nNicks = etaN[ index_ij][0];
  if( nNicks >= 1) {
    start = nicks[ etaN[ index_ij][1] + nNicks - 1]+1;
  }
  else {
    start = i+4;
  }

  for( d = start; d <= j; d++) {
    bp_penalty = 0.0;
    
    if( CanPair( seq[i], seq[ d]) == TRUE &&
       ( seq[i]) + ( seq[d]) == 5) {
         
         if( seq[i] != BASE_C && seq[d] != BASE_C) {
           bp_penalty = em->at_penalty;
         }

         extraTerms = EXP_FUNC( -(NickDangle( d+1,j,nicks, etaN,
                              FALSE, seq,seqlength, em) + 
                              bp_penalty)/(kB*em->temp_k) );

         if( em->dnarnacount == COUNT) 
           extraTerms = 1;

         Qs[ pf_ij] += Qb[ pf_index( i, d, seqlength) ] * 
           extraTerms;

         // ******************** 

         extraTerms =  ExplDangle( d+1, j, seq, seqlength, em) * 
           EXP_FUNC( -(bp_penalty + em->alpha_2 + 
                 em->alpha_3*(j-d))/(kB*em->temp_k) );

         if( em->dnarnacount == COUNT) 
           extraTerms = 1;
         Qms[ pf_ij] += Qb[ pf_index( i, d, seqlength) ] * 
           extraTerms;
    }
  }
}


/* ******************************* */
/* Q, Qm Recursions */
DEV
void MakeQ_Qm_N3( int i, int j, int seq[], int seqlength, 
                 DBL_TYPE *Q, DBL_TYPE *Qs, 
                 DBL_TYPE *Qms, DBL_TYPE *Qm,
                 int *nicks, int **etaN, energy_model_t *em) {
  int d;//left base of rightmost base pair.
  int pf_ij = pf_index( i, j, seqlength);

  DBL_TYPE extraTerms;
  Q[ pf_ij] = NickedEmptyQ( i, j, nicks, seq, seqlength, etaN, em);

  for( d = i; d <= j - 1; d++) {
    if( etaN[ EtaNIndex_same(d-0.5, seqlength)][0] == 0 || d == i ) {

      Q[ pf_ij] += Q[ pf_index(i, d-1, seqlength)] *
        Qs[ pf_index( d, j, seqlength)];

      if( em->dnarnacount == COUNT) 
        extraTerms = 1;
      else 
        extraTerms = ExplDangle( i, d-1, seq, seqlength, em) *
          EXP_FUNC( -(em->alpha_3)*(d-i)/(kB*em->temp_k));

      if( etaN[ EtaNIndex_same( d-0.5, seqlength)][0] == 0) { 
        //otherwise Qm not possible
       
        if( etaN[ EtaNIndex(i+0.5, d-0.5, seqlength)][0] == 0 ) {
          Qm[ pf_ij] += Qms[ pf_index( d, j, seqlength)] *
           extraTerms; //Single Pair
        }

        if( d >= i+2) {
         Qm[ pf_ij]+= Qm[ pf_index( i, d - 1, seqlength) ] *
           Qms[ pf_index( d, j, seqlength) ];
        }
      }
    }
  }
}

/* ******************************************* */

/* Functions in Q recursion */
// must be calculated after Qb, Qpk of same length

DEV void makeNewQx( int i, int j, int seq[], int seqlength, 
               int **etaN, DBL_TYPE Qb[], DBL_TYPE Qx[], energy_model_t *em) {
                 
  /*Determine the new entries of Qx(i,j,size) that are not extended 
  versions of Qx(i+1, j-1, size-2) */

  DBL_TYPE energy;
  int d, e; //Internal pair.(d, e will be restricted to special cases)

  int size, L1, L2; //size parameters: L1 + L2 = size, L1 = h-i-1, L2 = j-m-1

  //Add in all the cases that are not an extended version of an
  //extensible case.

  //Case 1:  L1 = 4, L2 >= 4;
  L1 = 4;
  d = i + L1 + 1;

  for( L2 = 4; L2 <= j - d - 2; L2++) {
    size = L1 + L2;
    e = j - L2 - 1;

    if( CanPair( seq[d], seq[e]) == TRUE &&
      (etaN[ EtaNIndex(i+0.5, d-0.5,seqlength)][0] == 0) &&
      (etaN[ EtaNIndex(e+0.5, j-0.5,seqlength)][0] == 0) ) {

        energy = asymmetryEfn( L1, L2, size, em) + InteriorMM( seq[e], seq[d],
            seq[e+1], seq[d-1], em);
        /*Exclude the i-j stacking energy here, just in case i-j 
        don't pair */

        if( em->dnarnacount == COUNT) 
          energy = 0;

        Qx[ fbixIndex( j-i, i, size, seqlength) ] += 
          EXP_FUNC(-energy/(kB*em->temp_k))*Qb[ pf_index(d, e, seqlength)];
    }
  }

  //Case 2  L1 > 4, L2 = 4
  L2 = 4;
  e = j - L2 -1;
  for( L1 = 5; L1 <= e-i-2; L1++) {
    size = L1 + L2;
    d = i + L1 + 1;

    if( CanPair( seq[d], seq[e]) == TRUE &&
      (etaN[ EtaNIndex(i+0.5, d-0.5,seqlength)][0] == 0) &&
      (etaN[ EtaNIndex(e+0.5, j-0.5,seqlength)][0] == 0) ) {

        energy = asymmetryEfn( L1, L2, size, em) + InteriorMM( seq[e], seq[d],
            seq[e+1], seq[d-1], em);
        /*Exclude the i-j stacking energy here, just in case i-j 
        don't pair */

        if( em->dnarnacount == COUNT) 
          energy = 0.0;

        Qx[ fbixIndex( j-i, i, size, seqlength)] +=
          EXP_FUNC(-energy/(kB*em->temp_k))*Qb[ pf_index(d, e, seqlength)];
    }
  }
}

/* ************************** */
DEV
void extendOldQx( int i, int j, int seqlength, DBL_TYPE Qx[], DBL_TYPE Qx_2[],
    energy_model_t *em) {
  /* Extends all entries of Qx */
  
  int size;
  DBL_TYPE oldSizeEnergy;
  DBL_TYPE newSizeEnergy;
  
  for( size = 8; size <= (j - i + 1) - 4; size++) {
    if( size <= 30) {
      oldSizeEnergy = em->loop37[ size - 1];
    }
    else {
      oldSizeEnergy = em->loop37[ 30 - 1];
      oldSizeEnergy += sizeLog(size, em);  //1.75*kB*TEMP_K*log( size/30.0); 
    }

    if( size + 2 <= 30) {
      newSizeEnergy = em->loop37[ size+2 - 1];
    }
    else {
      newSizeEnergy = em->loop37[ 30 - 1];
      newSizeEnergy += sizeLog (size+2, em); //1.75*kB*TEMP_K*log( (size+2)/30.0); 
    }

    if( em->dnarnacount == COUNT) 
      newSizeEnergy = oldSizeEnergy;
    
    Qx_2[ fbixIndex( j-i+2, i-1, size+2, seqlength)] = 
      Qx[ fbixIndex( j-i, i, size, seqlength)] * 
      EXP_FUNC( -(newSizeEnergy - oldSizeEnergy)/(kB*em->temp_k));
  }
}


/* ************************ */
DEV
DBL_TYPE SumExpInextensibleIL( int i, int j, int seq[], int seqlength, 
                              DBL_TYPE Qb[], int **etaN, energy_model_t *em) {
  /* This finds the minimum energy IL that has a special energy 
  calculation, i.e. small loops, bulge loops or GAIL case.  None of 
  these are allowed to be nicked
  */

  DBL_TYPE energy;
  DBL_TYPE sumexp = 0.0;
  int d, e; //Internal pair.(h, m will be restricted to special cases)  
  int L1, L2; //size parameters: L1 + L2 = size, L1 = h-i-1, L2 = j-m-1

  /* Consider "small" loops with special energy functions */

  for( L1 = 0; L1 <= 3; L1++) {
    d = i + L1 + 1;
    for( L2 = 0; L2 <= MIN( 3, j-d-2); L2++) {
      e = j - L2 - 1;

      if( CanPair( seq[d], seq[e]) == TRUE &&
         (etaN[ EtaNIndex(i+0.5, d-0.5,seqlength)][0] == 0) &&
         (etaN[ EtaNIndex(e+0.5, j-0.5,seqlength)][0] == 0) ) {

           energy = InteriorEnergy( i, j, d, e, seq, em);

           sumexp += EXP_FUNC( -energy/(kB*em->temp_k)) *
             Qb[ pf_index( d, e, seqlength)];
      }
    }
  }

  /* Next consider large bulges or large asymmetric loops */
  // Case 2a  L1 = 0,1,2,3, L2 >= 4;
  for( L1 = 0; L1 <= 3; L1++) {
    d = i + L1 + 1;
    for( L2 = 4; L2 <= j - d - 2; L2++) {
      e = j - L2 - 1;

      if( CanPair( seq[d], seq[e]) == TRUE &&
         (etaN[ EtaNIndex(i+0.5, d-0.5,seqlength)][0] == 0) &&
         (etaN[ EtaNIndex(e+0.5, j-0.5,seqlength)][0] == 0) ) {

           energy = InteriorEnergy( i, j, d, e, seq, em);

           sumexp += EXP_FUNC( -energy/(kB*em->temp_k)) *
             Qb[ pf_index( d, e, seqlength)]; 
      }
    }
  }

  // Case 2b L1 >= 4, L2 = 0,1,2,3;
  for( L2 = 0; L2 <= 3; L2++) {
    e = j - L2 - 1;
    for( L1 = 4;  L1 <= e - i - 2; L1++) {
      d = i + L1 + 1;

      if( CanPair( seq[d], seq[e]) == TRUE &&
         (etaN[ EtaNIndex(i+0.5, d-0.5,seqlength)][0] == 0) &&
         (etaN[ EtaNIndex(e+0.5, j-0.5,seqlength)][0] == 0) ) {

           energy = InteriorEnergy( i, j, d, e, seq, em);

           sumexp += EXP_FUNC( -energy/(kB*em->temp_k)) *
             Qb[ pf_index( d, e, seqlength)];
      }
    }
  }

  return sumexp;
}




