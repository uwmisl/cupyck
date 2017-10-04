/*
  ene.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 3/2006, Justin Bois 1/2007

  This file contains energy functions used for determining energies
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "pfuncUtilsHeader.h"
#include "DNAExternals.h"

DEV DBL_TYPE ExplDangleRaw( int i, int j, int seq[], int seqlength);

/* ************************************** */
DEV
DBL_TYPE HelixEnergy( int i, int j, int h, int m, energy_model_t *em) {
  // Calculate the energy of the helical region closed by pair
  // i-j and h-m.  Data from Zuker's mfold file stack.dgd

  int shift_ij; // Type of base pair
  int shift_hm; // Type of base pair

  shift_ij = GetMismatchShift( i, j);
  shift_hm = GetMismatchShift( h, m);

  if( shift_ij < 4 && shift_hm < 4) {
    return em->Stack[ ( i - 1)*6 + (h - 1) ];
  }

  if( shift_ij < 4 && shift_hm >= 4) {
    return em->Stack[ (i - 1)*6 + (h + 1) ];
  }

  if( shift_ij >= 4 && shift_hm < 4) {
    return em->Stack[ (i + 1)*6 + (h - 1) ];
  }

  if( shift_ij >= 4 && shift_hm >= 4) {
    return em->Stack[ (i + 1)*6 + (h + 1) ];
  }
  else {
    printf("Error in HelixEnergy!");
    assert(0);
    return NAD_INFINITY; // This never is returned
  }
}

// *******************************************************************
DEV
DBL_TYPE InteriorMM( char a, char b, char x, char y, energy_model_t *em) {
/*
  Interior Mismatch calculation

  This calculates the Mismatch interaction energies between positions
  1 -> 5' a x 3'
  2 -> 3' b y 5'
  Interactions energies taken from file tstacki2.dgd.
*/

  int cp_shift;
  DBL_TYPE energy = 0.0;

  cp_shift = GetMismatchShift( a, b );
  energy = em->MMEnergiesIL[ (4*(( x) - 1) + (( y) - 1) )*6 + cp_shift];

  return energy;
}

/* ********************************************** */

DEV
DBL_TYPE HairpinEnergy( int i, int j, int seq[], energy_model_t *em) {

  // This gives the energy of the hairpion closed by bases i and j
  DBL_TYPE energy;  //energy of hairpin

  int triloopnumber; //Index for specific triloop
  int tloopnumber; //index for tloops

  int size; //Doesn't include closing pair i-j

  int cp_shift; //Classification of base-pair for energy mismatch

  int polyC = TRUE;  //Is the hairpin a poly-C?
  int k;
  for( k = i+1; k < j; k++) {
    if( seq[k] != BASE_C) {
      polyC = FALSE;
      break;
    }
  }

  size = j - i - 1;

  if( size < 3) {
    return NAD_INFINITY;
  }

  if( CanPair( seq[i], seq[j]) == FALSE ) {
    return NAD_INFINITY;
  }

  if( size <= 30) {
    energy = em->loop37[ 60 + size - 1];
  }
  else {
    energy = em->loop37[ 60 + 30 - 1];
    energy += sizeLog (size, em); //1.75*kB*TEMP_K*LOG_FUNC( size/30.0);

    if( em->dnarnacount == COUNT) {
      energy = 0;
    }

  }

  if( size == 3) {
    //Get Triloop energy

    if( seq[i] != BASE_C && seq[j] != BASE_C) {
      energy += em->at_penalty;
    }

    triloopnumber = 256*(( seq[i]) - 1) +
      64*(( seq[i + 1]) - 1) +
      16*(( seq[i + 2]) - 1) +
      4*( ( seq[j - 1]) - 1) +
      1*( ( seq[j]) - 1);

    // 0 mismatch energy for triloops
    energy += em->triloop_energy[ triloopnumber];

    //Poly-C loop
    if( polyC == TRUE) {
      energy += em->polyc3;
    }
  }
  else if (size == 4) {

    tloopnumber = 1024*(( seq[i]) - 1) +
      256*(( seq[i + 1]) - 1) +
      64*( ( seq[i + 2]) - 1) +
      16*( ( seq[j - 2]) - 1) +
      4*(  ( seq[j - 1]) - 1) +
      1*(  ( seq[j])- 1);
    energy +=  em->tloop_energy[ tloopnumber];

    //Next do mismatches.
    cp_shift = GetMismatchShift( seq[i], seq[j]);

    energy += em->MMEnergiesHP[(4*(( seq[i + 1]) - 1) +
                            (( seq[j - 1]) - 1) )*6
                           + cp_shift];
    //Poly-C loop
    if( polyC == TRUE) {
      energy += em->polycslope*size + em->polycint;
    }
  }

  else if (size > 4) {
    // Calculate mismatch
    cp_shift = GetMismatchShift( seq[i], seq[j]);

    energy += em->MMEnergiesHP[(4*(( seq[i + 1]) - 1) +
                            (( seq[j - 1]) - 1) )*6
                           + cp_shift];

    //Poly-C loop
    if( polyC == TRUE) {
      energy += em->polycslope*size + em->polycint;
    }
  }
  return energy;
}



/* ****************************************** */
DEV
DBL_TYPE InteriorEnergy(  int i, int j, int h, int m, int seq[], energy_model_t
    * em) {
  return InteriorEnergyFull( i, j, h, m, seq, TRUE, em);
}

DEV
DBL_TYPE InteriorEnergyFull( int i, int j, int h, int m, int seq[],
                             int calcIJ, energy_model_t *em) {

  DBL_TYPE energy = 0.0;
  int L1, L2; //lengths of the 2 single stranded regions
  int size;
  int asymmetry;
  int cp_shift, ip_shift;  // For classifying basepairs

  if( em->dnarnacount == COUNT) return 0;
#ifdef DEBUG
  if( i >= h || h >= m || m >= j) {
    printf("Invalid boundary to interior loop! %d %d %d %d\n", i, h, m, j);
    assert(0);
  }
#endif

  L1 = h - i - 1;
  L2 = j - m - 1;
  size = L1 + L2;

  if( size == 0) { //Helical region
    energy = HelixEnergy( seq[i], seq[j], seq[h], seq[m], em);
  }

  else if ( L1*L2 == 0) { //Bulge
    if( size <= 30) {
      energy = em->loop37[ 30 + size - 1];
    }
    else {
      energy = em->loop37[ 30 + 30 - 1];
      energy += sizeLog (size, em); //1.75*kB*TEMP_K*LOG_FUNC( size/30.0);
    }

    if( L1 + L2 == 1 ) { //single bulge...treat as a stacked region
      energy += HelixEnergy( seq[i], seq[j], seq[h], seq[m], em);
      energy -= em->salt_correction;  // Correct for the extra salt correction
                                 // added from the HelixEnergy
    }
    else {
      // Next do AT_Penalty for no GC termination, assuming size >= 2
      if( seq[i] != BASE_C && seq[j] != BASE_C) {
        energy += em->at_penalty;
      }
      if( seq[h] != BASE_C && seq[m] != BASE_C) {
        energy += em->at_penalty;
      }
    }
  }
  else if ( L1 > 0 && L2 > 0) {
    asymmetry = abs( L1 - L2);
    if( asymmetry > 1 || size > 4) { //data not tabulated

      energy = asymmetryEfn( L1, L2, size, em);

      //Stacking Energy
      if( L1 > 1 && L2 > 1) { //Non-GAIL Version
        energy += InteriorMM( seq[m], seq[h], seq[m+1], seq[h-1], em);

        if( calcIJ == TRUE)
          energy += InteriorMM( seq[i], seq[j], seq[i+1], seq[j-1], em);
      }
      else if( L1 == 1 || L2 == 1) {// GAIL =>assume AA terminal mismatch
#ifndef NO_GAIL
        energy +=
          InteriorMM( seq[m], seq[h], BASE_A, BASE_A, em);
        if( calcIJ == TRUE)
          energy += InteriorMM( seq[i], seq[j], BASE_A, BASE_A, em);
#else
        energy += InteriorMM( seq[m], seq[h], seq[m+1], seq[h-1], em);
        if( calcIJ == TRUE)
          energy += InteriorMM( seq[i], seq[j], seq[i+1], seq[j-1], em)
#endif
      }
      else {
        printf("Error: Unclassified interior loop!\n");
        assert(0);
      }
    }
    else { //get tabulated data
      if( asymmetry == 0 && size == 2) {
        cp_shift = GetMismatchShift( seq[i], seq[j]);
        ip_shift = GetMismatchShift( seq[h], seq[m]);
        if (cp_shift==-1 || ip_shift==-1) return 0.0; //Wrongly called
        energy += em->IL_SInt2[ 96*cp_shift + 16*ip_shift +
                           4*(( seq[i+1]) - 1) +
                           (( seq[ j -1]) - 1) ];
      }
      else if( asymmetry == 0 && size == 4) {
        cp_shift = GetMismatchShift( seq[i], seq[j]);
        ip_shift = GetMismatchShift( seq[h], seq[m]);
        if (cp_shift==-1 || ip_shift==-1) return 0.0; //Wrongly called
        energy += em->IL_SInt4[ cp_shift*256*6 +  ip_shift*256 +
                           (4*(( seq[ i+1])  - 1) +
                            ( seq[ j - 1])   - 1)*16 +
                           (4*( ( seq[ i+2]) - 1) +
                            ( seq[ j - 2])   - 1) ];
      }
      else if( asymmetry == 1 && L1 == 1) {
        cp_shift = GetMismatchShift( seq[i], seq[j]);
        ip_shift = GetMismatchShift( seq[h], seq[m]);
        if (cp_shift==-1 || ip_shift==-1) return 0.0; //Wrongly called
        energy += em->IL_AsInt1x2[ cp_shift*4*24*4 +
                              (( seq[ j - 2]) - 1)*24*4 +
                              (( seq[ i + 1]) - 1)*24 +
                              4*ip_shift +
                              ((( seq[ j - 1]) - 1) % 4) ];
      }
      else if( asymmetry == 1 && L1 == 2) {
        cp_shift = GetMismatchShift( seq[j], seq[i]);
        ip_shift = GetMismatchShift( seq[m], seq[h]);
        if (cp_shift==-1 || ip_shift==-1) return 0.0; //Wrongly called
        //note reversed order of inputs above.
        //This is to comply with the format of asint1x2

        energy += em->IL_AsInt1x2[ ip_shift*4*24*4 +
                              (( seq[i + 1]) - 1)*24*4 +
                              (( seq[j - 1]) - 1)*24 +
                              4*cp_shift +
                              ((( seq[i + 2]) - 1) % 4) ];
      }
      else {
        printf("Error in tabulated Interior Loop!\n");
        assert(0);
      }
    }
  }
  else {
    printf("Improperly classified Interior Loop!\n");
    assert(0);
  }

  return energy;
}


/* ******************************** */
DEV
DBL_TYPE DangleEnergy( int i, int j, int seq[], int seqlength, energy_model_t
    *em) {
  //0 energy except for dangles

  DBL_TYPE dangle5 = 0;
  DBL_TYPE dangle3 = 0;
  int dangle_shift;

  if( em->dangletype != 2) {
    if( j == i - 1) {
      return 0;
    }
  }
  else if( j == i - 1 && (i == 0 || j == seqlength - 1) ) {
    return 0;
  }

  if( j == seqlength - 1) {
    dangle3 = 0;
  }
  else {
    int pt=GetPairType( seq[ j + 1]);
    if (pt==-1) {
      printf("i=%d j=%d seq[%d]=%d\n",i,j,j-1,seq[j-1]);
      assert(0);
    }
    dangle_shift = 3 - pt;
    dangle3 = em->dangle_energy[ 24 + dangle_shift*4 +
                            ( seq[ j]) - 1];
  }

  if( i == 0) {
    dangle5 = 0;
  }
  else {
    int pt=GetPairType( seq[ i - 1]);
    if (pt==-1) {
      printf("i=%d j=%d seq[%d]=%d\n",i,j,i-1,seq[i-1]);
      assert(0);
    }

    dangle_shift = pt;
    dangle5 = em->dangle_energy[ dangle_shift*4 +
                            ( seq[ i]) - 1];
  }

  if( em->dangletype != 2 && i == j && i != 0 && j != seqlength - 1) {
    return MIN(dangle3, dangle5 );
  }

  return dangle3 + dangle5;
}

/* ******************************** */
DEV
DBL_TYPE ExplDangleRaw( int i, int j, int seq[], int seqlength, energy_model_t
    *em) {
  //0 energy except for dangles

  DBL_TYPE dangle5 = 0;
  DBL_TYPE dangle3 = 0;
  int dangle_shift;

  if( (j == i - 1) || (j==-1 && i>0)) {
    return 1.0;
  }
  if( (j==-1 && i>0) || (j == i - 1 && (i == 0 || j == seqlength - 1)) ) {
    return 1.0;
  }

  if( j == seqlength - 1) {
    dangle3 = 0;
  }
  else {
    dangle_shift = 3 - GetPairType( seq[ j + 1]);
    dangle3 = em->dangle_energy[ 24 + dangle_shift*4 +
                            seq[ j] - 1];
  }

  if( i == 0) {
    dangle5 = 0;
  }
  else {
    dangle_shift = GetPairType( seq[i-1]);
    dangle5 = em->dangle_energy[ dangle_shift*4 +
                            seq[ i] - 1];
  }

  if(i == j && i != 0 && j != seqlength - 1) {
    if (em->dangletype == 2) return EXP_FUNC(-(dangle5 + dangle3)/(em->temp_k*kB));
    return EXP_FUNC( -MIN(dangle3, dangle5)/(em->temp_k*kB) );
  }

  return EXP_FUNC( -(dangle3 + dangle5)/(em->temp_k*kB) );
}

DEV
DBL_TYPE ExplDangle( int i, int j, int seq[], int seqlength, energy_model_t *em) {
  return ExplDangleRaw(i,j,seq,seqlength, em);
}


/* *********** */
DEV
DBL_TYPE NickDangle(int i, int j, const int *nicks, int **etaN, int hairpin,
                    int seq[], int seqlength, energy_model_t *em) {

  DBL_TYPE dangle5 = 0;
  DBL_TYPE dangle3 = 0;
  int dangle_shift;
  int nick;
  int nIndex;

  nick = -5;

  if( i != 0) { //if j == seqlength -1, this is still OK
    nIndex = EtaNIndex( i-0.5, j+0.5, seqlength);
  }
  else {
    nIndex = EtaNIndex( i+0.5, j+0.5, seqlength);
  }

  if( etaN[ nIndex][0] >= 2 ||
     ( etaN[ nIndex][0] == 1 && (i == 0 || j == seqlength -1)) ) {

       return NAD_INFINITY;
     }
  else if( etaN[ nIndex][0] >= 1) {
    nick = nicks[ etaN[ nIndex][1]];
  }

  if( em->dnarnacount == COUNT)
    return 0;

  if( j == i - 1) {
    return 0;
  }
  if( j == i - 1 && (i == 0 || j == seqlength - 1) ) {
    return 0;
  }

  if( j == seqlength - 1 || j == nick) {
    dangle3 = 0;
  }
  else {
    if( hairpin == FALSE) {
      dangle_shift = 3 - GetPairType( seq[ j + 1]);
    }
    else {
      dangle_shift = GetMismatchShift( seq[i-1], seq[j+1]);
    }

    dangle3 = em->dangle_energy[ 24 + dangle_shift*4 +
                            ( seq[ j]) - 1];
  }

  if( i == 0 || i-1 == nick) {
    dangle5 = 0;
  }
  else {
    if( hairpin == FALSE) {
      dangle_shift = GetPairType( seq[i-1]);
    }
    else {
      dangle_shift = GetMismatchShift( seq[i-1], seq[j+1]);
    }

    dangle5 = em->dangle_energy[ dangle_shift*4 +
                            ( seq[ i]) - 1];
  }

  if( nick >= i-1 && nick <= j) {
    return dangle3 + dangle5;
  }
  else {
    if( j > i || j == seqlength - 1 || i == 0) {
      return dangle3 + dangle5;
    }
    if(j == i && i != 0 && j != seqlength - 1) {
      if (em->dangletype == 2) return dangle3 + dangle5;
      return MIN(dangle3, dangle5);
    }
    //j == i-1 already handled above
  }
  printf("Error with function: NickDangle\n");
  assert(0);
  return -1; //Error!  This should never happen
}

/* ************** */
DEV
DBL_TYPE NickedEmptyQ( int i, int j, const int nicks[], int seq[],
                      int seqlength, int **etaN, energy_model_t *em) {

  if( j <= i || etaN[ EtaNIndex( i+0.5, j-0.5, seqlength)][0] == 0) {
    return EXP_FUNC( -1*NickDangle(i, j, nicks, etaN,
                              FALSE, seq, seqlength, em)
                /(kB*em->temp_k));
  }
  else { //disconnected
    return 0;
  }

}

/* ********* */
DEV
DBL_TYPE ExplInternal( int i, int j, int h, int m, int seq[], energy_model_t *em) {
  // Calculates E^(-energy/RT) of interior loop closed by i-j and h-m

  DBL_TYPE energy = InteriorEnergy( i, j, h, m, seq, em);
  if( energy == NAD_INFINITY) {
    return 0.0;
  }
  return EXP_FUNC( - energy/( kB*em->temp_k));
}

DEV
DBL_TYPE sizeLog(int size, energy_model_t *em){
  return 1.75*kB*em->temp_k*LOG_FUNC(size/30.0);
}


/* ******* */
DEV
DBL_TYPE asymmetryEfn( int L1, int L2, int size, energy_model_t *em) {

  int asymmetry_index;
  DBL_TYPE energy;
  int asymmetry = abs( L1 - L2);

  //Loop Size Energy
  if( size <= 30) {
    energy = em->loop37[ size - 1];
  }
  else {
    energy = em->loop37[ 30 - 1];
    energy += sizeLog(size, em);
  }

  //Asymmetry rountine copied from efn.f in Zuker's mfold package.
  asymmetry_index = 4;

  if( L1 < asymmetry_index) {
    asymmetry_index = L1;
  }

  if( L2 < asymmetry_index) {
    asymmetry_index = L2;
  }

  if( asymmetry*em->asymmetry_penalty[ asymmetry_index - 1] < em->max_asymmetry ) {
    energy += asymmetry*em->asymmetry_penalty[ asymmetry_index - 1];
  }
  else {
    energy += em->max_asymmetry; // MAX asymmetry penalty
  }
  return energy;
}



/* ********************** */
