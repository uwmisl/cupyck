/*
  pfuncUtils.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 7/2001, Justin Bois 1/2007

  Generic utility functions used in partition function calculations.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <ctype.h>

#include "pfuncUtilsHeader.h" //contains functions and structures

/* ********************************************** */
//Added 08/21/2001



/* ******************************************** */
int pf_index_old( int i, int j, int N) {
  /* Calculate index for partition function array 
  N = sequence length */

  static int seqlength  = -1;
  static int storeValues[(MAXSEQLENGTH)*(MAXSEQLENGTH)];
  static int ss2;

  int ind,value;

  if( seqlength == N) {
    ind = i*N+j;
    if( j != i-1) {
      value=storeValues[ind];
      if( value != -1) {
        return value;
      }
      else storeValues[ind] = value = ss2 - (N-i)*(N-i+1)/2+(j-i);
      return value;
    }
    else {
      return ss2 + i;
    }
  }
  else {
    seqlength = N;
    ss2 = N*(N+1)/2;
    int maxN=N*N;
    for( ind = 0; ind < maxN; ind++) storeValues[ind] = -1;
    if( j != i-1) {
      ind =i*N+j;
      storeValues[ ind] = ss2 - (N-i)*(N-i+1)/2+(j-i);
      return storeValues[ ind];
    }
    else {
      return ss2 + i;
    }
  }

  /*
  if( j == i - 1) {
  return N*(N+1)/2 + i;
  }
  
  if( i < 0 || j > N - 1 || j < i) {
  printf("Illegal partition function index!\n i = %d, j = %d, ind = %d\n",
  i, j, N*(N+1)/2 - (N - i)*( N - i + 1)/2 + (j - i));
  exit(1);
  }
  
  return N*(N+1)/2 - (N - i)*( N - i + 1)/2 + (j - i);
  //return (N+1-i)*(N-i)/ 2 + N+1-j;
  */
} 



/* **************************************************** */
int Base2int( char base) {
  /* Convert A, C, G, T, U to 1,2,3,4 respectively */
  if( base == 'A') {
    return BASE_A;
  }
  else if ( base == 'C') {
    return BASE_C;
  }
  else if (base == 'G') {
    return BASE_G;
  }
  else if (base == 'T' || base =='U') {
    return BASE_U;
  }
  else if (base == '+') {
    return STRAND_PLUS;
  }
  else {
    fprintf(stderr, "Error in Converting base %c!\n", base);
    exit(1);
    return NAD_INFINITY; // never returns this
  }
}

/* **************************************************** */
char Int2base( int base) {
  /* Convert 1,2,3,4,5 to A, C, G, T, U respectively */
  if( base == BASE_A) {
    return 'A';
  }
  else if ( base == BASE_C) {
    return 'C';
  }
  else if (base == BASE_G) {
    return 'G';
  }
  else if (base == BASE_T) {
    return 'T';
  }
  else if (base == BASE_U) {
    return 'U';
  }
  else if (base == STRAND_PLUS) {
    return '+';
  }
  else {
    fprintf(stderr, "Error in Converting base %d!\n", base);
    exit(1);
    return '\0'; // never returns this
  }
}


/* *********************************************************** */ 
int convertSeq(char seqchar[], int seqnum[], int seqlength){
  int i;
  for (i=0; i<seqlength; i++)
    seqnum[i]=Base2int(toupper(seqchar[i]));
  seqnum[i]=-1;
  return -1;
}

/* *********************************************************** */ 
int printSeqNum(int seqnum[]){
  int i;
  for (i=0; seqnum[i]>0 && seqnum[i]<6; i++)
    printf("%c",Int2base(seqnum[i]));
  printf("\n");
  return -1;
}

/* *********************************************************** */ 
DEV
int GetMismatchShift( int base1, int base2) {
  /* base1 and base2 are basepaired. the returned value is needed to 
  index energy arrays
  */
  int shift;
  
  if( base1 == BASE_A) { /* this is for retrieving 
  proper mismatch energy */
    shift = 0;
  }
  else if( base1 == BASE_C) {
    shift = 1;
  }
  else if( base1 == BASE_G && base2 == BASE_C) {
    shift = 2;
  }
  else if( base1 == BASE_G && (base2 == BASE_T|| base2 == BASE_U) ) {
    shift = 4;
  }
  else if( (base1 == BASE_T|| base1 == BASE_U) && base2 == BASE_A) {
    shift = 3;
  }
  else if( (base1 == BASE_T|| base1 == BASE_U) && base2 == BASE_G) {
    shift = 5;
  }
  else {
//    printf("Error in GetMismatchShift. %c and %c don't pair!\n", 
//           Int2base (base1), Int2base (base2));
    shift=-1; // No longer exit...
  }
  return shift;
}

/* ************************************************** */
DEV
int GetPairType( int b) { //assume pair of b is the watson crick pair
  int shift;
  
  if( b == BASE_A) { 
    shift = 0;
  }
  else if( b == BASE_C) {
    shift = 1;
  }
  else if( b == BASE_G) {
    shift = 2;
  }
  else if( b == BASE_T|| b == BASE_U) {
    shift = 3;
  }
  else {
    printf("Error in GetPairType: %d!\n",b);
    return -1;
  }
  return shift;
}


/* ******************************************** */
int fbixIndexOld( int d, int i, int size, int N ) {
  if( d < 0 || i < 0 || i + d >= N) {
    fprintf(stderr,  "Error in fbixIndex %d %d %d %d\n", d, i, size, N);
  }
  if( i*(d-1) + size > (N-d)*(d-1) ) {
    fprintf(stderr,   "Error in fbixIndex! %d > %d, %d %d %d %d,\n", i*(d-1) + size,
           (N-1-d)*(d-1), d,i,size,N);
  }

  return i*(d-1) + size;
}


/* ****** */
int fileExists( char* file) {
  FILE *fp;
  int exists = TRUE;

  fp = fopen( file, "r");
  if( fp == NULL) exists = FALSE;
  else fclose(fp);

  return exists;
}

/* ******** */

