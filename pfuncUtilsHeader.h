/*
  pfuncUtilsHeader.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 7/2006, Justin Bois 1/2007

  Note that DBL_TYPE is defined in pfuncUtilsConstants.h, and
  indicates what type of floating point variables should be used
  (float, double, long double)

  See below for descriptions of each function
*/

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#define DEV __device__
#define GLB __global__

#include "pfuncUtilsConstants.h"
#include "runtime_constants.h"
#include "physical_constants.h"
#include "utilsHeader.h"
#include "DNAExternals.h"

/* ************************ */
// Macros

// C compiler needs this! DO NOT REMOVE!  JZ
#ifndef MIN
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#endif
//This macro returns the min of x,y, regardless of type
#ifndef MAX
#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#endif

/* ****************************** */

//fold struct describes a secondary structure and sequence
typedef struct{
  int *seq; //the sequence
  int seqlength; //sequence length
  int *pairs; //array indicating what is paired with what (see nsStarPairsOrParens)
  int *pknots; //similar to pairs, but indicates ends of pseudoknots
  int *fixedBases; //this is used in design code only, and restricts the identity of a postion
  int  *isNicked; //indicates if a position is right before the end of a strand
  int nStrands; //number of strands in a multi-stranded complex
} fold;

//paramter sets.  DNA, RNA are mfold 2.3 parameters that allow for temperature dependence
//RNA37 is the mfold3.0 parameter set which is only good at 37C.  COUNT sets energies to zero,
//so that the "partition function" is simply a count of the number of possible structures
enum { DNA, RNA, RNA37, USE_SPECIFIED_PARAMETERS_FILE, COUNT};
enum { FALSE, TRUE};


// From utils.c in shared directory:
  /*
     Calculates the number of moles of water per liter at temperature
     (T in degrees C).
  */
double WaterDensity(double T);


/* ******************************************************************************** */
/*                          BEGIN FUNCTIONS FROM PF.C                               */
/* ******************************************************************************** */
/*
   pfuncFull: Calculates the partition function.
   Arguments:
   InputSeq is the sequence, including '+' characters to demarcate strand breaks.

   Complexity = 3 is the multi-stranded, pseudoknot-free method  [ O(N^3)]
   Complexity = 5 is the single-stranded, pseudoknot algorithm   [ O(N^5)]

   naType is chosen from the enum list given above, and chooses the parameter set

   dangles = 0, no dangle energies
   dangles = 1, mfold treatment of dangles
   dangles = 2, dangle energies are always summed, regardless of nearby structures
   (same as the Vienna package -D2 option)
   temperature = the temperature in celsius
   calcPairs = 1 indicates that pair probabilities are calculated and stored in global pairPr
   calcPairs = 0 skips pair probability calculations

   Ignores the possibility of symmetry
*/

DBL_TYPE pfuncFull( int inputSeq[], int complexity, int naType, int dangles,
                    DBL_TYPE temperature, int calcPairs,
                    DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);


GLB
void pfuncFullWithSymHelper(DBL_TYPE *pf, int ** inputSeqs, int * seqlengths,
    int * nStrands, int * permSymmetries);


/* pfunc
   Calls pfuncFull, and assumes complexity = 3, DNA parameters, T = 37, dangles = 1,
   calcPairs = 1, [Na+] = 1.0, [Mg++] = 0.0, and short helix model for salt correction
*/
DBL_TYPE pfunc( int seq[]);


/* ******************************************************************************** */
/*                            END FUNCTIONS FROM PF.C                               */
/* ******************************************************************************** */



/* ******************************************************************************** */
/*                         BEGIN FUNCTIONS FROM INIT.C                              */
/* ******************************************************************************** */
/* getSequenceLength() counts the number of strands (stored as nStrands)  in a given sequence (seq)
   and returns the sequence length */
int getSequenceLength( char *seq, int *nStrands);

/* getSequenceLength() counts the number of strands (stored as nStrands)  in a given sequence (seq)
   and returns the sequence length */
int getSequenceLengthInt( int seq[], int *nStrands);

/* processMultiSeqence() copies input sequence (containing strandPLUS members) to seq, but without the
   strandPLUS members.  The location of the strand breaks, strandPLUS, are stored in the nicks array */
DEV void processMultiSequence( int inputSeq[], int seqlength, int nStrands,
                           int seq[], int nicks[]);

//Allocates Q and sets the values to zero.
DEV void InitLDoublesMatrix( DBL_TYPE **Q, int size, char name[]);

//Sets Q to all zero
void ClearLDoublesMatrix(DBL_TYPE **Q, int size, char name[]);

//Memory management for "fast" interior loops subroutine
DEV void manageQx( DBL_TYPE **Qx, DBL_TYPE **Qx_1,
               DBL_TYPE **Qx_2, int len, int seqlength);


// Salt correction
DBL_TYPE computeSaltCorrection(DBL_TYPE sodiumConc, DBL_TYPE magenesiumConc,
                                 int useLongHelix);

//Load energy parameters.  Global variable DNARNACOUNT determines parameter set
void LoadEnergies(energy_model_t *energy_model, DBL_TYPE temp_k);
void setParametersToZero(energy_model_t *energy_model);

//Set Q[ pf_index(i, i-1, seqlength)] = 1;
DEV void nonZeroInit( DBL_TYPE Q[], int seq[], int seqlength, energy_model_t *em);

/*InitEtaN() initializes the etaN array.  etaN[ EtaNIndex(i,j,seqlength)][0] is the
  number of nicks between i and  j (i,j in [0.5, 1.5, 2.5, ...]).
  etaN[ EtaNIndex(i,j,seqlength)][1] is the index of the leftmost nick
  in the region between i and j, i.e. nicks[ EtaNIndex...] is the position
  of the leftmost nick between i and j.
*/
DEV void InitEtaN( int **etaN, const int *nicks, int seqlength);
int EtaNIndex_old( float i, float j, int seqlength);

/* ******************************************************************************** */
/*                           END FUNCTIONS FROM INIT.C                              */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                     BEGIN FUNCTIONS FROM PFUNCUTILS.C                            */
/* ******************************************************************************** */
//pf_index calculates the array index for a Q-type array
int pf_index_old( int i, int j, int N);
#define pf_index(i,j,N) ((j)==(i)-1?(N)*((N)+1)/2 + (i) : ((i)*(N)+(j)-(i)*(1+(i))/2))
#define pf_index_same(i,N) ((i)*(N)-(i)*((i)-1)/2)
#define EtaNIndex(i,j,N) pf_index((int)(i),(int)(j),N)
#define EtaNIndex_same(i,N) pf_index_same((int)(i),N)
//((int)(j)==(int)(i)-1?(int)((N)*((N)+1)/2 + (int)(i)) : ((int)(i)*(N)+(j)-(i)*(1+(i))/2))
//gap_index calculates the array index of a "gap" matrix.
int gap_index( int h, int r, int m, int s, int seqlength);

//fbixIndex computes the array index for a Qx/Fx array (fast i loops)
int fbixIndexOld( int d, int i, int size, int N );
#ifndef DEBUG
#define fbixIndex(d, i, size, N ) ((i)*((d)-1) + (size))
#else
#define fbixIndex(d, i, size, N ) fbixIndexOld(d, i, size, N )
#endif

//converts A->1, C->2, G->3, T/U->4
int Base2int( char base);
char Int2base( int base);
// Converts sequence of letters to sequence of
// returns -1 on success, index of faulty base letter on error
int convertSeq(char *seqchar, int* seqnum, int seqlength);
int printSeqNum(int seqnum[]);

//converts a pair to an index (For energy calculations)
DEV int GetMismatchShift( int base1, int base2);

//Returns the pair type (similar to GetMismatchShift), except assumes
//a Watson Crick (non-Wobble) pair
DEV int GetPairType( int b);

//returns TRUE if two bases can form a pair.
#ifdef NOGU
#define CanPair(i,j) ((i)+(j)==5 ?TRUE : FALSE)
#else
#define CanPair(i,j) ((i)+(j)==5 || (i)+(j)==7?TRUE : FALSE)
#endif

//check if a file exists
int fileExists( char*);
/* ******************************************************************************** */
/*                       END FUNCTIONS FROM PFUNCUTILS.C                            */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                       BEGIN FUNCTIONS FROM ENERGY.C                              */
/* ******************************************************************************** */
//Nearest neighbor energy of two consecutive base pairs
DEV DBL_TYPE HelixEnergy( int i, int j, int h, int m, energy_model_t *em);

//interior mismatch energy
DEV DBL_TYPE InteriorMM( char a, char b, char x, char y, energy_model_t *em);

//hairpin energy
DEV DBL_TYPE HairpinEnergy( int i, int j, int seq[], energy_model_t *em);

//interior loop energy
DEV DBL_TYPE InteriorEnergy(  int i, int j, int h, int m, int seq[], energy_model_t *em);

//interior loop energy, but allows for the exclusion of the i-j mismatch
//contribution (for use with fast i loop calculations)
DEV DBL_TYPE InteriorEnergyFull( int i, int j, int h, int m, int seq[], int, energy_model_t *em);

//Calculates the dangle energy of a subsequence (i,j), assuming
//i-1, j+1 are paired (unless near a strand break).
//DangleEnergy miscalculates the dangles of nearby wobble pairs
DEV DBL_TYPE DangleEnergy( int i, int j, int seq[], int seqlength, energy_model_t *em);

//Calculates exp(-(dangle energy)/RT) and
//exp( -(interior loop energy)/RT), respectively
DEV DBL_TYPE ExplDangle( int i, int j, int seq[], int seqlength, energy_model_t *em);
DEV DBL_TYPE ExplInternal( int i, int j, int h, int m, int seq[], energy_model_t *em);

//NickDangle calculates the dangle energy, taking into account the effects
//of strand breaks (nicks).  If hairpin == TRUE, then this region is a nicked hairpin
//and may be closed by a wobble pair
DEV DBL_TYPE NickDangle(  int i, int j, const int *nicks, int **etaN, int hairpin,
                      int seq[], int seqlength, energy_model_t *em);

/* Computes the energy of an exterior loop with no secondary structure,
   and returns either exp( -energy/RT) or simply energy*/
DEV DBL_TYPE NickedEmptyQ( int i, int j, const int nicks[], int seq[],
                       int seqlength, int **etaN, energy_model_t *em);

// Lookup table for the function 1.75*kB*TEMP_K*LOG_FUNC( size/30.0)
DEV DBL_TYPE sizeLog(int size, energy_model_t *em);

//Computes the contribution of asymmetry and size to a large interior loop
DEV DBL_TYPE asymmetryEfn( int L1, int L2, int size, energy_model_t *em);
/* ******************************************************************************** */
/*                         END FUNCTIONS FROM ENERGY.C                              */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                         BEGIN FUNCTIONS FROM SUMEXP.C                            */
/* ******************************************************************************** */
//Hairpin energy (exp)
DEV DBL_TYPE ExplHairpin( int i, int j, int seq[], int seqlength, int **etaN, energy_model_t *em);

//Calculates the contribution to the partition function of multiloops (non-nicked)
DEV DBL_TYPE SumExpMultiloops( int i, int j, int seq[],
                           DBL_TYPE *Qms, DBL_TYPE *Qm, int seqlength,
                           int **etaN, energy_model_t *em);
//Calculates the contribution of exterior loops
DEV DBL_TYPE SumExpExteriorLoop( int i,int j, int seq[], int seqlength,
                             DBL_TYPE *Q,
                             int *nicks, int **etaN, energy_model_t *em);

//Computes Qs, Qms  (pairs in exterior loops, multi loops)
DEV void MakeQs_Qms( int i, int j, int seq[], int seqlength,
                 DBL_TYPE *Qs, DBL_TYPE *Qms, DBL_TYPE *Qb,
                 int *nicks, int **etaN, energy_model_t *em);

//Computes Q, Qm for complexity = 3 algorithm
DEV void MakeQ_Qm_N3( int i, int j, int seq[], int seqlength,
                  DBL_TYPE *Q, DBL_TYPE *Qs,
                  DBL_TYPE *Qms, DBL_TYPE *Qm,
                  int *nicks, int **etaN, energy_model_t *em);

//Efficiently calculates the contribution of large interior loops
DEV void fastILoops( int i, int j, int L, int seqlength, int seq[],
                 int **etaN,
                 DBL_TYPE *Qb, DBL_TYPE *Qx, DBL_TYPE *Qx_2, energy_model_t *em);


//makeNewQx creates new "extensible" base cases for the interval i,j.
DEV void makeNewQx( int i, int j, int seq[], int seqlength,
                int **etaN, DBL_TYPE Qb[], DBL_TYPE Qx[], energy_model_t *em);
//extendOldQx extends Qx for the i-1, j+1 case
DEV void extendOldQx( int i, int j, int seqlength,
                  DBL_TYPE Qx[], DBL_TYPE Qx_2[], energy_model_t *em);

//Directly calculates the contribution of small interior loops
DEV DBL_TYPE SumExpInextensibleIL( int i, int j, int seq[], int seqlength,
                               DBL_TYPE Qb[],  int **etaN, energy_model_t *em);
/* ******************************************************************************** */
/*                           END FUNCTIONS FROM SUMEXP.C                            */
/* ******************************************************************************** */



#endif

