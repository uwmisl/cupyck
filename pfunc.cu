/*  
    pfunc.c is part of the NUPACK software suite
    Copyright (c) 2007 Caltech. All rights reserved.
    Coded by: Robert Dirks, 5/2005 and Justin Bois 1/2007   
    
    The purpose of this program is to calculate the partition function
    of all possible secondary structures of a given strand (or
    strands) of DNA/RNA.
    
    If there are multiple strands, then this algorithm will calculate
    the partition function for a single circular permutation of those
    strands, assuming each strand to be distinguishable from the
    others.  This algorithm is described in our SIAM Review paper
    published in 2007 (Dirks, Bois, Schaeffer, Winfree, Pierce).  The
    time complexity of this algorithm is O(N^3), where N is the total
    sequence of all the strands involved.
    
    For a single strand, the algorithm can be expanded to allow for
    the simplest kinds of pseudoknots, but with computational
    complexity of O(N^5), and storage complexity of O(N^4).  This
    algorithm is described in (Dirks, Pierce JCC 24:1664-77, 2003),
    and (Dirks, Pierce, JCC 25:1295-1304, 2004).  The Qp portion of
    the algorithm has been recently expanded from the paper to allow
    for gap spanning regions containing a single base pair.
    
    pfunc.c will compile as a stand alone executable to calculate
    partition functions.  The algorithms can also be compiled as a
    static library to be used by other code.
    
    Default pfunc.c usage: The default setings include the mfold2.3
    RNA parameter set at 37C.  Pseudoknots are enabled with the
    -pseudo flag unless the input sequence has multiple strands.
    These default parameters can be changed by modifying the call to
    pfuncFull.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <signal.h>


#include "pfuncUtilsHeader.h"
#include "DNAExternals.h"

int readMultiInput(char *inputFile, char ***seqs) {
  FILE *fin = fopen(inputFile, "r");
  if (!fin) {
    fprintf(stderr, "unable to open file\n");
    exit(1);
  }
  int nlines = 0;
  int rc = 0;
  while (rc != EOF) {
    rc = fgetc(fin);
    if (rc == '\n') {
      ++nlines;
    }
  }

  *seqs = (char**)malloc(nlines * sizeof(char*));
  rewind(fin);

  size_t len = MAXSEQLENGTH;
  for(int i = 0; i < nlines; ++i) {
    (*seqs)[i] = (char*)malloc(MAXSEQLENGTH);
    int nc = getline(&((*seqs)[i]), &len, fin);
    if (len > MAXSEQLENGTH) {
      fprintf(stderr, "seqs longer than MAXSEQLENGTH");
      exit(1);
    }
    (*seqs)[i][nc-1] = '\0';
  }
  fclose(fin);
  return nlines;
}

/* ************************************************ */

int main( int argc, char *argv[] ) {
  
  char seq[ MAXSEQLENGTH];
  int seqNum[ MAXSEQLENGTH+1];
  
  DBL_TYPE pf;

  int vs;
  int tmpLength;
  char inputFile[ MAXLINE];
  int inputFileSpecified;
  
  strcpy( inputFile, "");
  
  inputFileSpecified = ReadCommandLineNPK( argc, argv, inputFile);
  
  energy_model_t energy_model;
  energy_model.temp_k = TEMP_K;
  energy_model.dangletype = DANGLETYPE;
  energy_model.dnarnacount = DNARNACOUNT;
  LoadEnergies(&energy_model, TEMP_K);

  cudaMemcpyToSymbol(ENERGIES, &energy_model, sizeof(energy_model_t));

  if(NupackShowHelp) {
    printf("Usage: pfunc [OPTIONS] PREFIX\n");
    printf("Calculate the partition function of the input sequence.\n");
    printf("Example: pfunc -multi -T 25 -material dna example\n");
    PrintNupackThermoHelp();
    PrintNupackUtilitiesHelp();
    exit(1);
  }

  header( argc, argv, "pfunc","screen");
  if( !inputFileSpecified || 
      !ReadInputFile( inputFile, seq, &vs, NULL, NULL, NULL) ) {
       if (inputFileSpecified == 0) getUserInput( seq, &vs, NULL, NULL);
       else abort();
  }
  
  printInputs( argc, argv, seq, vs, NULL, NULL,"screen");
  
  //calculate partition function, without pairs info
  tmpLength = strlen( seq);
  convertSeq(seq, seqNum, tmpLength);

  pf = pfuncFullWithSym(seqNum, vs);

  printf("%s\n%s Free energy (kcal/mol) and partition function:\n",
	 COMMENT_STRING,COMMENT_STRING);

  printf("%.8Le\n",-1*(kB*TEMP_K)*logl( (long double) pf));
  printf( "%12.14Le\n", (long double) pf); 
  
  return 0;
}
/* ****** */













