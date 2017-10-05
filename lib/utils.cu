/*
  utils.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois, 1/2007, except where noted.

  UTILS.C 

  General utility functions for processing arrays (e.g., dot products,
  matrix multiplication, etc.), numbers (e.g., factorials, random
  number generation), strings, etc.  There is one physical function in
  here, which is to compute the density of water at a given
  temperature.

  For use with NUPACK.
*/

#include "utilsHeader.h"

/* ******************************************************************************** */
double WaterDensity(double T) {
  /* 
     Calculates the number of moles of water per liter at temperature
     (T in degrees C).

     Density of water calculated using data from:
     Tanaka M., Girard, G., Davis, R., Peuto A.,
     Bignell, N.   Recommended table for the denisty
     of water..., Metrologia, 2001, 38, 301-309
  */
   double a1 = -3.983035;
   double a2 = 301.797;
   double a3 = 522528.9;
   double a4 = 69.34881;
   double a5 = 999.974950;


   return a5 * (1 - (T+a1)*(T+a1)*(T+a2)/a3/(T+a4)) / 18.0152;

}
/* ******************************************************************************** */



