#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include <Rinternals.h>
#include "R_ext/Applic.h"

/*
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation (a copy of the GNU
 *  General Public License is available at
 *  http://www.r-project.org/Licenses/
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 */

/*
 * Grid.Binning is designed to redistributed the weight along a
 * defined grid.
 *
 *  ngrid is the number of grid points, which equals to the bin number
 *  plus one.
 *
 *  Last updated: Nov 2, 2012
 */
void GridBinning(double *x, double *w, int *nx,
		 double *xlo, double *bw, int *ngrid,
		 int *truncate, int *linbin, double *gcounts)
{
  int i, li, m=ngrid[0], n=nx[0];
  double binwidth = bw[0], lxi, rem, a=xlo[0];
  
  for(i=0; i<m; i++) gcounts[i] = 0.0;
  
  for(i=0; i<n; i++){
    lxi = (x[i] - a)/binwidth;
    li = (int) lxi;
    if(linbin[0] == 1)
      rem = lxi - li;
    else
      rem = 0.0;

    if((li>0) && (li<m-1)){
      gcounts[li] += (1.0-rem) * w[i];
      gcounts[li+1] += rem * w[i];
    }
    
    if((li<=0)&&(truncate==0))
      gcounts[0] += w[i];
    if((li>=m-1)&&(truncate==0)&&(linbin[0] == 1))
      gcounts[m-1] += w[i];
    if((li>=m-1)&&(truncate==0)&&(linbin[0] == 0))
      gcounts[m-2] += w[i];
  }
}

