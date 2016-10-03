/*
 Copyright (C) 2013 X. Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

 $Id$
*/
/*ET20160914: 
The original file found at http://www.tddft.org/svn/octopus/trunk/src/grid/hilbert.c
with some small modifications to compile under c++
*/

//ET20160914: one line removed
//#include <config.h>

#include <assert.h>

/* This code is based on the implementation of the Hilbert curves
   presented in:

   J. Skilling, Programming the Hilbert curve, AIP Conf. Proc. 707, 381 (2004); http://dx.doi.org/10.1063/1.1751381

*/

//ET20160914: put functions in namespace Hilbert
namespace hilbert_c {

void InttoTranspose(const int dim, const long long int h, int * x)
{/* the code uses some funny way of storing the bits */
	int idir, ibit, ifbit;

	for(idir = 0; idir < dim; idir++)
		x[idir] = 0;

	ifbit = 0;
	for(ibit = 0; ibit < 21; ibit++) {
		for(idir = dim - 1; idir >= 0; idir--) {
			x[idir] += (((h>>ifbit)&1)<<ibit);
			ifbit++;
		}
	}
}

void TransposetoAxes
  ( int* X // position
  , int b  // #bits
  , int n  // #dimensions
  )
{
	  int N = 2 << (b-1), P, Q, t;
	  int i;

 // Gray decode by H ^ (H/2)
	t = X[n-1] >> 1;
	for(i = n - 1; i > 0; i--) {
		X[i] ^= X[i-1];
	}
	X[0] ^= t;

 // Undo excess work
	for( Q = 2; Q != N; Q <<= 1 ) {
		P = Q - 1;
		for( i = n-1; i >= 0 ; i-- ) {
			if( X[i] & Q ) {// invert
				X[0] ^= P;
			} else {// exchange
				t = (X[0]^X[i]) & P;
				X[0] ^= t;
				X[i] ^= t;
			}
		}
    }
}

void AxestoTranspose
  ( int    * X	// I O  position   [n]
  , int      b	// I    # bits
  , int      n	// I    dimension
  )
{
	int P, Q, t; // originally coord_t instead of int
	int i;

 // Inverse undo
	for( Q = 1 << (b - 1); Q > 1; Q >>= 1 ) {
		P = Q - 1;
		if( X[0] & Q )	   			// invert
			X[0] ^= P;
		for( i = 1; i < n; i++ ) {
			if( X[i] & Q )
				X[0] ^= P;			// invert
			else { 					// exchange
				t = (X[0] ^ X[i]) & P;
				X[0] ^= t;
				X[i] ^= t;
			}
		}
	}

 // Gray encode (inverse of decode)
	for( i = 1; i < n; i++ ) {
		X[i] ^= X[i-1];
	}
	t = X[n-1];
	for( i = 1; i < b; i <<= 1 ) {
		X[n-1] ^= X[n-1] >> i;
	}
	t ^= X[n-1];
	for( i = n-2; i >= 0; i-- ) {
		X[i] ^= t;
	}
}

}//namespace hilbert_c
