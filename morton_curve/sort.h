//===============================================================================
// Copyright (C) 2011 by Huy T. Vo, Polytechnic Institute Of New York University 
// Contact email: hvo@poly.edu   -  URL: http://vgc.poly.edu/projects/meshlayout 
// This code is for educational purpose only. No guarantees that it's useful for 
// any other purpose. This code is part of the space-filling curve mesh layout.  
//===============================================================================
#ifndef SORT_H
#define SORT_H

#include "octcode.h"

typedef struct {
  OCTCODE code;
  unsigned index;
} OctNode;

#define MIN_FOR_RADIX 32
#define GET_BYTE(a, b) ((a >> (b * 8)) & 0xFF)
#define SWAP(x, y) { OctNode z(x); x = y; y = z; }

inline void insertionSort(OctNode *a, int n)
{
  unsigned i, j;
  OctNode k;
  for (i=0; i!=n; i++, a[j]=k)
	for (j=i, k=a[j]; j && k.code<a[j-1].code; j--) {
	  a[j] = a[j-1];
	}
}

static void inplaceRadixSortByte(OctNode *a, int n, int byte)
{
  if (n<MIN_FOR_RADIX)
	insertionSort(a, n);
  else {
	unsigned i, k, end;
	OctNode j;
	
	unsigned count[256] = {};
	for (i=0; i<n; i++)
	  count[GET_BYTE(a[i].code, byte)]++;
	 
	unsigned bucket[256];
	bucket[0] = 0;
	for (i=1; i<256; i++) {
	  bucket[i] = bucket[i-1] + count[i-1];	  
	}
	
	for (i=0; i<256; i++) {
	  end = (i>0?bucket[i-1]:0) + count[i];
	  for (; bucket[i]<end; bucket[i]++) {
		j = a[bucket[i]];
		while ((k=GET_BYTE(j.code, byte))!=i) {
		  register unsigned xx = bucket[k];
		  bucket[k]++;
		  SWAP(j, a[xx]);
		}
		a[bucket[i]] = j;
	  }
	}

	if (byte-->0) {
	  for (i=0; i<256; i++)
		if (count[i]>0)
		  inplaceRadixSortByte(a + bucket[i]-count[i], count[i], byte);
	}
  }
}

inline void inplaceRadixSort(OctNode *a, int n)
{
  inplaceRadixSortByte(a, n, sizeof(OCTCODE)-1);
}

#endif
