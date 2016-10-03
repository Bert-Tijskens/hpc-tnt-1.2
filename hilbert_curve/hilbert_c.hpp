namespace hilbert_c
{//=============================================================================
 // this namespace wraps the original functions in
 // http://www.tddft.org/svn/octopus/trunk/src/grid/hilbert.c
 //=============================================================================

	void InttoTranspose(const int dim, const long long int h, int * x);
	
	void TransposetoAxes
	  ( int* X // position
	  , int b  // #bits
	  , int n  // #dimensions
	  );
	  
	void AxestoTranspose
	  ( int    * X	// I O  position   [n]
	  , int      b	// I    # bits
	  , int      n	// I    dimension
	  );

 //=============================================================================
}// namespace hilbert_c
