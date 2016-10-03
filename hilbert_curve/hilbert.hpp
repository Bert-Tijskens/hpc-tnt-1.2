/*
 * hilbert.hpp provides a more practical interface to the functions in 
 *   http://www.tddft.org/svn/octopus/trunk/src/grid/hilbert.c 
 * queued to molecular dyynamics simulations.
 *
 *  Created on: 14 Sep 2016
 *      Author: etijskens
 */

#ifndef HILBERT_HPP_
#define HILBERT_HPP_

#include <cstddef> // size_t
#include <cmath> // floor
//#include <iostream>
#include <string>

#define DIMS 3
#define BITS 7
 // i,j,k can vary between 0 and 2**BITS - 1
 // h can vary between 0 and (2**BITS)**3 - 1
#ifdef VALIDATE
 // Check validity of cell and hilbert indices on input
    #define VALIDATE_IJK(i) validate_cell_index(i);
    #define VALIDATE_H(h) validate_hilbert_index(h);
#else
    #define VALIDATE_IJK(i)
    #define VALIDATE_H(h)
#endif

namespace hilbert
{//-------------------------------------------------------------------------------------------------
    typedef long long int HilbertIndex_t;
    typedef  unsigned int I_t;
     // couldn't get this to work with size_t, because I could not find out which numpy type
     // corresponds to size_t.
    static const int MAXIJK = (1 << BITS);
    static const int MAXH   = (1 << (BITS*3));
    void validate_cell_index(int const i);
    inline bool is_valid_ijk(int const* ijk) {
        return ijk[0]<MAXIJK &&    ijk[1]<MAXIJK &&    ijk[2]<MAXIJK
         && 0<=ijk[0]        && 0<=ijk[1]        && 0<=ijk[2];
    }
    inline bool is_valid_i(int const i) { return 0<=i && i<MAXIJK; }
    inline bool is_valid_h(int const h) { return 0<=h && h<MAXH; }

    void validate_hilbert_index(HilbertIndex_t const i);
 //-------------------------------------------------------------------------------------------------
 // return string with current constants and valid ranges of cell and Hilbert indices.
    std::string info();
 // return the largest valid cell index
    inline HilbertIndex_t cell_index_limit() { return MAXIJK; }
 // return the largest Hilbert index
    inline HilbertIndex_t hilbert_index_limit() { return MAXH; }
 // test if the library was compiled with VALIDATE defined
    inline bool is_validating() {
#ifdef VALIDATE
        return true;
#else
        return false;
#endif
    }
 //-------------------------------------------------------------------------------------------------
 // Hilbert index to 3D cell index
    void h2ijk( const long long int h
              , int* ijk // array
                         // Beware: ijk is modified!!! DO NOT use if you want to keep the values of i,j,k
              );
 //-------------------------------------------------------------------------------------------------
 // 3D cell index to Hilbert index
 	HilbertIndex_t ijk2h( int* ijk ); // ijk is array

 	inline HilbertIndex_t
 	ijk2h( int const i, int const j, int const k )
    {
        VALIDATE_IJK(i)
        VALIDATE_IJK(j)
        VALIDATE_IJK(k)
        int ijk[3] = {i,j,k};
        return ijk2h(ijk);
    }
 //-------------------------------------------------------------------------------------------------
 // 3D coordinates to Hilbert index
	inline HilbertIndex_t
	xyzw2h( double const x, double const y, double const z, double const w)
    {
        return ijk2h( static_cast<int>(floor(x/w))
                    , static_cast<int>(floor(y/w))
                    , static_cast<int>(floor(z/w)) );
    }
 //-------------------------------------------------------------------------------------------------
 // 3D coordinates to cell index and Hilbert index
    inline void
    xyzw2ijkh_float64
      ( double const    x   // input
      , double const    y   // input
      , double const    z   // input
      , double const    w   // input
      , int           & i   // output
      , int           & j   // output
      , int           & k   // output
      , HilbertIndex_t& h   // output
      )
    {
        i = static_cast<int>( floor(x/w) );
        j = static_cast<int>( floor(y/w) );
        k = static_cast<int>( floor(z/w) );
        h = ijk2h(i,j,k);
    }
    inline void
    xyzw2ijkh_float32
      ( float const     x   // input
      , float const     y   // input
      , float const     z   // input
      , float const     w   // input
      , int           & i   // output
      , int           & j   // output
      , int           & k   // output
      , HilbertIndex_t& h   // output
      )
    {
        i = static_cast<int>( floor(x/w) );
        j = static_cast<int>( floor(y/w) );
        k = static_cast<int>( floor(z/w) );
        h = ijk2h(i,j,k);
    }

 //-------------------------------------------------------------------------------------------------
    void insertion_sort
      ( I_t      const   n // (input) array length of h and I
      , HilbertIndex_t * h // (input and output) array of hilbert indices, sorted on output
      , I_t            * I // (output) array of new positions I[i] of i-th element
      );
     // Use insertion sort compute the new position I[i] of each h[i]
     // to sort any other array A in the same order:
     //     for( I_t i=0; i<n; i++ )
     //         Asorted[i] = A[[I[i]];
 //-------------------------------------------------------------------------------------------------
 // create a sorted copy of <unsorted>
    template <class T>
    void reorder( I_t const n, I_t const* I, T const * told, T* tnew )
    {
        for( I_t i=0; i<n; i++ ) {
            tnew[i] = told[I[i]];
        }
    }
 //-------------------------------------------------------------------------------------------------
#ifdef OBSOLETE
    template <class T>
    void reorder_inplace( I_t const n, T* t, I_t const* I )
    {// this algorithm is not general enough! It only works if all elements are visited.
     // it does not work if an element is already at the correct position
     // or when a subset is rotated.
     // cfr http://stackoverflow.com/questions/7365814/in-place-array-reordering
     // There it is stated that this can only be done by "destroying" the index array
     // this obviously makes it useless for our purpose.


//        std::cout<<"n="<<n<<std::endl;
//        for(int j=0;j<n;++j) {
//            std::cout<<I[j]<<' ';
//        }   std::cout<<std::endl;
        I_t i = I[0];
//        std::cout<<"i="<<i<<std::endl;
     // copy the element that eventually must become the first element:
        T t0 = t[i];
        while(i>0) {
//            for(int j=0;j<n;++j) {
//                std::cout<<t[j]<<' ';
//            }   std::cout<<std::endl;
            I_t Ii = I[i];
            t[i] = t[Ii];
            i = Ii;
        }
        t[0] = t0;
//        for(int j=0;j<n;++j) {
//            std::cout<<t[j]<<' ';
//        }   std::cout<<std::endl;
    }
#endif
 //-------------------------------------------------------------------------------------------------
}// namespace hilbert

#endif // HILBERT_HPP_ 
