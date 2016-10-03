/*
 * hilbert.cpp
 *
 *  Created on: 14 Sep 2016
 *      Author: etijskens
 */
 
#include "hilbert.hpp"
#include "hilbert_c.hpp"
#include <stdexcept>

#define DEBUG_IO
#ifdef DEBUG_IO
#  include <iostream>
#endif

namespace hilbert
{//-------------------------------------------------------------------------------------------------
    std::string info()
    {
        std::string sMAXIJK = std::to_string(MAXIJK);
        std::string s = std::string("libhilbert.a info")
                            .append("\n  DIMS = ").append( std::to_string(DIMS) ).append(" (number of spatial dimensions)")
                            .append("\n  BITS = ").append( std::to_string(BITS) ).append(" (").append(sMAXIJK).append("x").append(sMAXIJK).append("x").append(sMAXIJK).append(" cells)")
                            .append("\n  Valid cell indices i,j,k must satisfy: 0 <= i,j,k < ").append(sMAXIJK)
                            .append("\n  Valid hilbert indices h  must satisfy: 0 <=     h < ").append( std::to_string(MAXH) )
                            .append("\n  Compiled with cell and Hilbert index validation: ").append((is_validating()?"yes":"no"))
                            ;
        return s;
    }
 //-------------------------------------------------------------------------------------------------
    void 
    h2ijk( const HilbertIndex_t h, int* ijk )
    {
        VALIDATE_H(h)
        hilbert_c::InttoTranspose( DIMS, h, ijk );
        hilbert_c::TransposetoAxes(ijk,BITS,DIMS);
//        std::cout<<"hilbert::h2ijk "<<h<<"->"<<ijk[0]<<' '<<ijk[1]<<' '<<ijk[2]<<"->"<<ijk2h(ijk[0],ijk[1],ijk[2])<<std::endl;
    }
//-------------------------------------------------------------------------------------------------
    struct InvalidIndex : public std::logic_error {
        InvalidIndex( std::string const & what ) : std::logic_error(what) {}
    };
 //-------------------------------------------------------------------------------------------------
    void validate_cell_index(int const i) {
        if( i>=MAXIJK ) {
            std::string what = std::string("InvalidIndex\nCell index ")
                                   .append( std::to_string(i) ).append(" greater than or equal to limit (")
                                   .append( std::to_string(MAXIJK) ).append(").\n")
                                   .append("If this value should be valid, you must increase the value of BITS in hilbert.hpp and recompile the library.")
                                   ;
            throw InvalidIndex(what);
        }
    }
 //-------------------------------------------------------------------------------------------------
    void validate_hilbert_index(HilbertIndex_t const h) {
        if( h>=MAXH ) {
            std::string what = std::string("InvalidIndex\nHilbert index ")
                                   .append( std::to_string(h) ).append(" greater than or equal to limit (")
                                   .append( std::to_string(MAXH) ).append(").\n")
                                   .append("If this value should be valid, you must increase the value of BITS in hilbert.hpp and recompile the library.")
                                   ;
            throw InvalidIndex(what);
        }
    }
 //-------------------------------------------------------------------------------------------------
 //-------------------------------------------------------------------------------------------------
    HilbertIndex_t
    ijk2h( int* ijk ) // Beware: ijk is modified!!!
    {
        hilbert_c::AxestoTranspose( ijk, BITS, DIMS );
        HilbertIndex_t h = 0;
        for( int i=0;i<BITS;++i) {
            for( int j=0;j<DIMS;++j) {
                h += ((ijk[2-j]>>i & 1) << (DIMS*i+j));
            }
        }
        return h;
    }
 //-------------------------------------------------------------------------------------------------
 // insertion sort based on hilbert curve.
 // http://www.sorting-algorithms.com/nearly-sorted-initial-order
 // https://en.wikipedia.org/wiki/Insertion_sort
    void insertion_sort( I_t const n, HilbertIndex_t* h, I_t* I )
    {
        I[0] = 0;
        for( I_t i=1; i<n; ++i ) 
        {
            I[i] = i; // initialize I on the fly to [0,1,2,3,...]
            HilbertIndex_t htmp = h[i];
            I_t            itmp = i;// = I[i]
            I_t j = i;
            while( h[j-1] > htmp )
            {
                #ifdef DEBUG_IO
                std::cout<<j<<std::endl;
                #endif  
                
                h[j] = h[j-1];
                I[j] = I[j-1];
                if(j==0) break; 
                j -= 1;
            }
            h[j] = htmp;
            I[j] = itmp;
            
            #ifdef DEBUG_IO
            for(int i=0; i<8; ++i) {
                std::cout << h[i] << ' ';
            } std::cout << std::endl;
            #endif
       }
    }
 //-------------------------------------------------------------------------------------------------
}// namespace hilbert
