/*
 * hilbert_curve_py.cpp
 *
 *  Created on: 15 Sep 2016
 *      Author: etijskens
 */

#include <iostream>
#include "numpy_boost_python.hpp"
#include "hilbert.hpp"

// convert 3D positions (numpy arrays) to hilbert indices 
void xyzw2h_float64
  ( numpy_boost<double,1> x
  , numpy_boost<double,1> y
  , numpy_boost<double,1> z
  , double w
  , numpy_boost<hilbert::HilbertIndex_t,1> h 
  )
{
    hilbert::I_t const n = x.shape()[0];
    for( hilbert::I_t i=0; i<n; i++ ) {
        h[i] = hilbert::xyzw2h(x[i],y[i],z[i],w);
    }   
}

void xyzw2h_float32
  ( numpy_boost<float,1> x
  , numpy_boost<float,1> y
  , numpy_boost<float,1> z
  ,             float    w
  , numpy_boost<hilbert::HilbertIndex_t,1> h
  )
{
    hilbert::I_t const n = x.shape()[0];
    for( hilbert::I_t i=0; i<n; i++ ) {
        h[i] = hilbert::xyzw2h(x[i],y[i],z[i],w);
    }
}

void xyzw2ijkh_float64
  ( numpy_boost<double,1> x // input
  , numpy_boost<double,1> y // input
  , numpy_boost<double,1> z // input
  , double                w // input
  , numpy_boost<int                    ,1> i // output
  , numpy_boost<int                    ,1> j // output
  , numpy_boost<int                    ,1> k // output
  , numpy_boost<hilbert::HilbertIndex_t,1> h // output
  )
{
    hilbert::I_t const n = x.shape()[0];

    for( hilbert::I_t ip=0; ip<n; ip++ ) {
        hilbert::xyzw2ijkh_float64( x[ip], y[ip], z[ip], w
                                  , i[ip], j[ip], k[ip], h[ip] );
    }
}

void xyzw2ijkh_float32
  ( numpy_boost<float,1> x // input
  , numpy_boost<float,1> y // input
  , numpy_boost<float,1> z // input
  ,             float    w // input
  , numpy_boost<int                    ,1> i // output
  , numpy_boost<int                    ,1> j // output
  , numpy_boost<int                    ,1> k // output
  , numpy_boost<hilbert::HilbertIndex_t,1> h // output
  )
{
    hilbert::I_t const n = x.shape()[0];
    for( hilbert::I_t ip=0; ip<n; ip++ ) {
        hilbert::xyzw2ijkh_float32( x[ip], y[ip], z[ip], w
                                  , i[ip], j[ip], k[ip], h[ip] );
    }
}

void sort
  ( numpy_boost<hilbert::HilbertIndex_t,1> h
  , numpy_boost<hilbert::I_t           ,1> I  
  )
{
    hilbert::I_t const n = h.size();
    hilbert::insertion_sort( n, h.origin(), I.origin() );
}

void reorder_float64
  ( numpy_boost<hilbert::I_t,1> I  
  , numpy_boost<      double,1> told
  , numpy_boost<      double,1> tnew
  )
{
    hilbert::reorder( I.size(), I.origin(), told.origin(), tnew.origin() );
}

void reorder_float32
  ( numpy_boost<hilbert::I_t,1> I  
  , numpy_boost<       float,1> told
  , numpy_boost<       float,1> tnew
  )
{
    hilbert::reorder( I.size(), I.origin(), told.origin(), tnew.origin() );
}

void reorder_int32
  ( numpy_boost<hilbert::I_t,1> I  
  , numpy_boost<         int,1> told
  , numpy_boost<         int,1> tnew
  )
{
    hilbert::reorder( I.size(), I.origin(), told.origin(), tnew.origin() );
}

void reorder_uint32
  ( numpy_boost<hilbert::I_t,1> I  
  , numpy_boost<unsigned int,1> told
  , numpy_boost<unsigned int,1> tnew
  )
{
    hilbert::reorder( I.size(), I.origin(), told.origin(), tnew.origin() );
}

inline hilbert::HilbertIndex_t
ijk2h_1( numpy_boost<int,1> ijk )
{
    int const *ijk_data = ijk.origin();
    if( hilbert::is_valid_ijk(ijk_data) ) {
        hilbert::HilbertIndex_t h = hilbert::ijk2h(ijk_data[0],ijk_data[1],ijk_data[2]);
        return h;
    } else{
        return -1;
    }
}

inline bool
h2ijk_1( hilbert::HilbertIndex_t const h, numpy_boost<int,1> ijk )
{
    if( hilbert::is_valid_h(h) ) {
        hilbert::h2ijk(h,ijk.origin());
        return true;
    } else{
        return false;
    }
}

void
h2ijk( hilbert::HilbertIndex_t const h,numpy_boost<int,1> ijk )
{
    return hilbert::h2ijk(h,ijk.origin());
}

using namespace boost::python;

BOOST_PYTHON_MODULE(pyHilbertCpp)
{
 // Initialize the Numpy support
    IMPORT_ARRAY();

 // You must call this function inside of the BOOST_PYTHON_MODULE
 // init function to set up the type conversions.  It must be called
 // for every type and dimensionality of array you intend to return.
 // If you don't, the code will compile, but you will get error
 // messages at runtime.
    numpy_boost_python_register_type<hilbert::HilbertIndex_t,1>();
    numpy_boost_python_register_type<hilbert::I_t           ,1>();
//  numpy_boost_python_register_type<         unsigned int  ,1>(); //same as hilbert::I_t
    numpy_boost_python_register_type<         int           ,1>();
    numpy_boost_python_register_type<         double        ,1>();
    numpy_boost_python_register_type<         float         ,1>();
    
 // Exported functions:
    def("xyzw2h_float64", xyzw2h_float64 ); // 3D positions to hilbert index of the corresponding cell with width w.
    def("xyzw2h_float32", xyzw2h_float32 ); // 3D positions to hilbert index of the corresponding cell with width w.

    def("xyzw2ijkh_float32", xyzw2ijkh_float32 ); // 3D positions to cell indices and hilbert index of the corresponding cell with width w.
    def("xyzw2ijkh_float64", xyzw2ijkh_float64 ); // 3D positions to cell indices and hilbert index of the corresponding cell with width w.

    def("sort"  , sort   ); // sort the hilbert indices and produce reordering array.
    def("reorder_float32", reorder_float32);
    def("reorder_float64", reorder_float64);
    def("reorder_int32"  , reorder_int32  );
    def("reorder_uint32" , reorder_uint32 );

    def("ijk2h_1",ijk2h_1);
    def("h2ijk_1",h2ijk_1);
    def("info",hilbert::info);
    def("cell_index_limit",hilbert::cell_index_limit);
    def("hilbert_index_limit",hilbert::hilbert_index_limit);

    def("validate_cell_index",hilbert::validate_cell_index);
    def("validate_hilbert_index",hilbert::validate_hilbert_index);
    def("is_validating",hilbert::is_validating);
    def("is_valid_ijk",hilbert::is_valid_ijk);
    def("is_valid_i",hilbert::is_valid_i);
    def("is_valid_h",hilbert::is_valid_h);
}




