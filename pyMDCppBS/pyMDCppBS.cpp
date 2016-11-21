/*
 * pyMDCppBS.cpp
 *
 *  Created on: 26 Oct 2016
 *      Author: etijskens
 */
#include <boost/python.hpp>
using namespace boost::python;

#include <iostream>
#include "numpy_boost_python.hpp"

#include "md_bs.hpp"

float // cputime
compute_interactions_verlet_list_64_test
  ( numpy_boost<double,1> rx
  , numpy_boost<double,1> ry
  , numpy_boost<double,1> rz
  , numpy_boost<double,1> ax
  , numpy_boost<double,1> ay
  , numpy_boost<double,1> az
  , numpy_boost<int   ,1> verlet_n // number of pairs for each atom
  , numpy_boost<int   ,1> verlet_j // indices of second atom in each pair
  , int n_iter
  , int imax
  )
{
    return
    md::compute_interactions_verlet_list
      ( rx.origin()
      , ry.origin()
      , rz.origin()
      , ax.origin()
      , ay.origin()
      , az.origin()
      , verlet_n.origin()
      , verlet_j.origin()
      , static_cast<size_t>( rx      .shape()[0] ) // n_atoms
      , static_cast<size_t>( verlet_j.shape()[0] ) // n_pairs
      , static_cast<size_t>( n_iter )
      , static_cast<size_t>( imax )
      );
}

    float // returns cputime
    compute_interactions_verlet_list_64
      ( numpy_boost<double,1> rx
      , numpy_boost<double,1> ry
      , numpy_boost<double,1> rz
      , numpy_boost<double,1> ax
      , numpy_boost<double,1> ay
      , numpy_boost<double,1> az
      , numpy_boost<int   ,1> verlet_n // number of pairs for each atom
      , numpy_boost<int   ,1> verlet_j // indices of second atom in each pair
      , int n_iter
      )
    {
        return 
        md::compute_interactions_verlet_list
          ( rx.origin()
          , ry.origin()
          , rz.origin()
          , ax.origin()
          , ay.origin()
          , az.origin()
          , verlet_n.origin()
          , verlet_j.origin()
          , static_cast<size_t>( rx.shape()[0] ) // n_atoms
          , static_cast<size_t>( rx.shape()[0] ) // n_pairs
          , static_cast<size_t>( n_iter )
          );
    }

BOOST_PYTHON_MODULE(pyMDCppBS)
{
 // Initialize the Numpy support
    IMPORT_ARRAY();

 // You must call this function inside of the BOOST_PYTHON_MODULE
 // init function to set up the type conversions.  It must be called
 // for every type and dimensionality of array you intend to return.
 // If you don't, the code will compile, but you will get error
 // messages at runtime.
    numpy_boost_python_register_type<int   ,1>();
    numpy_boost_python_register_type<double,1>();
    numpy_boost_python_register_type<float ,1>();
    
 // Exported functions:
    def("compute_interactions_verlet_list_test",compute_interactions_verlet_list_64_test);
    def("compute_interactions_verlet_list"     ,compute_interactions_verlet_list_64);
}
