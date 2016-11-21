/*
 * hilbert_curve_py.cpp
 *
 *  Created on: 15 Sep 2016
 *      Author: etijskens
 */
#include <boost/python.hpp>
using namespace boost::python;

#include <iostream>
#include "numpy_boost_python.hpp"
#include <hilbert.hpp>

namespace helper
{
    template <class T>
    inline double sq(T const t ) {
        return t*t;
    }

    int E[13][2][3] =   { {{0,0,0},{1,0,0}}, //0
                          {{0,0,0},{0,1,0}},
                          {{0,0,0},{0,0,1}},
                          {{0,0,0},{1,1,0}},
                          {{0,0,0},{1,0,1}},
                          {{0,0,0},{0,1,1}},
                          {{0,0,0},{1,1,1}},
                          {{1,0,0},{0,1,0}}, //7
                          {{1,0,0},{0,0,1}},
                          {{1,0,0},{0,1,1}},
                          {{0,1,0},{0,0,1}}, //10
                          {{0,1,0},{1,0,1}},
                          {{0,0,1},{1,1,0}}  //12
                        };

    class VerletListBuilder
    {
    public:
     // Ctor
        VerletListBuilder
          ( numpy_boost<int   ,2>& verlet_list
          , numpy_boost<int   ,2>& hl_
          , numpy_boost<double,1> rx
          , numpy_boost<double,1> ry
          , numpy_boost<double,1> rz
          , numpy_boost<int   ,1> I
          , double                rcutoff2
          )
          : n_atoms_(verlet_list.shape()[1])
          , reserve_(verlet_list.shape()[0])
          , rx_(rx.origin())
          , ry_(ry.origin())
          , rz_(rz.origin())
          , I_ (( I.shape()[0]==0 ? nullptr : I.origin() ))
          , vl_(verlet_list)
          , hl_(hl_)
          , rcutoff2_(rcutoff2)
        {
            this->construct();
        }

        void add(int const i, int const j )
        {
            if(j<i) {
                this->vl_[0][i] += 1;
                int const n = this->vl_[0][i];
                this->vl_[n][i] = j;
            } else {
                this->vl_[0][j] += 1;
                int const n = this->vl_[0][j];
                this->vl_[n][j] = i;
            }
        }

        inline int ia2i( int ia ) const {
            return ( this->I_==0 ? ia           // this->I_ is nullptr
                                 : this->I_[ia] // this->I_ is nonzero pointer
                   );
        }

        void add_cell( hilbert::HilbertIndex_t h0 )
        {// update verlet lists for h0-h0 interactions (intra-cell)
            int const first_atom_in_h0 = this->hl_[h0][0];
            int const n_atoms_in_h0    = this->hl_[h0][1];
            for( int ia=first_atom_in_h0+1; ia<first_atom_in_h0+n_atoms_in_h0; ++ia )
            {
                int const i = ia2i(ia);
                double xi = this->rx_[i];
                double yi = this->ry_[i];
                double zi = this->rz_[i];
                for( int ja=first_atom_in_h0; ja<ia; ++ja )
                {
                    int const j = ia2i(ja);
                    double const r2 = sq(xi-this->rx_[j]) + sq(yi-this->ry_[j]) + sq(zi-this->rz_[j]);
                    if( r2<=this->rcutoff2_ )
                        this->add(i,j);
                }
            }
        }

        void add_cell_cell( hilbert::HilbertIndex_t h0, hilbert::HilbertIndex_t h1 )
        {// update verlet lists for h0-h1 interactions (h0!=h1)
            if( h0==h1 )
                std::cout<<"oops "<<h0<<"-"<<h1<<std::endl;
            assert(h0!=h1);
            int first_atom_in_h0 = this->hl_[h0][0];
            int n_atoms_in_h0    = this->hl_[h0][1];
            int first_atom_in_h1 = this->hl_[h1][0];
            int n_atoms_in_h1    = this->hl_[h1][1];
            for( int ia=first_atom_in_h0; ia<first_atom_in_h0+n_atoms_in_h0; ++ia )
            {
                int const i = ia2i(ia);
                double xi = this->rx_[i];
                double yi = this->ry_[i];
                double zi = this->rz_[i];
                for( int ja=first_atom_in_h1; ja<first_atom_in_h1+n_atoms_in_h1; ++ja )
                {
                    int const j = ia2i(ja);
                    double const r2 = sq(xi-this->rx_[j]) + sq(yi-this->ry_[j]) + sq(zi-this->rz_[j]);
                    if( r2<=this->rcutoff2_ ) {
                        this->add(i,j);
                    }
                }
            }
        }

        hilbert::HilbertIndex_t h_neighbour( int const ijk[3], int i0E, int i1E) {
            int ijkE[3];
            for( int c=0; c<3; ++c ) {
                ijkE[c] = ijk[c] + E[i0E][i1E][c];
            }
            return hilbert::ijk2h(ijkE);
        }

        void construct()
        {// loop over all cells
            int ijk00[3] = {0,0,0};
            int const nh = hl_.shape()[0];
            for( int h00=0; h00<nh; ++h00 )
            {
             // intra-cell
                this->add_cell(h00);

             // loop over neighbouring cells
                hilbert::h2ijk(h00,ijk00);
                int nb0[5] = {0,7,10,12,13};
                for( int inb0=0; inb0<4; ++inb0 )
                {
                    hilbert::HilbertIndex_t         h0E = this->h_neighbour(ijk00,nb0[inb0],0);
                    if( -1<h0E && h0E<nh ) {
                        for( int inb1=nb0[inb0]; inb1<nb0[inb0+1]; ++inb1 ) {
                            hilbert::HilbertIndex_t h1E = this->h_neighbour(ijk00,inb1,1);
                            if( -1<h1E && h1E<nh ) {
                                this->add_cell_cell(h0E,h1E);
                            }
                        }
                    }
                }
            }
        }
    private:
        int n_atoms_;
        int reserve_;
        double *rx_, *ry_, *rz_;
        int const * const I_;
        numpy_boost<int,2>& vl_;
        numpy_boost<int,2>& hl_;
        double rcutoff2_;
    };
}// namespace helper

void
build_verlet_list
  ( numpy_boost<int   ,2> verlet_list  // inout
  , numpy_boost<int   ,2> hilbert_list // in
  , numpy_boost<double,1> rx           // in
  , numpy_boost<double,1> ry           // in
  , numpy_boost<double,1> rz           // in
  , numpy_boost<int   ,1> I            // in
  , double                rcutoff2     // in
  )
{
    helper::VerletListBuilder vlb( verlet_list, hilbert_list, rx, ry, rz, I, rcutoff2 );
}

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
    numpy_boost_python_register_type<         int           ,2>();
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

    def("build_verlet_list",build_verlet_list);
}




