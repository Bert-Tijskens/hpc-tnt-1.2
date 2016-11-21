#include "md.hpp"
#include <iostream>

#define TEST
#ifdef  TEST
#   include <cmath>
#   include <cassert>
#   include <typeinfo>
#   include <vector>
#   include <boost/simd/memory/allocator.hpp>

    namespace bs = boost::simd;

    namespace test
    {//------------------------------------------------------------------------
        template
          < class T
          , class U = std::int32_t
          , class Abi=Vc::VectorAbi::Best<T>
          >
        void gather_scatter()
        {
            size_t n=8;
            typedef Vc::Vector<T> T_v;
            Vc::Memory<T_v> data(2*n);
            Vc::Memory<T_v> dat2(2*n);
            std::vector<U, bs::allocator<U>> idx (2*n);
            std::vector<U, bs::allocator<U>> idx2(2*n);
            for( size_t i=0; i<2*n; ++i ) {
                data[i] = T(i);
                dat2[i] = T(-1);
            }// data = { 1 2 ... 2*n }
            for( size_t i=0; i<n; ++i ) {
                idx [i  ] = T(2*i  );
                idx [i+n] = T(2*i+1);
                idx2[i  ] = T(2*i+1);
                idx2[i+n] = T(2*i  );
            }// idx  = { 0 2 4 ... 2*n   1 3 5 ... 2*n-1 }
             // idx2 = { 1 3 5 ... 2*n-1 0 2 4 ... 2*n   }

            constexpr size_t N = T_v::size();
            for( size_t i=0; i<2*n; i+=N ) {
                T_v p(&data[0],&idx[i]);            // gather
                std::cout<<p<<std::endl;
                for( size_t j=0; j<p.size(); ++j ) {//
                    dat2[idx2[i+j]] = p[j];         // scatter
                }                                   //
//                U_pack ip = bs::aligned_load<U_pack>( &idx [i] );
//                T_pack p  = bs::        load<T_pack>( &data[0], ip );
//                std::cout<<i<<"   p "<<p<<std::endl;
//                U_pack ip2= bs::aligned_load<U_pack>( &idx2[i] );
//                std::cout<<i<<" ip2 "<<ip2<<std::endl;
//                bs::store( p, &dat2[0], ip2 );
            }
            std::cout<<"[ "<<dat2[0];
            for( size_t i=1; i<2*n; ++i ) {
                std::cout<<", "<<dat2[i];
            }   std::cout<<" ]"<<std::endl;


        }
     //------------------------------------------------------------------------
        template
          < class T
          , class Abi=Vc::VectorAbi::Best<T>
          >
        void ljpot2()
        {
            std::cout<<"test::ljpot2<T="<<typeid(T).name()<<",Abi="<<typeid(Abi).name()<<"> ... "<<std::flush;

            typedef MD<T,Abi> md;

            typename md::T_v r2 = cst::one<T>();
            typename md::T_v vone = md::ljpot2(r2);
            assert( Vc::all_of( vone == cst::zero<T>() ) );

            typename md::T_v rmin = cst::rmin<T>();
            typename md::T_v vmin = md::ljpot2(rmin*rmin);
            assert( Vc::all_of( vmin == -cst::one<T>() ) );

            std::cout<<"done."<<std::endl;
        }
     //------------------------------------------------------------------------
    }// namespace test
#endif // TEST


#include <vector>
#include <chrono>

namespace cput
{//----------------------------------------------------------------------------
    template
      < class T
      , class Abi=Vc::VectorAbi::Best<T>
      >
    float ljpot2( size_t n, size_t niter=1, bool test=false )
    {
        typedef MD<T,Abi> md;
        typedef typename md::T_v T_v;
        Vc::Memory<T_v> r2(n);
        Vc::Memory<T_v> v(n);
     // initialize
        if( test ) { // unit test
            T_v rmin2 = ( cst::rmin<T>() * cst::rmin<T>() );
            for (size_t i = 0; i < r2.vectorsCount(); ++i) {
                r2.vector(i) = rmin2;
            }
        } else {
            for( size_t i = 0; i < r2.vectorsCount(); ++i) {
                r2.vector(i) = T_v::Random() + static_cast<T>(1);
            }
            std::cout<<"\nr2.vectorsCount()="<<r2.vectorsCount()
                     <<"\nr2.entriesCount()="<<r2.entriesCount()<<std::endl;
        }

        float cputime = std::numeric_limits<float>::max();
        for(size_t iter=0; iter<niter; ++iter )
        {
            auto t0 = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < r2.vectorsCount(); ++i) {
                v.vector(i) = md::ljpot2( r2.vector(i) );
            }
            auto t1 = std::chrono::high_resolution_clock::now();
            auto a = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
            std::cout<<"::"<<a<<std::endl;
            float cputime_iter = static_cast<float>( std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count() );
            std::cout<<"##"<<iter<<std::endl;
            std::cout<<"::"<<cputime_iter<<std::endl;
            std::cout<<"::"<<cputime<<"->";
            if( cputime_iter<cputime ) {
                cputime = cputime_iter;
            }
            std::cout<<"::"<<cputime<<std::endl;
        }
        if( test ) {
            for( size_t i=0; i<v.entriesCount(); ++i ) {
                T const& vv = v[i];
                assert( vv == -cst::one<T>() );
            }
        } else {
            std::cout<<cputime<<" µs"<<std::endl;
            std::cout<<( static_cast<float>(n)/cputime )*1000000<<" interactions/s"<<std::endl;
        }
        return cputime;
    }
 //----------------------------------------------------------------------------
    template
      < class T
      , class Abi=Vc::VectorAbi::Best<T>
      >
    float epot( size_t n, size_t niter=1, bool const test=false )
    {
        typedef MD<T,Abi> md;
        typedef typename md::T_v T_v;
        T w0, rmin;
        if( test ) {
            w0 = 0;
            rmin  = cst::rmin<T>();
        } else {
            w0 = static_cast<T>(-0.5);
        }
     // initialize
        T_v epot = T(); // = 0.0
        T_v x0 = w0;
        T_v y0 = w0;
        T_v z0 = w0;
        Vc::Memory<T_v> x(n);
        Vc::Memory<T_v> y(n);
        Vc::Memory<T_v> z(n);
        for (size_t i = 0; i < x.vectorsCount(); ++i) {
            if( test ) {
                x.vector(i) = rmin;
            } else {
                x.vector(i) = T_v::Random() + static_cast<T>(1);
            }
        }
        for (size_t i = 0; i < y.vectorsCount(); ++i) {
            if( test ) {
                y.vector(i) = T();
            } else {
                y.vector(i) = T_v::Random() + static_cast<T>(1);
            }
        }
        for (size_t i = 0; i < z.vectorsCount(); ++i) {
            if( test ) {
                z.vector(i) = T();
            } else {
                z.vector(i) = T_v::Random() + static_cast<T>(1);
            }
        }
        T_v r2;
        T epot0 = T();
        float cputime = std::numeric_limits<float>::max();
        for(size_t iter=0; iter<niter; ++iter )
        {
            auto t0 = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < x.vectorsCount(); ++i) {
                T_v dx = x.vector(i) - x0;
                r2 = dx*dx;
                T_v dy = y.vector(i) - y0;
                r2 += dy*dy;
                T_v dz = z.vector(i) - z0;
                r2 += dz*dz;
                epot += md::ljpot2( r2 );
            }
            for( size_t i=0; i<epot.size(); ++i ) {
                epot0 += epot[i];
            }
            auto t1 = std::chrono::high_resolution_clock::now();
            auto a = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
            std::cout<<"::"<<a<<std::endl;
            float cputime_iter = static_cast<float>( std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count() );
            std::cout<<"##"<<iter<<std::endl;
            std::cout<<"::"<<cputime_iter<<std::endl;
            std::cout<<"::"<<cputime<<"->";
            if( cputime_iter<cputime ) {
                cputime = cputime_iter;
            }
            std::cout<<"::"<<cputime<<std::endl;
        }
        if( test ) {
            assert( epot0 == -n*niter );
        } else {
            std::cout<<"\nx.vectorsCount()="<<x.vectorsCount()
                     <<"\nx.entriesCount()="<<x.entriesCount()<<std::endl;
            std::cout<<cputime<<" µs"<<std::endl;

            std::cout<<( static_cast<float>(n)/cputime)*1000000<<" interactions/s"<<std::endl;
        }
        return cputime;
    }
 //----------------------------------------------------------------------------
}// namespace cput

int main()
{
#ifdef TEST
    test::gather_scatter<float>();
    test::ljpot2<double,Vc::VectorAbi::Scalar>();
    test::ljpot2<double                      >();
    test::ljpot2<float ,Vc::VectorAbi::Scalar>();
    test::ljpot2<float                       >();
#endif // TEST

    cput::ljpot2<float >(8,true);
    cput::epot  <float >(8,true);
    cput::ljpot2<double>(8,true);
    cput::epot  <double>(8,true);

    size_t niter=100;
    size_t n=100000000;
    float cputime[4];

    cputime[0] = cput::ljpot2<float >(n,niter);
    cputime[1] = cput::ljpot2<double>(n,niter);
    cputime[2] = cput:: epot <float >(n,niter);
    cputime[3] = cput:: epot <double>(n,niter);

    for( size_t j=0; j<4; ++j ) {
        std::cout<<j<<"  "<<cputime[j]<<" "<<1e6*n/cputime[j]<<std::endl;
    }

    std::cout<<"vcMD main() done."<<std::endl;
    return 0;
}
