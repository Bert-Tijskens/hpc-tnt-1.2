#include "md.hpp"
#include <iostream>

#define TEST
#ifdef  TEST
#   include <cmath>
#   include <cassert>
#   include <typeinfo>
#   include <boost/simd/function/all.hpp>

    namespace test
    {//------------------------------------------------------------------------
        template <class T>
        struct LJpot2
        {
            static
            void test()
            {            
                static_assert(std::is_floating_point<T>::value,"oops");
                
                std::cout<<"test::ljpot2<T="<<typeid(T).name()<<"> ... "<<std::flush;
                
                T r2 = T(1);
                T vone = md::ljpot2(r2);
                assert( vone == T() );
    
                T rmin = cst::rmin<T>();
                T vmin = md::ljpot2(rmin*rmin);
                assert( vmin == T(-1) );
    
                std::cout<<"done."<<std::endl;
            }
        };
     //------------------------------------------------------------------------
        template <class T>
        struct LJpot2<bs::pack<T>>
        {
            static
            void test()
            {
                typedef bs::pack<T> pack_t;
                static_assert(std::is_floating_point<T>::value,"oops");
                
                std::cout<<"test::ljpot2<T="<<typeid(pack_t).name()<<"> ... "<<std::flush;
                
                pack_t r2( T(1) );
                pack_t vone = md::ljpot2(r2);
                assert( bs::all( vone == T() ) );
    
                pack_t rmin( cst::rmin<T>() );
                pack_t vmin = md::ljpot2(rmin*rmin);
                assert( bs::all( vmin == T(-1) ) );
    
                std::cout<<"done."<<std::endl;
            }
        };
     //------------------------------------------------------------------------
    }// namespace test
#endif // TEST

#include <boost/simd/function/aligned_load.hpp>
#include <boost/simd/function/aligned_store.hpp>
#include <boost/simd/function/sum.hpp>

#include <boost/simd/memory/allocator.hpp>

// too slow
// #include <random>
// std::default_random_engine dre;
#include "fastrand.hpp"

#include <vector>
#include <chrono>
 //------------------------------------------------------------------------
    template <class T>
    struct CPUTime
    {//------------------------------------------------------------------------
     //------------------------------------------------------------------------
        static
        float ljpot2( size_t n, size_t niter=1, bool test=false)
        {            
            static_assert(std::is_floating_point<T>::value,"oops");
            
            typedef bs::pack<T> pack_t;
            std::cout<<"CPUTime<T="<<typeid(T).name()<<">::ljpot2( n="<<n<<", bool test="<<test<<" )"<<std::flush;
            std::vector<T, bs::allocator<T>> r2(n+pack_t::static_size)
                                           , v (n+pack_t::static_size);
         // initialize
            if( test ) { // unit test
                T rmin2( cst::rmin<T>() * cst::rmin<T>() );
                for (size_t i = 0; i < r2.size(); ++i) {
                    r2[i] = rmin2;
                }
            } else {
//                 std::uniform_real_distribution<T> rnd;
//                 for( size_t i = 0; i < r2.size(); ++i) {
//                     r2[i] = rnd(dre) + T(1);
//                 }
                for( size_t i = 0; i < r2.size(); ++i) {
                    T rnd = fast_uniform_real<T>();
//                     std::cout<<rnd<<std::endl;
                    r2[i] = rnd + T(1);
                }
            }
    
            float cputime = std::numeric_limits<float>::max();
            for(size_t iter=0; iter<niter; ++iter )
            {
                auto t0 = std::chrono::high_resolution_clock::now();
                for (size_t i=0; i<n; i+=pack_t::static_size ) {
                    auto r2_pack = bs::aligned_load<pack_t>( &r2[i] );
                    pack_t v_pack = md::ljpot2( r2_pack );
                    bs::aligned_store( v_pack, &v[i] );
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
                for( size_t i=0; i<n; ++i ) {
//                     std::cout<<i<<": "<<r2[i]<<"->"<<v[i]<<std::endl;
                    T const& vv = v[i];
                    assert( vv == T(-1) );
                }
            } else {
                std::cout<<cputime<<" µs"<<std::endl;
                std::cout<<( static_cast<float>(n)/cputime )*1e6<<" interactions/s"<<std::endl;
            }
            std::cout<<"done."<<std::endl;
            return cputime;
        }
     //------------------------------------------------------------------------
        static
        float epot( size_t n, size_t niter=1, bool test=false)
        {            
            static_assert(std::is_floating_point<T>::value,"oops");
            
            typedef bs::pack<T> pack_t;
            std::cout<<"CPUTime<T="<<typeid(T).name()<<">::epot( n="<<n<<", bool test="<<test<<" )"<<std::flush;
            T w0, rmin;
            if( test ) {
                w0 = 0;
                rmin  = cst::rmin<T>();
            } else {
                w0 = static_cast<T>(-0.5);
            }
         // initialize
            T zero = T(); // = 0.0
            pack_t ep( zero );
            pack_t x0( w0 );
            pack_t y0( w0 );
            pack_t z0( w0 );
            std::vector<T, bs::allocator<T>> x(n+pack_t::static_size)
                                           , y(n+pack_t::static_size)
                                           , z(n+pack_t::static_size);
            for( size_t i = 0; i < x.size(); ++i ) {
                if( test ) {
                    x[i] = rmin;
                    y[i] = T(0);
                    z[i] = T(0);                    
                } else {
                    T rnd = fast_uniform_real<T>();
                    x[i]  = rnd + T(1);
                    rnd   = fast_uniform_real<T>();
                    y[i]  = rnd + T(1);
                    rnd   = fast_uniform_real<T>();
                    z[i]  = rnd + T(1);
                }
            }
            T epot0 = T(0);
            float cputime = std::numeric_limits<float>::max();
            for(size_t iter=0; iter<niter; ++iter )
            {
                auto t0 = std::chrono::high_resolution_clock::now();
                for (size_t i=0; i<n; i+=pack_t::static_size ) {
                    auto  xi = bs::aligned_load<pack_t>( &x[i] );
                    auto  yi = bs::aligned_load<pack_t>( &y[i] );
                    auto  zi = bs::aligned_load<pack_t>( &z[i] );
                    auto dx = xi - x0;
                    auto dy = yi - y0;
                    auto dz = zi - z0;
                    ep += md::ljpot2( dx*dx + dy*dy + dz*dz );
                }
                epot0 += bs::sum( ep );
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
            std::cout<<"epot = "<<epot0<<std::endl;
            if( test ) {
                assert( epot0 == -n );
            } else {
                std::cout<<cputime<<" µs"<<std::endl;    
                std::cout<<( static_cast<float>(n)/cputime)*1000000<<" interactions/s"<<std::endl;
            }    
            std::cout<<"done."<<std::endl;
            return cputime;
        }
     //------------------------------------------------------------------------
    };
 //------------------------------------------------------------------------

#include <iostream>
using namespace std;

int main()
{
#ifdef TEST
    test::LJpot2<         double >::test();
    test::LJpot2<bs::pack<double>>::test();
    test::LJpot2<         float  >::test();
    test::LJpot2<bs::pack<float >>::test();
#endif // TEST

    CPUTime<float >::ljpot2(8,true);
    CPUTime<double>::ljpot2(8,true);
    CPUTime<float >:: epot (8,true);
    CPUTime<double>:: epot (8,true);

    size_t niter=100;
    size_t n = 100000000;
    float cputime[4];
    cputime[0] = CPUTime<float >::ljpot2(n,niter);
    cputime[1] = CPUTime<double>::ljpot2(n,niter);
    cputime[2] = CPUTime<float >:: epot (n,niter);
    cputime[3] = CPUTime<double>:: epot (n,niter);
    for( size_t j=0; j<4; ++j ) {
        std::cout<<j<<"  "<<cputime[j]<<" "<<1e6*n/cputime[j]<<std::endl;
    }
	return 0;
}
