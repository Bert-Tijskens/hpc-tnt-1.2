#include <iostream>
#include <vector>

#include <boost/simd/memory/allocator.hpp>
#include <boost/simd/function/load.hpp>
#include <boost/simd/function/aligned_load.hpp>
#include <boost/simd/function/aligned_store.hpp>

#include <typeinfo>
#include <boost/simd/pack.hpp>
#include <boost/simd/function/rec.hpp>
namespace bs = boost::simd;

namespace test {
    template <class T, class U>
    struct Gather
    {
        static void test()
        {
            std::cout<<"test::Gather<T="<<typeid(T).name()<<">::test() ... "<<std::endl;

            size_t n = 8;
            std::vector<T, bs::allocator<T>> data(2*n);
            std::vector<T, bs::allocator<T>> dat2(2*n);
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
            }// idx = { 0 2 4 ... 2*n   1 3 5 ... 2*n-1 }
             // idx2= { 1 3 5 ... 2*n-1 0 2 4 ... 2*n   }
            constexpr size_t N = bs::pack<T>::static_size;
            typedef bs::pack<T,N> T_pack;
            typedef bs::pack<U,N> U_pack;
            for( size_t i=0; i<2*n; i+=U_pack::static_size ) {
                U_pack ip = bs::aligned_load<U_pack>( &idx [i] );
                T_pack p  = bs::        load<T_pack>( &data[0], ip );
                std::cout<<i<<"   p "<<p<<std::endl;
                U_pack ip2= bs::aligned_load<U_pack>( &idx2[i] );
                std::cout<<i<<" ip2 "<<ip2<<std::endl;
                bs::store( p, &dat2[0], ip2 );
            }
            std::cout<<"[ "<<dat2[0];
            for( size_t i=1; i<2*n; ++i ) {
                std::cout<<", "<<dat2[i];
            }   std::cout<<" ]"<<std::endl;
        }
    };
}// namespace test

int main()
{
    test::Gather<float ,std::int32_t>::test();
    test::Gather<double,std::int32_t>::test();
    return 0;
}
