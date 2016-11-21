#ifndef MD_HPP
#define MD_HPP

#include <boost/simd/pack.hpp>
#include <boost/simd/function/rec.hpp>
namespace bs = boost::simd;

#include <type_traits>

namespace helper
{
 // an alias for the Vc reciprocal
    template <class T>
    T
    reciprocal( T t ) {
        static_assert(std::is_floating_point<T>::value,"oops");
        return T(1)/t;
    }
    
 // specialization for T=float, in which case the Vc reciprocal is too inaccurate 
 // it is replaced by 1/t, which is more expensive, but accurate.
    template <class T>
    bs::pack<T>
    reciprocal( bs::pack<T> t ) {
        static_assert(std::is_floating_point<T>::value,"oops");
        return bs::rec(t);
    }
}

namespace cst
{//--------------------------------------------------------------------
    template <class T>
    T rmin() {
        T r = std::pow( static_cast<T>(2), static_cast<T>(1)/static_cast<T>(6) );
        return r;
    }
 //--------------------------------------------------------------------
}// namespace cst

namespace md
{//----------------------------------------------------------------------------
    template <class T>
    T ljpot2( T r2 )
    {   
        T const r2inv = helper::reciprocal(r2);
        T const r6inv = r2inv*r2inv*r2inv;
        T const epot = T(4)*r6inv*( r6inv-T(1) );
        return epot;
    }
 //============================================================================   
}// namespace md

#endif // MD_HPP
