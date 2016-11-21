#ifndef MD_HPP
#define MD_HPP

#include <Vc/Vc>

namespace helper
{
 // an alias for the Vc reciprocal
    template <class T, class Abi>
    Vc::Vector<T,Abi> 
    reciprocal( Vc::Vector<T,Abi> t ) {
        return Vc::reciprocal(t);
    }
    
 // specialization for T=float, in which case the Vc reciprocal is too inaccurate 
 // it is replaced by 1/t, which is more expensive, but accurate.
    template <class Abi>
    Vc::Vector<float,Abi> 
    reciprocal( Vc::Vector<float,Abi> t ) {
        return static_cast<float>(1) / t;
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
    template <class T>
    T zero() {
        T t = T();
        return t;
    }
//--------------------------------------------------------------------
    template <class T>
    T one() {
        T t = T(1);
        return t;
    }
 //--------------------------------------------------------------------
}// namespace cst

template
  < class T
  , class Abi=Vc::VectorAbi::Best<T>
  >
struct MD
{//----------------------------------------------------------------------------
    typedef Vc::Vector<T,Abi> T_v; 
 //----------------------------------------------------------------------------
    static inline
    T_v ljpot2( T_v r2 )
    {   
        T_v const r2inv = helper::reciprocal(r2);
        T_v const r6inv = r2inv*r2inv*r2inv;
        T_v const epot = T(4)*r6inv*( r6inv - T(1) );
        return epot;
    }
 //============================================================================   
};

#endif // MD_HPP