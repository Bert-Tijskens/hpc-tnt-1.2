#ifndef MD_HPP
#define MD_HPP

#include <Vc/Vc>

namespace md
{//============================================================================
    namespace scalar 
    {//========================================================================
        template <class T>
        inline T ljpot2( T const& r2 )
        {
            T r2inv = 1/r2;
            T r6inv = r2inv*r2inv*r2inv;
            T epot = 4.f*r6inv*(r6inv-1.f);
            return epot;
        }
     //========================================================================
    }
 //----------------------------------------------------------------------------
    template <class T>
    Vc::Vector<T> ljpot2( Vc::Vector<T> r2 )
    {
        typedef Vc::Vector<T> Tv;    
        Tv const r2inv = 1./r2;
        Tv const r6inv = r2inv*r2inv*r2inv;
        Tv const epot = 4.f*r6inv*(r6inv-1.f);
        return epot;
    }
 //============================================================================   
}// namespace md

#endif // MD_HPP