/*
 * md_bs.hpp
 *
 *  Created on: 26 Oct 2016
 *      Author: etijskens
 */

#ifndef MD_BS_HPP_
#define MD_BS_HPP_
#include <boost/simd/pack.hpp>
#include <boost/simd/function/rec.hpp>
#include <boost/simd/function/sum.hpp>
namespace bs = boost::simd;

#include <boost/multi_array.hpp>

#include <type_traits>
#include <chrono>

//#define VERBOSE_DEBUG

namespace md
{//-----------------------------------------------------------------------------
    namespace helper
    {//-----------------------------------------------------------------------------
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
     //-----------------------------------------------------------------------------
        template <class T>
        T
        sum( bs::pack<T> t,size_t n ) {
            static_assert(std::is_floating_point<T>::value,"oops");
            T s(0);
            for(size_t i=0; i<n; ++i ) {
                s += t[i];
            }
            return s;
        }
     //-----------------------------------------------------------------------------
    }// namespace helper
 //-----------------------------------------------------------------------------
    template <class T> // T can be a floating point type or a bs::pack<some_floating_point_type>
    T
    ljpot2
      ( T const r2 // squared distance
      )
    {   
        T const r2inv = helper::reciprocal(r2);
        T const r6inv = r2inv*r2inv*r2inv;
        T const epot = T(4)*r6inv*( r6inv-T(1) );
        return epot;
    }
 //-----------------------------------------------------------------------------
    template <class T> // T can be a floating point type or a bs::pack<some_floating_point_type>
    T
    lj_force_factor2
      ( T const r2 // squared distance
      )
    {// multiply inter-particle distance vector with this factor to get the force
        T const ri2 = helper::reciprocal( r2 );
        T const ri6 = ri2*ri2*ri2;
        T const ff = T(48)*ri2*ri6*( ri6 - T(0.5) );
        return ff;
    }
 //-----------------------------------------------------------------------------
    template <class T> // T can be a floating point type
    float // cputime 
    compute_interactions_verlet_list
      ( T const* rx
      , T const* ry
      , T const* rz
      , T* ax
      , T* ay
      , T* az
      , int const* verlet_n // number of pairs for each atom
      , int const* verlet_j // indices of second atom in each pair
      , size_t n_atoms
      , size_t n_pairs
      , size_t n_iter
      , int imax=-1
      )
    {
        size_t n_ATOMS = ( imax==-1 ? n_atoms : imax+1 );

        auto t0 = std::chrono::high_resolution_clock::now();
        typedef bs::pack<T> T_v;
        size_t constexpr N = T_v::static_size;
        typedef bs::pack<std::int32_t,N> I_v;
        for(size_t iter=0; iter<n_iter; ++iter )
        {
            size_t k0=0;
            for( size_t i=1; i<n_ATOMS; ++i ) //outer loop
            {
              #ifdef VERBOSE_DEBUG
                std::cout<<"\nol: i="<<i<<std::endl;
              #endif//VERBOSE_DEBUG

                size_t n_pairs_i = verlet_n[i];
                size_t k = k0;
                if( n_pairs_i > 0 )
                {
                    T_v rxi( rx[i] );
                    T_v ryi( ry[i] );
                    T_v rzi( rz[i] );
                  #ifdef VERBOSE_DEBUG
                    std::cout<<"ol: rxi = "<<rxi<<std::endl;
                    std::cout<<"ol: n_pairs_i="<<n_pairs_i<<std::endl;
                  #endif//VERBOSE_DEBUG
                    if( n_pairs_i >= N)
                    {// loop with full simd vectors
                        for( ; k<=k0+n_pairs_i-N; k+=N )
                        {
                          #ifdef VERBOSE_DEBUG
                            std::cout<<"ol: k="<<k<<std::endl;
                            std::cout<<"il: verlet_j ["
                                    <<verlet_j[k  ]<<' '
                                    <<verlet_j[k+1]<<' '
                                    <<verlet_j[k+2]<<' '
                                    <<verlet_j[k+3]<<"]["
                                    <<rx[verlet_j[k  ]]<<' '
                                    <<rx[verlet_j[k+1]]<<' '
                                    <<rx[verlet_j[k+2]]<<' '
                                    <<rx[verlet_j[k+3]]<<']'<<std::endl;
                          #endif//VERBOSE_DEBUG
                            I_v ip = bs::load<I_v>( &verlet_j[k] );
                            T_v dx = bs::load<T_v>( rx, ip );
                          #ifdef VERBOSE_DEBUG
                            std::cout<<"il: dx = "<<dx<<std::endl;
                          #endif//VERBOSE_DEBUG
                            dx-= rxi;
                          #ifdef VERBOSE_DEBUG
                            std::cout<<"il: dx = "<<dx<<std::endl;
                          #endif//VERBOSE_DEBUG
                            T_v dy = bs::load<T_v>( ry, ip );
                            dy-= ryi;
                            T_v dz = bs::load<T_v>( rz, ip );
                            dz-= rzi;

                            T_v aij = lj_force_factor2( dx*dx + dy*dy + dz*dz );

                            dx *= aij;
                            dy *= aij;
                            dz *= aij;
                         // update particle i acceleration
                            ax[i] += bs::sum(dx);
                            ay[i] += bs::sum(dy);
                            az[i] += bs::sum(dz);
                         // update particle ja acceleration
                            for( size_t l=0; l<N; ++l ) {
                                size_t const j = verlet_j[k+l];
                                ax[j] -= dx[l];
                                ay[j] -= dy[l];
                                az[j] -= dz[l];
                            }
                        }// loop over verlet pairs of atom i
                      #ifdef VERBOSE_DEBUG
                        std::cout<<"il: ax = [ "<<ax[0];
                        for(size_t a=1; a<n_atoms; ++a ) {
                            std::cout<<' '<<ax[a];
                        }   std::cout<<" ]"<<std::endl;
                      #endif//VERBOSE_DEBUG
                    }
                    size_t const n = n_pairs_i % N;
                    k0 += N*(n_pairs_i/N);
                    k = k0;
                    if( n>0 )
                    {// remainder part, not a full simd vector (using mask)
                        I_v ip = bs::load<I_v>( &verlet_j[k] );
                        T_v dx = bs::load<T_v>( rx, ip );
                      #ifdef VERBOSE_DEBUG
                        std::cout<<"il2: n="<<n<<std::endl;
                        std::cout<<"il2: k="<<k<<std::endl;
                        std::cout<<"il2: verlet_j ["
                                <<verlet_j[k  ]<<' '
                                <<verlet_j[k+1]<<' '
                                <<verlet_j[k+2]<<' '
                                <<verlet_j[k+3]<<"]["
                                <<rx[verlet_j[k  ]]<<' '
                                <<rx[verlet_j[k+1]]<<' '
                                <<rx[verlet_j[k+2]]<<' '
                                <<rx[verlet_j[k+3]]<<']'<<std::endl;
                        std::cout<<"il2: dx = "<<dx<<std::endl;
                      #endif//VERBOSE_DEBUG
                        dx-= rxi;
                      #ifdef VERBOSE_DEBUG
                        std::cout<<"il2: dx = "<<dx<<std::endl;
                      #endif//VERBOSE_DEBUG
                        T_v dy = bs::load<T_v>( ry, ip );
                        dy-= ryi;
                        T_v dz = bs::load<T_v>( rz, ip );
                        dz-= rzi;

                        T_v aij = lj_force_factor2( dx*dx + dy*dy + dz*dz );

                        dx *= aij;
                        dy *= aij;
                        dz *= aij;
                      #ifdef VERBOSE_DEBUG
                        std::cout<<"il2: dx = "<<dx<<std::endl;
                        std::cout<<"il2: sum(dx,n) = "<<helper::sum(dx,n)<<std::endl;
                        std::cout<<"il2: updating a[i] i="<<i<<std::endl;
                      #endif//VERBOSE_DEBUG
                     // update particle i acceleration
                        ax[i] += helper::sum(dx,n);
                        ay[i] += helper::sum(dy,n);
                        az[i] += helper::sum(dz,n);
                     // update particle j acceleration
                        for( size_t l=0; l<n; ++l ) {
                            size_t const j = verlet_j[k+l];
                          #ifdef VERBOSE_DEBUG
                            std::cout<<"il2: updating a[j] j="<<j<<std::endl;
                          #endif//VERBOSE_DEBUG
                            ax[j] -= dx[l];
                            ay[j] -= dy[l];
                            az[j] -= dz[l];
                        }
                      #ifdef VERBOSE_DEBUG
                        std::cout<<"il2: ax=[ "<<ax[0];
                        for(size_t a=1; a<n_atoms; ++a ) {
                            std::cout<<' '<<ax[a];
                        }   std::cout<<" ]"<<std::endl;
                      #endif//VERBOSE_DEBUG
                        k0 += N;
                    }
                }
            }// loop over atoms i
        }// loop over iterations
        auto t1 = std::chrono::high_resolution_clock::now();
        return static_cast<float>( std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count() );
    }
 //----------------------------------------------------------------------------
/**
    function compute_interactions_verlet_linear(n_atoms,rx,ry,rz,ax,ay,az,verlet_linear,m,n_iter)
      ! Computes the lennard-jones forces between all atom pairs in the linear verlet list and update
      ! the accelerations.
      ! It is assumed that the accelerations are appropriately initialized (typically by zeros).
      ! Repeat n_iter times (for accurate timing results)
      ! Return the cputime
        integer, parameter :: wp = kind(1.d0) ! double precision identifier
        integer                    ,intent(in)    :: n_atoms,m
        real(wp),dimension(n_atoms),intent(in)    :: rx,ry,rz ! positions
        real(wp),dimension(n_atoms),intent(inout) :: ax,ay,az ! accelerations
        integer ,dimension(m)      ,intent(in)    :: verlet_linear
        integer                    ,intent(in)    :: n_iter   ! number of iterations, usually 1, >1 to get meaningfull results for cputime measurement

        real :: compute_interactions_verlet_linear  ! the cpu time consumed by this subprogram

        real(wp)  :: aij,dx,dy,dz
        real :: cpustart,cpustop
        integer :: iter,ia,ja,ia_pairs,j,k

        call cpu_time(cpustart)
        do iter=1,n_iter
        j=1
        do ia=1,n_atoms
            ia_pairs = verlet_linear(j)
            !DIR$ SIMD
            do k = j+1,j+ia_pairs
                ja = verlet_linear(k)+1 ! +1 since fortran starts counting from 1 !
                dx = rx(ja)-rx(ia)
                dy = ry(ja)-ry(ia)
                dz = rz(ja)-rz(ia)
                aij = lj_force_factor2( dx**2 + dy**2 + dz**2 )
              ! update particle ia acceleration
                ax(ia) = ax(ia) + aij*dx
                ay(ia) = ay(ia) + aij*dy
                az(ia) = az(ia) + aij*dz
              ! update particle ja acceleration
                ax(ja) = ax(ja) - aij*dx
                ay(ja) = ay(ja) - aij*dy
                az(ja) = az(ja) - aij*dz
!                print *,ia,ja,ax(ia),ay(ia),az(ia),ax(ja),ay(ja),az(ja)
            enddo
            j = j + 1 + ia_pairs
        enddo
        enddo
        call cpu_time(cpustop)
        compute_interactions_verlet_linear = cpustop-cpustart
!        print *,cpustart,cpustop,cpustop-cpustart
    end function
**/
 //-----------------------------------------------------------------------------   
}// namespace md

#endif /* MD_BS_HPP_ */
