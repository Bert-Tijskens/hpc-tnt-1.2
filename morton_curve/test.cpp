#include "octcode.h"
#include <iostream>
#include <cassert>

struct BoundingBox {
    BoundingBox( float a=1.0)
      : lwr_({0,0,0})
      , upr_({a,a,a})
    {}
    float const* lwr() const { return lwr_; }
    float const* upr() const { return upr_; }
    
private:
    float lwr_[3];       
    float upr_[3];       
};

float hp[8][3] = { {0,0,1}
                 , {0,0,0}
                 , {1,0,0}
                 , {1,0,1}
                 , {1,1,1}
                 , {1,1,0}
                 , {0,1,0}
                 , {0,1,1}
                 };
                 
float zp[8][3] = { {0,0,0}
                 , {1,0,0}
                 , {0,1,0}
                 , {1,1,0}
                 , {0,0,1}
                 , {1,0,1}
                 , {0,1,1}
                 , {1,1,1}
                 };
                 
#include <limits> 
  
void test_n(int n)
{
    std::cout<<"test_n(int n="<<n<<")"<<std::endl;
    BoundingBox bb;
    OctCode oc( bb.lwr(), bb.upr() );

    int level = n-2;
    int n3 = n*n*n;
    float d = 1./n;
    float offset = 0.5*d;
    OCTCODE *hh = new OCTCODE[n3];
    for( int i=0; i<n3; ++i) {
        hh[i] = std::numeric_limits<OCTCODE>::max();
    }
    float p[3] = {offset,offset,offset};    
    for( int i=0; i<n; ++i ) {
        p[0] = offset+i*d;
    for( int j=0; j<n; ++j ) {
        p[1] = offset+j*d;
    for( int k=0; k<n; ++k ) {
        p[2] = offset+k*d;
        OCTCODE h = oc.hilbert(p,level);
        std::cout<<'['<<p[0]<<','<<p[1]<<','<<p[2]<<"] -> h="<<h<<std::endl;
        assert(h<n3);
        hh[h] = h;
    }}}
    for( int i=0; i<n3; ++i) {
        assert(hh[i]==i);
    }   
    delete[] hh;
}

int main()
{
    BoundingBox bb;
    OctCode oc( bb.lwr(), bb.upr() );
    
    for( int i=0; i<8; ++i ) {
        float const* p = hp[i];
        OCTCODE h = oc.hilbert(p,0);
        std::cout<<'['<<p[0]<<','<<p[1]<<','<<p[2]<<"] -> h="<<h<<std::endl;
        assert(h==i);        
    }
    std::cout<<std::endl;
    
    for( int i=0; i<8; ++i ) {
        float const* p = zp[i];
        OCTCODE z = oc.zOrder(p,0);
        std::cout<<'['<<p[0]<<','<<p[1]<<','<<p[2]<<"] -> z="<<z<<std::endl;        
        assert(z==i);        
    }
    std::cout<<std::endl;
    
    int n=2;
    for( int i=0; i<3; ++i ) {
        test_n(n);
        n *= 2;
    }
    
    std::cout<<"done."<<std::endl;    
}