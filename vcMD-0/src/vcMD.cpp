
#include "md.hpp"
#include <iostream>

#define MD_TEST
#ifdef MD_TEST
# include <cmath>
# include <cassert>
#endif

int main()
{
#ifdef MD_TEST
    double r2 = 1.;
    assert(md::scalar::ljpot2(r2)==0.);

    double rmin = std::pow(2.d,1.d/6.d);
    assert(md::scalar::ljpot2(rmin*rmin)==-1.d);

    std::cout<<"MD_TEST done."<<std::endl;
#endif
    return 0;
}
