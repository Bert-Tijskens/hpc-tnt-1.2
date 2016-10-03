//============================================================================
// Name        : hilbert_curve_test.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <hilbert.hpp>

template <class T>
void print( char const * s, T const *t ) {
    std::cout << s << ": ";
    for(int i=0; i<8; ++i) {
        std::cout << t[i] << ' ';
    }   std::cout << std::endl;
}

int main()
{
    hilbert::HilbertIndex_t h[8] = {6,5,3,1,8,7,2,4};
    size_t I[8]= {100,100,100,100,100,100,100,100};
    print("h", h);
    print("I", I);
    int H[8] = {6,5,3,1,8,7,2,4};
    print("H", H);
    
    hilbert::insertion_sort(8,h,I);
    print("h", h);
    print("I", I);    
	
	int HH[8];
    hilbert::reorder(8,H,HH,I);
    print("HH", HH);

    hilbert::reorder_inplace(8,H,I);
    print("H ", H);

	return 0;
}
