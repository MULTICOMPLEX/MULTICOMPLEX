

#include <iostream>
#include <chrono>

#include "/include/multicomplex.hpp"


void webMain()
{
	MX0 x(0.4, -0.5);
	MX5 y;
	sh(y, x);
		
	auto begin = std::chrono::steady_clock::now();	
	
	auto d = dv(sin(sqrt(y)));
	
	auto end = std::chrono::steady_clock::now();
	
	std::cout << "Derivatives, using Multicomplex Numbers, WebAssembly" << std::endl; 
	
	std::cout << "\nduration d^5/dx^5(sin(sqrt(y))) " <<
	std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " uS\n" << std::endl;
	
	std::cout << "d^5/dx^5(sin(sqrt(y))), x = 0.4 - 0.5i = -16.4786 - 18.9303 I =  ";	
	std::cout << d << std::endl << std::endl;
		
}
