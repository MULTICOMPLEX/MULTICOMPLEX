

#include <iostream>
#include <chrono>

#include "/include/multicomplex.hpp"
#include <iomanip>

void webMain()
{
	
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(12);
	
	
	MX0 x(0.4, -0.5);
	MX5 y;
	sh(y, x);
		
	auto begin = std::chrono::steady_clock::now();	
	
	auto d = dv(sin(sqrt(y)));
	
	auto end = std::chrono::steady_clock::now(); 
	
	std::cout << "d^5/dx^5(sin(sqrt(x))), x = 0.4 - 0.5i =  ";	
	std::cout << d << std::endl << std::endl;
	
	std::cout << "\nduration d^5/dx^5(sin(sqrt(x))) : " <<
	std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " uS\n" << std::endl;
		
}
