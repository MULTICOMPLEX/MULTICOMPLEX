#include <iostream> 
#include <fstream>    
#include <vector>
#include <chrono>
#include <map>
#include <sstream>
#include <Python.h>
#include <dh.hpp>
 
void PyRun_SimpleStringStd(const std::string&);
constexpr const double pi = 3.141592653589793238462643383279502884197169399375105820974;

void goto_url(const std::string&);
void goto_wiki(const std::string&);

template<typename T>
T Cartomancy(const T Range);
void random_wiki();
void geo_wiki();
std::pair<double,double> get_random_location(const double& x0, const double& y0, const int& radius);