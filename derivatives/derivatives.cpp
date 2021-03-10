
#include <array>

#include "include/multicomplex.hpp"

#include <cheerp/client.h> //Misc client side stuff
#include <cheerp/clientlib.h> //Complete DOM/HTML5 interface

using namespace client;
using namespace cheerp;

MX0 x(0.4, -0.5);
MX0 d;
MX5 y(0,0);

std::stringstream ss;

std::array<double, 4> va {double(pow(2,64))-1,2,3,4};


[[cheerp::genericjs]] void outputNumberOfElementsToTheConsole()
{
        double number = document.getElementsByTagName("*")->get_length();
        console.log("Live elements = ", number);
}

[[cheerp::genericjs]] int domOutput(const char* str)
{
    client::String* s = new client::String(str);
    client::console.log(s);

    return s->get_length();
}


//This function will be called only after the DOM is fully loaded

class [[cheerp::jsexport]] [[cheerp::genericjs]]  Graphics
{

private:
	
	// This method is the handler for requestAnimationFrame. The browser will call this
	// in sync with its graphics loop, usually at 60 fps.
	
		
	static void rafHandler()
	{
		mainLoop();
		client::requestAnimationFrame(cheerp::Callback(rafHandler));
	}
	
	
public:
	
	Graphics()
	{
		
	};
	
	inline auto update_array(double r, double i)
    {					
		
		static client::Float64Array * vec = MakeTypedArray<TypedArrayForPointerType<double>::type>(&va, va.size() * sizeof(double));
		
		MX0 x;
		
		x.real = r;
		x.imag = i;
		
		//sh(y, x);
		
		//auto d = dv(gamma(sin(y)));
		auto d = gamma(x);
		
		va[0] = d.real;
		va[1] = d.imag;
		
		return vec;
     
    }
	
	inline auto normalize(double val, double min, double max, int N) 
	{
		double delta = max - min;
        return (val - min) / delta;
	}
	
		
	inline auto conformal_map(int formula)
    {					
		int N = 200;
		
		
		MX0 mc;
		
		mcdv mcdv;
		
	
		int tel = 0;
		int total_index = 0;
		
		constexpr double grid_spacing = 0.25;
		
		constexpr double max = 1.0;
		
		constexpr int n_grid_lines = (2*max)/grid_spacing + 1;
		
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array;
		
		static client::Float64Array * vec = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&complex_array, complex_array.size() * sizeof(double));
		
		
		auto func1 = [formula](const auto & z0)   { 
			
			MX0 i;
			i.real = 0;
			i.imag = 1;
			
			MX1 z1;
			MX2 z2;
			MX3 z3;
			MX4 z4;
			
				 if(formula>12 && formula <25) sh(z1, z0);
			else if(formula>24 && formula <37) sh(z2, z0);
			else if(formula>36 && formula <48) sh(z3, z0);
			else if(formula>47 && formula <59) sh(z4, z0);
			else if(formula>58) sh(z2, z0);
			
			if(formula == 1) return z0;
			
			else if(formula == 2) return log(z0);
			else if(formula == 3) return sqrt(z0);
			else if(formula == 4) return sin(z0);
			else if(formula == 5) return cos(z0);
			else if(formula == 6) return tan(z0);
			else if(formula == 7) return sinh(z0);
			else if(formula == 8) return cosh(z0);
			else if(formula == 9) return tanh(z0);
			else if(formula == 10) return exp(z0);
			else if(formula == 11) return gamma(z0);
			else if(formula == 12) return LambertW(0,z0);
			
			else if(formula == 13) return dv((z1));
			else if(formula == 14) return dv(log(z1));
			else if(formula == 15) return dv(sqrt(z1));
			else if(formula == 16) return dv(sin(z1));
			else if(formula == 17) return dv(cos(z1));
			else if(formula == 18) return dv(tan(z1));
			else if(formula == 19) return dv(sinh(z1));
			else if(formula == 20) return dv(cosh(z1));
			else if(formula == 21) return dv(tanh(z1));
			else if(formula == 22) return dv(exp(z1));
			else if(formula == 23) return dv(gamma(z1));
			else if(formula == 24) return dv(LambertW(0,z1));
			
			else if(formula == 25) return dv((z2));
			else if(formula == 26) return dv(log(z2));
			else if(formula == 27) return dv(sqrt(z2));
			else if(formula == 28) return dv(sin(z2));
			else if(formula == 29) return dv(cos(z2));
			else if(formula == 30) return dv(tan(z2));
			else if(formula == 31) return dv(sinh(z2));
			else if(formula == 32) return dv(cosh(z2));
			else if(formula == 33) return dv(tanh(z2));
			else if(formula == 34) return dv(exp(z2));
			else if(formula == 35) return dv(gamma(z2));
			else if(formula == 36) return dv(LambertW(0,z2));
			
			else if(formula == 37) return dv(log(z3));
			else if(formula == 38) return dv(sqrt(z3));
			else if(formula == 39) return dv(sin(z3));
			else if(formula == 40) return dv(cos(z3));
			else if(formula == 41) return dv(tan(z3));
			else if(formula == 42) return dv(sinh(z3));
			else if(formula == 43) return dv(cosh(z3));
			else if(formula == 44) return dv(tanh(z3));
			else if(formula == 45) return dv(exp(z3));
			else if(formula == 46) return dv(gamma(z3));
			else if(formula == 47) return dv(LambertW(0,z3));
			
			else if(formula == 48) return dv(log(z4));
			else if(formula == 49) return dv(sqrt(z4));
			else if(formula == 50) return dv(sin(z4));
			else if(formula == 51) return dv(cos(z4));
			else if(formula == 52) return dv(tan(z4));
			else if(formula == 53) return dv(sinh(z4));
			else if(formula == 54) return dv(cosh(z4));
			else if(formula == 55) return dv(tanh(z4));
			else if(formula == 56) return dv(exp(z4));
			else if(formula == 57) return dv(gamma(z4));
			else if(formula == 58) return dv(LambertW(0,z4));
			
			else if(formula == 59) return dv(gamma(log(z2))) * i;
			else if(formula == 60) return dv(gamma(sqrt(z2))) * i;
			else if(formula == 61) return dv(gamma(sin(z2))) * i;
			else if(formula == 62) return dv(gamma(cos(z2))) * i;
			else if(formula == 63) return dv(gamma(tan(z2))) * i;
			else if(formula == 64) return dv(gamma(sinh(z2))) * i;
			else if(formula == 65) return dv(gamma(cosh(z2))) * i;
			else if(formula == 66) return dv(gamma(tanh(z2))) * i;
			else if(formula == 67) return dv(gamma(exp(z2))) * i;
			else if(formula == 68) return dv(gamma(LambertW(0,z2))) * i;
			
			
			else return dv(gamma(log(z2)));
		};
		
			
				
		for(double y = -max; y <= max; y+=grid_spacing)
		{
			for(int i = 0; i < N; i++)
			{		
				double t = (i-(max*100))/100.;
				
				mc.real = t;
				mc.imag = y;
		
				
				//mcdv.sh<0>(derivative, mc);
		
				MX0 d;
				
				d = func1(mc);
				
				
				//auto d = mcdv.dv<0>( gamma(log(derivative)) )  * r; 
				
			
				auto index = i+ tel * N;
				
				complex_array[index]             = d.real;
				complex_array[N*n_grid_lines + N*n_grid_lines + index] = d.imag;
				
				/////
				mc.real = y;
				mc.imag = t;
		
				
				//mcdv.sh<0>(derivative, mc);
				//d = mcdv.dv<0>( gamma(log(derivative)) )  * r; 
							
				
				d = func1(mc);
				
				
				complex_array[N*n_grid_lines + index]             = d.real;
				complex_array[N*n_grid_lines + N*n_grid_lines + N*n_grid_lines + index] = d.imag;
				
			}
			
			tel++;
		}

		return vec;
     
    }
	
	
	static void loadCallback()
	{
		initialize();
	}
	
	static void initialize()
	{

		ss.setf(std::ios::fixed, std::ios::floatfield);
		ss.precision(12);
		
		
		//domOutput("Hi from loadCallback!");
       
		client::requestAnimationFrame(cheerp::Callback(rafHandler));
		
		//auto inputValue = static_cast<HTMLInputElement*>(document.getElementById("myRange1"))->get_value();
		//console.log("inputSlider = ", inputValue);
		
		//domOutput("Bye from loadCallback!");
		
	}	
	
	static void mainLoop()
	{
		//Retrieve the <body> element
		
		static HTMLElement * body;
		
		body = document.get_body();
		
		//Create a new elements
		//static HTMLElement * h1 = document.createElement("h1");
	
		 static client::Element* titleElement = 
			client::document.getElementById("pagetitle");
		
		//h1->setAttribute("style", "font-size: 25px;" "font-family: Consolas MS;");
		

		auto slider1 = static_cast<HTMLInputElement*>(document.getElementById("myRange1"));
		auto slider2 = static_cast<HTMLInputElement*>(document.getElementById("myRange2"));
		
		auto output1 = static_cast<HTMLInputElement*>(document.getElementById("demo1"));
		auto output2 = static_cast<HTMLInputElement*>(document.getElementById("demo2"));
		
		auto v1 = slider1->get_value();
		//output1->set_textContent( v1 );
		
		auto v2 = slider2->get_value();
		//output2->set_textContent( v2 );
		
		auto i1 = parseInt( v1 );
		auto i2 = parseInt( v2 ); 		
		
		x.real = (i1/1000.0);
		x.imag = (i2/1000.0);
		
		ss << "gamma(z) = " << x << " = ";
		
		sh(y, x);
		
		//d = dv(gamma(sin(y)));
		//ss << d;
		ss << gamma(x);
		
		//Add the new elements to the <body>
		titleElement->set_textContent( ss.str().c_str() );
		
		ss.clear();
		ss.str("");
		
		//body->appendChild( h1 );
		//body->appendChild( slider1 );

	}
	
	static void init()
	{	
		document.addEventListener("DOMContentLoaded",cheerp::Callback(loadCallback));
	}	
	
};


class [[cheerp::jsexport]] [[cheerp::genericjs]] JsStruct
{
private:
        float a;
        int b;
		
		client::Float64Array * vec;
		
public:
        JsStruct(float _a):a(_a),b(0)
        {
            //client::console.log("Instance created");
			vec = MakeTypedArray<TypedArrayForPointerType<double>::type>(&va, va.size() * sizeof(double));
        }
		
		inline auto test()
        {					
			return vec;
            //client::console.log("Invoked test", a, b++);
        }
		
		int factorial(int n)
		{

			if (n < 2)
					return 1;
			return n * factorial(n-1);
		}
	
};


MX2 sx;
 
void webMain()
{
	//outputNumberOfElementsToTheConsole();
	
	//Graphics::init();

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(10);
	

	std::cout << std::endl;

	sx.random(-10,10);

	std::cout << "multicomplex roots, 2 x^2 - 10 x + 5 = 0" << std::endl <<  std::endl;
	
	std::cout << "initial : " << sx << std::endl << std::endl;
		
	
	//auto fx1 = [](const auto& x) { return Wilkinsons_polynomial(x, 20); };
	auto fx1 = [](const auto& x) { return 2 * pow(x, 2) - 10 * x + 5; };

	root(fx1, sx, 40);
	
	std::cout << std::endl;
	
	auto begin = std::chrono::steady_clock::now();	
	
	
	sh(y, x);
	
	auto d = dv(sin(sqrt(y)));
	
	auto end = std::chrono::steady_clock::now(); 
	
	std::cout << "d^5/dz^5(sin(sqrt(z))), z = 0.4 - 0.5i = ";	
	std::cout << d << std::endl << std::endl;
	
	std::cout << "duration : " << int(
	std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) << " uS\n" << std::endl;
	
	
	Graphics::initialize();//sliders
	
}




