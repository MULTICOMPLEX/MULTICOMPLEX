
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
	
	
	auto conformal_map(int formula, int mc_index)
    {					
		int N = 200;
		
		MX1 mc;
		
		int tel = 0;
		int total_index = 0;
		
		constexpr double grid_spacing = 0.25;
		
		constexpr double max = 1.0;
		
		constexpr int n_grid_lines = (2*max)/grid_spacing + 1;
		
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array1; //bicomplex real 
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array2; //bicomplex imag
		
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array10;//complex
		
		static client::Float64Array * vec1 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&complex_array1, complex_array1.size() * sizeof(double));
			
		static client::Float64Array * vec2 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&complex_array2, complex_array2.size() * sizeof(double));
			
		static client::Float64Array * vec3 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&complex_array10, complex_array10.size() * sizeof(double));
		
	
		for(double y = -max; y <= max; y+=grid_spacing)
		{
			for(int i = 0; i < N; i++)
			{		
				double t = (i-(max*100))/100.;
				
				mc.real.real = t;
				mc.real.imag = y;
				
				mc.imag.real = t;
				mc.imag.imag = y;
		
				MX1 d;
				
				auto index = i+ tel * N;
				
				//bicomplex
				if(mc_index==0){
				
				d = function(mc, formula);
				complex_array1[index]             = d.real.real;
				complex_array1[N*n_grid_lines + N*n_grid_lines + index] = d.real.imag;
				}
				
				//bicomplex
				if(mc_index==1){
				
				d = function(mc, formula);
				complex_array2[index]             = d.imag.real;
				complex_array2[N*n_grid_lines + N*n_grid_lines + index] = d.imag.imag;
				}
				
				//complex
				if(mc_index==10){
				
				mc.imag.real = 0;
				mc.imag.imag = 0;
				d = function(mc, formula);
				
				complex_array10[index]             = d.real.real;
				complex_array10[N*n_grid_lines + N*n_grid_lines + index] = d.real.imag;
				}
				
				/////
				mc.real.real = y;
				mc.real.imag = t;
				
				mc.imag.real = y;
				mc.imag.imag = t;
							

				//bicomplex
				if(mc_index==0){
					d = function(mc, formula);
				complex_array1[N*n_grid_lines + index]             = d.real.real;
				complex_array1[N*n_grid_lines + N*n_grid_lines + N*n_grid_lines + index] = d.real.imag;
				}
				
				//bicomplex
				if(mc_index==1){
					d = function(mc, formula);
				complex_array2[N*n_grid_lines + index]             = d.imag.real;
				complex_array2[N*n_grid_lines + N*n_grid_lines + N*n_grid_lines + index] = d.imag.imag;
				}
				
				//complex
				if(mc_index==10){
				
				mc.imag.real = 0;
				mc.imag.imag = 0;
				d = function(mc, formula);
				
				complex_array10[N*n_grid_lines + index]             = d.real.real;
				complex_array10[N*n_grid_lines + N*n_grid_lines + N*n_grid_lines + index] = d.real.imag;
				}
				
			}
			
			tel++;
		}

		if(mc_index==0)
			return vec1;
		else if(mc_index==1)
			return vec2;
		else 
			return vec3;
     
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




