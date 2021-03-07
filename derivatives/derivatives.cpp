
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
	
		
	inline auto conformal_map(int N)
    {					
		[[cheerp::wasm]] static std::array<double, 200*5*2 *2> complex_array;
		
		static client::Float64Array * vec = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&complex_array, complex_array.size() * sizeof(double));
		
		MX0 mc;
		MX2 derivative;
		
		MX0 r;
		r.real = 0;
		r.imag = 1;
		
		int tel = 0;
		int total_index = 0;
		
		auto func = [](const auto& z) { return dv( gamma(log(z)) ); };
		
		for(double y = -1; y <= 1; y+=0.5)
		{
			for(int i = 0; i < N; i++)
			{		
				double t = (i-100)/100.;
				
				mc.real = t;
				mc.imag = y;
		
				sh(derivative, mc);
		
				auto d = func(derivative) * r;
				
				//auto d = gamma(mc);
				
				auto index = i+ tel * N;
				
				complex_array[index]             = d.real;
				complex_array[N*5 + N*5 + index] = d.imag;
				
				/////
				mc.real = y;
				mc.imag = t;
		
				sh(derivative, mc);
		
				d = func(derivative) * r;
				
				complex_array[N*5 + index]             = d.real;
				complex_array[N*5 + N*5 + N*5 + index] = d.imag;
				
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




