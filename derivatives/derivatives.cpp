
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

class [[cheerp::jsexport]] [[cheerp::genericjs]] Graphics
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
	
	inline auto update_array()
    {					
		static client::Float64Array * vec = MakeTypedArray<TypedArrayForPointerType<double>::type>(&va, va.size() * sizeof(double));
		
		sh(y, x);
		
		auto d = dv(sin(y));
		
		va[0] = d.real;
		va[1] = d.imag;
		
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
		static HTMLElement * h1 = document.createElement("h1");
		
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
		
		x.real = 2 * pi * (i1/1000.0);
		x.imag = 2 * pi * (i2/1000.0);
		
		ss << "x = " << x << " = ";
		
		sh(y, x);
		
		d = dv(sin((y)));
		ss << d;
		
		//Add the new elements to the <body>
		h1->set_textContent( ss.str().c_str() );
		
		ss.clear();
		ss.str("");
		
		body->appendChild( h1 );

		//body->appendChild( slider1 );

		//body->appendChild( slider2 );

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

	std::cout << "multicomplex roots : 2 x^2 - 10 x + 5" << std::endl <<  std::endl;
	
	std::cout << "initial : " << sx << std::endl << std::endl;
		
	
	//auto fx1 = [](const auto& x) { return Wilkinsons_polynomial(x, 20); };
	auto fx1 = [](const auto& x) { return 2 * pow(x, 2) - 10 * x + 5; };

	root(fx1, sx, 40);
	
	std::cout << std::endl;
	
	auto begin = std::chrono::steady_clock::now();	
	
	
	sh(y, x);
	
	auto d = dv(sin(sqrt(y)));
	
	auto end = std::chrono::steady_clock::now(); 
	
	std::cout << "d^5/dx^5(sin(sqrt(x))), x = 0.4 - 0.5i = ";	
	std::cout << d << std::endl << std::endl;
	
	std::cout << "duration : " << int(
	std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) << " uS\n" << std::endl;
	
	
	Graphics::initialize();//sliders
	
}




