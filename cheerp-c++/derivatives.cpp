

#include <iostream>
#include <chrono>


#include "/include/multicomplex.hpp"
#include <iomanip>

#include <cheerp/client.h> //Misc client side stuff
#include <cheerp/clientlib.h> //Complete DOM/HTML5 interface

using namespace client;
using namespace cheerp;

MX0 x(0.4, -0.5);
MX5 y;

#include <sstream>

std::stringstream ss;

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

[[cheerp::genericjs]] void setupInputAndDisplay()
{
	
	//Retrieve the <body> element
	HTMLElement * body = document.get_body();
	
	//Create a range input element <input type="range">
	HTMLInputElement * inputSlider = static_cast<HTMLInputElement*>(document.createElement("input") );
	inputSlider->setAttribute("type", "range");
	
	//This sets the default value
	inputSlider->setAttribute("value", "0.5");
	
	//Create a new <h1> element
	HTMLElement * slider = document.createElement("h1");

	
	//use a C++11 lambda to capture the variables we need
	auto cb = [slider, inputSlider]() -> void { 
	
	//x.real += .1; 
	String * text = inputSlider->get_value();
	slider->set_textContent( text );
	
	};
	
	//Call the lambda to set the slider to the initial value
	cb();
	
	//Set up the handler for the input event. Use a C++11 lambda to capture the variables we need
	slider->addEventListener("input", cheerp::Callback(cb));
	
	//Add the new elements to the <body>
	body->appendChild( slider );
	body->appendChild( inputSlider );
}

//This function will be called only after the DOM is fully loaded

class [[cheerp::genericjs]] Graphics
{
	
	// This method is the handler for requestAnimationFrame. The browser will call this
	// in sync with its graphics loop, usually at 60 fps.
	static void rafHandler()
	{
		mainLoop();
		client::requestAnimationFrame(cheerp::Callback(rafHandler));
	}
	
public:
	
	
	
	static void initialize()
	{
        
		sh(y, x);
		
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

		static HTMLElement * body = document.get_body();
		
		auto slider1 = static_cast<HTMLInputElement*>(document.getElementById("myRange1"));
		auto output1 = static_cast<HTMLInputElement*>(document.getElementById("demo1"));
		
		
		
		auto v = slider1->get_value();
		output1->set_textContent( v );
	
		//Create a new <h4> element
		static HTMLElement * aa = document.createElement("h4");
		
		aa->setAttribute("style", "font-size: 25px;" "font-family: Consolas MS;");
		
		auto i = parseInt( v ); 
		x.real = i/100.0;
		
		sh(y, x);
		
		ss << "x = " << x << "result = ";
		
		auto d = dv(sin(sqrt(y)));
		ss << d;
		
		
		aa->set_textContent( ss.str().c_str() );
		
		ss.clear();
		ss.str("");
		
		//Add the new elements to the <body>
		body->appendChild( aa );
		
		//HTMLElement * lbr = document.createElement("<br />");
		//body->appendChild( lbr );
		
		body->appendChild( slider1 );
		

	}
	
};

 
void webMain()
{
	//outputNumberOfElementsToTheConsole();
	
	
	//client::document.addEventListener("DOMContentLoaded",cheerp::Callback(Graphics::loadCallback));
	
	Graphics::initialize();
	
	
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(12);

	auto begin = std::chrono::steady_clock::now();	
	
	auto d = dv(sin(sqrt(y)));
	
	auto end = std::chrono::steady_clock::now(); 
	
	std::cout << "d^5/dx^5(sin(sqrt(x))), x = 0.4 - 0.5i = ";	
	std::cout << d << std::endl << std::endl;
	
	std::cout << "\nduration d^5/dx^5(sin(sqrt(x))) : " << int(
	std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) << " uS\n" << std::endl;
	
}




