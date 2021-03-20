
#include "multicomplex.hpp"

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
	
	
	
	
	
	client::Float64Array * conformal_map(int formula, int mc_index)
    {					
		int N = 200;
		
		int tel = 0;
		int total_index = 0;
		
		constexpr double grid_spacing = 0.25;
		
		constexpr double max = 1.0;
		
		constexpr int n_grid_lines = (2*max)/grid_spacing + 1;
		
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array10;//complex
		
		std::array<double, 94> sphere;
		
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array1; //bicomplex real 
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array2; //bicomplex imag
		
		
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array3; //tricomplex real.real 
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array4; //tricomplex real.imag 
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array5; //tricomplex imag.real  
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array6; //tricomplex imag.imag 
		
		static client::Float64Array * vec10 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&complex_array10, complex_array10.size() * sizeof(double));
		
		////
		
		static client::Float64Array * vec1 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&complex_array1, complex_array1.size() * sizeof(double));
			
		static client::Float64Array * vec2 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&complex_array2, complex_array2.size() * sizeof(double));
			
		
		static client::Float64Array * vec3 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&complex_array3, complex_array3.size() * sizeof(double));
			
		static client::Float64Array * vec4 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&complex_array4, complex_array4.size() * sizeof(double));
			
		static client::Float64Array * vec5 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&complex_array5, complex_array5.size() * sizeof(double));
			
		static client::Float64Array * vec6 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&complex_array6, complex_array6.size() * sizeof(double));
		
		////
		
		std::array<double, 1242> X;
		std::array<double, 1242> Y;
		std::array<double, 1242> Z;
		
		static client::Float64Array * vec11 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&X, X.size() * sizeof(double));
		static client::Float64Array * vec12 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&Y, Y.size() * sizeof(double));
		static client::Float64Array * vec13 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&Z, Z.size() * sizeof(double));
		
	
	// draw unit sphere points (r=1 center=(0,0,0)) ... your rays directions
	
	if(mc_index==11)
	{
		int ia,na,ib,nb;
		int i = 0;
	
		double x,y,z,r;
		double a,b,da,db;
		na=32;                                  // number of slices
		da=pi/double(na-1);                   // latitude angle step
		
		MX0 d;
		MX0 mc;
		
		for (a=-half_pi,ia=0;ia<na;ia++,a+=da) // slice sphere to circles in xy planes
		{
			r=std::cos(a);                           // radius of actual circle in xy plane
			z=std::sin(a);                           // height of actual circle in xy plane
			nb=std::ceil(two_pi*r/da);
			db=two_pi/double(nb);             // longitude angle step
			if ((ia==0)||(ia==na-1)) { nb=1; db=0.0; }  // handle edge cases
			for (b=0.0,ib=0;ib<nb;ib++,b+=db)   // cut circle to vertexes
			{
				x=r*std::cos(b);                     // compute x,y of vertex
				y=r*std::sin(b);

				mc.real = x;
				mc.imag = y;
				
				d = function(mc, formula);	
				
			// this just draw the ray direction (x,y,z) as line in OpenGL
			// so you can ignore this
			// instead add the ray cast of yours
			double w=pi;
       
			//glColor3f(1.0,1.0,1.0); glVertex3d(x,y,z);
			//glColor3f(0.0,0.0,0.0); glVertex3d(w*x,w*y,w*z);
				//X[i] = x * w;
				//Y[i] = y * w;
				
				X[i] = d.real * w;
				Y[i] = d.imag * w;
				
				Z[i] = z * w;
				i++;
			}
		}
		 //std::cout << i << std::endl;
	}

	
		if(mc_index==0 || mc_index==2 || mc_index==10)
		{
		for(double y = -max; y <= max; y+=grid_spacing)
		{
			for(int i = 0; i < N; i++)
			{		
				double t = (i-(max*100))/100.;

				auto index = i+ tel * N;
				
				
				//complex
				if(mc_index==10){
				
				MX0 mc;
				mc.real = t*half_pi;
				mc.imag = y*half_pi;
				
				MX0 d;
				d = function(mc, formula);
				
				complex_array10[index]             = d.real;
				complex_array10[N*n_grid_lines + N*n_grid_lines + index] = d.imag;
				
				mc.real = y*half_pi;
				mc.imag = t*half_pi;
				
				d = function(mc, formula);
				complex_array10[N*n_grid_lines + index]             = d.real;
				complex_array10[N*n_grid_lines + N*n_grid_lines + N*n_grid_lines + index] = d.imag;
				
				}
				
				
				//bicomplex
				if(mc_index==0){

				MX1 mc;
				mc.real.real = t*half_pi;
				mc.real.imag = y*half_pi;
				
				mc.imag.real = t*half_pi;
				mc.imag.imag = y*half_pi;
				
				MX1 d;
				d = function(mc, formula);
				
				complex_array1[index]             = d.real.real;
				complex_array1[N*n_grid_lines + N*n_grid_lines + index] = d.real.imag;
				
				complex_array2[index]             = d.imag.real;
				complex_array2[N*n_grid_lines + N*n_grid_lines + index] = d.imag.imag;
				
				mc.real.real = y*half_pi;
				mc.real.imag = t*half_pi;
				
				mc.imag.real = y*half_pi;
				mc.imag.imag = t*half_pi;
				
				d = function(mc, formula);
				
				complex_array1[N*n_grid_lines + index]             = d.real.real;
				complex_array1[N*n_grid_lines + N*n_grid_lines + N*n_grid_lines + index] = d.real.imag;
				
				
				complex_array2[N*n_grid_lines + index]             = d.imag.real;
				complex_array2[N*n_grid_lines + N*n_grid_lines + N*n_grid_lines + index] = d.imag.imag;
				
				}

				
				//tricomplex
				if(mc_index==2){
				
				MX2 mc;
				mc.real.real.real = t*half_pi;
				mc.real.real.imag = y*half_pi;
				
				mc.real.imag.real = t*half_pi;
				mc.real.imag.imag = y*half_pi;
				
				
				mc.imag.real.real = t*half_pi;
				mc.imag.real.imag = y*half_pi;
				
				mc.imag.imag.real = t*half_pi;
				mc.imag.imag.imag = y*half_pi;
				
				MX2 d;
				d = function(mc, formula);
				
				complex_array3[index]             						= d.real.real.real;
				complex_array3[N*n_grid_lines + N*n_grid_lines + index] = d.real.real.imag;
				
				complex_array4[index]             						= d.real.imag.real;
				complex_array4[N*n_grid_lines + N*n_grid_lines + index] = d.real.imag.imag;
				
				complex_array5[index]             						= d.imag.real.real;
				complex_array5[N*n_grid_lines + N*n_grid_lines + index] = d.imag.real.imag;
				
				complex_array6[index]             						= d.imag.imag.real;
				complex_array6[N*n_grid_lines + N*n_grid_lines + index] = d.imag.imag.imag;
				
				mc.real.real.real = y*half_pi;
				mc.real.real.imag = t*half_pi;
				
				mc.real.imag.real = y*half_pi;
				mc.real.imag.imag = t*half_pi;
				
				mc.imag.real.real = y*half_pi;
				mc.imag.real.imag = t*half_pi;
				
				mc.imag.imag.real = y*half_pi;
				mc.imag.imag.imag = t*half_pi;
				
				d = function(mc, formula);
				
				complex_array3[N*n_grid_lines + index]             						 = d.real.real.real;
				complex_array3[N*n_grid_lines + N*n_grid_lines + N*n_grid_lines + index] = d.real.real.imag;
				
				complex_array4[N*n_grid_lines + index]             						 = d.real.imag.real;
				complex_array4[N*n_grid_lines + N*n_grid_lines + N*n_grid_lines + index] = d.real.imag.imag;
				
				complex_array5[N*n_grid_lines + index]             						 = d.imag.real.real;
				complex_array5[N*n_grid_lines + N*n_grid_lines + N*n_grid_lines + index] = d.imag.real.imag;
				
				complex_array6[N*n_grid_lines + index]            						 = d.imag.imag.real;
				complex_array6[N*n_grid_lines + N*n_grid_lines + N*n_grid_lines + index] = d.imag.imag.imag;
				
				}
				
			}	
			
			tel++;
		}}
		
		
		if(mc_index==0)
			return vec1;
		else if(mc_index==1)
			return vec2;
			
		else if(mc_index==2)
			return vec3;
		else if(mc_index==3)
			return vec4;
		else if(mc_index==4)
			return vec5;
		else if(mc_index==5)
			return vec6;
		
		else if(mc_index==10)
			return vec10;//complex
			
		else if(mc_index==11)
			return vec11;//sphere x
		else if(mc_index==12)
			return vec12;//sphere y
		else 
			return vec13;//sphere z
     
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




