
#include "multicomplex.hpp"

//#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <cheerp/client.h> //Misc client side stuff
#include <cheerp/clientlib.h> //Complete DOM/HTML5 interface

#include "AssociatedLegendre.hpp"

template <typename elem, int order, typename T> 
elem get_realY
(
	const T m, 
	const multicomplex<elem,order>& Y
)
{
	elem realY;
		
	if(m < 0){
		realY = abs(pow(-1,m) * root_two * Y.imag);}
						
	if(m == 0){
		realY = abs(Y.real);}
					
	if(m > 0){
		realY = abs(pow(-1,m) * root_two * Y.real);}
		
	return realY;
}

double sps(double theta, double phi)
{
	MX0 I(0,phi), i(1,1);//1 + i
	
	double a;
	
	a = (half * (1 - cos(sqrt(3./(two_pi)) * exp(I) * sin(theta)) ) ).real * 2 * root_two;//(sin(SphericalHarmonicY[1, 1])^2)
	return a;
}

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
		
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array1; //bicomplex real 
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array2; //bicomplex imag
		
		
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array3; //tricomplex real.real 
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array4; //tricomplex real.imag 
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array5; //tricomplex imag.real  
		std::array<double, int(2*max*100)*n_grid_lines*4> complex_array6; //tricomplex imag.imag 
		
		static client::Float64Array * vec10 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&complex_array10, complex_array10.size() * sizeof(double));
		
		////
		
		
		std::array<double, 10201*4> Yx;
		std::array<double, 10201*4> Yy;
		std::array<double, 10201*4> Yz;
		
		static client::Float64Array * vec11 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&Yx, Yx.size() * sizeof(double));
		static client::Float64Array * vec12 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&Yy, Yy.size() * sizeof(double));
		static client::Float64Array * vec13 = 
			MakeTypedArray<TypedArrayForPointerType<double>::type>(&Yz, Yz.size() * sizeof(double));
		
		
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
	

		// draw unit sphere points (r=1 center=(0,0,0)) ... your rays directions
		if(mc_index==11)
		{
			MX0 Y1, Y2;
			
			double realY1, realY2;
			const double w=pi;
			int i = 0;
			
			unsigned int l1,l2;
			l1 = 2; //s,p,d,f,g,h,i
					//0,1,2,3,4,5,6
			l2 = 3; 
			
			int m1,m2; //6,5,4,3,2,1,0,-1,-2,-3,-4,-5,-6
			m1 = 1;
			m2 = 0;
			
			double v1=0,v2=0;
			
			AssociatedLegendre al1(l1, abs(m1)); 
			AssociatedLegendre al2(l2, abs(m2));
			
			bool real_spherical_harmonics = true;
			
			MX0 THETA, PHI;
			MX1 THETA2, PHI2;
			mcdv mcdv;
			
			for (double phi=0; phi < two_pi; phi += two_pi/400.)
			{		
				for (double theta=0; theta < pi; theta += pi/100.)   
				{				
					
					THETA.real = theta;
					//THETA = function(THETA, formula);
					PHI.real = phi;
					
					
					mcdv.sh<0>(THETA2, THETA);
					mcdv.sh<0>(PHI2, PHI);
					
					//Y1 = al1.SphericalHarmonic(THETA, PHI);
					Y1 = mcdv.dv<0>(al1.SphericalHarmonic(THETA2, PHI2));
					
					Y2 = al2.SphericalHarmonic(THETA, PHI);
					
					Y1 = function(Y1, formula);
					//Y2 = function(Y2, formula);
					
					Y1 = pow(Y1,3.5);
					
					if(real_spherical_harmonics) {
						realY1 = get_realY(m1, Y1);	
						realY2 = get_realY(m1, Y2);
					}
					
					else {
						realY1 = abs(Y1); realY2 = abs(Y2); 
					 
					} //complex 
					
					double x,y,z;
					x = std::sin(theta) * std::sin(phi);
					y = std::sin(theta) * std::cos(phi);
					z = std::cos(theta);

					//v1 += (pow(x * realY1,2) + pow(y * realY1,2) + pow(z * realY1,2));
					//v2 += (pow(x * realY2,2) + pow(y * realY2,2) + pow(z * realY2,2));
					
					Yx[i] = x * w * realY1;
					Yy[i] = y * w * realY1;
					Yz[i] = z * w * realY1;
					
					//Yx[i+10201] = x * w * realY2+1;
					//Yy[i+10201] = y * w * realY2;
					//Yz[i+10201] = z * w * realY2;
					

					i++;
				}
			}
		//std::cout << "i = " << i << std::endl; 
		
		//std::cout << "v1    = " << v1 << std::endl;
		//std::cout << "v2/v1 = " << v2/v1 << std::endl; 		
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

	AssociatedLegendre al(2,1);
	MX0 theta,phi;

	mcdv mcdv;
	
	theta.real = 0.97;
	phi.real = 0.34;
	
	MX0 lp;
	lp.real = 0.34;
	lp.imag = -1.1298;
	
	auto ap = al.LegendreP(lp,lp,lp);
	std::cout << "LegendreP[0.34-1.1298i,0.34-1.1298i,0.34-1.1298i] = " << ap << std::endl;
	//LegendreP[n,m,x]

	auto SphericalH = al.SphericalHarmonic(theta, phi);
	std::cout << "SphericalHarmonicY[2, 1, 0.97, 0.34] = " << SphericalH << std::endl;
	
	MX1 THETA, PHI;
	PHI.real = phi;
	mcdv.sh<0>(THETA, theta);			
	SphericalH = mcdv.dv<0>(al.SphericalHarmonic(THETA, PHI));
	std::cout << "D[SphericalHarmonicY(2, 1, theta, 0.34), theta],theta = 0.97 = " << SphericalH << std::endl;
	
	THETA.clear();
	THETA.real = theta;
	mcdv.sh<0>(PHI, phi);		
	SphericalH = mcdv.dv<0>(al.SphericalHarmonic(THETA, PHI));
	std::cout << "D[SphericalHarmonicY(2, 1, 0.97, phi), phi],phi = 0.34 = " << SphericalH << std::endl << std::endl;
	 
	 //D[SphericalHarmonicY(2, 1, theta, 0.34), theta],theta = 0.97 = (+ 0.2628322396 + 0.0929734559*i1)
	 //D[SphericalHarmonicY(2, 1, 0.97, phi), phi],phi = 0.34 = (+ 0.1201370976 - 0.3396227679*i1)
	 
	//MX0 x;
	//x.real = 0.7;
	//x.imag = 0.8;
	//std::cout << "csc(0.7 + 0.8i) = " << csc(x) << std::endl;
	
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




