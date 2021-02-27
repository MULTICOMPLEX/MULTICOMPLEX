
#include <cheerp/client.h> //Misc client side stuff
#include <cheerp/clientlib.h> //Complete DOM/HTML5 interface

using namespace client;
using namespace cheerp;


int [[cheerp::jsexport]] [[cheerp::genericjs]] factorial(int n)
{
        if (n < 2)
                return 1;
        return n * factorial(n-1);
}

[[cheerp::jsexport]] int add(int a, int b) {
  return a + b;
}

[[cheerp::jsexport]] int sub(int a, int b) {
  return a - b;
}

class [[cheerp::jsexport]] JsStruct
{
private:
        float a;
        int b;
public:
        JsStruct(float _a, int _b):a(_a),b(_b)
        {
				client::console.log("Instance created");
        }
        void test()
        {
                client::console.log("Invoked test");
        }
		
		int factorial(int n)
		{
			if (n < 2)
                return 1;
			return n * factorial(n-1);
		}
};

void loadCallback()
{
 HTMLElement* body=document.get_body();
 auto a = document.getElementById("testId");
 }

void webMain()
{
	
 document.addEventListener("DOMContentLoaded",Callback(loadCallback));
 
}










