// cheerpWorker.cpp: Code that run inside the worker
#include <cheerp/clientlib.h>
#include <cheerp/client.h>
#include <string>

using namespace client;

class [[cheerp::jsexport]] [[cheerp::genericjs]] Graphics
{

private:
	
	// This method is the handler for requestAnimationFrame. The browser will call this
	// in sync with its graphics loop, usually at 60 fps.
	
	
public:
	
	Graphics()
	{
		
	};

	static void initialize()
	{
		addEventListener("message", cheerp::Callback([](MessageEvent* e) {
                               postMessage(e->get_data());
                               postMessage(e->get_data());
                               }));
		
	}	

	
};

void webMain()
{
       
	   Graphics::initialize();//sliders 
	   
}