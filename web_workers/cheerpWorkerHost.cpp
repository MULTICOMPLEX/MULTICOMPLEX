#include <cheerp/client.h>
#include <cheerp/clientlib.h>

using namespace client;



class Graphics2
{

private:
	
	// This method is the handler for requestAnimationFrame. The browser will call this
	// in sync with its graphics loop, usually at 60 fps.
	
	
public:
	
	Graphics2()
	{
		
	};

	static void initialize()
	{
		Worker* w = new Worker("cheerpWorker.js");
		
		w->addEventListener("message", cheerp::Callback([](MessageEvent* e) {
                                        console.log((String*)(e->get_data())); }));
        w->postMessage("Hello World");
		
	}	

	
};


void webMain()
{
        Graphics2::initialize();//sliders
        
}