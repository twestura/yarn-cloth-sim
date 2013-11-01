#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"

using namespace ci;
using namespace ci::app;
using namespace std;

class VisualizerApp : public AppNative {
  public:
	void setup();
	void mouseDown( MouseEvent event );	
	void update();
	void draw();
};

void VisualizerApp::setup()
{
}

void VisualizerApp::mouseDown( MouseEvent event )
{
}

void VisualizerApp::update()
{
}

void VisualizerApp::draw()
{
	// clear out the window with black
	gl::clear( Color( 0, 0, 0 ) ); 
}

CINDER_APP_NATIVE( VisualizerApp, RendererGl )
