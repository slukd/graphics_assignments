#include "InputManager.h"
// #include "display.h"
#include "game.h"
#include "glm\glm.hpp"

int main(int argc,char *argv[])
{
	const int DISPLAY_WIDTH = 800;
	const int DISPLAY_HEIGHT = 800;
	const float CAMERA_ANGLE = 45.0f;
	const float NEAR = 1.0f;
	const float FAR = 100.0f;

	Game *scn = new Game(CAMERA_ANGLE,(float)DISPLAY_WIDTH/DISPLAY_HEIGHT,NEAR,FAR);
	
	Display display(DISPLAY_WIDTH, DISPLAY_HEIGHT, "OpenGL");
	
	//Init(display);
	
	scn->Init();

	display.SetScene(scn);
	display.AddKeyCallBack(scn->process_key_input);
	display.AddMouseCallBacks(scn->handle_mouse_input, scn->handle_scroll, scn->handle_cursor_movement);

	while(!display.CloseWindow())
	{
		scn->Draw(1,0,scn->BACK,true,false);
		scn->Motion();
		display.SwapBuffers();
		display.PollEvents();		
	}
	delete scn;
	return 0;
}
