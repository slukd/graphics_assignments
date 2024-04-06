#include "game.h"
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <thread>
#include<chrono>

using namespace std;

//constants:
const double pi = 3.14159265358979323846;

//parameters:
float RUBIKS_CUBE_SIZE = 3.f;
float CUBE_SIZE = 2.f;
float ROTATION_ANGLE = 45;
float WHOLE_CUBE_ROTATION_ANGLE = 45;
map<Face, float> angles_rotated_relative;
map<Face, float> angles_rotated_absolute;
glm::vec3 cube_x_axis(1, 0, 0);	
glm::vec3 cube_y_axis(0, 1, 0);	
glm::vec3 cube_z_axis(0, 0, 1);

static void printMat(const glm::mat4 mat)
{
	cout<<" matrix:"<<endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout<< mat[j][i]<<" ";
		cout<<endl;
	}
}

Game::Game() : Scene()
{
}

Game::Game(float angle ,float relationWH, float near1, float far1) : Scene(angle,relationWH,near1,far1)
{ 	
}

void Game::Init()
{		

	AddShader("../res/shaders/pickingShader");	
	AddShader("../res/shaders/basicShader");
	
	AddTexture("../res/textures/plane.png",false);

	rubicks_cube = make_cube();

	int shape_indx = 0;
	for (int z = 0; z < RUBIKS_CUBE_SIZE; z++) {
		for (int y = 0; y < RUBIKS_CUBE_SIZE; y++) {
			for (int x = 0; x < RUBIKS_CUBE_SIZE; x++) {
				if (x == 0 && y == 0 && z == 0) {
					AddShape(Cube, -1, TRIANGLES);
					SetShapeTex(0, 0);
				} else {
					AddShapeCopy(0, -1, TRIANGLES);
				}
				rubicks_cube[x][y][z] = shapes[shape_indx];
				float offset = CUBE_SIZE * (RUBIKS_CUBE_SIZE-1) / 2;
				shapes[shape_indx]->MyTranslate(glm::vec3(x*CUBE_SIZE - offset, y*CUBE_SIZE - offset, -1*(z*CUBE_SIZE - offset)), 0);

				shape_indx++;
			}
		}
	}

	char axis_lst[] = {'x', 'y', 'z'}; 
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < RUBIKS_CUBE_SIZE; j++) {
			Face face = make_tuple(axis_lst[i], j);
			angles_rotated_relative[face] = 0;
			angles_rotated_absolute[face] = 0;
		}
	}

	pickedShape = 0;
	
	MoveCamera(0,zTranslate,45);
	MoveCamera(0,xTranslate,8);
	MoveCamera(0,yTranslate,8);

	pickedShape = -1;

	//ReadPixel(); //uncomment when you are reading from the z-buffer
}

void Game::Update(const glm::mat4 &MVP,const glm::mat4 &Model,const int  shaderIndx)
{
	Shader *s = shaders[shaderIndx];
	int r = ((pickedShape+1) & 0x000000FF) >>  0;
	int g = ((pickedShape+1) & 0x0000FF00) >>  8;
	int b = ((pickedShape+1) & 0x00FF0000) >> 16;
	s->Bind();
	s->SetUniformMat4f("MVP", MVP);
	s->SetUniformMat4f("Normal",Model);
	s->SetUniform4f("lightDirection", 0.0f , 0.0f, -1.0f, 0.0f);
	if(shaderIndx == 0)
		s->SetUniform4f("lightColor",r/255.0f, g/255.0f, b/255.0f,1.0f);
	else 
		s->SetUniform4f("lightColor",0.7f,0.8f,0.1f,1.0f);
	s->Unbind();
}

void Game::WhenRotate(float angle_x, float angle_y)
{
	float angle_x_radians = (angle_x/360) * 2 * pi;
	float angle_y_radians = (angle_y/360) * 2 * pi;
	float angle_z_radians = atan2(sin(angle_x_radians) * sin(angle_y_radians), cos(angle_x_radians)); //used euler angles to calculate the angle of rotation around the z axis
	float angle_z = (angle_z_radians/(2 * pi)) * 360;
	
	for (auto const& x: angles_rotated_absolute) {
		if (get<0>(x.first) == 'x')
			angles_rotated_absolute[x.first] += angle_x;
		else if (get<0>(x.first) == 'y')
			angles_rotated_absolute[x.first] += angle_y;
		else if (get<0>(x.first) == 'z'){
			angles_rotated_absolute[x.first] += angle_z;			
		}
	}

	bool should_switch_axes_x = true; 
	for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
		Face face_x = make_tuple('x', i);

		float sign_x = 0.f;
		if (angle_x != 0)
			sign_x = (angle_x/abs(angle_x));

		bool cond = sign_x != 0 && (angles_rotated_absolute[face_x] >= 45 || angles_rotated_absolute[face_x] <= -45);
		should_switch_axes_x &= cond;
		if (cond) {
			vector<vector<vector<Shape*>>> new_rubicks_cube = rubicks_cube;
			rotate_data_structure('x', i, new_rubicks_cube, angle_x);
			rubicks_cube = new_rubicks_cube;
			angles_rotated_absolute[face_x] -= sign_x * 90;
		}
	}

	bool should_switch_axes_y = true;
	for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
		Face face_y = make_tuple('y', i);

		float sign_y = 0.f;
		if (angle_y != 0)
			sign_y = (angle_y/abs(angle_y));

		bool cond = sign_y != 0 && (angles_rotated_absolute[face_y] >= 45 || angles_rotated_absolute[face_y] <= -45);
		should_switch_axes_y &= cond;
		
		if (cond) {
			vector<vector<vector<Shape*>>> new_rubicks_cube = rubicks_cube;
			rotate_data_structure('y', i, new_rubicks_cube, angle_y);
			rubicks_cube = new_rubicks_cube;
			angles_rotated_absolute[face_y] -= sign_y * 90;
		}
	}

	bool should_switch_axes_z = true; 
	for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
		Face face_z = make_tuple('z', i);

		float sign_z = 0.f;
		if (angle_z != 0)
			sign_z = (angle_z/abs(angle_z));

		bool cond = sign_z != 0 && (angles_rotated_absolute[face_z] >= 45 || angles_rotated_absolute[face_z] <= -45);
		should_switch_axes_z &= cond;
		
		if (cond) {
			vector<vector<vector<Shape*>>> new_rubicks_cube = rubicks_cube;
			rotate_data_structure('z', i, new_rubicks_cube, angle_z);
			rubicks_cube = new_rubicks_cube;
			angles_rotated_absolute[face_z] -= sign_z * 90;
		}
	}

	if (should_switch_axes_x) {
		vector<vector<vector<Shape*>>> new_rubicks_cube = rubicks_cube;
		switch_cube_axes('x', angle_x);
	}
	if (should_switch_axes_y) {
		vector<vector<vector<Shape*>>> new_rubicks_cube = rubicks_cube;
		switch_cube_axes('y', angle_y);
	}
	if (should_switch_axes_z) {
		vector<vector<vector<Shape*>>> new_rubicks_cube = rubicks_cube;
		switch_cube_axes('z', angle_z);
	}
}

void Game::WhenTranslate()
{
	
}

void Game::Motion()
{
	if(isActive)
	{
	}
}

void Game::my_key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {

	Game *scn = (Game*)glfwGetWindowUserPointer(window);
	
	if(action == GLFW_PRESS)
	{
		Shape* cube;
		switch (key)
		{			
			case GLFW_KEY_UP:
				WHOLE_CUBE_ROTATION_ANGLE = -1 * abs(WHOLE_CUBE_ROTATION_ANGLE);
				for (auto const& x: angles_rotated_absolute) {
					if (get<0>(x.first) == 'x')
						angles_rotated_absolute[x.first] += WHOLE_CUBE_ROTATION_ANGLE;
				}
				scn->rotate_cube(WHOLE_CUBE_ROTATION_ANGLE, glm::vec3(1, 0, 0)); 
				//rotate 45 deg counter-clockwise around the "real world's" x-axis
				break;
		
			case GLFW_KEY_DOWN:
				WHOLE_CUBE_ROTATION_ANGLE = abs(WHOLE_CUBE_ROTATION_ANGLE);
				for (auto const& x: angles_rotated_absolute) {
					if (get<0>(x.first) == 'x')
						angles_rotated_absolute[x.first] += WHOLE_CUBE_ROTATION_ANGLE;
				}
				scn->rotate_cube(WHOLE_CUBE_ROTATION_ANGLE, glm::vec3(1, 0, 0)); 
				//rotate 45 deg clockwise around the "real world's" x-axis
				break;
			
			case GLFW_KEY_RIGHT:
				WHOLE_CUBE_ROTATION_ANGLE = abs(WHOLE_CUBE_ROTATION_ANGLE);
				for (auto const& x: angles_rotated_absolute) {
					if (get<0>(x.first) == 'y')
						angles_rotated_absolute[x.first] += WHOLE_CUBE_ROTATION_ANGLE;
				}
				scn->rotate_cube(WHOLE_CUBE_ROTATION_ANGLE, glm::vec3(0, 1, 0)); 
				//rotate 45 deg clockwise around the "real world's" y-axis
				break;

			case GLFW_KEY_LEFT:
				WHOLE_CUBE_ROTATION_ANGLE = -1 * abs(WHOLE_CUBE_ROTATION_ANGLE);
				for (auto const& x: angles_rotated_absolute) {
					if (get<0>(x.first) == 'y')
						angles_rotated_absolute[x.first] += WHOLE_CUBE_ROTATION_ANGLE;
				}
				scn->rotate_cube(WHOLE_CUBE_ROTATION_ANGLE, glm::vec3(0, 1, 0)); 
				//rotate 45 deg counter-clockwise around the "real world's" y-axis
				break;

			case GLFW_KEY_C:
				WHOLE_CUBE_ROTATION_ANGLE = -1 * abs(WHOLE_CUBE_ROTATION_ANGLE);
				for (auto const& x: angles_rotated_absolute) {
					if (get<0>(x.first) == 'z')
						angles_rotated_absolute[x.first] += WHOLE_CUBE_ROTATION_ANGLE;
				}
				scn->rotate_cube(WHOLE_CUBE_ROTATION_ANGLE, glm::vec3(0, 0, 1)); 
				//rotate 45 deg clockwise around the "real world's" y-axis
				break;

			case GLFW_KEY_X:
				WHOLE_CUBE_ROTATION_ANGLE = abs(WHOLE_CUBE_ROTATION_ANGLE);
				for (auto const& x: angles_rotated_absolute) {
					if (get<0>(x.first) == 'z')
						angles_rotated_absolute[x.first] += WHOLE_CUBE_ROTATION_ANGLE;
				}
				scn->rotate_cube(WHOLE_CUBE_ROTATION_ANGLE, glm::vec3(0, 0, 1)); 
				//rotate 45 deg clockwise around the "real world's" y-axis
				break;
			case GLFW_KEY_R:
				angles_rotated_relative[make_tuple('x', RUBIKS_CUBE_SIZE-1)] += ROTATION_ANGLE;
				scn->rotate_face(ROTATION_ANGLE, cube_x_axis, RUBIKS_CUBE_SIZE-1, angles_rotated_relative);
				break;
			case GLFW_KEY_L:
				angles_rotated_relative[make_tuple('x', 0)] += ROTATION_ANGLE;
				scn->rotate_face(ROTATION_ANGLE, cube_x_axis, 0, angles_rotated_relative);
				break;
			case GLFW_KEY_U:
				angles_rotated_relative[make_tuple('y', RUBIKS_CUBE_SIZE-1)] += ROTATION_ANGLE;
				scn->rotate_face(ROTATION_ANGLE, cube_y_axis, RUBIKS_CUBE_SIZE-1, angles_rotated_relative);
				break;
			case GLFW_KEY_D:
				angles_rotated_relative[make_tuple('y', 0)] += ROTATION_ANGLE;
				scn->rotate_face(ROTATION_ANGLE, cube_y_axis, 0, angles_rotated_relative);
				break;
			case GLFW_KEY_B:
				angles_rotated_relative[make_tuple('z', RUBIKS_CUBE_SIZE-1)] += ROTATION_ANGLE;
				scn->rotate_face(ROTATION_ANGLE, cube_z_axis, RUBIKS_CUBE_SIZE-1, angles_rotated_relative);
				break;
			case GLFW_KEY_F:
				angles_rotated_relative[make_tuple('z', 0)] += ROTATION_ANGLE;
				scn->rotate_face(ROTATION_ANGLE, cube_z_axis, 0, angles_rotated_relative);
				break;
			case GLFW_KEY_SPACE:
				ROTATION_ANGLE *= -1;
				break;
			case GLFW_KEY_Z:
				ROTATION_ANGLE /= 2;
				WHOLE_CUBE_ROTATION_ANGLE /= 2;
				break;
			case GLFW_KEY_A:
				if (ROTATION_ANGLE <= 90)
					ROTATION_ANGLE *= 2;
				if (WHOLE_CUBE_ROTATION_ANGLE <= 90)
					WHOLE_CUBE_ROTATION_ANGLE *= 2;
				break;
			case GLFW_KEY_G:
				scn->generate(window);
				break;
			case GLFW_KEY_S:
				scn->solve(window);
				break;
			
		default:
			break;
		}
	}
}

void Game::my_mouse_callback(GLFWwindow* window,int button, int action, int mods) 
{
	if(action == GLFW_PRESS )
	{
		Game *scn = (Game*)glfwGetWindowUserPointer(window);
		double x2,y2;
		glfwGetCursorPos(window,&x2,&y2);
		scn->Picking((int)x2,(int)y2);
	}
}

void Game::my_cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
	Game *scn = (Game*)glfwGetWindowUserPointer(window);
	scn->UpdatePosition((float)xpos,(float)ypos);

	// moving cube
	if(glfwGetMouseButton(window,GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
	{
		scn->MouseProccessing(GLFW_MOUSE_BUTTON_RIGHT);
	}

	// rotating cube
	else if(glfwGetMouseButton(window,GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
	{
		scn->MouseProccessing(GLFW_MOUSE_BUTTON_LEFT);
	}

}

void Game::my_scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	Game *scn = (Game*)glfwGetWindowUserPointer(window);
	scn->MoveCamera(0, zTranslate, yoffset);
	
}

void Game::generate(GLFWwindow* window) {
	vector<int> keys{GLFW_KEY_R, GLFW_KEY_L, GLFW_KEY_U, GLFW_KEY_D, GLFW_KEY_B, GLFW_KEY_F};
	std::srand(time(0));
	int num_of_actions = rand() % 10 + 11; // in range 10 - 20

	float prev_rotation_angle = ROTATION_ANGLE;
	ROTATION_ANGLE = 45;
	for (int i = 0; i < num_of_actions; i++) {
		int key_index = rand() % keys.size();
		float sign = pow(-1, (rand() % 2 + 1)); // clockwise or counter-clockwise

		ROTATION_ANGLE *= sign;
		my_key_callback(window, keys[key_index], 0, GLFW_PRESS, 0);
		Draw(1,0,BACK,true,false);
		glfwSwapBuffers(window);
		std::this_thread::sleep_for(std::chrono::milliseconds(500));
		my_key_callback(window, keys[key_index], 0, GLFW_PRESS, 0);
		Draw(1,0,BACK,true,false);
		glfwSwapBuffers(window);
		std::this_thread::sleep_for(std::chrono::milliseconds(500));
	}
	ROTATION_ANGLE = prev_rotation_angle;
}

void Game::solve(GLFWwindow *window) {
	vector<vector<vector<Shape*>>> rubicks_cube_copy = rubicks_cube;
	map<vector<vector<vector<Shape*>>>, tuple<int, tuple<char, int, float, int>>> memo;
	vector<tuple<char, int, float, int>> steps{
		make_tuple('x', (int)RUBIKS_CUBE_SIZE-1, 90.f, GLFW_KEY_R), 
		make_tuple('x', 0, 90.f, GLFW_KEY_L), 
		make_tuple('y', (int)RUBIKS_CUBE_SIZE-1, 90.f, GLFW_KEY_U), 
		make_tuple('y', 0, 90.f, GLFW_KEY_D), 
		make_tuple('z', (int)RUBIKS_CUBE_SIZE-1, 90.f, GLFW_KEY_B), 
		make_tuple('z', 0, 90.f, GLFW_KEY_F),
		make_tuple('x', (int)RUBIKS_CUBE_SIZE-1, -90.f, GLFW_KEY_R), 
		make_tuple('x', 0, -90.f, GLFW_KEY_L), 
		make_tuple('y', (int)RUBIKS_CUBE_SIZE-1, -90.f, GLFW_KEY_U), 
		make_tuple('y', 0, -90.f, GLFW_KEY_D), 
		make_tuple('z', (int)RUBIKS_CUBE_SIZE-1, -90.f, GLFW_KEY_B), 
		make_tuple('z', 0, -90.f, GLFW_KEY_F)
	};

	vector<vector<vector<Shape*>>> prev_rubicks_cube = rubicks_cube;
	recursive_solve(window, memo, steps);
	rubicks_cube = prev_rubicks_cube;

	float prev_rotation_angle = ROTATION_ANGLE;
	int num_of_steps_left = 1;
	while (num_of_steps_left > 0) {
		printf("in while after recursive call\n\n");

		num_of_steps_left = get<0>(memo[rubicks_cube]);
		tuple<char, int, float, int> step = get<1>(memo[rubicks_cube]);
		cout << "num_of_steps: " << num_of_steps_left << "step: " << get<3>(step) << ", angle: " << get<2>(step) << "\n\n";
		
		ROTATION_ANGLE = get<2>(step)/2;

		my_key_callback(window, get<3>(step), 0, GLFW_PRESS, 0);
		Draw(1,0,BACK,true,false);
		glfwSwapBuffers(window);
		std::this_thread::sleep_for(std::chrono::milliseconds(500));
		my_key_callback(window, get<3>(step), 0, GLFW_PRESS, 0);
		Draw(1,0,BACK,true,false);
		glfwSwapBuffers(window);
		std::this_thread::sleep_for(std::chrono::milliseconds(500));
	}
	ROTATION_ANGLE = prev_rotation_angle;
}

int Game::recursive_solve(GLFWwindow* window, map<vector<vector<vector<Shape*>>>, tuple<int, tuple<char, int, float, int>>>& memo, vector<tuple<char, int, float, int>>& steps) {
	bool is_solved = true;
	vector<int> indexes = {0, (int)RUBIKS_CUBE_SIZE-1};
	vector<vector<vector<Shape*>>> prev_rubicks_cube(rubicks_cube);
	// vector<glm::vec3> axes = {glm::vec3(1, 0, 0), glm::vec3(0, 1, 0), glm::vec3(0, 0, 1)};
	printf("before fors\n\n");
	// for (int x : indexes) {
	for (int x : indexes) {
		glm::vec4 curr_x_axis = (rubicks_cube[x][0][0]->get_rot()) * glm::vec4(1, 0, 0, 1);
		for (int y = 0; y < RUBIKS_CUBE_SIZE && is_solved; y++) {
			for (int z = 0; z < RUBIKS_CUBE_SIZE && is_solved; z++) {
				Shape* cube = rubicks_cube[x][y][z];
				glm::vec4 x_axis = cube->get_rot() * glm::vec4(1, 0, 0, 1);
				if (!glm::all(glm::equal(curr_x_axis, x_axis))) {
					is_solved = false;
				}
			}
		}
		cout << "x axis, index: " << x << "   is_solved: " << is_solved << endl << endl; 
	}
	for (int y : indexes) {
		glm::vec4 curr_y_axis = (rubicks_cube[0][y][0]->get_rot()) * glm::vec4(0, 1, 0, 1);
		for (int z = 0; z < RUBIKS_CUBE_SIZE && is_solved; z++) {
			for (int x = 0; x < RUBIKS_CUBE_SIZE && is_solved; x++) {
				Shape* cube = rubicks_cube[x][y][z];
				glm::vec4 y_axis = cube->get_rot() * glm::vec4(0, 1, 0, 1);
				if (!glm::all(glm::equal(curr_y_axis, y_axis))) {
					is_solved = false;
				}
			}
		}
		cout << "y axis, index: " << y << "   is_solved: " << is_solved << "\n\n"; 
	}
	for (int z : indexes) {
		glm::vec4 curr_z_axis = (rubicks_cube[0][0][z]->get_rot()) * glm::vec4(0, 0, 1, 1);
		for (int y = 0; y < RUBIKS_CUBE_SIZE && is_solved; y++) {
			for (int x = 0; x < RUBIKS_CUBE_SIZE && is_solved; x++) {
				Shape* cube = rubicks_cube[x][y][z];
				glm::vec4 z_axis = cube->get_rot() * glm::vec4(0, 0, 1, 1);
				if (!glm::all(glm::equal(curr_z_axis, z_axis))) {
					is_solved = false;
				}
			}
		}
		cout << "z axis, index: " << z << "   is_solved: " << is_solved << "\n\n"; 
	}

	printf("before fors, before the end\n\n");


	if (is_solved) {
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
		
		memo[rubicks_cube] = make_tuple(0, make_tuple('x', 0, 0, 0));
		printf("after the end1\n\n");
		rubicks_cube = prev_rubicks_cube;
		return 0;
	} else if (memo.count(rubicks_cube) > 0) {
		// printf("accessing map, if2\n\n");
		int output = get<0>(memo[rubicks_cube]); 
		printf("after the end2\n\n");
		rubicks_cube = prev_rubicks_cube;
		return output;
	} else {
		printf("accessing map, if31\n\n");
		memo[rubicks_cube] = make_tuple(INT_MAX, make_tuple('x', 0, 0, 0));
		printf("after memo\n\n");
		int best_choice = INT_MAX;
		tuple<char, int, float, int> best_step;

		for (tuple<char, int, float, int> step: steps) {
			// vector<vector<vector<Shape*>>> rubicks_cube_copy1 = rubicks_cube_copy;
			printf("before rotating data structure\n\n");
			//rotate_data_structure2(get<0>(step), get<1>(step), new_rubicks_cube, rubicks_cube_copy, get<2>(step));
			vector<vector<vector<Shape*>>>* new_rubicks_cube = new vector<vector<vector<Shape*>>>(rubicks_cube);
			rotate_data_structure(get<0>(step), get<1>(step), *new_rubicks_cube, get<2>(step));
			rubicks_cube = *new_rubicks_cube;

			printf("after rotating data structure, before recursive call\n\n");
			
			int curr_choice = 1 + recursive_solve(window, memo, steps);
			printf("after recursive call\n\n");

			if (curr_choice < best_choice) {
				best_choice = curr_choice;
				best_step = step;
			}

			delete(new_rubicks_cube);
		}

		printf("accessing map, if32\n\n");

		memo[rubicks_cube] = make_tuple(best_choice, best_step);
		cout << "best_step: " << get<3>(best_step) << ", angle: " << get<2>(best_step) << "\n\n";
		printf("after the end3\n\n");
		rubicks_cube = prev_rubicks_cube;
		return best_choice;
	}

}

void Game::rotate_cube_axes(glm::vec3 axis, float angle)
{
	glm::mat4 rotation_mat = glm::rotate(glm::mat4(1), angle, axis);
	cube_x_axis = glm::vec3(rotation_mat * glm::vec4(cube_x_axis, 1));
	cube_y_axis = glm::vec3(rotation_mat * glm::vec4(cube_y_axis, 1));
	cube_z_axis = glm::vec3(rotation_mat * glm::vec4(cube_z_axis, 1));		
}

void Game::switch_cube_axes(char axis, float angle)
{
	float sign = angle/abs(angle);
		glm::vec3 tmp;
		if (axis == 'x') {
			tmp = cube_y_axis;
			cube_y_axis = -sign * cube_z_axis;
			cube_z_axis = sign * tmp;

			for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
				Face face_y = make_tuple('y', i);
				Face face_z = make_tuple('z', i);
				int prev_y_angle = angles_rotated_absolute[face_y];
				angles_rotated_absolute[face_y] = -sign * angles_rotated_absolute[face_z];
				angles_rotated_absolute[face_z] = sign * prev_y_angle;
			}

		} else if (axis == 'y') {
			tmp = cube_z_axis;
			cube_z_axis = -sign * cube_x_axis;
			cube_x_axis = sign * tmp;

			for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
				Face face_z = make_tuple('z', i);
				Face face_x = make_tuple('x', i);
				int prev_z_angle = angles_rotated_absolute[face_z];
				angles_rotated_absolute[face_z] = -sign * angles_rotated_absolute[face_x];
				angles_rotated_absolute[face_x] = sign * prev_z_angle;
			}

		} else {
			tmp = cube_x_axis;
			cube_x_axis = -sign * cube_y_axis;
			cube_y_axis = sign * tmp;

			for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
				Face face_x = make_tuple('x', i);
				Face face_y = make_tuple('y', i);
				int prev_x_angle = angles_rotated_absolute[face_x];
				angles_rotated_absolute[face_x] = -sign * angles_rotated_absolute[face_y];
				angles_rotated_absolute[face_y] = sign * prev_x_angle;
			}

		}
}

void Game::rotate_cube(float angle, glm::vec3 axis)
{
	float prev_angle = angle;
	if (angle > 90) {
		angle = 90;
	} else if (angle < -90) {
		angle = -90;
	}

	char rotation_axis = 'z';
	if (glm::all(glm::equal(axis, glm::vec3(1, 0, 0)))) {
		rotation_axis = 'x';
	} else if (glm::all(glm::equal(axis, glm::vec3(0, 1, 0)))) {
		rotation_axis = 'y';
	}
	
	bool rotated_data_structure = true;
	for (int z = 0; z < RUBIKS_CUBE_SIZE; z++) {
		rotated_data_structure &= rotate_face(angle, axis, z, angles_rotated_absolute);
	}

	rotate_cube_axes(axis, WHOLE_CUBE_ROTATION_ANGLE);

	//if a more than 45 deg turn has been completed, switch the faces
	if (rotated_data_structure) {
		switch_cube_axes(rotation_axis, angle);
	}
	
	if (prev_angle > 90) {
		rotate_cube(prev_angle - 90, axis);
	} else if (prev_angle < -90) {
		rotate_cube(prev_angle + 90, axis);
	}
}

vector<vector<vector<Shape*>>> Game::make_cube()
{
	vector<Shape*> line(RUBIKS_CUBE_SIZE);
	vector<vector<Shape*>> face(RUBIKS_CUBE_SIZE, line);
	vector<vector<vector<Shape*>>> cube(RUBIKS_CUBE_SIZE, face);

	return cube;
}

bool Game::rotate_face(float angle, glm::vec3 axis, int index, map<Face, float>& angles_rotated)
{
	float prev_angle = angle;
	if (angle > 90)
		angle = 90;
	else if (angle < -90)
		angle = -90;
	
	char rotation_axis = 'z';
	if ((angles_rotated == angles_rotated_relative &&  glm::all(glm::equal(axis, cube_x_axis))) || (angles_rotated == angles_rotated_absolute && glm::all(glm::equal(axis, glm::vec3(1, 0, 0))))) {
		rotation_axis = 'x';
	} else if ((angles_rotated == angles_rotated_relative &&  glm::all(glm::equal(axis,cube_y_axis))) || (angles_rotated == angles_rotated_absolute &&  glm::all(glm::equal(axis, glm::vec3(0, 1, 0))))) {
		rotation_axis = 'y';
	} else if ((angles_rotated == angles_rotated_relative &&  glm::all(glm::equal(axis, cube_z_axis))) || (angles_rotated == angles_rotated_absolute &&  glm::all(glm::equal(axis, glm::vec3(0, 0, 1))))) {
		rotation_axis = 'z';
	} else {
		cout << "rotate_face: wrong axis. axis sent: " << axis[0] << ", " << axis[1] << ", " << axis[2] << "\n\n";
		return false;
	}
	
	for (int j = 0; j < RUBIKS_CUBE_SIZE; j++) {
		for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
			Shape* cube = rubicks_cube[i][j][index];
			if (rotation_axis == 'x')
				cube = rubicks_cube[index][i][j];
			else if (rotation_axis == 'y')
				cube = rubicks_cube[i][index][j];
						
			glm::mat4 previous_trans = cube->get_trans();
			glm::mat4 trans_mat1 = -previous_trans; 
			glm::vec3 trans_vec1(trans_mat1[3][0], trans_mat1[3][1], trans_mat1[3][2]); //the translation matrix from the cube's position to the center of the rubick's cube
			cube->MyTranslate(trans_vec1, 0); //translate to the center of the rubick's cube

			cube->MyRotate(angle, axis, 1); //rotate the cube around the real world's axis (mode 1)

			glm::mat4 trans_mat2 = glm::rotate(glm::mat4(1), angle, axis) * -trans_mat1; //take the reverse of trans_mat1 and rotate it
			glm::vec3 trans_vec2(trans_mat2[3][0], trans_mat2[3][1], trans_mat2[3][2]);							
			cube->MyTranslate(trans_vec2, 0); //translate to the final location
		}
	}

	Face face = make_tuple(rotation_axis, index);
	bool rotated_data_strucrture = angles_rotated[face] >= 45 || angles_rotated[face] <= -45;
	if (rotated_data_strucrture) {
		vector<vector<vector<Shape*>>> new_rubicks_cube = rubicks_cube;
		rotate_data_structure(rotation_axis, index, new_rubicks_cube, angle);
		rubicks_cube = new_rubicks_cube;
		angles_rotated[face] -= (angle/abs(angle)) * 90;
	}

	if (prev_angle > 90) {
		return rotate_face(prev_angle - 90, axis, index, angles_rotated);
	} else if (prev_angle < -90) {
		return rotate_face(prev_angle + 90, axis, index, angles_rotated);
	}
	
	return rotated_data_strucrture;
}

void Game::rotate_data_structure(char axis, int index, vector<vector<vector<Shape*>>>& new_rubicks_cube, float angle)
{
	for (int j = 0; j < RUBIKS_CUBE_SIZE; j++) {
		for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
			if (axis == 'x') {
				Shape* cube = rubicks_cube[index][i][j];
				if (angle < 0)
					new_rubicks_cube[index][RUBIKS_CUBE_SIZE-1 - j][i] = cube;
				else
					new_rubicks_cube[index][j][RUBIKS_CUBE_SIZE-1 - i] = cube;
			} else if (axis == 'y') {
				Shape* cube = rubicks_cube[i][index][j];
				if (angle > 0)
					new_rubicks_cube[RUBIKS_CUBE_SIZE-1 - j][index][i] = cube;
				else 
					new_rubicks_cube[j][index][RUBIKS_CUBE_SIZE-1 - i] = cube;
			} else {
				Shape* cube = rubicks_cube[i][j][index];
				if (angle > 0)
					new_rubicks_cube[RUBIKS_CUBE_SIZE-1 - j][i][index] = cube;
				else
					new_rubicks_cube[j][RUBIKS_CUBE_SIZE-1 - i][index] = cube;
			} 
		}
	}
}

Game::~Game(void)
{
}
