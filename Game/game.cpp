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
float SINGLE_ROTATION_ANGLE = 45;
float TOTAL_ROTATION_ANGLE = 45;
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

	// Populate the rubicks_cube structure and initialize it with Shapes
rubicks_cube = make_cube();

int shape_index = 0;
for (int z = 0; z < RUBIKS_CUBE_SIZE; z++) {
    for (int y = 0; y < RUBIKS_CUBE_SIZE; y++) {
        for (int x = 0; x < RUBIKS_CUBE_SIZE; x++) {
            // Create a new shape for the first cube and set its texture
            if (x == 0 && y == 0 && z == 0) {
                AddShape(Cube, -1, TRIANGLES);
                SetShapeTex(0, 0);
            } else {
                // Create a copy of the first shape for all other cubes
                AddShapeCopy(0, -1, TRIANGLES);
            }
            // Assign the shape to the current position in the cube array
            rubicks_cube[x][y][z] = shapes[shape_index];
            // Calculate the offset to center the cube in the scene
            float offset = CUBE_SIZE * (RUBIKS_CUBE_SIZE - 1) / 2;
            shapes[shape_index]->MyTranslate(glm::vec3(x * CUBE_SIZE - offset, y * CUBE_SIZE - offset, -1 * (z * CUBE_SIZE - offset)), 0);

            shape_index++;
        }
    }
}

// Initialize rotation angles for each axis and layer of the cube
char axis_list[] = {'x', 'y', 'z'}; 
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < RUBIKS_CUBE_SIZE; j++) {
        Face face = make_tuple(axis_list[i], j);
        angles_rotated_relative[face] = 0; // Initialize relative rotation angles
        angles_rotated_absolute[face] = 0; // Initialize absolute rotation angles
    }
}

pickedShape = 0; // Initialize pickedShape for selection purposes

// Adjust camera positions to ensure the cube is visible and well-centered
MoveCamera(0, zTranslate, 45); // Adjust camera along z-axis
MoveCamera(0, xTranslate, 8);  // Adjust camera along x-axis
MoveCamera(0, yTranslate, 8);  // Adjust camera along y-axis

pickedShape = -1; // Reset picked shape after adjustments

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
	// Convert input angles from degrees to radians
    float angle_x_radians = (angle_x / 360) * 2 * pi;
    float angle_y_radians = (angle_y / 360) * 2 * pi;

    // Calculate the rotation angle around the z-axis using Euler angles
    float angle_z_radians = atan2(sin(angle_x_radians) * cos(angle_y_radians), cos(angle_x_radians));
    float angle_z = (angle_z_radians / (2 * pi)) * 360;

    // Update the absolute rotation angles for x, y, z based on input
	for (auto const& x: angles_rotated_absolute) {
		if (get<0>(x.first) == 'x')
			angles_rotated_absolute[x.first] += angle_x;
		else if (get<0>(x.first) == 'y')
			angles_rotated_absolute[x.first] += angle_y;
		else if (get<0>(x.first) == 'z'){
			angles_rotated_absolute[x.first] += angle_z;			
		}
	}

	bool need_to_switch_axis_x = true; 
	for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
		Face face_x = make_tuple('x', i);

		float sign_x = 0.f;
		if (angle_x != 0)
			sign_x = (angle_x/abs(angle_x));

		bool cond = sign_x != 0 && (angles_rotated_absolute[face_x] >= 45 || angles_rotated_absolute[face_x] <= -45);
		need_to_switch_axis_x &= cond;
		if (cond) {
			vector<vector<vector<Shape*>>> new_rubicks_cube = rubicks_cube;
			rotate_data_structure('x', i, new_rubicks_cube, angle_x);
			rubicks_cube = new_rubicks_cube;
			angles_rotated_absolute[face_x] -= sign_x * 90;
		}
	}

	bool need_to_switch_axis_y = true;
	for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
		Face face_y = make_tuple('y', i);

		float sign_y = 0.f;
		if (angle_y != 0)
			sign_y = (angle_y/abs(angle_y));

		bool cond = sign_y != 0 && (angles_rotated_absolute[face_y] >= 45 || angles_rotated_absolute[face_y] <= -45);
		need_to_switch_axis_y &= cond;
		
		if (cond) {
			vector<vector<vector<Shape*>>> new_rubicks_cube = rubicks_cube;
			rotate_data_structure('y', i, new_rubicks_cube, angle_y);
			rubicks_cube = new_rubicks_cube;
			angles_rotated_absolute[face_y] -= sign_y * 90;
		}
	}

	bool need_to_switch_axis_z = true; 
	for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
		Face face_z = make_tuple('z', i);

		float sign_z = 0.f;
		if (angle_z != 0)
			sign_z = (angle_z/abs(angle_z));

		bool cond = sign_z != 0 && (angles_rotated_absolute[face_z] >= 45 || angles_rotated_absolute[face_z] <= -45);
		need_to_switch_axis_z &= cond;
		
		if (cond) {
			vector<vector<vector<Shape*>>> new_rubicks_cube = rubicks_cube;
			rotate_data_structure('z', i, new_rubicks_cube, angle_z);
			rubicks_cube = new_rubicks_cube;
			angles_rotated_absolute[face_z] -= sign_z * 90;
		}
	}

	if (need_to_switch_axis_x) {
		vector<vector<vector<Shape*>>> new_rubicks_cube = rubicks_cube;
		switch_cube_axes('x', angle_x);
	}
	if (need_to_switch_axis_y) {
		vector<vector<vector<Shape*>>> new_rubicks_cube = rubicks_cube;
		switch_cube_axes('y', angle_y);
	}
	if (need_to_switch_axis_z) {
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
				TOTAL_ROTATION_ANGLE = -1 * abs(TOTAL_ROTATION_ANGLE);
				for (auto const& x: angles_rotated_absolute) {
					if (get<0>(x.first) == 'x')
						angles_rotated_absolute[x.first] += TOTAL_ROTATION_ANGLE;
				}
				scn->rotate_cube(TOTAL_ROTATION_ANGLE, glm::vec3(1, 0, 0)); 
				//rotate 45 deg counter-clockwise around the "real world's" x-axis
				break;
		
			case GLFW_KEY_DOWN:
				TOTAL_ROTATION_ANGLE = abs(TOTAL_ROTATION_ANGLE);
				for (auto const& x: angles_rotated_absolute) {
					if (get<0>(x.first) == 'x')
						angles_rotated_absolute[x.first] += TOTAL_ROTATION_ANGLE;
				}
				scn->rotate_cube(TOTAL_ROTATION_ANGLE, glm::vec3(1, 0, 0)); 
				//rotate 45 deg clockwise around the "real world's" x-axis
				break;
			
			case GLFW_KEY_RIGHT:
				TOTAL_ROTATION_ANGLE = abs(TOTAL_ROTATION_ANGLE);
				for (auto const& x: angles_rotated_absolute) {
					if (get<0>(x.first) == 'y')
						angles_rotated_absolute[x.first] += TOTAL_ROTATION_ANGLE;
				}
				scn->rotate_cube(TOTAL_ROTATION_ANGLE, glm::vec3(0, 1, 0)); 
				//rotate 45 deg clockwise around the "real world's" y-axis
				break;

			case GLFW_KEY_LEFT:
				TOTAL_ROTATION_ANGLE = -1 * abs(TOTAL_ROTATION_ANGLE);
				for (auto const& x: angles_rotated_absolute) {
					if (get<0>(x.first) == 'y')
						angles_rotated_absolute[x.first] += TOTAL_ROTATION_ANGLE;
				}
				scn->rotate_cube(TOTAL_ROTATION_ANGLE, glm::vec3(0, 1, 0)); 
				//rotate 45 deg counter-clockwise around the "real world's" y-axis
				break;

			case GLFW_KEY_C:
				TOTAL_ROTATION_ANGLE = -1 * abs(TOTAL_ROTATION_ANGLE);
				for (auto const& x: angles_rotated_absolute) {
					if (get<0>(x.first) == 'z')
						angles_rotated_absolute[x.first] += TOTAL_ROTATION_ANGLE;
				}
				scn->rotate_cube(TOTAL_ROTATION_ANGLE, glm::vec3(0, 0, 1)); 
				//rotate 45 deg clockwise around the "real world's" y-axis
				break;

			case GLFW_KEY_X:
				TOTAL_ROTATION_ANGLE = abs(TOTAL_ROTATION_ANGLE);
				for (auto const& x: angles_rotated_absolute) {
					if (get<0>(x.first) == 'z')
						angles_rotated_absolute[x.first] += TOTAL_ROTATION_ANGLE;
				}
				scn->rotate_cube(TOTAL_ROTATION_ANGLE, glm::vec3(0, 0, 1)); 
				//rotate 45 deg clockwise around the "real world's" y-axis
				break;
			case GLFW_KEY_R:
				angles_rotated_relative[make_tuple('x', RUBIKS_CUBE_SIZE-1)] += SINGLE_ROTATION_ANGLE;
				scn->rotate_face(SINGLE_ROTATION_ANGLE, cube_x_axis, RUBIKS_CUBE_SIZE-1, angles_rotated_relative);
				break;
			case GLFW_KEY_L:
				angles_rotated_relative[make_tuple('x', 0)] += SINGLE_ROTATION_ANGLE;
				scn->rotate_face(SINGLE_ROTATION_ANGLE, cube_x_axis, 0, angles_rotated_relative);
				break;
			case GLFW_KEY_U:
				angles_rotated_relative[make_tuple('y', RUBIKS_CUBE_SIZE-1)] += SINGLE_ROTATION_ANGLE;
				scn->rotate_face(SINGLE_ROTATION_ANGLE, cube_y_axis, RUBIKS_CUBE_SIZE-1, angles_rotated_relative);
				break;
			case GLFW_KEY_D:
				angles_rotated_relative[make_tuple('y', 0)] += SINGLE_ROTATION_ANGLE;
				scn->rotate_face(SINGLE_ROTATION_ANGLE, cube_y_axis, 0, angles_rotated_relative);
				break;
			case GLFW_KEY_B:
				angles_rotated_relative[make_tuple('z', RUBIKS_CUBE_SIZE-1)] += SINGLE_ROTATION_ANGLE;
				scn->rotate_face(SINGLE_ROTATION_ANGLE, cube_z_axis, RUBIKS_CUBE_SIZE-1, angles_rotated_relative);
				break;
			case GLFW_KEY_F:
				angles_rotated_relative[make_tuple('z', 0)] += SINGLE_ROTATION_ANGLE;
				scn->rotate_face(SINGLE_ROTATION_ANGLE, cube_z_axis, 0, angles_rotated_relative);
				break;
			case GLFW_KEY_SPACE:
				SINGLE_ROTATION_ANGLE *= -1;
				break;
			case GLFW_KEY_Z:
				SINGLE_ROTATION_ANGLE /= 2;
				TOTAL_ROTATION_ANGLE /= 2;
				break;
			case GLFW_KEY_A:
				if (SINGLE_ROTATION_ANGLE <= 90)
					SINGLE_ROTATION_ANGLE *= 2;
				if (TOTAL_ROTATION_ANGLE <= 90)
					TOTAL_ROTATION_ANGLE *= 2;
				break;
			case GLFW_KEY_G:
				scn->generate(window);
				break;

		default:
			break;
		}
	}
}

void Game::my_mouse_callback(GLFWwindow* window, int button, int action, int mods) 
{
    // Only proceed if the mouse button was pressed
    if (action == GLFW_PRESS)
    {
        // Retrieve the Game instance associated with the GLFW window
        Game *scn = static_cast<Game*>(glfwGetWindowUserPointer(window));
        
        // Obtain the current cursor position
        double x2, y2;
        glfwGetCursorPos(window, &x2, &y2);
        
        // Perform picking operation at the cursor location
        // The cursor position is cast to integers to conform with the expected input types
        scn->Picking(static_cast<int>(x2), static_cast<int>(y2));
    }
}

void Game::my_cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    // Retrieve the Game instance from the GLFW window's user pointer
    Game *scn = static_cast<Game*>(glfwGetWindowUserPointer(window));
    
    // Update the position with the current cursor coordinates
    scn->UpdatePosition(static_cast<float>(xpos), static_cast<float>(ypos));

    // Check the state of the mouse buttons to determine the type of processing required
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
    {
        // Process right mouse button action (typically used for moving objects)
        scn->MouseProccessing(GLFW_MOUSE_BUTTON_RIGHT);
    }
    else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
    {
        // Process left mouse button action (typically used for rotating objects)
        scn->MouseProccessing(GLFW_MOUSE_BUTTON_LEFT);
    }
}

void Game::my_scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    // Retrieve the Game instance from the GLFW window's user pointer
    Game *scn = static_cast<Game*>(glfwGetWindowUserPointer(window));
    
    const float zTranslate = 1.0f; // Example definition, adjust based on actual usage

    // Move the camera based on vertical scroll input (yoffset),
    // adjusting the camera's position or zoom level.
    scn->MoveCamera(0, zTranslate * yoffset, 0);
}


void Game::generate(GLFWwindow* window) {
    // Define the keys corresponding to cube rotations
    std::vector<int> keys{GLFW_KEY_R, GLFW_KEY_L, GLFW_KEY_U, GLFW_KEY_D, GLFW_KEY_B, GLFW_KEY_F};

    // Seed the random number generator
    std::srand(static_cast<unsigned int>(time(0)));

    // Generate a random number of actions between 10 and 20
    int num_of_actions = std::rand() % 10 + 10;

    // Store the current rotation angle to restore later
    float prev_rotation_angle = SINGLE_ROTATION_ANGLE;
    SINGLE_ROTATION_ANGLE = 45;  // Set a fixed rotation angle for the generation sequence

    for (int i = 0; i < num_of_actions; i++) {
        int key_index = std::rand() % keys.size();  // Choose a random key from the set
        float sign = pow(-1, std::rand() % 2 + 1);  // Determine the direction of rotation

        SINGLE_ROTATION_ANGLE *= sign;  // Apply the direction to the rotation angle
        my_key_callback(window, keys[key_index], 0, GLFW_PRESS, 0);  // Simulate key press

        // Draw and swap buffers with a delay to visually represent the rotation
        Draw(1, 0, BACK, true, false);
        glfwSwapBuffers(window);
        std::this_thread::sleep_for(std::chrono::milliseconds(500));

        // Optionally repeat or adjust additional actions, like further draws or state updates
        // (currently repeating the key press and draw for demonstration, adjust as needed)
        my_key_callback(window, keys[key_index], 0, GLFW_PRESS, 0);
        Draw(1, 0, BACK, true, false);
        glfwSwapBuffers(window);
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }

    // Restore the original rotation angle after the sequence is complete
    SINGLE_ROTATION_ANGLE = prev_rotation_angle;
}

void Game::rotate_cube_axes(glm::vec3 axis, float angle)
{
    // Create a rotation matrix using the provided axis and angle
    glm::mat4 rotation_matrix = glm::rotate(glm::mat4(1.0f), glm::radians(angle), axis);

    // Apply the rotation matrix to each of the cube's primary axes
    // glm::vec4 is used to convert the vec3 into vec4 with a homogeneous coordinate for proper matrix multiplication
    cube_x_axis = glm::vec3(rotation_matrix * glm::vec4(cube_x_axis, 0.0f));
    cube_y_axis = glm::vec3(rotation_matrix * glm::vec4(cube_y_axis, 0.0f));
    cube_z_axis = glm::vec3(rotation_matrix * glm::vec4(cube_z_axis, 0.0f));		
}

void Game::switch_cube_axes(char axis, float angle)
{
    // Determine the sign of the rotation angle
    float sign = angle / std::abs(angle);
    glm::vec3 tmp;

    // Rotate axes and update corresponding rotation angles in the absolute map
    if (axis == 'x') {
        // Switch Y and Z axes
        tmp = cube_y_axis;
        cube_y_axis = -sign * cube_z_axis;
        cube_z_axis = sign * tmp;

        // Update angles for Y and Z faces
        for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
            Face face_y = std::make_tuple('y', i);
            Face face_z = std::make_tuple('z', i);
            float prev_y_angle = angles_rotated_absolute[face_y];
            angles_rotated_absolute[face_y] = -sign * angles_rotated_absolute[face_z];
            angles_rotated_absolute[face_z] = sign * prev_y_angle;
        }

    } else if (axis == 'y') {
        // Switch Z and X axes
        tmp = cube_z_axis;
        cube_z_axis = -sign * cube_x_axis;
        cube_x_axis = sign * tmp;

        // Update angles for Z and X faces
        for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
            Face face_z = std::make_tuple('z', i);
            Face face_x = std::make_tuple('x', i);
            float prev_z_angle = angles_rotated_absolute[face_z];
            angles_rotated_absolute[face_z] = -sign * angles_rotated_absolute[face_x];
            angles_rotated_absolute[face_x] = sign * prev_z_angle;
        }

    } else if (axis == 'z') {
        // Switch X and Y axes
        tmp = cube_x_axis;
        cube_x_axis = -sign * cube_y_axis;
        cube_y_axis = sign * tmp;

        // Update angles for X and Y faces
        for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
            Face face_x = std::make_tuple('x', i);
            Face face_y = std::make_tuple('y', i);
            float prev_x_angle = angles_rotated_absolute[face_x];
            angles_rotated_absolute[face_x] = -sign * angles_rotated_absolute[face_y];
            angles_rotated_absolute[face_y] = sign * prev_x_angle;
        }
    }
}

void Game::rotate_cube(float angle, glm::vec3 axis)
{
    float prev_angle = angle; // Store the original angle for later use
    
    // Limit the angle to a maximum of 90 degrees or a minimum of -90 degrees
    if (angle > 90) {
        angle = 90;
    } else if (angle < -90) {
        angle = -90;
    }

    // Determine the rotation axis based on the unit vector 'axis'
    char rotation_axis = 'z'; // Default rotation axis is z
    if (glm::all(glm::equal(axis, glm::vec3(1, 0, 0)))) {
        rotation_axis = 'x'; // Rotation around x-axis
    } else if (glm::all(glm::equal(axis, glm::vec3(0, 1, 0)))) {
        rotation_axis = 'y'; // Rotation around y-axis
    }
    
    // Rotate each face of the cube and check if all rotations were successful
    bool rotated_data_structure = true;
    for (int z = 0; z < RUBIKS_CUBE_SIZE; z++) {
        rotated_data_structure &= rotate_face(angle, axis, z, angles_rotated_absolute);
    }

    // Apply accumulated total rotation to the cube's axes
    rotate_cube_axes(axis, TOTAL_ROTATION_ANGLE);

    // Switch the cube's faces if a turn of more than 45 degrees has been completed
    if (rotated_data_structure) {
        switch_cube_axes(rotation_axis, angle);
    }
    
    // Handle cases where the original angle exceeded 90 degrees
    if (prev_angle > 90) {
        rotate_cube(prev_angle - 90, axis); // Adjust for the remainder if angle was above 90
    } else if (prev_angle < -90) {
        rotate_cube(prev_angle + 90, axis); // Adjust for the remainder if angle was below -90
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
	// Clamp the angle to ensure it doesn't exceed +/- 90 degrees for each operation
float prev_angle = angle;
if (angle > 90)
    angle = 90;
else if (angle < -90)
    angle = -90;

// Determine the rotation axis based on the provided axis and the comparison with predefined axes
char rotation_axis = 'z'; // Default rotation axis
if ((angles_rotated == angles_rotated_relative && glm::all(glm::equal(axis, cube_x_axis))) || 
    (angles_rotated == angles_rotated_absolute && glm::all(glm::equal(axis, glm::vec3(1, 0, 0))))) {
    rotation_axis = 'x';
} else if ((angles_rotated == angles_rotated_relative && glm::all(glm::equal(axis, cube_y_axis))) || 
           (angles_rotated == angles_rotated_absolute && glm::all(glm::equal(axis, glm::vec3(0, 1, 0))))) {
    rotation_axis = 'y';
} else if ((angles_rotated == angles_rotated_relative && glm::all(glm::equal(axis, cube_z_axis))) || 
           (angles_rotated == angles_rotated_absolute && glm::all(glm::equal(axis, glm::vec3(0, 0, 1))))) {
    rotation_axis = 'z';
} else {
    return false;
}

// Iterate through the specified face of the Rubik's cube to apply transformations
for (int j = 0; j < RUBIKS_CUBE_SIZE; j++) {
    for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
        // Select the appropriate cube element based on the rotation axis
        Shape* cube = (rotation_axis == 'x') ? rubicks_cube[index][i][j] :
                      (rotation_axis == 'y') ? rubicks_cube[i][index][j] : 
                                               rubicks_cube[i][j][index];

        // Get the current transformation matrix of the cube
        glm::mat4 previous_trans = cube->get_trans();
        // Create a translation vector that moves the cube to the origin
        glm::vec3 trans_vec1(-previous_trans[3][0], -previous_trans[3][1], -previous_trans[3][2]);
        cube->MyTranslate(trans_vec1, 0);

        // Rotate the cube around the real world's axis
        cube->MyRotate(angle, axis, 1);

        // Calculate the new translation matrix post-rotation and move the cube back
        glm::mat4 trans_mat2 = glm::rotate(glm::mat4(1.0f), angle, axis) * glm::translate(glm::mat4(1.0f), -trans_vec1);
        glm::vec3 trans_vec2(trans_mat2[3][0], trans_mat2[3][1], trans_mat2[3][2]);
        cube->MyTranslate(trans_vec2, 0);
    }
}

// Update the angles in the data structure and check for complete rotations
Face face = make_tuple(rotation_axis, index);
bool to_rotate_cube_ds = angles_rotated[face] >= 45 || angles_rotated[face] <= -45;
if (to_rotate_cube_ds) {
    vector<vector<vector<Shape*>>> new_rubicks_cube = rubicks_cube;
rotate_data_structure(rotation_axis, index, new_rubicks_cube, angle);
rubicks_cube = new_rubicks_cube;
angles_rotated[face] -= (angle / std::abs(angle)) * 90; // Adjust the angle post-rotation
}

// Recursively adjust for angles greater than 90 degrees in either direction
if (prev_angle > 90) {
return rotate_face(prev_angle - 90, axis, index, angles_rotated);
} else if (prev_angle < -90) {
return rotate_face(prev_angle + 90, axis, index, angles_rotated);
}

return to_rotate_cube_ds; // Return whether a significant rotation has been made

}

void Game::rotate_data_structure(char axis, int index, vector<vector<vector<Shape*>>>& new_rubicks_cube, float angle)
{
// Loop through each layer of the Rubik's cube on the specified axis
for (int j = 0; j < RUBIKS_CUBE_SIZE; j++) {
    for (int i = 0; i < RUBIKS_CUBE_SIZE; i++) {
        if (axis == 'x') {
            // Access the element at the current layer and indices
            Shape* cube = rubicks_cube[index][i][j];
            // Rotate the layer clockwise or counterclockwise around the x-axis
            if (angle < 0)
                new_rubicks_cube[index][RUBIKS_CUBE_SIZE - 1 - j][i] = cube;
            else
                new_rubicks_cube[index][j][RUBIKS_CUBE_SIZE - 1 - i] = cube;
        } else if (axis == 'y') {
            // Access the element at the current layer and indices
            Shape* cube = rubicks_cube[i][index][j];
            // Rotate the layer clockwise or counterclockwise around the y-axis
            if (angle > 0)
                new_rubicks_cube[RUBIKS_CUBE_SIZE - 1 - j][index][i] = cube;
            else 
                new_rubicks_cube[j][index][RUBIKS_CUBE_SIZE - 1 - i] = cube;
        } else {
            // Access the element at the current layer and indices
            Shape* cube = rubicks_cube[i][j][index];
            // Rotate the layer clockwise or counterclockwise around the z-axis
            if (angle > 0)
                new_rubicks_cube[RUBIKS_CUBE_SIZE - 1 - j][i][index] = cube;
            else
                new_rubicks_cube[j][RUBIKS_CUBE_SIZE - 1 - i][index] = cube;
        }
    }
}
}


Game::~Game(void)
{
}
