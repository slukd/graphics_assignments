#include "game.h"
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <thread>
#include <chrono>

using namespace std;

// constants:
const double pi = 3.14159265358979323846;

// parameters:
float CUBE_SIZE = 2.f;
float RUBIKS_CUBE_SIZE = 3.f;
float SINGLE_ROTATION_ANGLE = 45;
float TOTAL_ROTATION_ANGLE = 45;

map<Face, float> angles_rotated_relative;
map<Face, float> angles_rotated_absolute;
glm::vec3 cube_x_axis(1, 0, 0);
glm::vec3 cube_y_axis(0, 1, 0);
glm::vec3 cube_z_axis(0, 0, 1);

static void printMat(const glm::mat4 mat)
{
	cout << " matrix:" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << mat[j][i] << " ";
		cout << endl;
	}
}

Game::Game() : Scene()
{
}

Game::Game(float angle, float relationWH, float near1, float far1) : Scene(angle, relationWH, near1, far1)
{
}

void Game::Init()
{

	AddShader("../res/shaders/pickingShader");
	AddShader("../res/shaders/basicShader");

	AddTexture("../res/textures/plane.png", false);

	// Populate the rubicks_cube structure and initialize it with Shapes
	rubicks_cube = make_cube();

	int shape_index = 0;
	for (int z = 0; z < RUBIKS_CUBE_SIZE; z++)
	{
		for (int y = 0; y < RUBIKS_CUBE_SIZE; y++)
		{
			for (int x = 0; x < RUBIKS_CUBE_SIZE; x++)
			{
				// Create a new shape for the first cube and set its texture
				if (x == 0 && y == 0 && z == 0)
				{
					AddShape(Cube, -1, TRIANGLES);
					SetShapeTex(0, 0);
				}
				else
				{
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
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < RUBIKS_CUBE_SIZE; j++)
		{
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

	// ReadPixel(); //uncomment when you are reading from the z-buffer
}

void Game::Update(const glm::mat4 &MVP, const glm::mat4 &Model, const int shaderIndx)
{
	Shader *s = shaders[shaderIndx];
	int r = ((pickedShape + 1) & 0x000000FF) >> 0;
	int g = ((pickedShape + 1) & 0x0000FF00) >> 8;
	int b = ((pickedShape + 1) & 0x00FF0000) >> 16;
	s->Bind();
	s->SetUniformMat4f("MVP", MVP);
	s->SetUniformMat4f("Normal", Model);
	s->SetUniform4f("lightDirection", 0.0f, 0.0f, -1.0f, 0.0f);
	if (shaderIndx == 0)
		s->SetUniform4f("lightColor", r / 255.0f, g / 255.0f, b / 255.0f, 1.0f);
	else
		s->SetUniform4f("lightColor", 0.7f, 0.8f, 0.1f, 1.0f);
	s->Unbind();
}

vector<vector<vector<Shape *>>> Game::make_cube()
{
	vector<Shape *> line(RUBIKS_CUBE_SIZE);
	vector<vector<Shape *>> face(RUBIKS_CUBE_SIZE, line);
	vector<vector<vector<Shape *>>> cube(RUBIKS_CUBE_SIZE, face);

	return cube;
}

void Game::WhenTranslate()
{
}

void Game::Motion()
{
	if (isActive)
	{
	}
}

void Game::process_key_input(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	// Retrieve the Game object from GLFW window user pointer
	Game *scn = static_cast<Game *>(glfwGetWindowUserPointer(window));

	// Check if the key action is a key press
	if (action == GLFW_PRESS)
	{
		// Declaration of variables
		Shape *cube;
		float rotation_direction = 1.0f; // Default rotation direction
		int axis_index = -1;			 // Default axis index

		// Determine the rotation direction and axis index based on the key pressed
		switch (key)
		{
		case GLFW_KEY_UP:
			rotation_direction = -1.0f; // Counter-clockwise rotation
			axis_index = 0;				// X-axis
			break;

		case GLFW_KEY_DOWN:
			rotation_direction = 1.0f; // Clockwise rotation
			axis_index = 0;			   // X-axis
			break;

		case GLFW_KEY_RIGHT:
			rotation_direction = 1.0f; // Clockwise rotation
			axis_index = 1;			   // Y-axis
			break;

		case GLFW_KEY_LEFT:
			rotation_direction = -1.0f; // Counter-clockwise rotation
			axis_index = 1;				// Y-axis
			break;

		case GLFW_KEY_C:
			rotation_direction = -1.0f; // Counter-clockwise rotation
			axis_index = 2;				// Z-axis
			break;

		case GLFW_KEY_X:
			rotation_direction = 1.0f; // Clockwise rotation
			axis_index = 2;			   // Z-axis
			break;

		case GLFW_KEY_R:
			// Update relative rotation angle for rotating the face on the X-axis
			angles_rotated_relative[make_tuple('x', RUBIKS_CUBE_SIZE - 1)] += SINGLE_ROTATION_ANGLE;
			scn->rotate_face(SINGLE_ROTATION_ANGLE, cube_x_axis, RUBIKS_CUBE_SIZE - 1, angles_rotated_relative);
			return;

		case GLFW_KEY_L:
			// Update relative rotation angle for rotating the face on the X-axis
			angles_rotated_relative[make_tuple('x', 0)] += SINGLE_ROTATION_ANGLE;
			scn->rotate_face(SINGLE_ROTATION_ANGLE, cube_x_axis, 0, angles_rotated_relative);
			return;

		case GLFW_KEY_U:
			// Update relative rotation angle for rotating the face on the Y-axis
			angles_rotated_relative[make_tuple('y', RUBIKS_CUBE_SIZE - 1)] += SINGLE_ROTATION_ANGLE;
			scn->rotate_face(SINGLE_ROTATION_ANGLE, cube_y_axis, RUBIKS_CUBE_SIZE - 1, angles_rotated_relative);
			return;

		case GLFW_KEY_D:
			// Update relative rotation angle for rotating the face on the Y-axis
			angles_rotated_relative[make_tuple('y', 0)] += SINGLE_ROTATION_ANGLE;
			scn->rotate_face(SINGLE_ROTATION_ANGLE, cube_y_axis, 0, angles_rotated_relative);
			return;

		case GLFW_KEY_B:
			// Update relative rotation angle for rotating the face on the Z-axis
			angles_rotated_relative[make_tuple('z', RUBIKS_CUBE_SIZE - 1)] += SINGLE_ROTATION_ANGLE;
			scn->rotate_face(SINGLE_ROTATION_ANGLE, cube_z_axis, RUBIKS_CUBE_SIZE - 1, angles_rotated_relative);
			return;

		case GLFW_KEY_F:
			// Update relative rotation angle for rotating the face on the Z-axis
			angles_rotated_relative[make_tuple('z', 0)] += SINGLE_ROTATION_ANGLE;
			scn->rotate_face(SINGLE_ROTATION_ANGLE, cube_z_axis, 0, angles_rotated_relative);
			return;

		case GLFW_KEY_SPACE:
			// Invert the rotation direction
			SINGLE_ROTATION_ANGLE *= -1;
			return;

		case GLFW_KEY_Z:
			// Halve the rotation angles
			SINGLE_ROTATION_ANGLE /= 2;
			TOTAL_ROTATION_ANGLE /= 2;
			return;

		case GLFW_KEY_A:
			// Double the rotation angles if they are less than or equal to 90 degrees
			if (SINGLE_ROTATION_ANGLE <= 90)
				SINGLE_ROTATION_ANGLE *= 2;
			if (TOTAL_ROTATION_ANGLE <= 90)
				TOTAL_ROTATION_ANGLE *= 2;
			return;

		case GLFW_KEY_G:
			// Generate new game state
			scn->generate(window);
			return;

		default:
			break;
		}

		// Apply rotation if axis index is valid
		if (axis_index >= 0)
		{
			// Update absolute rotation angles based on the axis and direction
			for (auto const &x : angles_rotated_absolute)
			{
				if (get<0>(x.first) == 'x' && axis_index == 0)
					angles_rotated_absolute[x.first] += rotation_direction * TOTAL_ROTATION_ANGLE;
				else if (get<0>(x.first) == 'y' && axis_index == 1)
					angles_rotated_absolute[x.first] += rotation_direction * TOTAL_ROTATION_ANGLE;
				else if (get<0>(x.first) == 'z' && axis_index == 2)
					angles_rotated_absolute[x.first] += rotation_direction * TOTAL_ROTATION_ANGLE;
			}

			// Perform cube rotation based on the axis and direction
			switch (axis_index)
			{
			case 0:
				scn->rotate_cube(rotation_direction * TOTAL_ROTATION_ANGLE, glm::vec3(1, 0, 0)); // Rotate around X-axis
				break;
			case 1:
				scn->rotate_cube(rotation_direction * TOTAL_ROTATION_ANGLE, glm::vec3(0, 1, 0)); // Rotate around Y-axis
				break;
			case 2:
				scn->rotate_cube(rotation_direction * TOTAL_ROTATION_ANGLE, glm::vec3(0, 0, 1)); // Rotate around Z-axis
				break;
			default:
				break;
			}
		}
	}
}

void Game::update_rotation_angles(float angle_x, float angle_y)
{
	// Convert input angles from degrees to radians
	float angle_x_radians = (angle_x / 360) * 2 * pi;
	float angle_y_radians = (angle_y / 360) * 2 * pi;

	// Calculate the rotation angle around the z-axis using Euler angles
	float angle_z_radians = atan2(sin(angle_x_radians) * cos(angle_y_radians), cos(angle_x_radians));
	float angle_z = (angle_z_radians / (2 * pi)) * 360;

	// Update the absolute rotation angles for x, y, z based on input
	for (auto &entry : angles_rotated_absolute)
	{
		char axis = get<0>(entry.first);
		if (axis == 'x')
			entry.second += angle_x;
		else if (axis == 'y')
			entry.second += angle_y;
		else if (axis == 'z')
			entry.second += angle_z;
	}

	// Process rotation updates and potential axis switching
	process_axis_rotation('x', angle_x);
	process_axis_rotation('y', angle_y);
	process_axis_rotation('z', angle_z);
}

void Game::process_axis_rotation(char axis, float angle)
{
	bool needs_switch = true;
	for (int i = 0; i < RUBIKS_CUBE_SIZE; i++)
	{
		Face face = make_tuple(axis, i);
		float sign = angle != 0 ? (angle / abs(angle)) : 0;

		bool condition = sign != 0 && (angles_rotated_absolute[face] >= 45 || angles_rotated_absolute[face] <= -45);
		needs_switch &= condition;

		if (condition)
		{
			std::vector<std::vector<std::vector<Shape *>>> new_rubicks_cube = rubicks_cube;
			rotate_data_structure(axis, i, new_rubicks_cube, angle);
			rubicks_cube = new_rubicks_cube;
			angles_rotated_absolute[face] -= sign * 90;
		}
	}

	if (needs_switch)
	{
		switch_cube_axis(axis, angle);
	}
}

void Game::handle_mouse_input(GLFWwindow *window, int button, int action, int mods)
{
	// Only proceed if the mouse button was pressed
	if (action == GLFW_PRESS)
	{
		// Retrieve the Game instance associated with the GLFW window
		Game *game_instance = static_cast<Game *>(glfwGetWindowUserPointer(window));

		// Obtain the current cursor position
		double cursor_x, cursor_y;
		glfwGetCursorPos(window, &cursor_x, &cursor_y);

		// Print the cursor position for demonstration
		std::cout << "Mouse clicked at position: (" << cursor_x << ", " << cursor_y << ")" << std::endl;

		// Perform a simulated action using the cursor position
		game_instance->Picking(static_cast<int>(cursor_x), static_cast<int>(cursor_y));
	}
}

void Game::handle_cursor_movement(GLFWwindow *window, double xpos, double ypos)
{
	// Retrieve the Game instance from the GLFW window's user pointer
	Game *game_instance = static_cast<Game *>(glfwGetWindowUserPointer(window));

	// Update the cursor position with the current cursor coordinates
	game_instance->UpdatePosition(static_cast<float>(xpos), static_cast<float>(ypos));

	// Check the state of the mouse buttons to determine the type of processing required
	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
	{
		// Process right mouse button action (typically used for moving objects)
		game_instance->MouseProccessing(GLFW_MOUSE_BUTTON_RIGHT);
	}
	else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
	{
		// Process left mouse button action (typically used for rotating objects)
		game_instance->MouseProccessing(GLFW_MOUSE_BUTTON_LEFT);
	}
}

void Game::handle_scroll(GLFWwindow *window, double xoffset, double yoffset)
{
	// Retrieve the Game instance from the GLFW window's user pointer
	Game *scn = static_cast<Game *>(glfwGetWindowUserPointer(window));

	const float zTranslate = 1.0f; // Example definition, adjust based on actual usage

	// Move the camera based on vertical scroll input (yoffset),
	// adjusting the camera's position or zoom level.
	scn->MoveCamera(0, zTranslate * yoffset, 0);
}

void Game::rotate_cube_axis(glm::vec3 axis, float angle)
{
	// Create a rotation matrix using the provided axis and angle
	glm::mat4 rotation_matrix = glm::rotate(glm::mat4(1.0f), glm::radians(angle), axis);

	// Apply the rotation matrix to each of the cube's primary axis
	// glm::vec4 is used to convert the vec3 into vec4 with a homogeneous coordinate for proper matrix multiplication
	cube_x_axis = glm::vec3(rotation_matrix * glm::vec4(cube_x_axis, 0.0f));
	cube_y_axis = glm::vec3(rotation_matrix * glm::vec4(cube_y_axis, 0.0f));
	cube_z_axis = glm::vec3(rotation_matrix * glm::vec4(cube_z_axis, 0.0f));
}

// Helper function to update angles for the switched axes
void update_angles(map<Face, float> &angles_rotated_absolute, char axis, float sign)
{
	for (int i = 0; i < RUBIKS_CUBE_SIZE; i++)
	{
		Face face1 = std::make_tuple(axis, i);
		Face face2 = std::make_tuple(axis == 'x' ? 'y' : (axis == 'y' ? 'z' : 'x'), i);
		float prev_angle = angles_rotated_absolute[face1];
		angles_rotated_absolute[face1] = -sign * angles_rotated_absolute[face2];
		angles_rotated_absolute[face2] = sign * prev_angle;
	}
}

// Helper function to switch cube axes
void switch_axes(glm::vec3 &axis1, glm::vec3 &axis2, float sign)
{
	glm::vec3 tmp = axis1;
	axis1 = -sign * axis2;
	axis2 = sign * tmp;
}

void Game::switch_cube_axis(char axis, float angle)
{
	// Determine the sign of the rotation angle
	float sign = angle / std::abs(angle);

	// Perform axis switching and angle updates based on the specified axis
	if (axis == 'x')
	{
		switch_axes(cube_y_axis, cube_z_axis, sign);
		update_angles(angles_rotated_absolute, 'y', sign);
	}
	else if (axis == 'y')
	{
		switch_axes(cube_z_axis, cube_x_axis, sign);
		update_angles(angles_rotated_absolute, 'z', sign);
	}
	else if (axis == 'z')
	{
		switch_axes(cube_x_axis, cube_y_axis, sign);
		update_angles(angles_rotated_absolute, 'x', sign);
	}
}

void Game::rotate_cube(float angle, glm::vec3 axis)
{
	// Store the original angle for later use
	float prev_angle = angle;

	// Clamp the angle to a maximum of 90 degrees or a minimum of -90 degrees
	angle = std::max(std::min(angle, 90.0f), -90.0f);

	// Determine the rotation axis based on the unit vector 'axis'
	char rotation_axis = 'z'; // Default rotation axis is z
	if (glm::all(glm::equal(axis, glm::vec3(1, 0, 0))))
	{
		rotation_axis = 'x'; // Rotation around x-axis
	}
	else if (glm::all(glm::equal(axis, glm::vec3(0, 1, 0))))
	{
		rotation_axis = 'y'; // Rotation around y-axis
	}

	// Rotate each face of the cube and check if all rotations were successful
	bool rotated_data_structure = true;
	for (int z = 0; z < RUBIKS_CUBE_SIZE; z++)
	{
		rotated_data_structure &= rotate_face(angle, axis, z, angles_rotated_absolute);
	}

	// Apply accumulated total rotation to the cube's axis
	rotate_cube_axis(axis, TOTAL_ROTATION_ANGLE);

	// Switch the cube's faces if a turn of more than 45 degrees has been completed
	if (rotated_data_structure)
	{
		// Switch the cube's faces based on the rotation axis
		glm::vec3 tmp;
		if (rotation_axis == 'x')
		{
			tmp = cube_y_axis;
			cube_y_axis = -angle / std::abs(angle) * cube_z_axis;
			cube_z_axis = angle / std::abs(angle) * tmp;

			for (int i = 0; i < RUBIKS_CUBE_SIZE; i++)
			{
				Face face_y = std::make_tuple('y', i);
				Face face_z = std::make_tuple('z', i);
				float prev_y_angle = angles_rotated_absolute[face_y];
				angles_rotated_absolute[face_y] = -angle / std::abs(angle) * angles_rotated_absolute[face_z];
				angles_rotated_absolute[face_z] = angle / std::abs(angle) * prev_y_angle;
			}
		}
		else if (rotation_axis == 'y')
		{
			tmp = cube_z_axis;
			cube_z_axis = -angle / std::abs(angle) * cube_x_axis;
			cube_x_axis = angle / std::abs(angle) * tmp;

			for (int i = 0; i < RUBIKS_CUBE_SIZE; i++)
			{
				Face face_z = std::make_tuple('z', i);
				Face face_x = std::make_tuple('x', i);
				float prev_z_angle = angles_rotated_absolute[face_z];
				angles_rotated_absolute[face_z] = -angle / std::abs(angle) * angles_rotated_absolute[face_x];
				angles_rotated_absolute[face_x] = angle / std::abs(angle) * prev_z_angle;
			}
		}
		else if (rotation_axis == 'z')
		{
			tmp = cube_x_axis;
			cube_x_axis = -angle / std::abs(angle) * cube_y_axis;
			cube_y_axis = angle / std::abs(angle) * tmp;

			for (int i = 0; i < RUBIKS_CUBE_SIZE; i++)
			{
				Face face_x = std::make_tuple('x', i);
				Face face_y = std::make_tuple('y', i);
				float prev_x_angle = angles_rotated_absolute[face_x];
				angles_rotated_absolute[face_x] = -angle / std::abs(angle) * angles_rotated_absolute[face_y];
				angles_rotated_absolute[face_y] = angle / std::abs(angle) * prev_x_angle;
			}
		}
	}

	// Handle cases where the original angle exceeded 90 degrees
	if (prev_angle > 90)
	{
		rotate_cube(prev_angle - 90, axis); // Adjust for the remainder if angle was above 90
	}
	else if (prev_angle < -90)
	{
		rotate_cube(prev_angle + 90, axis); // Adjust for the remainder if angle was below -90
	}
}

bool Game::setup_and_determine_rotation(float &angle, glm::vec3 &axis, int index, map<Face, float> &angles_rotated, char &rotation_axis)
{
	// Clamp the angle to ensure it doesn't exceed +/- 90 degrees for each operation
	angle = std::max(std::min(angle, 90.0f), -90.0f);

	// Determine the rotation axis based on the provided axis
	if ((angles_rotated == angles_rotated_relative && glm::all(glm::equal(axis, cube_x_axis))) ||
		(angles_rotated == angles_rotated_absolute && glm::all(glm::equal(axis, glm::vec3(1, 0, 0)))))
	{
		rotation_axis = 'x';
	}
	else if ((angles_rotated == angles_rotated_relative && glm::all(glm::equal(axis, cube_y_axis))) ||
			 (angles_rotated == angles_rotated_absolute && glm::all(glm::equal(axis, glm::vec3(0, 1, 0)))))
	{
		rotation_axis = 'y';
	}
	else if ((angles_rotated == angles_rotated_relative && glm::all(glm::equal(axis, cube_z_axis))) ||
			 (angles_rotated == angles_rotated_absolute && glm::all(glm::equal(axis, glm::vec3(0, 0, 1)))))
	{
		rotation_axis = 'z';
	}
	else
	{
		return false; // Axis not recognized
	}

	return true;
}

bool Game::perform_rotation(float angle, glm::vec3 axis, int index, map<Face, float> &angles_rotated, char rotation_axis)
{
	// Rotate the specified face of the Rubik's cube
	for (int j = 0; j < RUBIKS_CUBE_SIZE; j++)
	{
		for (int i = 0; i < RUBIKS_CUBE_SIZE; i++)
		{
			Shape *cube = (rotation_axis == 'x') ? rubicks_cube[index][i][j] : (rotation_axis == 'y') ? rubicks_cube[i][index][j]
																									  : rubicks_cube[i][j][index];

			glm::mat4 previous_trans = cube->get_trans();
			glm::vec3 trans_vec1(-previous_trans[3][0], -previous_trans[3][1], -previous_trans[3][2]);
			cube->MyTranslate(trans_vec1, 0);
			cube->MyRotate(angle, axis, 1);
			glm::mat4 trans_mat2 = glm::rotate(glm::mat4(1.0f), angle, axis) * glm::translate(glm::mat4(1.0f), -trans_vec1);
			glm::vec3 trans_vec2(trans_mat2[3][0], trans_mat2[3][1], trans_mat2[3][2]);
			cube->MyTranslate(trans_vec2, 0);
		}
	}

	// Update rotation angles and cube structure
	Face face = std::make_tuple(rotation_axis, index);
	if (angles_rotated[face] >= 45 || angles_rotated[face] <= -45)
	{
		vector<vector<vector<Shape *>>> new_cube = rubicks_cube;
		rotate_data_structure(rotation_axis, index, new_cube, angle);
		rubicks_cube = new_cube;
		angles_rotated[face] -= (angle / std::abs(angle)) * 90; // Adjust the angle post-rotation
	}

	return true;
}

bool Game::rotate_face(float angle, glm::vec3 axis, int index, map<Face, float> &angles_rotated)
{
	char rotation_axis;
	if (!setup_and_determine_rotation(angle, axis, index, angles_rotated, rotation_axis))
	{
		return false; // Early exit if axis is invalid
	}

	if (!perform_rotation(angle, axis, index, angles_rotated, rotation_axis))
	{
		return false; // Handle potential issues during rotation
	}

	// Recursive adjustment for angles exceeding 90 degrees in either direction
	float prev_angle = angle; // Store original angle for recursion
	if (prev_angle > 90)
	{
		return rotate_face(prev_angle - 90, axis, index, angles_rotated);
	}
	else if (prev_angle < -90)
	{
		return rotate_face(prev_angle + 90, axis, index, angles_rotated);
	}

	return true;
}

void Game::rotate_piece(char axis, int layerIndex, int row, int col, float angle, vector<vector<vector<Shape *>>> &newCube)
{
	Shape *cubePiece = nullptr;
	int targetRow, targetCol;

	if (axis == 'x')
	{
		cubePiece = rubicks_cube[layerIndex][row][col];
		targetRow = angle < 0 ? RUBIKS_CUBE_SIZE - 1 - col : col;
		targetCol = angle < 0 ? row : RUBIKS_CUBE_SIZE - 1 - row;
		newCube[layerIndex][targetRow][targetCol] = cubePiece;
	}
	else if (axis == 'y')
	{
		cubePiece = rubicks_cube[row][layerIndex][col];
		targetRow = angle > 0 ? RUBIKS_CUBE_SIZE - 1 - col : col;
		targetCol = angle > 0 ? row : RUBIKS_CUBE_SIZE - 1 - row;
		newCube[targetRow][layerIndex][targetCol] = cubePiece;
	}
	else if (axis == 'z')
	{
		cubePiece = rubicks_cube[row][col][layerIndex];
		targetRow = angle > 0 ? RUBIKS_CUBE_SIZE - 1 - col : col;
		targetCol = angle > 0 ? row : RUBIKS_CUBE_SIZE - 1 - row;
		newCube[targetRow][targetCol][layerIndex] = cubePiece;
	}
}

void Game::rotate_data_structure(char axis, int layerIndex, vector<vector<vector<Shape *>>> &newCube, float angle)
{
	// Iterate over each position in the specified cube layer to apply rotation
	for (int row = 0; row < RUBIKS_CUBE_SIZE; row++)
	{
		for (int col = 0; col < RUBIKS_CUBE_SIZE; col++)
		{
			// Delegate the rotation of each cube piece to the rotate_piece function
			rotate_piece(axis, layerIndex, row, col, angle, newCube);
		}
	}
}

void Game::generate(GLFWwindow *window)
{
	// Define the keys corresponding to cube rotations
	std::vector<int> keys{GLFW_KEY_R, GLFW_KEY_L, GLFW_KEY_U, GLFW_KEY_D, GLFW_KEY_B, GLFW_KEY_F};

	// Seed the random number generator
	std::srand(static_cast<unsigned int>(time(0)));

	// Generate a random number of actions between 10 and 20
	int num_of_actions = std::rand() % 10 + 10;

	// Store the current rotation angle to restore later
	float prev_rotation_angle = SINGLE_ROTATION_ANGLE;
	SINGLE_ROTATION_ANGLE = 45; // Set a fixed rotation angle for the generation sequence

	for (int i = 0; i < num_of_actions; i++)
	{
		// Choose a random key from the set
		int key_index = std::rand() % keys.size();
		// Determine the direction of rotation (either 1 or -1)
		float sign = pow(-1, std::rand() % 2);

		// Apply the direction to the rotation angle
		float rotation_angle = SINGLE_ROTATION_ANGLE * sign;
		// Simulate key press with the random rotation angle
		process_key_input(window, keys[key_index], 0, GLFW_PRESS, 0);
		// Draw and swap buffers with a delay to visually represent the rotation
		Draw(1, 0, BACK, true, false);
		glfwSwapBuffers(window);
		std::this_thread::sleep_for(std::chrono::milliseconds(500));

		// Optionally repeat or adjust additional actions, like further draws or state updates
		// This line simulates additional input handling that might be needed
		process_key_input(window, keys[key_index], 0, GLFW_PRESS, 0);
		Draw(1, 0, BACK, true, false);
		glfwSwapBuffers(window);
		std::this_thread::sleep_for(std::chrono::milliseconds(500));
	}

	// Restore the original rotation angle after the sequence is complete
	SINGLE_ROTATION_ANGLE = prev_rotation_angle;
}

Game::~Game(void)
{
}
