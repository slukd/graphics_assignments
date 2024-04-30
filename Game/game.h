#pragma once
#include "scene.h"
#include "display.h"
#include <map>

using namespace std;

typedef tuple<char, int> Face;
typedef tuple<vector<vector<vector<Shape *>>>, int, float> CubeState;

class Game : public Scene
{
public:
	Game();
	Game(float angle, float relationWH, float near, float far);
	void Init();
	void Update(const glm::mat4 &MVP, const glm::mat4 &Model, const int shaderIndx);
	void ControlPointUpdate();
	void update_rotation_angles(float angle_x, float angle_y);
	void process_axis_rotation(char axis, float angle);
	void WhenTranslate();
	void Motion();
	static void process_key_input(GLFWwindow *window, int key, int scancode, int action, int mods);
	static void handle_mouse_input(GLFWwindow *window, int button, int action, int mods);
	static void handle_cursor_movement(GLFWwindow *window, double xpos, double ypos);
	static void handle_scroll(GLFWwindow *window, double xoffset, double yoffset);
	void rotate_cube_axis(glm::vec3 axis, float angle);
	void switch_cube_axis(char axis, float angle);
	void update_face_rotation(std::map<Face, float> &angles_rotated_absolute, char axis1, char axis2, int size, float sign);
	float clamp_angle(float angle);
	void rotate_cube(float angle, glm::vec3 axis);
	vector<vector<vector<Shape *>>> make_cube();
	bool rotate_face(float angle, glm::vec3 axis, int index, map<Face, float> &angles_rotated);
	bool perform_rotation(float angle, glm::vec3 axis, int index, map<Face, float> &angles_rotated, char rotation_axis);
	bool setup_and_determine_rotation(float &angle, glm::vec3 &axis, int index, map<Face, float> &angles_rotated, char &rotation_axis);
	void rotate_data_structure(char axis, int index, vector<vector<vector<Shape *>>> &new_rubicks_cube, float angle);
	void rotate_piece(char axis, int index, int i, int j, float angle, vector<vector<vector<Shape *>>> &new_rubicks_cube);
	void simulate_key_rotation(GLFWwindow *window, int key, float angle);
	void generate(GLFWwindow *window);
	~Game(void);

private:
	vector<vector<vector<Shape *>>> rubicks_cube;
};
