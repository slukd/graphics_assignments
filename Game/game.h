#pragma once
#include "scene.h"
#include "display.h"
#include <map>

using namespace std;

typedef tuple<char, int> Face;
typedef tuple<vector<vector<vector<Shape*>>>, int, float> CubeState;

class Game : public Scene
{
public:
	
	Game();
	Game(float angle,float relationWH,float near, float far);
	void Init();
	void Update(const glm::mat4 &MVP,const glm::mat4 &Model,const int  shaderIndx);
	void ControlPointUpdate();
	void WhenRotate(float angle_x, float angle_y);
	void WhenTranslate();
	void Motion();
	static void my_key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
    static void my_mouse_callback(GLFWwindow *window, int button, int action, int mods);
    static void my_cursor_position_callback(GLFWwindow *window, double xpos, double ypos);
    static void my_scroll_callback(GLFWwindow *window, double xoffset, double yoffset);
    void generate(GLFWwindow* window);
    void solve(GLFWwindow *window);
	int recursive_solve(GLFWwindow* window, map<vector<vector<vector<Shape*>>>, tuple<int, tuple<char, int, float, int>>>& memo, vector<tuple<char, int, float, int>>& steps);
	void rotate_cube_axes(glm::vec3 axis, float angle);
    void switch_cube_axes(char axis, float angle);
    void rotate_cube(float angle, glm::vec3 axis);
    vector<vector<vector<Shape*>>> make_cube();
    bool rotate_face(float angle, glm::vec3 axis, int index, map<Face, float>& angles_rotated);
	void rotate_data_structure(char axis, int index, vector<vector<vector<Shape*>>>& new_rubicks_cube, float angle);
	~Game(void);

private:
	vector<vector<vector<Shape*>>> rubicks_cube;
};

