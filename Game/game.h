#pragma once
#include "scene.h"
#include <map>

class Game : public Scene
{
public:
	
	Game();
	Game(float angle,float relationWH,float near, float far);
	void Init();
	void Update(const glm::mat4 &MVP,const glm::mat4 &Model,const int  shaderIndx);
	void ControlPointUpdate();
	void WhenRotate();
	void WhenTranslate();
	void Motion();
	~Game(void);

private:
	void ray_tracing(std::string &scene_path, unsigned char *data);
    std::vector<glm::vec3> anti_aliasing(int i, int j);
    glm::vec3 get_pixel_coordinates(int i, int j);
    glm::vec4 send_ray(glm::vec3 origin, glm::vec3 direction, int previous_intersecting_shape_index, int num_of_call);
    void parse_scene(std::string &scene_path);
    std::vector<glm::vec3> check_shape_intersection(int shape_index, glm::vec3 origin, glm::vec3 direction, int num_of_call);
    bool check_light_intersection(int light_index, int intersecting_shape_index, glm::vec3 intersection_point);
	glm::vec4 diffuse(glm::vec3 origin, glm::vec3 intersection_point, int shape_index, int light_index);
	glm::vec4 specular(glm::vec3 origin, glm::vec3 intersection_point, int shape_index, int light_index);
	int split(const std::string& txt, std::vector<std::string>& strs, char delimeter);

};

