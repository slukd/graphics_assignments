#include <tuple>
#include <string>
#include <glm/gtc/matrix_transform.hpp>

class Light
{
public:
	glm::vec3 direction;
	glm::vec4 intensity;
	glm::vec3 location; //if it's directional, location is (INFINITY,INFINITY,INFINITY)
	float cos_of_angle; //if it's directional, angle is INFINITY


	Light(glm::vec4 direction);
	void set_intensity(glm::vec4 intensity);
	void set_location(glm::vec4 location);
	~Light(void);
};

