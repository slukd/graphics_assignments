#include <tuple>
#include <string>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/mat3x3.hpp>

class MyShape
{
public:
	glm::vec4 coordinates;
	glm::vec4 color;
	std::string o_r_t;
	float shininess;

	MyShape(std::string type, glm::vec4 coordinates);
	void set_color_and_shininess(glm::vec4 color);
	glm::vec3 get_normal(glm::vec3 intersection_point);
	~MyShape(void);
};

