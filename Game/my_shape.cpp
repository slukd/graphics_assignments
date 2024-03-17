#include "my_shape.h"

MyShape::MyShape(std::string type, glm::vec4 coordinates)
{
	this->coordinates = coordinates;
	this->o_r_t = type;
}

void MyShape::set_color_and_shininess(glm::vec4 color)
{
	glm::vec3 rgb = glm::vec3(color);
	this->color = glm::vec4(rgb, 0.f);
	this->shininess = color[3];
}

glm::vec3 MyShape::get_normal(glm::vec3 intersection_point)
{
	if (this->coordinates[3] < 0) {
		//plane
		return glm::normalize(glm::vec3(this->coordinates));
	} else {
		return glm::normalize(intersection_point - glm::vec3(this->coordinates));
	}
}


MyShape::~MyShape(void)
{
}
