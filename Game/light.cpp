#include "light.h"

Light::Light(glm::vec4 direction)
{
	this->direction = glm::normalize(glm::vec3(direction));
	this->location = glm::vec3(INFINITY, INFINITY, INFINITY);
	this->cos_of_angle = INFINITY;
}

void Light::set_intensity(glm::vec4 intensity)
{
	this->intensity = intensity;
}

void Light::set_location(glm::vec4 location)
{
	this->location = glm::vec3(location);
	this->cos_of_angle = location[3];
}

Light::~Light(void)
{
}
