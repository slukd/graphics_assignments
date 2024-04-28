#include "game.h"
#include <iostream>
#include <fstream>
#include <glm/gtc/matrix_transform.hpp>
#include "my_shape.h"
#include "light.h"

// constant parameters:
glm::vec4 REFLECTIVITY(0.7, 0.7, 0.7, 0.f);
int MAX_RECURSION_DEPTH = 5;
float IMAGE_RESOLUTION = 256.f;

// global variables:
std::vector<MyShape> scene_shapes;
std::vector<Light> scene_lights;
glm::vec3 camera;
glm::vec4 ambient_light;
int DATASIZE;
float WIDTH;
float HEIGHT;
int NUMOFCOLORS;
float PIXELHEIGHT;
float PIXELWIDTH;

static void printMat(const glm::mat4 mat)
{
	std::cout << " matrix:" << std::endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			std::cout << mat[j][i] << " ";
		std::cout << std::endl;
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

	WIDTH = IMAGE_RESOLUTION;
	HEIGHT = IMAGE_RESOLUTION;
	PIXELHEIGHT = 2.f / HEIGHT;
	PIXELWIDTH = 2.f / WIDTH;
	NUMOFCOLORS = 4;
	DATASIZE = WIDTH * HEIGHT * NUMOFCOLORS;
	unsigned char *data = new unsigned char[DATASIZE];

	for (int i = 0; i < DATASIZE; i++)
		data[0] = 0;

	std::string path = "../res/scenes/scene4.txt";
	ray_tracing(path, data);
	AddTexture(WIDTH, HEIGHT, data);

	AddShape(Plane, -1, TRIANGLES);

	pickedShape = 0;

	SetShapeTex(0, 0);
	MoveCamera(0, zTranslate, 10);
	pickedShape = -1;

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

void Game::WhenRotate()
{
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

Game::~Game(void)
{
}

// Performs ray tracing to generate an image from the scene configuration.
void Game::ray_tracing(std::string &scene_path, unsigned char *data)
{
	// Parse the scene from a file to set up objects and lighting.
	parse_scene(scene_path);

	glm::vec4 color(0.f); // Initialize the color accumulator to black.

	// Iterate over each pixel in the image.
	for (int i = 0; i < HEIGHT; ++i)
	{
		for (int j = 0; j < WIDTH; ++j)
		{

			// Reset color for the current pixel.
			color = glm::vec4(0.f);

			// Get sub-pixel coordinates for anti-aliasing.
			std::vector<glm::vec3> pixel_coordinates = anti_aliasing(i, j);

			// Accumulate color by sending rays through each sub-pixel coordinate.
			for (const auto &coord : pixel_coordinates)
			{
				glm::vec3 ray_direction = glm::normalize(coord - camera);
				color += send_ray(camera, ray_direction, -1, 0);
			}

			// Average the accumulated color by the number of samples.
			color /= static_cast<float>(pixel_coordinates.size());

			for (int k = 0; k < 4; k++)
			{
				data[i * (int)WIDTH * 4 + j * 4 + k] = color[k];
			}
		}
	}
}

// Calculates sub-pixel coordinates for anti-aliasing within a pixel.
// Anti-aliasing is achieved by sampling at several points within each pixel and averaging the results.
std::vector<glm::vec3> Game::anti_aliasing(int i, int j)
{
	// Adjust the calculation for the Y-coordinate to flip the image vertically.
	// This change ensures the top of the image corresponds to the top of the viewport.
	float top_left_corner_of_pixel_y = 1.f - (2.f * i) / HEIGHT;
	float top_left_corner_of_pixel_x = -1.f + (2.f * j) / WIDTH;

	std::vector<glm::vec3> pixel_coordinates; // Stores the sub-pixel coordinates.

	// Sampling points within the pixel to perform anti-aliasing.
	for (int subY = 1; subY <= 2; ++subY)
	{
		for (int subX = 1; subX <= 2; ++subX)
		{
			// Adjust sub-pixel's Y-position calculation to align with the flipped Y-axis.
			float pixel_y = top_left_corner_of_pixel_y - (PIXELHEIGHT * subY / 3.f) + PIXELHEIGHT / 2.f;
			float pixel_x = top_left_corner_of_pixel_x + (PIXELWIDTH * subX / 3.f) - PIXELWIDTH / 2.f;

			// Add the calculated sub-pixel coordinate to the list.
			pixel_coordinates.emplace_back(pixel_x, pixel_y, 0.f);
		}
	}

	return pixel_coordinates;
}

glm::vec4 Game::send_ray(glm::vec3 ray_origin, glm::vec3 ray_direction, int previous_intersecting_shape_index, int num_of_call)
{
	// Terminate if the maximum recursion depth has been reached to avoid infinite loops.
	if (num_of_call == MAX_RECURSION_DEPTH)
		return glm::vec4(0.f, 0.f, 0.f, 0.f);

	// Initialize variables to track the closest intersection point.
	glm::vec3 closest_hit_point(-INFINITY, -INFINITY, -INFINITY);
	int intersecting_shape_index = -1;
	float dist_to_intersection = INFINITY;

	// Loop through all the shapes in the scene to find intersections with the ray.
	for (int i = 0; i < scene_shapes.size(); i++)
	{
		// Skip the shape if it was the last intersected shape to prevent self-intersection.
		if (i != previous_intersecting_shape_index)
		{
			glm::vec3 new_intersection_point = check_shape_intersection(i, ray_origin, ray_direction, 0)[0];
			float new_dist = glm::length(new_intersection_point - ray_origin);
			// Check if this is the closest intersection so far.
			if (new_dist < dist_to_intersection)
			{
				closest_hit_point = new_intersection_point;
				intersecting_shape_index = i;
				dist_to_intersection = new_dist;
			}
		}
	}

	// Initialize color to black, which will accumulate the calculated color.
	glm::vec4 color(0.f, 0.f, 0.f, 0.f);

	// Check if an intersection was found and if we are not at the last recursion depth.
	if (intersecting_shape_index != -1 && num_of_call < MAX_RECURSION_DEPTH - 1)
	{
		MyShape shape = scene_shapes[intersecting_shape_index];
		// Correct the normal if the intersection is with the inside of a shape.
		glm::vec3 N = shape.get_normal(closest_hit_point);
		if (shape.coordinates[3] < 0 && glm::dot(N, ray_direction) < 0)
			N = N * -1.f;

		// Reflective material handling.
		if (shape.o_r_t == "r")
		{
			glm::vec3 R = glm::reflect(ray_direction, N);
			// Recursively send a new ray in the reflection direction.
			return color += send_ray(closest_hit_point, R, intersecting_shape_index, num_of_call + 1);
		}
		// Transparent material handling.
		else if (shape.o_r_t == "t")
		{
			glm::vec3 refracted_direction;
			// Check if the intersection is with a plane or a sphere to handle refraction.
			if (shape.coordinates[3] < 0)
			{
				// Handle refraction for a plane.
				refracted_direction = glm::normalize(glm::refract(ray_direction, N, 1.f / 1.5f));
				return color += send_ray(closest_hit_point, refracted_direction, intersecting_shape_index, num_of_call + 1);
			}
			else
			{
				// Handle refraction for a sphere.
				refracted_direction = glm::normalize(glm::refract(ray_direction, N, 1.f / 1.5f));
				glm::vec3 second_intersection_point = check_shape_intersection(intersecting_shape_index, closest_hit_point, refracted_direction, num_of_call + 1)[1];
				N = shape.get_normal(second_intersection_point) * -1.f;
				refracted_direction = glm::normalize(glm::refract(refracted_direction, N, 1.5f));
				// If it's the first recursion call, add lighting effects.
				if (num_of_call == 0)
				{
					for (int i = 0; i < scene_lights.size(); i++)
					{
						if (check_light_intersection(i, intersecting_shape_index, closest_hit_point))
						{
							color += specular(ray_origin, closest_hit_point, intersecting_shape_index, i);
						}
					}
					color += shape.color * ambient_light; // Add ambient light to the color.
				}
				return color += send_ray(second_intersection_point, refracted_direction, intersecting_shape_index, num_of_call + 1);
			}
		}
	}

	// If an intersection was found, but we're at the last recursion depth, handle local illumination.
	if (intersecting_shape_index != -1)
	{
		MyShape shape = scene_shapes[intersecting_shape_index];

		// Calculate local illumination effects from each light source in the scene.
		for (int i = 0; i < scene_lights.size(); i++)
		{
			if (check_light_intersection(i, intersecting_shape_index, closest_hit_point))
			{
				color += diffuse(ray_origin, closest_hit_point, intersecting_shape_index, i);
				color += specular(ray_origin, closest_hit_point, intersecting_shape_index, i);
			}
		}

		// Average the color by the number of lights and add the shape's own color.
		color /= static_cast<float>(scene_lights.size());
		color += shape.color * ambient_light;
	}

	// Return the calculated color, scaled to the range [0, 255].
	return color * 255.f;
}

// Defines a method to check for intersections between a ray and a shape within the game environment
std::vector<glm::vec3> Game::check_shape_intersection(int shape_index, glm::vec3 origin, glm::vec3 direction, int num_of_call)
{
	// Retrieves the shape from the scene based on its index
	MyShape shape = scene_shapes[shape_index];
	// Initializes intersection points with values that represent no intersection
	glm::vec3 intersection_point1(-INFINITY, -INFINITY, -INFINITY);
	glm::vec3 intersection_point2(-INFINITY, -INFINITY, -INFINITY);
	// Vector to store the intersection points
	std::vector<glm::vec3> intersection_points;
	// Vector to store the normal at the intersection point
	glm::vec3 N;

	// Checks if the shape is a sphere based on the fourth coordinate value
	if (shape.coordinates[3] > 0)
	{
		// Sphere logic
		glm::vec3 O = glm::vec3(shape.coordinates); // Sphere center
		glm::vec3 L = O - origin;					// Vector from ray origin to sphere center
		float r = shape.coordinates[3];				// Sphere radius
		float tm = glm::dot(L, direction);			// Distance along ray to point closest to sphere center
		float d2 = glm::dot(L, L) - tm * tm;		// Correct way to calculate squared distance from sphere center to ray

		if (d2 <= r * r)
		{
			// If d2 is less than or equal to squared radius, there are intersection points
			float th = sqrt(r * r - d2); // Half distance between intersection points
			float t1 = tm - th;			 // Distance to first intersection point
			float t2 = tm + th;			 // Distance to second intersection point

			intersection_point1 = origin + t1 * direction; // First intersection point
			intersection_point2 = origin + t2 * direction; // Second intersection point

			if (glm::length(origin - intersection_point1) > glm::length(origin - intersection_point2))
			{
				std::swap(intersection_point1, intersection_point2); // Ensure intersection_point1 is the closer one
			}

			// No change required here, as this part is correct
			N = glm::normalize(intersection_point1 - O);
		}
	}
	else
	{
		// Calculates intersection with a plane
		N = glm::normalize(glm::vec3(shape.coordinates)); // Plane normal
														  // Ensures the normal is facing towards the ray origin
		if (glm::dot(N, direction) < 0)
			N = -N;
		float N_dot_r_d = glm::dot(N, direction);
		// If N_dot_r_d is positive, there is an intersection
		if (N_dot_r_d > 0)
		{
			float v0 = -(glm::dot(N, origin) + shape.coordinates[3]);
			float t = v0 / N_dot_r_d; // Distance along ray to intersection
									  // If t is positive, calculates the intersection point
			if (t >= 0)
				intersection_point1 = origin + t * direction;
		}
	}
	// Adds the calculated intersection points to the vector
	intersection_points.push_back(intersection_point1);
	intersection_points.push_back(intersection_point2);

	// Returns the intersection points
	return intersection_points;
}

// Checks if the light from a specific source reaches an intersection point without being obstructed
bool Game::check_light_intersection(int light_index, int intersecting_shape_index, glm::vec3 intersection_point)
{
	// Retrieves the light from the scene based on its index
	Light light = scene_lights[light_index];

	// Initializes a flag to determine if the light is directional or if the intersection point is illuminated by the spotlight
	bool is_directional_or_intersection_point_is_in_the_spotlight = false;
	// Calculates the direction from the light to the intersection point
	glm::vec3 direction_from_light = glm::normalize(light.direction);
	// Checks if the light is a spotlight
	if (light.cos_of_angle != INFINITY)
	{
		direction_from_light = glm::normalize(intersection_point - light.location);
		// Determines if the intersection point is within the spotlight's cone
		if (glm::dot(direction_from_light, glm::normalize(light.direction)) > light.cos_of_angle)
		{
			is_directional_or_intersection_point_is_in_the_spotlight = true;
		}
	}
	else
	{
		// For directional lights, this flag is always true
		is_directional_or_intersection_point_is_in_the_spotlight = true;
	}

	// If the intersection point is illuminated by the light
	if (is_directional_or_intersection_point_is_in_the_spotlight)
	{
		// Checks for obstructions between the light and the intersection point
		for (int i = 0; i < scene_shapes.size(); i++)
		{
			if (i != intersecting_shape_index && scene_shapes[i].coordinates[3] > 0)
			{
				// Checks for intersections with other shapes in the scene
				std::vector<glm::vec3> light_intersection_points = check_shape_intersection(i, intersection_point, -direction_from_light, 0);
				// If an intersection is found, determines if it obstructs the light
				if (light_intersection_points[0][0] != -INFINITY || light_intersection_points[1][0] != -INFINITY)
				{
					return false; // Light is obstructed
				}
			}
		}
		return true; // Light reaches the intersection point unobstructed
	}

	return false; // Default case if light does not reach the intersection point
}

// Calculates the diffuse component of lighting at an intersection point
glm::vec4 Game::diffuse(glm::vec3 origin, glm::vec3 intersection_point, int shape_index, int light_index)
{
	// Initializes the diffuse color to black
	glm::vec4 diffuse_color(0.f, 0.f, 0.f, 0.f);
	// Retrieves the shape based on index
	MyShape shape = scene_shapes[shape_index];
	// Initializes the normal at the intersection point
	glm::vec3 N = glm::normalize(glm::vec3(shape.coordinates));
	// Adjusts the normal based on the shape type
	if (shape.coordinates[3] > 0)
	{
		// For spheres, calculates the normal at the intersection point
		glm::vec3 O = glm::vec3(shape.coordinates);
		N = glm::normalize(intersection_point - O);
	}
	else
	{
		// For planes, ensures the normal is oriented correctly
		if (glm::dot(N, (intersection_point - origin)) > 0)
			N = -N;
	}

	// Retrieves the light based on index
	Light light = scene_lights[light_index];
	// Calculates the direction of light at the intersection point
	glm::vec3 Li = glm::normalize(-light.direction);
	// Adjusts the light direction for spotlights
	if (light.cos_of_angle != INFINITY)
	{
		Li = glm::normalize(light.location - intersection_point);
	}
	// Retrieves the shape's color
	glm::vec4 shape_color = shape.color;
	// Modifies the shape color for checkerboard pattern on planes
	if (shape.coordinates[3] < 0)
	{
		bool cond = (int)(1.5 * std::abs(intersection_point[0])) % 2 == (int)(1.5 * std::abs(intersection_point[1])) % 2;
		if (intersection_point[0] < 0)
			cond = !cond;
		if (intersection_point[1] < 0)
			cond = !cond;
		if (cond)
			shape_color *= 0.5;
	}
	// Calculates the final diffuse color
	diffuse_color += shape_color * std::max(0.f, glm::dot(N, Li)) * light.intensity;

	// Returns the calculated diffuse color
	return diffuse_color;
}

// Calculates the specular component of lighting at an intersection point
glm::vec4 Game::specular(glm::vec3 origin, glm::vec3 intersection_point, int shape_index, int light_index)
{
	// If the origin and intersection point are the same, returns white (to avoid division by zero later)
	if (glm::length(origin - intersection_point) == 0)
	{
		return glm::vec4(1.f, 1.f, 1.f, 1.f);
	}

	// Retrieves the shape based on index
	MyShape shape = scene_shapes[shape_index];
	// Initializes the normal at the intersection point
	glm::vec3 N = glm::normalize(glm::vec3(shape.coordinates));
	// Initializes the specular color to black
	glm::vec4 specular_color(0.f, 0.f, 0.f, 0.f);
	// Adjusts the normal based on the shape type
	if (shape.coordinates[3] > 0)
	{
		// For spheres, calculates the normal at the intersection point
		glm::vec3 O = glm::vec3(shape.coordinates);
		N = glm::normalize(intersection_point - O);
	}
	else
	{
		// For planes, ensures the normal is oriented correctly
		if (glm::dot(N, (intersection_point - origin)) > 0)
			N = -N;
	}

	// Retrieves the light based on index
	Light light = scene_lights[light_index];
	// Calculates the direction from light to the intersection point
	glm::vec3 direction_from_light = light.direction;
	// Adjusts the direction for spotlights
	if (light.cos_of_angle != INFINITY)
	{
		direction_from_light = glm::normalize(intersection_point - light.location);
	}

	// Calculates the view vector and the reflection vector
	glm::vec3 V = glm::normalize(origin - intersection_point);
	glm::vec3 Ri = glm::normalize(glm::reflect(direction_from_light, N));

	// Adds the specular component to the color based on the angle between the view and reflection vectors
	specular_color += REFLECTIVITY * pow(std::max(0.f, glm::dot(V, Ri)), shape.shininess) * light.intensity;

	// Returns the calculated specular color
	return specular_color;
}

// Splits a string by a delimiter and returns a vector of the substrings
int Game::split(const std::string &txt, std::vector<std::string> &strs, char delimeter)
{
	// Finds the first occurrence of the delimiter
	size_t pos = txt.find(delimeter);
	// Initializes the starting position for the split
	size_t initialPos = 0;
	// Clears the vector to store the results
	strs.clear();

	// Loops through the string to find and add substrings
	while (pos != std::string::npos)
	{
		strs.push_back(txt.substr(initialPos, pos - initialPos));
		initialPos = pos + 1;
		pos = txt.find(delimeter, initialPos);
	}

	// Adds the last substring to the vector
	strs.push_back(txt.substr(initialPos, std::min(pos, txt.size()) - initialPos + 1));

	// Returns the number of substrings found
	return strs.size();
}

// Parses a scene from a file and initializes the game environment accordingly
void Game::parse_scene(std::string &scene_path)
{
	// String to store each line read from the file
	std::string line;
	// Opens the scene file
	std::ifstream myfile;
	myfile.open(scene_path);

	// Vector to hold the parts of each line after splitting
	std::vector<std::string> splitted;
	// Variables to keep track of instructions for setting color and intensity
	int colorInstruction = 0;
	int intensityInstruction = 0;
	// Variables for spotlight handling
	int spotlightIndex = 0;
	std::vector<int> spotlightsIndexes;

	// Checks if the file was successfully opened
	if (myfile.is_open())
	{
		// Reads the file line by line
		while (getline(myfile, line))
		{
			// Splits each line by spaces
			split(line, splitted, ' ');

			// Processes the line based on the first element (command)
			if (splitted[0] == "e")
			{
				// Sets the camera position
				camera = glm::vec3(std::stof(splitted[1]), std::stof(splitted[2]), std::stof(splitted[3]));
			}
			else if (splitted[0] == "a")
			{
				// Sets the ambient light
				ambient_light = glm::vec4(std::stof(splitted[1]), std::stof(splitted[2]),
										  std::stof(splitted[3]), std::stof(splitted[4]));
			}
			else if (splitted[0] == "o" || splitted[0] == "t" || splitted[0] == "r")
			{
				// Adds a new shape to the scene
				MyShape shape = MyShape(splitted[0], glm::vec4(std::stof(splitted[1]), std::stof(splitted[2]),
															   std::stof(splitted[3]), std::stof(splitted[4])));
				scene_shapes.push_back(shape);
			}
			else if (splitted[0] == "c")
			{
				// Sets the color and shininess of the last added shape
				scene_shapes[colorInstruction].set_color_and_shininess(glm::vec4(std::stof(splitted[1]), std::stof(splitted[2]),
																				 std::stof(splitted[3]), std::stof(splitted[4])));
				colorInstruction += 1;
			}
			else if (splitted[0] == "d")
			{
				// Adds a new light to the scene
				Light light = Light(glm::vec4(std::stof(splitted[1]), std::stof(splitted[2]),
											  std::stof(splitted[3]), std::stof(splitted[4])));
				// Checks if the light is a spotlight
				if (splitted[4] == "1.0")
					spotlightsIndexes.push_back(spotlightIndex);
				spotlightIndex += 1;
				scene_lights.push_back(light);
			}
			else if (splitted[0] == "p")
			{
				// Sets the location of the first spotlight and removes it from the tracking list
				scene_lights[spotlightsIndexes[0]].set_location(glm::vec4(std::stof(splitted[1]), std::stof(splitted[2]),
																		  std::stof(splitted[3]), std::stof(splitted[4])));
				spotlightsIndexes.erase(spotlightsIndexes.begin());
			}
			else if (splitted[0] == "i")
			{
				// Sets the intensity of the light
				scene_lights[intensityInstruction].set_intensity(glm::vec4(std::stof(splitted[1]), std::stof(splitted[2]),
																		   std::stof(splitted[3]), std::stof(splitted[4])));
				intensityInstruction += 1;
			}
		}
		// Closes the file
		myfile.close();
	}
	else
	{
		// Prints an error message if the file cannot be opened
		printf("Unable to open file\n");
	}
}
