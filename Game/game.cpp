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
glm::vec3 camera_pos;
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
				glm::vec3 ray_direction = glm::normalize(coord - camera_pos);
				color += send_ray(camera_pos, ray_direction, -1, 0);
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
						if (is_light_reaching_intersection(i, intersecting_shape_index, closest_hit_point))
						{
							color += calculate_specular(ray_origin, closest_hit_point, intersecting_shape_index, i);
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
			if (is_light_reaching_intersection(i, intersecting_shape_index, closest_hit_point))
			{
				color += calculate_diffuse_component(ray_origin, closest_hit_point, intersecting_shape_index, i);
				color += calculate_specular(ray_origin, closest_hit_point, intersecting_shape_index, i);
			}
		}

		// Average the color by the number of lights and add the shape's own color.
		color /= static_cast<float>(scene_lights.size());
		color += shape.color * ambient_light;
	}

	// Return the calculated color, scaled to the range [0, 255].
	return color * 255.f;
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

glm::vec4 Game::calculate_diffuse_component(glm::vec3 observer, glm::vec3 point, int shape_index, int light_index)
{
	glm::vec4 diffuse_color(0.f, 0.f, 0.f, 0.f); // Initialize diffuse color to black.

	MyShape current_shape = scene_shapes[shape_index]; // Retrieve the shape.

	glm::vec3 normal_at_point = glm::normalize(glm::vec3(current_shape.coordinates)); // Initialize the normal at the intersection point.

	if (current_shape.coordinates[3] > 0) // Adjust the normal for spheres.
	{
		glm::vec3 sphere_center = glm::vec3(current_shape.coordinates);
		normal_at_point = glm::normalize(point - sphere_center);
	}
	else // For planes, ensure the normal is oriented correctly.
	{
		if (glm::dot(normal_at_point, (point - observer)) > 0)
			normal_at_point = -normal_at_point;
	}

	Light current_light = scene_lights[light_index]; // Retrieve the light.

	glm::vec3 light_direction = glm::normalize(-current_light.direction); // Calculate the direction of light at the intersection point.

	if (current_light.cos_of_angle != INFINITY) // Adjust light direction for spotlights.
	{
		light_direction = glm::normalize(current_light.location - point);
	}

	glm::vec4 shape_color = current_shape.color; // Retrieve the shape's color.

	// Apply checkerboard pattern for planes.
	if (current_shape.coordinates[3] < 0)
	{
		bool checkered_condition = (int)(1.5 * std::abs(point[0])) % 2 == (int)(1.5 * std::abs(point[1])) % 2;
		if (point[0] < 0)
			checkered_condition = !checkered_condition;
		if (point[1] < 0)
			checkered_condition = !checkered_condition;
		if (checkered_condition)
			shape_color *= 0.5;
	}

	// Calculate the final diffuse color.
	diffuse_color += shape_color * std::max(0.f, glm::dot(normal_at_point, light_direction)) * current_light.intensity;

	return diffuse_color; // Return the calculated diffuse color.
}

bool Game::is_light_reaching_intersection(int light_index, int intersecting_shape_index, glm::vec3 intersection_point)
{
	Light current_light = scene_lights[light_index]; // Retrieve the light.

	bool is_illuminating_intersection = false; // Initialize flag to determine if the light reaches the intersection point.

	glm::vec3 light_direction = glm::normalize(current_light.direction); // Calculate direction from the light.

	if (current_light.cos_of_angle != INFINITY) // Check if the light is a spotlight.
	{
		light_direction = glm::normalize(intersection_point - current_light.location); // Adjust direction for spotlight.
		if (glm::dot(light_direction, glm::normalize(current_light.direction)) > current_light.cos_of_angle)
		{
			is_illuminating_intersection = true; // Intersection point is within spotlight's cone.
		}
	}
	else
	{
		is_illuminating_intersection = true; // For directional lights, intersection point is always illuminated.
	}

	if (is_illuminating_intersection)
	{
		for (int i = 0; i < scene_shapes.size(); i++)
		{
			if (i != intersecting_shape_index && scene_shapes[i].coordinates[3] > 0) // Check for obstructions with other shapes.
			{
				std::vector<glm::vec3> intersections = check_shape_intersection(i, intersection_point, -light_direction, 0);
				if (intersections[0][0] != -INFINITY || intersections[1][0] != -INFINITY) // Check if intersection occurs.
				{
					return false; // Light is obstructed.
				}
			}
		}
		return true; // Light reaches intersection point unobstructed.
	}

	return false; // Default case if light does not reach intersection point.
}

// Calculates the specular component of lighting at an intersection point
glm::vec4 Game::calculate_specular(glm::vec3 origin, glm::vec3 intersection_point, int shape_index, int light_index)
{
	// If the origin and intersection point are the same, return white to avoid division by zero later
	if (glm::length(origin - intersection_point) == 0)
	{
		return glm::vec4(1.f, 1.f, 1.f, 1.f);
	}

	MyShape current_shape = scene_shapes[shape_index]; // Retrieve the shape.

	glm::vec3 normal_at_intersection = glm::normalize(glm::vec3(current_shape.coordinates)); // Calculate normal.

	glm::vec4 specular_color(0.f, 0.f, 0.f, 0.f); // Initialize specular color.

	if (current_shape.coordinates[3] > 0) // Check if shape is a sphere.
	{
		glm::vec3 sphere_center = glm::vec3(current_shape.coordinates);
		normal_at_intersection = glm::normalize(intersection_point - sphere_center); // Calculate normal at intersection.
	}
	else
	{
		if (glm::dot(normal_at_intersection, (intersection_point - origin)) > 0)
			normal_at_intersection = -normal_at_intersection; // Ensure correct orientation for planes.
	}

	Light current_light = scene_lights[light_index]; // Retrieve the light.

	glm::vec3 light_direction = current_light.direction; // Calculate direction from light.

	if (current_light.cos_of_angle != INFINITY) // Check if light is a spotlight.
	{
		light_direction = glm::normalize(intersection_point - current_light.location); // Adjust direction for spotlight.
	}

	glm::vec3 view_vector = glm::normalize(origin - intersection_point);								 // Calculate view vector.
	glm::vec3 reflection_vector = glm::normalize(glm::reflect(light_direction, normal_at_intersection)); // Calculate reflection vector.

	specular_color += REFLECTIVITY * pow(std::max(0.f, glm::dot(view_vector, reflection_vector)), current_shape.shininess) * current_light.intensity; // Add specular component.

	return specular_color; // Return calculated specular color.
}

// Splits a string by a delimiter and returns a vector of the substrings
int Game::custom_split(const std::string &text, std::vector<std::string> &tokens, char delimiter)
{
	// Find the first occurrence of the delimiter
	size_t position = text.find(delimiter);
	// Initialize the starting position for the split
	size_t initial_position = 0;
	// Clear the vector to store the results
	tokens.clear();

	// Loop through the string to find and add substrings
	while (position != std::string::npos)
	{
		tokens.push_back(text.substr(initial_position, position - initial_position));
		initial_position = position + 1;
		position = text.find(delimiter, initial_position);
	}

	// Add the last substring to the vector
	tokens.push_back(text.substr(initial_position));

	// Return the number of substrings found
	return tokens.size();
}

// Parses a scene from a file and initializes the game environment accordingly
void Game::parse_scene(const std::string &scene_path)
{
	// String to store each line read from the file
	std::string line;
	// Opens the scene file
	std::ifstream file;
	file.open(scene_path);

	// Vector to hold the parts of each line after splitting
	std::vector<std::string> parts;
	// Variables to keep track of instructions for setting color and intensity
	int color_instruction = 0;
	int intensity_instruction = 0;
	// Variables for spotlight handling
	int spotlight_index = 0;
	std::vector<int> spotlight_indexes;

	// Checks if the file was successfully opened
	if (file.is_open())
	{
		// Reads the file line by line
		while (std::getline(file, line))
		{
			// Splits each line by spaces
			custom_split(line, parts, ' ');

			// Processes the line based on the first element (command)
			if (parts[0] == "e")
			{
				// Sets the camera position
				camera_pos = glm::vec3(std::stof(parts[1]), std::stof(parts[2]), std::stof(parts[3]));
			}
			else if (parts[0] == "a")
			{
				// Sets the ambient light
				ambient_light = glm::vec4(std::stof(parts[1]), std::stof(parts[2]), std::stof(parts[3]), std::stof(parts[4]));
			}
			else if (parts[0] == "o" || parts[0] == "t" || parts[0] == "r")
			{
				// Adds a new shape to the scene
				MyShape shape = MyShape(parts[0], glm::vec4(std::stof(parts[1]), std::stof(parts[2]), std::stof(parts[3]), std::stof(parts[4])));
				scene_shapes.push_back(shape);
			}
			else if (parts[0] == "c")
			{
				// Sets the color and shininess of the last added shape
				scene_shapes[color_instruction].set_color_and_shininess(glm::vec4(std::stof(parts[1]), std::stof(parts[2]), std::stof(parts[3]), std::stof(parts[4])));
				color_instruction += 1;
			}
			else if (parts[0] == "d")
			{
				// Adds a new light to the scene
				Light light = Light(glm::vec4(std::stof(parts[1]), std::stof(parts[2]), std::stof(parts[3]), std::stof(parts[4])));
				// Checks if the light is a spotlight
				if (parts[4] == "1.0")
					spotlight_indexes.push_back(spotlight_index);
				spotlight_index += 1;
				scene_lights.push_back(light);
			}
			else if (parts[0] == "p")
			{
				// Sets the location of the first spotlight and removes it from the tracking list
				scene_lights[spotlight_indexes[0]].set_location(glm::vec4(std::stof(parts[1]), std::stof(parts[2]), std::stof(parts[3]), std::stof(parts[4])));
				spotlight_indexes.erase(spotlight_indexes.begin());
			}
			else if (parts[0] == "i")
			{
				// Sets the intensity of the light
				scene_lights[intensity_instruction].set_intensity(glm::vec4(std::stof(parts[1]), std::stof(parts[2]), std::stof(parts[3]), std::stof(parts[4])));
				intensity_instruction += 1;
			}
		}
		// Closes the file
		file.close();
	}
	else
	{
		// Prints an error message if the file cannot be opened
		printf("Unable to open file\n");
	}
}
