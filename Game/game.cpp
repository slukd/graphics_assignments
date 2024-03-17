#include "game.h"
#include <iostream>
#include <fstream>
#include <glm/gtc/matrix_transform.hpp>
#include "my_shape.h"
#include "light.h"


//constant parameters:
glm::vec4 SPECULARVALUE(0.7, 0.7, 0.7, 0.f);
int MAXDEPTH = 5; 
float RESOLUTION = 256.f;


//global variables:
std::vector<MyShape> my_shapes; 
std::vector<Light> lights; 
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
	std::cout<<" matrix:"<<std::endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			std::cout<< mat[j][i]<<" ";
		std::cout<<std::endl;
	}
}

Game::Game() : Scene()
{
}

Game::Game(float angle ,float relationWH, float near1, float far1) : Scene(angle,relationWH,near1,far1)
{ 	
}

void Game::Init()
{		

	AddShader("../res/shaders/pickingShader");	
	AddShader("../res/shaders/basicShader");
	
	WIDTH = RESOLUTION;
	HEIGHT = RESOLUTION;
	PIXELHEIGHT = 2.f / HEIGHT;
	PIXELWIDTH = 2.f / WIDTH;
	NUMOFCOLORS = 4;
	DATASIZE = WIDTH * HEIGHT * NUMOFCOLORS;
	unsigned char* data = new unsigned char[DATASIZE];

	for (int i = 0; i < DATASIZE; i++)
		data[0] = 0;

	std::string path = "../res/scenes/scene4.txt";
	ray_tracing(path, data);
	AddTexture(WIDTH, HEIGHT, data);

	AddShape(Plane,-1,TRIANGLES);
	
	pickedShape = 0;
	
	SetShapeTex(0,0);
	MoveCamera(0,zTranslate,10);
	pickedShape = -1;
	
	//ReadPixel(); //uncomment when you are reading from the z-buffer
}

void Game::Update(const glm::mat4 &MVP,const glm::mat4 &Model,const int  shaderIndx)
{
	Shader *s = shaders[shaderIndx];
	int r = ((pickedShape+1) & 0x000000FF) >>  0;
	int g = ((pickedShape+1) & 0x0000FF00) >>  8;
	int b = ((pickedShape+1) & 0x00FF0000) >> 16;
	s->Bind();
	s->SetUniformMat4f("MVP", MVP);
	s->SetUniformMat4f("Normal",Model);
	s->SetUniform4f("lightDirection", 0.0f , 0.0f, -1.0f, 0.0f);
	if(shaderIndx == 0)
		s->SetUniform4f("lightColor",r/255.0f, g/255.0f, b/255.0f,1.0f);
	else 
		s->SetUniform4f("lightColor",0.7f,0.8f,0.1f,1.0f);
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
	if(isActive)
	{
	}
}

Game::~Game(void)
{
}

void Game::ray_tracing(std::string& scene_path, unsigned char* data)
{
	parse_scene(scene_path);

	glm::vec4 color = glm::vec4(0.f, 0.f, 0.f, 0.f);
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			
			std::vector<glm::vec3> pixel_coordinates = anti_aliasing(i, j);
			for (int k = 0; k < pixel_coordinates.size(); k++) {
				glm::vec3 ray_direction = glm::normalize(pixel_coordinates[k] - camera);
				color += send_ray(camera, ray_direction, -1, 0);
			}
			color /= ((float)pixel_coordinates.size() + 1.f);

			for (int k = 0; k < 4; k++) {
				data[i*(int)WIDTH*4 + j*4 + k] = color[k];
			}
		}
	}
}

std::vector<glm::vec3> Game::anti_aliasing(int i, int j)
{
	float top_left_corner_of_pixel_y = (((float)i / HEIGHT) * 2.f - 1.f) * -1.f;
	float top_left_corner_of_pixel_x = ((float)j / (WIDTH)) * 2.f - 1.f;

	std::vector<glm::vec3> pixel_coordinates;

	float pixel_y = top_left_corner_of_pixel_y - PIXELHEIGHT/3.f;
	float pixel_x = top_left_corner_of_pixel_x + PIXELWIDTH/3.f;
	glm::vec3 coordinates(pixel_x, pixel_y, 0.f);
	pixel_coordinates.push_back(coordinates);

	pixel_y = top_left_corner_of_pixel_y - 2.f*PIXELHEIGHT/3.f;
	pixel_x = top_left_corner_of_pixel_x + PIXELWIDTH/3.f;
	coordinates = glm::vec3(pixel_x, pixel_y, 0.f);
	pixel_coordinates.push_back(coordinates);

	pixel_y = top_left_corner_of_pixel_y - PIXELHEIGHT/3.f;
	pixel_x = top_left_corner_of_pixel_x + 2.f*PIXELWIDTH/3.f;
	coordinates = glm::vec3(pixel_x, pixel_y, 0.f);
	pixel_coordinates.push_back(coordinates);

	pixel_y = top_left_corner_of_pixel_y - 2.f*PIXELHEIGHT/3.f;
	pixel_x = top_left_corner_of_pixel_x + 2.f*PIXELWIDTH/3.f;
	coordinates = glm::vec3(pixel_x, pixel_y, 0.f);
	pixel_coordinates.push_back(coordinates);
	
	return pixel_coordinates;
}

glm::vec3 Game::get_pixel_coordinates(int i, int j)
{
	float top_left_corner_of_pixel_y = (((float)i / HEIGHT) * 2.f - 1.f) * -1.f;
	float top_left_corner_of_pixel_x = ((float)j / (WIDTH)) * 2.f - 1.f;

	float center_of_pixel_y = top_left_corner_of_pixel_y - PIXELHEIGHT/2.f;
	float center_of_pixel_x = top_left_corner_of_pixel_x + PIXELWIDTH/2.f;

	glm::vec3 pixel_coordinates(center_of_pixel_x, center_of_pixel_y, 0.f);
	return pixel_coordinates;
}

glm::vec4 Game::send_ray(glm::vec3 origin, glm::vec3 direction, int previous_intersecting_shape_index, int num_of_call)
{	
	if (num_of_call == MAXDEPTH)
		return glm::vec4(0.f, 0.f, 0.f, 0.f);
	
	glm::vec3 intersection_point(-INFINITY, -INFINITY, -INFINITY);
	int intersecting_shape_index = -1;
	float dist_to_intersection = INFINITY;
	for (int i = 0; i < my_shapes.size(); i++) {
		if (i != previous_intersecting_shape_index) {
			glm::vec3 new_intersection_point = check_shape_intersection(i, origin, direction, 0)[0];
			float new_dist = glm::length(new_intersection_point - origin);
			if (new_dist < dist_to_intersection) {
				intersection_point = new_intersection_point;
				intersecting_shape_index = i;
				dist_to_intersection = new_dist;
			}
		}
	}

	glm::vec4 color(0.f, 0.f, 0.f, 0.f);
	if (intersecting_shape_index != -1 && num_of_call < MAXDEPTH - 1) {
		MyShape shape = my_shapes[intersecting_shape_index];
		if (shape.o_r_t == "r") {
			glm::vec3 N = shape.get_normal(intersection_point);

			if (shape.coordinates[3] < 0 && glm::dot(N, direction) < 0)
				N = N * -1.f;
			glm::vec3 R = glm::reflect(direction, N);
			return color += send_ray(intersection_point, R, intersecting_shape_index, num_of_call + 1);
		} else if (shape.o_r_t == "t") {
			glm::vec3 N = shape.get_normal(intersection_point);
			if (shape.coordinates[3] < 0 && glm::dot(N, direction) < 0)
				N = N * -1.f;

			if (shape.coordinates[3] < 0) {
				//plane
				glm::vec3 refracted_direction = glm::normalize(glm::refract(direction, N, 1.f/1.5f));
				return color += send_ray(intersection_point, refracted_direction, intersecting_shape_index, num_of_call + 1);
			} else {
				//sphere
				glm::vec3 refracted_direction = glm::normalize(glm::refract(direction, N, 1.f/1.5f));
				glm::vec3 second_intersection_point = check_shape_intersection(intersecting_shape_index, intersection_point, refracted_direction, num_of_call + 1)[1];
				N = shape.get_normal(second_intersection_point) * -1.f;
				refracted_direction = glm::normalize(glm::refract(refracted_direction, N, 1.5f));
				if (num_of_call == 0) {
					for (int i = 0; i < lights.size(); i++) {
						if (check_light_intersection(i, intersecting_shape_index, intersection_point)) {
							color += specular(origin, intersection_point, intersecting_shape_index, i);
						}
					}
					color += shape.color * ambient_light;
				}

				return color += send_ray(second_intersection_point, refracted_direction, intersecting_shape_index, num_of_call + 1);
			}
		}
	}
	
	if (intersecting_shape_index != -1) {
		MyShape shape = my_shapes[intersecting_shape_index];

		for (int i = 0; i < lights.size(); i++) {
			if (check_light_intersection(i, intersecting_shape_index, intersection_point)) {
				color += diffuse(origin, intersection_point, intersecting_shape_index, i);
				color += specular(origin, intersection_point, intersecting_shape_index, i);
			}
		}

		color /= (float)lights.size();
		color += shape.color * ambient_light;
	}

	return color*255.f;
}

std::vector<glm::vec3> Game::check_shape_intersection(int shape_index, glm::vec3 origin, glm::vec3 direction, int num_of_call)
{
	MyShape shape = my_shapes[shape_index];
	glm::vec3 intersection_point1(-INFINITY, -INFINITY, -INFINITY);
	glm::vec3 intersection_point2(-INFINITY, -INFINITY, -INFINITY);
	std::vector<glm::vec3> intersection_points;
	glm::vec3 N;

	if (shape.coordinates[3] > 0) {
		//shape is a sphere:
		glm::vec3 O = glm::vec3(shape.coordinates);
		glm::vec3 L = (O - origin);
		float r = shape.coordinates[3];
		float tm = glm::dot(L, direction);
		float d2 = pow(glm::length(L), 2) - pow(tm, 2);

		if (d2 <= pow(r, 2)) {

			float th = sqrt(pow(r, 2) - d2);
			float t1 = tm - th;
			float t2 = tm + th;

			intersection_point1 = origin + t1 * direction;
			intersection_point2 = origin + t2 * direction;

			float dist1 = glm::length(origin - intersection_point1);
			float dist2 = glm::length(origin - intersection_point2);

			if (dist1 > dist2) {
				glm::vec3 tmp = intersection_point2;
				intersection_point2 = intersection_point1;
				intersection_point1 = tmp;
			}

			N = glm::normalize(intersection_point1 - O);			
		}
	} else {
		//shape is a plane:
		N = glm::normalize(glm::vec3(shape.coordinates));
		if (glm::dot(N, direction) < 0)
			N = N * -1.f;
		float N_dot_r_d = glm::dot(N, direction);
		if (N_dot_r_d > 0) {
			float v0 = -((glm::dot(N, origin)) + shape.coordinates[3]);
			float t = v0 / N_dot_r_d;
			if (t >= 0)
				intersection_point1 = origin + t * direction;
		}
	}
	intersection_points.push_back(intersection_point1);
	intersection_points.push_back(intersection_point2);
	
	return intersection_points;
}

bool Game::check_light_intersection(int light_index, int intersecting_shape_index, glm::vec3 intersection_point)
{
	Light light = lights[light_index];

	bool is_directional_or_intersection_point_is_in_the_spotlight = false;
	glm::vec3 direction_from_light = glm::normalize(light.direction);
	if (light.cos_of_angle != INFINITY) {
		direction_from_light = glm::normalize(intersection_point - light.location);
		if (glm::dot(direction_from_light, glm::normalize(light.direction)) > light.cos_of_angle) {
			is_directional_or_intersection_point_is_in_the_spotlight = true;
		}
	} else {
		is_directional_or_intersection_point_is_in_the_spotlight = true;
	}

	if (is_directional_or_intersection_point_is_in_the_spotlight) {
		for (int i = 0; i < my_shapes.size(); i++) {
			if (i != intersecting_shape_index && my_shapes[i].coordinates[3] > 0) {
				// fursther intersection point
				std::vector<glm::vec3> light_intersection_points = check_shape_intersection(i, intersection_point, -direction_from_light, 0);
				if (light_intersection_points[0][0] != -INFINITY) {
					glm::vec3 V = glm::normalize(intersection_point - light_intersection_points[0]);
					if (glm::dot(V, direction_from_light) > 0)
						return false;
				}
				if (light_intersection_points[1][0] != -INFINITY) {
					glm::vec3 V = glm::normalize(intersection_point - light_intersection_points[1]);
					if (glm::dot(V, direction_from_light) > 0)
						return false;
				}

			}
		}
		return true;	
	}

	return false;
}

glm::vec4 Game::diffuse(glm::vec3 origin, glm::vec3 intersection_point, int shape_index, int light_index)
{
	glm::vec4 diffuse_color(0.f, 0.f, 0.f, 0.f);
	MyShape shape = my_shapes[shape_index];
	glm::vec3 N = glm::normalize(glm::vec3(shape.coordinates));
	if (shape.coordinates[3] > 0) {
		//shape is a sphere
		glm::vec3 O = glm::vec3(shape.coordinates);
		N = glm::normalize(intersection_point - O);
	}
	else {
		if (glm::dot(N, (intersection_point - origin)) > 0)
			N = N * -1.f;
	}

	Light light = lights[light_index];
	glm::vec3 Li = glm::normalize(-light.direction);
	if (light.cos_of_angle != INFINITY) {
		//light is a spotlight
		Li = glm::normalize(light.location - intersection_point);
	}
	glm::vec4 shape_color = shape.color;
	if (shape.coordinates[3] < 0) {
		bool cond = (int)(1.5 * std::abs(intersection_point[0])) % 2 == (int)(1.5 * std::abs(intersection_point[1])) % 2;
		if (intersection_point[0] < 0)
			cond = !cond;

		if (intersection_point[1] < 0)
			cond = !cond;

		if (cond)
			shape_color *= 0.5;
	}
	diffuse_color += shape_color * std::max(0.f, glm::dot(N, Li)) * light.intensity;

	return diffuse_color;
}

glm::vec4 Game::specular(glm::vec3 origin, glm::vec3 intersection_point, int shape_index, int light_index)
{
	if (glm::length(origin - intersection_point) == 0) {
		return glm::vec4(1.f, 1.f, 1.f, 1.f);
	}
	
	MyShape shape = my_shapes[shape_index];
	glm::vec3 N = glm::normalize(glm::vec3(shape.coordinates));
	glm::vec4 specular_color(0.f, 0.f, 0.f, 0.f);
	if (shape.coordinates[3] > 0) {
		//shape is a sphere
		glm::vec3 O = glm::vec3(shape.coordinates);
		N = glm::normalize(intersection_point - O);
	}
	else {
		if (glm::dot(N, (intersection_point - origin)) > 0)
			N = N * -1.f;
	}

	Light light = lights[light_index];
	glm::vec3 direction_from_light = light.direction;
	if (light.cos_of_angle != INFINITY) {
		//light is a spotlight
		direction_from_light = glm::normalize(intersection_point - light.location);
	}

	glm::vec3 V = glm::normalize(origin - intersection_point);
	glm::vec3 Ri = glm::normalize(glm::reflect(direction_from_light, N));
	
	specular_color += SPECULARVALUE * pow(std::max(0.f, glm::dot(V, Ri)), shape.shininess) * light.intensity;

	return specular_color;
}

int Game::split(const std::string& txt, std::vector<std::string>& strs, char delimeter) {
	size_t pos = txt.find(delimeter);
	size_t initialPos = 0;
	strs.clear();

	// Decompose statement
	while (pos != std::string::npos) {
		strs.push_back(txt.substr(initialPos, pos - initialPos));
		initialPos = pos + 1;

		pos = txt.find(delimeter, initialPos);
	}

	// Add the last one
	strs.push_back(txt.substr(initialPos, std::min(pos, txt.size()) - initialPos + 1));

	return strs.size();
}


void Game::parse_scene(std::string& scene_path)
{
	std::string line;
	std::ifstream myfile;
	myfile.open(scene_path);

	std::vector<std::string> splitted;
	int colorInstruction = 0;
	int intensityInstruction = 0;
	int spotlightIndex = 0;
	std::vector<int> spotlightsIndexes;

	if (myfile.is_open()) {
		while (getline(myfile, line)) {
			split(line, splitted, ' ');

			if (splitted[0] == "e") {
				camera = glm::vec3(std::stof(splitted[1]), std::stof(splitted[2]), std::stof(splitted[3]));
			}

			else if (splitted[0] == "a") {
				ambient_light = glm::vec4(std::stof(splitted[1]), std::stof(splitted[2]),
					std::stof(splitted[3]), std::stof(splitted[4]));
			}

			else if (splitted[0] == "o" || splitted[0] == "t" || splitted[0] == "r") {

				MyShape shape = MyShape(splitted[0], glm::vec4(std::stof(splitted[1]), std::stof(splitted[2]),
					std::stof(splitted[3]), std::stof(splitted[4])));
				my_shapes.push_back(shape);

			}

			else if (splitted[0] == "c") {
				my_shapes[colorInstruction].set_color_and_shininess(glm::vec4(std::stof(splitted[1]), std::stof(splitted[2]),
					std::stof(splitted[3]), std::stof(splitted[4])));
				colorInstruction += 1;

			}

			else if (splitted[0] == "d") {
				Light light = Light(glm::vec4(std::stof(splitted[1]), std::stof(splitted[2]),
					std::stof(splitted[3]), std::stof(splitted[4])));
				if (splitted[4] == "1.0")
					spotlightsIndexes.push_back(spotlightIndex);

				spotlightIndex += 1;
				lights.push_back(light);
			}

			else if (splitted[0] == "p") {
				lights[spotlightsIndexes[0]].set_location(glm::vec4(std::stof(splitted[1]), std::stof(splitted[2]),
					std::stof(splitted[3]), std::stof(splitted[4])));
				spotlightsIndexes.erase(spotlightsIndexes.begin());

			}

			else if (splitted[0] == "i") {
				lights[intensityInstruction].set_intensity(glm::vec4(std::stof(splitted[1]), std::stof(splitted[2]),
					std::stof(splitted[3]), std::stof(splitted[4])));
				intensityInstruction += 1;
			}
		}

		myfile.close();
	}

	else {
		printf("Unable to open file\n");
	}
}
