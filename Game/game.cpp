

#define _USE_MATH_DEFINES

#include "game.h"
#include <iostream>
#include <fstream>
#include <glm/gtc/matrix_transform.hpp>
#include "stb_image.h"
#include <cmath>
#include <map>

// constants

int WIDTH;
int HEIGHT;
int COMP;
const double SIGMA = 1.175;
int DATA_SIZE;
const int NORMAL = 0;
const int GRAYSCALE = 1;
const int BLACK_AND_WHITE = 2;
const int THRESHOLD_UPPER_BOUND = 40;
const int THRESHOLD_LOWER_BOUND = 18;
int NUM_OF_COLORS = 4; // 4 is (R, G ,B, alpha)

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

	AddShape(Plane, -1, TRIANGLES);

	// Load the image file for texture mapping
	std::string imageFilePath = "../res/textures/lena256.jpg";
	unsigned char *imageData = stbi_load(imageFilePath.c_str(), &WIDTH, &HEIGHT, &COMP, NUM_OF_COLORS);
	if (imageData == NULL)
	{
		std::cerr << "ERROR: Unable to load the image file: " << imageFilePath << std::endl;
		exit(1);
	}
	DATA_SIZE = WIDTH * HEIGHT * NUM_OF_COLORS;
	AddTexture(imageFilePath, false);

	// Apply edge detection to the image
	unsigned char *edgeDetectionData = new unsigned char[DATA_SIZE];
	edge_detection(imageData, edgeDetectionData);
	AddTexture(WIDTH, HEIGHT, edgeDetectionData);
	print_matrix(edgeDetectionData, "../build/assignment/img4.txt", BLACK_AND_WHITE);
	delete[] edgeDetectionData;

	// Apply halftone effect to the image
	unsigned char *halftoneData = new unsigned char[DATA_SIZE * 4];
	halftone(imageData, halftoneData);
	AddTexture(2 * WIDTH, 2 * HEIGHT, halftoneData);
	print_matrix(halftoneData, "../build/assignment/img5.txt", BLACK_AND_WHITE);
	delete[] halftoneData;

	// Apply Floyd-Steinberg dithering to the image
	unsigned char *floydSteinbergData = new unsigned char[DATA_SIZE];
	floyd_steinberg(imageData, floydSteinbergData);
	AddTexture(WIDTH, HEIGHT, floydSteinbergData);
	print_matrix(floydSteinbergData, "../build/assignment/img6.txt", GRAYSCALE);
	delete[] floydSteinbergData;

	MoveCamera(0, zTranslate, 10);
	pickedShape = -1;

	stbi_image_free(imageData);

	// ReadPixel(); //uncomment when you are reading from the z-buffer
}

void Game::edge_detection(unsigned char *data, unsigned char *new_data)
{
	unsigned char *smoothed_data = new unsigned char[DATA_SIZE];
	unsigned char *compute_gradient_data = new unsigned char[DATA_SIZE];
	float *pixel_theta_data = new float[DATA_SIZE];
	unsigned char *non_max_suppression_data = new unsigned char[DATA_SIZE];

	smoothing(data, smoothed_data);

	compute_gradient(smoothed_data, compute_gradient_data, pixel_theta_data);

	non_max_suppression(compute_gradient_data, non_max_suppression_data, pixel_theta_data);

	hysteresis(non_max_suppression_data, new_data);

	delete[] smoothed_data;
	delete[] compute_gradient_data;
	delete[] pixel_theta_data;
	delete[] non_max_suppression_data;
}

void Game::convolution(unsigned char *data, unsigned char *new_data, std::vector<std::vector<float>> &kernel)
{
	const int kernel_size = kernel.size();
	const int middle = (kernel_size - 1) / 2;

	// generate matrix (in size of kernel) to be multiplied by kernel
	std::vector<std::vector<float>> matrix(kernel_size, std::vector<float>(kernel_size, 0));

	for (int data_i = 0; data_i < DATA_SIZE; data_i++)
	{

		int data_row_as_3d_mat = data_i / (NUM_OF_COLORS * WIDTH);
		int data_column_as_3d_mat = (data_i % (NUM_OF_COLORS * WIDTH)) / NUM_OF_COLORS;
		int data_color_as_3d_mat = data_i % NUM_OF_COLORS;

		for (int matrix_i = 0; matrix_i < kernel_size; matrix_i++)
		{

			int row = data_row_as_3d_mat - middle + matrix_i;
			for (int matrix_j = 0; matrix_j < kernel_size; matrix_j++)
			{

				int column = data_column_as_3d_mat - middle + matrix_j;
				if (row >= 0 && column >= 0 && row < HEIGHT && column < WIDTH)
				{
					if (DATA_SIZE <= row * WIDTH * NUM_OF_COLORS + column * NUM_OF_COLORS + data_color_as_3d_mat)
					{
						continue;
					}
					matrix[matrix_i][matrix_j] = data[row * WIDTH * NUM_OF_COLORS + column * NUM_OF_COLORS + data_color_as_3d_mat];
				}
			}
		}

		// multiply matrix[x][y] with kernel[x][y] and sum up to "new pixel"
		float new_pixel = 0;
		for (int matrix_i = 0; matrix_i < kernel_size; matrix_i++)
		{
			for (int matrix_j = 0; matrix_j < kernel_size; matrix_j++)
			{
				new_pixel += matrix[matrix_i][matrix_j] * kernel[matrix_i][matrix_j];
			}
		}

		// update
		new_data[data_i] = (char)(abs((int)new_pixel));
	}
}

// This function performs image smoothing using a Gaussian filter
void Game::smoothing(unsigned char *imageData, unsigned char *smoothedData)
{
	// Define the kernel radius
	const int kernel_r = 2;

	// Calculate the kernel size
	const int kernel_size = 2 * kernel_r + 1;

	// Initialize a 2D vector to store the Gaussian filter kernel
	std::vector<std::vector<float>> gaussian_kernel(kernel_size, std::vector<float>(kernel_size, 0));

	// Generate the Gaussian filter kernel
	for (int i = 0; i < kernel_size; i++)
	{
		for (int j = 0; j < kernel_size; j++)
		{
			// Calculate the Gaussian value at each point in the kernel
			float normalization_factor = 1 / (2.0 * M_PI * pow(SIGMA, 4));
			gaussian_kernel[i][j] = normalization_factor * exp(-((pow(i - kernel_r, 2) + pow(j - kernel_r, 2)) / (2.0 * pow(SIGMA, 2))));
		}
	}

	// Pass a reference to the Gaussian kernel for convolution
	std::vector<std::vector<float>> &gaussian_kernel_ref = gaussian_kernel;

	// Apply convolution with the Gaussian filter
	convolution(imageData, smoothedData, gaussian_kernel_ref);
}

// This function computes the gradient magnitude and direction of an image
void Game::compute_gradient(unsigned char *image_data, unsigned char *gradient_magnitude, float *gradient_direction)
{
	// Allocate memory for storing the derivatives along the x and y directions
	unsigned char *gradient_x = new unsigned char[DATA_SIZE];
	unsigned char *gradient_y = new unsigned char[DATA_SIZE];

	// Sobel filters for computing gradient along x and y directions
	std::vector<std::vector<float>> sobel_x{
		{1, 0, -1},
		{2, 0, -2},
		{1, 0, -1}};
	std::vector<std::vector<float>> sobel_y{
		{1, 2, 1},
		{0, 0, 0},
		{-1, -2, -1}};

	// Create references to the Sobel filters for convolution
	std::vector<std::vector<float>> &sobel_x_ref = sobel_x;
	std::vector<std::vector<float>> &sobel_y_ref = sobel_y;

	// Convolve the input image with the Sobel filters to compute the gradients
	convolution(image_data, gradient_x, sobel_x_ref);
	convolution(image_data, gradient_y, sobel_y_ref);

	// Compute gradient magnitude and direction for each pixel
	for (int i = 0; i < DATA_SIZE; i++)
	{
		// Compute gradient magnitude
		gradient_magnitude[i] = sqrt((pow(gradient_x[i], 2) + pow(gradient_y[i], 2)));

		// Compute gradient direction
		if (gradient_x[i] == 0)
			gradient_direction[i] = M_PI / 2;
		else
			gradient_direction[i] = atan(static_cast<float>(gradient_y[i]) / gradient_x[i]);
	}

	// Free dynamically allocated memory
	delete[] gradient_x;
	delete[] gradient_y;
}

// Function to determine the index of a neighboring pixel
int Game::determine_neighbor_index(int current_row, int current_column, int current_pixel, int row_offset, int column_offset)
{
	// Calculate the new row, column, and pixel values
	int new_row = current_row + row_offset;
	int new_column = current_column + column_offset;
	int new_pixel = current_pixel;

	// Calculate the index of the neighboring pixel in the flattened array
	int neighbor_index = new_row * WIDTH * NUM_OF_COLORS + new_column * NUM_OF_COLORS + new_pixel;

	// Check if the new row and column are within the image boundaries
	if ((0 <= new_row && new_row < HEIGHT) && (0 <= new_column && new_column < WIDTH))
	{
		return neighbor_index; // Return the index of the neighboring pixel
	}

	return -1; // Return -1 if the neighbor is out of bounds
}

// This function performs non-maximum suppression on an image
void Game::non_max_suppression(unsigned char *data, unsigned char *new_data, float *pixel_theta)
{
	int neighbor_index;

	// Define angle steps for different theta ranges
	std::map<int, std::vector<int>> angle_step_map = {
		{8, {-1, 0, 1, 0}},	 // Horizontal gradients
		{6, {-1, 1, 1, -1}}, // 45 and 135 degrees
		{4, {0, -1, 0, 1}},	 // 90 degrees
		{2, {-1, -1, 1, 1}}, // 135 and 45 degrees
		{0, {-1, 0, 1, 0}}	 // 0 degrees
	};

	// Loop through each pixel in the image
	for (int i = 0; i < DATA_SIZE; i++)
	{
		// Calculate 3D matrix coordinates from pixel index
		int row_as_3d_mat = i / (NUM_OF_COLORS * WIDTH);
		int column_as_3d_mat = (i % (NUM_OF_COLORS * WIDTH)) / NUM_OF_COLORS;
		int color_as_3d_mat = i % NUM_OF_COLORS;

		// Convert theta to an integer representing its range
		float current_theta = pixel_theta[i];
		int index_in_angle_step = ((current_theta + (M_PI / 2)) / (M_PI / 8) + 1);
		index_in_angle_step -= (index_in_angle_step % 2);

		// Find angle steps for the current theta range
		std::vector<int> angle_step;
		auto it = angle_step_map.find(index_in_angle_step);
		if (it != angle_step_map.end())
		{
			angle_step = it->second;
		}
		else
		{
			// If angle steps are not found, set default steps
			angle_step = {0, 0};
			continue;
		}

		// Compare edge strength of current pixel with neighboring pixels in positive and negative gradient directions
		neighbor_index = determine_neighbor_index(row_as_3d_mat, column_as_3d_mat, color_as_3d_mat, angle_step[0], angle_step[1]);
		if (neighbor_index != -1 && data[i] < data[neighbor_index])
		{
			// Suppress non-maximum pixels
			new_data[i] = 0;
			continue;
		}
		neighbor_index = determine_neighbor_index(row_as_3d_mat, column_as_3d_mat, color_as_3d_mat, angle_step[2], angle_step[3]);
		if (neighbor_index != -1 && data[i] < data[neighbor_index])
		{
			// Suppress non-maximum pixels
			new_data[i] = 0;
			continue;
		}
		new_data[i] = data[i];
	}
}

// Apply thresholding to detect edges in the image
void Game::apply_threshold(unsigned char *source_data, unsigned char *output_data, std::vector<std::vector<int>> &neighbor_offsets, int pixel_index)
{
	const int EDGE_VALUE = 255;
	const int NON_EDGE_VALUE = 0;

	// Extract row, column, and color information from the pixel index
	int row_idx_3d = pixel_index / (NUM_OF_COLORS * WIDTH);
	int col_idx_3d = (pixel_index % (NUM_OF_COLORS * WIDTH)) / NUM_OF_COLORS;
	int color_idx_3d = pixel_index % NUM_OF_COLORS;

	// Check if the pixel value exceeds the upper threshold
	if (source_data[pixel_index] > THRESHOLD_UPPER_BOUND)
	{
		// If so, mark it as an edge pixel
		output_data[pixel_index] = EDGE_VALUE;
		return;
	}

	// Check if the pixel value falls below the lower threshold
	else if (source_data[pixel_index] < THRESHOLD_LOWER_BOUND)
	{
		// If so, mark it as a non-edge pixel
		output_data[pixel_index] = NON_EDGE_VALUE;
		return;
	}

	bool weak_to_strong = false;

	// Check neighboring pixels for strong edges
	for (std::vector<int> &neighbor : neighbor_offsets)
	{
		// Determine the index of the neighboring pixel
		int neighbor_index = determine_neighbor_index(row_idx_3d, col_idx_3d, color_idx_3d, neighbor[0], neighbor[1]);

		// If a neighboring pixel is a strong edge, mark the current pixel as an edge pixel
		if (neighbor_index != -1 && source_data[neighbor_index] == EDGE_VALUE)
		{
			output_data[pixel_index] = EDGE_VALUE;
			weak_to_strong = true;
			return;
		}
	}

	// If no strong edge is found among neighbors, mark the current pixel as a non-edge pixel
	if (!weak_to_strong)
		output_data[pixel_index] = NON_EDGE_VALUE;
}

// This function performs edge detection using hysteresis
void Game::hysteresis(unsigned char *source_image, unsigned char *edge_image)
{

	// Allocate memory for edge maps in different directions
	unsigned char *diagonal_edge_1 = new unsigned char[DATA_SIZE];
	unsigned char *diagonal_edge_2 = new unsigned char[DATA_SIZE];
	unsigned char *horizontal_edge = new unsigned char[DATA_SIZE];
	unsigned char *vertical_edge = new unsigned char[DATA_SIZE];

	// Define neighbor directions for edge detection
	std::vector<std::vector<int>> neighbor_directions{
		{-1, 1},  // top-right
		{1, -1},  // bottom-left
		{0, -1},  // left
		{0, 1},	  // right
		{-1, -1}, // top-left
		{1, 1},	  // bottom-right
		{-1, 0},  // top
		{1, 0}	  // bottom
	};

	std::vector<std::vector<int>> &neighbor_directions_ref = neighbor_directions;

	// Apply thresholding for edge detection in different directions
	for (int i = 0; i < DATA_SIZE; i++)
	{
		// Detect diagonal edge (top-left to bottom-right)
		apply_threshold(source_image, diagonal_edge_1, neighbor_directions_ref, i);

		// Detect diagonal edge (bottom-right to top-left)
		apply_threshold(source_image, diagonal_edge_2, neighbor_directions_ref, (DATA_SIZE - 1) - i);
	}

	// Traverse the image and apply thresholding for edge detection in remaining directions
	int row = 0;
	while (row < HEIGHT)
	{
		for (int column = (WIDTH * NUM_OF_COLORS) - 1; column >= 0; column--)
		{
			// Calculate index in 1D array
			int index = row * WIDTH * NUM_OF_COLORS + column;

			// Detect horizontal edge
			apply_threshold(source_image, horizontal_edge, neighbor_directions_ref, index);

			// Detect vertical edge
			apply_threshold(source_image, vertical_edge, neighbor_directions_ref, (DATA_SIZE - 1) - index);
		}

		row += 1;
	}

	// Combine edge information from all directions
	for (int i = 0; i < DATA_SIZE; i++)
	{
		edge_image[i] = diagonal_edge_1[i] + diagonal_edge_2[i] + horizontal_edge[i] + vertical_edge[i];
	}

	// Free dynamically allocated memory
	delete[] diagonal_edge_1;
	delete[] diagonal_edge_2;
	delete[] horizontal_edge;
	delete[] vertical_edge;
}

// Apply halftone pattern to a single pixel
void Game::apply_halftone_to_pixel(unsigned char *source_data, unsigned char *output_data, int pixel_index, const std::vector<std::vector<unsigned char>> &halftone_patterns)
{
	// Calculate the row and column numbers in the source and output data
	int row_num_in_source_data = pixel_index / (NUM_OF_COLORS * WIDTH);
	int column_num_in_source_data = (pixel_index % (NUM_OF_COLORS * WIDTH)) / NUM_OF_COLORS;
	int row_num_in_output_data = 2 * row_num_in_source_data;
	int column_num_in_output_data = 2 * column_num_in_source_data;

	// Iterate over each color channel
	for (int color = 0; color < NUM_OF_COLORS; color++)
	{
		// Calculate the index in the output data for the current color channel
		int output_index = NUM_OF_COLORS * 2 * WIDTH * row_num_in_output_data + NUM_OF_COLORS * column_num_in_output_data + color;

		// Calculate the index of the halftone pattern for the current color channel
		int halftone_pattern_index = ((float)source_data[pixel_index + color] / 256) * 5;

		// Retrieve the halftone pattern for the current color channel
		const std::vector<unsigned char> &halftone_pattern = halftone_patterns[halftone_pattern_index];

		// Apply the halftone pattern to the output data
		output_data[output_index] = halftone_pattern[0];
		output_data[output_index + NUM_OF_COLORS] = halftone_pattern[1];
		output_data[output_index + NUM_OF_COLORS * 2 * WIDTH] = halftone_pattern[2];
		output_data[output_index + NUM_OF_COLORS * 2 * WIDTH + NUM_OF_COLORS] = halftone_pattern[3];
	}
}

// Apply halftone pattern to the entire image
void Game::halftone(unsigned char *source_data, unsigned char *output_data)
{
	// Define halftone patterns
	const std::vector<std::vector<unsigned char>> halftone_patterns{
		{0, 0, 0, 0},
		{0, 0, 255, 0},
		{0, 255, 255, 0},
		{0, 255, 255, 255},
		{255, 255, 255, 255},
	};

	// Iterate over each pixel in the image
	for (int i = 0; i < DATA_SIZE; i += NUM_OF_COLORS)
	{
		// Apply halftone pattern to the pixel
		apply_halftone_to_pixel(source_data, output_data, i, halftone_patterns);
	}
}

void Game::floyd_steinberg_pixel(std::vector<std::vector<float>> &output_values, int row_index, int column_index, std::vector<float> &color_palette)
{
	// Get the original color value of the pixel
	float original_value = output_values[row_index][column_index];

	// Determine the closest color from the palette for the original value
	float new_value = color_palette[(int)(original_value / 16)];

	// Calculate the error between the original value and the new value
	float error = original_value - new_value;

	// Assign the new value to the pixel
	output_values[row_index][column_index] = new_value;

	// Propagate the error to neighboring pixels using the Floyd-Steinberg weights
	if (column_index + NUM_OF_COLORS < WIDTH * NUM_OF_COLORS)
		output_values[row_index][column_index + NUM_OF_COLORS] += error * 7 / 16;

	if (row_index + 1 < HEIGHT)
	{
		if (column_index - NUM_OF_COLORS >= 0)
			output_values[row_index + 1][column_index - NUM_OF_COLORS] += error * 3 / 16;

		output_values[row_index + 1][column_index] += error * 5 / 16;

		if (column_index + NUM_OF_COLORS < WIDTH * NUM_OF_COLORS)
			output_values[row_index + 1][column_index + NUM_OF_COLORS] += error * 1 / 16;
	}
}

// Apply Floyd-Steinberg dithering algorithm to the entire image
void Game::floyd_steinberg(unsigned char *source_data, unsigned char *output_data)
{
	// Matrix representing the float values of the output data
	std::vector<std::vector<float>> output_float_values(HEIGHT, std::vector<float>(WIDTH * NUM_OF_COLORS));

	// Define the color palette with 16 color options
	std::vector<float> color_palette(16);
	for (int i = 0; i < 16; i++)
	{
		color_palette[i] = ((float)i) * 256 / 16;
	}

	// Apply Floyd-Steinberg dithering to each pixel in the image
	for (int i = 0; i < DATA_SIZE; i++)
	{
		int row_index = i / (NUM_OF_COLORS * WIDTH);
		int column_index = (i % (NUM_OF_COLORS * WIDTH)) / NUM_OF_COLORS;

		// Assign the original color to the pixel
		output_float_values[row_index][column_index + (i % NUM_OF_COLORS)] = (float)((int)source_data[i]);

		// Apply Floyd-Steinberg dithering to the pixel
		floyd_steinberg_pixel(output_float_values, row_index, column_index + (i % NUM_OF_COLORS), color_palette);

		// Convert the pixel's new color to unsigned char and assign it to the output data
		output_data[i] = (unsigned char)((int)output_float_values[row_index][column_index + (i % NUM_OF_COLORS)]);
	}
}

void Game::print_matrix(unsigned char *data, const std::string file_name, int type_of_image = NORMAL)
{
	// Constants for modulo operations
	const int MOD_VALUES[] = {256, 16, 2};

	// Open file for writing
	std::ofstream output_file(file_name);

	// Loop through each row of pixels
	for (int row = 0; row < HEIGHT; ++row)
	{
		// Loop through each pixel in the row
		for (int col = 0; col < WIDTH * NUM_OF_COLORS; col += NUM_OF_COLORS)
		{
			// Extract color components (RGBA)
			int red = (int)data[row * WIDTH * NUM_OF_COLORS + col] % MOD_VALUES[type_of_image];
			int green = (int)data[row * WIDTH * NUM_OF_COLORS + col + 1] % MOD_VALUES[type_of_image];
			int blue = (int)data[row * WIDTH * NUM_OF_COLORS + col + 2] % MOD_VALUES[type_of_image];
			int alpha = (int)data[row * WIDTH * NUM_OF_COLORS + col + 3] % MOD_VALUES[type_of_image];

			// Determine if it's the last pixel of the last row
			std::string delimiter = ",";
			if (row == HEIGHT - 1 && col == WIDTH * NUM_OF_COLORS - 4)
			{
				delimiter = "";
			}

			// Write color components to file
			output_file << red << "," << green << "," << blue << "," << alpha << delimiter;
		}
		// Add newline after each row
		output_file << "\n\n";
	}

	// Close the file
	output_file.close();
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
