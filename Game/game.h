#pragma once
#include "scene.h"
#include <functional>

class Game : public Scene
{
public:
    Game();
    Game(float angle, float relationWH, float near, float far);
    ~Game(void);

    void Init();
    void Update(const glm::mat4& MVP, const glm::mat4& Model, const int shaderIndx);
    void WhenRotate();
    void WhenTranslate();
    void Motion();

private:
    void print_matrix(unsigned char* data, const std::string file_name, int type_of_image);
    void edge_detection(unsigned char* data, unsigned char* new_data);
    void smoothing(unsigned char* data, unsigned char* new_data);
    void halftone(unsigned char* data, unsigned char* new_data);
    void floyd_steinberg(unsigned char* data, unsigned char* new_data);
    void non_max_suppression(unsigned char* data, unsigned char* new_data, float* pixel_theta);
    void compute_gradient(unsigned char* data, unsigned char* new_data, float* pixel_theta);
    void hysteresis(unsigned char* data, unsigned char* new_data);
    void convolution(unsigned char* data, unsigned char* new_data, std::vector<std::vector<float>>& kernel);
    void apply_threshold(unsigned char* data, unsigned char* new_data, std::vector<std::vector<int>>& neighbors_ref, int index);
    void floyd_steinberg_pixel(std::vector<std::vector<float>>& new_data_float_values, int row_num, int column_num, std::vector<float>& colors);
    void apply_halftone_to_pixel(unsigned char* data, unsigned char* new_data, int pixel_num, const std::vector<std::vector<unsigned char>>& halftone_patterns);
    int determine_neighbor_index(int x, int y, int z, int up_down, int left_right);
};
