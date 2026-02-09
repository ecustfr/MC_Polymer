#ifndef UTILS_H
#define UTILS_H

#include <random>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <array>
// Record 2D data
void sample_item(std::ofstream &file, bool time_or_not, const int time, double **data, const int row, const int col);

// Record 1D data
void sample_item(std::ofstream &file, bool time_or_not, const int time, double *data, const int row);

// Record single data
void sample_item(std::ofstream &file, bool time_or_not, const int time, double &data);

void generate_random_unit_vector(double *added_vec , std::mt19937 &gen);

double dis2_no_period(const double *x, const double *y);

double dis_period(double x, double y, double box);

int RouteWheel(double *P_try , int k_max , std::mt19937 &gen);

double sum_1d(const double *data, const int length);

void add_random_unit_vector(double *added_vec , std::mt19937 &gen);
bool acc_bias_or_not( double*  new_weight , double old_weight , int k_max , std::mt19937 &gen);

void read_line(std::ifstream &infile, const std::string label, int max_index, double *array);

inline double SQR(double x) {return x*x;}

inline double distance_array(const std::array<double, 3>& p1, const std::array<double, 3>& p2) 
{
    double dx = p1[0] - p2[0];
    double dy = p1[1] - p2[1];
    double dz = p1[2] - p2[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Find last monomer position for ring polymer
extern "C" void Find_Last(double *r1, double *r2, double **target, int kmax, double diameter);

// Calculate path probability P_road
double calculate_P_road(const double* begin_r, const double* end_r, double diameter, int try_num, double box_size_xy);

// Table version of calculate_P_road class
class P_Road_Lookup_Table {
private:
    static const int MAX_R = 10000;
    static const int MAX_TRY_NUM = 12;
    double lookup_table[MAX_TRY_NUM + 1][MAX_R + 1];
    bool initialized;

public:
    P_Road_Lookup_Table();
    void initialize();
    double get_P_road(double R, int try_num);
};

// Global table instance
extern P_Road_Lookup_Table g_p_road_table;

#endif