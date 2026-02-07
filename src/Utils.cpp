#define _USE_MATH_DEFINES
#include "Utils.h"

#include <random>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>


// 记录二维数据
void sample_item(std::ofstream &file, bool time_or_not, const int time, double **data, const int row, const int col) // 是否记录时间
{

    if (file.is_open())
    {
        if (time_or_not)
        {
            file << "Time: " << time << std::endl;
        }
        for (int par_index = 0; par_index < row; par_index++)
        {
            for (int j = 0; j < col; j++)
            {
                file << data[par_index][j] << " ";
            }
            file << std::endl;
        }
    }
}
// 记录一维数据
void sample_item(std::ofstream &file, bool time_or_not, const int time, double *data, const int row)
{
    if (file.is_open())
    {
        if (time_or_not)
        {
            file << "Time: " << time << std::endl;
        }
        for (int par_index = 0; par_index < row; par_index++)
        {
            file << data[par_index] << " ";
        }
        file << std::endl;
    }
}
// 记录一个数据
void sample_item(std::ofstream &file, bool time_or_not, const int time, double &data)
{
    if (file.is_open())
    {
        file << data << std::endl;
    }
}

void generate_random_unit_vector(double *vec,  std::mt19937 &gen)
{
    std::uniform_real_distribution<double> real_dis(-1.0,1.0);
    double len_sq;
    do
    {
        vec[0] = -real_dis(gen);
        vec[1] = -real_dis(gen);
        len_sq = vec[0]*vec[0] + vec[1]*vec[1];
    } while (len_sq > 1);

    double ranh = 2*sqrt(1 - len_sq);
    vec[0] = vec[0] * ranh;
    vec[1] = vec[1] * ranh;
    vec[2] = 1-2*len_sq;
}

double dis2_no_period(const double *x, const double *y)
{
    double dis2 = 0.0; 
    for(int i =0 ;i<3 ;i++)
    {
        dis2 += SQR(x[i] - y[i]);
    }
    return dis2;
}

double dis_period(double x, double y, double box)
{
    return (x-y) - round((x-y)/box)*box;
}

int RouteWheel(double *P_try , int k_max, std::mt19937 &gen)
{
    std::uniform_real_distribution<double> uni_dis(0.0,1.0);

    double sum =0.0;
    for(int i =0 ; i< k_max ; ++i) 
    {
        sum += P_try[i];
    }
    double cutoff = uni_dis(gen)*sum;
    for(int i =0 ;i <k_max ; i++)
    {
        cutoff -= P_try[i];
        if(cutoff <= 0)
        {
            return i;
        }
    }
    return -1;
}

double sum_1d(const double *data, const int length)
{
    double sum =0.0;
    for( int i = 0; i<length ; i++)
    {
        sum += data[i];
    }
    return sum;
}

void add_random_unit_vector(double *added_vec , std::mt19937 &gen)
{
    std::uniform_real_distribution<double> real_dis(-1.0,1.0); 
    double vec[3];
    double len_sq;
    do
    {
        vec[0] = real_dis(gen);
        vec[1] = real_dis(gen);
        len_sq = vec[0]*vec[0] + vec[1]*vec[1];
        /* code */
    } while (len_sq > 1);

    double ranh = 2*sqrt(1-len_sq);
    vec[0] = vec[0]*ranh;
    vec[1] = vec[1]*ranh;
    vec[2] = 1-2*len_sq;

    added_vec[0] += vec[0];
    added_vec[1] += vec[1];
    added_vec[2] += vec[2];
    
}

bool acc_bias_or_not( double*  new_weight , double old_weight , int k_max ,  std::mt19937 &gen)
{

    std::uniform_real_distribution<double> uni_dis(0.0,1.0);

    double new_weight_sum =0.0;
    double old_weight_sum = old_weight;
    
    for( int i =0 ; i< k_max ; i++)
    {
        new_weight_sum += new_weight[i];
    }
    
    for( int i =0 ; i< k_max-1 ; i++)
    {
        old_weight_sum += new_weight[i + k_max];
    }

    if( new_weight_sum / old_weight_sum > uni_dis(gen) )
    {
        return true;
    }
    else
    {
        return false;
    }
}

void read_line(std::ifstream &infile, const std::string label, int max_index, double *array) // 从文件的一行读取数据，文件流，标签，最大读取数量，存储数组
{
    double temp;
    std::string temp_label;
    int i = 0;
    if (infile.is_open())
    {
        infile >> temp_label;
        if (temp_label != label)
        {
            std::cout << "Expected label: " << label << ", but found: " << temp_label << std::endl;
            return;
        }

        while (infile >> temp)
        {
            array[i] = temp;
            i++;
            if (i >= max_index)
            {break;}
        }
        // infile >> temp_label; // read a blank lineS
    }
    else
    {
        std::cout << "the file is not open!" << std::endl;
    }
    // for test 
    if(max_index ==1 )
    {
        std::cout<< label << " = " << array[0] <<"\n";
    }
    else
    {
        for( int j=0 ; j< max_index ; j++)
        {
            std::cout<< label << "[" << j << "] = " << array[j] <<"  ";
        }
        std::cout << "\n";
    }
        
};



// 计算路径概率 P_road
double calculate_P_road(const double* begin_r, const double* end_r, double diameter, int try_num, double box_size_xy) 
{
    // 计算两点之间的距离
    double dx = begin_r[0] - end_r[0];
    double dy = begin_r[1] - end_r[1];
    double dz = begin_r[2] - end_r[2];
    double R = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
    
    // 使用打表版本计算路径概率
    return g_p_road_table.get_P_road(R, try_num);
}

// P_Road_Lookup_Table 类的实现

// 静态常量定义
const int P_Road_Lookup_Table::MAX_R;
const int P_Road_Lookup_Table::MAX_TRY_NUM;

// 构造函数
P_Road_Lookup_Table::P_Road_Lookup_Table() : initialized(false) {
    initialize();
}

// 初始化打表
void P_Road_Lookup_Table::initialize() {
    if (initialized) {
        return;
    }
    
    // 初始化打表
    for (int try_num = 1; try_num <= MAX_TRY_NUM; try_num++) {
        for (int r_int = 0; r_int <= MAX_R; r_int++) {
            double R = r_int / 1000.0; // 转换为实际距离
            double prob = 0.0;
            
            // 复制 calculate_P_road 中的计算逻辑
            switch (try_num) {
                case 2:
                    prob = (1 + (2 - R > 0 ? 2 : 0)) / 16 / M_PI / R;
                    break;
                case 3:
                    prob = ((3 - R) + fabs(R - 3)) / (32 * M_PI * R);
                    break;
                case 4:
                    prob = ((8 - 3 * R) * R - (R - 4) * fabs(R - 4) + 4 * (R - 2) * fabs(R - 2)) / (128 * M_PI * R);
                    break;
                case 5:
                    prob = (30 * R - 6 * pow(R, 3) + pow(fabs(R - 5), 3) - 5 * pow(fabs(-3 + R), 3) + 10 * pow(fabs(R - 1), 3)) / (768 * M_PI * R);
                    break;
                case 6:
                    prob = (2 * R * (96 + pow(R, 2) * (-24 + 5 * R)) - pow(R - 6, 3) * fabs(R - 6) + 6 * pow(R - 4, 3) * fabs(R - 4) - 15 * pow(-2 + R, 3) * fabs(-2 + R)) / (6144 * M_PI * R);
                    break;
                case 7:
                    prob = (20 * R * (77 - 14 * pow(R, 2) + pow(R, 4)) + pow(fabs(R - 7), 5) - 7 * pow(fabs(R - 5), 5) + 21 * pow(fabs(R - 3), 5) - 35 * pow(fabs(R - 1), 5)) / (61440 * M_PI * R);
                    break;
                case 8:
                    prob = - (1.0 / (737280 * M_PI * R)) * (70 * pow(R, 6) - 56 * pow(2 + R, 6) + 28 * pow(4 + R, 6) - 8 * pow(6 + R, 6) + pow(8 + R, 6) + 56 * pow(-2 + R, 6) * (R < 2 ? 1 : -1) - 28 * pow(-4 + R, 6) * (R < 4 ? 1 : -1) + 8 * pow(-6 + R, 6) * (R < 6 ? 1 : -1) - pow(-8 + R, 6) * (R < 8 ? 1 : -1));
                    break;
                case 9:
                    prob = 1.0 / (10321920 * M_PI * R) * (pow(fabs(-9 + R), 7) - 9 * pow(fabs(-7 + R), 7) + 2.0 * (91035 * R - 13545 * pow(R, 3) + 945 * pow(R, 5) - 35 * pow(R, 7) + 18 * pow(fabs(R - 5), 7) - 42 * pow(fabs(R - 3), 7) + 10 * pow(fabs(R - 1), 7)));
                    break;
                case 10:
                    prob = 1.0 / (165150720 * M_PI * R) * (252 * pow(R, 8) - 210 * pow(2 + R, 8) + 120 * pow(4 + R, 8) - 45 * pow(6 + R, 8) + 10 * pow(8 + R, 8) - pow(10 + R, 8) + 210 * pow(-2 + R, 8) * (R < 2 ? 1 : -1) - 120 * pow(-4 + R, 8) * (R < 4 ? 1 : -1) + 45 * pow(-6 + R, 8) * (R < 6 ? 1 : -1) - 10 * pow(-8 + R, 8) * (R < 8 ? 1 : -1) + pow(-10 + R, 8) * (R < 10 ? 1 : -1));
                    break;
                case 11:
                    prob = 1.0 / (2972712960 * M_PI * R) * (39415068 * R - 4904592 * pow(R, 3) + 293832 * pow(R, 5) - 11088 * pow(R, 7) + 252 * pow(R, 9) + pow(fabs(-11 + R), 9) - 11 * pow(fabs(-9 + R), 9) + 55 * pow(fabs(-7 + R), 9) - 165 * pow(fabs(-5 + R), 9) + 330 * pow(fabs(-3 + R), 9) - 462 * pow(fabs(-1 + R), 9));
                    break;
                case 12:
                    prob = -1.0 / (59454259200 * M_PI * R) * (924 * pow(R, 10) - 792 * pow(2 + R, 10) + 495 * pow(4 + R, 10) - 220 * pow(6 + R, 10) + 66 * pow(8 + R, 10) - 12 * pow(10 + R, 10) + pow(12 + R, 10) + 792 * pow(-2 + R, 10) * (R < 2 ? 1 : -1) - 495 * pow(-4 + R, 10) * (R < 4 ? 1 : -1) + 220 * pow(-6 + R, 10) * (R < 6 ? 1 : -1) - 66 * pow(-8 + R, 10) * (R < 8 ? 1 : -1) + 12 * pow(-10 + R, 10) * (R < 10 ? 1 : -1) - pow(-12 + R, 10) * (R < 12 ? 1 : -1));
                    break;
                default:
                    prob = 1.0 / (4 * M_PI * R);
                    break;
            }
            
            if (prob < 0) {
                prob = 0;
            }
            
            lookup_table[try_num][r_int] = prob;
        }
    }
    
    initialized = true;
}

// 获取路径概率
double P_Road_Lookup_Table::get_P_road(double R, int try_num) {
    if (try_num < 1 || try_num > MAX_TRY_NUM) {
        return 1.0 / (4 * M_PI * R);
    }
    
    int r_int = std::min(static_cast<int>(R * 1000), MAX_R);
    return lookup_table[try_num][r_int];
}

// 全局打表实例
P_Road_Lookup_Table g_p_road_table;


