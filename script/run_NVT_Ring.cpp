#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept> // 不清楚干什么的
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <vector>
#include <array>
#include <cstdlib> //中文
#include <format>

#include "MuVT_MC_RingPolymer.h"
#include "Utils.h"
#include "Data_Process.h"
#include "Ensemble_Data_Type.h"




int main(int argc, char *argv[])
{
    /*
    parameter file & data file reading
    */
    // read from input

    // std::string polymer_style(argv[1]); // linear or ring
    
    std::string polymer_topology = "ring";

    std::string input_root_address = "D:\\Desktop\\NewMuVTPolymer\\input\\";
    //std::string input_root_address(argv[2]);

    std::string output_root_address ="D:\\Desktop\\NewMuVTPolymer\\output\\";
    
    //std::string output_root_address(argv[3]);

    std::string parameter_file = "par_ring.dat";

    // std::string parameter_file(argv[4]); // parameter_file

    // std::string parameter_file =  "D:\\Desktop\\NewMuVTPolymer\\input\\par_ring.dat"; 
    
    // std::string configuration_address = "D:\\Desktop\\NewMuVTPolymer\\input\\TrefoilN64M40H20.000000rho0.100000.dat"; 

    std::string configuration_address = "test_1.dat";
    // std::string configuration_address(argv[5]); // initial configuration address

    configuration_address = input_root_address  +configuration_address; // + "init_config_Z\\" 
    
    parameter_file = input_root_address + parameter_file;

    // parameter processing ---------------------------------------------------------------------------------------------
    // def parameter
    double mu_b; 

    double box_xy;

    double M;

    double init_N;

    double dz;

    double rho_b;

    double H;

    // read parameter from file
    std::ifstream infile(parameter_file);
    std::string label;

    // read index
    // mu_b -> M -> init_N -> box_xy -> H -> dz -> An[4] -> Vlin_par1[4] -> Vlin_par2[4] -> xn1[4] -> xn2[4] -> phi_n[4]

    read_line(infile, "M", 1, &M);
    read_line(infile, "mu_b", 1, &mu_b);
    read_line(infile, "rho_bulk", 1, &rho_b);

    read_line(infile, "N_init", 1, &init_N);
    read_line(infile, "H", 1, &H);
    // read_line(infile, "H_permit", 1, &H_permit);
    read_line(infile, "box_xy", 1, &box_xy);
    read_line(infile, "dz", 1, &dz);

    std::cout << "Parameter reading finished." << std::endl;
    
    //---------------------------------------------------------------------------------------------------------------

    // simulation super-parameters
    double rc = 1.0;

    int max_N = round(2 * box_xy * box_xy * H / M); // permit max N

    // simulation length

    int sample_interval = 10;
    int sample_time = 20000;

    int sample_block = 8;

    int block = 0;


    std::string char_common = std::format("RHO{:.2f}H{:.2f}", rho_b, H);

    std::string rho_sample_address =  output_root_address +"rho_profile_" +char_common + "_.dat" ; // "test_sample_ring.dat";
    std::string trace_address = output_root_address + "trace_" + char_common + "_.dat" ; // 
    std::string rg_address = output_root_address + "RG_" + char_common + "_.dat" ; // 
    // mc simulation go

//---------------------------------------------------------------------------------------------------------------------------------------
    //函数参数列表  MuVT_MC_LinearPolymer(std::string configuration_address, const double mu_b, const int M, double int N_bulk,double rho_b, double box_xy, double H,double rcut, int max_N);
    MuVT_MC_RingPolymer mc_sys(
        configuration_address, 
        mu_b, 
        static_cast<int>(M), 
        static_cast<int>(init_N),
        rho_b ,  
        box_xy, 
        H,
        rc, 
        static_cast<int>(max_N));
    
    mc_sys.init_second();

    
    Data_Process rho_profile_proc (block , rho_sample_address);
    Data_Process trace_proc( block , trace_address );
    Data_Process RG_proc( block , rg_address );

    double ave_constant = dz * box_xy * box_xy;
    
    Ensemble_Data_Profile rho_profile("rho_profile", static_cast<int>(H/dz+1), block, ave_constant);

    // Ensemble_Data N;
    
    Data_Trace_Corr RG_z( "RG" , static_cast<int>(init_N));

    double *z  = new double[static_cast<int>(init_N)];
    double *rg = new double[static_cast<int>(init_N)];
    

    

    


//----------------------------------------------------------------------------------------------------------------------------------------
    
    // 调用封装的NVT模拟函数

    
for(int block = 0 ; block < sample_block  ;block ++  )
{
    rho_profile_proc.begin_block();
    trace_proc.begin_block();
    RG_proc.begin_block();
    auto start = std::chrono::high_resolution_clock::now();
    
    // 主模拟循环
    for(int stp = 0 ; stp < sample_time ; stp++)
    {
        // 执行NVT模拟步骤
        mc_sys.mc_one_step_NVT();
        
        // 每隔一定步数记录数据
        if(stp % sample_interval == 0)
        {
            // 计算并记录密度分布
            rho_profile_proc.cal_rho_seg_profile(mc_sys, dz, rho_profile.profile_data);
            rho_profile.sample_time++;
            
            // 记录轨迹
            trace_proc.record_trace(mc_sys, stp/sample_interval);

            // 记录
            RG_proc.get_RG_for_all_polymer(mc_sys , rg);
            RG_proc.get_rc_z_for_all_polymer(mc_sys, z);
            RG_z.sample_one_frame(z ,  rg , init_N);       
            RG_proc.pipe_to_file(RG_z,stp/sample_interval);

        }
    }

// ------------------
    mc_sys.end_block(block);
    
    rho_profile.get_average(); // 计算
    rho_profile_proc.pipe_to_file(rho_profile);
    
    
    rho_profile.reset();
    
    rho_profile_proc.end_block(block+1) ;

    RG_proc.end_block(block+1);




// -------------------
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "block("<<block<<") " <<"computer running time: " << duration.count() / 1000.0/1000.0 << " s" << std::endl;
}

    return 0;
}