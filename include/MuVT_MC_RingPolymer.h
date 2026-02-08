#ifndef MUVT_MC_RINGPOLYMER_H
#define MUVT_MC_RINGPOLYMER_H

// Must include parent class header
#include "MuVT_MC_LinearPolymer.h"

class MuVT_MC_RingPolymer : public MuVT_MC_LinearPolymer
{
public:
    // Constructor: pass parameters directly to parent class
    MuVT_MC_RingPolymer(std::string configuration,
                        const double mu_b,
                        const int M,
                        int init_N,
                        double rho_b,
                        double box_xy,
                        double H,
                        double rcut,
                        int max_N);

    // Destructor
    ~MuVT_MC_RingPolymer() override;

    // --- Override core logic ---
    void init_second() override;

    // 1. Check configuration: need to check connection between 0 and M-1
    bool check_configure_validity() override;

    // 2. Rotation move: ring has no end rotation, all are middle rotation
    void rot_polymer_move(int polymer_index) override;

    // 3. Middle monomer rotation: can directly inherit parent class implementation
    // void rot_mid_move(int monomer_index, int polymer_index) override;

    // 4. Insert move: last monomer must close the ring
    void insert_move(int k_max) override;

    // 5. Delete move: ring polymer deletion logic
    void delete_move(int k_max, int delete_polymer_index) override;

    // 6. Recursive insert move: use insert_recursive_ring method
    void insert_move_recursive_ring(int k_max);

    void build_topology_map() override;

    bool check_collision_except_monomer( const double* r_try , int monomer_index); // Do not detect collision with monomer_index
    bool insert_one_monomer(std::vector<std::array<double, 3>> &r_new, double *W, int monomer_index, int k_max); // Insert single monomer
    void delete_one_monomer(std::vector<std::array<double, 3>> &r_delete, double *W, int monomer_index, int k_max); // Delete single monomer

    double get_W_insert_ring(int k_max);
    double get_G_insert_ring(double insert_z, int k_max);
    bool insert_recursive_ring(int next_idx, int parent_idx, double &Z_eff, std::vector<std::array<double, 3>> &r_new, std::vector<int> &is_inserted, int k_max);

    // Calculate insertion weight for polymer, used for chemical potential calculation
    // double calculate_insertion_weight(int k_max) override;
};
#endif