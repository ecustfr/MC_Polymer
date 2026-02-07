#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

#include "MuVT_MC_LinearPolymer.h"
#include "MuVT_MC_RingPolymer.h"

namespace py = pybind11;

PYBIND11_MODULE(pymcpolymer, m) {
    m.doc() = "Monte Carlo Polymer Simulation Module";
    m.attr("__version__") = "0.1.0";

    // Add a simple function for testing
    m.def("add", [](int i, int j) {
        return i + j;
    }, "A simple function that adds two numbers");

    // Bind MuVT_MC_LinearPolymer class
    py::class_<MuVT_MC_LinearPolymer>(m, "MuVT_MC_LinearPolymer")
        .def(py::init<std::string, double, int, int, double, double, double, double, int>(),
             py::arg("configuration"),
             py::arg("mu_b"),
             py::arg("M"),
             py::arg("init_N"),
             py::arg("rho_b"),
             py::arg("box_xy"),
             py::arg("H"),
             py::arg("rcut"),
             py::arg("max_N"))
        .def("init_second", &MuVT_MC_LinearPolymer::init_second)
        .def("set_sim_parameters", &MuVT_MC_LinearPolymer::set_sim_parameters,
             py::arg("EPS_TRANS"),
             py::arg("ROT_RATIO"),
             py::arg("K_MAX"))
        .def("print_all_parameters", &MuVT_MC_LinearPolymer::print_all_parameters)
        .def("end_block", &MuVT_MC_LinearPolymer::end_block,
             py::arg("block"))
        .def("reset_mc_record", &MuVT_MC_LinearPolymer::reset_mc_record)
        .def("mc_one_step_NVT", &MuVT_MC_LinearPolymer::mc_one_step_NVT)
        .def("mc_one_step_MuVT", &MuVT_MC_LinearPolymer::mc_one_step_MuVT)
        // Bind getter methods
        .def("get_MN_now", &MuVT_MC_LinearPolymer::get_MN_now)
        .def("get_N_now", &MuVT_MC_LinearPolymer::get_N_now)
        .def("get_rho_now", &MuVT_MC_LinearPolymer::get_rho_now)
        .def("get_M", &MuVT_MC_LinearPolymer::get_M)
        .def("get_H", &MuVT_MC_LinearPolymer::get_H)
        .def("get_box_xy", &MuVT_MC_LinearPolymer::get_box_xy)
        .def("get_max_N", &MuVT_MC_LinearPolymer::get_max_N)
        // Add acceptance rate getter methods
        .def("get_trans_acceptance", &MuVT_MC_LinearPolymer::get_trans_acceptance)
        .def("get_rot_acceptance", &MuVT_MC_LinearPolymer::get_rot_acceptance)
        .def("get_insert_acceptance", &MuVT_MC_LinearPolymer::get_insert_acceptance)
        .def("get_delete_acceptance", &MuVT_MC_LinearPolymer::get_delete_acceptance)
        // Expose r_total as numpy array
        .def_property_readonly("r_total", [](const MuVT_MC_LinearPolymer &self) {
            int rows = self.get_MN_now();
            int cols = 3;
            
            // Create numpy array
            py::array_t<double> result({rows, cols});
            auto result_buf = result.mutable_unchecked<2>();
            
            // Get r_total pointer
            const double** r_total = self.get_r_total();
            
            // Copy data to numpy array
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    result_buf(i, j) = r_total[i][j];
                }
            }
            
            return result;
        })
        // Add new method bindings
        .def("insert_move", &MuVT_MC_LinearPolymer::insert_move,
             py::arg("k_max"),
             "Insert a polymer with k_max candidate positions")
        .def("delete_move", &MuVT_MC_LinearPolymer::delete_move,
             py::arg("k_max"),
             py::arg("delete_index"),
             "Delete a polymer with k_max candidate positions")
        .def("get_change_volumn", &MuVT_MC_LinearPolymer::get_change_volumn,
             py::arg("expan_or_shrink"),
             py::arg("mc_cond_eps_volumn"),
             "Volume change and overlap detection")
        .def("get_W_insert", &MuVT_MC_LinearPolymer::get_W_insert,
             py::arg("k_max"),
             "Calculate insertion weight for polymer")
        .def("get_G_insert", &MuVT_MC_LinearPolymer::get_G_insert,
             py::arg("insert_z"),
             py::arg("first_insert_index"),
             py::arg("k_max"),
             "Calculate insertion weight with specified z position and first insert index")

        .def("get_Wz_insert", &MuVT_MC_LinearPolymer::get_Wz_insert,
             py::arg("insert_z"),
             py::arg("first_insert_index"),
             py::arg("k_max"),
             "Calculate insertion weight at specific z position")

        .def("insert_recursive", &MuVT_MC_LinearPolymer::insert_recursive,
             py::arg("next_idx"),
             py::arg("parent_idx"),
             py::arg("step"),
             py::arg("total_W"),
             py::arg("r_new"),
             py::arg("is_inserted"),
             py::arg("k_max"),
             "Recursively insert monomers")
        // External potential related method bindings
        .def("set_external_potential", py::overload_cast<const std::string&>(&MuVT_MC_LinearPolymer::set_external_potential),
             py::arg("name"),
             "Set external potential by name")
        .def("set_external_potential", py::overload_cast<std::function<double(const double)>, const std::string&>(&MuVT_MC_LinearPolymer::set_external_potential),
             py::arg("potential"),
             py::arg("name"),
             "Set external potential by function")
        .def("Vext_pot", &MuVT_MC_LinearPolymer::Vext_pot,
             py::arg("z"),
             "Calculate external potential at position z")
        .def("Vext_bz", &MuVT_MC_LinearPolymer::Vext_bz,
             py::arg("pos"),
             "Calculate Boltzmann factor for external potential")
        .def("export_potential_table", &MuVT_MC_LinearPolymer::export_potential_table,
             py::arg("filename"),
             "Export potential table to file")
             
        .def("get_potential_table", &MuVT_MC_LinearPolymer::get_potential_table,
             "Get potential table values")
        .def("get_potential_table_z", &MuVT_MC_LinearPolymer::get_potential_table_z,
             "Get z coordinates corresponding to potential table values");


    // Bind MuVT_MC_RingPolymer class (inherits from MuVT_MC_LinearPolymer)
    py::class_<MuVT_MC_RingPolymer, MuVT_MC_LinearPolymer>(m, "MuVT_MC_RingPolymer")
        .def(py::init<std::string, double, int, int, double, double, double, double, int>(),
             py::arg("configuration"),
             py::arg("mu_b"),
             py::arg("M"),
             py::arg("init_N"),
             py::arg("rho_b"),
             py::arg("box_xy"),
             py::arg("H"),
             py::arg("rcut"),
             py::arg("max_N"));

    // Overridden MC move methods are automatically inherited, no need to rebind

    
}