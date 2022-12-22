#ifndef SIMULATION_CLASS
#define SIMULATION_CLASS

#include "crystal.h"
#include "hamiltonian.h"
#include "density_matrix.h"
#include <string>

class Simulation {

public:

    Simulation(std::string input_filename, std::string output_filename);

    // Function that is called by the main file and that runs the actual
    // simulation
    void simulate_dynamics(std::string name,
                           Crystal &crystal,
                           Hamiltonian &hamiltonian,
                           Density_matrix &sigma);

private:
    int configuration;
    int continuation;

    bool tree_pruning;

    double timestep;        // Timestep in ms

    int number_of_steps;    // Number of steps for the simulation

    int sampling_interval;
    // Data are written every such number of steps. By default, data
    // is written at every step.

    double omega_r;    // Spinning frequency

    // The result of the time evolution is written to a .dat file.  The
    // extension for the .dat file tells the highest spin order included
    // in the reduced Liouville space.
    std::string dat_extension;

    int number_of_crystallites;    // Number of crystallites taken into
    // account in the powder averaging scheme

    double **orientations;        // 2D array storing the orientations
    // and weights for powder averaging. The first dimension stores the index of
    // the orientation, and the second dimension stores alpha, beta, gamma and
    // weight in this order.

    std::string angle_set;    // Name of the angle set used for powder
    // averaging. By default, a single crystallite is considered, with the B
    // field along the z axis in  the frame used to define the structure.

    // The nuclei initially polarised can be
    // specified by a sequence of 4 indices.  The first three specify
    // which unit cell the nuclei is in, while the fourth gives the index
    // of the nuclei within its unit cell.  For instance, if
    // polarised_nuclei has the value 1, 2, 3, 4, then the unit cell
    // containing the nuclei is in the first position along the e vector, the
    // second position along the b vector, the third position along the
    // c vector; and the nuclei is the 4th (relevant) nuclei listed in the text
    // file which gives the crystal structure of the unit cell.
    // The 2D array polarised_nuclei stores for each nuclei the four integers
    // describe above, and in addition, as a fifth element, the index of
    // the nuclei with respect to the supercell
    std::vector <std::vector<int>> polarised_nuclei;


    // Function that opens an input file to extract the values for the member
    // variables of a simulations.
    void read_run_parameters(std::string input_filename);

    // Function that opens an output file and write the parameters, as a
    // summary of the simulation and to check that io is handled properly.
    void write_run_parameters(std::string output_filename);

    void read_number_of_crystallites(std::string orientations_filename);

    // Function that opens the file storing the orientations and weight for
    // powder averaging
    void read_orientations(std::string orientations_filename);


    // Function called by the public function simulate_dynamics if only
    // a single crystallite is included in the simulation
    void simulate_crystallite(std::vector<int>& polarised_spins,
                              std::string filename_prefix,
                              Crystal &crystal,
                              Hamiltonian &hamiltonian,
                              Density_matrix &sigma);

    // Function called by the public function simulate_dynamics if
    // a powder average is desired.
    void simulate_powder_average(std::vector<int>& polarised_spins,
                                 std::string filename_prefix,
                                 Crystal &crystal,
                                 Hamiltonian &hamiltonian,
                                 Density_matrix &sigma);

};

#endif
