#include "constants.h"
#include "tensor.h"
#include "crystal.h"
#include "simulation.h"
#include "hamiltonian.h"
#include "density_matrix.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <omp.h>


using namespace std;


// The main program first creates four objects: a Simulation, a Crystal, a
// Hamiltonian and a Density_matrix. Creation of these objects involves reading
// the input file, and writing some information in the output file.
// Once the objects are created, the numerical integration of the LVN equation
// can be performed.

// Usage is ssd <name>, with the input file named <name>.ssdi  
int main(int argc, char *argv[]) {


    if (argc != 2) {
        cerr << "\nError:  the format for calling this program is: ssd name\n";
        exit(1);
    }


    string name = argv[1];
    string input_filename = name + ".ssdi"; // name of the input file
    string output_filename = name + ".ssdo"; // name of the output file

    // A simulation is first constructed. The constructor reads the parameters
    // of the simulation from the input file, and reproduces them in the
    // output file.
    Simulation simulation(input_filename, output_filename);


    // The Crystal constructor then reads the input file to extract the
    // coordinates of the nuclei.
    // It also builds the simulation cell, and writes some information to the
    // output file.
    Crystal crystal(input_filename, output_filename);


    // The hamiltonian can then be created. By default, the orientation of the
    // Crystal and rotor frame are the same.
    // The A_20 component of the dipolar coupling hamiltonian is written for
    // each pair of spin in the output file.
    //int number_threads = omp_get_max_threads();
    //Hamiltonian *hamiltonian = (*Hamiltonian) malloc(sizeof(Hamiltonian) * number_threads);
    Hamiltonian hamiltonian(crystal, output_filename, crystal.S);


    // Finally the density matrix can be created.
    Density_matrix sigma(crystal.number_of_nuclei, crystal);

    // Simulate the dynamics and write to .dat files using "name" to
    // name the files
    printf("Starting simulation\n");
    simulation.simulate_dynamics(name, crystal, hamiltonian, sigma);

    return (0);

}
