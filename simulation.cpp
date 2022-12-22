#include "constants.h"
#include "crystal.h"
#include "hamiltonian.h"
#include "density_matrix.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <omp.h>

#define PSD 0
#define SATURATE 1
#define INVERT 2
#define SATURATEZ 3

//#define ZQT

using namespace std;


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


Simulation::Simulation(string input_filename, string output_filename) {

    // Extension added to the .dat file
    dat_extension = EXTENSION;

    // Some variables are initialised to impossible values for error checking
    timestep = -1.0;
    number_of_steps = -1;
    omega_r = -1.0;

    tree_pruning = false;

    configuration = PSD; // default to not saturating.

    // Some variables are set to a default value
    sampling_interval = 1;
    angle_set = string("none");
    continuation = 0;
    // Read the parameters for the run
    read_run_parameters(input_filename);


    // Check that values were read
    if (timestep < 0) {
        cerr << "\nError reading timestep.\n";
        exit(1);
    }

    if (number_of_steps < 0) {
        cerr << "\nError reading number of timesteps.\n";
        exit(1);
    }

    if (omega_r < 0) {
        cerr << "\nError reading spinning frequency.\n";
        exit(1);
    }

    // Check to see whether values were read into polarised_nuclei
    for (unsigned int i = 0; i < polarised_nuclei.size(); i++) {
        for (int j = 0; j < 4; j++) {
            if (polarised_nuclei[i][j] < 0) {
                cerr << "\nError in reading polarised nuclei indices." << endl;
                exit(1);
            }
        }
    }

    // Continue the initialisation. If no powder averaging is required, the
    // number of crystallites is set equal to 1. If a powder
    // averaging scheme is used, the number of orientations is read from the
    // .ang file
    if (angle_set == string("none")) {
        number_of_crystallites = 1;
    } else {
        string orientations_filename = angle_set + ".ang";
        read_number_of_crystallites(orientations_filename);
    }

    // Allocate memory for the orientations
    try {

        orientations = new double *[number_of_crystallites];
        for (int i = 0; i < number_of_crystallites; i++)
            orientations[i] = new double[4];

    }
    catch (bad_alloc) {
        cerr << "\nMemory could not be allocated for the orientations\n";
        exit(1);
    }

    // If no powder averaging is required, the only orientation has a rotor
    // frame identical to the crystal frame. Otherwise, the orientations are
    // read from the ang file and their validity is assessed.
    if (angle_set == string("none")) {

        for (int j = 0; j < 3; j++) {
            orientations[0][j] = 0.0;
        }
        orientations[0][3] = 1.0;

    } else {
        for (int i = 0; i < number_of_crystallites; i++) {
            for (int j = 0; j < 4; j++) {
                orientations[i][j] = -1.0;
            }
        }

        string orientations_filename = angle_set + ".ang";
        read_orientations(orientations_filename);

        for (int i = 0; i < number_of_crystallites; i++) {
            if ((orientations[i][0] < 0.0) || (orientations[i][0] > TWO_PI + 0.1)) {
                cerr << "\naplha must belong to [0, 2 * Pi]\n";
                exit(1);
            }
            if ((orientations[i][1] < 0.0) || (orientations[i][1] > PI + 0.1)) {
                cerr << "\nbeta must belong [0, Pi]\n";
                exit(1);
            }
            if ((orientations[i][2] < 0.0) || (orientations[i][2] > TWO_PI + 0.1)) {
                cerr << "\ngamma must belong to [0, 2 * Pi]\n";
                exit(1);
            }
            if (orientations[i][3] < 0.0) {
                cerr << "\nThe weight of an orientation must be strictly positive\n";
                exit(1);
            }

        } // for

    }   // else


    // Write a summary of the simulation in an output file
    write_run_parameters(output_filename);

} // Constructor

void Simulation::read_run_parameters(string input_filename) {

    // open an ifstream to read the input file
    ifstream input_filestream(input_filename.c_str());

    // Check that the file stream opened normally
    if (input_filestream.fail()) {
        cerr << "\nInput file " << input_filename
             << " could not be opened in main().\n";
        exit(1);
    }

    continuation = 0;
    // The input file is read line by line. For each line the first
    // token is compared to several possibilites to check whether the
    // line contains a parameter. If it does, the parameter is read, and
    // the corresponding variable is then written to the output file
    while (!input_filestream.eof()) {

        // The text is read line by line. Each line is temporarily stored in
        // a string
        string line = "";
        getline(input_filestream, line);

        // Creates a stream from the string line_srtream
        istringstream line_stream(line);

        string token = "";             // A token to be read from the string "line."

        // Read the first token from the line, with tokens separated by whitespace.
        line_stream >> token;


        // Check for the expected tokens, starting with a check for comments
        if (token[0] == '#')
            continue;            // Skip to the next line of text

        else if (token == "_timestep_in_ms") {

            // If the first token is "_timestep_in_ms" then read timestep
            line_stream >> timestep;
            if (line_stream.fail()) {
                cerr << "\nError extracting the timestep from the input file's line of text.\n";
                exit(1);
            }
        } else if (token == "_configuration") {
            string temp;
            line_stream >> temp;
            if (temp == "PSD") { configuration = PSD; }
            else if (temp == "SATURATE") { configuration = SATURATE; }
            else if (temp == "INVERT") { configuration = INVERT; }
            else if (temp == "SATURATEZ") { configuration = SATURATEZ; }
            else {
                cerr << "\nError did not understand configuration '" << temp <<"'\n";
                exit(1);
            }
        } else if (token == "_continuation") {
            line_stream >> continuation;
        } else if (token == "_tree_pruning") {
            tree_pruning = true;
        } else if (token == "_number_of_steps") {

            // If line begins with "_number_of_steps" then read number_of_steps
            line_stream >> number_of_steps;
            if (line_stream.fail()) {
                cerr << "\nError extracting number_of_text from the input file's line of text.\n";
                exit(1);
            }
        } else if (token == "_sampling_interval_in_steps") {

            // If line begins with "_sampling_interval_in_steps" then read
            // the sampling interval
            line_stream >> sampling_interval;
            if (line_stream.fail()) {
                cerr << "\nError extracting sampling_interval from the input file's line of text.\n";
                exit(1);
            }
        } else if (token == "_spinning_speed_in_kHz") {

            line_stream >> omega_r;
            if (line_stream.fail()) {
                cerr << "\nError extracting the spinning speed from the input file's line of text.\n";
                exit(1);
            }
        } else if (token == "_angle_set") {

            // If the first token is "_angle_set" then read the name of the angle set
            line_stream >> angle_set;
            if (line_stream.fail()) {
                cerr << "\nError extracting the name of the angle set from the input file's line of text.\n";
                exit(1);
            }
        } else if (token == "_polarised_nuclei") {

            // When a polarised nuclei is given in the input file, the size of
            // the polarised_nuclei 2D vector is increased, by adding a vector of
            // size 5, initialised to impossible values.
            std::vector<int> vec(5, -1);
            polarised_nuclei.push_back(vec);

            // The index of the newly created vector is used to assign the values
            // that are read from the input file
            int i = polarised_nuclei.size() - 1;

            // Read the coordinates of a polarised nuclei as a set of 4 indices:
            // (a_cell,b_cell,c_cell,index_in_unit_cell)
            for (int j = 0; j < 4; j++)
                line_stream >> polarised_nuclei[i][j]; // stores index j

        } else
            continue;            // None of the expected tokens were found.

    } // while

    input_filestream.close();

} // read_run_parameter


void Simulation::write_run_parameters(string output_filename) {

    // Open an ofstream to store a description of the simulation
    ofstream output_filestream(output_filename.c_str());

    // Check that the file stream opened normally
    if (output_filestream.fail()) {
        cerr << "\nOutput file " << output_filename << " could not be opened.\n";
        exit(1);
    }


    output_filestream << "_timestep_in_ms"
                      << "\t\t\t\t" << timestep << endl;

    output_filestream << "_number_of_steps"
                      << "\t\t\t" << number_of_steps << endl;

    output_filestream << "_sampling_interval_in_steps"
                      << "\t\t" << sampling_interval << endl
                      << endl << endl;

    output_filestream << "_spinning_speed_in_kHz"
                      << "\t\t\t" << omega_r << endl;


    // Convert the frequency in radian/ms
    omega_r *= TWO_PI;

    output_filestream << "_angle_set"
                      << "\t\t\t\t" << angle_set << endl
                      << endl << endl;


    for (unsigned int i = 0; i < polarised_nuclei.size(); i++) {

        output_filestream << "_polarised_nuclei\t\t\t"
                          << polarised_nuclei[i][0] << " "
                          << polarised_nuclei[i][1] << " "
                          << polarised_nuclei[i][2] << " "
                          << polarised_nuclei[i][3] << " " << endl;
    }
    output_filestream << endl << endl;

    output_filestream.close();

} // write_run_parameter


// This function reads the number of crystallites and also does error
// checking for the input file.  Blank lines or lines containing only
// whitespace are acceptable, but any non-whitespace characters must
// be interpretable as 4 real numbers separated by white space.
void Simulation::read_number_of_crystallites(string orientations_filename) {

    // open an ifstream
    ifstream orientations_filestream(orientations_filename.c_str());

    // Check that the file stream opened normally
    if (orientations_filestream.fail()) {
        cerr << "\nThe file " << orientations_filename
             << " could not be opened in read_number_of_crystallites().\n";
        exit(1);
    }

    // The file is read line by line to count the number of crystallites,
    // stored in crystallites_count.  Each line of text is temporarily
    // stored in the string "line."
    string line;

    string token;   // A token to be read from the string "line."

    int crystallites_count = 0;


    while (!orientations_filestream.eof()) {

        // Reset the strings "line" and "token," for cases where the
        // string operations fail.
        token = "";
        line = "";

        getline(orientations_filestream, line);

        // Creates a stream from the string line_srtream
        istringstream line_stream(line);

        // Read the first token from the line, with tokens separated by whitespace.
        line_stream >> token;


        // If the token contains no characters, skip to the next line
        if (token.length() == 0)
            continue;
        else {
            // Check that there are 4 numbers in the line and nothing else.

            // Return to the beginning of the stream
            line_stream.seekg(0, ios::beg);

            for (int i = 0; i < 4; i++) {
                double dummy = 0;
                line_stream >> dummy;
                if (line_stream.fail()) {
                    cerr << "\nError reading angles or weighting for a crystallite. 2\n";
                    exit(1);
                } // check for failure to input a double.

            }      // loop through inputting 4 numbers

            // Check to see whether there is junk left on the line
            token = "";
            line_stream >> token;    // This should be empty.
            if (token.length() != 0) {
                cerr << "\nError:  a line of the .ang file contains characters"
                     << endl << "which can't be interpreted as 4 numbers separated by whitespace.\n"
                     << endl;
                exit(1);
            }

            // This line of the .ang file passed the error test
            crystallites_count++;

        } // else

    }   // while


    // Set the varianble num_spin_unit equal to the number of spins
    // in the unit cell
    number_of_crystallites = crystallites_count;

    orientations_filestream.close();

}


void Simulation::read_orientations(string orientations_filename) {

    // open an ifstream
    ifstream orientations_filestream(orientations_filename.c_str());

    // Check that the file stream opened normally
    if (orientations_filestream.fail()) {
        cerr << "\nThe file " << orientations_filename
             << " could not be opened in read_orientations().\n";
        exit(1);
    }

    // The file is read line by line a first time to count the number of
    // orientations.  Each line of text is temporarily
    // stored in the string "line."
    string line;

    string token;   // A token to be read from the string "line."


    // The variable orientation_number will be used as an index to keep track
    // of the spin whose coordinates are being read.
    int orientation_number = 0;

    while (!orientations_filestream.eof()) {    // scans the file line by line

        // Reset the strings "line" and "token," for cases where the
        // string operations fail.
        token = "";
        line = "";

        getline(orientations_filestream, line);

        // Creates a stream from the string line_srtream
        istringstream line_stream(line);

        // Read the first token from the line, with tokens separated by whitespace.
        line_stream >> token;

        // If the token contains no characters, skip to the next line
        if (token.length() == 0)
            continue;
        else {

            // Check that there are 4 numbers in the line and nothing else.

            // Return to the beginning of the stream
            line_stream.seekg(0, ios::beg);

            for (int i = 0; i < 4; i++) {
                line_stream >> orientations[orientation_number][i];

                if (line_stream.fail()) {
                    cerr << "\nError reading angles or weighting for a crystallite.\n";
                    exit(1);
                } // check for failure to input a double.

            }      // loop through the 4 numbers on the line.

            orientation_number++;

        }   // else

    }    // while


    // close the input file
    orientations_filestream.close();

} // read orientations


void Simulation::simulate_dynamics(string name,
                                   Crystal &crystal,
                                   Hamiltonian &hamiltonian,
                                   Density_matrix &sigma) {

    // The global indices of the polarised_nuclei have to be calculated
    for (unsigned int i = 0; i < polarised_nuclei.size(); i++) {
        polarised_nuclei[i][4] = crystal.calculate_nucleus_index(polarised_nuclei[i][0], polarised_nuclei[i][1],
                                                                 polarised_nuclei[i][2], polarised_nuclei[i][3]);
    }


    // Loop over the initial density matrices. Each initial
    // density matrix is treated independently
    std::vector<int> polarised_spins;
    string filename_prefix = name;
    for (unsigned int i = 0; i < polarised_nuclei.size(); i++) {

        // The index of the nucleus initially polarised is read
        int polarised_spin = polarised_nuclei[i][4];

        // The names of the output files requires a string that contains the index
        // of the nuclei initially excited
        // The numbering starts at 1 for the user, and at 0 for the program, so
        // 1 has to be added to polarised_spin for the name
        stringstream polarised_spin_buffer;
        polarised_spin_buffer << polarised_spin + 1;
        string polarised_spin_string;
        polarised_spin_string = polarised_spin_buffer.str();

        polarised_spins.push_back(polarised_spin);
        filename_prefix += "_" + polarised_spin_string;
    }

    // Depending on the options for the simulation, a different function is
    // called
    if (number_of_crystallites == 1) // single orientation
        simulate_crystallite(polarised_spins,
                             filename_prefix,
                             crystal,
                             hamiltonian,
                             sigma);
    else   // powder average
        simulate_powder_average(polarised_spins,
                                filename_prefix,
                                crystal,
                                hamiltonian,
                                sigma);


}   // simulate_dynamics


void Simulation::simulate_crystallite(std::vector<int>& polarised_spins,
                                      string filename_prefix,
                                      Crystal &crystal,
                                      Hamiltonian &hamiltonian,
                                      Density_matrix &sigma) {
    return;
    // Create a set of Euler angles corresponding to the only
    // orientation
    Euler euler_C_to_R(orientations[0][0],
                       orientations[0][1],
                       orientations[0][2]);

    // The hamiltonian is updated for this orientation
    hamiltonian.rotate_crystallite(euler_C_to_R, crystal);

    // The density matrix is reset, and initialized by setting the value
    // of z to 1 for the polarised nuclei
    sigma.reset(crystal);
    for (auto ps : polarised_spins)
        sigma.z[ps] = complex1;

    // The result of the time evolution is written to a .dat file
    string dat_filename = filename_prefix + dat_extension + ".dat";
    ofstream dat_filestream(dat_filename.c_str());

    // Check that the file stream opened normally
    if (dat_filestream.fail()) {
        cerr << "\nDat file " << dat_filename << " could not be opened.\n";
        exit(1);
    }

    // Perform the numerical integration of the LVN equation
    // A counter is necessary to know when to write the data. The initial
    // density matrix is always written, so the counter is initialised to
    // its limit value
    int write_counter = sampling_interval;

    for (int step = 0; step < number_of_steps; step++) {

        if (write_counter == sampling_interval) {

            write_counter = 0;

            dat_filestream << scientific << setprecision(6) << setw(14)
                           << step * timestep;

            for (int i = 0; i < crystal.number_of_nuclei; i++) {

                // Individual magnetisation are printed
                dat_filestream << scientific << setprecision(6) << setw(14)
                               << sigma.z[i].real();

            } // i loop
            dat_filestream << "\n";

        } // if write_counter

        write_counter++;

        // If the spinning speed is non-zero, the sample is spun at the magic
        // angle. Otherwise the Rotor and Lab frame stay the same
        if (abs(omega_r) > 0) {
            hamiltonian.spin_at_magic_angle(omega_r * step * timestep, crystal);
        }

        // evolve the density matrix
        sigma.evolve_ST2(hamiltonian, timestep, crystal);

    } // step loop

    dat_filestream.close();

} // simulate_crystallite


void Simulation::simulate_powder_average(std::vector<int>& polarised_spins,
                                         string filename_prefix,
                                         Crystal &crystal,
                                         Hamiltonian &hamiltonian,
                                         Density_matrix &sigma) {

    // In the case where only powder average is calculated, the observables
    // cannot be printed at each step. They have to be stored in an array
    // and printed at the end of the simulation.
    double ***observables;
    double **observables_cumulate;

    double ***br_pm;
    double ***br_mp;
    double ***br_zz;

    complex<double> ****ZQTp;
    complex<double> ***ZQTp_cum;

    int i, j, k; // for loop indices
    int number_threads = omp_get_max_threads();
    printf("%d threads\n", number_threads);
    try {

        observables = new double **[crystal.number_of_nuclei];
        observables_cumulate = new double *[crystal.number_of_nuclei];
#ifdef ZQT
        ZQTp = new complex<double> ***[crystal.number_of_nuclei];
        ZQTp_cum = new complex<double> **[crystal.number_of_nuclei];
#endif
        for (i = 0; i < crystal.number_of_nuclei; i++) {
#ifdef ZQT
            ZQTp[i] = new complex<double> **[crystal.number_of_nuclei];
            ZQTp_cum[i] = new complex<double> *[crystal.number_of_nuclei];
            //for (j = 0; j < crystal.number_of_nuclei; j++) {
            for (auto j : crystal.nearest_neighbours[i]){
                ZQTp[i][j] = new complex<double> *[number_of_steps];
                ZQTp_cum[i][j] = new complex<double> [number_of_steps];
                for (k = 0; k < number_of_steps; k++) {
                    ZQTp[i][j][k] = new complex<double>[number_threads];
                }

            }
#endif

            observables[i] = new double*[number_threads];
            for (j = 0; j < number_threads; j++) {
                observables[i][j] = new double[number_of_steps];
            }

            observables_cumulate[i] = new double[number_of_steps];
        }

        if (tree_pruning) {
            br_pm = new double **[number_threads];
            br_mp = new double **[number_threads];
            br_zz = new double **[number_threads];
            for (i = 0; i < number_threads; i++) {
                br_pm[i] = new double *[crystal.number_of_nuclei];
                br_mp[i] = new double *[crystal.number_of_nuclei];
                br_zz[i] = new double *[crystal.number_of_nuclei];
                for (j = 0; j < crystal.number_of_nuclei; j++) {
                    br_pm[i][j] = new double[j];
                    br_mp[i][j] = new double[j];
                    br_zz[i][j] = new double[j];
                }
            }
        }


    }

    // observables[nuclei][thread][step]
    catch (bad_alloc) {
        cerr << "\nMemory could not be allocated for the observables\n";
        exit(1);
    }
    for (j = 0; j < crystal.number_of_nuclei; j++) {
        for (i = 0; i < number_threads; i++) {
            for (k = 0; k < number_of_steps; k++) {
                observables[j][i][k] = 0.0;
                observables_cumulate[j][k] = 0.0;
            }
        }
    }

    if (tree_pruning) {
        for (i = 0; i < number_threads; i++) {
            for (j = 0; j < crystal.number_of_nuclei; j++) {
                for (k = 0; k < j; k++) {
                    br_pm[i][j][k] = 0.0;
                    br_mp[i][j][k] = 0.0;
                    br_zz[i][j][k] = 0.0;
                }
            }
        }
    }

#ifdef ZQT
    for (i = 0; i < crystal.number_of_nuclei; i++) {
        //for (j = 0; j < crystal.number_of_nuclei; j++) {
        for (auto j : crystal.nearest_neighbours[i]) {
            for (k = 0; k < number_of_steps; k++) {
                ZQTp_cum[i][j][k] = complex0;
                for (int l = 0; l < number_threads; l++) {
                    ZQTp[i][j][k][l] = complex0;
                }
            }

        }
    }
#endif

    // The outer loop runs over the set of orientations. For each
    // orientation, the result of the simulation are stored in the
    // observables array
 #pragma omp parallel for shared(observables, ZQTp, br_pm, br_mp, br_zz)
    for (int crystallite = 0; crystallite < number_of_crystallites;
         crystallite++) {
        int mythread = omp_get_thread_num();
        printf("My thread is %d, crystallite %d -> %d\n", mythread, crystallite, configuration);
        /* need local copies of
         *
                                         Hamiltonian &hamiltonian,
                                         Density_matrix &sigma
                                         */

        //Hamiltonian local_hamiltonian = ;
        Density_matrix local_sigma = Density_matrix(sigma.number_of_spins, crystal);
        string output_filename_hamiltonian = filename_prefix + std::to_string(mythread) + ".ssdo";
        Hamiltonian local_hamiltonian = Hamiltonian(crystal, output_filename_hamiltonian, crystal.S);

        // hamiltonian(crystal, output_filename);
        // Create a set of Euler angles corresponding to the current
        // orientation
        Euler euler_C_to_R(orientations[crystallite][0],
                           orientations[crystallite][1],
                           orientations[crystallite][2]);

        // The hamiltonian is updated for the orientation under study
        local_hamiltonian.rotate_crystallite(euler_C_to_R, crystal);

        // The density matrix is reset, and initialized by setting the value
        // of z to 1 for the polarised nuclei
        local_sigma.reset(crystal);
        int i;
        /* SATURATION
         *  change this to set all to be complex1.
         *  then set sigma.z[polarised_spin] = 0  - eg this is now saturated spin.
         */

        if (configuration == SATURATE || configuration == SATURATEZ) {
            for (i = 0; i < local_sigma.number_of_spins; i++) {
                local_sigma.z[i] = complex1;
            }

            for (auto ps : polarised_spins)
                local_sigma.z[ps] = complex0;
        } else if (configuration == PSD) {
            for (auto ps : polarised_spins)
                local_sigma.z[ps] = complex1;
        } else if (configuration == INVERT) {
            for (i = 0; i < local_sigma.number_of_spins; i++) {
                local_sigma.z[i] = complex1;
            }
            for (auto ps : polarised_spins)
                local_sigma.z[ps] = -complex1;
        }
        if (continuation != 0) {
            local_sigma.read_continuation(filename_prefix + "_" + std::to_string(crystallite) + "_");
        }
        // Perform the numerical integration of the LVN equation
        for (int step = 0; step < number_of_steps; step++) {

            // Individual magnetisation are stored in the observables array,
            // with the appropriate weight obtained from the angle set
            for (i = 0; i < crystal.number_of_nuclei; i++) {
                observables[i][mythread][step] +=
                        local_sigma.z[i].real() * orientations[crystallite][3];
#ifdef ZQT
                
               for (auto j : crystal.nearest_neighbours[i]) {
                    ZQTp[i][j][step][mythread] += local_sigma.pm[i + local_sigma.number_of_spins * j] * orientations[crystallite][3];
                }
#endif
            } // I loop

            // If the spinning speed is non-zero, the sample is spun at the magic
            // angle. Otherwise the Rotor and Lab frame stay the same
            if (abs(omega_r) > 0) {
                local_hamiltonian.spin_at_magic_angle(omega_r * step * timestep, crystal);
            }

            // evolve the density matrix
            local_sigma.evolve_ST1(local_hamiltonian, timestep, crystal);
            local_sigma.evolve_CS(crystal, local_hamiltonian, timestep);
            if (configuration == SATURATE)
                local_sigma.saturate(polarised_spins);
            else if (configuration == SATURATEZ) {
                for (auto spin: polarised_spins)
                    local_sigma.z[spin] = complex0;
            }

            //if (step % 1000)
            //    local_sigma.write_continuation(filename_prefix + "_" + std::to_string(crystallite) + "_");
            //printf("%d elements in map\n", local_sigma.ppmm.size());
        } // step loop

        if (tree_pruning)
            local_sigma.store_br(br_zz[mythread], br_pm[mythread], br_mp[mythread], crystal);

    } // orientation loop
    printf("I finish the loop\n");


    // bring all the observables together
    for (i = 0 ; i< crystal.number_of_nuclei; i++) {
        for (j = 0; j < number_of_steps; j++) {
            for (k = 0; k < number_threads; k++) {
                observables_cumulate[i][j] += observables[i][k][j];
                //for (int l = 0; l < crystal.number_of_nuclei; l++) {
#ifdef ZQT
                for (auto l : crystal.nearest_neighbours[i]) {
                    ZQTp_cum[i][l][j] += ZQTp[i][l][j][k];
                }
#endif

            }
        }
    }

    // The result of the time evolution is written to a .dat file
    string dat_filename = filename_prefix + dat_extension + ".dat";
    ofstream dat_filestream(dat_filename.c_str(), (continuation != 0) ? ios::app : ios::out);


    // Check that the file stream opened normally
    if (dat_filestream.fail()) {
        cerr << "\nDat file " << dat_filename << " could not be opened.\n";
        exit(1);
    }

    // A counter is necessary to knwow when to write the data. The initial
    // density matrix is always written, so the counter is initialise to
    // its limit value
    int write_counter = sampling_interval;
    for (int step = 0; step < number_of_steps; step++) {
        if (write_counter == sampling_interval) {
            write_counter = 0;
            dat_filestream << scientific << setprecision(6) << setw(14)
                           << step * timestep;
            for (int nuclei = 0; nuclei < crystal.number_of_nuclei; nuclei++) {
                dat_filestream << scientific << setprecision(6) << setw(14)
                               << observables_cumulate[nuclei][step];

            } // i loop
            dat_filestream << "\n";

        }    // if write_counter
        write_counter++;
    } // step loop
    dat_filestream.close();

    if (tree_pruning) {
        /* if tree pruning, then output a file with all of the br values (to determine which spins to include!) */
        string tp_filename = filename_prefix + ".tp";
        ofstream tp_filestream(tp_filename.c_str());
        for (j = 0; j < crystal.number_of_nuclei; j++) {
            for (k = 0; k < j; k++) {
                double zz = 0, pm = 0, mp = 0;
                for (i = 0; i < number_threads; i++) {
                    zz += br_zz[i][j][k];
                    pm += br_pm[i][j][k];
                    mp += br_mp[i][j][k];
                }
                tp_filestream << setw(8) << j;
                tp_filestream << setw(8) << k;
                tp_filestream << scientific << setprecision(6) << setw(14)
                              << zz;
                tp_filestream << scientific << setprecision(6) << setw(14)
                              << pm;
                tp_filestream << scientific << setprecision(6) << setw(14)
                              << mp;
                tp_filestream << endl;
            }
        }
    }

#ifdef ZQT
    /* If in ZQT mode, output the I+I- FID */
    for (int nuclei = 0 ; nuclei < crystal.number_of_nuclei; nuclei++) {
        string zqt_filename = filename_prefix + "_" + std::to_string(nuclei) + ".zqt";
        cout << zqt_filename << endl;
        ofstream zqt_filestream(zqt_filename.c_str());
        if (zqt_filestream.fail()) {
            cerr << "\nZQT file " << zqt_filename << " could not be opened.\n";
            exit(1);
        }
        write_counter = sampling_interval;
        zqt_filestream << setw(14) << "# timestep";
        for (auto nuclei2 : crystal.nearest_neighbours[nuclei]) {
            zqt_filestream << setw(14) << std::to_string(nuclei2) + "_real";
            zqt_filestream << setw(14) << std::to_string(nuclei2) + "_imag";
        }
        zqt_filestream << "\n";

        for (int step = 0; step < number_of_steps; step++) {
            if (write_counter == sampling_interval) {
                write_counter = 0;

                zqt_filestream << scientific << setprecision(6) << setw(14) << step * timestep;
                for (auto nuclei2 : crystal.nearest_neighbours[nuclei]) {
                    zqt_filestream << scientific << setprecision(6) << setw(14)
                                   << ZQTp_cum[nuclei][nuclei2][step].real();
                    zqt_filestream << scientific << setprecision(6) << setw(14)
                                   << ZQTp_cum[nuclei][nuclei2][step].imag();
                }
                zqt_filestream << "\n";

            }    // if write_counter
            write_counter++;
        } // step loop

        zqt_filestream.close();
    }
#endif
}    // simulate_powder_average
