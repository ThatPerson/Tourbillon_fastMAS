#include "constants.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <set>
#include<bits/stdc++.h>

#define NMAX 16 
// maximum number of neighbours in simulation (eg Nmax in the paper).

using namespace std;


// The Crystal class contains geometrical information about a system 
// of interest. Functions are available to read the size and content of a 
// unit cell from an input file, and to build a supercell. The constructor 
// needs to be given the name of a properly formatted input file, along with 
// the name of an output file where some of the information is written.
class Crystal {

public:

    int number_of_nuclei;        // Number of nuclei in the crystal

    double Rz, Rpm;

    // Keyword for the use of periodic boundary conditions. By default PBC are
    // turned off
    bool periodic_boundary_conditions;

    double ***nuclei_separation;    // Array of vector separations between pairs of
    // nuclei, given in spherical coordinates. The first dimension corresponds to
    // rho, theta and phi for 0, 1 and 2 respectively. The second and third
    // dimensions are the indices of the nuclei
    // By convention, the cell vector a is along x, b is in the xy plane, and
    // (a, b, c) is direct.

    double *csiso_ppm;
    double *csiso_ppm_unit;
    double *csaniso_ppm, *csaniso_ppm_unit;
    double *cseta, *cseta_unit;
    double *csalpha, *csalpha_unit, *csbeta, *csbeta_unit, *csgamma, *csgamma_unit, *csavg, *csavg_unit;
    double methyl_rotation_val;
    double field_strength;
    // Constructor function
    Crystal(std::string input_filename, std::string output_filename);

    // Function that calculates the index of a nucleus in the full spin system,
    // given the position of the unit cell replica in which the spin is, and the
    // index of the spin in the initial unit cell
    int calculate_nucleus_index(int cell_index_a,
                                int cell_index_b,
                                int cell_index_c,
                                int nucleus_index_in_unit_cell);

    std::set<int> *nearest_neighbours;
    // Set containing the nearest neighbour spins to each spin.
    double S;


private:
    int get_nearest(int i, std::set<int> &results);
    int read_nearest(std::set<int> *results);
    // Informations about the unit cell: its dimensions, the number of nuclei it
    // contains, and the fractional cartesian coordinates of the nuclei
    double unit_cell_len[3];    // a, b and c length
    bool nearest_CS;

    string nn_file;

    double cell_ang[3];    // alpha, beta and gamma angles

    int number_of_nuclei_unit;    // Number of nuclei in unit cell

    double **coord_unit;        // 2D array containing the fractional
    // coordinates of the nuclei in the unit cell. coord_unit[i][j] stores the
    // Ith coordinate of the Jth nucleus, where a, b, c, correspond to I=0, 1, 2,
    // respectively.

    // Informations about the supercell: its dimensions, the number of nuclei it
    // contains, and the coordinates of the nuclei
    int number_unit_cells[3];        // Number of cells along a, b and c
    // By default, the simulation cell and the unit cell are identical

    double cell_len[3];    // Lattice length along a, b and c

    double **coord;        // 2D array storing the fractional coordinates,
    // defined with respect to the supercell vectors


    double **cart_coord;        // 2D array storing Cartesian coordinates,
    // only for printing in the output file

    double **H;            // Change-of-basis matrix from
    // the Cartesian basis to the (a, b, c) basis. It contains in its ith column
    // the Cartesian coordinates of the ith cell vector


    // Functions called by the constructor.
    void read_number_of_nuclei(std::string input_filename);

    void read_crystal_parameters(std::string input_filename);

    void write_crystal_parameters(std::string output_filename);

    void calculate_H();

    void generate_supercell();

    void calculate_nuclei_separation();


    // The following function is called by the constructor only when
    // periodic boundary conditions are used.
    // Function that wraps the simulation cell, by selecting, for each nucleus,
    // the periodic image that has fractional coordinates between 0 and 1.
    void wrap_cell();

};


/* Identifies the nearest spins to a given spin, i */
int Crystal::get_nearest(int i, std::set<int> &results) {
    std::vector<std::pair <double, int>> distances;
    /* pair, double is the distance between spins, int is the spin index */
    int j;
    double xsep, ysep, zsep;
    double dist;
    for (j = 0; j < number_of_nuclei; j++) {
        if (i == j)
            continue;

        dist = nearest_CS ? fabs(csiso_ppm[i] - csiso_ppm[j]) : nuclei_separation[0][i][j];

        distances.push_back(make_pair(dist, j));
    }

    sort(distances.begin(), distances.end());

    int nmax = (NMAX > distances.size()) ? distances.size() : NMAX;
    // if we have fewer spins than we have NMAX, then use the number of spins
    printf("%d\n", nmax);
    for (j = 0; j < nmax; j++) {
        if (distances[j].second > i) 
            /* greater than is as in Perras, less than is as in Dumez.
               don't want to double count - have to make a decision to either keep the pairs
               so that the neighbours all have greater index (e.g., a 3 - 4 pair would have {4}
               in the neighbour set of 3 but {3} would not be in the neighbour set of 4).*/
            results.insert(distances[j].second);

    }
    printf("%d: {", i+1);
    for (auto m : results) {
        printf("%d, ", m+1);
    }
    printf("}\n");

    return 1;
}

int Crystal::read_nearest(std::set<int> *results) {
    ifstream input_filestream(nn_file.c_str());
    if (input_filestream.fail()) {
        cerr << "\nInput file " << nn_file
             << " could not be opened in read_nearest().\n";
        exit(1);
    }

    string line, token;
    int i, j;
    while (!input_filestream.eof()) {
        token = "";
        line = "";
        getline(input_filestream, line);
        cout << line << endl;
        istringstream line_stream(line);
        if (line.length() == 0)
            continue;
        line_stream >> i;
        line_stream >> j;
        if (j > i) {
            if (i < number_of_nuclei && j < number_of_nuclei) {
                //if ((csiso_ppm[i] - csiso_ppm[j]) > 0.001) // if the spins aren't magnetically equivalent
                    results[i].insert(j);

            } else {
                cerr << "Unknown spins " << i << " and " << j << endl;
            }
        } else {
            cerr << "Misordered spins " << i << " and " << j << endl;
        }
        /*cout << line << endl;
        if (token.length() == 0)
            continue;
        try {
            i = stoi(token);
            line_stream >> token;
            if (token.length() == 0)
                continue;
            j = stoi(token);
            if (j > i) {
                if (i < number_of_nuclei && j < number_of_nuclei) {
                    results[i].insert(j);
                } else
                    cout << "Unknown spins " << i << " and " << j << endl;
            } else
                cout << "Malformed NN file" << endl;
        } catch (...) {
            cout << "Could not parse" << endl << "\t" << line << endl;
            break; // not sure why it keeps going at EOF
        }*/

    }
    input_filestream.close();
  //  exit(-1);
    return 1;

}
// If a filename is passed to the constructor function, it first reads the
// number of nuclei and allocates memory for the arrays coord_unit, coord,
// nuclei_separation. It then calls successively four functions:
// read_crystal_parameters, which reads the coordinates of the nuclei
// in the unit cell; generate_supercell, which builds the
// supercell by periodically repeating the unit cell, and fills the array
// coord; calculate_nuclei_separation, which calculates the
// vector separation between each pair of nuclei, 
// optionally using periodic boundary conditions 
// and a minimum image convention; write_crystal_parameters, that appends some
// of this information to the output file.
// The objet thus constructed is a crystal, ready for use to initialise
// a hamiltonian.
Crystal::Crystal(string input_filename, string output_filename) {
    nearest_CS = false;
    read_number_of_nuclei(input_filename);
    nn_file = "";
    S = 1; // dipolar coupling averaging
    Rz = 0;
    Rpm = 0;
    methyl_rotation_val = 1;
    try {

        // Allocates memory for a 3*number_of_nuclei_unit 2D array
        coord_unit = new double *[3];
        for (int i = 0; i < 3; i++)
            coord_unit[i] = new double[number_of_nuclei_unit];
        csiso_ppm_unit = new double[number_of_nuclei_unit];
        csaniso_ppm_unit = new double[number_of_nuclei_unit];
        cseta_unit = new double[number_of_nuclei_unit];
        csalpha_unit = new double[number_of_nuclei_unit];
        csbeta_unit = new double[number_of_nuclei_unit];
        csgamma_unit = new double[number_of_nuclei_unit];
        csavg_unit = new double[number_of_nuclei_unit]; // unused
    }
    catch (bad_alloc) {
        cerr << "\nMemory could not be allocated for geometrical information about the unit cell\n";
        exit(1);
    }


    // Initialise some arrays to 0
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < number_of_nuclei_unit; j++) {
            coord_unit[i][j] = 0.0;
        }
    }


    // Set values of parameters which will be read to impossible values
    // for error-checking
    for (int i = 0; i < 3; i++) {
        unit_cell_len[i] = -1.0;
        cell_ang[i] = -1.0;
    }

    // By default, the simulation cell and the unit cell are identical
    for (int i = 0; i < 3; i++) {
        number_unit_cells[i] = 1;
    }

    // Initialise the keywords to their default value
    periodic_boundary_conditions = false;
    field_strength = 700;
    // the function read_crystal_parameters reads the parameter file,
    // obtaining values for unit_cell_len[], coord_unit[][],
    // number_unit_cells[]
    read_crystal_parameters(input_filename);

    // check to see if values were read into unit_cell_len
    for (int i = 0; i < 3; i++) {
        if ((unit_cell_len[i] < 0) || (cell_ang[i] < 0)) {
            cerr << "\nError in reading values of cell dimensions." << endl;
            exit(1);
        } // if
    }   // for

    // calculate the number of nuclei
    number_of_nuclei = number_of_nuclei_unit * number_unit_cells[0] * number_unit_cells[1] * number_unit_cells[2];

    // Allocate memory for the cartesian coordinates of the nuclei the
    // distances between the nuclei expressed in spherical coordinates,
    // and the coupling constants.
    try {

        coord = new double *[3];
        cart_coord = new double *[3];
        for (int i = 0; i < 3; i++) {
            coord[i] = new double[number_of_nuclei];
            cart_coord[i] = new double[number_of_nuclei];
        }
        csiso_ppm = new double[number_of_nuclei];
        csaniso_ppm = new double[number_of_nuclei];
        cseta = new double[number_of_nuclei];
        csalpha = new double[number_of_nuclei];
        csbeta = new double[number_of_nuclei];
        csgamma = new double[number_of_nuclei];
        nuclei_separation = new double **[3];
        for (int i = 0; i < 3; i++) {
            nuclei_separation[i] = new double *[number_of_nuclei];
            for (int j = 0; j < number_of_nuclei; j++)
                nuclei_separation[i][j] = new double[number_of_nuclei];
        }

        H = new double *[3];
        for (int i = 0; i < 3; i++)
            H[i] = new double[3];

    }
    catch (bad_alloc) {
        cerr << "\nMemory could not be allocated for geometrical information about the supercell\n";
        exit(1);
    }


    // all the array elements must be set to zero
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < number_of_nuclei; j++) {
            coord[i][j] = 0.0;
            cart_coord[i][j] = 0.0;
            for (int k = 0; k < number_of_nuclei; k++) {
                nuclei_separation[i][j][k] = 0.0;
            }    // k loop
        }    // j loop
    }    // i loop

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            H[i][j] = 0.0;
        }
    }

    // calculate and store all nuclei coordinates
    generate_supercell();

    // calculate the matrix for the change of basis
    calculate_H();

    // If periodic boundary conditions are used, the system first has to be
    // wrapped inside the simulation cell
    if (periodic_boundary_conditions == true)
        wrap_cell();

    // calculate and store all coupling constants
    calculate_nuclei_separation();

    nearest_neighbours = new std::set<int>[number_of_nuclei];
    if (nn_file != "") {
        read_nearest(nearest_neighbours);
    } else {
        for (int i = 0; i < number_of_nuclei; i++) {
            get_nearest(i, nearest_neighbours[i]);
        }
    }
    // append some of the parameters to the output file
    write_crystal_parameters(output_filename);

} // constructor that receives filename


int Crystal::calculate_nucleus_index(int cell_index_a,
                                     int cell_index_b,
                                     int cell_index_c,
                                     int nucleus_index_in_unit_cell) {
    int nucleus_index_in_supercell;

    // The convention for indexing a nucleus in the supercell is as follows:
    //  index = index_unit + num_nuclei_xy * k
    //                     + num_nuclei_x * j
    //                     + number_of_nuclei_unit * i
    // where index_unit is the index of the nuclei in the original unit cell,
    // num_nuclei_xy is the number of nuclei in a set of unit cells with
    // identical z, and num_nuclei_x is the number of nuclei in a set of
    // unit cells with identical z and y.

    nucleus_index_in_supercell = nucleus_index_in_unit_cell - 1
                                 + number_of_nuclei_unit * number_unit_cells[0] * number_unit_cells[1] *
                                   (cell_index_c - 1)
                                 + number_of_nuclei_unit * number_unit_cells[0] * (cell_index_b - 1)
                                 + number_of_nuclei_unit * (cell_index_a - 1);

    // Verify that the indices in the correct range of values.
    if (nucleus_index_in_supercell < 0 || nucleus_index_in_supercell >= number_of_nuclei) {
        cerr << "\nError: A specified nucleus index is out of range."
             << endl;
        exit(1);
    }

    return (nucleus_index_in_supercell);

} // calculate_index


void Crystal::read_number_of_nuclei(string input_filename) {

    // open an ifstream
    ifstream input_filestream(input_filename.c_str());

    // Check that the file stream opened normally
    if (input_filestream.fail()) {
        cerr << "\nInput file " << input_filename
             << " could not be opened in read_crystal_parameters().\n";
        exit(1);
    }

    // The file is read line by line.  Each line of text is temporarily
    // stored in the string "line."
    string line;

    string token;   // A token to be read from the string "line."

    int spin_count = 0;


    while (!input_filestream.eof()) {

        // Reset the strings "line" and "token," for cases where the
        // string operations fail.
        token = "";
        line = "";

        getline(input_filestream, line);

        // creates a stream from the string filestream
        istringstream line_stream(line);

        // Read the first token from the line, with tokens separated by whitespace.
        line_stream >> token;

        // Check to see whether the token begins with 'H'

        // If the token contains no characters, skip to the next line
        if (token.length() == 0)
            continue;

        else if (token[0] == 'H')
            spin_count++;    // counts the number of proton in the unit cell
    }   // while


    // Set the varianble num_spin_unit equal to the number of spins
    // in the unit cell
    number_of_nuclei_unit = spin_count;

    // close the input file
    input_filestream.close();

} // read_number_of_nuclei


// The function read_crystal_parameters() reads the crystal parameters from
// a file:  unit_cell_len[], cell_ang[], coord_unit[][], number_cell[].
// This function assumes that the name of the file to read is in the
// member variable filename.
void Crystal::read_crystal_parameters(string input_filename) {

    // open an ifstream
    ifstream input_filestream(input_filename.c_str());

    // Check that the file stream opened normally
    if (input_filestream.fail()) {
        cerr << "\nInput file " << input_filename
             << " could not be opened in read_crystal_parameters().\n";
        exit(1);
    }

    // The file is read line by line.  Each line of text is temporarily
    // stored in the string "line."
    string line;

    string token;   // A token to be read from the string "line."


    // The variable spin_number will be used as an index to keep track
    // of the spin whose coordinates are being read.
    int spin_number = 0;

    while (!input_filestream.eof()) {    // scans the file line by line

        // Reset the strings "line" and "token," for cases where the
        // string operations fail.
        token = "";
        line = "";

        getline(input_filestream, line);

        // Creates a stream from the string line_srtream
        istringstream line_stream(line);

        // Read the first token from the line, with tokens separated by whitespace.
        line_stream >> token;

        // If the token contains no characters, skip to the next line
        if (token.length() == 0)
            continue;

        else if (token.substr(0, 13) == "_cell_length_") {

            // If token begins with "_cell_length_" then read cell dimensions
            switch (token[13]) {
                case 'a':
                    line_stream >> unit_cell_len[0];
                    break;
                case 'b':
                    line_stream >> unit_cell_len[1];
                    break;
                case 'c':
                    line_stream >> unit_cell_len[2];
                    break;
                default:
                    cerr << "\n Error in reading cell lengths" << endl;
                    exit(1);
            }    // switch
            if (line_stream.fail()) {
                cerr << "\nError in reading cell angles\n";
                exit(1);
            }
        } else if (token.substr(0, 12) == "_cell_angle_") {

            // If token begins with "_cell_angle_" then read cell dimensions
            switch (token[12]) {
                case 'a':
                    line_stream >> cell_ang[0];
                    break;
                case 'b':
                    line_stream >> cell_ang[1];
                    break;
                case 'g':
                    line_stream >> cell_ang[2];
                    break;
                default:
                    cerr << "\n Error in reading cell angles" << endl;
                    exit(1);
            }    // switch
            if (line_stream.fail()) {
                cerr << "\nError in reading cell lenghts\n";
                exit(1);
            }
        } else if (token.substr(0, 17) == "_number_of_cells_") {

            // If the token begins with "_number_of_cells_" then read the number of
            // cells along a given translation vector.
            switch (token[17]) {
                case 'a':
                    line_stream >> number_unit_cells[0];
                    break;
                case 'b':
                    line_stream >> number_unit_cells[1];
                    break;
                case 'c':
                    line_stream >> number_unit_cells[2];
                    break;
                default:
                    cerr << "\n Error in reading number of unit cells" << endl;
                    exit(1);
            }    // switch
            if (line_stream.fail()) {
                cerr << "\nError in reading  number of unit cells\n";
                exit(1);
            }
        } else if (token[0] == 'H') {
            for (int i = 0; i < 3; i++) {
                line_stream >> coord_unit[i][spin_number];
            }
            line_stream >> csiso_ppm_unit[spin_number];
            line_stream >> csaniso_ppm_unit[spin_number];
            line_stream >> cseta_unit[spin_number];
            line_stream >> csalpha_unit[spin_number];
            line_stream >> csbeta_unit[spin_number];
            line_stream >> csgamma_unit[spin_number];
            line_stream >> csavg_unit[spin_number];
            if (line_stream.fail()) {
                cerr << "\nError in reading atomic information\n";
                exit(1);
            }
            spin_number++;
        } else if (token == "_nearest_CS") {
            nearest_CS = true;// if
        }
        else if (token == "_field_strength") {
            line_stream >> field_strength;
        }
        else if (token == "_nearest_neighbour_file") {
            line_stream >> nn_file;
        }
        else if (token == "_methyl_rotation") {
            line_stream >> methyl_rotation_val;
        }
        else if (token == "_periodic_boundary_conditions") {
            periodic_boundary_conditions = true;
        }   // if
        else if (token == "_Rz") {
            line_stream >> Rz;
	    Rz /= 1000.; // Rz is given in s-1, but we want it in ms-1.
        }
        else if (token == "_Rpm") {
            line_stream >> Rpm;
	    Rpm /= 1000.0;
        }
        else if (token == "_dipolar_averaging") {
            line_stream >> S;
        }

        else
            continue;            // none of the expected token were found

    }    // while


    // close the input file
    input_filestream.close();

} // function read_crystal_parameters


void Crystal::write_crystal_parameters(string output_filename) {

    // Open an ofstream in append mode
    ofstream output_filestream(output_filename.c_str(), ios_base::app);

    // Check that the file stream opened normally
    if (output_filestream.fail()) {
        cerr << "\nOutput file " << output_filename << " could not be opened.\n";
        exit(1);
    }

    // Lattice parameters read from the input file
    output_filestream << "_cell_length_a\t\t\t\t"
                      << fixed << setprecision(3)
                      << unit_cell_len[0]
                      << endl
                      << "_cell_length_b\t\t\t\t"
                      << fixed << setprecision(3)
                      << unit_cell_len[1]
                      << endl
                      << "_cell_length_c\t\t\t\t"
                      << fixed << setprecision(3)
                      << unit_cell_len[2]
                      << endl << endl;

    // Lattice parameters read from the input file
    output_filestream << "_cell_angles_alpha\t\t\t"
                      << fixed << setprecision(3)
                      << cell_ang[0]
                      << endl
                      << "_cell_angles_beta\t\t\t"
                      << fixed << setprecision(3)
                      << cell_ang[1]
                      << endl
                      << "_cell_angles_gamma\t\t\t"
                      << fixed << setprecision(3)
                      << cell_ang[2]
                      << endl << endl;

    // number of unit cells in the simulation cell, read from the input file
    output_filestream << "_number_of_cells_a_direction\t\t"
                      << number_unit_cells[0] << endl
                      << "_number_of_cells_b_direction\t\t"
                      << number_unit_cells[1] << endl
                      << "_number_of_cells_c_direction\t\t"
                      << number_unit_cells[2] << endl
                      << endl;

    output_filestream << "_field" << "\t\t\t" << field_strength << endl;

    output_filestream << "_Rz\t\t\t" << Rz << endl;
    output_filestream << "_Rpm\t\t\t" << Rpm << endl;



    // Number of nuclei in the unit cell, as read from the input file, and
    // number of nuclei in the simulation cell
    output_filestream << "_number_of_nuclei_in_unit_cell\t\t"
                      << number_of_nuclei_unit << endl
                      << "_number_of_nuclei\t\t\t"
                      << number_of_nuclei << endl
                      << endl;

    // Keyword for PBC
    if (periodic_boundary_conditions == true) {
        output_filestream << "_periodic_boundary_conditions\t\tyes" << endl
                          << endl;
    } else {
        output_filestream << "_periodic_boundary_conditions\t\tno" << endl
                          << endl;
    }

    // Cartesian coordinates of the nuclei, which have to be calulated first

    for (int i = 0; i < number_of_nuclei; i++) {
        for (int axis = 0; axis < 3; axis++) {
            cart_coord[axis][i] = coord[0][i] * H[axis][0]
                                  + coord[1][i] * H[axis][1]
                                  + coord[2][i] * H[axis][2];
        }
    }

    output_filestream << "nearest_CS" << (nearest_CS) ? "true" : "false";
    output_filestream << "Cartesian coordinates of the nuclei\n";
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < number_of_nuclei; j++) {
            output_filestream << setw(6) << fixed << setprecision(2) <<
                              cart_coord[i][j] << " ";
        }
        output_filestream << "\n";
    }
    output_filestream << "\n";

    // Rho for each nuclei pair, as calculated by Crystal
    output_filestream << "nuclei_separation: rho, in angstrom\n";
    for (int i = 0; i < number_of_nuclei; i++) {
        for (int j = 0; j < number_of_nuclei; j++) {
            output_filestream << setw(6) << fixed << setprecision(2) <<
                              nuclei_separation[0][i][j] << " ";
        }
        output_filestream << "\n";
    }
    output_filestream << "\n";

    // Theta for each nuclei pair, as calculated by Crystal
    output_filestream << "nuclei_separation: theta, in degrees\n";
    for (int i = 0; i < number_of_nuclei; i++) {
        for (int j = 0; j < number_of_nuclei; j++) {
            output_filestream << setw(7) << fixed << setprecision(2) <<
                              nuclei_separation[1][i][j] / PI * 180 << " ";
        }
        output_filestream << "\n";
    }
    output_filestream << "\n";

    // Phi for each nuclei pair, as calculated by Crystal
    output_filestream << "nuclei_separation: phi, in degrees\n";
    for (int i = 0; i < number_of_nuclei; i++) {
        for (int j = 0; j < number_of_nuclei; j++) {
            output_filestream << setw(7) << fixed << setprecision(2) <<
                              nuclei_separation[2][i][j] / PI * 180 << " ";
        }
        output_filestream << "\n";
    }
    output_filestream << "\n";

    output_filestream << "nearest_neighbours\n";
    for (int i = 0; i < number_of_nuclei; i++) {
        output_filestream << setw(5) << fixed << i << "\t\t{";
        for (auto m : nearest_neighbours[i]) {
            output_filestream << m << " ";
        }
        output_filestream << "}" << endl;
    }
    output_filestream << endl;

    output_filestream << "dipolar_averaging\n";
    output_filestream << setw(5) << fixed << S << endl;

    output_filestream.close();

}


// The function generate_supercell requires the variables number_unit_cells[], 
// number_of_nuclei_unit, and coord_unit[][] to have been previously
// initialised. It then generates an array coord[][] of fractional
// coordinates for all nuclei in the supercell.
// Memory must have been allocated for coord.
// Once the function has been called, the array coord contains 
// the fractional coordinates of the nuclei for the whole supercell (i.e. the
// whole set of periodically repeated unit cells).
void Crystal::generate_supercell() {

    // In order to calculate the index of a nuclei in the crystal, it
    // is useful to know the number of nuclei in a set of unit cells with
    // identical z, here named num_nuclei_xy, and the number of nuclei in a
    // set of unit cells with identical z and y, here named num_nuclei_x.
    int num_nuclei_xy = number_of_nuclei_unit * number_unit_cells[0] * number_unit_cells[1];
    int num_nuclei_x = number_of_nuclei_unit * number_unit_cells[0];


    // The crystal, which consists of the periodic repetition of a unit cell, is
    // explored unit cell by unit cell. Within each repetition of the unit cell
    // an index is given to each nuclei:
    //  index = index_unit + num_nuclei_xy * k
    //                     + num_nuclei_x * j
    //                     + number_of_nuclei_unit * i
    // where index_unit is the index of the nuclei in the original unit cell

    // Loop over all the cells in the supercell
    for (int k = 0; k < number_unit_cells[2]; k++) {
        for (int j = 0; j < number_unit_cells[1]; j++) {
            for (int i = 0; i < number_unit_cells[0]; i++) {

                // The indices of the nuclei in the image cell are shifted
                int shift = num_nuclei_xy * k + num_nuclei_x * j
                            + number_of_nuclei_unit * i;

                for (int index_unit = 0; index_unit < number_of_nuclei_unit; index_unit++) {

                    // The index is obtained by a simple shift
                    int index = index_unit + shift;

                    // the fractional coordinates with respect to the unit cell vectors
                    // are obtained by adding an integer
                    coord[0][index] = coord_unit[0][index_unit] + i;
                    coord[1][index] = coord_unit[1][index_unit] + j;
                    coord[2][index] = coord_unit[2][index_unit] + k;
                    csiso_ppm[index] = csiso_ppm_unit[index_unit];
                    csaniso_ppm[index] = csaniso_ppm_unit[index_unit];
                    cseta[index] = cseta_unit[index_unit];
                    csalpha[index] = csalpha_unit[index_unit];
                    csbeta[index] = csbeta_unit[index_unit];
                    csgamma[index] = csgamma_unit[index_unit];
                }


            }      // i loop
        }      // j loop
    }      // k loop


    // Finally, calculate the dimension of the supercell
    for (int i = 0; i < 3; i++)
        cell_len[i] = number_unit_cells[i] * unit_cell_len[i];


    // The fractional coordinates with respect to the supercell vectors are
    // obtained by dividing the current values by the number of unit cells
    // in the supercell
    for (int nucleus = 0; nucleus < number_of_nuclei; nucleus++) {
        for (int cell_vector = 0; cell_vector < 3; cell_vector++) {
            coord[cell_vector][nucleus] = coord[cell_vector][nucleus]
                                          / (number_unit_cells[cell_vector] * 1.0);
        }
    }


} // generate_supercell


void Crystal::calculate_H() {

    // The change-of-basis matrix from Cartesian to cell basis consists of
    // the decomposition of the cell vectors in the Cartesian basis.
    // A convention has to be chosen for the orientation of the cell vector
    // with respect to the Cartesian vector:
    // . a is along Ox
    // . b is in the xOy plane
    // A third constraint is needed (there are three angles and three length for
    // the cell vectors, and H has nine elements); the positive root squared is
    // taken for the z component of a

    // A few quantities are needed to calculate the element of H:
    double cosalpha = cos(cell_ang[0] * PI / 180.0);
    double cosbeta = cos(cell_ang[1] * PI / 180.0);
    double cosgamma = cos(cell_ang[2] * PI / 180.0);
    double singamma = sin(cell_ang[2] * PI / 180.0);

    // Components of a
    H[0][0] = cell_len[0];
    H[1][0] = 0.0;
    H[2][0] = 0.0;

    // Components of b
    H[0][1] = cell_len[1] * cosgamma;
    H[1][1] = cell_len[1] * singamma;
    H[2][1] = 0.0;

    // Components of c
    H[0][2] = cell_len[2] * cosbeta;
    H[1][2] = cell_len[2] * (cosalpha - cosbeta * cosgamma) / singamma;
    H[2][2] = sqrt(cell_len[2] * cell_len[2]
                   - H[0][2] * H[0][2]
                   - H[1][2] * H[1][2]);


} // Calculate H


// The function calculate_nuclei_separation assumes that the array of
// nuclei coordinates coord has been generated. It calculates and stores
// in spherical coordinates all vectors from a nucleus i to a nucleus j.
// If periodic boundary conditions are used, for each nucleus i the 
// internuclear vector between i and each other nucleus j is calculated by 
// considering the periodic image of j that is inside a unit cell centred 
// on i.
void Crystal::calculate_nuclei_separation() {

    // loop over all nuclei i, and all nuclei j with j > i
    for (int i = 0; i < number_of_nuclei; i++) {
        for (int j = i + 1; j < number_of_nuclei; j++) {

            // The variable delta contains the vector between nucleus i and
            // the periodic image of j which is the closest, expressed
            // in the cell basis:
            // { xj-xi; yj-yi; zj-zi }
            double delta_frac[3] = {0.0, 0.0, 0.0};

            // To start with, consider the distance inside the initial system
            for (int cell_vector = 0; cell_vector < 3; cell_vector++) {
                delta_frac[cell_vector] = coord[cell_vector][j]
                                          - coord[cell_vector][i];
            }


            // Under PBC with the minimum image convention,
            // the distance is calculated
            // between i and the periodic image of j that is inside a unit cell
            // centered on i. This image is found by considering the three
            // components of the internuclear vector in the cell basis. For each
            // component, if the value is >= 0.5 or < -0.5, j is outside
            // the unit cell centered on i and has to be replaced by one
            // of its periodic image.
            if (periodic_boundary_conditions == true) {

                for (int cell_vector = 0; cell_vector < 3; cell_vector++) {

                    if (coord[cell_vector][j] - coord[cell_vector][i] >= 0.5)
                        delta_frac[cell_vector] = (coord[cell_vector][j] - 1)
                                                  - coord[cell_vector][i];

                    else if (coord[cell_vector][j] - coord[cell_vector][i] < -0.5)
                        delta_frac[cell_vector] = (coord[cell_vector][j] + 1)
                                                  - coord[cell_vector][i];
                }
            }


            // Now, the internuclear vector has to be converted from fractional
            // to cartesian coordinates, in order to calculate the spherical
            // coordinates
            double delta[3] = {0.0, 0.0, 0.0};

            for (int axis = 0; axis < 3; axis++) {
                delta[axis] = delta_frac[0] * H[axis][0]
                              + delta_frac[1] * H[axis][1]
                              + delta_frac[2] * H[axis][2];
            }


            double rho = sqrt(delta[0] * delta[0]
                              + delta[1] * delta[1]
                              + delta[2] * delta[2]);

            nuclei_separation[0][i][j] = rho;
            nuclei_separation[1][i][j] = acos(delta[2] / rho);
            nuclei_separation[2][i][j] = atan2(delta[1], delta[0]);

            // In spherical coordinates, phi must belong to [ 0 , TWO_PI ]
            if (nuclei_separation[2][i][j] < 0)
                nuclei_separation[2][i][j] += TWO_PI;


            // The lower triangle matrix, where i > j, can be obtained by calculating
            // the components of Rji, in spherical coordinates, from those of Rij
            nuclei_separation[0][j][i] = nuclei_separation[0][i][j];
            nuclei_separation[1][j][i] = PI - nuclei_separation[1][i][j];

            // In spherical coordinates, phi must belong to [ 0 , TWO_PI ]
            if (nuclei_separation[2][i][j] >= PI)
                nuclei_separation[2][j][i] = nuclei_separation[2][i][j] - PI;
            else
                nuclei_separation[2][j][i] = nuclei_separation[2][i][j] + PI;

        } // j loop
    }   // i loop

} // calculate_nuclei_separation


// This function wraps the system inside the simulation cell, by selecting,
// for each nucleus, the periodic image that has fractional coordinates 
// between 0 and 1
void Crystal::wrap_cell() {

    for (int nucleus = 0; nucleus < number_of_nuclei; nucleus++) {
        for (int cell_vector = 0; cell_vector < 3; cell_vector++) {

            double unwrap_coord = coord[cell_vector][nucleus];

            // If the fractional coordinate is larger than 1 or smaller than 0,
            // it should be replaced by its initial value minus the largest integer
            // that is not larger
            if ((unwrap_coord >= 1) || (unwrap_coord < 0))
                coord[cell_vector][nucleus] = unwrap_coord - floor(unwrap_coord);

        }
    }

} // wrap_cell
