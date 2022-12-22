#ifndef CRYSTAL_CLASS
#define CRYSTAL_CLASS


#include <string>
#include <vector>
#include <set>


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
    double S;
private:
    int get_nearest(int i, std::set<int> &results);
    int read_nearest(std::set<int> *results);

    // Informations about the unit cell: its dimensions, the number of nuclei it
    // contains, and the fractional cartesian coordinates of the nuclei
    double unit_cell_len[3];    // a, b and c length
    bool nearest_CS;
    std::string nn_file;

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

#endif 
