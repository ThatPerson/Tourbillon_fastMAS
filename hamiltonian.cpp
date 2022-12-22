#include "constants.h"
#include "tensor.h"
#include "crystal.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <complex>


using namespace std;


// The class Hamiltonian stores the spatial part of the
// dipolar coupling Hamiltonian for the system of interest. 
// The dipolar coupling hamiltonian can be calculated
// in several frames: from the principal axis system (calculated with 
// geometrical information about the crystal) to the lab frame.
class Hamiltonian {


public:

    // omega_D_L_20[i][j] stores the 20 component of the dipolar coupling
    // between i and j expressed in the lab frame, in rad/ms. This component,
    // which corresponds to the secular part of the hamiltonian, is purely
    // real.
    // This quantity is the coeffficient of the following operator:
    // 2 * I_izI_jz - (1/2)(I_i+I_j- + I_i-I_j+)
    double **omega_D_L_20;
    double *omega_CSA_L_20;
    double *omega_CSA_L_2p2;
    double *omega_CSA_L_2m2;
    // Constructor
    Hamiltonian(const Crystal &crystal,
                std::string output_filename,
                const double S);


    // This function performs a rotation from the Rotor frame to the Lab
    // frame, for a rotor set at the magic angle, and rotated around its
    // axis by omega_r_t
    // omega_r_t has to be given in radian
    void spin_at_magic_angle(double omega_r_t, const Crystal &crystal);

    // This functions performs a rotation from the Crystal frame to the Rotor
    // frame
    void rotate_crystallite(const Euler &euler_C_to_R, const Crystal &crystal);


private:

    int number_of_spins;        // Number of spins in the system

    double gamma;            // Gyromagnetic ratio of the nuclei under study

    // Arrays storing, for each pair of spin, the irreducible spherical tensor
    // describing the spatial part of the dipolar interaction in a given frame:
    // P: principal axis system
    // C: crystal frame
    // R: rotor frame
    // L: laboratory frame
    R2_IST **omega_D_P;
    R2_IST **omega_D_C;
    R2_IST **omega_D_R;
    R2_IST **omega_D_L;

    R2_IST *omega_CSA_P;
    R2_IST *omega_CSA_C;
    R2_IST *omega_CSA_R;
    R2_IST *omega_CSA_L;

    // Functions called by the constuctor to initialise the hamiltonian
    void calculate_omega_D(const Crystal &crystal, const double S);

    // Appends some of the information calculated by the constructor to the
    // output file
    void write_omega_D(std::string output_filename);

};


// The constructor function allocates memory for the arrays and initialises
// them to zero. It then calls functions using information from a Crystal
// to calculate the hamiltonian in the lab frame. Finally it
// appends some of the information to the output file
Hamiltonian::Hamiltonian(const Crystal &crystal,
                         string output_filename,
                         const double S) {

    number_of_spins = crystal.number_of_nuclei;

    gamma = GAMMA_H;


    // Allocates memory
    try {
        omega_D_L_20 = new double *[number_of_spins];
        omega_D_P = new R2_IST *[number_of_spins];
        omega_D_C = new R2_IST *[number_of_spins];
        omega_D_R = new R2_IST *[number_of_spins];
        omega_D_L = new R2_IST *[number_of_spins];
        omega_CSA_P = new R2_IST [number_of_spins];
        omega_CSA_C = new R2_IST [number_of_spins];
        omega_CSA_R = new R2_IST [number_of_spins];
        omega_CSA_L = new R2_IST [number_of_spins];
        omega_CSA_L_20 = new double [number_of_spins];
        omega_CSA_L_2p2 = new double [number_of_spins];
        omega_CSA_L_2m2 = new double [number_of_spins];

        for (int i = 0; i < number_of_spins; i++) {
            omega_D_L_20[i] = new double[number_of_spins];
            omega_D_P[i] = new R2_IST[number_of_spins];
            omega_D_C[i] = new R2_IST[number_of_spins];
            omega_D_R[i] = new R2_IST[number_of_spins];
            omega_D_L[i] = new R2_IST[number_of_spins];
        }

    }
    catch (bad_alloc) {
        cerr << "\nMemory could not be allocated for the Hamiltonian\n";
        exit(1);
    }


    // Initialises the arrays to zero
    for (int i = 0; i < number_of_spins; i++) {
        for (int j = 0; j < number_of_spins; j++) {
            omega_D_L_20[i][j] = 0.0;
            omega_D_P[i][j] = R2_IST();
            omega_D_C[i][j] = R2_IST();
            omega_D_R[i][j] = R2_IST();
            omega_D_L[i][j] = R2_IST();
        } // j loop
    }   // i loop


    // Calculate the elements of the hamiltonian in several frames
    calculate_omega_D(crystal, S);

    // Appends some values to the output file
    write_omega_D(output_filename);
}


// This function calculates the elements of the irreducible tensor for
// dipolar coupling between each pair of spin in several frames:
void Hamiltonian::calculate_omega_D(const Crystal &crystal, const double S) {

    // This expression assumes distances are in angstrom, and yields coupling
    // constants in 10^3 radian
    double const_term = (-1.0) * MU_0 * S * gamma * gamma * H_BAR * 1.0e30\
 / (2 * TWO_PI) / 1000.0;

	/* : Only loop over auto j : crystal.neighbours[i] */
    for (int i = 0; i < number_of_spins; i++) {
        /* parameters are in crystal.csaniso, crystal.cseta, crystal.csalpha, crystal.csbeta, crystal.csgamma
         * First; set components of omega_CSA_P.
         * Then; convert from P to C frame.
         * Then; convert from C to R
         * finally: convert from R to L.
         */

        omega_CSA_P[i].components[0] = -sqrt(2/3.) * (1/2.) * crystal.cseta[i] * crystal.csaniso_ppm[i];
        omega_CSA_P[i].components[4] = -sqrt(2/3.) * (1/2.) * crystal.cseta[i] * crystal.csaniso_ppm[i];
        omega_CSA_P[i].components[2] = crystal.csaniso_ppm[i];

        Euler euler_P_to_C(crystal.csalpha[i], crystal.csbeta[i], crystal.csgamma[i]);
        omega_CSA_C[i] = R2_IST(omega_CSA_P[i], euler_P_to_C);
        omega_CSA_R[i] = R2_IST(omega_CSA_C[i]);
        omega_CSA_L[i] = R2_IST(omega_CSA_R[i]);
        omega_CSA_L_20[i] = omega_D_L[i]->components[2].real();
        omega_CSA_L_2p2[i] = omega_D_L[i]->components[4].real();
        omega_CSA_L_2m2[i] = omega_D_L[i]->components[0].real();

        for (auto j : crystal.nearest_neighbours[i]) {
       // for (int j = 0; j < number_of_spins; j++) {

            // If i and j are equal nothing has to be done
            if (i == j)
                continue;

            // The vector separation between two spins can be expressed in spherical
            // coordinates
            double rho = crystal.nuclei_separation[0][i][j];
            double theta = crystal.nuclei_separation[1][i][j];
            double phi = crystal.nuclei_separation[2][i][j];

            double rot_avg = 1;

            /* if the two spins are in the same methyl/NH3 group, then we average the coupling.
             * this will be the case if the two spins are < 1.8 A in distance, and their CS are all identical.
             */
            if (rho < 1.8) {
                if ((crystal.csiso_ppm[i] - crystal.csiso_ppm[j]) < 0.01 &&
                        (crystal.csaniso_ppm[i] - crystal.csaniso_ppm[j]) < 0.01 &&
                        (crystal.csalpha[i] - crystal.csalpha[j]) < 0.01 &&
                        (crystal.csbeta[i] - crystal.csbeta[j]) < 0.01 &&
                        (crystal.csgamma[i] - crystal.csgamma[j]) < 0.01) {
                    rot_avg = crystal.methyl_rotation_val;
                    cout << "Methyl rotation between " << i << " and " << j << endl;
                }
            }

            // In the PAS, the 20 component of omega_D is the only non-zero
            // component (axially symmetric anisotropic tensor)
            omega_D_P[i][j].components[2] = rot_avg * const_term / pow(rho, 3);


            // The euler angle to go from the PAS to the Crystal frame can be
            // expressed as a function of the vector separation between the
            // spins
            double beta = theta;
            double gamma = PI - phi;

            // By convention, gamma belongs to [0, TWO_PI]
            if (gamma < 0.0)
                gamma += TWO_PI;

            // Define the euler angles and perform the rotation
            Euler euler_P_to_C(0.0, beta, gamma);
            omega_D_C[i][j] = R2_IST(omega_D_P[i][j], euler_P_to_C);

            // The rotation from the Crystal frame to the rotor frame depends
            // on each crystallites. By default, the orientation of the Crystal
            // and Rotor frame are identical
            omega_D_R[i][j] = R2_IST(omega_D_C[i][j]);


            // By default, the rotor is aligned along the magnetic field,
            // and the sample is not spinning, i.e. the rotor frame and the
            // lab frame are the same.
            omega_D_L[i][j] = R2_IST(omega_D_R[i][j]);

            omega_D_L_20[i][j] = omega_D_L[i][j].components[2].real();

        }    // j loop
    }    // i loop

}    // calculate_omega_D() function


// This function performs a rotation from the Rotor frame to the Lab frame,
// assuming the rotor is set at the magic angle, and rotated around its
// axis by an angle omega_r_t
void Hamiltonian::spin_at_magic_angle(double omega_r_t, const Crystal &crystal) {
	
    double magic_angle = acos(sqrt(1 / 3.0));

    Euler euler_R_to_L(-1.0 * omega_r_t, magic_angle, 0.0);

    for (int i = 0; i < number_of_spins; i++) {
        omega_CSA_L[i] = R2_IST(omega_CSA_R[i], euler_R_to_L);
        omega_CSA_L_20[i] = omega_CSA_L[i].components[2].real();
        omega_CSA_L_2p2[i] = omega_CSA_L[i].components[4].real();
        omega_CSA_L_2m2[i] = omega_CSA_L[i].components[0].real();
        for (auto j: crystal.nearest_neighbours[i]) {
        //for (int j = 0; j < number_of_spins; j++) {

            // If i and j are equal nothing has to be done
            if (i == j)
                continue;
            else {
                omega_D_L[i][j] = R2_IST(omega_D_R[i][j], euler_R_to_L);

                omega_D_L_20[i][j] = omega_D_L[i][j].components[2].real();

            }    // else command

        }    // j loop
    }    // i loop

}


// This function performs a rotation from the Crystal frame to the Rotor frame,
// given a set of Euler angle
void Hamiltonian::rotate_crystallite(const Euler &euler_C_to_R, const Crystal &crystal) {
	/* : Restrict to only j in neighbour set of i */
    for (int i = 0; i < number_of_spins; i++) {
        omega_CSA_R[i] = R2_IST(omega_CSA_C[i], euler_C_to_R);
        omega_CSA_L[i] = R2_IST(omega_CSA_R[i]);
        omega_CSA_L_20[i] = omega_CSA_L[i].components[2].real();
        omega_CSA_L_2p2[i] = omega_CSA_L[i].components[4].real();
        omega_CSA_L_2m2[i] = omega_CSA_L[i].components[0].real();
        for (auto j : crystal.nearest_neighbours[i]) {
        //for (int j = 0; j < number_of_spins; j++) {

            // If i and j are equal nothing has to be done
            if (i == j)
                continue;
            else {

                omega_D_R[i][j] = R2_IST(omega_D_C[i][j], euler_C_to_R);

                // By default, the rotor is aligned along the magnetic field,
                // and the sample is not spinning, i.e. the rotor frame and the
                // lab frame are the same.
                omega_D_L[i][j] = R2_IST(omega_D_R[i][j]);

                omega_D_L_20[i][j] = omega_D_L[i][j].components[2].real();

            }    // else command

        }    // j loop
    }    // i loop

}


void Hamiltonian::write_omega_D(string output_filename) {

    // Open an ofstream in append mode
    ofstream output_filestream(output_filename.c_str(), ios_base::app);

    // Check that the file stream opened normally
    if (output_filestream.fail()) {
        cerr << "\nOutput file " << output_filename << " could not be opened.\n";
        exit(1);
    }


    // 20 component of the dipolar coupling hamiltonian in the principal
    // axis system
    output_filestream << "omega_D_P, in kHz\n";
    for (int i = 0; i < number_of_spins; i++) {
        for (int j = 0; j < number_of_spins; j++) {

            output_filestream << setw(6) << fixed << setprecision(2) <<
                              omega_D_P[i][j].components[2].real() / (TWO_PI) << " ";

        }
        output_filestream << "\n";
    }
    output_filestream << "\n";

    output_filestream.close();

}
