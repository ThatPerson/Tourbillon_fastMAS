#ifndef HAMILTONIAN_CLASS
#define HAMILTONIAN_CLASS


#include "tensor.h"
#include "crystal.h"
#include <string>
#include <complex>


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

#endif
