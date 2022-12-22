#ifndef DENSITY_MATRIX_CLASS
#define DENSITY_MATRIX_CLASS


#include "hamiltonian.h"
#include <complex>
#include <vector>
#include <unordered_map>


// A Density_matrix stores, in a series of multidimensional arrays, the 
// coefficient b_r of the expansion of the density matrix in a basis set of
// product operators: sigma = sum_r( b_r * B_r ). 
// The following notation is used for the product operators: 
// z : I_z, p : I_+ * (1/sqrt(2)), m = I_- * (1/sqrt(2))
// The product operators have a 2^(q-1) normalisation factor, where q is the 
// number of spins; 
// e.g. zpm[i][j][k] stand for 4 * I_iz * (I_j+ / sqrt(2)) * (I_k- / sqrt(2))

// A Density_matrix has an evolve member function, that perfoms a numerical 
// integration of the LVN equation using a Suzuki-Trotter scheme
class Density_matrix {

public:

    std::unordered_map<int, std::complex<double>> z;    // Element [i] is the coefficient of I_iz

    std::unordered_map<int, std::complex<double>> pm;       // Element [i][j] is the coefficient of
    // I_i+ * I_j-

    std::unordered_map<int, std::complex<double>> zz;    // Element [i][j] is the coefficient of
    // 2 * I_iz * I_jz

    std::unordered_map<int, std::complex<double>> zpm;    // Element [i][j][k] is the coefficient of
    // 2 * I_iz * I_j+ * I_k-

    std::unordered_map<int, std::complex<double>> zzz;    // Element [i][j][k] is the coefficient of
    // 4 * I_iz * I_jz * I_kz.

    std::unordered_map<int, std::complex<double>> ppmm;    // Element [i][j][k][l] is the
    //coefficient of 2 * I_i+ * I_j+ * I_k- * I_l-

    std::unordered_map<int, std::complex<double>> zzpm;    // Element [i][j][k][l] is the
    // coefficient of 4 * I_iz * I_jz * I_k+ * I_l-
    std::vector<std::complex<double> *> saturate_z;


    //std::complex<double> ****zzzz;    // Element [i][j][k][l] is the
    // coefficient of 8 * I_iz * I_jz * I_kz * I_lz



    // Constructor function
    Density_matrix(int number_of_nuclei, Crystal &crystal);

    // Performs the numerical integration of the LVN equation with a
    // Suzuki-Trotter algorithm of order 1 or 2.
    void evolve_ST1(const Hamiltonian &hamiltonian, double timestep, const Crystal &crystal);
    void evolve_CS(const Crystal &crystal, const Hamiltonian &hamiltonian, double timestep);
    void evolve_ST2(const Hamiltonian &hamiltonian, double timestep, const Crystal &crystal);

    // sets all the elements to zero
    void reset(const Crystal &crystal);
    int number_of_spins;        // Number of spins
    void write_continuation(std::string filename);
    void read_continuation(std::string filename);
    void saturate(std::vector<int>& polarised_spins);

    void store_br(double **br_zz, double **br_pm, double **br_mp, const Crystal &crystal);

private:

    int get_hash(std::vector<int> spins);
    // The following functions are called by the main evolve function. update_R
    // calculates the matrix elements for the rotation matrices, and evolve
    // implements their action, for the dipolar Hamiltonian between spins i and j
    void update_RHH(const Hamiltonian &hamiltonian,
                    double timestep,
                    int i,
                    int j);

    void evolve_HH(int i, int j, const Crystal &crystal, const Hamiltonian &hamiltonian, double timestep);


    // The evolution makes repeated use of the same matrix elements
    double RHHzqc_11, RHHzqc_21;
    std::complex<double> RHHzqc_31, RHHzqc_41;
    double RHHpoqc_11, RHHpoqc_21;
    std::complex<double> RHHpoqc_31, RHHpoqc_41;
    double RHHmoqc_11, RHHmoqc_21;
    std::complex<double> RHHmoqc_31, RHHmoqc_41;

    // The following functions multiply the vector X = [a;b;c;d]
    // by a 4 by 4 rotation matrix
    void RHHzqcX(std::complex<double> &a, std::complex<double> &b,
                 std::complex<double> &c, std::complex<double> &d);

    void RHHpoqcX(std::complex<double> &a, std::complex<double> &b,
                  std::complex<double> &c, std::complex<double> &d);

    void RHHmoqcX(std::complex<double> &a, std::complex<double> &b,
                  std::complex<double> &c, std::complex<double> &d);

};

#endif
